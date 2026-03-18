/*
 * Copyright 2025-2026 Daniel Cederberg
 * Copyright 2026 [Your Name]
 *
 * This file is part of the PSQP project (QP Presolver).
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/* 
 * QR decomposition support: P = Q + R*R^T
 * 
 * This file contains the implementation for initializing presolver with
 * P decomposed as Q + R*R^T. The Q and R matrices are stored separately
 * and shrinked independently during presolve.
 */

#include "PSQP_API.h"
#include "Activity.h"
#include "Constraints.h"
#include "Debugger.h"
#include "Locks.h"
#include "Matrix.h"
#include "Memory_wrapper.h"
#include "Numerics.h"
#include "Problem.h"
#include "PSQP_stats.h"
#include "State.h"
#include "Timer.h"
#include "Workspace.h"
#include "glbopts.h"
#include "psqp_thread.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* Forward declaration from Presolver.c */
struct PresolveStats *init_stats(size_t n_rows, size_t n_cols, size_t nnz);

// Helper struct for parallel initialization work
typedef struct
{
    Matrix *A;
    Work *work;
    size_t n_cols;
    size_t n_rows;
    const double *lbs;
    const double *ubs;
    double *lhs_copy;
    double *rhs_copy;
    Bound *bounds;
    ColTag *col_tags;
    RowTag *row_tags;
    Lock *locks;
    Activity *activities;
    int *row_sizes;
} ParallelInitDataQR;

// thread function for initializing data
static void *init_thread_func_qr(void *arg)
{
    ParallelInitDataQR *data = (ParallelInitDataQR *) arg;

    data->row_tags = new_rowtags(data->lhs_copy, data->rhs_copy, data->n_rows);

    for (int i = 0; i < data->n_cols; i++)
    {
        data->bounds[i].lb = data->lbs[i];
        data->bounds[i].ub = data->ubs[i];

        if (IS_NEG_INF(data->lbs[i]))
        {
            UPDATE_TAG(data->col_tags[i], C_TAG_LB_INF);
        }

        if (IS_POS_INF(data->ubs[i]))
        {
            UPDATE_TAG(data->col_tags[i], C_TAG_UB_INF);
        }
    }

    data->locks = new_locks(data->A, data->row_tags);
    data->activities = new_activities(data->A, data->col_tags, data->bounds);
    count_rows(data->A, data->row_sizes);

    return NULL;
}

typedef struct clean_up_scope_qr
{
    Matrix *A, *AT;
    double *lhs_copy, *rhs_copy, *c_copy;
    int *row_sizes, *col_sizes;
    RowTag *row_tags;
    Lock *locks;
    Activity *activities;
    State *data;
    Constraints *constraints;
    Presolver *presolver;
    Objective *obj;
    ColTag *col_tags;
    Bound *bounds;
    Work *work;
} clean_up_scope_qr;

static void presolver_clean_up_qr(clean_up_scope_qr scope)
{
    free_matrix(scope.A);
    free_matrix(scope.AT);
    PS_FREE(scope.lhs_copy);
    PS_FREE(scope.rhs_copy);
    PS_FREE(scope.c_copy);
    PS_FREE(scope.row_sizes);
    PS_FREE(scope.col_sizes);
    PS_FREE(scope.row_tags);
    PS_FREE(scope.col_tags);
    PS_FREE(scope.bounds);
    free_activities(scope.activities);
    free_locks(scope.locks);
    free_work(scope.work);
    free_state(scope.data);
    free_constraints(scope.constraints);
    free_objective(scope.obj);
    free_presolver(scope.presolver);
}

/* Initialize presolver with P = Q + R*R^T decomposition
 * 
 * This function creates a presolver that stores Q and R separately.
 * The Q and R matrices are shrinked independently during presolve.
 */
Presolver *new_qp_presolver_qr(const double *Ax, const int *Ai, const int *Ap,
                               size_t m, size_t n, size_t nnz, const double *lhs,
                               const double *rhs, const double *lbs, const double *ubs,
                               const double *c,
                               const double *Qx, const int *Qi, const int *Qp, size_t Qnnz,
                               const double *Rx, const int *Ri, const int *Rp, size_t Rnnz,
                               size_t k, const Settings *stgs)
{
    Timer timer;
    clock_gettime(CLOCK_MONOTONIC, &timer.start);
    Matrix *A = NULL, *AT = NULL;
    double *lhs_copy = NULL, *rhs_copy = NULL, *c_copy = NULL;
    int *row_sizes = NULL, *col_sizes = NULL;
    RowTag *row_tags = NULL;
    Lock *locks = NULL;
    Activity *activities = NULL;
    State *data = NULL;
    Constraints *constraints = NULL;
    Presolver *presolver = NULL;
    Objective *obj = NULL;
    ColTag *col_tags = NULL;
    Bound *bounds = NULL;
    Work *work = NULL;
    QuadTermQR *quad_qr = NULL;

    /* Validate input: at least one of Q or R must be provided */
    if ((Qx == NULL || Qnnz == 0) && (Rx == NULL || Rnnz == 0))
    {
        return NULL;
    }

    //  ---------------------------------------------------------------------------
    //   Copy data and allocate memory.
    //  ---------------------------------------------------------------------------
    size_t n_rows = m;
    size_t n_cols = n;
    lhs_copy = (double *) ps_malloc(n_rows, sizeof(double));
    rhs_copy = (double *) ps_malloc(n_rows, sizeof(double));
    c_copy = (double *) ps_malloc(n_cols, sizeof(double));
    col_tags = (ColTag *) ps_calloc(n_cols, sizeof(ColTag));
    bounds = (Bound *) ps_malloc(n_cols, sizeof(Bound));
    work = new_work(n_rows, n_cols);
    row_sizes = (int *) ps_malloc(n_rows, sizeof(int));
    col_sizes = (int *) ps_malloc(n_cols, sizeof(int));

    if (!lhs_copy || !rhs_copy || !c_copy || !col_tags || !bounds || !work ||
        !row_sizes || !col_sizes)
    {
        goto cleanup;
    }

    memcpy(lhs_copy, lhs, n_rows * sizeof(double));
    memcpy(rhs_copy, rhs, n_rows * sizeof(double));
    memcpy(c_copy, c, n_cols * sizeof(double));

    // ---------------------------------------------------------------------------
    //  Create quadratic term QR if provided
    // ---------------------------------------------------------------------------
    quad_qr = quad_term_qr_new(Qx, Qi, Qp, Qnnz, Rx, Ri, Rp, Rnnz, n_cols, k);
    if (!quad_qr) goto cleanup;

    // ---------------------------------------------------------------------------
    //  Build bounds, row tags, A and AT.
    // ---------------------------------------------------------------------------
    A = matrix_new_no_extra_space(Ax, Ai, Ap, n_rows, n_cols, nnz);
    if (!A) goto cleanup;

    ps_thread_t thread_id;
    ParallelInitDataQR parallel_data = {A,    work,     n_cols,   n_rows,   lbs,
                                        ubs,  lhs_copy, rhs_copy, bounds,   col_tags,
                                        NULL, NULL,     NULL,     row_sizes};

    ps_thread_create(&thread_id, NULL, init_thread_func_qr, &parallel_data);

    // Main thread: Transpose A and count rows
    AT = transpose(A, work->iwork_n_cols);
    if (!AT)
    {
        ps_thread_join(&thread_id, NULL);
        goto cleanup;
    }
    count_rows(AT, col_sizes);

    // sync threads
    ps_thread_join(&thread_id, NULL);

    row_tags = parallel_data.row_tags;
    if (!row_tags) goto cleanup;
    locks = parallel_data.locks;
    activities = parallel_data.activities;
    if (!locks || !activities) goto cleanup;

    // ---------------------------------------------------------------------------
    //  Initialize internal data and constraints
    // ---------------------------------------------------------------------------
    data = new_state(row_sizes, col_sizes, locks, n_rows, n_cols, activities, work,
                     row_tags);

    if (!data) goto cleanup;
    constraints =
        constraints_new(A, AT, lhs_copy, rhs_copy, bounds, data, row_tags, col_tags);
    if (!constraints) goto cleanup;

    // ---------------------------------------------------------------------------
    //             Allocate the actual presolver
    // ---------------------------------------------------------------------------
    presolver = (Presolver *) ps_malloc(1, sizeof(Presolver));
    obj = objective_new_qr(c_copy, quad_qr);  /* QP with QR decomposition */
    if (!presolver || !obj) goto cleanup;
    presolver->stgs = stgs;
    presolver->prob = new_problem(constraints, obj);
    presolver->stats = init_stats(A->m, A->n, nnz);
    presolver->stats->nnz_removed_trivial = nnz - (size_t) A->nnz; /* explicit 0's */
    presolver->reduced_prob =
        (PresolvedProblem *) ps_calloc(1, sizeof(PresolvedProblem));
    DEBUG(run_debugger(constraints, false));

    // ---------------------------------------------------------------------------
    //           Allocate space for returning the solution
    // ---------------------------------------------------------------------------
    presolver->sol = (Solution *) ps_malloc(1, sizeof(Solution));
    if (!presolver->sol) goto cleanup;
    presolver->sol->x = (double *) ps_malloc(n_cols, sizeof(double));
    presolver->sol->y = (double *) ps_malloc(n_rows, sizeof(double));
    presolver->sol->z = (double *) ps_malloc(n_cols, sizeof(double));
    presolver->sol->dim_x = n_cols;
    presolver->sol->dim_y = n_rows;
    if (!presolver->sol->x || !presolver->sol->y || !presolver->sol->z)
    {
        goto cleanup;
    }

    clock_gettime(CLOCK_MONOTONIC, &timer.end);
    presolver->stats->time_init = GET_ELAPSED_SECONDS(timer);

    return presolver;

cleanup:
{
    struct clean_up_scope_qr scope = {
        A,         AT,       lhs_copy, rhs_copy,   c_copy, row_sizes,
        col_sizes, row_tags, locks,    activities, data,   constraints,
        presolver, obj,      col_tags, bounds,     work};
    presolver_clean_up_qr(scope);
    free_quad_term_qr(quad_qr);
}
    return NULL;
}
