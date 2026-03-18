/*
 * Copyright 2025-2026 Daniel Cederberg
 *
 * This file is part of the PSLP project (LP Presolver).
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

#include "Problem.h"
#include "Constraints.h"
#include "Matrix.h"
#include "Memory_wrapper.h"
#include "State.h"
#include "Workspace.h"
#include "utils.h"

Problem *new_problem(Constraints *constraints, Objective *obj)
{
    Problem *problem = (Problem *) ps_malloc(1, sizeof(Problem));
    RETURN_PTR_IF_NULL(problem, NULL);

    problem->constraints = constraints;
    problem->obj = obj;

    return problem;
}

void free_problem(Problem *problem)
{
    RETURN_IF_NULL(problem);
    RETURN_IF_NULL(problem->constraints);
    State *data = problem->constraints->state;

    if (data)
    {
        free_work(data->work);
        PS_FREE(data->col_locks);
        PS_FREE(data->activities);
        PS_FREE(data->row_sizes);
        PS_FREE(data->col_sizes);
        free_state(data);
    }

    free_constraints(problem->constraints);
    free_objective(problem->obj);
    PS_FREE(problem);
}

/* Helper function to transpose a matrix in CSR format */
static void transpose_matrix(const double *Ax, const int *Ai, const int *Ap,
                             double **ATx, int **ATi, int **ATp,
                             size_t n_rows, size_t n_cols, size_t nnz)
{
    *ATx = (double *) ps_malloc(nnz, sizeof(double));
    *ATi = (int *) ps_malloc(nnz, sizeof(int));
    *ATp = (int *) ps_malloc(n_cols + 1, sizeof(int));
    
    if (!*ATx || !*ATi || !*ATp) return;
    
    /* Count elements per column */
    for (size_t i = 0; i <= n_cols; i++) (*ATp)[i] = 0;
    for (size_t i = 0; i < nnz; i++) (*ATp)[Ai[i] + 1]++;
    
    /* Cumulative sum */
    for (size_t i = 0; i < n_cols; i++) (*ATp)[i + 1] += (*ATp)[i];
    
    /* Fill values */
    int *col_counts = (int *) ps_malloc(n_cols, sizeof(int));
    for (size_t i = 0; i < n_cols; i++) col_counts[i] = 0;
    
    for (size_t i = 0; i < n_rows; i++) {
        for (int j = Ap[i]; j < Ap[i + 1]; j++) {
            int col = Ai[j];
            int idx = (*ATp)[col] + col_counts[col];
            (*ATi)[idx] = (int) i;
            (*ATx)[idx] = Ax[j];
            col_counts[col]++;
        }
    }
    PS_FREE(col_counts);
}

QuadTerm *quad_term_new(const double *Px, const int *Pi, const int *Pp, size_t nnz,
                        size_t n)
{
    if (nnz == 0 || Px == NULL)
    {
        return NULL;
    }

    QuadTerm *quad = (QuadTerm *) ps_malloc(1, sizeof(QuadTerm));
    RETURN_PTR_IF_NULL(quad, NULL);

    quad->nnz = nnz;
    quad->has_quad = true;

    /* Allocate and copy P matrix data */
    quad->Px = (double *) ps_malloc(nnz, sizeof(double));
    quad->Pi = (int *) ps_malloc(nnz, sizeof(int));
    quad->Pp = (int *) ps_malloc(n + 1, sizeof(int));

    if (!quad->Px || !quad->Pi || !quad->Pp)
    {
        free_quad_term(quad);
        return NULL;
    }

    memcpy(quad->Px, Px, nnz * sizeof(double));
    memcpy(quad->Pi, Pi, nnz * sizeof(int));
    memcpy(quad->Pp, Pp, (n + 1) * sizeof(int));
    
    /* Create transpose for efficient column access */
    transpose_matrix(Px, Pi, Pp, &quad->PTx, &quad->PTi, &quad->PTp, n, n, nnz);

    return quad;
}

void free_quad_term(QuadTerm *quad)
{
    RETURN_IF_NULL(quad);
    PS_FREE(quad->Px);
    PS_FREE(quad->Pi);
    PS_FREE(quad->Pp);
    PS_FREE(quad->PTx);
    PS_FREE(quad->PTi);
    PS_FREE(quad->PTp);
    PS_FREE(quad);
}

Objective *objective_new(double *c, QuadTerm *quad)
{
    Objective *obj = (Objective *) ps_malloc(1, sizeof(Objective));
    RETURN_PTR_IF_NULL(obj, NULL);
    obj->c = c;
    obj->offset = 0.0;
    obj->quad = quad;
    obj->quad_qr = NULL;

    return obj;
}

Objective *objective_new_qr(double *c, QuadTermQR *quad_qr)
{
    Objective *obj = (Objective *) ps_malloc(1, sizeof(Objective));
    RETURN_PTR_IF_NULL(obj, NULL);
    obj->c = c;
    obj->offset = 0.0;
    obj->quad = NULL;
    obj->quad_qr = quad_qr;

    return obj;
}

void free_objective(Objective *obj)
{
    RETURN_IF_NULL(obj);
    PS_FREE(obj->c);
    free_quad_term(obj->quad);
    free_quad_term_qr(obj->quad_qr);
    PS_FREE(obj);
}

void objective_shrink(double *c, int *map, size_t len)
{
    dPtr_shrink(c, map, len);
}

void fix_var_in_obj(Objective *obj, int col, double value)
{
    obj->offset += obj->c[col] * value;
}

void fix_var_in_obj_qp(Objective *obj, int col, double value)
{
    /* Linear term contribution */
    obj->offset += obj->c[col] * value;

    /* Handle QR format: P = Q + R*R^T */
    if (obj->quad_qr != NULL && obj->quad_qr->has_quad)
    {
        fix_var_in_obj_qp_qr(obj->quad_qr, col, value, &obj->offset, obj->c);
        return;
    }

    /* Quadratic term contribution: (1/2) * P_ii * value^2 + sum_{j!=i} P_ij * value * x_j */
    if (obj->quad == NULL || !obj->quad->has_quad)
    {
        return;
    }

    QuadTerm *quad = obj->quad;
    
    /* Use PT (transpose of P) to efficiently access column 'col' of P */
    /* PT stores P in CSC format, so PT column 'col' contains all non-zeros in P row 'col' */
    if (quad->PTp != NULL) {
        int col_start = quad->PTp[col];
        int col_end = quad->PTp[col + 1];
        
        for (int idx = col_start; idx < col_end; ++idx) {
            int row = quad->PTi[idx];  /* row in original P matrix */
            double pval = quad->PTx[idx];
            
            if (row == col) {
                /* Diagonal term: (1/2) * P_ii * value^2 */
                obj->offset += 0.5 * pval * value * value;
            } else {
                /* Off-diagonal term: P_ij * value * x_j
                 * Add to c[j]: P_ij * value (symmetric, so also affects row j) */
                obj->c[row] += pval * value;
            }
        }
    }
}

/* Shrink quadratic term P by removing fixed/substituted columns
 * This function filters out rows/columns corresponding to deleted variables
 * and re-indexes the remaining entries */
static void quad_term_shrink(QuadTerm *quad, int *col_map, size_t n_cols_old)
{
    if (quad == NULL || !quad->has_quad || quad->nnz == 0)
    {
        return;
    }
    
    /* Count new columns (active columns) */
    size_t n_cols_new = 0;
    for (size_t i = 0; i < n_cols_old; i++)
    {
        if (col_map[i] >= 0)
        {
            n_cols_new++;
        }
    }
    
    /* If no columns were removed, nothing to do */
    if (n_cols_new == n_cols_old)
    {
        return;
    }
    
    if (n_cols_new == 0)
    {
        /* All columns removed - clear P matrix */
        PS_FREE(quad->Px);
        PS_FREE(quad->Pi);
        PS_FREE(quad->Pp);
        PS_FREE(quad->PTx);
        PS_FREE(quad->PTi);
        PS_FREE(quad->PTp);
        quad->nnz = 0;
        quad->has_quad = false;
        return;
    }
    
    /* First pass: count new nnz (only keeping entries where both row and col are active) */
    size_t new_nnz = 0;
    for (size_t i = 0; i < n_cols_old; i++)
    {
        if (col_map[i] < 0) continue;  /* Skip deleted rows */
        
        int row_start = quad->Pp[i];
        int row_end = quad->Pp[i + 1];
        for (int idx = row_start; idx < row_end; idx++)
        {
            int col = quad->Pi[idx];
            if (col_map[col] >= 0)  /* Keep only if column is also active */
            {
                new_nnz++;
            }
        }
    }
    
    if (new_nnz == 0)
    {
        /* No quadratic terms left */
        PS_FREE(quad->Px);
        PS_FREE(quad->Pi);
        PS_FREE(quad->Pp);
        PS_FREE(quad->PTx);
        PS_FREE(quad->PTi);
        PS_FREE(quad->PTp);
        quad->nnz = 0;
        quad->has_quad = false;
        return;
    }
    
    /* Allocate new arrays */
    double *new_Px = (double *) ps_malloc(new_nnz, sizeof(double));
    int *new_Pi = (int *) ps_malloc(new_nnz, sizeof(int));
    int *new_Pp = (int *) ps_malloc(n_cols_new + 1, sizeof(int));
    
    if (!new_Px || !new_Pi || !new_Pp)
    {
        PS_FREE(new_Px);
        PS_FREE(new_Pi);
        PS_FREE(new_Pp);
        return;
    }
    
    /* Second pass: fill new P matrix */
    new_Pp[0] = 0;
    size_t idx_new = 0;
    
    for (size_t i = 0; i < n_cols_old; i++)
    {
        if (col_map[i] < 0) continue;  /* Skip deleted rows */
        
        int new_row = col_map[i];
        int row_start = quad->Pp[i];
        int row_end = quad->Pp[i + 1];
        
        for (int idx = row_start; idx < row_end; idx++)
        {
            int col = quad->Pi[idx];
            if (col_map[col] >= 0)  /* Keep only if column is also active */
            {
                new_Px[idx_new] = quad->Px[idx];
                new_Pi[idx_new] = col_map[col];  /* New column index */
                idx_new++;
            }
        }
        new_Pp[new_row + 1] = (int) idx_new;
    }
    
    /* Free old arrays */
    PS_FREE(quad->Px);
    PS_FREE(quad->Pi);
    PS_FREE(quad->Pp);
    PS_FREE(quad->PTx);
    PS_FREE(quad->PTi);
    PS_FREE(quad->PTp);
    
    /* Assign new arrays */
    quad->Px = new_Px;
    quad->Pi = new_Pi;
    quad->Pp = new_Pp;
    quad->nnz = new_nnz;
    
    /* Rebuild transpose */
    transpose_matrix(quad->Px, quad->Pi, quad->Pp, 
                     &quad->PTx, &quad->PTi, &quad->PTp,
                     n_cols_new, n_cols_new, new_nnz);
}

void problem_clean(Problem *prob, bool remove_all)
{
    Mapping *maps = prob->constraints->state->work->mappings;
    size_t n_cols_old = prob->constraints->n;
    size_t n_rows_old = prob->constraints->m;
    constraints_clean(prob->constraints, maps, remove_all);
    clean_state(prob->constraints->state, maps, n_rows_old, n_cols_old);
    objective_shrink(prob->obj->c, maps->cols, n_cols_old);
    quad_term_shrink(prob->obj->quad, maps->cols, n_cols_old);
    quad_term_qr_shrink(prob->obj->quad_qr, maps->cols, n_cols_old);
}

void sub_var_in_obj(Objective *obj, const double *vals, const int *cols, int len,
                    int k, double aik, double rhs)
{
    if (obj->c[k] == 0.0)
    {
        return;
    }

    double ratio = obj->c[k] / aik;
    for (int i = 0; i < len; ++i)
    {
        obj->c[cols[i]] -= ratio * vals[i];
    }

    obj->offset += rhs * ratio;
}

void sub_var_in_obj_qp(Objective *obj, const double *vals, const int *cols, int len,
                       int k, double aik, double rhs)
{
    /* Substitution for QP with general P matrix is NOT implemented.
     * 
     * Reason: Variable substitution x_k = (rhs - sum a_j*x_j) / a_k in QP
     * requires complex updates to the P matrix to maintain mathematical equivalence.
     * When substituting, the quadratic term (1/2) x^T P x produces:
     * - New constant terms (contribute to offset)
     * - New linear terms (modify c)
     * - New quadratic terms between remaining variables (modify P entries)
     * 
     * The cross-terms in P require rebuilding the sparse matrix structure,
     * which is computationally expensive and complex.
     * 
     * Instead of approximate (incorrect) updates, we use a conservative strategy:
     * - Check if variable k has any quadratic terms using has_quadratic_terms()
     * - If yes, skip the substitution (return UNCHANGED)
     * - If no (variable k is "QP-free"), use standard LP substitution
     * 
     * This ensures mathematical equivalence while supporting QP presolve
     * for many practical cases (e.g., when substituted variables only appear
     * linearly in the objective).
     */
    
    /* Handle linear part using standard LP substitution.
     * This is only called when has_quadratic_terms() returned false,
     * so we know variable k has no quadratic associations.
     */
    sub_var_in_obj(obj, vals, cols, len, k, aik, rhs);
}
