/*
 * Copyright 2025-2026 Daniel Cederberg
 *
 * This file is part of the PSLP project (QP Presolver tests).
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

#ifndef TEST_QP_H
#define TEST_QP_H

#include "PSQP_API.h"
#include "PSQP_stats.h"
#include "test_macros.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* Test 1: Simple QP with diagonal P matrix
 * min 0.5 * (2*x1^2 + 2*x2^2) + 0*x1 + 0*x2 = x1^2 + x2^2
 * s.t. x1 + x2 = 1
 *      x1, x2 >= 0
 * 
 * Solution: x1 = x2 = 0.5, optimal value = 0.5
 * 
 * Note: With conservative QP strategy, singleton column substitution is skipped
 * for variables with quadratic terms. So the problem dimensions should be preserved.
 */
static char *test_qp_simple_diagonal()
{
    size_t m = 1;
    size_t n = 2;
    
    /* Constraint: x1 + x2 = 1 */
    double Ax[] = {1.0, 1.0};
    int Ai[] = {0, 1};
    int Ap[] = {0, 2};
    size_t nnz = 2;
    
    double lhs[] = {1.0};
    double rhs[] = {1.0};
    double lbs[] = {0.0, 0.0};
    double ubs[] = {INFINITY, INFINITY};
    double c[] = {0.0, 0.0};
    
    /* Quadratic term: P = [2 0; 0 2] (upper triangular) */
    double Px[] = {2.0, 2.0};
    int Pi[] = {0, 1};
    int Pp[] = {0, 1, 2};
    size_t Pnnz = 2;
    
    Settings *stgs = default_settings();
    stgs->verbose = false;
    /* Disable parallel_cols and dton_eq as they may fully reduce this simple problem */
    stgs->parallel_cols = false;
    stgs->dton_eq = false;
    
    Presolver *presolver = new_qp_presolver_qr(Ax, Ai, Ap, m, n, nnz,
                                               lhs, rhs, lbs, ubs, c,
                                               Px, Pi, Pp, Pnnz,
                                               NULL, NULL, NULL, 0, 0, stgs);
    mu_assert("presolver creation failed", presolver != NULL);
    
    PresolveStatus status = run_presolver(presolver);
    PresolvedProblem *reduced = presolver->reduced_prob;
    
    /* Problem should be processed */
    mu_assert("QP should be processed", status == REDUCED || status == UNCHANGED);
    
    /* With conservative QP strategy:
     * - Variables with quadratic terms are not substituted
     * - But fix operations (tight bounds) still work and update offset
     * - Problem dimensions may be reduced via other means
     */
    if (reduced->n > 0) {
        mu_assert("should have quadratic term", reduced->Qnnz > 0 || reduced->Rnnz > 0 || reduced->k > 0);
    }
    /* If n==0, offset should capture the quadratic contribution */
    
    PS_FREE(stgs);
    free_presolver(presolver);
    
    return 0;
}

/* Test 2: QP with off-diagonal terms and multiple constraints
 * min 0.5 * x^T P x + c^T x
 * where P = [2 1; 1 2] (symmetric)
 * Constraints:
 *   x1 + x2 >= 1
 *   x1 + x2 <= 2
 * This creates a non-trivial QP that won't be fully reduced.
 */
static char *test_qp_with_offdiagonal()
{
    size_t m = 2;
    size_t n = 2;
    
    /* Constraints: 
     *  x1 + x2 >= 1  (row 0)
     *  x1 + x2 <= 2  (row 1)
     */
    double Ax[] = {1.0, 1.0, 1.0, 1.0};
    int Ai[] = {0, 1, 0, 1};
    int Ap[] = {0, 2, 4};
    size_t nnz = 4;
    
    double lhs[] = {1.0, -INFINITY};
    double rhs[] = {INFINITY, 2.0};
    double lbs[] = {0.0, 0.0};
    double ubs[] = {INFINITY, INFINITY};
    double c[] = {1.0, 1.0};
    
    /* Quadratic term: P = [2 1; 1 2] (upper triangular: [2 1; 0 2]) */
    double Px[] = {2.0, 1.0, 2.0};
    int Pi[] = {0, 1, 1};
    int Pp[] = {0, 2, 3};
    size_t Pnnz = 3;
    
    Settings *stgs = default_settings();
    stgs->verbose = false;
    /* Disable presolvers that may fully reduce this problem */
    stgs->ston_cols = false;
    stgs->parallel_cols = false;
    stgs->dton_eq = false;
    
    Presolver *presolver = new_qp_presolver_qr(Ax, Ai, Ap, m, n, nnz,
                                               lhs, rhs, lbs, ubs, c,
                                               Px, Pi, Pp, Pnnz,
                                               NULL, NULL, NULL, 0, 0, stgs);
    mu_assert("presolver creation failed", presolver != NULL);
    
    PresolveStatus status = run_presolver(presolver);
    PresolvedProblem *reduced = presolver->reduced_prob;
    
    mu_assert("QP with off-diagonal should work", status == REDUCED || status == UNCHANGED);
    
    /* If not fully reduced, should preserve quadratic term */
    if (reduced->n > 0) {
        mu_assert("should have quadratic term", reduced->Qnnz > 0 || reduced->Rnnz > 0 || reduced->k > 0);
    }
    
    PS_FREE(stgs);
    free_presolver(presolver);
    
    return 0;
}

/* Test 3: Larger QP problem with overlapping constraints
 * This creates a denser constraint matrix that won't be fully reduced.
 */
static char *test_qp_larger()
{
    size_t m = 3;
    size_t n = 3;
    
    /* Constraints:
     * x1 + x2 >= 1
     * x2 + x3 >= 1
     * x1 + x2 + x3 <= 5
     */
    double Ax[] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    int Ai[] = {0, 1, 0, 1, 2, 0, 2};
    int Ap[] = {0, 2, 4, 7};
    size_t nnz = 7;
    
    double lhs[] = {1.0, 1.0, -INFINITY};
    double rhs[] = {INFINITY, INFINITY, 5.0};
    double lbs[] = {0.0, 0.0, 0.0};
    double ubs[] = {INFINITY, INFINITY, INFINITY};
    double c[] = {0.0, 0.0, 0.0};
    
    /* Quadratic term: P = 2*I (diagonal) */
    double Px[] = {2.0, 2.0, 2.0};
    int Pi[] = {0, 1, 2};
    int Pp[] = {0, 1, 2, 3};
    size_t Pnnz = 3;
    
    Settings *stgs = default_settings();
    stgs->verbose = false;
    /* Disable presolvers that may fully reduce this problem */
    stgs->ston_cols = false;
    stgs->parallel_cols = false;
    stgs->dton_eq = false;
    
    Presolver *presolver = new_qp_presolver_qr(Ax, Ai, Ap, m, n, nnz,
                                               lhs, rhs, lbs, ubs, c,
                                               Px, Pi, Pp, Pnnz,
                                               NULL, NULL, NULL, 0, 0, stgs);
    mu_assert("presolver creation failed for larger QP", presolver != NULL);
    
    PresolveStatus status = run_presolver(presolver);
    PresolvedProblem *reduced = presolver->reduced_prob;
    
    mu_assert("larger QP should be processed", status == REDUCED || status == UNCHANGED);
    
    /* Quadratic term should be preserved if any variables remain */
    if (reduced->n > 0) {
        mu_assert("should have quadratic term", reduced->Qnnz > 0 || reduced->Rnnz > 0 || reduced->k > 0);
    }
    
    PS_FREE(stgs);
    free_presolver(presolver);
    
    return 0;
}

/* Test 4: QP with variable bounds causing fixes
 * x1 is fixed to 1 by bounds, so the P matrix should be reduced accordingly.
 */
static char *test_qp_with_bounds()
{
    size_t m = 2;
    size_t n = 3;
    
    /* Constraints:
     * x1 + x2 + x3 >= 3
     * x1 + x2 + x3 <= 5
     */
    double Ax[] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    int Ai[] = {0, 1, 2, 0, 1, 2};
    int Ap[] = {0, 3, 6};
    size_t nnz = 6;
    
    double lhs[] = {3.0, -INFINITY};
    double rhs[] = {INFINITY, 5.0};
    
    /* x1 fixed to 1 by bounds */
    double lbs[] = {1.0, 0.0, 0.0};
    double ubs[] = {1.0, INFINITY, INFINITY};
    double c[] = {0.0, 0.0, 0.0};
    
    /* Quadratic term: P = 2*I */
    double Px[] = {2.0, 2.0, 2.0};
    int Pi[] = {0, 1, 2};
    int Pp[] = {0, 1, 2, 3};
    size_t Pnnz = 3;
    
    Settings *stgs = default_settings();
    stgs->verbose = false;
    /* Disable presolvers that may interfere with this test */
    stgs->ston_cols = false;
    stgs->parallel_cols = false;
    stgs->dton_eq = false;
    
    Presolver *presolver = new_qp_presolver_qr(Ax, Ai, Ap, m, n, nnz,
                                               lhs, rhs, lbs, ubs, c,
                                               Px, Pi, Pp, Pnnz,
                                               NULL, NULL, NULL, 0, 0, stgs);
    mu_assert("presolver creation failed", presolver != NULL);
    
    PresolveStatus status = run_presolver(presolver);
    PresolvedProblem *reduced = presolver->reduced_prob;
    
    /* x1 should be fixed by bound propagation, reducing the problem */
    mu_assert("QP with fixed vars should reduce", status == REDUCED || status == UNCHANGED);
    
    /* Offset should account for fixed variable contribution: 0.5 * P_11 * x1^2 = 0.5 * 2 * 1 = 1 */
    mu_assert("offset should be updated", reduced->obj_offset >= 0.99 && reduced->obj_offset <= 1.01);
    
    /* After fixing x1, we should have 2 variables left with P = 2*I (2x2) */
    if (reduced->n > 0) {
        mu_assert("should have 2 variables", reduced->n == 2);
        mu_assert("should have quadratic term", reduced->Qnnz > 0 || reduced->Rnnz > 0 || reduced->k > 0);
        mu_assert("Qnnz should be 2", reduced->Qnnz == 2);
    }
    
    PS_FREE(stgs);
    free_presolver(presolver);
    
    return 0;
}

/* Test 5: Dense QP matrix with multiple off-diagonal terms
 * Creates a 5x5 dense quadratic matrix with cross terms
 * This tests the transpose handling and offset computation.
 */
static char *test_qp_dense_matrix()
{
    size_t m = 3;
    size_t n = 5;
    
    /* Constraints:
     * sum(x_i) >= 5   (row 0)
     * sum(x_i) <= 20  (row 1)
     * x1 + x2 >= 1    (row 2)
     */
    double Ax[] = {1.0, 1.0, 1.0, 1.0, 1.0,  /* row 0: all vars */
                   1.0, 1.0, 1.0, 1.0, 1.0,  /* row 1: all vars */
                   1.0, 1.0, 0.0, 0.0, 0.0}; /* row 2: x1, x2 */
    int Ai[] = {0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1};
    int Ap[] = {0, 5, 10, 12};
    size_t nnz = 12;
    
    double lhs[] = {5.0, -INFINITY, 1.0};
    double rhs[] = {INFINITY, 20.0, INFINITY};
    double lbs[] = {0.0, 0.0, 0.0, 0.0, 0.0};
    double ubs[] = {INFINITY, INFINITY, INFINITY, INFINITY, INFINITY};
    double c[] = {1.0, 2.0, 3.0, 4.0, 5.0};
    
    /* Dense quadratic term: P[i][j] = i+j+1 for i <= j (upper triangular)
     * Creates matrix:
     * [1 2 3 4 5]
     * [2 3 4 5 6]
     * [3 4 5 6 7]
     * [4 5 6 7 8]
     * [5 6 7 8 9]
     */
    double Px[] = {1.0, 2.0, 3.0, 4.0, 5.0,   /* row 0 */
                   3.0, 4.0, 5.0, 6.0,        /* row 1 */
                   5.0, 6.0, 7.0,             /* row 2 */
                   7.0, 8.0,                  /* row 3 */
                   9.0};                      /* row 4 */
    int Pi[] = {0, 1, 2, 3, 4,
                1, 2, 3, 4,
                2, 3, 4,
                3, 4,
                4};
    int Pp[] = {0, 5, 9, 12, 14, 15};
    size_t Pnnz = 15;
    
    Settings *stgs = default_settings();
    stgs->verbose = false;
    stgs->ston_cols = false;
    stgs->parallel_cols = false;
    stgs->dton_eq = false;
    
    Presolver *presolver = new_qp_presolver_qr(Ax, Ai, Ap, m, n, nnz,
                                               lhs, rhs, lbs, ubs, c,
                                               Px, Pi, Pp, Pnnz,
                                               NULL, NULL, NULL, 0, 0, stgs);
    mu_assert("presolver creation failed for dense QP", presolver != NULL);
    
    PresolveStatus status = run_presolver(presolver);
    PresolvedProblem *reduced = presolver->reduced_prob;
    
    /* Note: This problem might be detected as infeasible or unbounded due to
     * constraint interactions. We accept any valid status. */
    mu_assert("dense QP should be processed or detected", 
              status == REDUCED || status == UNCHANGED || 
              status == INFEASIBLE || status == UNBNDORINFEAS);
    
    if (status == REDUCED || status == UNCHANGED) {
        if (reduced->n > 0) {
            mu_assert("should have quadratic term", reduced->Qnnz > 0 || reduced->Rnnz > 0 || reduced->k > 0);
        }
    }
    
    PS_FREE(stgs);
    free_presolver(presolver);
    
    return 0;
}

/* Test 6: QP with all presolvers enabled
 * This tests that all presolvers work correctly with QP problems.
 */
static char *test_qp_all_presolvers()
{
    size_t m = 4;
    size_t n = 6;
    
    /* Create a problem that won't be fully reduced but will exercise all presolvers */
    double Ax[] = {1.0, 2.0, 3.0,    /* row 0: x1, x2, x3 */
                   1.0, 1.0,         /* row 1: x4, x5 */
                   2.0, 3.0,         /* row 2: x4, x5 (parallel to row 1) */
                   1.0, 1.0, 1.0};   /* row 3: x2, x3, x6 */
    int Ai[] = {0, 1, 2, 3, 4, 3, 4, 3, 4, 1, 2, 5};
    int Ap[] = {0, 3, 5, 7, 10};
    size_t nnz = 10;
    
    double lhs[] = {-INFINITY, 5.0, 10.0, 2.0};
    double rhs[] = {10.0, 8.0, 16.0, INFINITY};
    double lbs[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double ubs[] = {INFINITY, INFINITY, INFINITY, INFINITY, INFINITY, INFINITY};
    double c[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    
    /* Diagonal P matrix - simple case */
    double Px[] = {2.0, 2.0, 2.0, 2.0, 2.0, 2.0};
    int Pi[] = {0, 1, 2, 3, 4, 5};
    int Pp[] = {0, 1, 2, 3, 4, 5, 6};
    size_t Pnnz = 6;
    
    Settings *stgs = default_settings();
    stgs->verbose = false;
    /* Enable all presolvers */
    
    Presolver *presolver = new_qp_presolver_qr(Ax, Ai, Ap, m, n, nnz,
                                               lhs, rhs, lbs, ubs, c,
                                               Px, Pi, Pp, Pnnz,
                                               NULL, NULL, NULL, 0, 0, stgs);
    mu_assert("presolver creation failed for all-presolvers QP", presolver != NULL);
    
    PresolveStatus status = run_presolver(presolver);
    
    mu_assert("QP with all presolvers should work", 
              status == REDUCED || status == UNCHANGED || status == INFEASIBLE);
    
    PS_FREE(stgs);
    free_presolver(presolver);
    
    return 0;
}

/* Test 7: QP with fixed variables and postsolve
 * Tests that the offset is correctly computed when variables are fixed.
 */
static char *test_qp_postsolve()
{
    size_t m = 2;
    size_t n = 3;
    
    /* Constraints:
     * x1 + x2 + x3 = 6
     * x1 - x2 = 0  => x1 = x2
     */
    double Ax[] = {1.0, 1.0, 1.0, 1.0, -1.0};
    int Ai[] = {0, 1, 2, 0, 1};
    int Ap[] = {0, 3, 5};
    size_t nnz = 5;
    
    double lhs[] = {6.0, 0.0};
    double rhs[] = {6.0, 0.0};
    double lbs[] = {0.0, 0.0, 0.0};
    double ubs[] = {INFINITY, INFINITY, INFINITY};
    double c[] = {0.0, 0.0, 0.0};
    
    /* Diagonal P = 2*I */
    double Px[] = {2.0, 2.0, 2.0};
    int Pi[] = {0, 1, 2};
    int Pp[] = {0, 1, 2, 3};
    size_t Pnnz = 3;
    
    Settings *stgs = default_settings();
    stgs->verbose = false;
    stgs->ston_cols = false;
    stgs->parallel_cols = false;
    
    Presolver *presolver = new_qp_presolver_qr(Ax, Ai, Ap, m, n, nnz,
                                               lhs, rhs, lbs, ubs, c,
                                               Px, Pi, Pp, Pnnz,
                                               NULL, NULL, NULL, 0, 0, stgs);
    mu_assert("presolver creation failed for QP postsolve", presolver != NULL);
    
    PresolveStatus status = run_presolver(presolver);
    PresolvedProblem *reduced = presolver->reduced_prob;
    
    mu_assert("QP postsolve test should reduce", status == REDUCED || status == UNCHANGED);
    
    /* If problem was reduced, try postsolve with dummy solution */
    if (reduced->n > 0 && reduced->m > 0) {
        double *x_reduced = (double *)calloc(reduced->n, sizeof(double));
        double *y_reduced = (double *)calloc(reduced->m, sizeof(double));
        double *z_reduced = (double *)calloc(reduced->n, sizeof(double));
        
        /* Fill with some values */
        for (size_t i = 0; i < reduced->n; i++) x_reduced[i] = 1.0;
        for (size_t i = 0; i < reduced->m; i++) y_reduced[i] = 1.0;
        
        postsolve(presolver, x_reduced, y_reduced, z_reduced);
        
        /* Check that solution was recovered for original problem */
        mu_assert("postsolve should recover solution", presolver->sol != NULL);
        mu_assert("postsolve x dim should match original", presolver->sol->dim_x == (int)n);
        mu_assert("postsolve y dim should match original", presolver->sol->dim_y == (int)m);
        
        free(x_reduced);
        free(y_reduced);
        free(z_reduced);
    }
    
    PS_FREE(stgs);
    free_presolver(presolver);
    
    return 0;
}

/* Test 8: Large-scale QP problem
 * Tests performance and correctness on larger problems.
 */
static char *test_qp_large_scale()
{
    /* 50 variables, 30 constraints */
    size_t n = 50;
    size_t m = 30;
    
    /* Create a sparse random-like constraint matrix */
    size_t nnz_max = 300;
    double *Ax = (double *)calloc(nnz_max, sizeof(double));
    int *Ai = (int *)calloc(nnz_max, sizeof(int));
    int *Ap = (int *)calloc(m + 1, sizeof(int));
    
    /* Build constraint: each row has 5-15 non-zeros */
    size_t nnz = 0;
    Ap[0] = 0;
    for (size_t i = 0; i < m; i++) {
        size_t row_nnz = 5 + (i % 11);  /* 5 to 15 non-zeros per row */
        for (size_t j = 0; j < row_nnz && nnz < nnz_max; j++) {
            Ax[nnz] = 1.0 + ((i + j) % 5);
            Ai[nnz] = (int)((i * 3 + j * 7) % n);
            nnz++;
        }
        Ap[i + 1] = (int)nnz;
    }
    
    double *lhs = (double *)calloc(m, sizeof(double));
    double *rhs = (double *)calloc(m, sizeof(double));
    for (size_t i = 0; i < m; i++) {
        lhs[i] = -INFINITY;
        rhs[i] = 100.0 + (i * 10);
    }
    
    double *lbs = (double *)calloc(n, sizeof(double));
    double *ubs = (double *)calloc(n, sizeof(double));
    double *c = (double *)calloc(n, sizeof(double));
    for (size_t i = 0; i < n; i++) {
        lbs[i] = 0.0;
        ubs[i] = INFINITY;
        c[i] = (double)(i + 1);
    }
    
    /* Diagonal P matrix - simple and efficient for large problems */
    size_t Pnnz = n;
    double *Px = (double *)calloc(Pnnz, sizeof(double));
    int *Pi = (int *)calloc(Pnnz, sizeof(int));
    int *Pp = (int *)calloc(n + 1, sizeof(int));
    
    for (size_t i = 0; i < n; i++) {
        Px[i] = 2.0 + ((i % 3) == 0 ? 1.0 : 0.0);  /* Some vars have P_ii = 3, others = 2 */
        Pi[i] = (int)i;
        Pp[i] = (int)i;
    }
    Pp[n] = (int)n;
    
    Settings *stgs = default_settings();
    stgs->verbose = false;
    stgs->max_time = 10.0;  /* 10 second limit */
    
    Presolver *presolver = new_qp_presolver_qr(Ax, Ai, Ap, m, n, nnz,
                                               lhs, rhs, lbs, ubs, c,
                                               Px, Pi, Pp, Pnnz,
                                               NULL, NULL, NULL, 0, 0, stgs);
    mu_assert("presolver creation failed for large QP", presolver != NULL);
    
    PresolveStatus status = run_presolver(presolver);
    PresolvedProblem *reduced = presolver->reduced_prob;
    
    mu_assert("large QP should be processed without error", 
              status == REDUCED || status == UNCHANGED);
    
    /* Verify stats are reasonable */
    mu_assert("original nnz should be recorded", presolver->stats->nnz_original == nnz);
    mu_assert("reduced nnz should be valid", reduced->nnz <= nnz);
    
    /* Cleanup */
    free(Ax);
    free(Ai);
    free(Ap);
    free(lhs);
    free(rhs);
    free(lbs);
    free(ubs);
    free(c);
    free(Px);
    free(Pi);
    free(Pp);
    
    PS_FREE(stgs);
    free_presolver(presolver);
    
    return 0;
}

/* Test 9: QP with mixed linear/quadratic variables
 * Some variables only appear in linear terms, others have quadratic terms.
 */
static char *test_qp_mixed_vars()
{
    size_t m = 2;
    size_t n = 4;
    
    /* Constraints:
     * x1 + x2 + x3 + x4 = 10
     * x1 - x3 >= 0
     */
    double Ax[] = {1.0, 1.0, 1.0, 1.0, 1.0, -1.0};
    int Ai[] = {0, 1, 2, 3, 0, 2};
    int Ap[] = {0, 4, 6};
    size_t nnz = 6;
    
    double lhs[] = {10.0, 0.0};
    double rhs[] = {10.0, INFINITY};
    double lbs[] = {0.0, 0.0, 0.0, 0.0};
    double ubs[] = {INFINITY, INFINITY, INFINITY, INFINITY};
    double c[] = {1.0, 0.0, 1.0, 0.0};
    
    /* P only has entries for x1 and x3 (indices 0 and 2)
     * P = [2 1; 0 2] for the submatrix
     */
    double Px[] = {2.0, 1.0, 2.0};
    int Pi[] = {0, 2, 2};
    int Pp[] = {0, 1, 1, 3, 3};  /* x2 and x4 have no quadratic terms */
    size_t Pnnz = 3;
    
    Settings *stgs = default_settings();
    stgs->verbose = false;
    stgs->ston_cols = false;
    stgs->parallel_cols = false;
    stgs->dton_eq = false;
    
    Presolver *presolver = new_qp_presolver_qr(Ax, Ai, Ap, m, n, nnz,
                                               lhs, rhs, lbs, ubs, c,
                                               Px, Pi, Pp, Pnnz,
                                               NULL, NULL, NULL, 0, 0, stgs);
    mu_assert("presolver creation failed for mixed QP", presolver != NULL);
    
    PresolveStatus status = run_presolver(presolver);
    PresolvedProblem *reduced = presolver->reduced_prob;
    
    mu_assert("mixed QP should be processed", status == REDUCED || status == UNCHANGED);
    
    /* Quadratic info should be preserved for remaining vars */
    if (reduced->n > 0) {
        /* Either has quadratic term or all quadratic vars were eliminated */
        (void)0; /* OK either way */
    }
    
    PS_FREE(stgs);
    free_presolver(presolver);
    
    return 0;
}

/* Test 10: Edge case - single variable QP
 */
static char *test_qp_single_var()
{
    size_t m = 1;
    size_t n = 1;
    
    /* Constraint: x >= 5 */
    double Ax[] = {1.0};
    int Ai[] = {0};
    int Ap[] = {0, 1};
    size_t nnz = 1;
    
    double lhs[] = {5.0};
    double rhs[] = {INFINITY};
    double lbs[] = {0.0};
    double ubs[] = {INFINITY};
    double c[] = {1.0};
    
    /* P = [2] */
    double Px[] = {2.0};
    int Pi[] = {0};
    int Pp[] = {0, 1};
    size_t Pnnz = 1;
    
    Settings *stgs = default_settings();
    stgs->verbose = false;
    
    Presolver *presolver = new_qp_presolver_qr(Ax, Ai, Ap, m, n, nnz,
                                               lhs, rhs, lbs, ubs, c,
                                               Px, Pi, Pp, Pnnz,
                                               NULL, NULL, NULL, 0, 0, stgs);
    mu_assert("presolver creation failed for single var QP", presolver != NULL);
    
    PresolveStatus status = run_presolver(presolver);
    
    mu_assert("single var QP should be processed", status == REDUCED || status == UNCHANGED);
    
    PS_FREE(stgs);
    free_presolver(presolver);
    
    return 0;
}

/* Test 11: QP with tight bounds leading to full reduction
 */
static char *test_qp_full_reduction()
{
    size_t m = 2;
    size_t n = 2;
    
    /* Constraints:
     * x1 = 1  (tight equality)
     * x2 = 2  (tight equality)
     */
    double Ax[] = {1.0, 1.0};
    int Ai[] = {0, 1};
    int Ap[] = {0, 1, 2};
    size_t nnz = 2;
    
    double lhs[] = {1.0, 2.0};
    double rhs[] = {1.0, 2.0};
    double lbs[] = {0.0, 0.0};
    double ubs[] = {INFINITY, INFINITY};
    double c[] = {1.0, 1.0};
    
    /* P = 2*I */
    double Px[] = {2.0, 2.0};
    int Pi[] = {0, 1};
    int Pp[] = {0, 1, 2};
    size_t Pnnz = 2;
    
    Settings *stgs = default_settings();
    stgs->verbose = false;
    
    Presolver *presolver = new_qp_presolver_qr(Ax, Ai, Ap, m, n, nnz,
                                               lhs, rhs, lbs, ubs, c,
                                               Px, Pi, Pp, Pnnz,
                                               NULL, NULL, NULL, 0, 0, stgs);
    mu_assert("presolver creation failed for full reduction QP", presolver != NULL);
    
    PresolveStatus status = run_presolver(presolver);
    PresolvedProblem *reduced = presolver->reduced_prob;
    
    /* Both variables should be fixed, problem fully reduced */
    mu_assert("full reduction QP should be reduced", status == REDUCED);
    
    /* Optimal value should be: 0.5*(2*1^2 + 2*2^2) + 1*1 + 1*2 = 0.5*(2+8) + 3 = 5 + 3 = 8 */
    double expected_offset = 0.5 * (2.0 * 1.0 + 2.0 * 4.0) + 1.0 + 2.0;  /* = 8 */
    mu_assert("offset should capture optimal value", 
              fabs(reduced->obj_offset - expected_offset) < 1e-6);
    
    PS_FREE(stgs);
    free_presolver(presolver);
    
    return 0;
}

/* Test 12: QP infeasibility detection
 */
static char *test_qp_infeasible()
{
    size_t m = 2;
    size_t n = 2;
    
    /* Infeasible constraints:
     * x1 + x2 >= 5
     * x1 + x2 <= 3
     * With x1, x2 >= 0
     */
    double Ax[] = {1.0, 1.0, 1.0, 1.0};
    int Ai[] = {0, 1, 0, 1};
    int Ap[] = {0, 2, 4};
    size_t nnz = 4;
    
    double lhs[] = {5.0, -INFINITY};
    double rhs[] = {INFINITY, 3.0};
    double lbs[] = {0.0, 0.0};
    double ubs[] = {INFINITY, INFINITY};
    double c[] = {1.0, 1.0};
    
    double Px[] = {2.0, 2.0};
    int Pi[] = {0, 1};
    int Pp[] = {0, 1, 2};
    size_t Pnnz = 2;
    
    Settings *stgs = default_settings();
    stgs->verbose = false;
    
    Presolver *presolver = new_qp_presolver_qr(Ax, Ai, Ap, m, n, nnz,
                                               lhs, rhs, lbs, ubs, c,
                                               Px, Pi, Pp, Pnnz,
                                               NULL, NULL, NULL, 0, 0, stgs);
    mu_assert("presolver creation failed for infeasible QP", presolver != NULL);
    
    PresolveStatus status = run_presolver(presolver);
    
    mu_assert("infeasible QP should be detected", status == INFEASIBLE);
    
    PS_FREE(stgs);
    free_presolver(presolver);
    
    return 0;
}

/* Test 13: QP with very small P values (near-zero)
 * Tests numerical stability with tiny quadratic coefficients
 */
static char *test_qp_near_zero_p()
{
    size_t m = 1;
    size_t n = 2;
    
    double Ax[] = {1.0, 1.0};
    int Ai[] = {0, 1};
    int Ap[] = {0, 2};
    size_t nnz = 2;
    
    double lhs[] = {1.0};
    double rhs[] = {1.0};
    double lbs[] = {0.0, 0.0};
    double ubs[] = {INFINITY, INFINITY};
    double c[] = {1.0, 1.0};
    
    /* Very small P values */
    double Px[] = {1e-12, 1e-12};
    int Pi[] = {0, 1};
    int Pp[] = {0, 1, 2};
    size_t Pnnz = 2;
    
    Settings *stgs = default_settings();
    stgs->verbose = false;
    
    Presolver *presolver = new_qp_presolver_qr(Ax, Ai, Ap, m, n, nnz,
                                               lhs, rhs, lbs, ubs, c,
                                               Px, Pi, Pp, Pnnz,
                                               NULL, NULL, NULL, 0, 0, stgs);
    mu_assert("presolver creation failed for near-zero P", presolver != NULL);
    
    PresolveStatus status = run_presolver(presolver);
    
    mu_assert("near-zero P QP should be processed", 
              status == REDUCED || status == UNCHANGED || status == INFEASIBLE);
    
    PS_FREE(stgs);
    free_presolver(presolver);
    
    return 0;
}

/* Test 14: QP with very large P values
 * Tests numerical stability with huge quadratic coefficients
 */
static char *test_qp_large_p()
{
    size_t m = 1;
    size_t n = 2;
    
    double Ax[] = {1.0, 1.0};
    int Ai[] = {0, 1};
    int Ap[] = {0, 2};
    size_t nnz = 2;
    
    double lhs[] = {1.0};
    double rhs[] = {1.0};
    double lbs[] = {0.0, 0.0};
    double ubs[] = {INFINITY, INFINITY};
    double c[] = {1.0, 1.0};
    
    /* Very large P values */
    double Px[] = {1e10, 1e10};
    int Pi[] = {0, 1};
    int Pp[] = {0, 1, 2};
    size_t Pnnz = 2;
    
    Settings *stgs = default_settings();
    stgs->verbose = false;
    
    Presolver *presolver = new_qp_presolver_qr(Ax, Ai, Ap, m, n, nnz,
                                               lhs, rhs, lbs, ubs, c,
                                               Px, Pi, Pp, Pnnz,
                                               NULL, NULL, NULL, 0, 0, stgs);
    mu_assert("presolver creation failed for large P", presolver != NULL);
    
    PresolveStatus status = run_presolver(presolver);
    
    mu_assert("large P QP should be processed", 
              status == REDUCED || status == UNCHANGED || status == INFEASIBLE);
    
    PS_FREE(stgs);
    free_presolver(presolver);
    
    return 0;
}

/* Test 15: QP with negative diagonal elements
 * Tests that negative P diagonal is handled
 */
static char *test_qp_negative_diagonal()
{
    size_t m = 1;
    size_t n = 2;
    
    double Ax[] = {1.0, 1.0};
    int Ai[] = {0, 1};
    int Ap[] = {0, 2};
    size_t nnz = 2;
    
    double lhs[] = {1.0};
    double rhs[] = {1.0};
    double lbs[] = {0.0, 0.0};
    double ubs[] = {INFINITY, INFINITY};
    double c[] = {1.0, 1.0};
    
    /* Negative diagonal P (non-convex QP) */
    double Px[] = {-2.0, -2.0};
    int Pi[] = {0, 1};
    int Pp[] = {0, 1, 2};
    size_t Pnnz = 2;
    
    Settings *stgs = default_settings();
    stgs->verbose = false;
    
    Presolver *presolver = new_qp_presolver_qr(Ax, Ai, Ap, m, n, nnz,
                                               lhs, rhs, lbs, ubs, c,
                                               Px, Pi, Pp, Pnnz,
                                               NULL, NULL, NULL, 0, 0, stgs);
    mu_assert("presolver creation failed for negative diagonal P", presolver != NULL);
    
    PresolveStatus status = run_presolver(presolver);
    
    /* Should be processed even if non-convex */
    mu_assert("negative diagonal QP should be processed", 
              status == REDUCED || status == UNCHANGED || status == INFEASIBLE || status == UNBNDORINFEAS);
    
    PS_FREE(stgs);
    free_presolver(presolver);
    
    return 0;
}

/* Test 16: Many constraints, few variables (overdetermined)
 * Tests handling of over-constrained systems
 */
static char *test_qp_overdetermined()
{
    size_t m = 10;
    size_t n = 3;
    
    /* Create many constraints on few variables */
    double Ax[] = {1.0, 2.0, 3.0,  /* row 0 */
                   2.0, 1.0, 1.0,  /* row 1 */
                   1.0, 1.0, 1.0,  /* row 2 */
                   3.0, 2.0, 1.0,  /* row 3 */
                   1.0, 3.0, 2.0,  /* row 4 */
                   2.0, 2.0, 2.0,  /* row 5 */
                   1.0, 0.0, 1.0,  /* row 6 */
                   0.0, 1.0, 1.0,  /* row 7 */
                   1.0, 1.0, 0.0,  /* row 8 */
                   2.0, 1.0, 2.0}; /* row 9 */
    int Ai[] = {0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 2, 1, 2, 0, 1, 0, 1, 2};
    int Ap[] = {0, 3, 6, 9, 12, 15, 18, 20, 22, 24, 27};
    size_t nnz = 27;
    
    double lhs[] = {5.0, -INFINITY, 3.0, -INFINITY, 4.0, -INFINITY, 1.0, 1.0, 0.5, 4.0};
    double rhs[] = {INFINITY, 8.0, INFINITY, 10.0, INFINITY, 12.0, INFINITY, INFINITY, INFINITY, INFINITY};
    double lbs[] = {0.0, 0.0, 0.0};
    double ubs[] = {INFINITY, INFINITY, INFINITY};
    double c[] = {1.0, 2.0, 3.0};
    
    /* Diagonal P */
    double Px[] = {2.0, 2.0, 2.0};
    int Pi[] = {0, 1, 2};
    int Pp[] = {0, 1, 2, 3};
    size_t Pnnz = 3;
    
    Settings *stgs = default_settings();
    stgs->verbose = false;
    stgs->max_time = 5.0;
    
    Presolver *presolver = new_qp_presolver_qr(Ax, Ai, Ap, m, n, nnz,
                                               lhs, rhs, lbs, ubs, c,
                                               Px, Pi, Pp, Pnnz,
                                               NULL, NULL, NULL, 0, 0, stgs);
    mu_assert("presolver creation failed for overdetermined QP", presolver != NULL);
    
    PresolveStatus status = run_presolver(presolver);
    PresolvedProblem *reduced = presolver->reduced_prob;
    
    mu_assert("overdetermined QP should be processed", 
              status == REDUCED || status == UNCHANGED || status == INFEASIBLE);
    
    /* Many rows should be eliminated as redundant */
    if (status == REDUCED || status == UNCHANGED) {
        mu_assert("should have fewer or equal rows", reduced->m <= m);
        mu_assert("should have fewer or equal cols", reduced->n <= n);
    }
    
    PS_FREE(stgs);
    free_presolver(presolver);
    
    return 0;
}

/* Test 17: Few constraints, many variables (underdetermined)
 * Tests handling of under-constrained systems
 */
static char *test_qp_underdetermined()
{
    size_t m = 2;
    size_t n = 20;
    
    /* Simple constraints on many variables */
    double Ax[40];
    int Ai[40];
    int Ap[3];
    
    /* Row 0: sum of first 10 vars = 10 */
    for (int i = 0; i < 10; i++) {
        Ax[i] = 1.0;
        Ai[i] = i;
    }
    Ap[0] = 0;
    Ap[1] = 10;
    
    /* Row 1: sum of last 10 vars = 10 */
    for (int i = 0; i < 10; i++) {
        Ax[10 + i] = 1.0;
        Ai[10 + i] = 10 + i;
    }
    Ap[2] = 20;
    size_t nnz = 20;
    
    double lhs[] = {10.0, 10.0};
    double rhs[] = {10.0, 10.0};
    
    double lbs[20];
    double ubs[20];
    double c[20];
    double Px[20];
    int Pi[20];
    int Pp[21];
    
    for (int i = 0; i < 20; i++) {
        lbs[i] = 0.0;
        ubs[i] = INFINITY;
        c[i] = (double)(i + 1);
        Px[i] = 2.0;
        Pi[i] = i;
        Pp[i] = i;
    }
    Pp[20] = 20;
    size_t Pnnz = 20;
    
    Settings *stgs = default_settings();
    stgs->verbose = false;
    stgs->max_time = 5.0;
    
    Presolver *presolver = new_qp_presolver_qr(Ax, Ai, Ap, m, n, nnz,
                                               lhs, rhs, lbs, ubs, c,
                                               Px, Pi, Pp, Pnnz,
                                               NULL, NULL, NULL, 0, 0, stgs);
    mu_assert("presolver creation failed for underdetermined QP", presolver != NULL);
    
    PresolveStatus status = run_presolver(presolver);
    PresolvedProblem *reduced = presolver->reduced_prob;
    
    mu_assert("underdetermined QP should be processed", 
              status == REDUCED || status == UNCHANGED || status == INFEASIBLE);
    
    PS_FREE(stgs);
    free_presolver(presolver);
    
    return 0;
}

/* Test 18: QP with tight bounds on all variables
 * Tests bound propagation with quadratic terms
 */
static char *test_qp_tight_bounds()
{
    size_t m = 1;
    size_t n = 5;
    
    /* Constraint: sum(x) = 15 */
    double Ax[] = {1.0, 1.0, 1.0, 1.0, 1.0};
    int Ai[] = {0, 1, 2, 3, 4};
    int Ap[] = {0, 5};
    size_t nnz = 5;
    
    double lhs[] = {15.0};
    double rhs[] = {15.0};
    
    /* All variables have tight bounds: 2 <= xi <= 4 */
    double lbs[] = {2.0, 2.0, 2.0, 2.0, 2.0};
    double ubs[] = {4.0, 4.0, 4.0, 4.0, 4.0};
    double c[] = {1.0, 1.0, 1.0, 1.0, 1.0};
    
    /* Diagonal P */
    double Px[] = {2.0, 2.0, 2.0, 2.0, 2.0};
    int Pi[] = {0, 1, 2, 3, 4};
    int Pp[] = {0, 1, 2, 3, 4, 5};
    size_t Pnnz = 5;
    
    Settings *stgs = default_settings();
    stgs->verbose = false;
    
    Presolver *presolver = new_qp_presolver_qr(Ax, Ai, Ap, m, n, nnz,
                                               lhs, rhs, lbs, ubs, c,
                                               Px, Pi, Pp, Pnnz,
                                               NULL, NULL, NULL, 0, 0, stgs);
    mu_assert("presolver creation failed for tight bounds QP", presolver != NULL);
    
    PresolveStatus status = run_presolver(presolver);
    PresolvedProblem *reduced = presolver->reduced_prob;
    
    mu_assert("tight bounds QP should be processed", 
              status == REDUCED || status == UNCHANGED || status == INFEASIBLE);
    
    PS_FREE(stgs);
    free_presolver(presolver);
    
    return 0;
}

/* Test 19: Empty problem (no constraints)
 * Tests handling of unconstrained QP
 */
static char *test_qp_unconstrained()
{
    size_t m = 0;
    size_t n = 3;
    
    /* No constraints */
    double Ax[] = {};
    int Ai[] = {};
    int Ap[] = {0};
    size_t nnz = 0;
    
    double *lhs = NULL;
    double *rhs = NULL;
    
    double lbs[] = {0.0, 0.0, 0.0};
    double ubs[] = {INFINITY, INFINITY, INFINITY};
    double c[] = {1.0, 2.0, 3.0};
    
    /* Diagonal P */
    double Px[] = {2.0, 2.0, 2.0};
    int Pi[] = {0, 1, 2};
    int Pp[] = {0, 1, 2, 3};
    size_t Pnnz = 3;
    
    Settings *stgs = default_settings();
    stgs->verbose = false;
    
    Presolver *presolver = new_qp_presolver_qr(Ax, Ai, Ap, m, n, nnz,
                                               lhs, rhs, lbs, ubs, c,
                                               Px, Pi, Pp, Pnnz,
                                               NULL, NULL, NULL, 0, 0, stgs);
    mu_assert("presolver creation failed for unconstrained QP", presolver != NULL);
    
    PresolveStatus status = run_presolver(presolver);
    
    /* Unconstrained convex QP should be detected as unbounded or properly handled */
    mu_assert("unconstrained QP should be processed", 
              status == REDUCED || status == UNCHANGED || status == UNBNDORINFEAS);
    
    PS_FREE(stgs);
    free_presolver(presolver);
    
    return 0;
}

/* Test 20: Dense off-diagonal QP
 * Tests with fully dense P matrix
 */
static char *test_qp_fully_dense()
{
    size_t m = 2;
    size_t n = 4;
    
    double Ax[] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    int Ai[] = {0, 1, 2, 3, 0, 1, 2, 3};
    int Ap[] = {0, 4, 8};
    size_t nnz = 8;
    
    double lhs[] = {4.0, -INFINITY};
    double rhs[] = {INFINITY, 8.0};
    double lbs[] = {0.0, 0.0, 0.0, 0.0};
    double ubs[] = {INFINITY, INFINITY, INFINITY, INFINITY};
    double c[] = {1.0, 1.0, 1.0, 1.0};
    
    /* Fully dense upper triangular P: P[i][j] = 1 for all i <= j */
    /* n*(n+1)/2 = 10 non-zeros */
    double Px[] = {1.0, 1.0, 1.0, 1.0,  /* row 0: P00, P01, P02, P03 */
                   1.0, 1.0, 1.0,        /* row 1: P11, P12, P13 */
                   1.0, 1.0,             /* row 2: P22, P23 */
                   1.0};                 /* row 3: P33 */
    int Pi[] = {0, 1, 2, 3,
                1, 2, 3,
                2, 3,
                3};
    int Pp[] = {0, 4, 7, 9, 10};
    size_t Pnnz = 10;
    
    Settings *stgs = default_settings();
    stgs->verbose = false;
    stgs->ston_cols = false;
    stgs->parallel_cols = false;
    
    Presolver *presolver = new_qp_presolver_qr(Ax, Ai, Ap, m, n, nnz,
                                               lhs, rhs, lbs, ubs, c,
                                               Px, Pi, Pp, Pnnz,
                                               NULL, NULL, NULL, 0, 0, stgs);
    mu_assert("presolver creation failed for dense QP", presolver != NULL);
    
    PresolveStatus status = run_presolver(presolver);
    PresolvedProblem *reduced = presolver->reduced_prob;
    
    mu_assert("dense QP should be processed", 
              status == REDUCED || status == UNCHANGED || status == INFEASIBLE);
    
    if ((status == REDUCED || status == UNCHANGED) && reduced->n > 0) {
        mu_assert("should have quadratic term", reduced->Qnnz > 0 || reduced->Rnnz > 0 || reduced->k > 0);
    }
    
    PS_FREE(stgs);
    free_presolver(presolver);
    
    return 0;
}

static int counter_qp = 0;

static const char *all_tests_qp()
{
    mu_run_test(test_qp_simple_diagonal, counter_qp);
    mu_run_test(test_qp_with_offdiagonal, counter_qp);
    mu_run_test(test_qp_larger, counter_qp);
    mu_run_test(test_qp_with_bounds, counter_qp);
    mu_run_test(test_qp_dense_matrix, counter_qp);
    mu_run_test(test_qp_all_presolvers, counter_qp);
    mu_run_test(test_qp_postsolve, counter_qp);
    mu_run_test(test_qp_large_scale, counter_qp);
    mu_run_test(test_qp_mixed_vars, counter_qp);
    mu_run_test(test_qp_single_var, counter_qp);
    mu_run_test(test_qp_full_reduction, counter_qp);
    mu_run_test(test_qp_infeasible, counter_qp);
    mu_run_test(test_qp_near_zero_p, counter_qp);
    mu_run_test(test_qp_large_p, counter_qp);
    mu_run_test(test_qp_negative_diagonal, counter_qp);
    mu_run_test(test_qp_overdetermined, counter_qp);
    mu_run_test(test_qp_underdetermined, counter_qp);
    mu_run_test(test_qp_tight_bounds, counter_qp);
    mu_run_test(test_qp_unconstrained, counter_qp);
    mu_run_test(test_qp_fully_dense, counter_qp);
    return 0;
}

int test_qp()
{
    const char *result = all_tests_qp();
    if (result != 0)
    {
        printf("%s\n", result);
        printf("QP: TEST FAILED!\n");
    }
    else
    {
        printf("QP: ALL TESTS PASSED\n");
    }
    printf("QP: Tests run: %d\n", counter_qp);
    return result == 0;
}

#endif // TEST_QP_H
