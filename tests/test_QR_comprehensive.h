/*
 * Comprehensive tests for P = Q + R*R^T decomposition
 * Tests presolve and postsolve correctness
 */

#ifndef TEST_QR_COMPREHENSIVE_H
#define TEST_QR_COMPREHENSIVE_H

#include "PSQP_API.h"
#include "PSQP_sol.h"
#include "test_macros.h"
#include <math.h>
#include <string.h>

#define TOLERANCE 1e-6

/* Helper: Compute objective value for QP with P = Q + R*R^T */
static double compute_objective_qr(const double *x, const double *c, 
                                   const double *Qx, const int *Qi, const int *Qp,
                                   const double *Rx, const int *Ri, const int *Rp,
                                   size_t n, size_t k, double offset)
{
    double obj = offset;
    
    /* Linear term: c^T x */
    for (size_t i = 0; i < n; i++) {
        obj += c[i] * x[i];
    }
    
    /* Quadratic term: 0.5 * x^T Q x */
    if (Qx != NULL && Qp != NULL) {
        for (size_t i = 0; i < n; i++) {
            int row_start = Qp[i];
            int row_end = Qp[i + 1];
            for (int idx = row_start; idx < row_end; idx++) {
                int j = Qi[idx];
                double val = Qx[idx];
                if (i == j) {
                    obj += 0.5 * val * x[i] * x[i];
                } else {
                    obj += 0.5 * val * x[i] * x[j];
                }
            }
        }
    }
    
    /* Quadratic term: 0.5 * x^T R*R^T x = 0.5 * ||R^T x||^2 */
    if (Rx != NULL && Rp != NULL && k > 0) {
        /* Compute R^T x (k-dimensional vector) */
        double *RTx = (double *)calloc(k, sizeof(double));
        for (size_t i = 0; i < n; i++) {
            int row_start = Rp[i];
            int row_end = Rp[i + 1];
            for (int idx = row_start; idx < row_end; idx++) {
                int col = Ri[idx];
                RTx[col] += Rx[idx] * x[i];
            }
        }
        /* Add 0.5 * ||R^T x||^2 to objective */
        for (size_t j = 0; j < k; j++) {
            obj += 0.5 * RTx[j] * RTx[j];
        }
        free(RTx);
    }
    
    return obj;
}

/* Helper: Check if constraint is satisfied */
static int check_constraint(const double *x, const double *Ax, const int *Ai, 
                            const int *Ap, int row, double lhs, double rhs,
                            size_t n)
{
    double sum = 0.0;
    int row_start = Ap[row];
    int row_end = Ap[row + 1];
    
    for (int idx = row_start; idx < row_end; idx++) {
        int col = Ai[idx];
        sum += Ax[idx] * x[col];
    }
    
    if (!isinf(lhs) && sum < lhs - 1e-6) return 0;
    if (!isinf(rhs) && sum > rhs + 1e-6) return 0;
    return 1;
}

/* Test 1: End-to-end test with postsolve verification
 * Simple problem with known solution
 */
static char *test_qr_postsolve_simple()
{
    size_t m = 1;
    size_t n = 2;
    
    /* Constraint: x0 + x1 = 1 */
    double Ax[] = {1.0, 1.0};
    int Ai[] = {0, 1};
    int Ap[] = {0, 2};
    size_t nnz = 2;
    
    double lhs[] = {1.0};
    double rhs[] = {1.0};
    double lbs[] = {0.0, 0.0};
    double ubs[] = {INFINITY, INFINITY};
    double c[] = {0.0, 0.0};  /* No linear term */
    
    /* Q = 2*I, so quadratic term is 0.5 * (2*x0^2 + 2*x1^2) = x0^2 + x1^2 */
    double Qx[] = {2.0, 2.0};
    int Qi[] = {0, 1};
    int Qp[] = {0, 1, 2};
    size_t Qnnz = 2;
    
    /* R = [1; 0], so R*R^T adds [1, 0; 0, 0]
     * Total P = [3, 0; 0, 2]
     * Objective: 0.5 * (3*x0^2 + 2*x1^2) = 1.5*x0^2 + x1^2
     * With x0 + x1 = 1, optimal is x0 = 0.4, x1 = 0.6
     */
    double Rx[] = {1.0, 0.0};
    int Ri[] = {0, 0};
    int Rp[] = {0, 1, 1};
    size_t Rnnz = 2;
    size_t k = 1;
    
    Settings *stgs = default_settings();
    stgs->verbose = false;
    stgs->ston_cols = false;  /* Keep problem size for testing */
    stgs->dton_eq = false;
    
    Presolver *presolver = new_qp_presolver_qr(Ax, Ai, Ap, m, n, nnz,
                                               lhs, rhs, lbs, ubs, c,
                                               Qx, Qi, Qp, Qnnz,
                                               Rx, Ri, Rp, Rnnz, k, stgs);
    
    mu_assert("presolver creation failed", presolver != NULL);
    
    PresolveStatus status = run_presolver(presolver);
    PresolvedProblem *reduced = presolver->reduced_prob;
    
    mu_assert("presolve should work", status == REDUCED || status == UNCHANGED);
    
    /* Create a feasible solution for the reduced problem
     * Note: x_reduced must satisfy the reduced constraints */
    double *x_reduced = (double *)calloc(reduced->n, sizeof(double));
    double *y_reduced = (double *)calloc(reduced->m, sizeof(double));
    double *z_reduced = (double *)calloc(reduced->n, sizeof(double));
    
    /* Fill with values that sum to 1 (satisfy the constraint) */
    if (reduced->n > 0) {
        x_reduced[0] = 0.5;
        if (reduced->n > 1) {
            x_reduced[1] = 0.5;
        }
    }
    
    /* Run postsolve */
    postsolve(presolver, x_reduced, y_reduced, z_reduced);
    
    /* Get original solution */
    double *x_original = presolver->sol->x;
    
    /* Verify constraint is satisfied in original solution */
    mu_assert("constraint should be satisfied after postsolve", 
              check_constraint(x_original, Ax, Ai, Ap, 0, lhs[0], rhs[0], n));
    
    /* Verify bounds are satisfied */
    for (size_t i = 0; i < n; i++) {
        mu_assert("lower bound should be satisfied", x_original[i] >= lbs[i] - TOLERANCE);
        mu_assert("upper bound should be satisfied", x_original[i] <= ubs[i] + TOLERANCE);
    }
    
    free(x_reduced);
    free(y_reduced);
    free(z_reduced);
    free_settings(stgs);
    free_presolver(presolver);
    
    return 0;
}

/* Test 2: Variable fixing with postsolve
 * Tests that fixed variables are correctly recovered
 */
static char *test_qr_postsolve_fixed_vars()
{
    size_t m = 1;
    size_t n = 3;
    
    /* Constraint: x0 + x1 + x2 = 3 */
    double Ax[] = {1.0, 1.0, 1.0};
    int Ai[] = {0, 1, 2};
    int Ap[] = {0, 3};
    size_t nnz = 3;
    
    double lhs[] = {3.0};
    double rhs[] = {3.0};
    /* Fix x1 = 1.0 via tight bounds */
    double lbs[] = {0.0, 1.0, 0.0};
    double ubs[] = {INFINITY, 1.0, INFINITY};
    double c[] = {1.0, 2.0, 1.0};
    
    /* Q = diagonal */
    double Qx[] = {2.0, 3.0, 2.0};
    int Qi[] = {0, 1, 2};
    int Qp[] = {0, 1, 2, 3};
    size_t Qnnz = 3;
    
    /* R = rank-1 */
    double Rx[] = {1.0, 1.0, 1.0};
    int Ri[] = {0, 0, 0};
    int Rp[] = {0, 1, 2, 3};
    size_t Rnnz = 3;
    size_t k = 1;
    
    Settings *stgs = default_settings();
    stgs->verbose = false;
    stgs->ston_cols = false;
    stgs->dton_eq = false;
    
    Presolver *presolver = new_qp_presolver_qr(Ax, Ai, Ap, m, n, nnz,
                                               lhs, rhs, lbs, ubs, c,
                                               Qx, Qi, Qp, Qnnz,
                                               Rx, Ri, Rp, Rnnz, k, stgs);
    
    mu_assert("presolver creation failed", presolver != NULL);
    
    PresolveStatus status = run_presolver(presolver);
    mu_assert("presolve should work", status == REDUCED || status == UNCHANGED);
    
    /* Verify fixed variable is correctly identified */
    PresolvedProblem *reduced = presolver->reduced_prob;
    
    /* Create solution for reduced problem */
    double *x_reduced = (double *)calloc(reduced->n, sizeof(double));
    double *y_reduced = (double *)calloc(reduced->m, sizeof(double));
    double *z_reduced = (double *)calloc(reduced->n, sizeof(double));
    
    for (size_t i = 0; i < reduced->n; i++) {
        x_reduced[i] = 1.0;
    }
    
    /* Run postsolve */
    postsolve(presolver, x_reduced, y_reduced, z_reduced);
    
    /* Get original solution */
    double *x_original = presolver->sol->x;
    
    /* Verify x1 (fixed variable) is recovered correctly */
    mu_assert("fixed variable x1 should be 1.0", 
              fabs(x_original[1] - 1.0) < TOLERANCE);
    
    /* Verify constraint is satisfied */
    double sum = x_original[0] + x_original[1] + x_original[2];
    mu_assert("constraint sum should be 3", fabs(sum - 3.0) < TOLERANCE);
    
    free(x_reduced);
    free(y_reduced);
    free(z_reduced);
    free_settings(stgs);
    free_presolver(presolver);
    
    return 0;
}

/* Test 3: QR shrink verification
 * Tests that Q and R are correctly shrinked when variables are removed
 */
static char *test_qr_shrink_verification()
{
    size_t m = 1;
    size_t n = 4;
    
    /* Constraint: x0 + x1 + x2 + x3 >= 0 */
    double Ax[] = {1.0, 1.0, 1.0, 1.0};
    int Ai[] = {0, 1, 2, 3};
    int Ap[] = {0, 4};
    size_t nnz = 4;
    
    double lhs[] = {0.0};
    double rhs[] = {INFINITY};
    /* Fix x2 = 1.0 via tight bounds - this should remove x2 */
    double lbs[] = {0.0, 0.0, 1.0, 0.0};
    double ubs[] = {INFINITY, INFINITY, 1.0, INFINITY};
    double c[] = {1.0, 1.0, 1.0, 1.0};
    
    /* Q matrix with connections between all variables */
    /* Upper triangular: Q00, Q01, Q11, Q02, Q12, Q22, Q33 */
    double Qx[] = {2.0, 0.5, 2.0, 0.3, 0.4, 3.0, 2.0};
    int Qi[] = {0, 1, 1, 2, 2, 2, 3};
    int Qp[] = {0, 1, 3, 6, 7};
    size_t Qnnz = 7;
    
    /* R matrix: 4x2 */
    double Rx[] = {1.0, 0.5, 0.8, 0.3};  /* row 0, row 1, row 2, row 3 */
    int Ri[] = {0, 0, 1, 1};
    int Rp[] = {0, 1, 2, 3, 4};
    size_t Rnnz = 4;
    size_t k = 2;
    
    Settings *stgs = default_settings();
    stgs->verbose = false;
    stgs->ston_cols = false;
    stgs->dton_eq = false;
    
    Presolver *presolver = new_qp_presolver_qr(Ax, Ai, Ap, m, n, nnz,
                                               lhs, rhs, lbs, ubs, c,
                                               Qx, Qi, Qp, Qnnz,
                                               Rx, Ri, Rp, Rnnz, k, stgs);
    
    mu_assert("presolver creation failed", presolver != NULL);
    
    PresolveStatus status = run_presolver(presolver);
    PresolvedProblem *reduced = presolver->reduced_prob;
    
    mu_assert("presolve should work", status == REDUCED || status == UNCHANGED);
    
    /* Verify QR output if present */
    if (reduced->Qnnz > 0 || reduced->Rnnz > 0 || reduced->k > 0 && reduced->n > 0) {
        /* Q and R dimensions should match reduced problem size */
        mu_assert("Q should have n+1 row pointers", 
                  reduced->Qp != NULL && reduced->Qp[reduced->n] >= 0);
        mu_assert("R should have n+1 row pointers", 
                  reduced->Rp != NULL && reduced->Rp[reduced->n] >= 0);
        
        /* R should still have k columns */
        mu_assert("R should preserve k columns", reduced->k == k);
        
        /* Verify Q is upper triangular */
        for (size_t i = 0; i < reduced->n; i++) {
            int row_start = reduced->Qp[i];
            int row_end = reduced->Qp[i + 1];
            for (int idx = row_start; idx < row_end; idx++) {
                int col = reduced->Qi[idx];
                mu_assert("Q should be upper triangular (col >= row)", col >= (int)i);
            }
        }
    }
    
    free_settings(stgs);
    free_presolver(presolver);
    
    return 0;
}

/* Test 4: Empty Q (R only) with postsolve */
static char *test_qr_empty_q_postsolve()
{
    size_t m = 1;
    size_t n = 2;
    
    /* Constraint: x0 + x1 = 1 */
    double Ax[] = {1.0, 1.0};
    int Ai[] = {0, 1};
    int Ap[] = {0, 2};
    size_t nnz = 2;
    
    double lhs[] = {1.0};
    double rhs[] = {1.0};
    double lbs[] = {0.0, 0.0};
    double ubs[] = {INFINITY, INFINITY};
    double c[] = {0.0, 0.0};
    
    /* Q is empty (NULL) */
    /* R = [1; 1] (2x1), so P = R R^T = [[1, 1], [1, 1]] */
    double Rx[] = {1.0, 1.0};
    int Ri[] = {0, 0};
    int Rp[] = {0, 1, 2};
    size_t Rnnz = 2;
    size_t k = 1;
    
    Settings *stgs = default_settings();
    stgs->verbose = false;
    stgs->ston_cols = false;
    stgs->dton_eq = false;
    
    Presolver *presolver = new_qp_presolver_qr(Ax, Ai, Ap, m, n, nnz,
                                               lhs, rhs, lbs, ubs, c,
                                               NULL, NULL, NULL, 0,  /* No Q */
                                               Rx, Ri, Rp, Rnnz, k, stgs);
    
    mu_assert("presolver creation failed for empty Q", presolver != NULL);
    
    PresolveStatus status = run_presolver(presolver);
    mu_assert("presolve should work", status == REDUCED || status == UNCHANGED);
    
    /* Verify output has R but no Q */
    PresolvedProblem *reduced = presolver->reduced_prob;
    
    if (reduced->Qnnz > 0 || reduced->Rnnz > 0 || reduced->k > 0 && reduced->n > 0) {
        mu_assert("Qnnz should be 0", reduced->Qnnz == 0);
        mu_assert("Qx should be NULL", reduced->Qx == NULL);
        mu_assert("Rnnz should be > 0", reduced->Rnnz > 0);
        mu_assert("Rx should not be NULL", reduced->Rx != NULL);
    }
    
    /* Test postsolve */
    double *x_reduced = (double *)calloc(reduced->n, sizeof(double));
    double *y_reduced = (double *)calloc(reduced->m, sizeof(double));
    double *z_reduced = (double *)calloc(reduced->n, sizeof(double));
    
    for (size_t i = 0; i < reduced->n; i++) {
        x_reduced[i] = 0.5;
    }
    
    postsolve(presolver, x_reduced, y_reduced, z_reduced);
    
    /* Verify constraint */
    double *x_original = presolver->sol->x;
    double sum = x_original[0] + x_original[1];
    mu_assert("constraint should be satisfied", fabs(sum - 1.0) < TOLERANCE);
    
    free(x_reduced);
    free(y_reduced);
    free(z_reduced);
    free_settings(stgs);
    free_presolver(presolver);
    
    return 0;
}

/* Test 5: Empty R (Q only) with postsolve */
static char *test_qr_empty_r_postsolve()
{
    size_t m = 1;
    size_t n = 2;
    
    /* Constraint: x0 + x1 = 1 */
    double Ax[] = {1.0, 1.0};
    int Ai[] = {0, 1};
    int Ap[] = {0, 2};
    size_t nnz = 2;
    
    double lhs[] = {1.0};
    double rhs[] = {1.0};
    double lbs[] = {0.0, 0.0};
    double ubs[] = {INFINITY, INFINITY};
    double c[] = {0.0, 0.0};
    
    /* Q = 2*I, R is empty */
    double Qx[] = {2.0, 2.0};
    int Qi[] = {0, 1};
    int Qp[] = {0, 1, 2};
    size_t Qnnz = 2;
    
    Settings *stgs = default_settings();
    stgs->verbose = false;
    stgs->ston_cols = false;
    stgs->dton_eq = false;
    
    Presolver *presolver = new_qp_presolver_qr(Ax, Ai, Ap, m, n, nnz,
                                               lhs, rhs, lbs, ubs, c,
                                               Qx, Qi, Qp, Qnnz,
                                               NULL, NULL, NULL, 0, 0,  /* No R */
                                               stgs);
    
    mu_assert("presolver creation failed for empty R", presolver != NULL);
    
    PresolveStatus status = run_presolver(presolver);
    mu_assert("presolve should work", status == REDUCED || status == UNCHANGED);
    
    /* Verify output has Q but no R */
    PresolvedProblem *reduced = presolver->reduced_prob;
    
    if (reduced->Qnnz > 0 || reduced->Rnnz > 0 || reduced->k > 0 && reduced->n > 0) {
        mu_assert("Rnnz should be 0", reduced->Rnnz == 0);
        mu_assert("Rx should be NULL", reduced->Rx == NULL);
        mu_assert("Qnnz should be > 0", reduced->Qnnz > 0);
        mu_assert("Qx should not be NULL", reduced->Qx != NULL);
    }
    
    /* Test postsolve */
    double *x_reduced = (double *)calloc(reduced->n, sizeof(double));
    double *y_reduced = (double *)calloc(reduced->m, sizeof(double));
    double *z_reduced = (double *)calloc(reduced->n, sizeof(double));
    
    for (size_t i = 0; i < reduced->n; i++) {
        x_reduced[i] = 0.5;
    }
    
    postsolve(presolver, x_reduced, y_reduced, z_reduced);
    
    /* Verify constraint */
    double *x_original = presolver->sol->x;
    double sum = x_original[0] + x_original[1];
    mu_assert("constraint should be satisfied", fabs(sum - 1.0) < TOLERANCE);
    
    free(x_reduced);
    free(y_reduced);
    free(z_reduced);
    free_settings(stgs);
    free_presolver(presolver);
    
    return 0;
}

/* Test 6: Larger problem with multiple reductions */
static char *test_qr_large_problem()
{
    size_t n = 10;
    size_t m = 3;
    
    /* Create a larger problem with:
     * - sum(x) = 5
     * - x0 - x1 = 0
     * - x2 >= 1
     */
    double Ax[13] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,  /* row 0 */
                     1.0, -1.0,                                          /* row 1 */
                     1.0};                                               /* row 2 (x2 >= 1) */
    int Ai[13] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 1, 2};
    int Ap[4] = {0, 10, 12, 13};
    size_t nnz = 13;
    
    double lhs[3] = {5.0, 0.0, 1.0};
    double rhs[3] = {5.0, 0.0, INFINITY};
    double lbs[10] = {0};
    double ubs[10] = {INFINITY, INFINITY, INFINITY, INFINITY, INFINITY,
                      INFINITY, INFINITY, INFINITY, INFINITY, INFINITY};
    double c[10] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    
    /* Q: diagonal with some off-diagonal */
    double Qx[12] = {2, 0.5, 2, 0.3, 3, 0.4, 0.6, 4, 5, 6, 7, 8};
    int Qi[12] = {0, 1, 1, 2, 2, 3, 3, 3, 4, 5, 6, 7};
    int Qp[11] = {0, 1, 3, 5, 8, 9, 10, 11, 12, 12, 12};
    size_t Qnnz = 12;
    
    /* R: 10x3 low-rank */
    double Rx[10] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    int Ri[10] = {0, 1, 2, 0, 1, 2, 0, 1, 2, 0};
    int Rp[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    size_t Rnnz = 10;
    size_t k = 3;
    
    Settings *stgs = default_settings();
    stgs->verbose = false;
    
    Presolver *presolver = new_qp_presolver_qr(Ax, Ai, Ap, m, n, nnz,
                                               lhs, rhs, lbs, ubs, c,
                                               Qx, Qi, Qp, Qnnz,
                                               Rx, Ri, Rp, Rnnz, k, stgs);
    
    mu_assert("presolver creation failed for large problem", presolver != NULL);
    
    PresolveStatus status = run_presolver(presolver);
    mu_assert("large problem presolve should work", status == REDUCED || status == UNCHANGED);
    
    /* Test postsolve with random feasible solution */
    PresolvedProblem *reduced = presolver->reduced_prob;
    
    double *x_reduced = (double *)calloc(reduced->n, sizeof(double));
    double *y_reduced = (double *)calloc(reduced->m, sizeof(double));
    double *z_reduced = (double *)calloc(reduced->n, sizeof(double));
    
    /* Set feasible values within bounds */
    for (size_t i = 0; i < reduced->n; i++) {
        x_reduced[i] = reduced->lbs[i] + 0.5;
        if (isinf(x_reduced[i])) x_reduced[i] = 1.0;
    }
    
    postsolve(presolver, x_reduced, y_reduced, z_reduced);
    
    /* Verify bounds */
    double *x_original = presolver->sol->x;
    for (size_t i = 0; i < n; i++) {
        mu_assert("bound check after postsolve", x_original[i] >= lbs[i] - TOLERANCE);
    }
    
    free(x_reduced);
    free(y_reduced);
    free(z_reduced);
    free_settings(stgs);
    free_presolver(presolver);
    
    return 0;
}

static int counter_qr_comp = 0;

static const char *all_tests_qr_comprehensive()
{
    mu_run_test(test_qr_postsolve_simple, counter_qr_comp);
    mu_run_test(test_qr_postsolve_fixed_vars, counter_qr_comp);
    mu_run_test(test_qr_shrink_verification, counter_qr_comp);
    mu_run_test(test_qr_empty_q_postsolve, counter_qr_comp);
    mu_run_test(test_qr_empty_r_postsolve, counter_qr_comp);
    mu_run_test(test_qr_large_problem, counter_qr_comp);
    return 0;
}

int test_qr_comprehensive()
{
    const char *result = all_tests_qr_comprehensive();
    if (result != 0)
    {
        printf("%s\n", result);
        printf("QR Comprehensive: TEST FAILED!\n");
    }
    else
    {
        printf("QR Comprehensive: ALL TESTS PASSED\n");
    }
    printf("QR Comprehensive: Tests run: %d\n", counter_qr_comp);
    return result == 0;
}

#endif // TEST_QR_COMPREHENSIVE_H
