/*
 * End-to-end tests for QR decomposition
 * Simulates full workflow: presolve -> solve reduced -> postsolve -> verify
 */

#ifndef TEST_QR_END2END_H
#define TEST_QR_END2END_H

#include "PSQP_API.h"
#include "PSQP_sol.h"
#include "test_macros.h"
#include <math.h>
#include <string.h>

#define E2E_TOL 1e-6

/* Test 1: Portfolio optimization style problem */
static char *test_qr_portfolio_optimization()
{
    printf("\n  Testing portfolio optimization...\n");
    
    size_t n = 5;
    size_t m = 1;
    
    double Ax[] = {1.0, 1.0, 1.0, 1.0, 1.0};
    int Ai[] = {0, 1, 2, 3, 4};
    int Ap[] = {0, 5};
    size_t nnz = 5;
    
    double lhs[] = {1.0};
    double rhs[] = {1.0};
    double lbs[] = {0.0, 0.0, 0.0, 0.0, 0.0};
    double ubs[] = {INFINITY, INFINITY, INFINITY, INFINITY, INFINITY};
    double c[] = {-0.1, -0.12, -0.08, -0.15, -0.11};
    
    double Dx[] = {0.05, 0.04, 0.06, 0.03, 0.05};
    int Di[] = {0, 1, 2, 3, 4};
    int Dp[] = {0, 1, 2, 3, 4, 5};
    size_t Dnnz = 5;
    
    double Fx[] = {0.3, 0.4, 0.2, 0.5, 0.35, 0.1, -0.1, 0.2, 0.0, 0.15};
    int Fi[] = {0, 0, 0, 0, 0, 1, 1, 1, 1, 1};
    int Fp[] = {0, 5, 10};
    size_t Fnnz = 10;
    size_t k = 2;
    
    Settings *stgs = default_settings();
    stgs->verbose = false;
    
    Presolver *presolver = new_qp_presolver_qr(Ax, Ai, Ap, m, n, nnz,
                                               lhs, rhs, lbs, ubs, c,
                                               Dx, Di, Dp, Dnnz,
                                               Fx, Fi, Fp, Fnnz, k, stgs);
    
    mu_assert("presolver creation failed", presolver != NULL);
    
    PresolveStatus status = run_presolver(presolver);
    PresolvedProblem *reduced = presolver->reduced_prob;
    
    mu_assert("presolve should complete", status == REDUCED || status == UNCHANGED);
    
    printf("    Original: %zu vars, %zu constraints\n", n, m);
    printf("    Reduced:  %zu vars, %zu constraints\n", reduced->n, reduced->m);
    
    double *x_reduced = (double *)calloc(reduced->n, sizeof(double));
    double *y_reduced = (double *)calloc(reduced->m, sizeof(double));
    double *z_reduced = (double *)calloc(reduced->n, sizeof(double));
    
    for (size_t i = 0; i < reduced->n; i++) {
        x_reduced[i] = 1.0 / reduced->n;
    }
    
    postsolve(presolver, x_reduced, y_reduced, z_reduced);
    
    double *x_orig = presolver->sol->x;
    double sum = 0.0;
    for (size_t i = 0; i < n; i++) {
        sum += x_orig[i];
    }
    mu_assert("budget constraint should be satisfied", fabs(sum - 1.0) < E2E_TOL);
    
    for (size_t i = 0; i < n; i++) {
        mu_assert("x should be non-negative", x_orig[i] >= -E2E_TOL);
    }
    
    printf("    Budget constraint: sum(x) = %.6f OK\n", sum);
    
    free(x_reduced);
    free(y_reduced);
    free(z_reduced);
    free_settings(stgs);
    free_presolver(presolver);
    
    return 0;
}

/* Test 2: Equality constraints */
static char *test_qr_equality_constraints()
{
    printf("\n  Testing equality constraints...\n");
    
    size_t n = 4;
    size_t m = 2;
    
    double Ax[] = {1.0, 1.0, 1.0, 1.0, 1.0, -1.0};
    int Ai[] = {0, 1, 2, 3, 0, 1};
    int Ap[] = {0, 4, 6};
    size_t nnz = 6;
    
    double lhs[] = {4.0, 0.0};
    double rhs[] = {4.0, 0.0};
    double lbs[] = {0.0, 0.0, 0.0, 0.0};
    double ubs[] = {INFINITY, INFINITY, INFINITY, INFINITY};
    double c[] = {1.0, 1.0, 1.0, 1.0};
    
    double Qx[] = {1.0, 1.0, 2.0, 2.0};
    int Qi[] = {0, 1, 2, 3};
    int Qp[] = {0, 1, 2, 3, 4};
    size_t Qnnz = 4;
    
    double Rx[] = {0.5, 0.5, 0.5, 0.5};
    int Ri[] = {0, 0, 0, 0};
    int Rp[] = {0, 1, 2, 3, 4};
    size_t Rnnz = 4;
    size_t k = 1;
    
    Settings *stgs = default_settings();
    stgs->verbose = false;
    
    Presolver *presolver = new_qp_presolver_qr(Ax, Ai, Ap, m, n, nnz,
                                               lhs, rhs, lbs, ubs, c,
                                               Qx, Qi, Qp, Qnnz,
                                               Rx, Ri, Rp, Rnnz, k, stgs);
    
    mu_assert("presolver creation failed", presolver != NULL);
    
    PresolveStatus status = run_presolver(presolver);
    PresolvedProblem *reduced = presolver->reduced_prob;
    
    printf("    Original: %zu vars, %zu constraints\n", n, m);
    printf("    Reduced:  %zu vars, %zu constraints\n", reduced->n, reduced->m);
    
    double *x_reduced = (double *)calloc(reduced->n, sizeof(double));
    double *y_reduced = (double *)calloc(reduced->m, sizeof(double));
    double *z_reduced = (double *)calloc(reduced->n, sizeof(double));
    
    for (size_t i = 0; i < reduced->n; i++) {
        x_reduced[i] = 1.0;
    }
    
    postsolve(presolver, x_reduced, y_reduced, z_reduced);
    
    double *x_orig = presolver->sol->x;
    double sum = x_orig[0] + x_orig[1] + x_orig[2] + x_orig[3];
    double diff = x_orig[0] - x_orig[1];
    
    mu_assert("sum constraint should be satisfied", fabs(sum - 4.0) < E2E_TOL);
    mu_assert("x0 = x1 should be satisfied", fabs(diff) < E2E_TOL);
    
    printf("    Sum constraint: sum(x) = %.6f OK\n", sum);
    printf("    Equality constraint: x0 - x1 = %.6f OK\n", diff);
    
    free(x_reduced);
    free(y_reduced);
    free(z_reduced);
    free_settings(stgs);
    free_presolver(presolver);
    
    return 0;
}

/* Test 3: Inequality constraints */
static char *test_qr_inequality_constraints()
{
    printf("\n  Testing inequality constraints...\n");
    
    size_t n = 3;
    size_t m = 2;
    
    double Ax[] = {1.0, 1.0, 1.0, 1.0};
    int Ai[] = {0, 1, 2, 0};
    int Ap[] = {0, 3, 4};
    size_t nnz = 4;
    
    double lhs[] = {-INFINITY, 1.0};
    double rhs[] = {5.0, INFINITY};
    double lbs[] = {0.0, 0.0, 0.0};
    double ubs[] = {INFINITY, INFINITY, INFINITY};
    double c[] = {1.0, 2.0, 3.0};
    
    double Qx[] = {2.0, 2.0, 2.0};
    int Qi[] = {0, 1, 2};
    int Qp[] = {0, 1, 2, 3};
    size_t Qnnz = 3;
    
    Settings *stgs = default_settings();
    stgs->verbose = false;
    stgs->ston_cols = false;
    stgs->dton_eq = false;
    
    Presolver *presolver = new_qp_presolver_qr(Ax, Ai, Ap, m, n, nnz,
                                               lhs, rhs, lbs, ubs, c,
                                               Qx, Qi, Qp, Qnnz,
                                               NULL, NULL, NULL, 0, 0,
                                               stgs);
    
    mu_assert("presolver creation failed", presolver != NULL);
    
    PresolveStatus status = run_presolver(presolver);
    PresolvedProblem *reduced = presolver->reduced_prob;
    
    printf("    Original: %zu vars, %zu constraints\n", n, m);
    printf("    Reduced:  %zu vars, %zu constraints\n", reduced->n, reduced->m);
    
    double *x_reduced = (double *)calloc(reduced->n, sizeof(double));
    double *y_reduced = (double *)calloc(reduced->m, sizeof(double));
    double *z_reduced = (double *)calloc(reduced->n, sizeof(double));
    
    for (size_t i = 0; i < reduced->n; i++) {
        x_reduced[i] = 1.0;
    }
    
    postsolve(presolver, x_reduced, y_reduced, z_reduced);
    
    double *x_orig = presolver->sol->x;
    
    mu_assert("x0 should be >= 1", x_orig[0] >= 1.0 - E2E_TOL);
    double sum = x_orig[0] + x_orig[1] + x_orig[2];
    mu_assert("sum should be <= 5", sum <= 5.0 + E2E_TOL);
    
    printf("    Lower bound: x0 = %.6f >= 1 OK\n", x_orig[0]);
    printf("    Inequality: sum(x) = %.6f <= 5 OK\n", sum);
    
    free(x_reduced);
    free(y_reduced);
    free(z_reduced);
    free_settings(stgs);
    free_presolver(presolver);
    
    return 0;
}

/* Test 4: QR output validity */
static char *test_qr_output_validity()
{
    printf("\n  Testing QR output validity...\n");
    
    size_t n = 4;
    size_t m = 1;
    
    double Ax[] = {1.0, 1.0, 1.0, 1.0};
    int Ai[] = {0, 1, 2, 3};
    int Ap[] = {0, 4};
    size_t nnz = 4;
    
    double lhs[] = {1.0};
    double rhs[] = {1.0};
    double lbs[] = {0.0, 0.0, 0.0, 0.0};
    double ubs[] = {INFINITY, INFINITY, INFINITY, INFINITY};
    double c[] = {1.0, 1.0, 1.0, 1.0};
    
    double Qx[] = {2.0, 0.5, 2.0, 0.3, 2.0, 2.0};
    int Qi[] = {0, 1, 1, 2, 2, 3};
    int Qp[] = {0, 1, 3, 5, 6};
    size_t Qnnz = 6;
    
    double Rx[] = {1.0, 0.5, 0.3, 0.8, 0.2, 0.4};
    int Ri[] = {0, 1, 0, 1, 0, 1};
    int Rp[] = {0, 2, 4, 5, 6};
    size_t Rnnz = 6;
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
    
    run_presolver(presolver);
    PresolvedProblem *reduced = presolver->reduced_prob;
    
    if (reduced->Qnnz > 0 || reduced->Rnnz > 0 || reduced->k > 0 && reduced->n > 0) {
        mu_assert("Qp[0] should be 0", reduced->Qp[0] == 0);
        mu_assert("Qp[n] should equal Qnnz", reduced->Qp[reduced->n] == (int)reduced->Qnnz);
        mu_assert("Rp[0] should be 0", reduced->Rp[0] == 0);
        mu_assert("Rp[n] should equal Rnnz", reduced->Rp[reduced->n] == (int)reduced->Rnnz);
        mu_assert("k should be preserved", reduced->k == k);
        
        /* Check Q is upper triangular */
        for (size_t i = 0; i < reduced->n; i++) {
            for (int idx = reduced->Qp[i]; idx < reduced->Qp[i+1]; idx++) {
                mu_assert("Q should be upper triangular", reduced->Qi[idx] >= (int)i);
            }
        }
        
        printf("    Q: %zu nnz, %zu rows OK\n", reduced->Qnnz, reduced->n);
        printf("    R: %zu nnz, %zu rows, k=%zu OK\n", reduced->Rnnz, reduced->n, reduced->k);
    }
    
    free_settings(stgs);
    free_presolver(presolver);
    
    return 0;
}

/* Test 5: Infeasible QP with QR
 * Tests that infeasible problems are correctly detected
 */
static char *test_qr_infeasible()
{
    printf("\n  Testing infeasible QP detection...\n");
    
    size_t n = 2;
    size_t m = 2;
    
    /* Infeasible constraints:
     * x0 + x1 >= 5
     * x0 + x1 <= 3
     * With x0, x1 >= 0
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
    
    /* QR format quadratic term */
    double Qx[] = {2.0, 2.0};
    int Qi[] = {0, 1};
    int Qp[] = {0, 1, 2};
    size_t Qnnz = 2;
    
    double Rx[] = {1.0, 1.0};
    int Ri[] = {0, 0};
    int Rp[] = {0, 1, 1};
    size_t Rnnz = 2;
    size_t k = 1;
    
    Settings *stgs = default_settings();
    stgs->verbose = false;
    
    Presolver *presolver = new_qp_presolver_qr(Ax, Ai, Ap, m, n, nnz,
                                               lhs, rhs, lbs, ubs, c,
                                               Qx, Qi, Qp, Qnnz,
                                               Rx, Ri, Rp, Rnnz, k, stgs);
    
    mu_assert("presolver creation failed", presolver != NULL);
    
    PresolveStatus status = run_presolver(presolver);
    
    printf("    Status: %s\n", 
           status == INFEASIBLE ? "INFEASIBLE" :
           status == REDUCED ? "REDUCED" :
           status == UNCHANGED ? "UNCHANGED" : "OTHER");
    
    /* Should detect infeasibility */
    mu_assert("infeasible QR QP should be detected", 
              status == INFEASIBLE || status == UNBNDORINFEAS);
    
    free_settings(stgs);
    free_presolver(presolver);
    
    return 0;
}

/* Test 6: Infeasible bounds with QR
 * Tests that conflicting bounds are detected
 */
static char *test_qr_infeasible_bounds()
{
    printf("\n  Testing infeasible bounds...\n");
    
    size_t n = 2;
    size_t m = 1;
    
    /* Constraint: x0 + x1 = 1 */
    double Ax[] = {1.0, 1.0};
    int Ai[] = {0, 1};
    int Ap[] = {0, 2};
    size_t nnz = 2;
    
    double lhs[] = {1.0};
    double rhs[] = {1.0};
    /* Infeasible bounds: x0 >= 2 and x0 <= 0.5 */
    double lbs[] = {2.0, 0.0};
    double ubs[] = {0.5, INFINITY};
    double c[] = {1.0, 1.0};
    
    /* QR format */
    double Qx[] = {1.0, 1.0};
    int Qi[] = {0, 1};
    int Qp[] = {0, 1, 2};
    size_t Qnnz = 2;
    
    double Rx[] = {0.5, 0.5};
    int Ri[] = {0, 0};
    int Rp[] = {0, 1, 1};
    size_t Rnnz = 2;
    size_t k = 1;
    
    Settings *stgs = default_settings();
    stgs->verbose = false;
    
    Presolver *presolver = new_qp_presolver_qr(Ax, Ai, Ap, m, n, nnz,
                                               lhs, rhs, lbs, ubs, c,
                                               Qx, Qi, Qp, Qnnz,
                                               Rx, Ri, Rp, Rnnz, k, stgs);
    
    mu_assert("presolver creation failed", presolver != NULL);
    
    PresolveStatus status = run_presolver(presolver);
    
    printf("    Status: %s\n", 
           status == INFEASIBLE ? "INFEASIBLE" :
           status == REDUCED ? "REDUCED" :
           status == UNCHANGED ? "UNCHANGED" : "OTHER");
    
    /* Should detect infeasibility from conflicting bounds */
    mu_assert("infeasible bounds should be detected", 
              status == INFEASIBLE || status == UNBNDORINFEAS);
    
    free_settings(stgs);
    free_presolver(presolver);
    
    return 0;
}

static int counter_qr_e2e = 0;

static const char *all_tests_qr_end2end()
{
    printf("\n=== QR End-to-End Tests ===");
    mu_run_test(test_qr_portfolio_optimization, counter_qr_e2e);
    mu_run_test(test_qr_equality_constraints, counter_qr_e2e);
    mu_run_test(test_qr_inequality_constraints, counter_qr_e2e);
    mu_run_test(test_qr_output_validity, counter_qr_e2e);
    mu_run_test(test_qr_infeasible, counter_qr_e2e);
    mu_run_test(test_qr_infeasible_bounds, counter_qr_e2e);
    return 0;
}

int test_qr_end2end()
{
    const char *result = all_tests_qr_end2end();
    if (result != 0)
    {
        printf("%s\n", result);
        printf("QR End-to-End: TEST FAILED!\n");
    }
    else
    {
        printf("\nQR End-to-End: ALL TESTS PASSED\n");
    }
    printf("QR End-to-End: Tests run: %d\n", counter_qr_e2e);
    return result == 0;
}

#endif // TEST_QR_END2END_H
