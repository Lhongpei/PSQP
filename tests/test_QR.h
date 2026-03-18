/*
 * Test cases for P = Q + R*R^T decomposition
 */

#ifndef TEST_QR_H
#define TEST_QR_H

#include "PSQP_API.h"
#include "test_macros.h"
#include <math.h>

/* Test 1: Simple QR decomposition
 * P = Q + R*R^T where Q is diagonal and R is low-rank
 */
static char *test_qr_simple()
{
    size_t m = 1;
    size_t n = 3;
    
    /* Constraint: sum(x) = 1 */
    double Ax[] = {1.0, 1.0, 1.0};
    int Ai[] = {0, 1, 2};
    int Ap[] = {0, 3};
    size_t nnz = 3;
    
    double lhs[] = {1.0};
    double rhs[] = {1.0};
    double lbs[] = {0.0, 0.0, 0.0};
    double ubs[] = {INFINITY, INFINITY, INFINITY};
    double c[] = {0.0, 0.0, 0.0};
    
    /* Q = I (diagonal) */
    double Qx[] = {1.0, 1.0, 1.0};
    int Qi[] = {0, 1, 2};
    int Qp[] = {0, 1, 2, 3};
    size_t Qnnz = 3;
    
    /* R = [1, 0; 0, 1; 0, 0] (3x2, rank 2)
     * R*R^T = [1, 0, 0; 0, 1, 0; 0, 0, 0]
     * P = Q + R*R^T = [2, 0, 0; 0, 2, 0; 0, 0, 1]
     */
    double Rx[] = {1.0, 1.0};
    int Ri[] = {0, 1};
    int Rp[] = {0, 1, 2, 2};
    size_t Rnnz = 2;
    size_t k = 2;
    
    Settings *stgs = default_settings();
    stgs->verbose = false;
    
    Presolver *presolver = new_qp_presolver_qr(Ax, Ai, Ap, m, n, nnz,
                                               lhs, rhs, lbs, ubs, c,
                                               Qx, Qi, Qp, Qnnz,
                                               Rx, Ri, Rp, Rnnz, k, stgs);
    
    mu_assert("presolver creation failed for QR", presolver != NULL);
    
    PresolveStatus status = run_presolver(presolver);
    PresolvedProblem *reduced = presolver->reduced_prob;
    
    mu_assert("QR presolve should work", status == REDUCED || status == UNCHANGED);
    
    if (reduced->n > 0) {
        mu_assert("should have quadratic term", reduced->has_quadratic);
    }
    
    free_settings(stgs);
    free_presolver(presolver);
    
    return 0;
}

/* Test 2: Only Q (no R) */
static char *test_qr_q_only()
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
    double c[] = {0.0, 0.0};
    
    /* Q = 2*I, R = NULL */
    double Qx[] = {2.0, 2.0};
    int Qi[] = {0, 1};
    int Qp[] = {0, 1, 2};
    size_t Qnnz = 2;
    
    Settings *stgs = default_settings();
    stgs->verbose = false;
    
    Presolver *presolver = new_qp_presolver_qr(Ax, Ai, Ap, m, n, nnz,
                                               lhs, rhs, lbs, ubs, c,
                                               Qx, Qi, Qp, Qnnz,
                                               NULL, NULL, NULL, 0, 0, stgs);
    
    mu_assert("presolver creation failed for Q-only", presolver != NULL);
    
    PresolveStatus status = run_presolver(presolver);
    mu_assert("Q-only presolve should work", status == REDUCED || status == UNCHANGED);
    
    free_settings(stgs);
    free_presolver(presolver);
    
    return 0;
}

/* Test 3: Only R (no Q) */
static char *test_qr_r_only()
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
    double c[] = {0.0, 0.0};
    
    /* Q = NULL, R = [1; 1] (2x1)
     * R*R^T = [1, 1; 1, 1]
     */
    double Rx[] = {1.0, 1.0};
    int Ri[] = {0, 0};
    int Rp[] = {0, 1, 1};
    size_t Rnnz = 2;
    size_t k = 1;
    
    Settings *stgs = default_settings();
    stgs->verbose = false;
    
    Presolver *presolver = new_qp_presolver_qr(Ax, Ai, Ap, m, n, nnz,
                                               lhs, rhs, lbs, ubs, c,
                                               NULL, NULL, NULL, 0,
                                               Rx, Ri, Rp, Rnnz, k, stgs);
    
    mu_assert("presolver creation failed for R-only", presolver != NULL);
    
    PresolveStatus status = run_presolver(presolver);
    mu_assert("R-only presolve should work", status == REDUCED || status == UNCHANGED);
    
    free_settings(stgs);
    free_presolver(presolver);
    
    return 0;
}

static int counter_qr = 0;

static const char *all_tests_qr()
{
    mu_run_test(test_qr_simple, counter_qr);
    mu_run_test(test_qr_q_only, counter_qr);
    mu_run_test(test_qr_r_only, counter_qr);
    return 0;
}

int test_qr()
{
    const char *result = all_tests_qr();
    if (result != 0)
    {
        printf("%s\n", result);
        printf("QR: TEST FAILED!\n");
    }
    else
    {
        printf("QR: ALL TESTS PASSED\n");
    }
    printf("QR: Tests run: %d\n", counter_qr);
    return result == 0;
}

#endif // TEST_QR_H
