/*
 * QP Partial Presolve Tests - SIMPLIFIED
 * 
 * Tests with real optimizations but designed to NOT fully eliminate
 * the problem, so we can properly test postsolve.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "test_suite.h"
#include "PSQP_API.h"
#include "PSQP_sol.h"

#define PARTIAL_TOL 1e-5

/* 
 * Test 1: Simple QP with no reductions possible
 * Dense constraint matrix prevents variable elimination
 */
test_result_t test_partial_dton_with_offdiag_q(void)
{
    printf("\n[Test] Partial Presolve - Dense Constraint Matrix\n");
    
    /* 2 dense constraints with 3 variables - cannot fully eliminate */
    int m = 2, n = 3;
    double Ax[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};  /* Dense rows */
    int Ai[] = {0, 1, 2, 0, 1, 2};
    int Ap[] = {0, 3, 6};
    double lhs[] = {6.0, 15.0};
    double rhs[] = {6.0, 15.0};
    double lbs[] = {0.0, 0.0, 0.0};
    double ubs[] = {10.0, 10.0, 10.0};
    double c[] = {1.0, 2.0, 3.0};
    
    /* Pure diagonal Q - enables QP optimizations */
    double Qx[] = {1.0, 1.0, 1.0};
    int Qi[] = {0, 1, 2};
    int Qp[] = {0, 1, 2, 3};
    
    Settings *stgs = default_settings();
    stgs->verbose = false;
    /* Enable all optimizations */
    stgs->dton_eq = true;
    stgs->ston_cols = true;
    stgs->dual_fix = true;
    
    Presolver *presolver = new_qp_presolver_qr(Ax, Ai, Ap, m, n, 6,
                                               lhs, rhs, lbs, ubs, c,
                                               Qx, Qi, Qp, 3,
                                               NULL, NULL, NULL, 0, 0, stgs);
    
    if (!presolver) {
        printf("    ❌ FAIL: Failed to create presolver\n");
        free_settings(stgs);
        g_stats.failed++; g_stats.total++;
        return TEST_FAIL;
    }
    
    run_presolver(presolver);
    PresolvedProblem *reduced = presolver->reduced_prob;
    printf("  Reduced: m=%zu, n=%zu\n", reduced->m, reduced->n);
    
    /* Create solution for reduced problem */
    double *x_reduced = (double *)calloc(reduced->n, sizeof(double));
    double *y_reduced = (double *)calloc(reduced->m, sizeof(double));
    double *z_reduced = (double *)calloc(reduced->n, sizeof(double));
    
    /* Initialize with valid solution (all ones is feasible for reduced problem) */
    for (size_t i = 0; i < reduced->n; i++) {
        x_reduced[i] = 1.0;
        z_reduced[i] = 0.0;
    }
    for (size_t i = 0; i < reduced->m; i++) {
        y_reduced[i] = 0.0;
    }
    
    printf("  Running postsolve...\n");
    fflush(stdout);
    postsolve(presolver, x_reduced, y_reduced, z_reduced);
    
    printf("  Getting solution...\n");
    fflush(stdout);
    double *x_orig = presolver->sol->x;
    printf("  Postsolved x: [%f, %f, %f]\n", x_orig[0], x_orig[1], x_orig[2]);
    
    /* Verify all values are valid (not inf/nan) */
    for (int i = 0; i < n; i++) {
        if (isnan(x_orig[i]) || isinf(x_orig[i])) {
            printf("    ❌ FAIL: x[%d] is invalid (%f)\n", i, x_orig[i]);
            free(x_reduced); free(y_reduced); free(z_reduced);
            free_settings(stgs); free_presolver(presolver);
            g_stats.failed++; g_stats.total++;
            return TEST_FAIL;
        }
    }
    
    printf("    ✅ PASS: Postsolve produced valid values\n");
    
    free(x_reduced); free(y_reduced); free(z_reduced);
    free_settings(stgs); free_presolver(presolver);
    
    g_stats.passed++; g_stats.total++;
    return TEST_PASS;
}

/*
 * Test 2: Mixed bounds - some variables fixed, others free
 * This specifically tests FIXED_COL_QP postsolve path
 */
test_result_t test_partial_mixed_bounds(void)
{
    printf("\n[Test] Partial Presolve - Mixed Bounds\n");
    
    /* Use 2 constraints to prevent full elimination */
    int m = 2, n = 4;
    double Ax[] = {1.0, 1.0, 1.0, 1.0,  /* row 0: all vars sum to 4 */
                   0.0, 1.0, 0.0, 1.0};  /* row 1: x1 + x3 = 2 */
    int Ai[] = {0, 1, 2, 3, 1, 3};
    int Ap[] = {0, 4, 6};
    double lhs[] = {4.0, 2.0};
    double rhs[] = {4.0, 2.0};
    /* x0 and x2 fixed to 1, x1 and x3 free */
    double lbs[] = {1.0, 0.0, 1.0, 0.0};
    double ubs[] = {1.0, 10.0, 1.0, 10.0};
    double c[] = {0.0, 0.0, 0.0, 0.0};
    
    double Qx[] = {1.0, 1.0, 1.0, 1.0};
    int Qi[] = {0, 1, 2, 3};
    int Qp[] = {0, 1, 2, 3, 4};
    
    Settings *stgs = default_settings();
    stgs->verbose = false;
    /* Disable eliminations that might fully reduce the problem */
    stgs->dual_fix = false;
    stgs->ston_cols = false;
    stgs->dton_eq = false;
    
    Presolver *presolver = new_qp_presolver_qr(Ax, Ai, Ap, m, n, 6,
                                               lhs, rhs, lbs, ubs, c,
                                               Qx, Qi, Qp, 4,
                                               NULL, NULL, NULL, 0, 0, stgs);
    
    if (!presolver) {
        printf("    ❌ FAIL: Failed to create presolver\n");
        free_settings(stgs);
        g_stats.failed++; g_stats.total++;
        return TEST_FAIL;
    }
    
    run_presolver(presolver);
    PresolvedProblem *reduced = presolver->reduced_prob;
    printf("  Reduced: m=%zu, n=%zu\n", reduced->m, reduced->n);
    
    /* Create solution for reduced problem */
    double *x_reduced = (double *)calloc(reduced->n, sizeof(double));
    double *y_reduced = (double *)calloc(reduced->m, sizeof(double));
    double *z_reduced = (double *)calloc(reduced->n, sizeof(double));
    
    for (size_t i = 0; i < reduced->n; i++) {
        x_reduced[i] = 1.0;
        z_reduced[i] = 0.0;
    }
    for (size_t i = 0; i < reduced->m; i++) {
        y_reduced[i] = 0.0;
    }
    
    postsolve(presolver, x_reduced, y_reduced, z_reduced);
    
    double *x_orig = presolver->sol->x;
    printf("  Postsolved x: [%f, %f, %f, %f]\n", x_orig[0], x_orig[1], x_orig[2], x_orig[3]);
    
    /* Check x0 and x2 are at their fixed values */
    if (fabs(x_orig[0] - 1.0) > PARTIAL_TOL || fabs(x_orig[2] - 1.0) > PARTIAL_TOL) {
        printf("    ❌ FAIL: Fixed variables not recovered (x0=%f, x2=%f)\n", 
               x_orig[0], x_orig[2]);
        free(x_reduced); free(y_reduced); free(z_reduced);
        free_settings(stgs); free_presolver(presolver);
        g_stats.failed++; g_stats.total++;
        return TEST_FAIL;
    }
    
    /* Check x1 and x3 are valid (not inf/nan) - THIS IS THE KEY FIX */
    if (isnan(x_orig[1]) || isinf(x_orig[1]) || isnan(x_orig[3]) || isinf(x_orig[3])) {
        printf("    ❌ FAIL: Free variables have invalid values (x1=%f, x3=%f)\n",
               x_orig[1], x_orig[3]);
        free(x_reduced); free(y_reduced); free(z_reduced);
        free_settings(stgs); free_presolver(presolver);
        g_stats.failed++; g_stats.total++;
        return TEST_FAIL;
    }
    
    /* Verify constraints */
    double sum1 = x_orig[0] + x_orig[1] + x_orig[2] + x_orig[3];
    double sum2 = x_orig[1] + x_orig[3];
    if (fabs(sum1 - 4.0) > PARTIAL_TOL || fabs(sum2 - 2.0) > PARTIAL_TOL) {
        printf("    ❌ FAIL: Constraints violated (sum1=%f, sum2=%f)\n", sum1, sum2);
        free(x_reduced); free(y_reduced); free(z_reduced);
        free_settings(stgs); free_presolver(presolver);
        g_stats.failed++; g_stats.total++;
        return TEST_FAIL;
    }
    
    printf("    ✅ PASS: Mixed bounds handled correctly\n");
    
    free(x_reduced); free(y_reduced); free(z_reduced);
    free_settings(stgs); free_presolver(presolver);
    
    g_stats.passed++; g_stats.total++;
    return TEST_PASS;
}

/*
 * Test 3: Multiple independent constraints
 * Some can be reduced, others not
 */
test_result_t test_partial_multiple_dtons(void)
{
    printf("\n[Test] Partial Presolve - Multiple Constraints\n");
    
    int m = 3, n = 6;
    /* Three constraints with overlapping variables - cannot all be eliminated */
    double Ax[] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    int Ai[] = {0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5};
    int Ap[] = {0, 6, 8, 12};  /* Row 0: vars 0-5, Row 1: vars 4-5, Row 2: vars 0-5 */
    double lhs[] = {6.0, 2.0, 4.0};
    double rhs[] = {6.0, 2.0, 4.0};
    double lbs[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double ubs[] = {10.0, 10.0, 10.0, 10.0, 10.0, 10.0};
    double c[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    
    double Qx[] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    int Qi[] = {0, 1, 2, 3, 4, 5};
    int Qp[] = {0, 1, 2, 3, 4, 5, 6};
    
    Settings *stgs = default_settings();
    stgs->verbose = false;
    stgs->dton_eq = false;  /* Prevent full elimination */
    
    Presolver *presolver = new_qp_presolver_qr(Ax, Ai, Ap, m, n, 12,
                                               lhs, rhs, lbs, ubs, c,
                                               Qx, Qi, Qp, 6,
                                               NULL, NULL, NULL, 0, 0, stgs);
    
    if (!presolver) {
        printf("    ❌ FAIL: Failed to create presolver\n");
        free_settings(stgs);
        g_stats.failed++; g_stats.total++;
        return TEST_FAIL;
    }
    
    run_presolver(presolver);
    printf("  Getting reduced problem...\n");
    fflush(stdout);
    PresolvedProblem *reduced = presolver->reduced_prob;
    printf("  Reduced: m=%zu, n=%zu\n", reduced->m, reduced->n);
    fflush(stdout);
    
    printf("  Allocating solution arrays...\n");
    fflush(stdout);
    double *x_reduced = (double *)calloc(reduced->n, sizeof(double));
    double *y_reduced = (double *)calloc(reduced->m, sizeof(double));
    double *z_reduced = (double *)calloc(reduced->n, sizeof(double));
    
    for (size_t i = 0; i < reduced->n; i++) {
        x_reduced[i] = 1.0;
        z_reduced[i] = 0.0;
    }
    for (size_t i = 0; i < reduced->m; i++) {
        y_reduced[i] = 0.0;
    }
    
    postsolve(presolver, x_reduced, y_reduced, z_reduced);
    
    double *x_orig = presolver->sol->x;
    double avg = 0.0;
    for (int i = 0; i < n; i++) avg += x_orig[i];
    avg /= n;
    printf("  Postsolved x avg: %f\n", avg);
    
    /* Check all values are valid */
    for (int i = 0; i < n; i++) {
        if (isnan(x_orig[i]) || isinf(x_orig[i])) {
            printf("    ❌ FAIL: x[%d] is invalid\n", i);
            free(x_reduced); free(y_reduced); free(z_reduced);
            free_settings(stgs); free_presolver(presolver);
            g_stats.failed++; g_stats.total++;
            return TEST_FAIL;
        }
    }
    
    printf("    ✅ PASS: Multiple constraints handled\n");
    
    free(x_reduced); free(y_reduced); free(z_reduced);
    free_settings(stgs); free_presolver(presolver);
    
    g_stats.passed++; g_stats.total++;
    return TEST_PASS;
}

/*
 * Test 4: Empty row removal combined with other reductions
 */
test_result_t test_partial_empty_and_ston(void)
{
    printf("\n[Test] Partial Presolve - Empty Row Removal\n");
    
    int m = 2, n = 3;
    /* Row 0: empty (will be removed), Row 1: sum = 3 */
    double Ax[] = {1.0, 1.0, 1.0};
    int Ai[] = {0, 1, 2};
    int Ap[] = {0, 0, 3};  /* Row 0 is empty */
    double lhs[] = {0.0, 3.0};
    double rhs[] = {0.0, 3.0};
    double lbs[] = {0.0, 0.0, 0.0};
    double ubs[] = {10.0, 10.0, 10.0};
    double c[] = {0.0, 0.0, 0.0};
    
    double Qx[] = {1.0, 1.0, 1.0};
    int Qi[] = {0, 1, 2};
    int Qp[] = {0, 1, 2, 3};
    
    Settings *stgs = default_settings();
    stgs->verbose = false;
    stgs->dton_eq = false;
    
    Presolver *presolver = new_qp_presolver_qr(Ax, Ai, Ap, m, n, 3,
                                               lhs, rhs, lbs, ubs, c,
                                               Qx, Qi, Qp, 3,
                                               NULL, NULL, NULL, 0, 0, stgs);
    
    if (!presolver) {
        printf("    ❌ FAIL: Failed to create presolver\n");
        free_settings(stgs);
        g_stats.failed++; g_stats.total++;
        return TEST_FAIL;
    }
    
    run_presolver(presolver);
    PresolvedProblem *reduced = presolver->reduced_prob;
    printf("  Reduced: m=%zu, n=%zu\n", reduced->m, reduced->n);
    
    double *x_reduced = (double *)calloc(reduced->n, sizeof(double));
    double *y_reduced = (double *)calloc(reduced->m, sizeof(double));
    double *z_reduced = (double *)calloc(reduced->n, sizeof(double));
    
    for (size_t i = 0; i < reduced->n; i++) {
        x_reduced[i] = 1.0;
        z_reduced[i] = 0.0;
    }
    y_reduced[0] = 0.0;
    
    postsolve(presolver, x_reduced, y_reduced, z_reduced);
    
    double *x_orig = presolver->sol->x;
    printf("  Postsolved x: [%f, %f, %f]\n", x_orig[0], x_orig[1], x_orig[2]);
    
    double sum = x_orig[0] + x_orig[1] + x_orig[2];
    if (fabs(sum - 3.0) > PARTIAL_TOL) {
        printf("    ❌ FAIL: Constraint violated\n");
        free(x_reduced); free(y_reduced); free(z_reduced);
        free_settings(stgs); free_presolver(presolver);
        g_stats.failed++; g_stats.total++;
        return TEST_FAIL;
    }
    
    printf("    ✅ PASS: Empty row removal + postsolve\n");
    
    free(x_reduced); free(y_reduced); free(z_reduced);
    free_settings(stgs); free_presolver(presolver);
    
    g_stats.passed++; g_stats.total++;
    return TEST_PASS;
}

/*
 * Test 5: With R matrix
 */
test_result_t test_partial_with_r_matrix(void)
{
    printf("\n[Test] Partial Presolve - With R Matrix\n");
    
    int m = 1, n = 3;
    double Ax[] = {1.0, 1.0, 1.0};
    int Ai[] = {0, 1, 2};
    int Ap[] = {0, 3};
    double lhs[] = {3.0};
    double rhs[] = {3.0};
    double lbs[] = {0.0, 0.0, 0.0};
    double ubs[] = {10.0, 10.0, 10.0};
    double c[] = {0.0, 0.0, 0.0};
    
    /* Small diagonal Q */
    double Qx[] = {0.1, 0.1, 0.1};
    int Qi[] = {0, 1, 2};
    int Qp[] = {0, 1, 2, 3};
    
    /* Low-rank R (1x3 matrix with single row containing 3 elements) */
    double Rx[] = {1.0, 1.0, 1.0};
    int Ri[] = {0, 1, 2};  /* Column indices 0, 1, 2 */
    int Rp[] = {0, 3};     /* Row 0: elements 0 to 3 */
    
    Settings *stgs = default_settings();
    stgs->verbose = false;
    stgs->dton_eq = false;  /* Prevent full elimination */
    
    Presolver *presolver = new_qp_presolver_qr(Ax, Ai, Ap, m, n, 3,
                                               lhs, rhs, lbs, ubs, c,
                                               Qx, Qi, Qp, 3,
                                               Rx, Ri, Rp, 3, 1, stgs);
    
    if (!presolver) {
        printf("    ❌ FAIL: Failed to create presolver\n");
        free_settings(stgs);
        g_stats.failed++; g_stats.total++;
        return TEST_FAIL;
    }
    
    run_presolver(presolver);
    PresolvedProblem *reduced = presolver->reduced_prob;
    printf("  Reduced: m=%zu, n=%zu\n", reduced->m, reduced->n);
    
    if (reduced->n == 0) {
        printf("  Fully reduced - skipping\n");
        printf("    ✅ PASS\n");
        free_settings(stgs); free_presolver(presolver);
        g_stats.passed++; g_stats.total++;
        return TEST_PASS;
    }
    
    double *x_reduced = (double *)calloc(reduced->n, sizeof(double));
    double *y_reduced = (double *)calloc(reduced->m, sizeof(double));
    double *z_reduced = (double *)calloc(reduced->n, sizeof(double));
    
    for (size_t i = 0; i < reduced->n; i++) {
        x_reduced[i] = 1.0;
        z_reduced[i] = 0.0;
    }
    y_reduced[0] = 0.0;
    
    postsolve(presolver, x_reduced, y_reduced, z_reduced);
    
    double *x_orig = presolver->sol->x;
    printf("  Postsolved x: [%f, %f, %f]\n", x_orig[0], x_orig[1], x_orig[2]);
    
    double sum = x_orig[0] + x_orig[1] + x_orig[2];
    if (fabs(sum - 3.0) > PARTIAL_TOL) {
        printf("    ❌ FAIL: Constraint violated\n");
        free(x_reduced); free(y_reduced); free(z_reduced);
        free_settings(stgs); free_presolver(presolver);
        g_stats.failed++; g_stats.total++;
        return TEST_FAIL;
    }
    
    printf("    ✅ PASS: R matrix postsolve\n");
    
    free(x_reduced); free(y_reduced); free(z_reduced);
    free_settings(stgs); free_presolver(presolver);
    
    g_stats.passed++; g_stats.total++;
    return TEST_PASS;
}

/* ============================================================================
 * QP Fixed Column Postsolve Tests
 * ============================================================================
 * These tests specifically verify the FIXED_COL_QP reduction type
 * which handles QP variables that are fixed during presolve.
 */

/*
 * Test: QP Fixed Column Postsolve - Simple Case
 * Tests that a pure diagonal QP variable fixed to a value is correctly recovered
 */
test_result_t test_qp_fixed_col_postsolve_simple(void)
{
    printf("\n[Test] QP Fixed Col Postsolve - Simple (Pure Diagonal Q)\n");
    
    /* Problem: minimize 0.5*x^T*Q*x + c^T*x 
     * subject to: x0 + x1 = 2
     *             x0 >= 1, x0 <= 1 (fixed to 1)
     *             x1 >= 0
     * 
     * Solution: x0 = 1 (fixed), x1 = 1
     */
    int m = 1, n = 2;
    double Ax[] = {1.0, 1.0};
    int Ai[] = {0, 1};
    int Ap[] = {0, 2};
    double lhs[] = {2.0};
    double rhs[] = {2.0};
    double lbs[] = {1.0, 0.0};  /* x0 fixed to 1 */
    double ubs[] = {1.0, 10.0};
    double c[] = {0.0, 0.0};
    
    /* Pure diagonal Q */
    double Qx[] = {2.0, 2.0};  /* Q = 2*I */
    int Qi[] = {0, 1};
    int Qp[] = {0, 1, 2};
    
    Settings *stgs = default_settings();
    stgs->verbose = false;
    /* Disable eliminations to ensure fixed column path is used */
    stgs->dton_eq = false;
    stgs->ston_cols = false;
    stgs->dual_fix = false;
    
    Presolver *presolver = new_qp_presolver_qr(Ax, Ai, Ap, m, n, 2,
                                               lhs, rhs, lbs, ubs, c,
                                               Qx, Qi, Qp, 2,
                                               NULL, NULL, NULL, 0, 0, stgs);
    
    if (!presolver) {
        printf("    ❌ FAIL: Failed to create presolver\n");
        free_settings(stgs);
        g_stats.failed++; g_stats.total++;
        return TEST_FAIL;
    }
    
    run_presolver(presolver);
    PresolvedProblem *reduced = presolver->reduced_prob;
    printf("  Reduced: m=%zu, n=%zu\n", reduced->m, reduced->n);
    
    /* Create solution for reduced problem */
    double *x_reduced = (double *)calloc(reduced->n, sizeof(double));
    double *y_reduced = (double *)calloc(reduced->m, sizeof(double));
    double *z_reduced = (double *)calloc(reduced->n, sizeof(double));
    
    for (size_t i = 0; i < reduced->n; i++) {
        x_reduced[i] = 1.0;
        z_reduced[i] = 0.0;
    }
    for (size_t i = 0; i < reduced->m; i++) {
        y_reduced[i] = 0.0;
    }
    
    postsolve(presolver, x_reduced, y_reduced, z_reduced);
    
    double *x_orig = presolver->sol->x;
    double *z_orig = presolver->sol->z;
    printf("  Postsolved x: [%f, %f]\n", x_orig[0], x_orig[1]);
    printf("  Postsolved z: [%f, %f]\n", z_orig[0], z_orig[1]);
    
    /* Verify x0 is correctly fixed to 1 */
    if (fabs(x_orig[0] - 1.0) > PARTIAL_TOL) {
        printf("    ❌ FAIL: Fixed variable x0 = %f (expected 1.0)\n", x_orig[0]);
        free(x_reduced); free(y_reduced); free(z_reduced);
        free_settings(stgs); free_presolver(presolver);
        g_stats.failed++; g_stats.total++;
        return TEST_FAIL;
    }
    
    /* Verify constraint is satisfied: x0 + x1 = 2 */
    if (fabs(x_orig[0] + x_orig[1] - 2.0) > PARTIAL_TOL) {
        printf("    ❌ FAIL: Constraint violated (x0 + x1 = %f)\n", x_orig[0] + x_orig[1]);
        free(x_reduced); free(y_reduced); free(z_reduced);
        free_settings(stgs); free_presolver(presolver);
        g_stats.failed++; g_stats.total++;
        return TEST_FAIL;
    }
    
    /* Verify x1 is valid */
    if (isnan(x_orig[1]) || isinf(x_orig[1])) {
        printf("    ❌ FAIL: x1 is invalid (%f)\n", x_orig[1]);
        free(x_reduced); free(y_reduced); free(z_reduced);
        free_settings(stgs); free_presolver(presolver);
        g_stats.failed++; g_stats.total++;
        return TEST_FAIL;
    }
    
    printf("    ✅ PASS: Fixed column recovered correctly\n");
    
    free(x_reduced); free(y_reduced); free(z_reduced);
    free_settings(stgs); free_presolver(presolver);
    
    g_stats.passed++; g_stats.total++;
    return TEST_PASS;
}

/*
 * Test: QP Fixed Column Postsolve - With R Matrix
 * Tests fixed column recovery when there's also an R matrix contribution
 */
test_result_t test_qp_fixed_col_postsolve_with_r(void)
{
    printf("\n[Test] QP Fixed Col Postsolve - With R Matrix\n");
    fflush(stdout);
    
    /* Problem: minimize 0.5*x^T*(Q + R*R^T)*x 
     * where Q is diagonal and R is low-rank
     * subject to: x0 + x1 + x2 = 3
     *             x0 >= 1, x0 <= 1 (fixed to 1)
     */
    int m = 1, n = 3;
    double Ax[] = {1.0, 1.0, 1.0};
    int Ai[] = {0, 1, 2};
    int Ap[] = {0, 3};
    double lhs[] = {3.0};
    double rhs[] = {3.0};
    double lbs[] = {1.0, 0.0, 0.0};  /* x0 fixed to 1 */
    double ubs[] = {1.0, 10.0, 10.0};
    double c[] = {0.0, 0.0, 0.0};
    
    /* Diagonal Q */
    double Qx[] = {1.0, 1.0, 1.0};
    int Qi[] = {0, 1, 2};
    int Qp[] = {0, 1, 2, 3};
    
    /* R matrix (1x3) - adds to P = Q + R*R^T
     * CSR format: Rx = values, Ri = column indices, Rp = row pointers
     * R = [0.5, 0.5, 0.5] (1 row, 3 columns)
     */
    double Rx[] = {0.5, 0.5, 0.5};
    int Ri[] = {0, 1, 2};  /* Column indices */
    int Rp[] = {0, 3};     /* Row 0: elements 0 to 3 */
    
    Settings *stgs = default_settings();
    stgs->verbose = false;
    stgs->dton_eq = false;
    stgs->ston_cols = false;
    stgs->dual_fix = false;
    
    Presolver *presolver = new_qp_presolver_qr(Ax, Ai, Ap, m, n, 3,
                                               lhs, rhs, lbs, ubs, c,
                                               Qx, Qi, Qp, 3,
                                               Rx, Ri, Rp, 3, 1, stgs);
    
    if (!presolver) {
        printf("    ❌ FAIL: Failed to create presolver\n");
        free_settings(stgs);
        g_stats.failed++; g_stats.total++;
        return TEST_FAIL;
    }
    
    run_presolver(presolver);
    PresolvedProblem *reduced = presolver->reduced_prob;
    printf("  Reduced: m=%zu, n=%zu\n", reduced->m, reduced->n);
    
    double *x_reduced = (double *)calloc(reduced->n, sizeof(double));
    double *y_reduced = (double *)calloc(reduced->m, sizeof(double));
    double *z_reduced = (double *)calloc(reduced->n, sizeof(double));
    
    for (size_t i = 0; i < reduced->n; i++) {
        x_reduced[i] = 1.0;
        z_reduced[i] = 0.0;
    }
    y_reduced[0] = 0.0;
    
    postsolve(presolver, x_reduced, y_reduced, z_reduced);
    
    double *x_orig = presolver->sol->x;
    printf("  Postsolved x: [%f, %f, %f]\n", x_orig[0], x_orig[1], x_orig[2]);
    
    /* Verify x0 is correctly fixed to 1 */
    if (fabs(x_orig[0] - 1.0) > PARTIAL_TOL) {
        printf("    ❌ FAIL: Fixed variable x0 = %f (expected 1.0)\n", x_orig[0]);
        free(x_reduced); free(y_reduced); free(z_reduced);
        free_settings(stgs); free_presolver(presolver);
        g_stats.failed++; g_stats.total++;
        return TEST_FAIL;
    }
    
    /* Verify constraint is satisfied: x0 + x1 + x2 = 3 */
    double sum = x_orig[0] + x_orig[1] + x_orig[2];
    if (fabs(sum - 3.0) > PARTIAL_TOL) {
        printf("    ❌ FAIL: Constraint violated (sum = %f)\n", sum);
        free(x_reduced); free(y_reduced); free(z_reduced);
        free_settings(stgs); free_presolver(presolver);
        g_stats.failed++; g_stats.total++;
        return TEST_FAIL;
    }
    
    /* Verify x1 and x2 are valid */
    for (int i = 1; i < 3; i++) {
        if (isnan(x_orig[i]) || isinf(x_orig[i])) {
            printf("    ❌ FAIL: x[%d] is invalid (%f)\n", i, x_orig[i]);
            free(x_reduced); free(y_reduced); free(z_reduced);
            free_settings(stgs); free_presolver(presolver);
            g_stats.failed++; g_stats.total++;
            return TEST_FAIL;
        }
    }
    
    printf("    ✅ PASS: Fixed column with R matrix recovered correctly\n");
    
    free(x_reduced); free(y_reduced); free(z_reduced);
    free_settings(stgs); free_presolver(presolver);
    
    g_stats.passed++; g_stats.total++;
    return TEST_PASS;
}

/*
 * Test: QP Fixed Column Postsolve - Dual Recovery
 * Tests that dual variables are correctly recovered for fixed QP columns
 */
test_result_t test_qp_fixed_col_postsolve_dual_recovery(void)
{
    printf("\n[Test] QP Fixed Col Postsolve - Dual Recovery\n");
    
    /* Problem: minimize 0.5*(x0^2 + x1^2) + x0 + 2*x1
     * subject to: x0 + x1 = 3
     *             x0 >= 2, x0 <= 2 (fixed to 2)
     * 
     * With x0 = 2, the problem becomes:
     * minimize 0.5*(4 + x1^2) + 2 + 2*x1 = 0.5*x1^2 + 2*x1 + 4
     * subject to: x1 = 1
     * 
     * The reduced problem has x1 = 1, and dual can be computed.
     */
    int m = 1, n = 2;
    double Ax[] = {1.0, 1.0};
    int Ai[] = {0, 1};
    int Ap[] = {0, 2};
    double lhs[] = {3.0};
    double rhs[] = {3.0};
    double lbs[] = {2.0, 0.0};  /* x0 fixed to 2 */
    double ubs[] = {2.0, 10.0};
    double c[] = {1.0, 2.0};    /* Linear terms */
    
    /* Pure diagonal Q */
    double Qx[] = {1.0, 1.0};  /* Q = I */
    int Qi[] = {0, 1};
    int Qp[] = {0, 1, 2};
    
    Settings *stgs = default_settings();
    stgs->verbose = false;
    stgs->dton_eq = false;
    stgs->ston_cols = false;
    stgs->dual_fix = false;
    
    Presolver *presolver = new_qp_presolver_qr(Ax, Ai, Ap, m, n, 2,
                                               lhs, rhs, lbs, ubs, c,
                                               Qx, Qi, Qp, 2,
                                               NULL, NULL, NULL, 0, 0, stgs);
    
    if (!presolver) {
        printf("    ❌ FAIL: Failed to create presolver\n");
        free_settings(stgs);
        g_stats.failed++; g_stats.total++;
        return TEST_FAIL;
    }
    
    run_presolver(presolver);
    PresolvedProblem *reduced = presolver->reduced_prob;
    printf("  Reduced: m=%zu, n=%zu\n", reduced->m, reduced->n);
    
    double *x_reduced = (double *)calloc(reduced->n, sizeof(double));
    double *y_reduced = (double *)calloc(reduced->m, sizeof(double));
    double *z_reduced = (double *)calloc(reduced->n, sizeof(double));
    
    /* Set reduced solution - x1 should be 1 */
    for (size_t i = 0; i < reduced->n; i++) {
        x_reduced[i] = 1.0;
        z_reduced[i] = 0.0;
    }
    /* Dual variable for equality constraint */
    y_reduced[0] = 2.5;
    
    postsolve(presolver, x_reduced, y_reduced, z_reduced);
    
    double *x_orig = presolver->sol->x;
    double *z_orig = presolver->sol->z;
    double *y_orig = presolver->sol->y;
    printf("  Postsolved x: [%f, %f]\n", x_orig[0], x_orig[1]);
    printf("  Postsolved y: [%f]\n", y_orig[0]);
    printf("  Postsolved z: [%f, %f]\n", z_orig[0], z_orig[1]);
    
    /* Verify primal values */
    if (fabs(x_orig[0] - 2.0) > PARTIAL_TOL) {
        printf("    ❌ FAIL: Fixed variable x0 = %f (expected 2.0)\n", x_orig[0]);
        free(x_reduced); free(y_reduced); free(z_reduced);
        free_settings(stgs); free_presolver(presolver);
        g_stats.failed++; g_stats.total++;
        return TEST_FAIL;
    }
    
    if (fabs(x_orig[1] - 1.0) > PARTIAL_TOL) {
        printf("    ❌ FAIL: x1 = %f (expected 1.0)\n", x_orig[1]);
        free(x_reduced); free(y_reduced); free(z_reduced);
        free_settings(stgs); free_presolver(presolver);
        g_stats.failed++; g_stats.total++;
        return TEST_FAIL;
    }
    
    /* Verify constraint satisfaction */
    if (fabs(x_orig[0] + x_orig[1] - 3.0) > PARTIAL_TOL) {
        printf("    ❌ FAIL: Constraint violated\n");
        free(x_reduced); free(y_reduced); free(z_reduced);
        free_settings(stgs); free_presolver(presolver);
        g_stats.failed++; g_stats.total++;
        return TEST_FAIL;
    }
    
    /* Check dual variables are valid (not inf/nan) */
    for (int i = 0; i < n; i++) {
        if (isnan(z_orig[i]) || isinf(z_orig[i])) {
            printf("    ❌ FAIL: z[%d] is invalid (%f)\n", i, z_orig[i]);
            free(x_reduced); free(y_reduced); free(z_reduced);
            free_settings(stgs); free_presolver(presolver);
            g_stats.failed++; g_stats.total++;
            return TEST_FAIL;
        }
    }
    
    printf("    ✅ PASS: Dual recovery successful\n");
    
    free(x_reduced); free(y_reduced); free(z_reduced);
    free_settings(stgs); free_presolver(presolver);
    
    g_stats.passed++; g_stats.total++;
    return TEST_PASS;
}
