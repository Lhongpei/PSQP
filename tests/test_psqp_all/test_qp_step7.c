/*
 * QP Step 7 - Call has_only_q_diag
 */

#include <stdio.h>
#include <stdlib.h>
#include "test_suite.h"
#include "PSQP_API.h"
#include "Problem.h"

/* Forward declaration */
extern bool has_only_q_diag(const QuadTermQR *quad_qr, int k);

/* Test: Step 7 - Call has_only_q_diag */
test_result_t test_qp_step7(void)
{
    printf("\n[Test] QP - Step 7: Call has_only_q_diag\n");
    
    printf("  Step 7.1: Creating settings\n");
    Settings *stgs = default_settings();
    if (!stgs) {
        printf("    ❌ FAIL: Failed to create settings\n");
        g_stats.failed++; g_stats.total++;
        return TEST_FAIL;
    }
    stgs->verbose = false;
    
    printf("  Step 7.2: Preparing arrays\n");
    int n = 2, m = 1;
    double c[] = {1.0, 2.0};
    double lb[] = {0.0, 0.0};
    double ub[] = {10.0, 10.0};
    double Qx[] = {2.0, 3.0};
    int Qi[] = {0, 1};
    int Qp[] = {0, 1, 2};
    double A_x[] = {1.0, 1.0};
    int A_i[] = {0, 0};
    int A_p[] = {0, 2};
    double lhs[] = {3.0};
    double rhs[] = {3.0};
    
    printf("  Step 7.3: Creating presolver\n");
    Presolver *presolver = new_qp_presolver_qr(
        A_x, A_i, A_p, m, n, 2,
        lhs, rhs, lb, ub, c,
        Qx, Qi, Qp, 2,
        NULL, NULL, NULL, 0, 0,
        stgs
    );
    
    if (!presolver) {
        printf("    ❌ FAIL: Failed to create presolver\n");
        free_settings(stgs);
        g_stats.failed++; g_stats.total++;
        return TEST_FAIL;
    }
    printf("  Step 7.4: Presolver created\n");
    
    printf("  Step 7.5: Accessing quad_qr\n");
    QuadTermQR *qr = presolver->prob->obj->quad_qr;
    if (!qr) {
        printf("    ❌ FAIL: quad_qr is NULL\n");
        free_presolver(presolver);
        free_settings(stgs);
        g_stats.failed++; g_stats.total++;
        return TEST_FAIL;
    }
    printf("  Step 7.6: quad_qr OK\n");
    
    printf("  Step 7.7: Calling has_only_q_diag(qr, 0)\n");
    bool result0 = has_only_q_diag(qr, 0);
    printf("  Step 7.8: has_only_q_diag(qr, 0) = %d\n", result0);
    
    printf("  Step 7.9: Calling has_only_q_diag(qr, 1)\n");
    bool result1 = has_only_q_diag(qr, 1);
    printf("  Step 7.10: has_only_q_diag(qr, 1) = %d\n", result1);
    
    if (!result0 || !result1) {
        printf("    ❌ FAIL: Variables should be pure Q diagonal\n");
        free_presolver(presolver);
        free_settings(stgs);
        g_stats.failed++; g_stats.total++;
        return TEST_FAIL;
    }
    
    printf("  Step 7.11: About to free presolver\n");
    free_presolver(presolver);
    printf("  Step 7.12: About to free settings\n");
    free_settings(stgs);
    
    printf("    ✅ PASS: Step 7 completed\n");
    g_stats.passed++; g_stats.total++;
    return TEST_PASS;
}
