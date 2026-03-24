/*
 * QP Step 8 - Call run_presolver
 */

#include <stdio.h>
#include <stdlib.h>
#include "test_suite.h"
#include "PSQP_API.h"
#include "Problem.h"

/* Test: Step 8 - Call run_presolver */
test_result_t test_qp_step8(void)
{
    printf("\n[Test] QP - Step 8: Call run_presolver\n");
    
    printf("  Step 8.1: Creating settings\n");
    Settings *stgs = default_settings();
    if (!stgs) {
        printf("    ❌ FAIL: Failed to create settings\n");
        g_stats.failed++; g_stats.total++;
        return TEST_FAIL;
    }
    stgs->verbose = false;
    stgs->dton_eq = false;  /* Disable to avoid complex substitution */
    
    printf("  Step 8.2: Preparing arrays\n");
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
    
    printf("  Step 8.3: Creating presolver\n");
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
    printf("  Step 8.4: Presolver created\n");
    
    printf("  Step 8.5: Calling run_presolver\n");
    PresolveStatus status = run_presolver(presolver);
    printf("  Step 8.6: run_presolver returned status = %d\n", (int)status);
    
    printf("  Step 8.7: Checking reduced_prob\n");
    if (!presolver->reduced_prob) {
        printf("    ❌ FAIL: reduced_prob is NULL\n");
        free_presolver(presolver);
        free_settings(stgs);
        g_stats.failed++; g_stats.total++;
        return TEST_FAIL;
    }
    printf("  Step 8.8: reduced_prob OK\n");
    
    printf("  Step 8.9: About to free presolver\n");
    free_presolver(presolver);
    printf("  Step 8.10: About to free settings\n");
    free_settings(stgs);
    
    printf("    ✅ PASS: Step 8 completed\n");
    g_stats.passed++; g_stats.total++;
    return TEST_PASS;
}
