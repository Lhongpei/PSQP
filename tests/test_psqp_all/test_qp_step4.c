/*
 * QP Step 4 - Access presolver->prob
 */

#include <stdio.h>
#include <stdlib.h>
#include "test_suite.h"
#include "PSQP_API.h"
#include "Problem.h"

/* Test: Step 4 - Access presolver->prob */
test_result_t test_qp_step4(void)
{
    printf("\n[Test] QP - Step 4: Access presolver->prob\n");
    
    printf("  Step 4.1: Creating settings\n");
    Settings *stgs = default_settings();
    if (!stgs) {
        printf("    ❌ FAIL: Failed to create settings\n");
        g_stats.failed++; g_stats.total++;
        return TEST_FAIL;
    }
    stgs->verbose = false;
    
    printf("  Step 4.2: Preparing arrays\n");
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
    
    printf("  Step 4.3: Creating presolver\n");
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
    printf("  Step 4.4: Presolver created at %p\n", (void*)presolver);
    
    printf("  Step 4.5: About to access presolver->prob\n");
    struct Problem *prob = presolver->prob;
    printf("  Step 4.6: presolver->prob = %p\n", (void*)prob);
    
    if (!prob) {
        printf("    ❌ FAIL: presolver->prob is NULL\n");
        free_presolver(presolver);
        free_settings(stgs);
        g_stats.failed++; g_stats.total++;
        return TEST_FAIL;
    }
    
    printf("  Step 4.7: About to free presolver\n");
    free_presolver(presolver);
    printf("  Step 4.8: About to free settings\n");
    free_settings(stgs);
    
    printf("    ✅ PASS: Step 4 completed\n");
    g_stats.passed++; g_stats.total++;
    return TEST_PASS;
}
