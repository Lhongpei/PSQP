/*
 * QP Step 5 - Access presolver->prob->obj
 */

#include <stdio.h>
#include <stdlib.h>
#include "test_suite.h"
#include "PSQP_API.h"
#include "Problem.h"

/* Test: Step 5 - Access presolver->prob->obj */
test_result_t test_qp_step5(void)
{
    printf("\n[Test] QP - Step 5: Access presolver->prob->obj\n");
    
    printf("  Step 5.1: Creating settings\n");
    Settings *stgs = default_settings();
    if (!stgs) {
        printf("    ❌ FAIL: Failed to create settings\n");
        g_stats.failed++; g_stats.total++;
        return TEST_FAIL;
    }
    stgs->verbose = false;
    
    printf("  Step 5.2: Preparing arrays\n");
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
    
    printf("  Step 5.3: Creating presolver\n");
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
    printf("  Step 5.4: Presolver created\n");
    
    printf("  Step 5.5: About to access presolver->prob\n");
    struct Problem *prob = presolver->prob;
    printf("  Step 5.6: presolver->prob = %p\n", (void*)prob);
    
    if (!prob) {
        printf("    ❌ FAIL: presolver->prob is NULL\n");
        free_presolver(presolver);
        free_settings(stgs);
        g_stats.failed++; g_stats.total++;
        return TEST_FAIL;
    }
    
    printf("  Step 5.7: About to access presolver->prob->obj\n");
    Objective *obj = presolver->prob->obj;
    printf("  Step 5.8: presolver->prob->obj = %p\n", (void*)obj);
    
    if (!obj) {
        printf("    ❌ FAIL: presolver->prob->obj is NULL\n");
        free_presolver(presolver);
        free_settings(stgs);
        g_stats.failed++; g_stats.total++;
        return TEST_FAIL;
    }
    
    printf("  Step 5.9: About to free presolver\n");
    free_presolver(presolver);
    printf("  Step 5.10: About to free settings\n");
    free_settings(stgs);
    
    printf("    ✅ PASS: Step 5 completed\n");
    g_stats.passed++; g_stats.total++;
    return TEST_PASS;
}
