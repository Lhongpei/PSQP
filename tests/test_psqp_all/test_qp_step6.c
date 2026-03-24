/*
 * QP Step 6 - Access presolver->prob->obj->quad_qr
 */

#include <stdio.h>
#include <stdlib.h>
#include "test_suite.h"
#include "PSQP_API.h"
#include "Problem.h"

/* Test: Step 6 - Access presolver->prob->obj->quad_qr */
test_result_t test_qp_step6(void)
{
    printf("\n[Test] QP - Step 6: Access quad_qr\n");
    
    printf("  Step 6.1: Creating settings\n");
    Settings *stgs = default_settings();
    if (!stgs) {
        printf("    ❌ FAIL: Failed to create settings\n");
        g_stats.failed++; g_stats.total++;
        return TEST_FAIL;
    }
    stgs->verbose = false;
    
    printf("  Step 6.2: Preparing arrays\n");
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
    
    printf("  Step 6.3: Creating presolver\n");
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
    printf("  Step 6.4: Presolver created\n");
    
    printf("  Step 6.5: Accessing presolver->prob\n");
    if (!presolver->prob) {
        printf("    ❌ FAIL: presolver->prob is NULL\n");
        free_presolver(presolver);
        free_settings(stgs);
        g_stats.failed++; g_stats.total++;
        return TEST_FAIL;
    }
    printf("  Step 6.6: presolver->prob OK\n");
    
    printf("  Step 6.7: Accessing presolver->prob->obj\n");
    if (!presolver->prob->obj) {
        printf("    ❌ FAIL: presolver->prob->obj is NULL\n");
        free_presolver(presolver);
        free_settings(stgs);
        g_stats.failed++; g_stats.total++;
        return TEST_FAIL;
    }
    printf("  Step 6.8: presolver->prob->obj OK\n");
    
    printf("  Step 6.9: Accessing presolver->prob->obj->quad_qr\n");
    QuadTermQR *qr = presolver->prob->obj->quad_qr;
    printf("  Step 6.10: quad_qr = %p\n", (void*)qr);
    
    if (!qr) {
        printf("    ❌ FAIL: quad_qr is NULL\n");
        free_presolver(presolver);
        free_settings(stgs);
        g_stats.failed++; g_stats.total++;
        return TEST_FAIL;
    }
    
    printf("  Step 6.11: Checking quad_qr fields\n");
    printf("    qr->n = %zu\n", qr->n);
    printf("    qr->Qnnz = %zu\n", qr->Qnnz);
    printf("    qr->has_quad = %d\n", qr->has_quad);
    
    printf("  Step 6.12: About to free presolver\n");
    free_presolver(presolver);
    printf("  Step 6.13: About to free settings\n");
    free_settings(stgs);
    
    printf("    ✅ PASS: Step 6 completed\n");
    g_stats.passed++; g_stats.total++;
    return TEST_PASS;
}
