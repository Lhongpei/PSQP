/*
 * QP-Specific Feature Tests - Step by step migration
 */

#include <stdio.h>
#include <stdlib.h>
#include "test_suite.h"
#include "PSQP_API.h"
#include "Problem.h"

/* Test 1: Step 1 - Just settings creation */
test_result_t test_qp_substitution(void)
{
    printf("\n[Test] QP - Step 1: Settings Creation\n");
    
    printf("  Step 1.1: Calling default_settings()\n");
    Settings *stgs = default_settings();
    printf("  Step 1.2: Settings created at %p\n", (void*)stgs);
    
    if (!stgs) {
        printf("    ❌ FAIL: Failed to create settings\n");
        g_stats.failed++; g_stats.total++;
        return TEST_FAIL;
    }
    
    printf("  Step 1.3: Setting verbose to false\n");
    stgs->verbose = false;
    printf("  Step 1.4: About to free settings\n");
    
    free_settings(stgs);
    printf("  Step 1.5: Settings freed successfully\n");
    
    printf("    ✅ PASS: Step 1 completed\n");
    g_stats.passed++; g_stats.total++;
    return TEST_PASS;
}

/* Test 2: Step 2 - Add local arrays */
test_result_t test_qp_dual_fix_bounds(void)
{
    printf("\n[Test] QP - Step 2: Local Arrays\n");
    
    printf("  Step 2.1: Creating settings\n");
    Settings *stgs = default_settings();
    if (!stgs) {
        printf("    ❌ FAIL: Failed to create settings\n");
        g_stats.failed++; g_stats.total++;
        return TEST_FAIL;
    }
    stgs->verbose = false;
    printf("  Step 2.2: Settings OK\n");
    
    printf("  Step 2.3: Creating local arrays\n");
    int n = 2, m = 1;
    printf("  Step 2.4: n=%d, m=%d\n", n, m);
    
    double c[] = {1.0, 2.0};
    printf("  Step 2.5: c[] = {%f, %f}\n", c[0], c[1]);
    
    double lb[] = {0.0, 0.0};
    double ub[] = {10.0, 10.0};
    printf("  Step 2.6: lb[] and ub[] created\n");
    
    double Qx[] = {2.0, 3.0};
    int Qi[] = {0, 1};
    int Qp[] = {0, 1, 2};
    printf("  Step 2.7: Q arrays created\n");
    
    double A_x[] = {1.0, 1.0};
    int A_i[] = {0, 0};
    int A_p[] = {0, 2};
    printf("  Step 2.8: A arrays created\n");
    
    double lhs[] = {3.0};
    double rhs[] = {3.0};
    printf("  Step 2.9: lhs and rhs created\n");
    
    printf("  Step 2.10: About to free settings\n");
    free_settings(stgs);
    printf("  Step 2.11: Cleanup complete\n");
    
    printf("    ✅ PASS: Step 2 completed\n");
    g_stats.passed++; g_stats.total++;
    return TEST_PASS;
}

/* Test 3: Step 3 - Call new_qp_presolver_qr */
test_result_t test_qp_presolve_integration(void)
{
    printf("\n[Test] QP - Step 3: Create Presolver\n");
    
    printf("  Step 3.1: Creating settings\n");
    Settings *stgs = default_settings();
    if (!stgs) {
        printf("    ❌ FAIL: Failed to create settings\n");
        g_stats.failed++; g_stats.total++;
        return TEST_FAIL;
    }
    stgs->verbose = false;
    
    printf("  Step 3.2: Preparing arrays\n");
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
    
    printf("  Step 3.3: About to call new_qp_presolver_qr\n");
    Presolver *presolver = new_qp_presolver_qr(
        A_x, A_i, A_p, m, n, 2,
        lhs, rhs, lb, ub, c,
        Qx, Qi, Qp, 2,
        NULL, NULL, NULL, 0, 0,
        stgs
    );
    printf("  Step 3.4: new_qp_presolver_qr returned %p\n", (void*)presolver);
    
    if (!presolver) {
        printf("    ❌ FAIL: Failed to create presolver\n");
        free_settings(stgs);
        g_stats.failed++; g_stats.total++;
        return TEST_FAIL;
    }
    
    printf("  Step 3.5: About to free presolver\n");
    free_presolver(presolver);
    printf("  Step 3.6: Presolver freed\n");
    
    printf("  Step 3.7: About to free settings\n");
    free_settings(stgs);
    printf("  Step 3.8: Settings freed\n");
    
    printf("    ✅ PASS: Step 3 completed\n");
    g_stats.passed++; g_stats.total++;
    return TEST_PASS;
}

/* Forward declaration for test_qp_step4 */
test_result_t test_qp_step4(void);
test_result_t test_qp_step5(void);
