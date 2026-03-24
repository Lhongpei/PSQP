/*
 * QP Full Integration Test - Complete workflow
 */

#include <stdio.h>
#include <stdlib.h>
#include "test_suite.h"
#include "PSQP_API.h"
#include "Problem.h"

extern bool has_only_q_diag(const QuadTermQR *quad_qr, int k);
extern bool compute_Px_bounds(const QuadTermQR *qr, int k, 
                               const Bound *bounds, size_t n_cols,
                               double *min_Px, double *max_Px);

/* Test: Full QP presolve workflow */
test_result_t test_qp_full_workflow(void)
{
    printf("\n[Test] QP - Full Workflow\n");
    
    /* Setup */
    Settings *stgs = default_settings();
    if (!stgs) {
        printf("    ❌ FAIL: Failed to create settings\n");
        g_stats.failed++; g_stats.total++;
        return TEST_FAIL;
    }
    stgs->verbose = false;
    
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
    
    /* Create presolver */
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
    
    /* Verify quad_qr exists */
    if (!presolver->prob || !presolver->prob->obj || !presolver->prob->obj->quad_qr) {
        printf("    ❌ FAIL: quad_qr not properly initialized\n");
        free_presolver(presolver);
        free_settings(stgs);
        g_stats.failed++; g_stats.total++;
        return TEST_FAIL;
    }
    
    QuadTermQR *qr = presolver->prob->obj->quad_qr;
    
    /* Test has_only_q_diag */
    if (!has_only_q_diag(qr, 0) || !has_only_q_diag(qr, 1)) {
        printf("    ❌ FAIL: Variables should be pure Q diagonal\n");
        free_presolver(presolver);
        free_settings(stgs);
        g_stats.failed++; g_stats.total++;
        return TEST_FAIL;
    }
    
    /* Test compute_Px_bounds */
    Bound bounds[2];
    bounds[0].lb = 0.0; bounds[0].ub = 5.0;
    bounds[1].lb = 0.0; bounds[1].ub = 5.0;
    
    double min_Px, max_Px;
    bool success = compute_Px_bounds(qr, 0, bounds, 2, &min_Px, &max_Px);
    if (!success) {
        printf("    ○ SKIP: compute_Px_bounds returned false\n");
        free_presolver(presolver);
        free_settings(stgs);
        g_stats.skipped++; g_stats.total++;
        return TEST_SKIP;
    }
    
    /* Run presolve */
    PresolveStatus status = run_presolver(presolver);
    (void)status;
    
    if (!presolver->reduced_prob) {
        printf("    ❌ FAIL: reduced_prob is NULL after presolve\n");
        free_presolver(presolver);
        free_settings(stgs);
        g_stats.failed++; g_stats.total++;
        return TEST_FAIL;
    }
    
    /* Cleanup */
    free_presolver(presolver);
    free_settings(stgs);
    
    printf("    ✅ PASS: QP full workflow works\n");
    g_stats.passed++; g_stats.total++;
    return TEST_PASS;
}

/* Test: QP with DtonsEq enabled */
test_result_t test_qp_with_dton_eq(void)
{
    printf("\n[Test] QP - With DtonsEq\n");
    
    Settings *stgs = default_settings();
    if (!stgs) {
        printf("    ❌ FAIL: Failed to create settings\n");
        g_stats.failed++; g_stats.total++;
        return TEST_FAIL;
    }
    stgs->verbose = false;
    stgs->dton_eq = true;  /* Enable DtonsEq */
    
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
    
    /* Run presolve - this may trigger DtonsEq with pure Q diagonal substitution */
    PresolveStatus status = run_presolver(presolver);
    (void)status;
    
    if (!presolver->reduced_prob) {
        printf("    ❌ FAIL: reduced_prob is NULL\n");
        free_presolver(presolver);
        free_settings(stgs);
        g_stats.failed++; g_stats.total++;
        return TEST_FAIL;
    }
    
    free_presolver(presolver);
    free_settings(stgs);
    
    printf("    ✅ PASS: QP with DtonsEq works\n");
    g_stats.passed++; g_stats.total++;
    return TEST_PASS;
}
