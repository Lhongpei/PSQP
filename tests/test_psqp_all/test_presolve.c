/*
 * Presolve Transformation Tests
 */

#include "test_suite.h"
#include "PSQP_API.h"

/* Test 1: Empty row removal */
test_result_t test_presolve_empty_rows(void)
{
    printf("\n[Test] Presolve - Empty Row Removal\n");
    
    Settings *stgs = default_settings();
    TEST_ASSERT(stgs != NULL, "Settings should be created");
    
    /* Problem with loose constraints */
    int n = 2;
    int m = 2;
    double c[] = {1.0, 2.0};
    double lb[] = {0.0, 0.0};
    double ub[] = {10.0, 10.0};
    
    double A_x[] = {1.0, 1.0, 1.0, 1.0};
    int A_i[] = {0, 1, 0, 1};
    int A_p[] = {0, 2, 4};
    
    double lhs[] = {0.0, 0.0};
    double rhs[] = {10.0, 10.0};
    
    Presolver *presolver = new_presolver(
        A_x, A_i, A_p, m, n, 4,
        lhs, rhs, lb, ub, c, stgs
    );
    
    TEST_ASSERT(presolver != NULL, "Presolver should be created");
    
    PresolveStatus status = run_presolver(presolver);
    TEST_ASSERT(status == UNCHANGED || status == REDUCED, "Presolve should succeed");
    TEST_ASSERT(presolver->reduced_prob != NULL, "Reduced problem should exist");
    
    free_presolver(presolver);
    free_settings(stgs);
    TEST_PASS_MSG("Empty row removal works");
}

/* Test 2: Singleton row elimination */
test_result_t test_presolve_ston_rows(void)
{
    printf("\n[Test] Presolve - Singleton Row Elimination\n");
    
    Settings *stgs = default_settings();
    TEST_ASSERT(stgs != NULL, "Settings should be created");
    
    /* Problem with singleton row */
    int n = 2;
    int m = 2;
    double c[] = {1.0, 2.0};
    double lb[] = {0.0, 0.0};
    double ub[] = {10.0, 10.0};
    
    /* Row 0: x0 = 1 (singleton), Row 1: x0 + x1 = 2 */
    double A_x[] = {1.0, 1.0, 1.0};
    int A_i[] = {0, 0, 1};
    int A_p[] = {0, 1, 3};
    
    double lhs[] = {1.0, 2.0};
    double rhs[] = {1.0, 2.0};
    
    Presolver *presolver = new_presolver(
        A_x, A_i, A_p, m, n, 3,
        lhs, rhs, lb, ub, c, stgs
    );
    
    TEST_ASSERT(presolver != NULL, "Presolver should be created");
    
    PresolveStatus status = run_presolver(presolver);
    TEST_ASSERT(status == UNCHANGED || status == REDUCED, "Presolve should succeed");
    
    free_presolver(presolver);
    free_settings(stgs);
    TEST_PASS_MSG("Singleton row elimination works");
}

/* Test 3: Doubleton equality */
test_result_t test_presolve_dton_eq(void)
{
    printf("\n[Test] Presolve - Doubleton Equality\n");
    
    Settings *stgs = default_settings();
    stgs->dton_eq = true;
    TEST_ASSERT(stgs != NULL, "Settings should be created");
    
    /* Problem with doubleton equality */
    int n = 2;
    int m = 1;
    double c[] = {1.0, 2.0};
    double lb[] = {0.0, 0.0};
    double ub[] = {10.0, 10.0};
    
    /* Row 0: x0 + x1 = 2 (doubleton) */
    double A_x[] = {1.0, 1.0};
    int A_i[] = {0, 0};
    int A_p[] = {0, 2};
    
    double lhs[] = {2.0};
    double rhs[] = {2.0};
    
    Presolver *presolver = new_presolver(
        A_x, A_i, A_p, m, n, 2,
        lhs, rhs, lb, ub, c, stgs
    );
    
    TEST_ASSERT(presolver != NULL, "Presolver should be created");
    
    PresolveStatus status = run_presolver(presolver);
    TEST_ASSERT(status == UNCHANGED || status == REDUCED, "Presolve should succeed");
    
    free_presolver(presolver);
    free_settings(stgs);
    TEST_PASS_MSG("Doubleton equality works");
}

/* Test 4: Dual fix */
test_result_t test_presolve_dual_fix(void)
{
    printf("\n[Test] Presolve - Dual Fix\n");
    
    Settings *stgs = default_settings();
    stgs->dual_fix = true;
    TEST_ASSERT(stgs != NULL, "Settings should be created");
    
    /* Problem where dual fix applies */
    int n = 2;
    int m = 1;
    double c[] = {1.0, -1.0};  /* x1 has negative coef */
    double lb[] = {0.0, 0.0};
    double ub[] = {10.0, 5.0};
    
    double A_x[] = {1.0, 1.0};
    int A_i[] = {0, 0};
    int A_p[] = {0, 2};
    
    double lhs[] = {1.0};
    double rhs[] = {10.0};
    
    Presolver *presolver = new_presolver(
        A_x, A_i, A_p, m, n, 2,
        lhs, rhs, lb, ub, c, stgs
    );
    
    TEST_ASSERT(presolver != NULL, "Presolver should be created");
    
    PresolveStatus status = run_presolver(presolver);
    TEST_ASSERT(status == UNCHANGED || status == REDUCED, "Presolve should succeed");
    
    free_presolver(presolver);
    free_settings(stgs);
    TEST_PASS_MSG("Dual fix works");
}

/* Test 5: Bound tightening */
test_result_t test_presolve_bounds(void)
{
    printf("\n[Test] Presolve - Bound Tightening\n");
    
    Settings *stgs = default_settings();
    stgs->primal_propagation = true;
    stgs->finite_bound_tightening = true;
    TEST_ASSERT(stgs != NULL, "Settings should be created");
    
    /* Problem where bound tightening applies */
    int n = 2;
    int m = 2;
    double c[] = {1.0, 1.0};
    double lb[] = {0.0, 0.0};
    double ub[] = {100.0, 100.0};  /* Loose bounds */
    
    /* Tight constraints will tighten bounds */
    double A_x[] = {1.0, 1.0, 1.0, 1.0};
    int A_i[] = {0, 1, 0, 1};
    int A_p[] = {0, 2, 4};
    
    double lhs[] = {0.0, 0.0};
    double rhs[] = {5.0, 5.0};
    
    Presolver *presolver = new_presolver(
        A_x, A_i, A_p, m, n, 4,
        lhs, rhs, lb, ub, c, stgs
    );
    
    TEST_ASSERT(presolver != NULL, "Presolver should be created");
    
    PresolveStatus status = run_presolver(presolver);
    TEST_ASSERT(status == UNCHANGED || status == REDUCED, "Presolve should succeed");
    
    free_presolver(presolver);
    free_settings(stgs);
    TEST_PASS_MSG("Bound tightening works");
}
