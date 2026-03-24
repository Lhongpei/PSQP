/*
 * Problem Construction Tests
 */

#include "test_suite.h"
#include "PSQP_API.h"

/* Test 1: Create LP problem */
test_result_t test_problem_create_lp(void)
{
    printf("\n[Test] Problem Creation - LP\n");
    
    /* Simple LP: min c^T x s.t. Ax = b, lb <= x <= ub */
    int n = 3;
    int m = 2;
    
    double c[] = {1.0, 2.0, 3.0};
    double lb[] = {0.0, 0.0, 0.0};
    double ub[] = {10.0, 10.0, 10.0};
    
    /* A = [1 1 0; 0 1 1] */
    double A_x[] = {1.0, 1.0, 1.0, 1.0};
    int A_i[] = {0, 1, 0, 1};
    int A_p[] = {0, 2, 4};
    
    double lhs[] = {1.0, 2.0};
    double rhs[] = {1.0, 2.0};
    
    Settings *stgs = default_settings();
    TEST_ASSERT(stgs != NULL, "Settings should be created");
    
    Presolver *presolver = new_presolver(
        A_x, A_i, A_p, m, n, 4,
        lhs, rhs, lb, ub, c, stgs
    );
    
    TEST_ASSERT(presolver != NULL, "LP presolver should be created");
    TEST_ASSERT(presolver->prob != NULL, "Problem should exist");
    
    free_presolver(presolver);
    free_settings(stgs);
    TEST_PASS_MSG("LP problem creation works");
}

/* Test 2: Create QP problem */
test_result_t test_problem_create_qp(void)
{
    printf("\n[Test] Problem Creation - QP\n");
    
    int n = 2;
    int m = 1;
    
    double c[] = {1.0, 2.0};
    double lb[] = {0.0, 0.0};
    double ub[] = {10.0, 10.0};
    
    /* P = Q + R*R^T: Q = diag(2, 3), R = [1, 1] */
    double Qx[] = {2.0, 3.0};
    int Qi[] = {0, 1};
    int Qp[] = {0, 1, 2};
    double Rx[] = {1.0, 1.0};
    int Ri[] = {0, 1};
    int Rp[] = {0, 2};
    
    double A_x[] = {1.0, 1.0};
    int A_i[] = {0, 0};
    int A_p[] = {0, 2};
    
    double lhs[] = {1.0};
    double rhs[] = {1.0};
    
    Settings *stgs = default_settings();
    TEST_ASSERT(stgs != NULL, "Settings should be created");
    
    Presolver *presolver = new_qp_presolver_qr(
        A_x, A_i, A_p, m, n, 2,
        lhs, rhs, lb, ub, c,
        Qx, Qi, Qp, 2,
        Rx, Ri, Rp, 2, 1,
        stgs
    );
    
    TEST_ASSERT(presolver != NULL, "QP presolver should be created");
    TEST_ASSERT(presolver->prob != NULL, "Problem should exist");
    TEST_ASSERT(presolver->prob->obj->quad_qr != NULL, "QR should exist");
    
    free_presolver(presolver);
    free_settings(stgs);
    TEST_PASS_MSG("QP problem creation works");
}

/* Test 3: Problem with settings modification */
test_result_t test_problem_modify(void)
{
    printf("\n[Test] Problem Settings Modification\n");
    
    Settings *stgs = default_settings();
    TEST_ASSERT(stgs != NULL, "Settings should be created");
    
    /* Modify settings */
    stgs->verbose = false;
    stgs->dual_fix = true;
    stgs->max_time = 60.0;
    
    /* Create a simple problem */
    int n = 2;
    int m = 1;
    double c[] = {1.0, 2.0};
    double lb[] = {0.0, 0.0};
    double ub[] = {10.0, 10.0};
    double A_x[] = {1.0, 1.0};
    int A_i[] = {0, 0};
    int A_p[] = {0, 2};
    double lhs[] = {1.0};
    double rhs[] = {1.0};
    
    Presolver *presolver = new_presolver(
        A_x, A_i, A_p, m, n, 2,
        lhs, rhs, lb, ub, c, stgs
    );
    
    TEST_ASSERT(presolver != NULL, "Presolver should be created");
    TEST_ASSERT(presolver->stgs->dual_fix == true, "Settings should be applied");
    
    free_presolver(presolver);
    free_settings(stgs);
    TEST_PASS_MSG("Problem settings modification works");
}
