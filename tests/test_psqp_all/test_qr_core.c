/*
 * QR Core Functionality Tests
 */

#include "test_suite.h"

/* Test 1: Basic QR construction */
test_result_t test_qr_basic(void)
{
    printf("\n[Test] QR Basic Construction\n");
    
    /* Test 1.1: Empty QR */
    {
        QuadTermQR *qr = quad_term_qr_new(NULL, NULL, NULL, 0, NULL, NULL, NULL, 0, 0, 0);
        TEST_ASSERT(qr == NULL, "Empty QR should return NULL");
    }
    
    /* Test 1.2: Q-only QR */
    {
        double Qx[] = {2.0, 3.0};
        int Qi[] = {0, 1};
        int Qp[] = {0, 1, 2};
        QuadTermQR *qr = quad_term_qr_new(Qx, Qi, Qp, 2, NULL, NULL, NULL, 0, 2, 0);
        TEST_ASSERT(qr != NULL, "Q-only QR should be created");
        TEST_ASSERT(qr->n == 2, "n should be 2");
        TEST_ASSERT(qr->k == 0, "k should be 0 (no R)");
        TEST_ASSERT(qr->Qnnz == 2, "Qnnz should be 2");
        free_quad_term_qr(qr);
    }
    
    /* Test 1.3: R-only QR */
    {
        double Rx[] = {1.0, 1.0};
        int Ri[] = {0, 1};
        int Rp[] = {0, 2};
        QuadTermQR *qr = quad_term_qr_new(NULL, NULL, NULL, 0, Rx, Ri, Rp, 2, 2, 1);
        TEST_ASSERT(qr != NULL, "R-only QR should be created");
        TEST_ASSERT(qr->k == 1, "k should be 1");
        free_quad_term_qr(qr);
    }
    
    /* Test 1.4: Full QR (Q + R with multi-element row)
     * Note: Single-element R rows get auto-extracted to Q diagonal
     */
    {
        double Qx[] = {2.0};
        int Qi[] = {0};
        int Qp[] = {0, 1, 1};
        /* Use 2 rows in R with 2 elements each so they don't get extracted */
        double Rx[] = {1.0, 1.0, 1.0, 1.0};
        int Ri[] = {0, 1, 0, 1};
        int Rp[] = {0, 2, 4};
        QuadTermQR *qr = quad_term_qr_new(Qx, Qi, Qp, 1, Rx, Ri, Rp, 4, 2, 2);
        TEST_ASSERT(qr != NULL, "Full QR should be created");
        TEST_ASSERT(qr->Qnnz >= 1, "Qnnz should be at least 1");
        TEST_ASSERT(qr->k >= 1, "k should be at least 1");
        TEST_ASSERT(qr->RTp != NULL, "R^T should be built");
        free_quad_term_qr(qr);
    }
    
    TEST_PASS_MSG("QR basic construction works");
}

/* Test 2: Collinear row merging */
test_result_t test_qr_collinear_merge(void)
{
    printf("\n[Test] QR Collinear Row Merging\n");
    
    /* Create QR with collinear rows: r2 = 2*r1 */
    double Qx[] = {2.0, 4.0, 6.0};
    int Qi[] = {0, 1, 2};
    int Qp[] = {0, 1, 2, 3};
    double Rx[] = {1.0, 1.0, 2.0, 2.0};  /* row0=[1,1], row1=[2,2]=2*row0 */
    int Ri[] = {0, 1, 0, 1};
    int Rp[] = {0, 2, 4};
    
    QuadTermQR *qr = quad_term_qr_new(Qx, Qi, Qp, 3, Rx, Ri, Rp, 4, 3, 2);
    TEST_ASSERT(qr != NULL, "QR should be created");
    TEST_ASSERT(qr->k == 1, "Rows should be merged to k=1");
    TEST_ASSERT(qr->Rnnz == 2, "Rnnz should be 2 after merge");
    
    /* Check merged values: sqrt(1^2 + 2^2) = sqrt(5) */
    double expected = sqrt(5.0);
    TEST_ASSERT_NEAR(qr->Rx[0], expected, 1e-10, "Merged row[0] incorrect");
    TEST_ASSERT_NEAR(qr->Rx[1], expected, 1e-10, "Merged row[1] incorrect");
    
    free_quad_term_qr(qr);
    TEST_PASS_MSG("Collinear row merging works");
}

/* Test 3: Isolated quadratic variable detection */
test_result_t test_qr_isolated_vars(void)
{
    printf("\n[Test] QR Isolated Variable Detection\n");
    
    
    /* Create QR:
     * - Q: diagonal with [2.0, 3.0]
     * - R: 1 row with 2 elements (affects both vars, won't be auto-extracted)
     */
    double Qx[] = {2.0, 3.0};
    int Qi[] = {0, 1};
    int Qp[] = {0, 1, 2};
    /* R affects both vars (multi-element, won't be extracted) */
    double Rx[] = {1.0, 1.0};
    int Ri[] = {0, 1};
    int Rp[] = {0, 2};
    
    QuadTermQR *qr = quad_term_qr_new(Qx, Qi, Qp, 2, Rx, Ri, Rp, 2, 2, 1);
    TEST_ASSERT(qr != NULL, "QR should be created");
    
    /* Both vars have R entries (R has 1 row with 2 elements) */
    TEST_ASSERT(!has_only_q_diag(qr, 0), "Var 0 should not be pure diagonal (has R)");
    TEST_ASSERT(!has_only_q_diag(qr, 1), "Var 1 should not be pure diagonal (has R)");
    
    /* Test pure Q diagonal case (no R) */
    QuadTermQR *qr2 = quad_term_qr_new(Qx, Qi, Qp, 2, NULL, NULL, NULL, 0, 2, 0);
    TEST_ASSERT(qr2 != NULL, "Q-only QR should be created");
    TEST_ASSERT(has_only_q_diag(qr2, 0), "Var 0 should be pure diagonal (Q only)");
    TEST_ASSERT(has_only_q_diag(qr2, 1), "Var 1 should be pure diagonal (Q only)");
    free_quad_term_qr(qr2);
    
    free_quad_term_qr(qr);
    TEST_PASS_MSG("Isolated variable detection works");
}

/* Test 4: DtonsEq with pure diagonal Q */
test_result_t test_qr_dton_diagonal(void)
{
    printf("\n[Test] QR DtonsEq Pure Diagonal Substitution\n");
    
    
    /* Create QR with pure diagonal Q for var 0, no R */
    double Qx[] = {2.0, 3.0};
    int Qi[] = {0, 1};
    int Qp[] = {0, 1, 2};
    
    QuadTermQR *qr = quad_term_qr_new(Qx, Qi, Qp, 2, NULL, NULL, NULL, 0, 2, 0);
    TEST_ASSERT(qr != NULL, "QR should be created");
    TEST_ASSERT(has_only_q_diag(qr, 0), "Var 0 should be pure diagonal");
    TEST_ASSERT(has_only_q_diag(qr, 1), "Var 1 should be pure diagonal");
    
    free_quad_term_qr(qr);
    TEST_PASS_MSG("DtonsEq pure diagonal check works");
}
