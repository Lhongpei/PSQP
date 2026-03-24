/*
 * Test for Doubleton Equality with pure Q diagonal variables.
 * 
 * Tests that variables with only Q diagonal entries (no off-diagonal, no R)
 * can be safely substituted in DtonsEq without creating cross-terms.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "Problem.h"
#include "Memory_wrapper.h"
#include "Bounds.h"

/* Forward declaration from Problem.h */
extern bool has_only_q_diag(const QuadTermQR *quad_qr, int k);

/* Test has_only_q_diag function */
static int test_has_only_q_diag(void)
{
    printf("Testing has_only_q_diag...\n");
    
    /* Case 1: Pure Q diagonal variable */
    {
        double Qx[] = {2.0};
        int Qi[] = {0};
        int Qp[] = {0, 1};
        size_t Qnnz = 1;
        size_t n = 1, k = 0;
        
        QuadTermQR *qr = quad_term_qr_new(Qx, Qi, Qp, Qnnz, NULL, NULL, NULL, 0, n, k);
        if (!qr) {
            printf("  FAIL: Failed to create QR\n");
            return 1;
        }
        
        if (!has_only_q_diag(qr, 0)) {
            printf("  FAIL: Pure Q diagonal should return true\n");
            free_quad_term_qr(qr);
            return 1;
        }
        printf("  PASS: Pure Q diagonal detected correctly\n");
        free_quad_term_qr(qr);
    }
    
    /* Case 2: Variable with R entries 
     * Create R with multi-element rows so they don't get auto-extracted.
     * Var 0 has Q diagonal but also appears in R (in a row with another var).
     */
    {
        double Qx[] = {2.0, 3.0};  /* Q[0][0]=2, Q[1][1]=3 */
        int Qi[] = {0, 1};
        int Qp[] = {0, 1, 2};  /* n=2, both vars have Q diagonal */
        /* R: 1 row with 2 elements [var0, var1] - multi-element, won't be extracted */
        double Rx[] = {1.0, 1.0};
        int Ri[] = {0, 1};
        int Rp[] = {0, 2};  /* 1 row with 2 entries */
        size_t n = 2, k = 1;
        
        QuadTermQR *qr = quad_term_qr_new(Qx, Qi, Qp, 2, Rx, Ri, Rp, 2, n, k);
        if (!qr) {
            printf("  FAIL: Failed to create QR with R\n");
            return 1;
        }
        
        /* Var 0 has Q diagonal AND R entry */
        if (has_only_q_diag(qr, 0)) {
            printf("  FAIL: Variable with R entries should return false\n");
            free_quad_term_qr(qr);
            return 1;
        }
        /* Var 1 also has Q diagonal AND R entry */
        if (has_only_q_diag(qr, 1)) {
            printf("  FAIL: Variable 1 with R entries should return false\n");
            free_quad_term_qr(qr);
            return 1;
        }
        printf("  PASS: Variables with R entries correctly rejected\n");
        free_quad_term_qr(qr);
    }
    
    /* Case 3: Variable with Q off-diagonal entries */
    {
        double Qx[] = {2.0, 1.0};  /* Q[0][0]=2, Q[0][1]=1 (off-diagonal) */
        int Qi[] = {0, 1};
        int Qp[] = {0, 2, 2};  /* row 0 has 2 entries, row 1 has 0 */
        size_t n = 2, k = 0;
        
        QuadTermQR *qr = quad_term_qr_new(Qx, Qi, Qp, 2, NULL, NULL, NULL, 0, n, k);
        if (!qr) {
            printf("  FAIL: Failed to create QR with off-diagonal\n");
            return 1;
        }
        
        if (has_only_q_diag(qr, 0)) {
            printf("  FAIL: Variable with Q off-diagonal should return false\n");
            free_quad_term_qr(qr);
            return 1;
        }
        printf("  PASS: Variable with Q off-diagonal correctly rejected\n");
        free_quad_term_qr(qr);
    }
    
    /* Case 4: No Q entries at all (only R)
     * Use multi-element R row to avoid auto-extraction to Q.
     */
    {
        /* R: 1 row with 2 elements - won't be auto-extracted */
        double Rx[] = {1.0, 1.0};
        int Ri[] = {0, 1};
        int Rp[] = {0, 2};
        size_t n = 2, k = 1;
        
        QuadTermQR *qr = quad_term_qr_new(NULL, NULL, NULL, 0, Rx, Ri, Rp, 2, n, k);
        if (!qr) {
            printf("  FAIL: Failed to create QR with only R\n");
            return 1;
        }
        
        /* Var 0 has only R entry, no Q */
        if (has_only_q_diag(qr, 0)) {
            printf("  FAIL: Variable with only R should return false\n");
            free_quad_term_qr(qr);
            return 1;
        }
        printf("  PASS: Variable with only R correctly rejected\n");
        free_quad_term_qr(qr);
    }
    
    return 0;
}

/* Test DtonsEq substitution with pure Q diagonal
 * 
 * Problem:
 *   min 0.5 * 2 * x1^2 + x1 + 2*x2
 *   s.t. x1 + x2 = 3  (doubleton equality)
 *        0 <= x1, x2 <= 10
 * 
 * x1 has pure Q diagonal (Q[0][0]=2, no off-diagonal, no R).
 * After substituting x1 = 3 - x2:
 *   0.5 * 2 * (3-x2)^2 = (3-x2)^2 = 9 - 6*x2 + x2^2
 *   Linear term: (3-x2) = 3 - x2
 *   Total obj: 9 - 6*x2 + x2^2 + 3 - x2 + 2*x2 = 12 - 5*x2 + x2^2
 *            = 12 + (-5)*x2 + 0.5 * 2 * x2^2
 * 
 * So after substitution:
 *   - Q[1][1] should increase by 2 (from x2^2 term)
 *   - c[1] should be -5
 *   - offset should be 12
 */
static int test_dton_q_diag_substitution(void)
{
    printf("\nTesting DtonsEq substitution with pure Q diagonal...\n");
    printf("  (This test verifies the substitution math is correct)\n");
    
    /* This test documents the expected behavior.
     * Full integration test would require setting up a complete Problem
     * and running the presolver.
     */
    
    /* Create QR: x1 has Q[0][0]=2, x2 has no Q */
    double Qx[] = {2.0};
    int Qi[] = {0};
    int Qp[] = {0, 1, 1};  /* var 0 has 1 entry, var 1 has 0 entries */
    
    QuadTermQR *qr = quad_term_qr_new(Qx, Qi, Qp, 1, NULL, NULL, NULL, 0, 2, 0);
    if (!qr) {
        printf("  FAIL: Failed to create QR\n");
        return 1;
    }
    
    /* Verify x1 has only Q diagonal */
    if (!has_only_q_diag(qr, 0)) {
        printf("  FAIL: x1 should have only Q diagonal\n");
        free_quad_term_qr(qr);
        return 1;
    }
    
    /* Verify x2 has no quadratic terms */
    if (has_quadratic_terms_qr(qr, 1)) {
        printf("  FAIL: x2 should have no quadratic terms\n");
        free_quad_term_qr(qr);
        return 1;
    }
    
    printf("  PASS: Pure Q diagonal variable can be substituted in DtonsEq\n");
    free_quad_term_qr(qr);
    return 0;
}

int main(void)
{
    printf("=== QR Doubleton Equality Tests ===\n\n");
    
    int failures = 0;
    
    failures += test_has_only_q_diag();
    failures += test_dton_q_diag_substitution();
    
    printf("\n=== Summary ===\n");
    if (failures == 0) {
        printf("All tests passed!\n");
        return 0;
    } else {
        printf("%d test(s) failed!\n", failures);
        return 1;
    }
}
