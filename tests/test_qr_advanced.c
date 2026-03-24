/*
 * Test for advanced QR optimization features:
 * 1. Collinear row detection and merging (integrated in quad_term_qr_new)
 * 2. Isolated quadratic variable fixing (fix_isolated_quadratic_vars_qp)
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "Problem.h"
#include "Memory_wrapper.h"
#include "Bounds.h"

/* Helper: Create a simple QuadTermQR for testing collinear row merging */
static QuadTermQR *create_test_qr_with_collinear_rows(void)
{
    /* Create QR with two collinear rows in R: r2 = 2*r1 */
    double Qx[] = {2.0, 4.0, 6.0};  /* Q diagonal */
    int Qi[] = {0, 1, 2};
    int Qp[] = {0, 1, 2, 3};
    size_t Qnnz = 3;
    
    /* R: row 0 = [1, 1, 0], row 1 = [2, 2, 0] = 2*row 0 */
    double Rx[] = {1.0, 1.0, 2.0, 2.0};
    int Ri[] = {0, 1, 0, 1};
    int Rp[] = {0, 2, 4};
    size_t Rnnz = 4;
    
    size_t n = 3;
    size_t k = 2;
    
    /* quad_term_qr_new will automatically call merge_collinear_rows_in_R */
    return quad_term_qr_new(Qx, Qi, Qp, Qnnz, Rx, Ri, Rp, Rnnz, n, k);
}

/* Test collinear row merging (now integrated in constructor) */
static int test_merge_collinear_rows(void)
{
    printf("Testing merge_collinear_rows_in_R (via quad_term_qr_new)...\n");
    
    QuadTermQR *qr = create_test_qr_with_collinear_rows();
    if (!qr) {
        printf("  FAIL: Failed to create test QR\n");
        return 1;
    }
    
    printf("  After constructor: k=%zu, Rnnz=%zu\n", qr->k, qr->Rnnz);
    
    /* Check results: k should be 1 (rows merged) */
    if (qr->k != 1) {
        printf("  FAIL: Expected k=1 after merging, got %zu\n", qr->k);
        free_quad_term_qr(qr);
        return 1;
    }
    
    if (qr->Rnnz != 2) {
        printf("  FAIL: Expected Rnnz=2 after merging, got %zu\n", qr->Rnnz);
        free_quad_term_qr(qr);
        return 1;
    }
    
    /* Check merged row values: r_new = sqrt(1^2 + 2^2) * [1, 1] = sqrt(5) * [1, 1] */
    double expected_scale = sqrt(5.0);
    if (fabs(qr->Rx[0] - expected_scale) > 1e-10 || 
        fabs(qr->Rx[1] - expected_scale) > 1e-10) {
        printf("  FAIL: Merged row values incorrect (expected %f, %f, got %f, %f)\n",
               expected_scale, expected_scale, qr->Rx[0], qr->Rx[1]);
        free_quad_term_qr(qr);
        return 1;
    }
    
    printf("  PASS: Collinear rows merged correctly in constructor\n");
    free_quad_term_qr(qr);
    return 0;
}

/* Test non-collinear rows are not merged */
static int test_no_false_merging(void)
{
    printf("\nTesting that non-collinear rows are not merged...\n");
    
    /* R: row 0 = [1, 1, 0], row 1 = [1, 2, 0] (not collinear) */
    double Qx[] = {2.0, 4.0, 6.0};
    int Qi[] = {0, 1, 2};
    int Qp[] = {0, 1, 2, 3};
    size_t Qnnz = 3;
    
    double Rx[] = {1.0, 1.0, 1.0, 2.0};
    int Ri[] = {0, 1, 0, 1};
    int Rp[] = {0, 2, 4};
    size_t Rnnz = 4;
    
    size_t n = 3;
    size_t k = 2;
    
    QuadTermQR *qr = quad_term_qr_new(Qx, Qi, Qp, Qnnz, Rx, Ri, Rp, Rnnz, n, k);
    if (!qr) {
        printf("  FAIL: Failed to create test QR\n");
        return 1;
    }
    
    if (qr->k != 2) {
        printf("  FAIL: Expected k=2 (no merging), got %zu\n", qr->k);
        free_quad_term_qr(qr);
        return 1;
    }
    
    printf("  PASS: Non-collinear rows correctly kept separate\n");
    free_quad_term_qr(qr);
    return 0;
}

/* Test single row in R (no merging possible) */
static int test_single_row_no_merge(void)
{
    printf("\nTesting single row in R (no merging possible)...\n");
    
    double Qx[] = {2.0};
    int Qi[] = {0};
    int Qp[] = {0, 1};
    size_t Qnnz = 1;
    
    double Rx[] = {1.0, 1.0};
    int Ri[] = {0, 1};
    int Rp[] = {0, 2};
    size_t Rnnz = 2;
    
    size_t n = 2;
    size_t k = 1;
    
    QuadTermQR *qr = quad_term_qr_new(Qx, Qi, Qp, Qnnz, Rx, Ri, Rp, Rnnz, n, k);
    if (!qr) {
        printf("  FAIL: Failed to create test QR\n");
        return 1;
    }
    
    if (qr->k != 1) {
        printf("  FAIL: Expected k=1, got %zu\n", qr->k);
        free_quad_term_qr(qr);
        return 1;
    }
    
    printf("  PASS: Single row preserved correctly\n");
    free_quad_term_qr(qr);
    return 0;
}

/* Test empty R (no merging) */
static int test_empty_r(void)
{
    printf("\nTesting empty R matrix...\n");
    
    double Qx[] = {2.0, 3.0};
    int Qi[] = {0, 1};
    int Qp[] = {0, 1, 2};
    size_t Qnnz = 2;
    
    size_t n = 2;
    size_t k = 0;
    
    QuadTermQR *qr = quad_term_qr_new(Qx, Qi, Qp, Qnnz, NULL, NULL, NULL, 0, n, k);
    if (!qr) {
        printf("  Note: quad_term_qr_new returned NULL for empty R (expected)\n");
        printf("  PASS: Empty R handled correctly\n");
        return 0;
    }
    
    if (qr->k != 0) {
        printf("  FAIL: Expected k=0, got %zu\n", qr->k);
        free_quad_term_qr(qr);
        return 1;
    }
    
    printf("  PASS: Empty R preserved correctly\n");
    free_quad_term_qr(qr);
    return 0;
}

int main(void)
{
    printf("=== Advanced QR Optimization Tests ===\n\n");
    
    int failures = 0;
    
    failures += test_merge_collinear_rows();
    failures += test_no_false_merging();
    failures += test_single_row_no_merge();
    failures += test_empty_r();
    
    printf("\n=== Summary ===\n");
    if (failures == 0) {
        printf("All tests passed!\n");
        return 0;
    } else {
        printf("%d test(s) failed!\n", failures);
        return 1;
    }
}
