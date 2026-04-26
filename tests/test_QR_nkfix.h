/*
 * Focused tests for Q + R*R^T correctness.
 *
 * These tests hit paths that the earlier QR tests avoided — in particular,
 * multi-variable fixing which forces the presolver to read Rp[i] for i up to
 * n-1 and to evaluate P_ii = Q_ii + sum_j R_ij^2 for every fixed variable.
 *
 * With the original bug (Rp allocated as k+1 instead of n+1 and R stored as
 * k×n internally despite being advertised as n×k), these tests would read
 * out-of-bounds memory and/or return a mathematically wrong offset. They are
 * designed to pass only when R is stored consistently as n×k CSR.
 */

#ifndef TEST_QR_NKFIX_H
#define TEST_QR_NKFIX_H

#include "PSQP_API.h"
#include "test_macros.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define NKTOL 1e-9

/* Reference: P = Q + R*R^T for R stored as n×k CSR.
 * Returns the full dense P[row][col]. */
static double ref_P(const double *Qx, const int *Qi, const int *Qp,
                    const double *Rx, const int *Ri, const int *Rp,
                    size_t n, int row, int col)
{
    double v = 0.0;
    if (Qx && Qp)
    {
        int a = (row <= col) ? row : col;
        int b = (row <= col) ? col : row;
        for (int idx = Qp[a]; idx < Qp[a + 1]; idx++)
        {
            if (Qi[idx] == b) { v += Qx[idx]; break; }
        }
    }
    if (Rx && Rp)
    {
        int r_start_i = Rp[row];
        int r_end_i = Rp[row + 1];
        int r_start_j = Rp[col];
        int r_end_j = Rp[col + 1];
        for (int a = r_start_i; a < r_end_i; a++)
        {
            int fa = Ri[a];
            double va = Rx[a];
            for (int b = r_start_j; b < r_end_j; b++)
            {
                if (Ri[b] == fa) { v += va * Rx[b]; }
            }
        }
    }
    (void)n;
    return v;
}

/* Test: Fix two variables in a 3-variable QP, verify offset and c update.
 *
 * Problem:
 *   min 0.5 x^T (Q + R R^T) x + c^T x
 *   s.t. x0 + x1 + x2 = 3,   0 <= x_i
 *   Fix x0 = 1.0 and x1 = 2.0 via tight bounds (only x2 remains free).
 *
 * Q = diag(1, 1, 1).
 * R = [[1, 0], [0, 1], [1, 1]]  (n=3, k=2).
 * => P = Q + R R^T = [[2, 0, 1], [0, 2, 1], [1, 1, 3]]
 *
 * Expected after fixing:
 *   offset = 0.5 * P00 * 1^2 + 0.5 * P11 * 2^2 + P01 * 1 * 2
 *          = 0.5 * 2 + 0.5 * 2 * 4 + 0 = 1 + 4 = 5
 *   c2 += P02 * 1 + P12 * 2 = 1 + 2 = 3, so c2 becomes 3 + original.
 *
 * If the presolver reads Rp[2] out of bounds (original bug) the result would
 * be undefined; with the fix, it must hit exactly the reference value. */
static char *test_qr_two_fixed_vars_offdiag_R()
{
    size_t m = 1;
    size_t n = 3;

    double Ax[] = {1.0, 1.0, 1.0};
    int Ai[] = {0, 1, 2};
    int Ap[] = {0, 3};
    size_t nnz = 3;

    double lhs[] = {3.0};
    double rhs[] = {3.0};
    double lbs[] = {1.0, 2.0, 0.0};
    double ubs[] = {1.0, 2.0, INFINITY};
    double c[] = {0.0, 0.0, 0.0};

    double Qx[] = {1.0, 1.0, 1.0};
    int Qi[] = {0, 1, 2};
    int Qp[] = {0, 1, 2, 3};
    size_t Qnnz = 3;

    /* R (3x2):
     * row 0: [1, 0]
     * row 1: [0, 1]
     * row 2: [1, 1]
     */
    double Rx[] = {1.0, 1.0, 1.0, 1.0};
    int Ri[] = {0, 1, 0, 1};
    int Rp[] = {0, 1, 2, 4};
    size_t Rnnz = 4;
    size_t k = 2;

    double expected_offset =
        0.5 * ref_P(Qx, Qi, Qp, Rx, Ri, Rp, n, 0, 0) * 1.0 * 1.0 +
        0.5 * ref_P(Qx, Qi, Qp, Rx, Ri, Rp, n, 1, 1) * 2.0 * 2.0 +
        ref_P(Qx, Qi, Qp, Rx, Ri, Rp, n, 0, 1) * 1.0 * 2.0;
    double expected_c2 =
        ref_P(Qx, Qi, Qp, Rx, Ri, Rp, n, 0, 2) * 1.0 +
        ref_P(Qx, Qi, Qp, Rx, Ri, Rp, n, 1, 2) * 2.0;

    Settings *stgs = default_settings();
    stgs->verbose = false;

    Presolver *presolver = new_qp_presolver_qr(Ax, Ai, Ap, m, n, nnz,
                                               lhs, rhs, lbs, ubs, c,
                                               Qx, Qi, Qp, Qnnz,
                                               Rx, Ri, Rp, Rnnz, k, stgs);
    mu_assert("presolver creation failed", presolver != NULL);

    PresolveStatus status = run_presolver(presolver);
    PresolvedProblem *reduced = presolver->reduced_prob;

    mu_assert("presolve must not fail", status == REDUCED || status == UNCHANGED);

    if (fabs(reduced->obj_offset - expected_offset) >= NKTOL)
    {
        printf("\n    two_fixed: reduced_offset=%.9f expected=%.9f diff=%.3e\n",
               reduced->obj_offset, expected_offset, reduced->obj_offset - expected_offset);
        printf("    reduced->n=%zu\n", reduced->n);
        if (reduced->n == 1)
            printf("    reduced->c[0]=%.9f (expected %.9f)\n", reduced->c[0], expected_c2);
    }
    mu_assert("offset mismatch after fixing two vars with R coupling",
              fabs(reduced->obj_offset - expected_offset) < NKTOL);

    /* Reduced problem should have at most 1 variable remaining (x2). */
    if (reduced->n == 1)
    {
        mu_assert("reduced c[0] should equal expected c2 update",
                  fabs(reduced->c[0] - expected_c2) < NKTOL);
    }

    free_settings(stgs);
    free_presolver(presolver);
    return 0;
}

/* Test: n > k stress test that touches Rp[i] for i > k.
 * With the original bug (Rp allocated size k+1), this reads out-of-bounds
 * memory for every variable with index >= k that gets fixed.
 *
 * Setup: n=4, k=1, R = [1; 1; 1; 1]. Q = diag(2). With the equality
 * constraint sum(x)=4 and x0=x1=x2 fixed to 1, the presolver also fixes
 * x3=1, so all four variables end up fixed.
 *
 * The objective at (1,1,1,1) is 0.5 * 1^T (Q + R R^T) 1 = 0.5 * 24 = 12.
 * That value must appear exactly in obj_offset; any off-by-one in Rp or any
 * missed R cross-term contribution would show up as a difference here. */
static char *test_qr_n_gt_k_multifix()
{
    size_t m = 1;
    size_t n = 4;

    double Ax[] = {1.0, 1.0, 1.0, 1.0};
    int Ai[] = {0, 1, 2, 3};
    int Ap[] = {0, 4};
    size_t nnz = 4;

    double lhs[] = {4.0};
    double rhs[] = {4.0};
    double lbs[] = {1.0, 1.0, 1.0, 0.0};
    double ubs[] = {1.0, 1.0, 1.0, INFINITY};
    double c[] = {0.0, 0.0, 0.0, 0.0};

    double Qx[] = {2.0, 2.0, 2.0, 2.0};
    int Qi[] = {0, 1, 2, 3};
    int Qp[] = {0, 1, 2, 3, 4};
    size_t Qnnz = 4;

    double Rx[] = {1.0, 1.0, 1.0, 1.0};
    int Ri[] = {0, 0, 0, 0};
    int Rp[] = {0, 1, 2, 3, 4};
    size_t Rnnz = 4;
    size_t k = 1;

    /* Compute expected offset: 0.5 x^T P x + c^T x at x=(1,1,1,1). */
    double final_x[4] = {1.0, 1.0, 1.0, 1.0};
    double expected_offset = 0.0;
    for (int i = 0; i < 4; i++)
    {
        expected_offset += c[i] * final_x[i];
        expected_offset += 0.5 * ref_P(Qx, Qi, Qp, Rx, Ri, Rp, n, i, i)
                           * final_x[i] * final_x[i];
        for (int j = i + 1; j < 4; j++)
        {
            expected_offset += ref_P(Qx, Qi, Qp, Rx, Ri, Rp, n, i, j)
                               * final_x[i] * final_x[j];
        }
    }

    Settings *stgs = default_settings();
    stgs->verbose = false;

    Presolver *presolver = new_qp_presolver_qr(Ax, Ai, Ap, m, n, nnz,
                                               lhs, rhs, lbs, ubs, c,
                                               Qx, Qi, Qp, Qnnz,
                                               Rx, Ri, Rp, Rnnz, k, stgs);
    mu_assert("presolver creation failed", presolver != NULL);

    PresolveStatus status = run_presolver(presolver);
    PresolvedProblem *reduced = presolver->reduced_prob;

    mu_assert("presolve must not fail", status == REDUCED || status == UNCHANGED);

    if (fabs(reduced->obj_offset - expected_offset) >= NKTOL)
    {
        printf("\n    n_gt_k_multifix: reduced_offset=%.9f expected=%.9f diff=%.3e\n",
               reduced->obj_offset, expected_offset, reduced->obj_offset - expected_offset);
        printf("    reduced->n=%zu reduced->m=%zu\n", reduced->n, reduced->m);
    }
    mu_assert("offset should equal 0.5 x^T P x at the implied solution",
              fabs(reduced->obj_offset - expected_offset) < NKTOL);

    free_settings(stgs);
    free_presolver(presolver);
    return 0;
}

/* Test: Output structure invariants.
 * After presolving a QP where no variables are eliminated, the output Rp
 * must have n+1 entries and Rp[n] must equal Rnnz. This directly verifies
 * that Rp is allocated and populated correctly for n×k storage. */
static char *test_qr_output_rp_size()
{
    size_t m = 1;
    size_t n = 5;

    double Ax[] = {1.0, 1.0, 1.0, 1.0, 1.0};
    int Ai[] = {0, 1, 2, 3, 4};
    int Ap[] = {0, 5};
    size_t nnz = 5;

    double lhs[] = {-INFINITY};
    double rhs[] = {INFINITY};
    double lbs[] = {0.0, 0.0, 0.0, 0.0, 0.0};
    double ubs[] = {INFINITY, INFINITY, INFINITY, INFINITY, INFINITY};
    double c[] = {1.0, 1.0, 1.0, 1.0, 1.0};

    double Qx[] = {1.0, 1.0, 1.0, 1.0, 1.0};
    int Qi[] = {0, 1, 2, 3, 4};
    int Qp[] = {0, 1, 2, 3, 4, 5};
    size_t Qnnz = 5;

    /* R is 5x2 with one nonzero per row */
    double Rx[] = {0.5, 0.4, 0.3, 0.2, 0.1};
    int Ri[] = {0, 1, 0, 1, 0};
    int Rp[] = {0, 1, 2, 3, 4, 5};
    size_t Rnnz = 5;
    size_t k = 2;

    Settings *stgs = default_settings();
    stgs->verbose = false;
    /* Disable aggressive reductions so we can inspect the output */
    stgs->ston_cols = false;
    stgs->dton_eq = false;
    stgs->parallel_rows = false;
    stgs->parallel_cols = false;
    stgs->primal_propagation = false;
    stgs->dual_fix = false;

    Presolver *presolver = new_qp_presolver_qr(Ax, Ai, Ap, m, n, nnz,
                                               lhs, rhs, lbs, ubs, c,
                                               Qx, Qi, Qp, Qnnz,
                                               Rx, Ri, Rp, Rnnz, k, stgs);
    mu_assert("presolver creation failed", presolver != NULL);
    run_presolver(presolver);
    PresolvedProblem *reduced = presolver->reduced_prob;

    if (reduced->Rnnz > 0)
    {
        mu_assert("Rp[0] must be 0", reduced->Rp[0] == 0);
        mu_assert("Rp[n] must equal Rnnz",
                  (size_t)reduced->Rp[reduced->n] == reduced->Rnnz);
        /* Monotone */
        for (size_t i = 0; i < reduced->n; i++)
        {
            mu_assert("Rp must be monotone nondecreasing",
                      reduced->Rp[i + 1] >= reduced->Rp[i]);
        }
        /* Column indices within [0, k) */
        for (size_t i = 0; i < reduced->Rnnz; i++)
        {
            mu_assert("Ri must be in [0, k)",
                      reduced->Ri[i] >= 0 && (size_t)reduced->Ri[i] < reduced->k);
        }
    }

    free_settings(stgs);
    free_presolver(presolver);
    return 0;
}

/* Test: a factor whose only non-zero loads a single variable must be
 * absorbed into the Q diagonal (rank-reduction optimization).
 *
 * R is 3×2:
 *   factor 0: R[0,0] = 2 (single-nnz factor → extract to Q[0,0] += 4)
 *   factor 1: R[0,1] = 1, R[1,1] = 1 (two-variable factor, stays)
 *
 * After extraction we expect k == 1 and Q[0,0] absorbed += 4.
 * The overall P matrix must still be the original Q + R R^T. */
static char *test_qr_extract_single_nnz_factor()
{
    size_t m = 1;
    size_t n = 3;

    double Ax[] = {1.0, 1.0, 1.0};
    int Ai[] = {0, 1, 2};
    int Ap[] = {0, 3};
    size_t nnz = 3;

    double lhs[] = {-INFINITY};
    double rhs[] = {INFINITY};
    double lbs[] = {0.0, 0.0, 0.0};
    double ubs[] = {INFINITY, INFINITY, INFINITY};
    double c[] = {0.0, 0.0, 0.0};

    double Qx[] = {1.0, 1.0, 1.0};
    int Qi[] = {0, 1, 2};
    int Qp[] = {0, 1, 2, 3};
    size_t Qnnz = 3;

    /* R (3x2): factor 0 single-nnz at var 0; factor 1 touches vars 0,1. */
    double Rx[] = {2.0, 1.0, 1.0};
    int Ri[] = {0, 1, 1};         /* col indices inside each row */
    int Rp[] = {0, 2, 3, 3};      /* row 0: two entries; row 1: one; row 2: none */
    size_t Rnnz = 3;
    size_t k = 2;

    Settings *stgs = default_settings();
    stgs->verbose = false;
    stgs->ston_cols = false;
    stgs->dton_eq = false;
    stgs->parallel_rows = false;
    stgs->parallel_cols = false;
    stgs->primal_propagation = false;
    stgs->dual_fix = false;

    Presolver *presolver = new_qp_presolver_qr(Ax, Ai, Ap, m, n, nnz,
                                               lhs, rhs, lbs, ubs, c,
                                               Qx, Qi, Qp, Qnnz,
                                               Rx, Ri, Rp, Rnnz, k, stgs);
    mu_assert("presolver creation failed", presolver != NULL);
    run_presolver(presolver);
    PresolvedProblem *reduced = presolver->reduced_prob;

    /* Expected diagonal values after extraction */
    double expected_diag0 = 1.0 + 2.0 * 2.0 + 1.0 * 1.0;  /* Q0 + R00^2 + R01^2 */
    double expected_diag1 = 1.0 + 1.0 * 1.0;              /* Q1 + R11^2 */
    double expected_diag2 = 1.0;                          /* Q2 */

    /* Compute dense P diagonal from reduced Q + reduced R R^T. */
    for (size_t v = 0; v < reduced->n; v++)
    {
        double p_vv = 0.0;
        if (reduced->Qp)
        {
            for (int idx = reduced->Qp[v]; idx < reduced->Qp[v + 1]; idx++)
            {
                if (reduced->Qi[idx] == (int)v) { p_vv += reduced->Qx[idx]; break; }
            }
        }
        if (reduced->Rp)
        {
            for (int idx = reduced->Rp[v]; idx < reduced->Rp[v + 1]; idx++)
            {
                p_vv += reduced->Rx[idx] * reduced->Rx[idx];
            }
        }
        double expect = (v == 0) ? expected_diag0
                       : (v == 1) ? expected_diag1 : expected_diag2;
        mu_assert("reduced P diagonal matches original",
                  fabs(p_vv - expect) < NKTOL);
    }

    /* After extraction the single-nnz factor is absorbed, so k drops to 1. */
    mu_assert("single-nnz factor should be extracted, k becomes 1",
              reduced->k <= 1);

    free_settings(stgs);
    free_presolver(presolver);
    return 0;
}

/* Test: two collinear factors should be merged, reducing k by 1 while
 * keeping R R^T numerically unchanged.
 *
 * R is 2×2 with col(0) = (1, 2)^T and col(1) = (2, 4)^T = 2 * col(0).
 * These two factors are collinear: R R^T from them equals
 *   col0 col0^T + col1 col1^T = (1 + 4) * col0 col0^T
 * After merging we expect k == 1 and the surviving factor scaled by sqrt(5).
 * Off-diagonal entries of P must match the pre-merge value. */
static char *test_qr_merge_collinear_factors()
{
    size_t m = 1;
    size_t n = 2;

    double Ax[] = {1.0, 1.0};
    int Ai[] = {0, 1};
    int Ap[] = {0, 2};
    size_t nnz = 2;

    double lhs[] = {-INFINITY};
    double rhs[] = {INFINITY};
    double lbs[] = {0.0, 0.0};
    double ubs[] = {INFINITY, INFINITY};
    double c[] = {0.0, 0.0};

    /* R (2x2): col0 = (1, 2), col1 = (2, 4). Collinear (col1 = 2*col0). */
    double Rx[] = {1.0, 2.0, 2.0, 4.0};
    int Ri[] = {0, 1, 0, 1};           /* factor indices in CSR row-major */
    int Rp[] = {0, 2, 4};              /* row 0: 2 entries; row 1: 2 entries */
    size_t Rnnz = 4;
    size_t k = 2;

    Settings *stgs = default_settings();
    stgs->verbose = false;
    stgs->ston_cols = false;
    stgs->dton_eq = false;
    stgs->parallel_rows = false;
    stgs->parallel_cols = false;
    stgs->primal_propagation = false;
    stgs->dual_fix = false;

    Presolver *presolver = new_qp_presolver_qr(Ax, Ai, Ap, m, n, nnz,
                                               lhs, rhs, lbs, ubs, c,
                                               NULL, NULL, NULL, 0,
                                               Rx, Ri, Rp, Rnnz, k, stgs);
    mu_assert("presolver creation failed", presolver != NULL);
    run_presolver(presolver);
    PresolvedProblem *reduced = presolver->reduced_prob;

    /* Original R R^T:
     *   P[0,0] = 1*1 + 2*2 = 5
     *   P[1,1] = 2*2 + 4*4 = 20
     *   P[0,1] = 1*2 + 2*4 = 10
     */
    double expected[2][2] = {{5.0, 10.0}, {10.0, 20.0}};

    mu_assert("k should reduce to 1 after merging collinear factors",
              reduced->k == 1 || reduced->k == 0);

    for (size_t i = 0; i < reduced->n; i++)
    {
        for (size_t j = 0; j < reduced->n; j++)
        {
            double p_ij = 0.0;
            if (reduced->Qp)
            {
                int a = i < j ? (int)i : (int)j;
                int b = i < j ? (int)j : (int)i;
                for (int idx = reduced->Qp[a]; idx < reduced->Qp[a + 1]; idx++)
                {
                    if (reduced->Qi[idx] == b) { p_ij += reduced->Qx[idx]; break; }
                }
            }
            if (reduced->Rp)
            {
                for (int ai = reduced->Rp[i]; ai < reduced->Rp[i + 1]; ai++)
                {
                    for (int bj = reduced->Rp[j]; bj < reduced->Rp[j + 1]; bj++)
                    {
                        if (reduced->Ri[ai] == reduced->Ri[bj])
                            p_ij += reduced->Rx[ai] * reduced->Rx[bj];
                    }
                }
            }
            mu_assert("merged P entry matches the pre-merge value",
                      fabs(p_ij - expected[i][j]) < 1e-7);
        }
    }

    free_settings(stgs);
    free_presolver(presolver);
    return 0;
}

/* Test: parallel_cols must NOT merge two variables that share quadratic
 * structure. This reproduces the bug where ston_cols + parallel_cols caused
 * spurious infeasibility for a QP with equality-constrained QP variables. */
static char *test_qr_parallel_cols_skips_qp_vars()
{
    size_t m = 1;
    size_t n = 3;

    /* x0 + x1 + x2 = 3 (equality) */
    double Ax[] = {1.0, 1.0, 1.0};
    int Ai[] = {0, 1, 2};
    int Ap[] = {0, 3};
    size_t nnz = 3;

    double lhs[] = {3.0};
    double rhs[] = {3.0};
    /* Fix x1 via tight bounds; x0, x2 remain parallel columns in the
     * residual row but both have quadratic structure. */
    double lbs[] = {0.0, 1.0, 0.0};
    double ubs[] = {INFINITY, 1.0, INFINITY};
    double c[] = {1.0, 2.0, 1.0};

    double Qx[] = {2.0, 3.0, 2.0};
    int Qi[] = {0, 1, 2};
    int Qp[] = {0, 1, 2, 3};
    size_t Qnnz = 3;

    double Rx[] = {1.0, 1.0, 1.0};
    int Ri[] = {0, 0, 0};
    int Rp[] = {0, 1, 2, 3};
    size_t Rnnz = 3;
    size_t k = 1;

    /* All defaults on: ston_cols, parallel_cols, primal_propagation, etc.
     * This is the exact configuration that previously reported INFEASIBLE. */
    Settings *stgs = default_settings();
    stgs->verbose = false;

    Presolver *presolver = new_qp_presolver_qr(Ax, Ai, Ap, m, n, nnz,
                                               lhs, rhs, lbs, ubs, c,
                                               Qx, Qi, Qp, Qnnz,
                                               Rx, Ri, Rp, Rnnz, k, stgs);
    mu_assert("presolver creation failed", presolver != NULL);

    PresolveStatus status = run_presolver(presolver);
    mu_assert("feasible QP must NOT be reported infeasible",
              status != INFEASIBLE && status != UNBNDORINFEAS);

    free_settings(stgs);
    free_presolver(presolver);
    return 0;
}

static int counter_qr_nkfix = 0;

static const char *all_tests_qr_nkfix()
{
    mu_run_test(test_qr_two_fixed_vars_offdiag_R, counter_qr_nkfix);
    mu_run_test(test_qr_n_gt_k_multifix, counter_qr_nkfix);
    mu_run_test(test_qr_output_rp_size, counter_qr_nkfix);
    mu_run_test(test_qr_extract_single_nnz_factor, counter_qr_nkfix);
    mu_run_test(test_qr_merge_collinear_factors, counter_qr_nkfix);
    mu_run_test(test_qr_parallel_cols_skips_qp_vars, counter_qr_nkfix);
    return 0;
}

int test_qr_nkfix()
{
    const char *result = all_tests_qr_nkfix();
    if (result != 0)
    {
        printf("%s\n", result);
        printf("QR nkfix: TEST FAILED!\n");
    }
    else
    {
        printf("QR nkfix: ALL TESTS PASSED\n");
    }
    printf("QR nkfix: Tests run: %d\n", counter_qr_nkfix);
    return result == 0;
}

#endif /* TEST_QR_NKFIX_H */
