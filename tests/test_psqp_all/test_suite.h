/*
 * PSQP Test Suite - Header File
 */

#ifndef TEST_SUITE_H
#define TEST_SUITE_H

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include "Problem.h"
#include "Memory_wrapper.h"
#include "Bounds.h"

/* Forward declarations for internal functions not in public headers */
extern bool has_only_q_diag(const QuadTermQR *quad_qr, int k);
extern bool compute_Px_bounds(const QuadTermQR *qr, int k, 
                               const Bound *bounds, size_t n_cols,
                               double *min_Px, double *max_Px);

/* Test result codes */
typedef enum {
    TEST_PASS = 0,
    TEST_FAIL = 1,
    TEST_SKIP = 2
} test_result_t;

/* Test statistics */
typedef struct {
    size_t total;
    size_t passed;
    size_t failed;
    size_t skipped;
} test_stats_t;

extern test_stats_t g_stats;

/* Test assertion macros */
#define TEST_ASSERT(cond, msg) do { \
    if (!(cond)) { \
        printf("    ❌ FAIL: %s (line %d)\n", msg, __LINE__); \
        g_stats.failed++; \
        g_stats.total++; \
        return TEST_FAIL; \
    } \
} while(0)

#define TEST_ASSERT_NEAR(a, b, tol, msg) do { \
    if (fabs((a) - (b)) > (tol)) { \
        printf("    ❌ FAIL: %s (expected %g, got %g, line %d)\n", msg, (double)(b), (double)(a), __LINE__); \
        g_stats.failed++; \
        g_stats.total++; \
        return TEST_FAIL; \
    } \
} while(0)

#define TEST_PASS_MSG(msg) do { \
    printf("    ✅ PASS: %s\n", msg); \
    g_stats.passed++; \
    g_stats.total++; \
    return TEST_PASS; \
} while(0)

#define TEST_SKIP_MSG(msg) do { \
    printf("    ○ SKIP: %s\n", msg); \
    g_stats.skipped++; \
    g_stats.total++; \
    return TEST_SKIP; \
} while(0)

/* Test function declarations - QR Core */
test_result_t test_qr_basic(void);
test_result_t test_qr_collinear_merge(void);
test_result_t test_qr_isolated_vars(void);
test_result_t test_qr_dton_diagonal(void);

/* Test function declarations - Problem Construction */
test_result_t test_problem_create_lp(void);
test_result_t test_problem_create_qp(void);
test_result_t test_problem_modify(void);

/* Test function declarations - Presolve */
test_result_t test_presolve_empty_rows(void);
test_result_t test_presolve_ston_rows(void);
test_result_t test_presolve_dton_eq(void);
test_result_t test_presolve_dual_fix(void);
test_result_t test_presolve_bounds(void);

/* Test function declarations - QP Specific */
test_result_t test_qp_substitution(void);
test_result_t test_qp_dual_fix_bounds(void);
test_result_t test_qp_presolve_integration(void);
test_result_t test_qp_step4(void);
test_result_t test_qp_step5(void);
test_result_t test_qp_step6(void);
test_result_t test_qp_step7(void);
test_result_t test_qp_step8(void);
test_result_t test_qp_full_workflow(void);
test_result_t test_qp_with_dton_eq(void);

/* Test function declarations - Edge Cases */
test_result_t test_edge_empty_problem(void);
test_result_t test_edge_infeasible(void);
test_result_t test_edge_unbounded(void);
test_result_t test_edge_numerical(void);
test_result_t test_edge_no_constraints(void);
test_result_t test_edge_single_var(void);
test_result_t test_edge_equality_only(void);
test_result_t test_edge_duplicate_constraints(void);

/* Test function declarations - QP Postsolve */
test_result_t test_qp_dton_postsolve_simple(void);
test_result_t test_qp_dton_postsolve_dual(void);
test_result_t test_qp_isolated_postsolve(void);
test_result_t test_qp_postsolve_end2end(void);

/* Test function declarations - KKT Verification */
test_result_t test_kkt_simple_qp(void);
test_result_t test_kkt_active_bounds(void);
test_result_t test_kkt_multiple_constraints(void);
test_result_t test_kkt_inequality(void);
test_result_t test_kkt_fixed_vars(void);

/* Test function declarations - Partial Presolve */
test_result_t test_partial_dton_with_offdiag_q(void);
test_result_t test_partial_mixed_bounds(void);
test_result_t test_partial_multiple_dtons(void);
test_result_t test_partial_empty_and_ston(void);
test_result_t test_partial_with_r_matrix(void);

/* Test function declarations - QP Fixed Column Postsolve */
test_result_t test_qp_fixed_col_postsolve_simple(void);
test_result_t test_qp_fixed_col_postsolve_with_r(void);
test_result_t test_qp_fixed_col_postsolve_dual_recovery(void);

#endif /* TEST_SUITE_H */
