/*
 * KKT verification after postsolve.
 *
 * For the original problem
 *   min 0.5 x^T P x + c^T x,  lhs <= Ax <= rhs,  lb <= x <= ub
 * PSQP's dual convention (from Postsolver.c) is:
 *   z = c + P x - A^T y
 * with z acting as the bound multiplier and y as the row multiplier.
 *
 * These tests construct a small QP, analytically solve the REDUCED problem
 * produced by PSQP, feed that reduced KKT point into postsolve(), and then
 * verify the full KKT conditions for the ORIGINAL problem:
 *   1. Primal feasibility (Ax bounds, variable bounds)
 *   2. Stationarity: c + Px - A^T y - z = 0
 *   3. Complementary slackness for active bounds
 *
 * If presolve or postsolve mishandles the QR decomposition (e.g. forgets
 * off-diagonal R R^T contributions), the stationarity residual on the
 * recovered z will be non-zero.
 */

#ifndef TEST_QR_KKT_H
#define TEST_QR_KKT_H

#include "PSQP_API.h"
#include "PSQP_sol.h"
#include "test_macros.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define KKT_TOL 1e-7

/* Dense-P helper: compute P[i][j] from (Q, R) in user coordinates
 * (Q upper-triangular CSR, R in n×k CSR). */
static double kkt_dense_P(const double *Qx, const int *Qi, const int *Qp,
                          const double *Rx, const int *Ri, const int *Rp,
                          int i, int j)
{
    double v = 0.0;
    if (Qx && Qp)
    {
        int a = i <= j ? i : j;
        int b = i <= j ? j : i;
        for (int idx = Qp[a]; idx < Qp[a + 1]; idx++)
        {
            if (Qi[idx] == b) { v += Qx[idx]; break; }
        }
    }
    if (Rx && Rp)
    {
        for (int ai = Rp[i]; ai < Rp[i + 1]; ai++)
        {
            for (int bj = Rp[j]; bj < Rp[j + 1]; bj++)
            {
                if (Ri[ai] == Ri[bj]) v += Rx[ai] * Rx[bj];
            }
        }
    }
    return v;
}

/* Verify KKT residuals for the original problem given (x, y, z). */
static char *kkt_assert(size_t n, size_t m,
                        const double *Ax, const int *Ai, const int *Ap,
                        const double *lhs, const double *rhs,
                        const double *lbs, const double *ubs,
                        const double *c,
                        const double *Qx, const int *Qi, const int *Qp,
                        const double *Rx, const int *Ri, const int *Rp,
                        const double *x, const double *y, const double *z)
{
    /* Primal feasibility: Ax */
    for (size_t i = 0; i < m; i++)
    {
        double axi = 0.0;
        for (int idx = Ap[i]; idx < Ap[i + 1]; idx++)
        {
            axi += Ax[idx] * x[Ai[idx]];
        }
        if (!isinf(lhs[i]))
            mu_assert("primal: Ax >= lhs violated", axi >= lhs[i] - KKT_TOL);
        if (!isinf(rhs[i]))
            mu_assert("primal: Ax <= rhs violated", axi <= rhs[i] + KKT_TOL);
    }

    /* Primal feasibility: bounds */
    for (size_t j = 0; j < n; j++)
    {
        if (!isinf(lbs[j]))
            mu_assert("primal: x >= lb violated", x[j] >= lbs[j] - KKT_TOL);
        if (!isinf(ubs[j]))
            mu_assert("primal: x <= ub violated", x[j] <= ubs[j] + KKT_TOL);
    }

    /* Stationarity: c + P x - A^T y - z = 0 for each variable j. */
    for (size_t j = 0; j < n; j++)
    {
        double stat = c[j];
        for (size_t i = 0; i < n; i++)
        {
            stat += kkt_dense_P(Qx, Qi, Qp, Rx, Ri, Rp, (int)j, (int)i) * x[i];
        }
        for (size_t i = 0; i < m; i++)
        {
            for (int idx = Ap[i]; idx < Ap[i + 1]; idx++)
            {
                if (Ai[idx] == (int)j) { stat -= Ax[idx] * y[i]; break; }
            }
        }
        stat -= z[j];
        if (fabs(stat) >= KKT_TOL)
        {
            fprintf(stderr, "  KKT stationarity violated at j=%zu: residual=%.3e\n",
                    j, stat);
            fprintf(stderr, "    x=[");
            for (size_t ii = 0; ii < n; ii++) fprintf(stderr, "%.6f ", x[ii]);
            fprintf(stderr, "] y=[");
            for (size_t ii = 0; ii < m; ii++) fprintf(stderr, "%.6f ", y[ii]);
            fprintf(stderr, "] z=[");
            for (size_t ii = 0; ii < n; ii++) fprintf(stderr, "%.6f ", z[ii]);
            fprintf(stderr, "]\n");
        }
        mu_assert("stationarity c + Px - A^T y - z != 0",
                  fabs(stat) < KKT_TOL);
    }

    /* Complementary slackness: z[j] * (bound slack) = 0.
     * When both bounds finite and equal (fixed var), z is unrestricted.
     * When lb active (x==lb): z must have the "right" sign for min problem (z >= 0).
     * When ub active (x==ub): z <= 0.
     * When strictly interior: z == 0. */
    for (size_t j = 0; j < n; j++)
    {
        bool lb_finite = !isinf(lbs[j]);
        bool ub_finite = !isinf(ubs[j]);
        bool lb_active = lb_finite && fabs(x[j] - lbs[j]) < KKT_TOL;
        bool ub_active = ub_finite && fabs(x[j] - ubs[j]) < KKT_TOL;

        if (lb_finite && ub_finite && fabs(ubs[j] - lbs[j]) < KKT_TOL)
        {
            /* Fixed variable: z is unrestricted. Skip sign test. */
            continue;
        }

        if (!lb_active && !ub_active)
        {
            mu_assert("z must be zero for strictly interior x",
                      fabs(z[j]) < KKT_TOL);
        }
        else if (lb_active && !ub_active)
        {
            mu_assert("z >= 0 when lb is active (min)",
                      z[j] >= -KKT_TOL);
        }
        else if (ub_active && !lb_active)
        {
            mu_assert("z <= 0 when ub is active (min)",
                      z[j] <= KKT_TOL);
        }
    }

    return 0;
}

/* ------------------------------------------------------------------
 * Test A: QP that gets fully determined by presolve (0 reduced vars).
 *
 *   min 0.5*(x0^2 + x1^2)
 *   s.t. x0 + x1 = 2, x0 = 1 (tight bounds), x1 >= 0
 *
 * Presolve fixes x0=1, then the singleton row forces x1=1.
 * Expected full KKT: x=(1,1), y=(1), z=(0,0).
 * (Because P*x = (1,1), A^T y = (1,1), no active bounds except the fixed
 * x0 whose z is unrestricted — so z=0 is a valid certificate.)
 * ------------------------------------------------------------------ */
static char *test_kkt_q_only_fix_plus_ston()
{
    size_t m = 1;
    size_t n = 2;

    double Ax[] = {1.0, 1.0};
    int Ai[] = {0, 1};
    int Ap[] = {0, 2};

    double lhs[] = {2.0};
    double rhs[] = {2.0};
    double lbs[] = {1.0, 0.0};
    double ubs[] = {1.0, INFINITY};
    double c[] = {0.0, 0.0};

    double Qx[] = {1.0, 1.0};
    int Qi[] = {0, 1};
    int Qp[] = {0, 1, 2};
    size_t Qnnz = 2;

    Settings *stgs = default_settings();
    stgs->verbose = false;

    Presolver *p = new_qp_presolver_qr(Ax, Ai, Ap, m, n, 2,
                                       lhs, rhs, lbs, ubs, c,
                                       Qx, Qi, Qp, Qnnz,
                                       NULL, NULL, NULL, 0, 0, stgs);
    mu_assert("presolver creation failed", p != NULL);
    PresolveStatus s = run_presolver(p);
    mu_assert("presolve must not fail", s != INFEASIBLE && s != UNBNDORINFEAS);

    PresolvedProblem *red = p->reduced_prob;
    /* Reduced may be 0-var 0-row. Allocate anyway. */
    double *xr = (double *)calloc(red->n > 0 ? red->n : 1, sizeof(double));
    double *yr = (double *)calloc(red->m > 0 ? red->m : 1, sizeof(double));
    double *zr = (double *)calloc(red->n > 0 ? red->n : 1, sizeof(double));

    postsolve(p, xr, yr, zr);

    char *err = kkt_assert(n, m, Ax, Ai, Ap, lhs, rhs, lbs, ubs, c,
                           Qx, Qi, Qp, NULL, NULL, NULL,
                           p->sol->x, p->sol->y, p->sol->z);

    free(xr); free(yr); free(zr);
    free_settings(stgs);
    free_presolver(p);
    return err;
}

/* ------------------------------------------------------------------
 * Test B: QP with R having off-diagonal R R^T contribution, fixed var.
 *
 *   min 0.5 x^T (Q + R R^T) x
 *   s.t. x0 + x1 = 2, x0 = 1 (tight bounds), x1 >= 0
 *
 * Q = diag(1,1), R = [1; 1]   =>  P = [[2, 1], [1, 2]].
 * Fixing x0 = 1 forces x1 = 1 via the constraint.
 *
 * Expected KKT:
 *   Stationarity at x1 (free): 2*1 + 1*1 - 1*y - z1 = 0. Interior x1=1 > 0 so z1=0 => y = 3.
 *   Stationarity at x0 (fixed): 2*1 + 1*1 - 1*y - z0 = 0 => z0 = 3 - y = 0.
 *
 * With the pre-fix postsolve (diagonal-only R R^T), z0 would be off by
 * the missed off-diagonal contribution R[0,:]·R[1,:]·x1.
 * ------------------------------------------------------------------ */
static char *test_kkt_R_offdiag_fix()
{
    size_t m = 1;
    size_t n = 2;

    double Ax[] = {1.0, 1.0};
    int Ai[] = {0, 1};
    int Ap[] = {0, 2};

    double lhs[] = {2.0};
    double rhs[] = {2.0};
    double lbs[] = {1.0, 0.0};
    double ubs[] = {1.0, INFINITY};
    double c[] = {0.0, 0.0};

    double Qx[] = {1.0, 1.0};
    int Qi[] = {0, 1};
    int Qp[] = {0, 1, 2};
    size_t Qnnz = 2;

    /* R = [1; 1] (2x1) */
    double Rx[] = {1.0, 1.0};
    int Ri[] = {0, 0};
    int Rp[] = {0, 1, 2};
    size_t Rnnz = 2;
    size_t k = 1;

    Settings *stgs = default_settings();
    stgs->verbose = false;

    Presolver *p = new_qp_presolver_qr(Ax, Ai, Ap, m, n, 2,
                                       lhs, rhs, lbs, ubs, c,
                                       Qx, Qi, Qp, Qnnz,
                                       Rx, Ri, Rp, Rnnz, k, stgs);
    mu_assert("presolver creation failed", p != NULL);
    run_presolver(p);

    PresolvedProblem *red = p->reduced_prob;
    double *xr = (double *)calloc(red->n > 0 ? red->n : 1, sizeof(double));
    double *yr = (double *)calloc(red->m > 0 ? red->m : 1, sizeof(double));
    double *zr = (double *)calloc(red->n > 0 ? red->n : 1, sizeof(double));

    postsolve(p, xr, yr, zr);

    char *err = kkt_assert(n, m, Ax, Ai, Ap, lhs, rhs, lbs, ubs, c,
                           Qx, Qi, Qp, Rx, Ri, Rp,
                           p->sol->x, p->sol->y, p->sol->z);

    free(xr); free(yr); free(zr);
    free_settings(stgs);
    free_presolver(p);
    return err;
}

/* ------------------------------------------------------------------
 * Test C: UNCHANGED presolve — presolver doesn't reduce anything, so
 * postsolve is just identity pass-through. Feed a known optimum.
 *
 *   min 0.5*(x0^2 + x1^2) + x0 + x1
 *   s.t. x0 + x1 >= -1, x_i free
 *
 * Unique optimum is x0=x1=-1 (interior, constraint inactive since sum=-2<=-1? wait lhs -1: -2 < -1, violates).
 * Let's use constraint x0+x1 >= -10 so constraint is inactive at (-1,-1).
 *
 * Analytical: min_x (0.5||x||^2 + 1^T x) without constraint gives x=(-1,-1).
 * Feed this as reduced solution. Postsolve should return same values with
 * z=0, y=0 (inactive constraint). Then verify KKT.
 * ------------------------------------------------------------------ */
static char *test_kkt_identity_passthrough()
{
    size_t m = 1;
    size_t n = 2;

    double Ax[] = {1.0, 1.0};
    int Ai[] = {0, 1};
    int Ap[] = {0, 2};

    double lhs[] = {-10.0};
    double rhs[] = {10.0};
    double lbs[] = {-INFINITY, -INFINITY};
    double ubs[] = {INFINITY, INFINITY};
    double c[] = {1.0, 1.0};

    double Qx[] = {1.0, 1.0};
    int Qi[] = {0, 1};
    int Qp[] = {0, 1, 2};
    size_t Qnnz = 2;

    Settings *stgs = default_settings();
    stgs->verbose = false;

    Presolver *p = new_qp_presolver_qr(Ax, Ai, Ap, m, n, 2,
                                       lhs, rhs, lbs, ubs, c,
                                       Qx, Qi, Qp, Qnnz,
                                       NULL, NULL, NULL, 0, 0, stgs);
    mu_assert("presolver creation failed", p != NULL);
    run_presolver(p);

    PresolvedProblem *red = p->reduced_prob;
    /* Provide the analytical optimum for the reduced problem.
     * If presolve didn't reduce anything, feed x=(-1,-1), y=0, z=0. */
    double *xr = (double *)calloc(red->n > 0 ? red->n : 1, sizeof(double));
    double *yr = (double *)calloc(red->m > 0 ? red->m : 1, sizeof(double));
    double *zr = (double *)calloc(red->n > 0 ? red->n : 1, sizeof(double));
    for (size_t i = 0; i < red->n; i++) xr[i] = -1.0;

    postsolve(p, xr, yr, zr);

    char *err = kkt_assert(n, m, Ax, Ai, Ap, lhs, rhs, lbs, ubs, c,
                           Qx, Qi, Qp, NULL, NULL, NULL,
                           p->sol->x, p->sol->y, p->sol->z);

    free(xr); free(yr); free(zr);
    free_settings(stgs);
    free_presolver(p);
    return err;
}

/* ------------------------------------------------------------------
 * Test D: QP with Q off-diagonal + R off-diagonal, sequential fixing.
 *
 *   min 0.5 x^T P x    where P = Q + R R^T
 *     Q = [[1, 0.5], [0.5, 2]]  (upper stored)
 *     R = [[1], [1]]     => R R^T = [[1,1],[1,1]]
 *     P = [[2, 1.5], [1.5, 3]]
 *   s.t. x0 + x1 = 2, x0 = 1 tight, x1 >= 0
 *
 * Fixing x0 = 1 forces x1 = 1 via singleton row.
 * Analytical KKT (primal interior for x1):
 *   stat at x1: 1.5 + 3 - y - z1 = 0, z1 = 0  ⇒  y = 4.5
 *   stat at x0: 2 + 1.5 - y - z0 = 0, z0 = -1.0
 *   (x0 is at tight bound, so z0 is unrestricted; sign is fine.)
 * ------------------------------------------------------------------ */
static char *test_kkt_Q_offdiag_plus_R_offdiag()
{
    size_t m = 1;
    size_t n = 2;

    double Ax[] = {1.0, 1.0};
    int Ai[] = {0, 1};
    int Ap[] = {0, 2};

    double lhs[] = {2.0};
    double rhs[] = {2.0};
    double lbs[] = {1.0, 0.0};
    double ubs[] = {1.0, INFINITY};
    double c[] = {0.0, 0.0};

    /* Q = [[1, 0.5], [0.5, 2]] upper triangular */
    double Qx[] = {1.0, 0.5, 2.0};
    int Qi[] = {0, 1, 1};
    int Qp[] = {0, 2, 3};
    size_t Qnnz = 3;

    double Rx[] = {1.0, 1.0};
    int Ri[] = {0, 0};
    int Rp[] = {0, 1, 2};
    size_t Rnnz = 2;
    size_t k = 1;

    Settings *stgs = default_settings();
    stgs->verbose = false;

    Presolver *p = new_qp_presolver_qr(Ax, Ai, Ap, m, n, 2,
                                       lhs, rhs, lbs, ubs, c,
                                       Qx, Qi, Qp, Qnnz,
                                       Rx, Ri, Rp, Rnnz, k, stgs);
    mu_assert("presolver creation failed", p != NULL);
    run_presolver(p);

    PresolvedProblem *red = p->reduced_prob;
    double *xr = (double *)calloc(red->n > 0 ? red->n : 1, sizeof(double));
    double *yr = (double *)calloc(red->m > 0 ? red->m : 1, sizeof(double));
    double *zr = (double *)calloc(red->n > 0 ? red->n : 1, sizeof(double));

    postsolve(p, xr, yr, zr);

    char *err = kkt_assert(n, m, Ax, Ai, Ap, lhs, rhs, lbs, ubs, c,
                           Qx, Qi, Qp, Rx, Ri, Rp,
                           p->sol->x, p->sol->y, p->sol->z);

    free(xr); free(yr); free(zr);
    free_settings(stgs);
    free_presolver(p);
    return err;
}

static int counter_qr_kkt = 0;

static const char *all_tests_qr_kkt()
{
    mu_run_test(test_kkt_identity_passthrough, counter_qr_kkt);
    mu_run_test(test_kkt_q_only_fix_plus_ston, counter_qr_kkt);
    mu_run_test(test_kkt_R_offdiag_fix, counter_qr_kkt);
    mu_run_test(test_kkt_Q_offdiag_plus_R_offdiag, counter_qr_kkt);
    return 0;
}

int test_qr_kkt(void)
{
    const char *result = all_tests_qr_kkt();
    if (result != 0)
    {
        printf("%s\n", result);
        printf("QR KKT: TEST FAILED!\n");
    }
    else
    {
        printf("QR KKT: ALL TESTS PASSED\n");
    }
    printf("QR KKT: Tests run: %d\n", counter_qr_kkt);
    return result == 0;
}

#endif /* TEST_QR_KKT_H */
