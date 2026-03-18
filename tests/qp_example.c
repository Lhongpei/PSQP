/*
 * QP Presolver Usage Example
 * 
 * This file demonstrates how to use PSLP's QP presolver API.
 * 
 * Problem format:
 *   minimize    (1/2) * x^T * P * x + c^T * x
 *   subject to  lhs <= A * x <= rhs
 *               lbs <= x <= ubs
 * 
 * P is stored in CSR format (upper triangular part only since P is symmetric)
 * A is stored in CSR format
 */

#include "PSQP_API.h"
#include "PSQP_stats.h"
#include "PSQP_sol.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Example 1: Simple portfolio optimization-like problem
 * 
 * minimize    (1/2) * x^T * P * x - c^T * x
 * subject to  sum(x) = 1
 *             x >= 0
 * 
 * where P represents the covariance matrix and c represents expected returns.
 */
void example_portfolio_qp()
{
    printf("\n=== Example 1: Portfolio-like QP ===\n");
    
    size_t n = 4;  /* 4 assets */
    size_t m = 1;  /* 1 constraint (budget) */
    
    /* Constraint: sum(x) = 1 */
    double Ax[] = {1.0, 1.0, 1.0, 1.0};
    int Ai[] = {0, 1, 2, 3};
    int Ap[] = {0, 4};
    size_t nnz = 4;
    
    double lhs[] = {1.0};
    double rhs[] = {1.0};
    double lbs[] = {0.0, 0.0, 0.0, 0.0};
    double ubs[] = {INFINITY, INFINITY, INFINITY, INFINITY};
    /* Maximize returns: minimize negative returns */
    double c[] = {-0.1, -0.15, -0.08, -0.12};  /* negative because we minimize */
    
    /* Covariance matrix P (symmetric, positive semidefinite)
     * P = [0.1  0.02 0.01 0.03
     *      0.02 0.15 0.02 0.04
     *      0.01 0.02 0.08 0.02
     *      0.03 0.04 0.02 0.12]
     * Upper triangular stored as:
     */
    double Px[] = {0.1, 0.02, 0.01, 0.03,   /* row 0 */
                   0.15, 0.02, 0.04,        /* row 1 */
                   0.08, 0.02,              /* row 2 */
                   0.12};                   /* row 3 */
    int Pi[] = {0, 1, 2, 3, 1, 2, 3, 2, 3, 3};
    int Pp[] = {0, 4, 7, 9, 10};
    size_t Pnnz = 10;
    
    /* Create settings */
    Settings *stgs = default_settings();
    stgs->verbose = true;  /* Print presolving info */
    
    /* Create QP presolver */
    Presolver *presolver = new_qp_presolver(Ax, Ai, Ap, m, n, nnz,
                                            lhs, rhs, lbs, ubs, c,
                                            Px, Pi, Pp, Pnnz, stgs);
    if (!presolver) {
        printf("ERROR: Failed to create presolver\n");
        return;
    }
    
    /* Run presolver */
    PresolveStatus status = run_presolver(presolver);
    
    printf("\nPresolve status: %s\n", 
           status == REDUCED ? "REDUCED" : 
           status == UNCHANGED ? "UNCHANGED" : 
           status == INFEASIBLE ? "INFEASIBLE" : 
           status == UNBNDORINFEAS ? "UNBOUNDED_OR_INFEASIBLE" : "UNKNOWN");
    
    /* Access reduced problem */
    PresolvedProblem *reduced = presolver->reduced_prob;
    printf("\nReduced problem size: %zu x %zu, nnz=%zu\n", 
           reduced->m, reduced->n, reduced->nnz);
    printf("Objective offset: %f\n", reduced->obj_offset);
    
    if (reduced->has_quadratic) {
        printf("Reduced P matrix: %zu non-zeros\n", reduced->Pnnz);
    }
    
    /* Cleanup */
    free_settings(stgs);
    free_presolver(presolver);
}

/* Example 2: QP with variable bounds causing fixes
 * Demonstrates how fixed variables affect the objective offset
 */
void example_fixed_variables()
{
    printf("\n=== Example 2: QP with Fixed Variables ===\n");
    
    size_t n = 3;
    size_t m = 1;
    
    /* Constraint: x1 + x2 + x3 = 6 */
    double Ax[] = {1.0, 1.0, 1.0};
    int Ai[] = {0, 1, 2};
    int Ap[] = {0, 3};
    size_t nnz = 3;
    
    double lhs[] = {6.0};
    double rhs[] = {6.0};
    /* x1 is fixed to 1 by bounds */
    double lbs[] = {1.0, 0.0, 0.0};
    double ubs[] = {1.0, INFINITY, INFINITY};
    double c[] = {0.0, 0.0, 0.0};
    
    /* P = 2*I (identity scaled by 2) */
    double Px[] = {2.0, 2.0, 2.0};
    int Pi[] = {0, 1, 2};
    int Pp[] = {0, 1, 2, 3};
    size_t Pnnz = 3;
    
    Settings *stgs = default_settings();
    stgs->verbose = true;
    
    Presolver *presolver = new_qp_presolver(Ax, Ai, Ap, m, n, nnz,
                                            lhs, rhs, lbs, ubs, c,
                                            Px, Pi, Pp, Pnnz, stgs);
    if (!presolver) {
        printf("ERROR: Failed to create presolver\n");
        return;
    }
    
    PresolveStatus status = run_presolver(presolver);
    PresolvedProblem *reduced = presolver->reduced_prob;
    
    printf("\nOriginal size: %zu x %zu\n", m, n);
    printf("Reduced size: %zu x %zu\n", reduced->m, reduced->n);
    
    /* The offset should include contribution from fixed x1:
     * 0.5 * P_00 * x1^2 = 0.5 * 2 * 1^2 = 1
     */
    printf("Objective offset: %f (expected contribution from x1: 1.0)\n", 
           reduced->obj_offset);
    
    free_settings(stgs);
    free_presolver(presolver);
}

/* Example 3: Postsolve - recovering original solution
 * Shows how to use postsolve to recover the solution to the original problem
 */
void example_postsolve()
{
    printf("\n=== Example 3: Postsolve Example ===\n");
    
    size_t n = 4;
    size_t m = 2;
    
    /* Constraints:
     * x1 + x2 + x3 + x4 = 10
     * x1 - x2 = 0  (implies x1 = x2)
     */
    double Ax[] = {1.0, 1.0, 1.0, 1.0, 1.0, -1.0};
    int Ai[] = {0, 1, 2, 3, 0, 1};
    int Ap[] = {0, 4, 6};
    size_t nnz = 6;
    
    double lhs[] = {10.0, 0.0};
    double rhs[] = {10.0, 0.0};
    double lbs[] = {0.0, 0.0, 0.0, 0.0};
    double ubs[] = {INFINITY, INFINITY, INFINITY, INFINITY};
    double c[] = {1.0, 2.0, 3.0, 4.0};
    
    /* Simple diagonal P */
    double Px[] = {2.0, 2.0, 2.0, 2.0};
    int Pi[] = {0, 1, 2, 3};
    int Pp[] = {0, 1, 2, 3, 4};
    size_t Pnnz = 4;
    
    Settings *stgs = default_settings();
    stgs->verbose = true;
    
    Presolver *presolver = new_qp_presolver(Ax, Ai, Ap, m, n, nnz,
                                            lhs, rhs, lbs, ubs, c,
                                            Px, Pi, Pp, Pnnz, stgs);
    if (!presolver) {
        printf("ERROR: Failed to create presolver\n");
        return;
    }
    
    PresolveStatus status = run_presolver(presolver);
    PresolvedProblem *reduced = presolver->reduced_prob;
    
    printf("\nReduced problem: %zu x %zu\n", reduced->m, reduced->n);
    
    /* Simulate solving the reduced problem (dummy solution) */
    if (reduced->n > 0 && reduced->m > 0) {
        /* In practice, you would call your QP solver here
         * For this example, we use dummy values
         */
        double *x_reduced = (double *)calloc(reduced->n, sizeof(double));
        double *y_reduced = (double *)calloc(reduced->m, sizeof(double));
        double *z_reduced = (double *)calloc(reduced->n, sizeof(double));
        
        /* Fill with dummy solution values */
        for (size_t i = 0; i < reduced->n; i++) {
            x_reduced[i] = 2.5;  /* Dummy value */
        }
        for (size_t i = 0; i < reduced->m; i++) {
            y_reduced[i] = 1.0;  /* Dummy multiplier */
        }
        
        /* Postsolve to recover original solution */
        postsolve(presolver, x_reduced, y_reduced, z_reduced);
        
        /* Access recovered solution */
        printf("\nRecovered solution to original problem:\n");
        printf("x = [");
        for (size_t i = 0; i < n; i++) {
            printf("%.4f%s", presolver->sol->x[i], i < n-1 ? ", " : "");
        }
        printf("]\n");
        
        printf("y = [");
        for (size_t i = 0; i < m; i++) {
            printf("%.4f%s", presolver->sol->y[i], i < m-1 ? ", " : "");
        }
        printf("]\n");
        
        free(x_reduced);
        free(y_reduced);
        free(z_reduced);
    } else {
        printf("Problem fully reduced - no postsolve needed\n");
    }
    
    free_settings(stgs);
    free_presolver(presolver);
}

/* Example 4: Using LP through QP interface
 * Demonstrates that you can use new_qp_presolver for LP problems too
 */
void example_lp_through_qp()
{
    printf("\n=== Example 4: LP through QP Interface ===\n");
    
    size_t n = 3;
    size_t m = 2;
    
    /* Constraints:
     * x1 + x2 + x3 = 6
     * x1 - x3 >= 1
     */
    double Ax[] = {1.0, 1.0, 1.0, 1.0, -1.0};
    int Ai[] = {0, 1, 2, 0, 2};
    int Ap[] = {0, 3, 5};
    size_t nnz = 5;
    
    double lhs[] = {6.0, 1.0};
    double rhs[] = {6.0, INFINITY};
    double lbs[] = {0.0, 0.0, 0.0};
    double ubs[] = {INFINITY, INFINITY, INFINITY};
    double c[] = {1.0, 2.0, 3.0};  /* minimize c^T * x */
    
    /* No P matrix - this is an LP */
    Settings *stgs = default_settings();
    stgs->verbose = true;
    
    /* Pass NULL for P matrix arguments */
    Presolver *presolver = new_qp_presolver(Ax, Ai, Ap, m, n, nnz,
                                            lhs, rhs, lbs, ubs, c,
                                            NULL, NULL, NULL, 0, stgs);
    if (!presolver) {
        printf("ERROR: Failed to create presolver\n");
        return;
    }
    
    PresolveStatus status = run_presolver(presolver);
    PresolvedProblem *reduced = presolver->reduced_prob;
    
    printf("\nReduced problem: %zu x %zu\n", reduced->m, reduced->n);
    printf("Has quadratic term: %s\n", reduced->has_quadratic ? "yes" : "no");
    
    free_settings(stgs);
    free_presolver(presolver);
}

int main()
{
    printf("========================================\n");
    printf("   PSLP QP Presolver Examples\n");
    printf("========================================\n");
    
    example_portfolio_qp();
    example_fixed_variables();
    example_postsolve();
    example_lp_through_qp();
    
    printf("\n========================================\n");
    printf("   All examples completed!\n");
    printf("========================================\n");
    
    return 0;
}
