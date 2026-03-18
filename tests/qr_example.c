/*
 * QR decomposition example: P = Q + R*R^T
 * 
 * This example demonstrates using PSQP with P decomposed as Q + R*R^T.
 * This is useful when P has a low-rank component.
 */

#include "PSQP_API.h"
#include "PSQP_stats.h"
#include "PSQP_sol.h"
#include <stdio.h>
#include <math.h>

void example_factor_model()
{
    printf("\n=== Factor Model Portfolio Optimization ===\n");
    printf("P = D + F*F^T where D is diagonal (idiosyncratic risk)\n");
    printf("               and F is factor loadings (systematic risk)\n\n");
    
    size_t n = 5;  /* 5 assets */
    size_t m = 1;  /* 1 constraint (budget) */
    
    /* Constraint: sum(x) = 1 */
    double Ax[] = {1.0, 1.0, 1.0, 1.0, 1.0};
    int Ai[] = {0, 1, 2, 3, 4};
    int Ap[] = {0, 5};
    size_t nnz = 5;
    
    double lhs[] = {1.0};
    double rhs[] = {1.0};
    double lbs[] = {0.0, 0.0, 0.0, 0.0, 0.0};
    double ubs[] = {INFINITY, INFINITY, INFINITY, INFINITY, INFINITY};
    double c[] = {-0.1, -0.12, -0.08, -0.15, -0.11};  /* Expected returns */
    
    /* D matrix: diagonal idiosyncratic variances */
    /* P_ii = D_ii + sum_k F_ik^2 */
    double Dx[] = {0.05, 0.04, 0.06, 0.03, 0.05};
    int Di[] = {0, 1, 2, 3, 4};
    int Dp[] = {0, 1, 2, 3, 4, 5};
    size_t Dnnz = 5;
    
    /* F matrix: 5 assets x 2 factors */
    /* Factor 1: market factor, Factor 2: size factor */
    double Fx[] = {0.3, 0.4, 0.2, 0.5, 0.35,   /* Factor 1 loadings */
                   0.1, -0.1, 0.2, 0.0, 0.15}; /* Factor 2 loadings */
    int Fi[] = {0, 1, 2, 3, 4, 0, 1, 2, 3, 4};
    int Fp[] = {0, 5, 10};
    size_t Fnnz = 10;
    size_t k = 2;  /* 2 factors */
    
    Settings *stgs = default_settings();
    stgs->verbose = true;
    
    /* Use QR decomposition API */
    Presolver *presolver = new_qp_presolver_qr(Ax, Ai, Ap, m, n, nnz,
                                               lhs, rhs, lbs, ubs, c,
                                               Dx, Di, Dp, Dnnz,
                                               Fx, Fi, Fp, Fnnz, k, stgs);
    
    if (!presolver) {
        printf("ERROR: Failed to create presolver\n");
        return;
    }
    
    PresolveStatus status = run_presolver(presolver);
    PresolvedProblem *reduced = presolver->reduced_prob;
    
    printf("\nPresolve status: %s\n", 
           status == REDUCED ? "REDUCED" : 
           status == UNCHANGED ? "UNCHANGED" : "OTHER");
    
    printf("Original: %zu vars, %zu constraints\n", n, m);
    printf("Reduced:  %zu vars, %zu constraints\n", reduced->n, reduced->m);
    
    if (reduced->has_quadratic) {
        printf("Quadratic term: %zu non-zeros in P\n", reduced->Pnnz);
    }
    
    /* Output Q and R separately */
    if (reduced->has_quad_qr) {
        printf("\n--- Q matrix (sparse symmetric) ---\n");
        printf("Qnnz = %zu\n", reduced->Qnnz);
        if (reduced->Qnnz > 0 && reduced->Qx) {
            printf("Q values: ");
            for (size_t i = 0; i < reduced->Qnnz; i++) {
                printf("%.2f ", reduced->Qx[i]);
            }
            printf("\n");
        }
        
        printf("\n--- R matrix (n x k) ---\n");
        printf("Rnnz = %zu, k = %zu\n", reduced->Rnnz, reduced->k);
        if (reduced->Rnnz > 0 && reduced->Rx) {
            printf("R values: ");
            for (size_t i = 0; i < reduced->Rnnz; i++) {
                printf("%.2f ", reduced->Rx[i]);
            }
            printf("\n");
        }
    }
    
    free_settings(stgs);
    free_presolver(presolver);
}

void example_low_rank_plus_diagonal()
{
    printf("\n=== Low-Rank Plus Diagonal Structure ===\n");
    printf("P = diag(d) + u*u^T where u is a sparse vector\n\n");
    
    size_t n = 4;
    size_t m = 1;
    
    /* Constraint: sum(x) = 1 */
    double Ax[] = {1.0, 1.0, 1.0, 1.0};
    int Ai[] = {0, 1, 2, 3};
    int Ap[] = {0, 4};
    size_t nnz = 4;
    
    double lhs[] = {1.0};
    double rhs[] = {1.0};
    double lbs[] = {0.0, 0.0, 0.0, 0.0};
    double ubs[] = {INFINITY, INFINITY, INFINITY, INFINITY};
    double c[] = {1.0, 1.0, 1.0, 1.0};
    
    /* Diagonal Q */
    double Qx[] = {2.0, 2.0, 2.0, 2.0};
    int Qi[] = {0, 1, 2, 3};
    int Qp[] = {0, 1, 2, 3, 4};
    size_t Qnnz = 4;
    
    /* Rank-1 component: u = [1, 1, 0, 0]^T
     * u*u^T adds [1, 1, 0, 0; 1, 1, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0]
     * So P = [3, 1, 0, 0; 1, 3, 0, 0; 0, 0, 2, 0; 0, 0, 0, 2]
     */
    double Rx[] = {1.0, 1.0};  /* Only first two elements non-zero */
    int Ri[] = {0, 1};
    int Rp[] = {0, 2, 2, 2, 2};  /* Each row of R has specific non-zeros */
    size_t Rnnz = 2;
    size_t k = 1;  /* Rank-1 */
    
    /* Note: Actually Rp should be [0,1,2,2,2] if R is 4x1
     * Let me fix this: R[0][0]=1, R[1][0]=1, R[2][0]=0, R[3][0]=0 */
    int Rp_correct[] = {0, 1, 2, 2, 2};
    
    Settings *stgs = default_settings();
    stgs->verbose = false;
    
    Presolver *presolver = new_qp_presolver_qr(Ax, Ai, Ap, m, n, nnz,
                                               lhs, rhs, lbs, ubs, c,
                                               Qx, Qi, Qp, Qnnz,
                                               Rx, Ri, Rp_correct, Rnnz, k, stgs);
    
    if (!presolver) {
        printf("ERROR: Failed to create presolver\n");
        return;
    }
    
    PresolveStatus status = run_presolver(presolver);
    PresolvedProblem *reduced = presolver->reduced_prob;
    
    printf("Status: %s\n", status == REDUCED ? "REDUCED" : "UNCHANGED");
    printf("P = Q + u*u^T where u = [1,1,0,0]^T\n");
    printf("Expected P = [3,1,0,0; 1,3,0,0; 0,0,2,0; 0,0,0,2]\n");
    
    /* Output Q and R separately */
    if (reduced->has_quad_qr) {
        printf("\n--- Q matrix after presolve ---\n");
        printf("Qnnz = %zu\n", reduced->Qnnz);
        if (reduced->Qnnz > 0 && reduced->Qx) {
            printf("Q values: ");
            for (size_t i = 0; i < reduced->Qnnz; i++) {
                printf("%.2f ", reduced->Qx[i]);
            }
            printf("\n");
        }
        
        printf("\n--- R matrix after presolve ---\n");
        printf("Rnnz = %zu, k = %zu\n", reduced->Rnnz, reduced->k);
        if (reduced->Rnnz > 0 && reduced->Rx) {
            printf("R values: ");
            for (size_t i = 0; i < reduced->Rnnz; i++) {
                printf("%.2f ", reduced->Rx[i]);
            }
            printf("\n");
        }
    }
    
    free_settings(stgs);
    free_presolver(presolver);
}

void example_qr_shrink_demo()
{
    printf("\n=== QR Shrink Demo: Variables Removed ===\n");
    printf("Demonstrates Q and R matrices being shrinked when variables are removed\n\n");
    
    size_t n = 5;
    size_t m = 1;
    
    /* Simple constraint: x0 + x1 + x2 + x3 + x4 >= 0 */
    double Ax[] = {1.0, 1.0, 1.0, 1.0, 1.0};
    int Ai[] = {0, 1, 2, 3, 4};
    int Ap[] = {0, 5};
    size_t nnz = 5;
    
    double lhs[] = {0.0};
    double rhs[] = {INFINITY};
    /* Fix x2 = 1.0 (tight bounds), leave others free */
    double lbs[] = {0.0, 0.0, 1.0, 0.0, 0.0};
    double ubs[] = {INFINITY, INFINITY, 1.0, INFINITY, INFINITY};
    double c[] = {1.0, 1.0, 1.0, 1.0, 1.0};
    
    /* Q matrix: diagonal with some off-diagonal terms */
    /* Q = diag(2, 2, 3, 3, 2) with Q[0,1] = Q[1,0] = 0.5 */
    double Qx[] = {2.0, 0.5, 2.0, 3.0, 3.0, 2.0};  /* Upper triangular: Q00, Q01, Q11, Q22, Q33, Q44 */
    int Qi[] = {0, 1, 1, 2, 3, 4};
    int Qp[] = {0, 1, 3, 4, 5, 6};
    size_t Qnnz = 6;
    
    /* R matrix: 5x2 low-rank component */
    /* R = [1 0; 1 0; 0 1; 0 1; 1 0] - first two vars share factor 1, next two share factor 2 */
    double Rx[] = {1.0,    /* Row 0: col 0 */
                   1.0,    /* Row 1: col 0 */
                   1.0,    /* Row 2: col 1 */
                   1.0,    /* Row 3: col 1 */
                   1.0};   /* Row 4: col 0 */
    int Ri[] = {0, 0, 1, 1, 0};
    int Rp[] = {0, 1, 2, 3, 4, 5};  /* Row pointers: 1 nnz per row */
    size_t Rnnz = 5;
    size_t k = 2;
    
    Settings *stgs = default_settings();
    stgs->verbose = true;
    stgs->ston_cols = false;  /* Disable singleton column removal to keep more variables */
    stgs->dton_eq = false;    /* Disable doubleton equality removal */
    
    Presolver *presolver = new_qp_presolver_qr(Ax, Ai, Ap, m, n, nnz,
                                               lhs, rhs, lbs, ubs, c,
                                               Qx, Qi, Qp, Qnnz,
                                               Rx, Ri, Rp, Rnnz, k, stgs);
    
    if (!presolver) {
        printf("ERROR: Failed to create presolver\n");
        return;
    }
    
    PresolveStatus status = run_presolver(presolver);
    PresolvedProblem *reduced = presolver->reduced_prob;
    
    printf("\n--- Results ---\n");
    printf("Original: %zu vars, %zu constraints\n", n, m);
    printf("Reduced:  %zu vars, %zu constraints\n", reduced->n, reduced->m);
    printf("Status: %s\n", status == REDUCED ? "REDUCED" : "UNCHANGED");
    
    /* Output Q and R separately */
    if (reduced->has_quad_qr) {
        printf("\n--- Q matrix after presolve (shrinked) ---\n");
        printf("Qnnz = %zu (original: %zu)\n", reduced->Qnnz, Qnnz);
        if (reduced->Qnnz > 0 && reduced->Qx) {
            printf("Q values: ");
            for (size_t i = 0; i < reduced->Qnnz; i++) {
                printf("%.2f ", reduced->Qx[i]);
            }
            printf("\nQ row pointers (Qp): ");
            for (size_t i = 0; i <= reduced->n; i++) {
                printf("%d ", reduced->Qp[i]);
            }
            printf("\n");
        }
        
        printf("\n--- R matrix after presolve (shrinked) ---\n");
        printf("Rnnz = %zu (original: %zu), k = %zu\n", reduced->Rnnz, Rnnz, reduced->k);
        if (reduced->Rnnz > 0 && reduced->Rx) {
            printf("R values: ");
            for (size_t i = 0; i < reduced->Rnnz; i++) {
                printf("%.2f ", reduced->Rx[i]);
            }
            printf("\nR row pointers (Rp): ");
            for (size_t i = 0; i <= reduced->n; i++) {
                printf("%d ", reduced->Rp[i]);
            }
            printf("\n");
        }
        
        /* Verify dimensions match */
        if (reduced->n > 0) {
            printf("\n--- Verification ---\n");
            printf("Q matrix row count (from Qp): %d (should match n=%zu)\n", 
                   reduced->Qp ? reduced->Qp[reduced->n] : -1, reduced->n);
            printf("R matrix row count (from Rp): %d (should match n=%zu)\n", 
                   reduced->Rp ? reduced->Rp[reduced->n] : -1, reduced->n);
        }
    } else {
        printf("\nNo QR data in reduced problem (all variables removed or no quadratic term)\n");
    }
    
    free_settings(stgs);
    free_presolver(presolver);
}

int main()
{
    printf("========================================\n");
    printf("   PSQP QR Decomposition Examples\n");
    printf("   P = Q + R*R^T\n");
    printf("========================================\n");
    
    example_factor_model();
    example_low_rank_plus_diagonal();
    example_qr_shrink_demo();
    
    printf("\n========================================\n");
    printf("   Examples completed!\n");
    printf("========================================\n");
    
    return 0;
}
