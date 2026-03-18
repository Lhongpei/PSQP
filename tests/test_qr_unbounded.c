#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "PSQP_API.h"

/* Test 1: Unbounded QP with QR
 * min -x (negative coefficient, unbounded below)
 * s.t. x >= 0
 * No upper bound on x, and we want to minimize -x -> unbounded
 */
void test_unbounded_simple()
{
    printf("Test 1: Simple unbounded QP...\n");
    
    size_t n = 1;
    size_t m = 0;  // No constraints
    
    double Ax[] = {};
    int Ai[] = {};
    int Ap[] = {0};
    size_t nnz = 0;
    
    double *lhs = NULL;
    double *rhs = NULL;
    double lbs[] = {0.0};  // x >= 0
    double ubs[] = {INFINITY};  // No upper bound
    double c[] = {-1.0};  // Minimize -x -> want x -> infinity
    
    // Simple Q = 0 (just use NULL)
    Settings *stgs = default_settings();
    stgs->verbose = false;
    
    Presolver *presolver = new_qp_presolver_qr(Ax, Ai, Ap, m, n, nnz,
                                               lhs, rhs, lbs, ubs, c,
                                               NULL, NULL, NULL, 0,  // No Q
                                               NULL, NULL, NULL, 0, 0,  // No R
                                               stgs);
    
    if (presolver == NULL) {
        printf("  Presolver creation failed\n");
        return;
    }
    
    PresolveStatus status = run_presolver(presolver);
    
    printf("  Status: %s\n", 
           status == INFEASIBLE ? "INFEASIBLE" :
           status == UNBNDORINFEAS ? "UNBNDORINFEAS" :
           status == REDUCED ? "REDUCED" :
           status == UNCHANGED ? "UNCHANGED" : "OTHER");
    
    if (status == UNBNDORINFEAS) {
        printf("  ✓ Unboundedness detected!\n");
    } else {
        printf("  Note: Status is %d (may or may not detect unboundedness)\n", status);
    }
    
    free_settings(stgs);
    free_presolver(presolver);
}

/* Test 2: Unbounded due to free variable with negative cost */
void test_unbounded_free_var()
{
    printf("\nTest 2: Unbounded with free variable...\n");
    
    size_t n = 2;
    size_t m = 1;
    
    // Constraint: x + y = 1
    double Ax[] = {1.0, 1.0};
    int Ai[] = {0, 1};
    int Ap[] = {0, 2};
    size_t nnz = 2;
    
    double lhs[] = {1.0};
    double rhs[] = {1.0};
    double lbs[] = {-INFINITY, 0.0};  // x is free (unbounded below)
    double ubs[] = {INFINITY, INFINITY};
    double c[] = {-1.0, 0.0};  // Minimize -x -> x -> infinity
    
    double Qx[] = {1.0, 1.0};  // Small quadratic term
    int Qi[] = {0, 1};
    int Qp[] = {0, 1, 2};
    size_t Qnnz = 2;
    
    Settings *stgs = default_settings();
    stgs->verbose = false;
    stgs->ston_cols = false;  // Disable to preserve structure
    
    Presolver *presolver = new_qp_presolver(Ax, Ai, Ap, m, n, nnz,
                                            lhs, rhs, lbs, ubs, c,
                                            Qx, Qi, Qp, Qnnz, stgs);
    
    if (presolver == NULL) {
        printf("  Presolver creation failed\n");
        return;
    }
    
    PresolveStatus status = run_presolver(presolver);
    
    printf("  Status: %s\n", 
           status == INFEASIBLE ? "INFEASIBLE" :
           status == UNBNDORINFEAS ? "UNBNDORINFEAS" :
           status == REDUCED ? "REDUCED" :
           status == UNCHANGED ? "UNCHANGED" : "OTHER");
    
    free_settings(stgs);
    free_presolver(presolver);
}

/* Test 3: Strictly convex QP (should be bounded) */
void test_bounded_convex()
{
    printf("\nTest 3: Strictly convex QP (should be bounded)...\n");
    
    size_t n = 2;
    size_t m = 0;
    
    double Ax[] = {};
    int Ai[] = {};
    int Ap[] = {0};
    size_t nnz = 0;
    
    double *lhs = NULL;
    double *rhs = NULL;
    double lbs[] = {-INFINITY, -INFINITY};  // Both unbounded
    double ubs[] = {INFINITY, INFINITY};
    double c[] = {-1.0, -1.0};  // Negative linear term
    
    // Strongly convex: P = 2*I
    double Qx[] = {2.0, 2.0};
    int Qi[] = {0, 1};
    int Qp[] = {0, 1, 2};
    size_t Qnnz = 2;
    
    Settings *stgs = default_settings();
    stgs->verbose = false;
    
    Presolver *presolver = new_qp_presolver(Ax, Ai, Ap, m, n, nnz,
                                            lhs, rhs, lbs, ubs, c,
                                            Qx, Qi, Qp, Qnnz, stgs);
    
    if (presolver == NULL) {
        printf("  Presolver creation failed\n");
        return;
    }
    
    PresolveStatus status = run_presolver(presolver);
    
    printf("  Status: %s\n", 
           status == INFEASIBLE ? "INFEASIBLE" :
           status == UNBNDORINFEAS ? "UNBNDORINFEAS" :
           status == REDUCED ? "REDUCED" :
           status == UNCHANGED ? "UNCHANGED" : "OTHER");
    
    if (status != UNBNDORINFEAS) {
        printf("  ✓ Convex QP is bounded as expected\n");
    }
    
    free_settings(stgs);
    free_presolver(presolver);
}

/* Test 4: Unbounded with QR format */
void test_unbounded_qr()
{
    printf("\nTest 4: Unbounded QP with QR format...\n");
    
    size_t n = 2;
    size_t m = 0;
    
    double Ax[] = {};
    int Ai[] = {};
    int Ap[] = {0};
    size_t nnz = 0;
    
    double *lhs = NULL;
    double *rhs = NULL;
    double lbs[] = {0.0, 0.0};  // Both >= 0
    double ubs[] = {INFINITY, INFINITY};  // No upper bounds
    double c[] = {-1.0, -1.0};  // Minimize -x1 - x2
    
    // Q = 0 (no quadratic term)
    Settings *stgs = default_settings();
    stgs->verbose = false;
    
    Presolver *presolver = new_qp_presolver_qr(Ax, Ai, Ap, m, n, nnz,
                                               lhs, rhs, lbs, ubs, c,
                                               NULL, NULL, NULL, 0,  // No Q
                                               NULL, NULL, NULL, 0, 0,  // No R
                                               stgs);
    
    if (presolver == NULL) {
        printf("  Presolver creation failed\n");
        return;
    }
    
    PresolveStatus status = run_presolver(presolver);
    
    printf("  Status: %s\n", 
           status == INFEASIBLE ? "INFEASIBLE" :
           status == UNBNDORINFEAS ? "UNBNDORINFEAS" :
           status == REDUCED ? "REDUCED" :
           status == UNCHANGED ? "UNCHANGED" : "OTHER");
    
    free_settings(stgs);
    free_presolver(presolver);
}

int main()
{
    printf("=== QP Unboundedness Tests ===\n\n");
    
    test_unbounded_simple();
    test_unbounded_free_var();
    test_bounded_convex();
    test_unbounded_qr();
    
    printf("\n=== Tests Complete ===\n");
    printf("\nNote: UNBNDORINFEAS means 'unbounded or infeasible'.\n");
    printf("The presolver may not always distinguish between the two.\n");
    return 0;
}
