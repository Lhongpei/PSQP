#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "PSQP_API.h"

/* Test 1: LP-like unbounded (no quadratic term) */
void test_unbounded_lp()
{
    printf("Test 1: LP-like unbounded problem...\n");
    
    size_t n = 2;
    size_t m = 1;
    
    // Constraint: x + y = 1
    double Ax[] = {1.0, 1.0};
    int Ai[] = {0, 1};
    int Ap[] = {0, 2};
    size_t nnz = 2;
    
    double lhs[] = {1.0};
    double rhs[] = {1.0};
    double lbs[] = {0.0, -INFINITY};  // y is unbounded below
    double ubs[] = {INFINITY, INFINITY};
    double c[] = {0.0, -1.0};  // Minimize -y -> y -> infinity
    
    // No P matrix (LP)
    Settings *stgs = default_settings();
    stgs->verbose = false;
    
    Presolver *presolver = new_presolver(Ax, Ai, Ap, m, n, nnz,
                                         lhs, rhs, lbs, ubs, c, stgs);
    
    if (presolver == NULL) {
        printf("  Presolver creation failed\n");
        return;
    }
    
    PresolveStatus status = run_presolver(presolver);
    
    printf("  Status: %s\n", 
           status == UNBNDORINFEAS ? "UNBNDORINFEAS" :
           status == INFEASIBLE ? "INFEASIBLE" :
           status == REDUCED ? "REDUCED" : "UNCHANGED");
    
    if (status == UNBNDORINFEAS) {
        printf("  ✓ Detected as unbounded/infeasible\n");
    }
    
    free_settings(stgs);
    free_presolver(presolver);
}

/* Test 2: QP with weak convexity - can still be unbounded */
void test_unbounded_qp_weak()
{
    printf("\nTest 2: QP with weak convexity...\n");
    
    size_t n = 2;
    size_t m = 1;
    
    // Constraint: x + y = 1
    double Ax[] = {1.0, 1.0};
    int Ai[] = {0, 1};
    int Ap[] = {0, 2};
    size_t nnz = 2;
    
    double lhs[] = {1.0};
    double rhs[] = {1.0};
    double lbs[] = {0.0, -INFINITY};  // y is free
    double ubs[] = {INFINITY, INFINITY};
    double c[] = {0.0, -1.0};  // Want y -> infinity
    
    // P = diag(0, 0) - no quadratic cost on y (semidefinite)
    // Actually use very small Q
    double Qx[] = {1e-8};  // Tiny Q for x only
    int Qi[] = {0};
    int Qp[] = {0, 1, 1};
    size_t Qnnz = 1;
    
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
           status == UNBNDORINFEAS ? "UNBNDORINFEAS" :
           status == INFEASIBLE ? "INFEASIBLE" :
           status == REDUCED ? "REDUCED" : "UNCHANGED");
    
    free_settings(stgs);
    free_presolver(presolver);
}

/* Test 3: QR format with strong convexity - should be bounded */
void test_bounded_qr()
{
    printf("\nTest 3: QR format with strong convexity (bounded)...\n");
    
    size_t n = 2;
    size_t m = 1;
    
    double Ax[] = {1.0, 1.0};
    int Ai[] = {0, 1};
    int Ap[] = {0, 2};
    size_t nnz = 2;
    
    double lhs[] = {1.0};
    double rhs[] = {1.0};
    double lbs[] = {-INFINITY, -INFINITY};  // Both free
    double ubs[] = {INFINITY, INFINITY};
    double c[] = {-1.0, -1.0};  // Negative cost
    
    // Strong convexity via QR: Q = I, R = small
    double Qx[] = {1.0, 1.0};  // Q = I
    int Qi[] = {0, 1};
    int Qp[] = {0, 1, 2};
    size_t Qnnz = 2;
    
    double Rx[] = {0.1, 0.1};  // Small R
    int Ri[] = {0, 0};
    int Rp[] = {0, 1, 2};
    size_t Rnnz = 2;
    size_t k = 1;
    
    Settings *stgs = default_settings();
    stgs->verbose = false;
    
    Presolver *presolver = new_qp_presolver_qr(Ax, Ai, Ap, m, n, nnz,
                                               lhs, rhs, lbs, ubs, c,
                                               Qx, Qi, Qp, Qnnz,
                                               Rx, Ri, Rp, Rnnz, k, stgs);
    
    if (presolver == NULL) {
        printf("  Presolver creation failed\n");
        return;
    }
    
    PresolveStatus status = run_presolver(presolver);
    
    printf("  Status: %s\n", 
           status == UNBNDORINFEAS ? "UNBNDORINFEAS" :
           status == INFEASIBLE ? "INFEASIBLE" :
           status == REDUCED ? "REDUCED" : "UNCHANGED");
    
    if (status != UNBNDORINFEAS) {
        printf("  ✓ Bounded QP detected correctly\n");
    }
    
    free_settings(stgs);
    free_presolver(presolver);
}

/* Test 4: Dual infeasibility detection */
void test_dual_infeasible()
{
    printf("\nTest 4: Testing dual infeasibility scenario...\n");
    
    size_t n = 1;
    size_t m = 0;
    
    // No constraints
    double Ax[] = {};
    int Ai[] = {};
    int Ap[] = {0};
    size_t nnz = 0;
    
    double *lhs = NULL;
    double *rhs = NULL;
    double lbs[] = {-INFINITY};  // Free variable
    double ubs[] = {INFINITY};
    double c[] = {1.0};  // Minimize x where x is free -> unbounded below
    
    Settings *stgs = default_settings();
    stgs->verbose = false;
    
    // Use presolver directly (which handles m=0)
    Presolver *presolver = new_presolver(Ax, Ai, Ap, 0, n, nnz,
                                         lhs, rhs, lbs, ubs, c, stgs);
    
    if (presolver == NULL) {
        printf("  Presolver creation failed (expected for m=0)\n");
        return;
    }
    
    PresolveStatus status = run_presolver(presolver);
    
    printf("  Status: %s\n", 
           status == UNBNDORINFEAS ? "UNBNDORINFEAS" :
           status == INFEASIBLE ? "INFEASIBLE" :
           status == REDUCED ? "REDUCED" : "UNCHANGED");
    
    free_settings(stgs);
    free_presolver(presolver);
}

int main()
{
    printf("=== QP Unboundedness Detection Tests ===\n\n");
    printf("Testing if presolver detects unbounded problems...\n\n");
    
    test_unbounded_lp();
    test_unbounded_qp_weak();
    test_bounded_qr();
    test_dual_infeasible();
    
    printf("\n=== Tests Complete ===\n\n");
    printf("Summary:\n");
    printf("- LP mode can detect unboundedness via dual/simplex checks\n");
    printf("- QP mode focuses on constraint feasibility\n");
    printf("- QR format inherits the same detection mechanisms\n");
    printf("- UNBNDORINFEAS = 'unbounded or infeasible' (may not distinguish)\n");
    return 0;
}
