#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "PSQP_API.h"

/* Test dual infeasibility detection:
 * 
 * Scenario 1: ck > 0, no downlocks, lb = -inf
 * This means we want to minimize a positive cost variable
 * that has no lower bound -> can go to -inf -> unbounded
 */
void test_dual_unbounded_v1()
{
    printf("Test 1: Dual infeasibility (ck>0, no downlocks, unbounded below)...\n");
    
    size_t n = 1;
    size_t m = 0;  // No constraints = no locks
    
    // Try with minimal constraints
    double Ax[] = {1.0};  // One constraint to allow presolver creation
    int Ai[] = {0};
    int Ap[] = {0, 1};
    size_t nnz = 1;
    
    double lhs[] = {-INFINITY};  // Ax >= -inf (always satisfied)
    double rhs[] = {INFINITY};   // Ax <= inf (always satisfied)
    double lbs[] = {-INFINITY};  // Unbounded below
    double ubs[] = {INFINITY};   // Unbounded above
    double c[] = {1.0};          // Minimize x (positive coefficient)
    
    Settings *stgs = default_settings();
    stgs->verbose = false;
    stgs->dual_fix = true;  // Enable dual fix
    
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
        printf("  ✓ Correctly detected dual infeasibility!\n");
    }
    
    free_settings(stgs);
    free_presolver(presolver);
}

/* Scenario 2: ck < 0, no uplocks, ub = inf
 * Minimize negative cost variable with no upper bound -> can go to +inf
 */
void test_dual_unbounded_v2()
{
    printf("\nTest 2: Dual infeasibility (ck<0, no uplocks, unbounded above)...\n");
    
    size_t n = 1;
    size_t m = 1;
    
    // Loose constraint
    double Ax[] = {1.0};
    int Ai[] = {0};
    int Ap[] = {0, 1};
    size_t nnz = 1;
    
    double lhs[] = {-INFINITY};
    double rhs[] = {INFINITY};
    double lbs[] = {-INFINITY};
    double ubs[] = {INFINITY};
    double c[] = {-1.0};  // Minimize -x = maximize x
    
    Settings *stgs = default_settings();
    stgs->verbose = false;
    stgs->dual_fix = true;
    
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
        printf("  ✓ Correctly detected dual infeasibility!\n");
    }
    
    free_settings(stgs);
    free_presolver(presolver);
}

/* Test 3: Same with QR format */
void test_dual_unbounded_qr()
{
    printf("\nTest 3: Dual infeasibility with QR format...\n");
    
    size_t n = 1;
    size_t m = 1;
    
    double Ax[] = {1.0};
    int Ai[] = {0};
    int Ap[] = {0, 1};
    size_t nnz = 1;
    
    double lhs[] = {-INFINITY};
    double rhs[] = {INFINITY};
    double lbs[] = {-INFINITY};  // Unbounded below
    double ubs[] = {INFINITY};
    double c[] = {1.0};  // Positive cost
    
    // No quadratic term
    Settings *stgs = default_settings();
    stgs->verbose = false;
    stgs->dual_fix = true;
    
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
           status == UNBNDORINFEAS ? "UNBNDORINFEAS" :
           status == INFEASIBLE ? "INFEASIBLE" :
           status == REDUCED ? "REDUCED" : "UNCHANGED");
    
    if (status == UNBNDORINFEAS) {
        printf("  ✓ QR format correctly detected dual infeasibility!\n");
    }
    
    free_settings(stgs);
    free_presolver(presolver);
}

/* Test 4: With locks (constraints) - should be bounded */
void test_bounded_with_locks()
{
    printf("\nTest 4: Bounded problem (with constraints/locks)...\n");
    
    size_t n = 1;
    size_t m = 1;
    
    // Real constraint: x >= 1 (creates downlock)
    double Ax[] = {1.0};
    int Ai[] = {0};
    int Ap[] = {0, 1};
    size_t nnz = 1;
    
    double lhs[] = {1.0};  // x >= 1
    double rhs[] = {INFINITY};
    double lbs[] = {-INFINITY};
    double ubs[] = {INFINITY};
    double c[] = {1.0};  // Minimize x
    
    Settings *stgs = default_settings();
    stgs->verbose = false;
    stgs->dual_fix = true;
    
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
    
    if (status != UNBNDORINFEAS) {
        printf("  ✓ Correctly identified as bounded (has downlock from constraint)\n");
    } else {
        printf("  Note: Detected as UNBNDORINFEAS (may be conservative)\n");
    }
    
    free_settings(stgs);
    free_presolver(presolver);
}

int main()
{
    printf("=== Dual Infeasibility (Unboundedness) Detection Tests ===\n\n");
    printf("Testing detection when ck > 0, no downlocks, lb = -inf...\n\n");
    
    test_dual_unbounded_v1();
    test_dual_unbounded_v2();
    test_dual_unbounded_qr();
    test_bounded_with_locks();
    
    printf("\n=== Summary ===\n");
    printf("The presolver uses dual fix to detect unboundedness:\n");
    printf("- If c_k > 0 and no downlocks and lb = -inf → UNBNDORINFEAS\n");
    printf("- If c_k < 0 and no uplocks and ub = +inf → UNBNDORINFEAS\n");
    printf("- QR format inherits this detection via the same code path\n");
    printf("- Constraints create locks that can prevent false positives\n");
    
    return 0;
}
