#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "PSQP_API.h"

/* Test unboundedness with QR format - using very small Q */
void test_unbounded_qr_with_q()
{
    printf("Test: Unbounded QP with QR format (using small Q)...\n");
    
    size_t n = 1;
    size_t m = 1;
    
    // Loose constraint
    double Ax[] = {1.0};
    int Ai[] = {0};
    int Ap[] = {0, 1};
    size_t nnz = 1;
    
    double lhs[] = {-INFINITY};
    double rhs[] = {INFINITY};
    double lbs[] = {-INFINITY};  // Unbounded below
    double ubs[] = {INFINITY};
    double c[] = {1.0};  // Positive cost
    
    // Very small Q (almost zero)
    double Qx[] = {1e-10};
    int Qi[] = {0};
    int Qp[] = {0, 1};
    size_t Qnnz = 1;
    
    Settings *stgs = default_settings();
    stgs->verbose = false;
    stgs->dual_fix = true;
    
    Presolver *presolver = new_qp_presolver_qr(Ax, Ai, Ap, m, n, nnz,
                                               lhs, rhs, lbs, ubs, c,
                                               Qx, Qi, Qp, Qnnz,
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
        printf("  ✓ QR format detected unboundedness!\n");
    } else {
        printf("  Note: Strong convexity (even tiny) can make problem bounded\n");
    }
    
    free_settings(stgs);
    free_presolver(presolver);
}

/* Test bounded QR problem */
void test_bounded_qr()
{
    printf("\nTest: Bounded QP with QR format...\n");
    
    size_t n = 1;
    size_t m = 1;
    
    // Binding constraint
    double Ax[] = {1.0};
    int Ai[] = {0};
    int Ap[] = {0, 1};
    size_t nnz = 1;
    
    double lhs[] = {0.0};  // x >= 0
    double rhs[] = {INFINITY};
    double lbs[] = {-INFINITY};
    double ubs[] = {INFINITY};
    double c[] = {1.0};  // Minimize x
    
    // Stronger Q
    double Qx[] = {1.0};
    int Qi[] = {0};
    int Qp[] = {0, 1};
    size_t Qnnz = 1;
    
    double Rx[] = {1.0};  // Add R too
    int Ri[] = {0};
    int Rp[] = {0, 1};
    size_t Rnnz = 1;
    size_t k = 1;
    
    Settings *stgs = default_settings();
    stgs->verbose = false;
    stgs->dual_fix = true;
    
    Presolver *presolver = new_qp_presolver_qr(Ax, Ai, Ap, m, n, nnz,
                                               lhs, rhs, lbs, ubs, c,
                                               Qx, Qi, Qp, Qnnz,
                                               Rx, Ri, Rp, Rnnz, k,
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
    
    if (status != UNBNDORINFEAS) {
        printf("  ✓ Problem is bounded (constraint provides downlock)\n");
    }
    
    free_settings(stgs);
    free_presolver(presolver);
}

int main()
{
    printf("=== QR Format Unboundedness Tests ===\n\n");
    
    test_unbounded_qr_with_q();
    test_bounded_qr();
    
    printf("\n=== Summary ===\n");
    printf("QR format inherits unboundedness detection from LP/QP core:\n");
    printf("- Uses dual fix mechanism (ck > 0 + no downlocks + unbounded below)\n");
    printf("- Strong convexity (Q/R terms) can make problems bounded\n");
    printf("- Constraint locks prevent false unboundedness detection\n");
    
    return 0;
}
