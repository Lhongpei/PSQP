#include <stdio.h>
#include <stdlib.h>
#include "PSQP_API.h"
#include "test_macros.h"

#define E2E_TOL 1e-6

static char *test_qr_infeasible()
{
    printf("Testing infeasible QP detection with QR...\n");
    
    size_t n = 2;
    size_t m = 2;
    
    /* Infeasible constraints:
     * x0 + x1 >= 5
     * x0 + x1 <= 3
     */
    double Ax[] = {1.0, 1.0, 1.0, 1.0};
    int Ai[] = {0, 1, 0, 1};
    int Ap[] = {0, 2, 4};
    size_t nnz = 4;
    
    double lhs[] = {5.0, -INFINITY};
    double rhs[] = {INFINITY, 3.0};
    double lbs[] = {0.0, 0.0};
    double ubs[] = {INFINITY, INFINITY};
    double c[] = {1.0, 1.0};
    
    /* QR format quadratic term */
    double Qx[] = {2.0, 2.0};
    int Qi[] = {0, 1};
    int Qp[] = {0, 1, 2};
    size_t Qnnz = 2;
    
    double Rx[] = {1.0, 1.0};
    int Ri[] = {0, 0};
    int Rp[] = {0, 1, 1};
    size_t Rnnz = 2;
    size_t k = 1;
    
    Settings *stgs = default_settings();
    stgs->verbose = false;
    
    Presolver *presolver = new_qp_presolver_qr(Ax, Ai, Ap, m, n, nnz,
                                               lhs, rhs, lbs, ubs, c,
                                               Qx, Qi, Qp, Qnnz,
                                               Rx, Ri, Rp, Rnnz, k, stgs);
    
    if (presolver == NULL) {
        return "presolver creation failed";
    }
    
    PresolveStatus status = run_presolver(presolver);
    
    printf("  Status: %s\n", 
           status == INFEASIBLE ? "INFEASIBLE" :
           status == REDUCED ? "REDUCED" :
           status == UNCHANGED ? "UNCHANGED" : 
           status == UNBNDORINFEAS ? "UNBNDORINFEAS" : "OTHER");
    
    int success = (status == INFEASIBLE || status == UNBNDORINFEAS);
    
    free_settings(stgs);
    free_presolver(presolver);
    
    if (!success) {
        return "infeasible QR QP should be detected";
    }
    
    return 0;
}

int main()
{
    printf("=== QR Infeasibility Test ===\n\n");
    
    char *result = test_qr_infeasible();
    
    if (result != 0) {
        printf("\nTEST FAILED: %s\n", result);
        return 1;
    }
    
    printf("\nTEST PASSED: Infeasibility correctly detected!\n");
    return 0;
}
