/* Comprehensive test for QR interface */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "PSQP/PSQP_API.h"
#include "PSQP/PSQP_sol.h"
#include "PSQP/PSQP_status.h"

int main() {
    printf("=== Comprehensive QR Test ===\n\n");
    
    Settings *stgs = default_settings();
    set_settings_true(stgs);
    stgs->verbose = false;
    
    int passed = 0, total = 0;
    
    // Test 1: R-only
    total++;
    printf("Test 1: R-only... ");
    {
        double Ax[] = {1.0, 1.0}, lhs[] = {1.0}, rhs[] = {1.0}, c[] = {0.0, 0.0};
        int Ai[] = {0, 0}, Ap[] = {0, 1, 2};
        double lbs[] = {0.0, 0.0}, ubs[] = {1e20, 1e20};
        double Rx[] = {1.0, 1.0}; int Ri[] = {0, 1}, Rp[] = {0, 2};
        
        Presolver *ps = new_qp_presolver_qr(Ax, Ai, Ap, 1, 2, 2, lhs, rhs, lbs, ubs, c,
                                            NULL, NULL, NULL, 0, Rx, Ri, Rp, 2, 1, stgs);
        if (ps && run_presolver(ps) != INFEASIBLE) passed++;
        else printf("FAILED ");
        if (ps) free_presolver(ps);
    }
    printf("OK\n");
    
    // Test 2: Q-only via unified interface
    total++;
    printf("Test 2: Q-only via unified... ");
    {
        double Ax[] = {1.0, 1.0}, lhs[] = {1.0}, rhs[] = {1.0}, c[] = {0.0, 0.0};
        int Ai[] = {0, 0}, Ap[] = {0, 1, 2};
        double lbs[] = {0.0, 0.0}, ubs[] = {1e20, 1e20};
        double Qx[] = {2.0, 2.0}; int Qi[] = {0, 1}, Qp[] = {0, 1, 2};
        
        Presolver *ps = new_qp_presolver(Ax, Ai, Ap, 1, 2, 2, lhs, rhs, lbs, ubs, c,
                                         Qx, Qi, Qp, 2, stgs);
        if (ps && run_presolver(ps) != INFEASIBLE) passed++;
        else printf("FAILED ");
        if (ps) free_presolver(ps);
    }
    printf("OK\n");
    
    // Test 3: Mixed Q+R
    total++;
    printf("Test 3: Mixed Q+R... ");
    {
        double Ax[] = {1.0, 1.0}, lhs[] = {1.0}, rhs[] = {1.0}, c[] = {0.0, 0.0};
        int Ai[] = {0, 0}, Ap[] = {0, 1, 2};
        double lbs[] = {0.0, 0.0}, ubs[] = {1e20, 1e20};
        double Qx[] = {2.0, 2.0}; int Qi[] = {0, 1}, Qp[] = {0, 1, 2};
        double Rx[] = {1.0, 1.0}; int Ri[] = {0, 1}, Rp[] = {0, 1, 2};
        
        Presolver *ps = new_qp_presolver_qr(Ax, Ai, Ap, 1, 2, 2, lhs, rhs, lbs, ubs, c,
                                            Qx, Qi, Qp, 2, Rx, Ri, Rp, 2, 2, stgs);
        if (ps && run_presolver(ps) != INFEASIBLE) passed++;
        else printf("FAILED ");
        if (ps) free_presolver(ps);
    }
    printf("OK\n");
    
    // Test 4: Pure LP
    total++;
    printf("Test 4: Pure LP... ");
    {
        double Ax[] = {1.0, 1.0}, lhs[] = {1.0}, rhs[] = {1.0}, c[] = {1.0, 1.0};
        int Ai[] = {0, 0}, Ap[] = {0, 1, 2};
        double lbs[] = {0.0, 0.0}, ubs[] = {1e20, 1e20};
        
        Presolver *ps = new_presolver(Ax, Ai, Ap, 1, 2, 2, lhs, rhs, lbs, ubs, c, stgs);
        if (ps && run_presolver(ps) != INFEASIBLE) passed++;
        else printf("FAILED ");
        if (ps) free_presolver(ps);
    }
    printf("OK\n");
    
    // Test 5: QP with mixed linear/quadratic vars (x1 quadratic, x2 linear)
    total++;
    printf("Test 5: Mixed linear/quadratic vars... ");
    {
        double Ax[] = {1.0, 1.0}, lhs[] = {1.0}, rhs[] = {1.0}, c[] = {0.0, 1.0};
        int Ai[] = {0, 0}, Ap[] = {0, 1, 2};
        double lbs[] = {0.0, 0.0}, ubs[] = {1e20, 1e20};
        double Qx[] = {1.0}; int Qi[] = {0}, Qp[] = {0, 1, 1};  // only x1 has Q
        
        Presolver *ps = new_qp_presolver(Ax, Ai, Ap, 1, 2, 2, lhs, rhs, lbs, ubs, c,
                                         Qx, Qi, Qp, 1, stgs);
        if (ps && run_presolver(ps) != INFEASIBLE) passed++;
        else printf("FAILED ");
        if (ps) free_presolver(ps);
    }
    printf("OK\n");
    
    free_settings(stgs);
    
    printf("\n=== Results: %d/%d passed ===\n", passed, total);
    return (passed == total) ? 0 : 1;
}
