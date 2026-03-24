/*
 * QP Debug Test - Step by step
 */

#include <stdio.h>
#include <stdlib.h>
#include "PSQP_API.h"
#include "Problem.h"

int main() {
    printf("=== QP Debug Test ===\n");
    
    printf("Step 1: Create settings\n");
    Settings *stgs = default_settings();
    if (!stgs) {
        printf("ERROR: Failed to create settings\n");
        return 1;
    }
    printf("  Settings at %p\n", (void*)stgs);
    stgs->verbose = false;
    
    printf("Step 2: Prepare problem data\n");
    int n = 2, m = 1;
    double c[] = {1.0, 2.0};
    double lb[] = {0.0, 0.0};
    double ub[] = {10.0, 10.0};
    double Qx[] = {2.0, 3.0};
    int Qi[] = {0, 1};
    int Qp[] = {0, 1, 2};
    double A_x[] = {1.0, 1.0};
    int A_i[] = {0, 0};
    int A_p[] = {0, 2};
    double lhs[] = {3.0};
    double rhs[] = {3.0};
    
    printf("Step 3: Create QP presolver\n");
    Presolver *presolver = new_qp_presolver_qr(
        A_x, A_i, A_p, m, n, 2,
        lhs, rhs, lb, ub, c,
        Qx, Qi, Qp, 2,
        NULL, NULL, NULL, 0, 0,
        stgs
    );
    
    if (!presolver) {
        printf("ERROR: Failed to create presolver\n");
        free_settings(stgs);
        return 1;
    }
    printf("  Presolver at %p\n", (void*)presolver);
    
    printf("Step 4: Check presolver fields\n");
    printf("  presolver->stgs = %p\n", (void*)presolver->stgs);
    printf("  presolver->prob = %p\n", (void*)presolver->prob);
    
    if (!presolver->prob) {
        printf("ERROR: presolver->prob is NULL\n");
        free_presolver(presolver);
        free_settings(stgs);
        return 1;
    }
    
    printf("Step 5: Check problem fields\n");
    printf("  presolver->prob->obj = %p\n", (void*)presolver->prob->obj);
    
    if (!presolver->prob->obj) {
        printf("ERROR: presolver->prob->obj is NULL\n");
        free_presolver(presolver);
        free_settings(stgs);
        return 1;
    }
    
    printf("Step 6: Check objective fields\n");
    printf("  presolver->prob->obj->quad_qr = %p\n", (void*)presolver->prob->obj->quad_qr);
    
    if (!presolver->prob->obj->quad_qr) {
        printf("ERROR: presolver->prob->obj->quad_qr is NULL\n");
        free_presolver(presolver);
        free_settings(stgs);
        return 1;
    }
    
    printf("Step 7: Check QR fields\n");
    QuadTermQR *qr = presolver->prob->obj->quad_qr;
    printf("  qr->n = %zu\n", qr->n);
    printf("  qr->Qnnz = %zu\n", qr->Qnnz);
    printf("  qr->has_quad = %d\n", qr->has_quad);
    
    printf("Step 8: Cleanup\n");
    free_presolver(presolver);
    free_settings(stgs);
    
    printf("SUCCESS!\n");
    return 0;
}
