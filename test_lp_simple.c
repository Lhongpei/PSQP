#include <stdio.h>
#include <stdlib.h>
#include "PSQP/PSQP_API.h"

int main() {
    // Problem: min x + 2y s.t. x + y = 1, x>=0, y>=0
    // Expected: x=1, y=0, but we fix both to bounds that make them equal
    // Let's try a simpler case: just fix one variable
    
    double c[] = {1.0, 2.0};
    double lhs[] = {1.0};
    double rhs[] = {1.0};
    double lb[] = {1.0, 0.0};  // x fixed at 1
    double ub[] = {1.0, 0.0};  // y fixed at 0
    
    // A = [1, 1] (CSR format)
    double Ax[] = {1.0, 1.0};
    int Ai[] = {0, 1};
    int Ap[] = {0, 2};
    
    size_t m = 1;
    size_t n = 2;
    size_t nnz = 2;
    
    Settings *stgs = default_settings();
    stgs->verbose = true;
    
    printf("Creating LP presolver...\n");
    Presolver *presolver = new_presolver(Ax, Ai, Ap, m, n, nnz, lhs, rhs, lb, ub, c, stgs);
    
    if (!presolver) {
        printf("Failed to create presolver\n");
        return 1;
    }
    
    printf("Running presolver...\n");
    PresolveStatus status = run_presolver(presolver);
    printf("Presolve status: %d\n", status);
    
    if (presolver->reduced_prob) {
        printf("Reduced problem: m=%zu, n=%zu\n", presolver->reduced_prob->m, presolver->reduced_prob->n);
    } else {
        printf("No reduced problem\n");
    }
    
    if (presolver->reduced_prob && presolver->reduced_prob->n == 0) {
        printf("n=0, calling postsolve...\n");
        double x_reduced[0];
        double y_reduced[1] = {0.0};
        double z_reduced[0];
        postsolve(presolver, x_reduced, y_reduced, z_reduced);
        printf("Postsolve done\n");
    }
    
    free_presolver(presolver);
    free_settings(stgs);
    
    printf("Done\n");
    return 0;
}
