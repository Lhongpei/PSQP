/* Debug test for inf value in postsolve */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "PSQP_API.h"
#include "PSQP_sol.h"

int main() {
    printf("Debug: Testing postsolve with fully reduced problem\n");
    
    int m = 1, n = 4;
    double Ax[] = {1.0, 1.0, 1.0, 1.0};
    int Ai[] = {0, 1, 2, 3};
    int Ap[] = {0, 4};
    double lhs[] = {4.0};
    double rhs[] = {4.0};
    double lbs[] = {1.0, 0.0, 1.0, 0.0};  /* x0=1, x2=1 fixed */
    double ubs[] = {1.0, 10.0, 1.0, 10.0};
    double c[] = {0.0, 0.0, 0.0, 0.0};
    
    double Qx[] = {1.0, 1.0, 1.0, 1.0};
    int Qi[] = {0, 1, 2, 3};
    int Qp[] = {0, 1, 2, 3, 4};
    
    Settings *stgs = default_settings();
    stgs->verbose = true;  /* Enable verbose to see what happens */
    
    Presolver *presolver = new_qp_presolver_qr(Ax, Ai, Ap, m, n, 4,
                                               lhs, rhs, lbs, ubs, c,
                                               Qx, Qi, Qp, 4,
                                               NULL, NULL, NULL, 0, 0, stgs);
    
    if (!presolver) {
        printf("Failed to create presolver\n");
        return 1;
    }
    
    PresolveStatus status = run_presolver(presolver);
    printf("\nPresolve status: %d\n", status);
    
    PresolvedProblem *reduced = presolver->reduced_prob;
    printf("Reduced: m=%zu, n=%zu\n", reduced->m, reduced->n);
    
    /* Even if reduced to 0, we should be able to postsolve */
    double *x_reduced = (double *)calloc(reduced->n ? reduced->n : 1, sizeof(double));
    double *y_reduced = (double *)calloc(reduced->m ? reduced->m : 1, sizeof(double));
    double *z_reduced = (double *)calloc(reduced->n ? reduced->n : 1, sizeof(double));
    
    postsolve(presolver, x_reduced, y_reduced, z_reduced);
    
    double *x_orig = presolver->sol->x;
    printf("\nPostsolved x: [%f, %f, %f, %f]\n", 
           x_orig[0], x_orig[1], x_orig[2], x_orig[3]);
    
    /* Check for inf/nan */
    for (int i = 0; i < n; i++) {
        if (isinf(x_orig[i])) {
            printf("ERROR: x[%d] = inf\n", i);
        }
        if (isnan(x_orig[i])) {
            printf("ERROR: x[%d] = nan\n", i);
        }
    }
    
    free(x_reduced); free(y_reduced); free(z_reduced);
    free_settings(stgs);
    free_presolver(presolver);
    
    return 0;
}
