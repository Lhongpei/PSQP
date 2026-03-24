/* Debug test for inf value in postsolve - detailed */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "PSQP_API.h"
#include "PSQP_sol.h"
#include "Postsolver.h"
#include "Postsolver.h"

int main() {
    printf("Debug: Testing postsolve with fully reduced problem\n\n");
    
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
    stgs->verbose = false;
    stgs->dual_fix = true;  /* Enable dual fix to trigger the bug */
    
    Presolver *presolver = new_qp_presolver_qr(Ax, Ai, Ap, m, n, 4,
                                               lhs, rhs, lbs, ubs, c,
                                               Qx, Qi, Qp, 4,
                                               NULL, NULL, NULL, 0, 0, stgs);
    
    if (!presolver) {
        printf("Failed to create presolver\n");
        return 1;
    }
    
    PresolveStatus status = run_presolver(presolver);
    printf("Presolve status: %d\n", status);
    
    PresolvedProblem *reduced = presolver->reduced_prob;
    printf("Reduced: m=%zu, n=%zu\n\n", reduced->m, reduced->n);
    
    /* Check col_map */
    Mapping *maps = presolver->prob->constraints->state->work->mappings;
    printf("Column mappings (col_map):\n");
    for (int i = 0; i < n; i++) {
        printf("  col_map[%d] = %d\n", i, maps->cols[i]);
    }
    printf("\n");
    
    /* Check reductions */
    PostsolveInfo *info = presolver->prob->constraints->state->postsolve_info;
    printf("Number of reductions: %zu\n", info->type->len);
    printf("Reductions:\n");
    
    const char* type_names[] = {
        "FIXED_COL", "FIXED_COL_INF", "SUB_COL", "PARALLEL_COL",
        "DELETED_ROW", "ADDED_ROW", "ADDED_ROWS", "LHS_CHANGE",
        "RHS_CHANGE", "EQ_TO_INEQ", "BOUND_CHANGE_NO_ROW",
        "BOUND_CHANGE_THE_ROW", "FIXED_COL_QP", "SUB_COL_QP"
    };
    
    for (size_t i = 0; i < info->type->len; i++) {
        int type = info->type->data[i];
        int start = info->starts->data[i];
        int end = info->starts->data[i+1];
        
        const char* name = "UNKNOWN";
        if (type >= 0 && type < 14) {
            /* Map bit flags to names */
            if (type == 0) name = "FIXED_COL";
            else if (type == 1) name = "FIXED_COL_INF";
            else if (type == 2) name = "SUB_COL";
            else if (type == 4) name = "PARALLEL_COL";
            else if (type == 8) name = "DELETED_ROW";
            else if (type == 16) name = "ADDED_ROW";
            else if (type == 32) name = "ADDED_ROWS";
            else if (type == 64) name = "LHS_CHANGE";
            else if (type == 128) name = "RHS_CHANGE";
            else if (type == 256) name = "EQ_TO_INEQ";
            else if (type == 512) name = "BOUND_CHANGE_NO_ROW";
            else if (type == 1024) name = "BOUND_CHANGE_THE_ROW";
            else if (type == 2048) name = "FIXED_COL_QP";
            else if (type == 4096) name = "SUB_COL_QP";
        }
        
        printf("  [%zu] type=%d (%s), indices=[", i, type, name);
        for (int j = start; j < end && j < start + 5; j++) {
            printf("%d", info->indices->data[j]);
            if (j < end - 1) printf(", ");
        }
        if (end - start > 5) printf("...");
        printf("]\n");
    }
    printf("\n");
    
    /* Even if reduced to 0, we should be able to postsolve */
    double *x_reduced = (double *)calloc(reduced->n ? reduced->n : 1, sizeof(double));
    double *y_reduced = (double *)calloc(reduced->m ? reduced->m : 1, sizeof(double));
    double *z_reduced = (double *)calloc(reduced->n ? reduced->n : 1, sizeof(double));
    
    postsolve(presolver, x_reduced, y_reduced, z_reduced);
    
    double *x_orig = presolver->sol->x;
    printf("Postsolved x: [%f, %f, %f, %f]\n", 
           x_orig[0], x_orig[1], x_orig[2], x_orig[3]);
    
    /* Check for inf/nan */
    int has_inf = 0;
    for (int i = 0; i < n; i++) {
        if (isinf(x_orig[i])) {
            printf("ERROR: x[%d] = inf\n", i);
            has_inf = 1;
        }
        if (isnan(x_orig[i])) {
            printf("ERROR: x[%d] = nan\n", i);
            has_inf = 1;
        }
    }
    
    if (!has_inf) {
        printf("SUCCESS: All values are finite\n");
    }
    
    free(x_reduced); free(y_reduced); free(z_reduced);
    free_settings(stgs);
    free_presolver(presolver);
    
    return has_inf;
}
