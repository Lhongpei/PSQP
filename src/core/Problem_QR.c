/*
 * Copyright 2025-2026 Daniel Cederberg
 * Copyright 2026 Hongpei Li
 *
 * This file is part of the PSQP project (QP Presolver).
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "Problem.h"
#include "Memory_wrapper.h"
#include "Bounds.h"
#include <math.h>
#include <string.h>

/* Helper: Transpose a sparse matrix in CSR format */
static void transpose_csr(const double *Ax, const int *Ai, const int *Ap,
                          double **ATx, int **ATi, int **ATp,
                          size_t n_rows, size_t n_cols, size_t nnz)
{
    *ATx = (double *) ps_malloc(nnz, sizeof(double));
    *ATi = (int *) ps_malloc(nnz, sizeof(int));
    *ATp = (int *) ps_malloc(n_cols + 1, sizeof(int));
    
    if (!*ATx || !*ATi || !*ATp) return;
    
    /* Count elements per column */
    for (size_t i = 0; i <= n_cols; i++) (*ATp)[i] = 0;
    for (size_t i = 0; i < nnz; i++) (*ATp)[Ai[i] + 1]++;
    
    /* Cumulative sum */
    for (size_t i = 0; i < n_cols; i++) (*ATp)[i + 1] += (*ATp)[i];
    
    /* Fill values */
    int *col_counts = (int *) ps_malloc(n_cols, sizeof(int));
    if (!col_counts) return;
    
    for (size_t i = 0; i < n_cols; i++) col_counts[i] = 0;
    
    for (size_t i = 0; i < n_rows; i++) {
        for (int j = Ap[i]; j < Ap[i + 1]; j++) {
            int col = Ai[j];
            int idx = (*ATp)[col] + col_counts[col];
            (*ATi)[idx] = (int) i;
            (*ATx)[idx] = Ax[j];
            col_counts[col]++;
        }
    }
    PS_FREE(col_counts);
}

/* ============================================================================
 * QuadTermQR: P = Q + R*R^T implementation
 * ============================================================================ */

QuadTermQR *quad_term_qr_new(const double *Qx, const int *Qi, const int *Qp, size_t Qnnz,
                             const double *Rx, const int *Ri, const int *Rp, size_t Rnnz,
                             size_t n, size_t k)
{
    if ((Qnnz == 0 || Qx == NULL) && (Rnnz == 0 || Rx == NULL))
    {
        return NULL;
    }

    QuadTermQR *qr = (QuadTermQR *) ps_malloc(1, sizeof(QuadTermQR));
    RETURN_PTR_IF_NULL(qr, NULL);
    
    memset(qr, 0, sizeof(QuadTermQR));
    qr->n = n;
    qr->k = k;
    qr->has_quad = true;

    /* Copy Q matrix */
    if (Qx != NULL && Qnnz > 0)
    {
        qr->Qnnz = Qnnz;
        qr->Qx = (double *) ps_malloc(Qnnz, sizeof(double));
        qr->Qi = (int *) ps_malloc(Qnnz, sizeof(int));
        qr->Qp = (int *) ps_malloc(n + 1, sizeof(int));
        
        if (!qr->Qx || !qr->Qi || !qr->Qp)
        {
            free_quad_term_qr(qr);
            return NULL;
        }
        
        memcpy(qr->Qx, Qx, Qnnz * sizeof(double));
        memcpy(qr->Qi, Qi, Qnnz * sizeof(int));
        memcpy(qr->Qp, Qp, (n + 1) * sizeof(int));
    }

    /* Copy R matrix. R is n×k CSR: user must pass Rp of size n+1.
     * We take Rp[n] as the authoritative nnz if it fits within Rnnz. */
    if (Rx != NULL && Rnnz > 0 && Rp != NULL)
    {
        qr->Rnnz = Rnnz;
        qr->Rx = (double *) ps_malloc(Rnnz, sizeof(double));
        qr->Ri = (int *) ps_malloc(Rnnz, sizeof(int));
        qr->Rp = (int *) ps_malloc(n + 1, sizeof(int));

        if (!qr->Rx || !qr->Ri || !qr->Rp)
        {
            free_quad_term_qr(qr);
            return NULL;
        }

        memcpy(qr->Rx, Rx, Rnnz * sizeof(double));
        memcpy(qr->Ri, Ri, Rnnz * sizeof(int));
        memcpy(qr->Rp, Rp, (n + 1) * sizeof(int));

        /* Reconcile Rnnz with Rp[n]: trust Rp as the authoritative structure. */
        size_t rp_nnz = (size_t)qr->Rp[n];
        if (rp_nnz < qr->Rnnz) qr->Rnnz = rp_nnz;

        /* Build R transpose (k×n) for efficient column access of R */
        transpose_csr(qr->Rx, qr->Ri, qr->Rp, &qr->RTx, &qr->RTi, &qr->RTp,
                      n, k, qr->Rnnz);
        if (!qr->RTx)
        {
            free_quad_term_qr(qr);
            return NULL;
        }

        /* Clean small entries */
        clean_zero_entries_qr(qr, 1e-13);

        /* Extract single-nnz factors (columns of R) into Q diagonal. */
        extract_single_nnz_rows_from_R(qr);

        /* Merge collinear factors (columns of R) if we still have >= 2. */
        if (qr->k > 1 && qr->Rnnz > 0)
        {
            merge_collinear_rows_in_R(qr, 1e-10);
        }
    }

    return qr;
}

void free_quad_term_qr(QuadTermQR *qr)
{
    RETURN_IF_NULL(qr);
    
    PS_FREE(qr->Qx);
    PS_FREE(qr->Qi);
    PS_FREE(qr->Qp);
    
    PS_FREE(qr->Rx);
    PS_FREE(qr->Ri);
    PS_FREE(qr->Rp);
    
    PS_FREE(qr->RTx);
    PS_FREE(qr->RTi);
    PS_FREE(qr->RTp);
    
    PS_FREE(qr);
}

/* Compute P_ii = Q_ii + sum_j R_ij^2 */
double compute_p_diag(const QuadTermQR *qr, int col)
{
    if (!qr || !qr->has_quad) return 0.0;
    
    double result = 0.0;
    
    /* Q contribution: Q_ii */
    if (qr->Qx != NULL && qr->Qp != NULL)
    {
        int q_start = qr->Qp[col];
        int q_end = qr->Qp[col + 1];
        for (int idx = q_start; idx < q_end; idx++)
        {
            if (qr->Qi[idx] == col)  /* Diagonal element */
            {
                result += qr->Qx[idx];
                break;
            }
        }
    }
    
    /* R*R^T contribution: sum_j R_ij^2 */
    if (qr->Rx != NULL && qr->Rp != NULL)
    {
        int r_start = qr->Rp[col];
        int r_end = qr->Rp[col + 1];
        for (int idx = r_start; idx < r_end; idx++)
        {
            double val = qr->Rx[idx];
            result += val * val;
        }
    }
    
    return result;
}

/* Add contribution from fixed variable to linear term c
 * For fixed x_i, we add to c[j]: P_ij * x_i = (Q_ij + sum_l R_il*R_jl) * x_i
 */
void add_fixed_var_contribution_qr(QuadTermQR *qr, int fixed_col, double value, double *c)
{
    if (!qr || !qr->has_quad || value == 0.0) return;
    
    /* Q contribution: for each j, add Q_{fixed_col,j} * value to c[j] */
    if (qr->Qx != NULL && qr->Qp != NULL)
    {
        int q_start = qr->Qp[fixed_col];
        int q_end = qr->Qp[fixed_col + 1];
        for (int idx = q_start; idx < q_end; idx++)
        {
            int j = qr->Qi[idx];
            double q_val = qr->Qx[idx];
            c[j] += q_val * value;
        }
    }
    
    /* R*R^T contribution: for each j, add (sum_l R_{fixed_col,l} * R_{j,l}) * value */
    if (qr->Rx != NULL && qr->Rp != NULL && qr->RTp != NULL)
    {
        /* Iterate over non-zeros in row fixed_col of R */
        int r_row_start = qr->Rp[fixed_col];
        int r_row_end = qr->Rp[fixed_col + 1];
        
        for (int idx_i = r_row_start; idx_i < r_row_end; idx_i++)
        {
            double r_il = qr->Rx[idx_i];
            int l = qr->Ri[idx_i];  /* column l in R */
            
            /* Now iterate over all j where R_{j,l} != 0 using RT */
            /* RT is k x n, stored CSR. Row l of RT is column l of R. */
            int rt_row_start = qr->RTp[l];
            int rt_row_end = qr->RTp[l + 1];
            
            for (int idx_j = rt_row_start; idx_j < rt_row_end; idx_j++)
            {
                int j = qr->RTi[idx_j];
                double r_jl = qr->RTx[idx_j];
                
                /* Add R_{fixed_col,l} * R_{j,l} * value to c[j] */
                c[j] += r_il * r_jl * value;
            }
        }
    }
}

/* Update offset when variable is fixed */
void fix_var_in_obj_qp_qr(QuadTermQR *qr, int col, double value, double *offset, double *c)
{
    if (!qr || !qr->has_quad || value == 0.0) return;
    
    /* Compute P_ii = Q_ii + sum_j R_ij^2 */
    double p_ii = compute_p_diag(qr, col);
    
    /* offset += 0.5 * P_ii * value^2 */
    *offset += 0.5 * p_ii * value * value;
    
    /* Update linear term c from off-diagonal contributions */
    add_fixed_var_contribution_qr(qr, col, value, c);
}

/* Check if variable has R terms */
bool has_r_terms(const QuadTermQR *qr, int col)
{
    if (!qr || !qr->has_quad || qr->Rp == NULL) return false;
    
    int r_start = qr->Rp[col];
    int r_end = qr->Rp[col + 1];
    return (r_end > r_start);
}

/* Shrink Q and R when variables are removed */
void quad_term_qr_shrink(QuadTermQR *qr, int *col_map, size_t n_cols_old)
{
    if (!qr || !qr->has_quad) return;
    
    /* Count new columns */
    size_t n_cols_new = 0;
    for (size_t i = 0; i < n_cols_old; i++)
    {
        if (col_map[i] >= 0) n_cols_new++;
    }
    
    if (n_cols_new == n_cols_old) return;  /* Nothing to do */
    
    qr->n = n_cols_new;
    
    if (n_cols_new == 0)
    {
        /* All columns removed */
        PS_FREE(qr->Qx); qr->Qx = NULL;
        PS_FREE(qr->Qi); qr->Qi = NULL;
        PS_FREE(qr->Qp); qr->Qp = NULL;
        PS_FREE(qr->Rx); qr->Rx = NULL;
        PS_FREE(qr->Ri); qr->Ri = NULL;
        PS_FREE(qr->Rp); qr->Rp = NULL;
        PS_FREE(qr->RTx); qr->RTx = NULL;
        PS_FREE(qr->RTi); qr->RTi = NULL;
        PS_FREE(qr->RTp); qr->RTp = NULL;
        qr->Qnnz = 0;
        qr->Rnnz = 0;
        qr->has_quad = false;
        return;
    }
    
    /* Shrink Q matrix: rebuild CSR with new indices */
    if (qr->Qx != NULL && qr->Qnnz > 0)
    {
        /* First pass: count new nnz (only keeping entries where both row and col are active) */
        size_t new_Qnnz = 0;
        for (size_t i = 0; i < n_cols_old; i++)
        {
            if (col_map[i] < 0) continue;  /* Skip deleted rows */
            
            int row_start = qr->Qp[i];
            int row_end = qr->Qp[i + 1];
            for (int idx = row_start; idx < row_end; idx++)
            {
                int col = qr->Qi[idx];
                if (col_map[col] >= 0)  /* Keep only if column is also active */
                {
                    new_Qnnz++;
                }
            }
        }
        
        if (new_Qnnz == 0)
        {
            /* No Q terms left */
            PS_FREE(qr->Qx);
            PS_FREE(qr->Qi);
            PS_FREE(qr->Qp);
            qr->Qnnz = 0;
        }
        else
        {
            /* Allocate new arrays */
            double *new_Qx = (double *) ps_malloc(new_Qnnz, sizeof(double));
            int *new_Qi = (int *) ps_malloc(new_Qnnz, sizeof(int));
            int *new_Qp = (int *) ps_malloc(n_cols_new + 1, sizeof(int));
            
            if (!new_Qx || !new_Qi || !new_Qp)
            {
                PS_FREE(new_Qx);
                PS_FREE(new_Qi);
                PS_FREE(new_Qp);
                return;
            }
            
            /* Second pass: fill new Q matrix */
            new_Qp[0] = 0;
            size_t idx_new = 0;
            
            for (size_t i = 0; i < n_cols_old; i++)
            {
                if (col_map[i] < 0) continue;  /* Skip deleted rows */
                
                int new_row = col_map[i];
                int row_start = qr->Qp[i];
                int row_end = qr->Qp[i + 1];
                
                for (int idx = row_start; idx < row_end; idx++)
                {
                    int col = qr->Qi[idx];
                    if (col_map[col] >= 0)  /* Keep only if column is also active */
                    {
                        new_Qx[idx_new] = qr->Qx[idx];
                        new_Qi[idx_new] = col_map[col];  /* New column index */
                        idx_new++;
                    }
                }
                new_Qp[new_row + 1] = (int) idx_new;
            }
            
            /* Free old arrays and assign new ones */
            PS_FREE(qr->Qx);
            PS_FREE(qr->Qi);
            PS_FREE(qr->Qp);
            qr->Qx = new_Qx;
            qr->Qi = new_Qi;
            qr->Qp = new_Qp;
            qr->Qnnz = new_Qnnz;
        }
    }
    
    /* Shrink R matrix: rebuild CSR with new row indices */
    if (qr->Rx != NULL && qr->Rnnz > 0)
    {
        /* First pass: count new nnz (only keeping entries where row is active) */
        size_t new_Rnnz = 0;
        for (size_t i = 0; i < n_cols_old; i++)
        {
            if (col_map[i] < 0) continue;  /* Skip deleted rows */
            
            int row_start = qr->Rp[i];
            int row_end = qr->Rp[i + 1];
            new_Rnnz += (size_t)(row_end - row_start);
        }
        
        if (new_Rnnz == 0)
        {
            /* No R terms left */
            PS_FREE(qr->Rx);
            PS_FREE(qr->Ri);
            PS_FREE(qr->Rp);
            PS_FREE(qr->RTx);
            PS_FREE(qr->RTi);
            PS_FREE(qr->RTp);
            qr->Rx = NULL; qr->Ri = NULL; qr->Rp = NULL;
            qr->RTx = NULL; qr->RTi = NULL; qr->RTp = NULL;
            qr->Rnnz = 0;
        }
        else
        {
            /* Allocate new arrays */
            double *new_Rx = (double *) ps_malloc(new_Rnnz, sizeof(double));
            int *new_Ri = (int *) ps_malloc(new_Rnnz, sizeof(int));
            int *new_Rp = (int *) ps_malloc(n_cols_new + 1, sizeof(int));
            
            if (!new_Rx || !new_Ri || !new_Rp)
            {
                PS_FREE(new_Rx);
                PS_FREE(new_Ri);
                PS_FREE(new_Rp);
                return;
            }
            
            /* Second pass: fill new R matrix */
            new_Rp[0] = 0;
            size_t idx_new = 0;
            
            for (size_t i = 0; i < n_cols_old; i++)
            {
                if (col_map[i] < 0) continue;  /* Skip deleted rows */
                
                int new_row = col_map[i];
                int row_start = qr->Rp[i];
                int row_end = qr->Rp[i + 1];
                
                for (int idx = row_start; idx < row_end; idx++)
                {
                    new_Rx[idx_new] = qr->Rx[idx];
                    new_Ri[idx_new] = qr->Ri[idx];  /* Column indices in R don't change */
                    idx_new++;
                }
                new_Rp[new_row + 1] = (int) idx_new;
            }
            
            /* Free old arrays and assign new ones */
            PS_FREE(qr->Rx);
            PS_FREE(qr->Ri);
            PS_FREE(qr->Rp);
            PS_FREE(qr->RTx);
            PS_FREE(qr->RTi);
            PS_FREE(qr->RTp);
            qr->Rx = new_Rx;
            qr->Ri = new_Ri;
            qr->Rp = new_Rp;
            qr->Rnnz = new_Rnnz;
            qr->RTx = NULL; qr->RTi = NULL; qr->RTp = NULL;
            
            /* Rebuild R transpose */
            transpose_csr(qr->Rx, qr->Ri, qr->Rp, &qr->RTx, &qr->RTi, &qr->RTp,
                          n_cols_new, qr->k, qr->Rnnz);
        }
    }
    
    /* Update has_quad flag */
    qr->has_quad = (qr->Qnnz > 0) || (qr->Rnnz > 0);
}

/* Remove small entries below tolerance from R matrix.
 * R is stored as n×k CSR (n rows = variables, k cols = factors). */
void clean_zero_entries_qr(QuadTermQR *qr, double tol)
{
    if (!qr || !qr->Rx || qr->Rnnz == 0) return;
    if (tol <= 0) tol = 1e-13;

    int n = (int)qr->n;
    int k = (int)qr->k;

    /* First pass: count entries above tolerance */
    size_t new_nnz = 0;
    for (size_t i = 0; i < qr->Rnnz; i++)
    {
        if (fabs(qr->Rx[i]) > tol)
            new_nnz++;
    }

    if (new_nnz == qr->Rnnz) return;  /* Nothing to remove */

    /* Allocate new arrays (n+1 row pointers for n×k CSR) */
    double *new_Rx = (double *)ps_malloc(new_nnz > 0 ? new_nnz : 1, sizeof(double));
    int *new_Ri = (int *)ps_malloc(new_nnz > 0 ? new_nnz : 1, sizeof(int));
    int *new_Rp = (int *)ps_malloc((size_t)n + 1, sizeof(int));

    if (!new_Rx || !new_Ri || !new_Rp)
    {
        PS_FREE(new_Rx);
        PS_FREE(new_Ri);
        PS_FREE(new_Rp);
        return;
    }

    /* Second pass: copy surviving entries (iterate all n rows) */
    new_Rp[0] = 0;
    size_t idx_new = 0;
    for (int row = 0; row < n; row++)
    {
        int row_start = qr->Rp[row];
        int row_end = qr->Rp[row + 1];

        for (int idx = row_start; idx < row_end; idx++)
        {
            if (fabs(qr->Rx[idx]) > tol)
            {
                new_Rx[idx_new] = qr->Rx[idx];
                new_Ri[idx_new] = qr->Ri[idx];
                idx_new++;
            }
        }
        new_Rp[row + 1] = (int)idx_new;
    }

    /* Free old arrays */
    PS_FREE(qr->Rx);
    PS_FREE(qr->Ri);
    PS_FREE(qr->Rp);
    PS_FREE(qr->RTx);
    PS_FREE(qr->RTi);
    PS_FREE(qr->RTp);

    qr->Rx = new_Rx;
    qr->Ri = new_Ri;
    qr->Rp = new_Rp;
    qr->Rnnz = new_nnz;
    qr->RTx = NULL;
    qr->RTi = NULL;
    qr->RTp = NULL;

    /* Rebuild transpose (k×n) if any entries remain */
    if (new_nnz > 0)
        transpose_csr(qr->Rx, qr->Ri, qr->Rp, &qr->RTx, &qr->RTi, &qr->RTp, (size_t)n, (size_t)k, new_nnz);
}

/* Extract rows of R that have only a single non-zero entry and move them to Q diagonal.
 * If the target column already has a Q diagonal entry, the values are ADDED.
 * Returns number of rows extracted.
 */
/* For n×k R, find factors (columns of R = rows of RT) that have a single
 * non-zero entry. If factor l has only R[v,l] = val, then factor l's
 * rank-1 contribution to R R^T is val^2 · e_v e_v^T (diagonal only),
 * so we add val^2 to Q[v,v] and remove that column from R (reducing k).
 *
 * Returns the number of factors extracted. */
size_t extract_single_nnz_rows_from_R(QuadTermQR *qr)
{
    if (!qr || !qr->Rx || qr->Rnnz == 0 || qr->k == 0) return 0;
    if (!qr->RTp || !qr->RTi || !qr->RTx) return 0;

    int n = (int)qr->n;
    int k = (int)qr->k;

    /* Count factors (RT rows) with exactly one non-zero. */
    int single_cnt = 0;
    for (int l = 0; l < k; l++)
    {
        if (qr->RTp[l + 1] - qr->RTp[l] == 1) single_cnt++;
    }
    if (single_cnt == 0) return 0;

    /* factor_keep[l] = new factor index, or -1 if extracted */
    int *factor_keep = (int *)ps_malloc((size_t)k, sizeof(int));
    if (!factor_keep) return 0;

    int new_k = 0;
    for (int l = 0; l < k; l++)
    {
        if (qr->RTp[l + 1] - qr->RTp[l] == 1) factor_keep[l] = -1;
        else                                  factor_keep[l] = new_k++;
    }
    int k_new = new_k;

    /* For each diagonal (variable), accumulate val^2 to add to Q[v,v]. */
    double *diag_add = (double *)ps_calloc((size_t)n, sizeof(double));
    if (!diag_add) { PS_FREE(factor_keep); return 0; }

    for (int l = 0; l < k; l++)
    {
        if (factor_keep[l] >= 0) continue;
        int rt_idx = qr->RTp[l];
        int v = qr->RTi[rt_idx];
        double val = qr->RTx[rt_idx];
        diag_add[v] += val * val;
    }

    /* Locate existing Q diagonal entries (by variable). */
    int *q_diag_idx = (int *)ps_malloc((size_t)n, sizeof(int));
    if (!q_diag_idx) { PS_FREE(factor_keep); PS_FREE(diag_add); return 0; }
    for (int i = 0; i < n; i++) q_diag_idx[i] = -1;

    if (qr->Qx && qr->Qp)
    {
        for (int v = 0; v < n; v++)
        {
            for (int idx = qr->Qp[v]; idx < qr->Qp[v + 1]; idx++)
            {
                if (qr->Qi[idx] == v) { q_diag_idx[v] = idx; break; }
            }
        }
    }

    /* Count variables that need a new diagonal slot. */
    int new_diag_needed = 0;
    for (int v = 0; v < n; v++)
    {
        if (diag_add[v] != 0.0 && q_diag_idx[v] < 0) new_diag_needed++;
    }

    size_t new_Qnnz = qr->Qnnz + (size_t)new_diag_needed;

    double *new_Qx = (double *)ps_malloc(new_Qnnz > 0 ? new_Qnnz : 1, sizeof(double));
    int    *new_Qi = (int *)   ps_malloc(new_Qnnz > 0 ? new_Qnnz : 1, sizeof(int));
    int    *new_Qp = (int *)   ps_malloc((size_t)n + 1, sizeof(int));
    if (!new_Qx || !new_Qi || !new_Qp)
    {
        PS_FREE(factor_keep); PS_FREE(diag_add); PS_FREE(q_diag_idx);
        PS_FREE(new_Qx); PS_FREE(new_Qi); PS_FREE(new_Qp);
        return 0;
    }

    /* Rebuild Q row by row. Q is upper-triangular CSR with entries sorted by
     * column index within each row; insert a new diagonal as the first entry
     * of that row when the diagonal slot didn't previously exist. */
    new_Qp[0] = 0;
    size_t q_new = 0;
    for (int v = 0; v < n; v++)
    {
        bool diag_existed = (q_diag_idx[v] >= 0);
        bool need_new_diag = (diag_add[v] != 0.0) && !diag_existed;

        int row_start = (qr->Qp ? qr->Qp[v] : 0);
        int row_end   = (qr->Qp ? qr->Qp[v + 1] : 0);

        if (need_new_diag)
        {
            /* Insert diagonal first (smallest col index). */
            new_Qx[q_new] = diag_add[v];
            new_Qi[q_new] = v;
            q_new++;
        }

        for (int idx = row_start; idx < row_end; idx++)
        {
            new_Qx[q_new] = qr->Qx[idx];
            new_Qi[q_new] = qr->Qi[idx];
            if (qr->Qi[idx] == v && diag_existed)
            {
                /* Augment existing diagonal with accumulated val^2. */
                new_Qx[q_new] += diag_add[v];
            }
            q_new++;
        }
        new_Qp[v + 1] = (int)q_new;
    }

    PS_FREE(qr->Qx); PS_FREE(qr->Qi); PS_FREE(qr->Qp);
    qr->Qx = new_Qx;
    qr->Qi = new_Qi;
    qr->Qp = new_Qp;
    qr->Qnnz = new_Qnnz;

    PS_FREE(q_diag_idx);
    PS_FREE(diag_add);

    /* Now rebuild R to drop extracted factor columns and remap remaining
     * factor indices. R is n×k CSR; we walk each variable row and keep only
     * entries whose factor is not extracted. */
    if (k_new == 0)
    {
        PS_FREE(qr->Rx); PS_FREE(qr->Ri); PS_FREE(qr->Rp);
        PS_FREE(qr->RTx); PS_FREE(qr->RTi); PS_FREE(qr->RTp);
        qr->Rx = NULL; qr->Ri = NULL; qr->Rp = NULL;
        qr->RTx = NULL; qr->RTi = NULL; qr->RTp = NULL;
        qr->Rnnz = 0;
        qr->k = 0;
    }
    else
    {
        /* Count surviving entries. */
        size_t new_Rnnz = 0;
        for (size_t i = 0; i < qr->Rnnz; i++)
        {
            if (factor_keep[qr->Ri[i]] >= 0) new_Rnnz++;
        }

        double *new_Rx = (double *)ps_malloc(new_Rnnz > 0 ? new_Rnnz : 1, sizeof(double));
        int    *new_Ri = (int *)   ps_malloc(new_Rnnz > 0 ? new_Rnnz : 1, sizeof(int));
        int    *new_Rp = (int *)   ps_malloc((size_t)n + 1, sizeof(int));
        if (!new_Rx || !new_Ri || !new_Rp)
        {
            PS_FREE(new_Rx); PS_FREE(new_Ri); PS_FREE(new_Rp);
            PS_FREE(factor_keep);
            return 0;
        }

        new_Rp[0] = 0;
        size_t idx_new = 0;
        for (int v = 0; v < n; v++)
        {
            int row_start = qr->Rp[v];
            int row_end   = qr->Rp[v + 1];
            for (int idx = row_start; idx < row_end; idx++)
            {
                int old_factor = qr->Ri[idx];
                int new_factor = factor_keep[old_factor];
                if (new_factor < 0) continue;  /* dropped */
                new_Rx[idx_new] = qr->Rx[idx];
                new_Ri[idx_new] = new_factor;
                idx_new++;
            }
            new_Rp[v + 1] = (int)idx_new;
        }

        PS_FREE(qr->Rx); PS_FREE(qr->Ri); PS_FREE(qr->Rp);
        PS_FREE(qr->RTx); PS_FREE(qr->RTi); PS_FREE(qr->RTp);
        qr->Rx = new_Rx;
        qr->Ri = new_Ri;
        qr->Rp = new_Rp;
        qr->Rnnz = new_Rnnz;
        qr->k = (size_t)k_new;
        qr->RTx = NULL; qr->RTi = NULL; qr->RTp = NULL;
        if (new_Rnnz > 0)
            transpose_csr(qr->Rx, qr->Ri, qr->Rp, &qr->RTx, &qr->RTi, &qr->RTp,
                          (size_t)n, (size_t)k_new, new_Rnnz);
    }

    PS_FREE(factor_keep);
    return (size_t)single_cnt;
}

/* Simple hash function for sparse row */
static uint64_t hash_sparse_row(const double *vals, const int *cols, int nnz, double scale, double tol)
{
    uint64_t h = 14695981039346656037ULL;
    const uint64_t FNV_PRIME = 1099511628211ULL;
    
    for (int i = 0; i < nnz; i++)
    {
        double normalized = (fabs(vals[i]) > tol) ? vals[i] * scale : 0.0;
        int col = cols[i];
        
        /* Hash column */
        h ^= (uint64_t)col;
        h *= FNV_PRIME;
        
        /* Hash value (as bits) */
        uint64_t val_bits;
        memcpy(&val_bits, &normalized, sizeof(double));
        h ^= val_bits;
        h *= FNV_PRIME;
    }
    
    return h;
}

/* Check if two sparse rows are collinear within tolerance */
static bool rows_are_collinear(const double *vals1, const int *cols1, int nnz1,
                               const double *vals2, const int *cols2, int nnz2,
                               double tol, double *alpha)
{
    if (nnz1 != nnz2) return false;
    if (nnz1 == 0) return false;
    
    /* Check same sparsity pattern */
    for (int i = 0; i < nnz1; i++)
    {
        if (cols1[i] != cols2[i])
            return false;
    }
    
    /* Compute ratio of first non-zero elements */
    int first_idx = 0;
    while (first_idx < nnz1 && fabs(vals1[first_idx]) < tol)
        first_idx++;
    if (first_idx >= nnz1) return false;
    
    double ratio = vals2[first_idx] / vals1[first_idx];
    
    /* Check if all elements match this ratio */
    for (int i = first_idx; i < nnz1; i++)
    {
        if (fabs(vals1[i]) < tol && fabs(vals2[i]) < tol)
            continue;
        if (fabs(vals1[i]) < tol || fabs(vals2[i]) < tol)
            return false;
        double computed = vals1[i] * ratio;
        if (fabs(computed - vals2[i]) > tol * fmax(fabs(computed), fabs(vals2[i])))
            return false;
    }
    
    *alpha = ratio;
    return true;
}

/* Merge collinear factors (columns of R = rows of RT) using hash-based
 * detection. Two factors l1, l2 are collinear if col(R, l1) = alpha * col(R, l2).
 * Their rank-1 contributions combine: col_l2 gets scaled by sqrt(1 + alpha^2)
 * and col_l1 is removed, reducing k.
 *
 * R is n×k CSR. We operate on RT (k×n) where row l of RT is column l of R. */
size_t merge_collinear_rows_in_R(QuadTermQR *qr, double tol)
{
    if (!qr || !qr->Rx || qr->Rnnz == 0 || qr->k <= 1) return 0;
    if (!qr->RTp || !qr->RTi || !qr->RTx) return 0;

    int n = (int)qr->n;
    int k = (int)qr->k;

    /* Hash buckets over factors (rows of RT). */
    typedef struct HashEntry {
        int row;
        struct HashEntry *next;
    } HashEntry;

    uint64_t num_buckets = (uint64_t)k * 2 + 1;
    HashEntry **buckets = (HashEntry **)ps_calloc((size_t)num_buckets, sizeof(HashEntry *));
    HashEntry *entries = (HashEntry *)ps_malloc((size_t)k, sizeof(HashEntry));

    if (!buckets || !entries)
    {
        PS_FREE(buckets); PS_FREE(entries);
        return 0;
    }

    for (int i = 0; i < k; i++) { entries[i].row = i; entries[i].next = NULL; }

    /* Hash each factor (RT row) by its normalized sparsity pattern + values. */
    for (int l = 0; l < k; l++)
    {
        int rt_start = qr->RTp[l];
        int rt_end   = qr->RTp[l + 1];
        int nnz      = rt_end - rt_start;
        if (nnz == 0) continue;

        double first_val = qr->RTx[rt_start];
        double scale = (fabs(first_val) > tol) ? 1.0 / fabs(first_val) : 1.0;

        uint64_t h = hash_sparse_row(&qr->RTx[rt_start], &qr->RTi[rt_start],
                                     nnz, scale, tol);
        int bucket = (int)(h % num_buckets);

        entries[l].next = buckets[bucket];
        buckets[bucket] = &entries[l];
    }

    int   *merge_with = (int *)   ps_malloc((size_t)k, sizeof(int));
    double *alpha_vals = (double *)ps_malloc((size_t)k, sizeof(double));
    if (!merge_with || !alpha_vals)
    {
        PS_FREE(buckets); PS_FREE(entries);
        PS_FREE(merge_with); PS_FREE(alpha_vals);
        return 0;
    }
    for (int i = 0; i < k; i++) { merge_with[i] = -1; alpha_vals[i] = 0.0; }

    for (uint64_t b = 0; b < num_buckets; b++)
    {
        HashEntry *e1 = buckets[b];
        while (e1)
        {
            int l1 = e1->row;
            if (merge_with[l1] >= 0) { e1 = e1->next; continue; }

            int s1 = qr->RTp[l1];
            int n1 = qr->RTp[l1 + 1] - s1;

            HashEntry *e2 = e1->next;
            while (e2)
            {
                int l2 = e2->row;
                if (merge_with[l2] >= 0) { e2 = e2->next; continue; }

                int s2 = qr->RTp[l2];
                int n2 = qr->RTp[l2 + 1] - s2;

                double alpha;
                if (rows_are_collinear(&qr->RTx[s1], &qr->RTi[s1], n1,
                                       &qr->RTx[s2], &qr->RTi[s2], n2,
                                       tol, &alpha))
                {
                    merge_with[l2] = l1;
                    alpha_vals[l2] = alpha;
                    break;
                }
                e2 = e2->next;
            }
            e1 = e1->next;
        }
    }

    /* Assign new factor indices to survivors. */
    int *new_factor_idx = (int *)ps_malloc((size_t)k, sizeof(int));
    if (!new_factor_idx)
    {
        PS_FREE(buckets); PS_FREE(entries);
        PS_FREE(merge_with); PS_FREE(alpha_vals);
        return 0;
    }

    int merged = 0;
    int survivor_idx = 0;
    for (int l = 0; l < k; l++)
    {
        if (merge_with[l] < 0) new_factor_idx[l] = survivor_idx++;
        else                   { new_factor_idx[l] = -1; merged++; }
    }

    if (merged == 0)
    {
        PS_FREE(buckets); PS_FREE(entries);
        PS_FREE(merge_with); PS_FREE(alpha_vals);
        PS_FREE(new_factor_idx);
        return 0;
    }

    int k_new = survivor_idx;

    /* Each survivor parent gets scaled by sqrt(1 + sum(alpha_i^2)) over its
     * merged children. Walk parent chains to their root. */
    double *factor_scales = (double *)ps_malloc((size_t)k, sizeof(double));
    if (!factor_scales)
    {
        PS_FREE(buckets); PS_FREE(entries);
        PS_FREE(merge_with); PS_FREE(alpha_vals);
        PS_FREE(new_factor_idx);
        return 0;
    }
    for (int i = 0; i < k; i++) factor_scales[i] = 1.0;

    for (int l = 0; l < k; l++)
    {
        if (merge_with[l] < 0) continue;
        int parent = merge_with[l];
        while (merge_with[parent] >= 0) parent = merge_with[parent];
        factor_scales[parent] += alpha_vals[l] * alpha_vals[l];
    }
    for (int l = 0; l < k; l++)
    {
        if (merge_with[l] < 0) factor_scales[l] = sqrt(factor_scales[l]);
    }

    /* Rebuild R: drop entries whose factor was merged away, and scale the
     * remaining entries by the scale of their (kept) parent factor. */
    size_t new_Rnnz = 0;
    for (size_t i = 0; i < qr->Rnnz; i++)
    {
        if (new_factor_idx[qr->Ri[i]] >= 0) new_Rnnz++;
    }

    double *new_Rx = (double *)ps_malloc(new_Rnnz > 0 ? new_Rnnz : 1, sizeof(double));
    int    *new_Ri = (int *)   ps_malloc(new_Rnnz > 0 ? new_Rnnz : 1, sizeof(int));
    int    *new_Rp = (int *)   ps_malloc((size_t)n + 1, sizeof(int));
    if (!new_Rx || !new_Ri || !new_Rp)
    {
        PS_FREE(buckets); PS_FREE(entries);
        PS_FREE(merge_with); PS_FREE(alpha_vals);
        PS_FREE(new_factor_idx); PS_FREE(factor_scales);
        PS_FREE(new_Rx); PS_FREE(new_Ri); PS_FREE(new_Rp);
        return 0;
    }

    new_Rp[0] = 0;
    size_t idx_new = 0;
    for (int v = 0; v < n; v++)
    {
        int rs = qr->Rp[v];
        int re = qr->Rp[v + 1];
        for (int idx = rs; idx < re; idx++)
        {
            int old_f = qr->Ri[idx];
            int new_f = new_factor_idx[old_f];
            if (new_f < 0) continue;
            new_Rx[idx_new] = qr->Rx[idx] * factor_scales[old_f];
            new_Ri[idx_new] = new_f;
            idx_new++;
        }
        new_Rp[v + 1] = (int)idx_new;
    }

    PS_FREE(qr->Rx); PS_FREE(qr->Ri); PS_FREE(qr->Rp);
    PS_FREE(qr->RTx); PS_FREE(qr->RTi); PS_FREE(qr->RTp);

    qr->Rx = new_Rx;
    qr->Ri = new_Ri;
    qr->Rp = new_Rp;
    qr->Rnnz = new_Rnnz;
    qr->k = (size_t)k_new;
    qr->RTx = NULL; qr->RTi = NULL; qr->RTp = NULL;

    if (new_Rnnz > 0)
        transpose_csr(qr->Rx, qr->Ri, qr->Rp, &qr->RTx, &qr->RTi, &qr->RTp,
                      (size_t)n, (size_t)k_new, new_Rnnz);

    PS_FREE(buckets); PS_FREE(entries);
    PS_FREE(merge_with); PS_FREE(alpha_vals);
    PS_FREE(new_factor_idx); PS_FREE(factor_scales);

    return (size_t)merged;
}

/* Check if variable k has any quadratic terms in QR decomposition (Q or R matrix).
 * Returns true if column k has non-zero entries in Q matrix or R matrix.
 */
bool has_quadratic_terms_qr(const QuadTermQR *quad_qr, int k)
{
    if (!quad_qr || !quad_qr->has_quad) return false;
    if (k < 0 || k >= quad_qr->n) return false;
    
    /* Check Q matrix - row k of Q (which is column k since Q is symmetric) */
    if (quad_qr->Qx != NULL && quad_qr->Qp != NULL)
    {
        int q_start = quad_qr->Qp[k];
        int q_end = quad_qr->Qp[k + 1];
        if (q_end > q_start) return true;  /* Has Q entries */
    }
    
    /* Check R matrix - row k of R (variable k's factor loadings).
     * R is n×k CSR so row pointer Rp has n+1 entries indexed by variable. */
    if (quad_qr->Rx != NULL && quad_qr->Rp != NULL)
    {
        int r_start = quad_qr->Rp[k];
        int r_end = quad_qr->Rp[k + 1];
        if (r_end > r_start) return true;  /* Has R entries */
    }

    return false;
}

/* Check if variable k has ONLY diagonal Q entry (no Q off-diagonal, no R entries).
 * This is useful for safe variable substitution in Doubleton Equality:
 * - Q[k][k] > 0 (diagonal entry exists)
 * - No other Q[k][j] or Q[j][k] entries (j != k)
 * - No R[:,k] entries (column k of R is empty)
 */
bool has_only_q_diag(const QuadTermQR *quad_qr, int k)
{
    if (!quad_qr || !quad_qr->has_quad) return false;
    if (k < 0 || k >= quad_qr->n) return false;
    
    /* Check if variable k has any R entries (row k of n×k R). */
    if (quad_qr->Rx != NULL && quad_qr->Rp != NULL)
    {
        int r_start = quad_qr->Rp[k];
        int r_end = quad_qr->Rp[k + 1];
        if (r_end > r_start) return false;  /* Has R entries */
    }
    
    /* Check Q matrix for column k */
    if (quad_qr->Qx != NULL && quad_qr->Qp != NULL)
    {
        /* Check if row k of Q (which is column k since Q is symmetric in CSR storage) 
         * has only the diagonal entry */
        int q_start = quad_qr->Qp[k];
        int q_end = quad_qr->Qp[k + 1];
        
        /* No Q entries at all for this variable */
        if (q_start == q_end) return false;
        
        /* Check if there's only the diagonal entry */
        if (q_end - q_start > 1) return false;  /* Has off-diagonal entries in row k */
        
        /* Verify the single entry is on the diagonal */
        if (quad_qr->Qi[q_start] != k) return false;  /* Not a diagonal entry */
        
        /* Also need to check if column k appears in other rows (j != k)
         * Since Q is stored in CSR, we need to scan all rows for column k */
        for (int row = 0; row < quad_qr->n; row++)
        {
            if (row == k) continue;  /* Skip diagonal row, already checked */
            
            int row_start = quad_qr->Qp[row];
            int row_end = quad_qr->Qp[row + 1];
            
            for (int idx = row_start; idx < row_end; idx++)
            {
                if (quad_qr->Qi[idx] == k)
                    return false;  /* Found off-diagonal reference to column k */
            }
        }
    }
    else
    {
        /* No Q matrix - can't have "only Q diagonal" if there's no Q at all
         * (unless we have no quadratic terms, but then has_quad would be false) */
        return false;
    }
    
    return true;
}

/* Compute bounds for (P*x)_k using QR representation (only for row k) */
bool compute_Px_bounds(const QuadTermQR *qr, int k,
                       const Bound *bounds, size_t n_cols,
                       double *min_Px, double *max_Px)
{
    if (!qr || !qr->has_quad || !bounds || !min_Px || !max_Px)
        return false;
    
    if (k < 0 || (size_t)k >= n_cols || n_cols != qr->n)
        return false;
    
    /* Initialize output */
    *min_Px = 0.0;
    *max_Px = 0.0;
    
    /* Add Q contribution for row k */
    if (qr->Qx != NULL && qr->Qp != NULL)
    {
        int q_start = qr->Qp[k];
        int q_end = qr->Qp[k + 1];
        for (int idx = q_start; idx < q_end; idx++)
        {
            int col = qr->Qi[idx];
            double qval = qr->Qx[idx];
            
            double lb = bounds[col].lb;
            double ub = bounds[col].ub;
            
            double contrib_lb = qval * ((qval >= 0) ? lb : ub);
            double contrib_ub = qval * ((qval >= 0) ? ub : lb);
            
            *min_Px += contrib_lb;
            *max_Px += contrib_ub;
        }
    }
    
    /* Add R*R^T contribution for row k */
    if (qr->Rx != NULL && qr->Rp != NULL && qr->RTp != NULL)
    {
        // int k_dim = (int)qr->k;
        
        /* For each non-zero R[k][l]: add R[k][l] * (R^T * x)_l */
        int r_start = qr->Rp[k];
        int r_end = qr->Rp[k + 1];
        
        for (int idx = r_start; idx < r_end; idx++)
        {
            int l = qr->Ri[idx];
            double r_val = qr->Rx[idx];
            
            /* Compute bounds for (R^T * x)_l = sum_i R[i][l] * x_i */
            double rtx_lb = 0.0;
            double rtx_ub = 0.0;
            
            int rt_start = qr->RTp[l];
            int rt_end = qr->RTp[l + 1];
            
            for (int rt_idx = rt_start; rt_idx < rt_end; rt_idx++)
            {
                int i = qr->RTi[rt_idx];
                double rt_val = qr->RTx[rt_idx];
                
                double lb = bounds[i].lb;
                double ub = bounds[i].ub;
                
                rtx_lb += rt_val * ((rt_val >= 0) ? lb : ub);
                rtx_ub += rt_val * ((rt_val >= 0) ? ub : lb);
            }
            
            /* Add contribution from R[k][l] * (R^T * x)_l */
            *min_Px += r_val * ((r_val >= 0) ? rtx_lb : rtx_ub);
            *max_Px += r_val * ((r_val >= 0) ? rtx_ub : rtx_lb);
        }
    }
    
    return true;
}
