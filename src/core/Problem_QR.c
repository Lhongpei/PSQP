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

    /* Copy R matrix */
    if (Rx != NULL && Rnnz > 0)
    {
        qr->Rnnz = Rnnz;
        qr->Rx = (double *) ps_malloc(Rnnz, sizeof(double));
        qr->Ri = (int *) ps_malloc(Rnnz, sizeof(int));
        /* R is k×n (k rows, n cols), CSR format: Rp has k+1 elements */
        qr->Rp = (int *) ps_malloc(k + 1, sizeof(int));
        
        if (!qr->Rx || !qr->Ri || !qr->Rp)
        {
            free_quad_term_qr(qr);
            return NULL;
        }
        
        memcpy(qr->Rx, Rx, Rnnz * sizeof(double));
        memcpy(qr->Ri, Ri, Rnnz * sizeof(int));
        memcpy(qr->Rp, Rp, (k + 1) * sizeof(int));
        
        /* Build R transpose for efficient access */
        /* R is k×n, so n_rows=k, n_cols=n */
        transpose_csr(Rx, Ri, Rp, &qr->RTx, &qr->RTi, &qr->RTp, k, n, Rnnz);
        if (!qr->RTx)
        {
            free_quad_term_qr(qr);
            return NULL;
        }
        
        /* Clean small entries */
        clean_zero_entries_qr(qr, 1e-13);
        
        /* Extract single-nnz rows to Q diagonal */
        extract_single_nnz_rows_from_R(qr);
        
        /* Merge collinear rows if k > 1 */
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
            new_Rnnz += (row_end - row_start);
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

/* Remove small entries below tolerance from R matrix */
void clean_zero_entries_qr(QuadTermQR *qr, double tol)
{
    if (!qr || !qr->Rx || qr->Rnnz == 0) return;
    if (tol <= 0) tol = 1e-13;
    
    int n = qr->n;
    int k = qr->k;
    
    /* First pass: count entries above tolerance */
    size_t new_nnz = 0;
    for (size_t i = 0; i < qr->Rnnz; i++)
    {
        if (fabs(qr->Rx[i]) > tol)
            new_nnz++;
    }
    
    if (new_nnz == qr->Rnnz) return;  /* Nothing to remove */
    
    /* Allocate new arrays */
    double *new_Rx = (double *)ps_malloc(new_nnz, sizeof(double));
    int *new_Ri = (int *)ps_malloc(new_nnz, sizeof(int));
    int *new_Rp = (int *)ps_malloc(k + 1, sizeof(int));
    
    if (!new_Rx || !new_Ri || !new_Rp)
    {
        PS_FREE(new_Rx);
        PS_FREE(new_Ri);
        PS_FREE(new_Rp);
        return;
    }
    
    /* Second pass: copy surviving entries */
    new_Rp[0] = 0;
    size_t idx_new = 0;
    for (int row = 0; row < k; row++)
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
    
    /* Rebuild transpose */
    if (new_nnz > 0)
        transpose_csr(qr->Rx, qr->Ri, qr->Rp, &qr->RTx, &qr->RTi, &qr->RTp, k, n, new_nnz);
}

/* Extract rows of R that have only a single non-zero entry and move them to Q diagonal.
 * Returns number of rows extracted.
 */
size_t extract_single_nnz_rows_from_R(QuadTermQR *qr)
{
    if (!qr || !qr->Rx || qr->Rnnz == 0 || qr->k <= 0) return 0;
    
    int n = qr->n;
    int k = qr->k;
    
    /* Count rows with single nnz */
    int single_nnz_rows = 0;
    for (int row = 0; row < k; row++)
    {
        int row_start = qr->Rp[row];
        int row_end = qr->Rp[row + 1];
        if (row_end - row_start == 1)
            single_nnz_rows++;
    }
    
    if (single_nnz_rows == 0) return 0;
    
    /* Calculate additional Q storage needed */
    size_t additional_q_nnz = single_nnz_rows;
    size_t new_Qnnz = qr->Qnnz + additional_q_nnz;
    
    /* Allocate new Q arrays */
    double *new_Qx = (double *)ps_malloc(new_Qnnz, sizeof(double));
    int *new_Qi = (int *)ps_malloc(new_Qnnz, sizeof(int));
    int *new_Qp = (int *)ps_malloc(n + 1, sizeof(int));
    
    if (!new_Qx || !new_Qi || !new_Qp)
    {
        PS_FREE(new_Qx);
        PS_FREE(new_Qi);
        PS_FREE(new_Qp);
        return 0;
    }
    
    /* Build row map: which rows to keep in R */
    int *row_keep = (int *)ps_malloc(k, sizeof(int));
    if (!row_keep)
    {
        PS_FREE(new_Qx);
        PS_FREE(new_Qi);
        PS_FREE(new_Qp);
        return 0;
    }
    
    int new_row_idx = 0;
    for (int row = 0; row < k; row++)
    {
        int row_start = qr->Rp[row];
        int row_end = qr->Rp[row + 1];
        if (row_end - row_start == 1)
            row_keep[row] = -1;  /* Mark for extraction */
        else
            row_keep[row] = new_row_idx++;
    }
    
    int k_new = new_row_idx;
    
    /* Copy existing Q entries */
    memcpy(new_Qx, qr->Qx, qr->Qnnz * sizeof(double));
    memcpy(new_Qi, qr->Qi, qr->Qnnz * sizeof(int));
    memcpy(new_Qp, qr->Qp, (n + 1) * sizeof(int));
    
    /* Add extracted rows to Q diagonal */
    size_t q_idx = qr->Qnnz;
    for (int row = 0; row < k; row++)
    {
        if (row_keep[row] < 0)
        {
            /* Extract this row */
            int idx = qr->Rp[row];
            double val = qr->Rx[idx];
            int col = qr->Ri[idx];
            
            /* Add val^2 to Q[col,col] */
            new_Qx[q_idx] = val * val;
            new_Qi[q_idx] = col;
            q_idx++;
        }
    }
    
    /* Sort Q entries by row (column index in this transpose representation) */
    /* Simple insertion sort for small arrays */
    for (size_t i = qr->Qnnz + 1; i < new_Qnnz; i++)
    {
        double val = new_Qx[i];
        int col = new_Qi[i];
        size_t j = i;
        while (j > qr->Qnnz && new_Qi[j-1] > col)
        {
            new_Qx[j] = new_Qx[j-1];
            new_Qi[j] = new_Qi[j-1];
            j--;
        }
        new_Qx[j] = val;
        new_Qi[j] = col;
    }
    
    /* Rebuild Qp */
    /* First count entries per row */
    int *row_counts = (int *)ps_calloc(n, sizeof(int));
    for (size_t i = 0; i < new_Qnnz; i++)
        row_counts[new_Qi[i]]++;
    
    new_Qp[0] = 0;
    for (int i = 0; i < n; i++)
        new_Qp[i+1] = new_Qp[i] + row_counts[i];
    
    PS_FREE(row_counts);
    
    /* Free old Q and assign new */
    PS_FREE(qr->Qx);
    PS_FREE(qr->Qi);
    PS_FREE(qr->Qp);
    qr->Qx = new_Qx;
    qr->Qi = new_Qi;
    qr->Qp = new_Qp;
    qr->Qnnz = new_Qnnz;
    
    /* Build new R with remaining rows */
    if (k_new == 0)
    {
        /* All rows extracted */
        PS_FREE(qr->Rx);
        PS_FREE(qr->Ri);
        PS_FREE(qr->Rp);
        PS_FREE(qr->RTx);
        PS_FREE(qr->RTi);
        PS_FREE(qr->RTp);
        qr->Rx = NULL; qr->Ri = NULL; qr->Rp = NULL;
        qr->RTx = NULL; qr->RTi = NULL; qr->RTp = NULL;
        qr->Rnnz = 0;
        qr->k = 0;
    }
    else
    {
        /* Count surviving R entries */
        size_t new_Rnnz = 0;
        for (int row = 0; row < k; row++)
        {
            if (row_keep[row] >= 0)
            {
                int row_start = qr->Rp[row];
                int row_end = qr->Rp[row + 1];
                new_Rnnz += (row_end - row_start);
            }
        }
        
        double *new_Rx = (double *)ps_malloc(new_Rnnz, sizeof(double));
        int *new_Ri = (int *)ps_malloc(new_Rnnz, sizeof(int));
        int *new_Rp = (int *)ps_malloc(k_new + 1, sizeof(int));
        
        if (!new_Rx || !new_Ri || !new_Rp)
        {
            PS_FREE(new_Rx);
            PS_FREE(new_Ri);
            PS_FREE(new_Rp);
            PS_FREE(row_keep);
            return 0;
        }
        
        new_Rp[0] = 0;
        size_t idx_new = 0;
        for (int row = 0; row < k; row++)
        {
            if (row_keep[row] >= 0)
            {
                int row_start = qr->Rp[row];
                int row_end = qr->Rp[row + 1];
                for (int idx = row_start; idx < row_end; idx++)
                {
                    new_Rx[idx_new] = qr->Rx[idx];
                    new_Ri[idx_new] = qr->Ri[idx];
                    idx_new++;
                }
                new_Rp[row_keep[row] + 1] = (int)idx_new;
            }
        }
        
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
        qr->k = k_new;
        qr->RTx = NULL;
        qr->RTi = NULL;
        qr->RTp = NULL;
        
        /* Rebuild transpose */
        if (new_Rnnz > 0)
            transpose_csr(qr->Rx, qr->Ri, qr->Rp, &qr->RTx, &qr->RTi, &qr->RTp, k_new, n, new_Rnnz);
    }
    
    PS_FREE(row_keep);
    return (size_t)single_nnz_rows;
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

/* Merge collinear rows in R using hash-based detection */
size_t merge_collinear_rows_in_R(QuadTermQR *qr, double tol)
{
    if (!qr || !qr->Rx || qr->Rnnz == 0 || qr->k <= 1) return 0;
    
    int n = qr->n;
    int k = qr->k;
    
    /* Build hash buckets */
    typedef struct HashEntry {
        int row;
        struct HashEntry *next;
    } HashEntry;
    
    int num_buckets = k * 2 + 1;
    HashEntry **buckets = (HashEntry **)ps_calloc(num_buckets, sizeof(HashEntry *));
    HashEntry *entries = (HashEntry *)ps_malloc(k, sizeof(HashEntry));
    
    if (!buckets || !entries)
    {
        PS_FREE(buckets);
        PS_FREE(entries);
        return 0;
    }
    
    for (int i = 0; i < k; i++)
    {
        entries[i].row = i;
        entries[i].next = NULL;
    }
    
    /* Hash each row and add to bucket */
    for (int row = 0; row < k; row++)
    {
        int row_start = qr->Rp[row];
        int row_end = qr->Rp[row + 1];
        int nnz = row_end - row_start;
        
        if (nnz == 0) continue;
        
        /* Normalize row by dividing by first element magnitude */
        double first_val = qr->Rx[row_start];
        double scale = (fabs(first_val) > tol) ? 1.0 / fabs(first_val) : 1.0;
        
        uint64_t h = hash_sparse_row(&qr->Rx[row_start], &qr->Ri[row_start], nnz, scale, tol);
        int bucket = (int)(h % num_buckets);
        
        entries[row].next = buckets[bucket];
        buckets[bucket] = &entries[row];
    }
    
    /* Find collinear pairs */
    int *merge_with = (int *)ps_malloc(k, sizeof(int));
    double *alpha_vals = (double *)ps_malloc(k, sizeof(double));
    
    if (!merge_with || !alpha_vals)
    {
        PS_FREE(buckets);
        PS_FREE(entries);
        PS_FREE(merge_with);
        PS_FREE(alpha_vals);
        return 0;
    }
    
    for (int i = 0; i < k; i++)
    {
        merge_with[i] = -1;
        alpha_vals[i] = 0.0;
    }
    
    /* For each bucket, check all pairs */
    for (int b = 0; b < num_buckets; b++)
    {
        HashEntry *e1 = buckets[b];
        while (e1)
        {
            int row1 = e1->row;
            if (merge_with[row1] >= 0)  /* Already merged */
            {
                e1 = e1->next;
                continue;
            }
            
            int r1_start = qr->Rp[row1];
            int r1_end = qr->Rp[row1 + 1];
            int nnz1 = r1_end - r1_start;
            
            HashEntry *e2 = e1->next;
            while (e2)
            {
                int row2 = e2->row;
                if (merge_with[row2] >= 0)  /* Already merged */
                {
                    e2 = e2->next;
                    continue;
                }
                
                int r2_start = qr->Rp[row2];
                int r2_end = qr->Rp[row2 + 1];
                int nnz2 = r2_end - r2_start;
                
                double alpha;
                if (rows_are_collinear(&qr->Rx[r1_start], &qr->Ri[r1_start], nnz1,
                                       &qr->Rx[r2_start], &qr->Ri[r2_start], nnz2,
                                       tol, &alpha))
                {
                    merge_with[row2] = row1;
                    alpha_vals[row2] = alpha;
                    break;
                }
                e2 = e2->next;
            }
            e1 = e1->next;
        }
    }
    
    /* Count merged rows and survivors */
    int merged = 0;
    int *new_row_idx = (int *)ps_malloc(k, sizeof(int));
    int survivor_idx = 0;
    
    for (int i = 0; i < k; i++)
    {
        if (merge_with[i] < 0)
            new_row_idx[i] = survivor_idx++;
        else
        {
            new_row_idx[i] = -1;
            merged++;
        }
    }
    
    if (merged == 0)
    {
        PS_FREE(buckets);
        PS_FREE(entries);
        PS_FREE(merge_with);
        PS_FREE(alpha_vals);
        PS_FREE(new_row_idx);
        return 0;
    }
    
    int k_new = survivor_idx;
    
    /* First pass: scale surviving rows and count new nnz */
    /* For row i that survives and has merged rows: new_R_i = sqrt(1 + sum alpha^2) * R_i */
    double *row_scales = (double *)ps_malloc(k, sizeof(double));
    for (int i = 0; i < k; i++)
        row_scales[i] = 1.0;
    
    /* Accumulate scale factors */
    for (int i = 0; i < k; i++)
    {
        if (merge_with[i] >= 0)
        {
            int parent = merge_with[i];
            while (merge_with[parent] >= 0)
                parent = merge_with[parent];
            row_scales[parent] += alpha_vals[i] * alpha_vals[i];
        }
    }
    
    /* Apply sqrt to get final scales */
    for (int i = 0; i < k; i++)
    {
        if (merge_with[i] < 0)
            row_scales[i] = sqrt(row_scales[i]);
    }
    
    /* Count new nnz (only surviving rows) */
    size_t new_Rnnz = 0;
    for (int row = 0; row < k; row++)
    {
        if (new_row_idx[row] >= 0)
        {
            int row_start = qr->Rp[row];
            int row_end = qr->Rp[row + 1];
            new_Rnnz += (row_end - row_start);
        }
    }
    
    /* Allocate new R arrays */
    double *new_Rx = (double *)ps_malloc(new_Rnnz, sizeof(double));
    int *new_Ri = (int *)ps_malloc(new_Rnnz, sizeof(int));
    int *new_Rp = (int *)ps_malloc(k_new + 1, sizeof(int));
    
    if (!new_Rx || !new_Ri || !new_Rp)
    {
        PS_FREE(buckets);
        PS_FREE(entries);
        PS_FREE(merge_with);
        PS_FREE(alpha_vals);
        PS_FREE(new_row_idx);
        PS_FREE(row_scales);
        PS_FREE(new_Rx);
        PS_FREE(new_Ri);
        PS_FREE(new_Rp);
        return 0;
    }
    
    /* Fill new R with scaled rows */
    new_Rp[0] = 0;
    size_t idx_new = 0;
    for (int row = 0; row < k; row++)
    {
        int new_idx = new_row_idx[row];
        if (new_idx < 0) continue;
        
        int row_start = qr->Rp[row];
        int row_end = qr->Rp[row + 1];
        double scale = row_scales[row];
        
        for (int idx = row_start; idx < row_end; idx++)
        {
            new_Rx[idx_new] = qr->Rx[idx] * scale;
            new_Ri[idx_new] = qr->Ri[idx];
            idx_new++;
        }
        new_Rp[new_idx + 1] = (int)idx_new;
    }
    
    /* Free old R and assign new */
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
    qr->k = k_new;
    qr->RTx = NULL;
    qr->RTi = NULL;
    qr->RTp = NULL;
    
    /* Rebuild transpose */
    if (new_Rnnz > 0)
        transpose_csr(qr->Rx, qr->Ri, qr->Rp, &qr->RTx, &qr->RTi, &qr->RTp, k_new, n, new_Rnnz);
    
    /* Cleanup */
    PS_FREE(buckets);
    PS_FREE(entries);
    PS_FREE(merge_with);
    PS_FREE(alpha_vals);
    PS_FREE(new_row_idx);
    PS_FREE(row_scales);
    
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
    
    /* Check R matrix - column k of R (row k of R transpose) */
    if (quad_qr->Rx != NULL && quad_qr->RTp != NULL)
    {
        int rt_start = quad_qr->RTp[k];
        int rt_end = quad_qr->RTp[k + 1];
        if (rt_end > rt_start) return true;  /* Has R entries */
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
    
    /* Check if column k of R has any entries (using R transpose) */
    if (quad_qr->Rx != NULL && quad_qr->RTp != NULL)
    {
        /* RT is k×n, row k of RT is column k of R */
        int rt_start = quad_qr->RTp[k];
        int rt_end = quad_qr->RTp[k + 1];
        if (rt_end > rt_start) return false;  /* Has R entries */
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

/* Compute bounds for P*x using QR representation */
bool compute_Px_bounds(const QuadTermQR *qr, int k,
                       const Bound *bounds, size_t n_cols,
                       double *min_Px, double *max_Px)
{
    (void)k;  /* Unused, kept for API compatibility */
    if (!qr || !qr->has_quad || !bounds || !min_Px || !max_Px)
        return false;
    
    int n = (int)n_cols;
    if (n != qr->n) return false;
    
    /* Extract lb/ub arrays from bounds */
    double *x_lb = (double *)ps_malloc(n, sizeof(double));
    double *x_ub = (double *)ps_malloc(n, sizeof(double));
    if (!x_lb || !x_ub)
    {
        PS_FREE(x_lb);
        PS_FREE(x_ub);
        return false;
    }
    
    for (int i = 0; i < n; i++)
    {
        x_lb[i] = bounds[i].lb;
        x_ub[i] = bounds[i].ub;
    }
    
    double *px_lb = min_Px;
    double *px_ub = max_Px;
    
    /* Initialize bounds */
    for (int i = 0; i < n; i++)
    {
        px_lb[i] = 0.0;
        px_ub[i] = 0.0;
    }
    
    /* Add Q contribution (diagonal) */
    if (qr->Qx != NULL && qr->Qp != NULL)
    {
        for (int row = 0; row < n; row++)
        {
            int q_start = qr->Qp[row];
            int q_end = qr->Qp[row + 1];
            for (int idx = q_start; idx < q_end; idx++)
            {
                int col = qr->Qi[idx];
                double qval = qr->Qx[idx];
                
                /* Contribution to row 'row' from column 'col' */
                double lb = qval * ((qval >= 0) ? x_lb[col] : x_ub[col]);
                double ub = qval * ((qval >= 0) ? x_ub[col] : x_lb[col]);
                
                px_lb[row] += lb;
                px_ub[row] += ub;
                
                /* Symmetric contribution if off-diagonal */
                if (col != row && qr->Qp != NULL)
                {
                    /* Add to column's row as well (symmetric Q) */
                    /* Actually Q is symmetric in QP, but stored in CSR format */
                    /* This is tricky - we assume Q is diagonal for now */
                }
            }
        }
    }
    
    /* Add R*R^T contribution */
    if (qr->Rx != NULL && qr->Rp != NULL && qr->RTp != NULL)
    {
        int k_dim = qr->k;
        
        /* For each row i: (R*R^T)_{i,:} * x = sum_l R_{i,l} * (R^T * x)_l */
        /* First compute R^T * x bounds for each column l of R (row l of R^T) */
        double *rtx_lb = (double *)ps_malloc(k_dim, sizeof(double));
        double *rtx_ub = (double *)ps_malloc(k_dim, sizeof(double));
        
        if (!rtx_lb || !rtx_ub)
        {
            PS_FREE(x_lb);
            PS_FREE(x_ub);
            PS_FREE(rtx_lb);
            PS_FREE(rtx_ub);
            return false;
        }
        
        /* Compute (R^T * x)_l bounds: sum_i R_{i,l} * x_i */
        for (int l = 0; l < k_dim; l++)
        {
            rtx_lb[l] = 0.0;
            rtx_ub[l] = 0.0;
            
            int rt_start = qr->RTp[l];
            int rt_end = qr->RTp[l + 1];
            
            for (int idx = rt_start; idx < rt_end; idx++)
            {
                int i = qr->RTi[idx];
                double r_val = qr->RTx[idx];
                
                double lb = r_val * ((r_val >= 0) ? x_lb[i] : x_ub[i]);
                double ub = r_val * ((r_val >= 0) ? x_ub[i] : x_lb[i]);
                
                rtx_lb[l] += lb;
                rtx_ub[l] += ub;
            }
        }
        
        /* Now compute R * (R^T * x) bounds */
        for (int row = 0; row < n; row++)
        {
            int r_start = qr->Rp[row];
            int r_end = qr->Rp[row + 1];
            
            for (int idx = r_start; idx < r_end; idx++)
            {
                int l = qr->Ri[idx];
                double r_val = qr->Rx[idx];
                
                double lb = r_val * ((r_val >= 0) ? rtx_lb[l] : rtx_ub[l]);
                double ub = r_val * ((r_val >= 0) ? rtx_ub[l] : rtx_lb[l]);
                
                px_lb[row] += lb;
                px_ub[row] += ub;
            }
        }
        
        PS_FREE(rtx_lb);
        PS_FREE(rtx_ub);
    }
    
    PS_FREE(x_lb);
    PS_FREE(x_ub);
    return true;
}
