/*
 * Copyright 2025-2026 Daniel Cederberg
 * Copyright 2026 [Your Name]
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
        qr->Rp = (int *) ps_malloc(n + 1, sizeof(int));
        
        if (!qr->Rx || !qr->Ri || !qr->Rp)
        {
            free_quad_term_qr(qr);
            return NULL;
        }
        
        memcpy(qr->Rx, Rx, Rnnz * sizeof(double));
        memcpy(qr->Ri, Ri, Rnnz * sizeof(int));
        memcpy(qr->Rp, Rp, (n + 1) * sizeof(int));
        
        /* Build R transpose for efficient access */
        transpose_csr(Rx, Ri, Rp, &qr->RTx, &qr->RTi, &qr->RTp, n, k, Rnnz);
        if (!qr->RTx)
        {
            free_quad_term_qr(qr);
            return NULL;
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
