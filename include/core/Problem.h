/*
 * Copyright 2025-2026 Daniel Cederberg
 *
 * This file is part of the PSLP project (LP Presolver).
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

#ifndef CORE_PROBLEM_H
#define CORE_PROBLEM_H

#include <stdbool.h>
#include <stddef.h>
#include "PSQP_status.h"
struct Constraints;

/* ============================================================================
 * QuadTermQR: P = Q + R*R^T decomposition
 * 
 * This representation is useful when P has a low-rank component R*R^T plus
 * a sparse component Q (usually diagonal).
 * 
 * Q: sparse symmetric matrix (usually diagonal) in CSR upper triangular format
 * R: sparse n x k matrix in CSR format
 * ============================================================================ */

typedef struct
{
    /* Q matrix: sparse symmetric, stored upper triangular */
    double *Qx;     /* values */
    int *Qi;        /* column indices */
    int *Qp;        /* row pointers */
    size_t Qnnz;    /* number of non-zeros */
    size_t n;       /* dimension (n x n) */
    
    /* R matrix: sparse n x k */
    double *Rx;     /* values */
    int *Ri;        /* column indices */
    int *Rp;        /* row pointers */
    size_t Rnnz;    /* number of non-zeros */
    size_t k;       /* number of columns in R */
    
    /* R transpose: k x n, for efficient column access of R */
    double *RTx;
    int *RTi;
    int *RTp;
    
    bool has_quad;  /* true if Q or R is non-empty */
} QuadTermQR;

typedef struct Objective
{
    double *c;
    double offset;
    QuadTermQR *quad_qr; /* quadratic term in QR format (P = Q + R*R^T), NULL for LP */
} Objective;

typedef struct Problem
{
    struct Constraints *constraints;
    Objective *obj;
} Problem;

/* Constructor and destructor */
Problem *new_problem(struct Constraints *constraints, Objective *obj);
void free_problem(Problem *problem);

/* Cleans the problem in the sense that all data associated with
   inactive rows and inactive columns is removed. Rows and columns
   are also re-indexed to take the removal of rows/columns into
   account. */
void problem_clean(Problem *problem, bool remove_all);

/* Constructor and destructor */
Objective *objective_new_qr(double *c, QuadTermQR *quad_qr);
void free_objective(Objective *obj);

/* Updates the offset when a variable is fixed */
void fix_var_in_obj(Objective *obj, int col, double value);

/* Updates the offset and P matrix when a variable is fixed (QP version) */
void fix_var_in_obj_qp(Objective *obj, int col, double value);

/* Substitutes variable 'k' from the objective using row 'i'
   (assumed to be an equality). Row 'i' is specified by
   (vals, cols, len, rhs) */
void sub_var_in_obj(Objective *obj, const double *vals, const int *cols, int len,
                    int k, double aik, double rhs);

/* Substitutes variable 'k' from the QP objective using row 'i' (QP version) */
void sub_var_in_obj_qp(Objective *obj, const double *vals, const int *cols, int len,
                       int k, double aik, double rhs);

/* Constructor and destructor for QR decomposition */
QuadTermQR *quad_term_qr_new(const double *Qx, const int *Qi, const int *Qp, size_t Qnnz,
                             const double *Rx, const int *Ri, const int *Rp, size_t Rnnz,
                             size_t n, size_t k);
void free_quad_term_qr(QuadTermQR *quad_qr);

/* Compute P_ii = Q_ii + sum_j R_ij^2 (diagonal element of P = Q + R*R^T) */
double compute_p_diag(const QuadTermQR *quad_qr, int col);

/* Compute contribution to linear term from fixed variable:
 * For fixed x_i, adds to c[j]: (Q_ij + sum_l R_il*R_jl) * x_i
 * This handles both Q and R*R^T contributions */
void add_fixed_var_contribution_qr(QuadTermQR *quad_qr, int fixed_col, double value, double *c);

/* Updates offset when variable is fixed (QR version):
 * offset += 0.5 * (Q_ii + sum_j R_ij^2) * x_i^2 */
void fix_var_in_obj_qp_qr(QuadTermQR *quad_qr, int col, double value, double *offset, double *c);

/* Shrink Q and R when variables are removed */
void quad_term_qr_shrink(QuadTermQR *quad_qr, int *col_map, size_t n_cols_old);

/* Check if variable has quadratic terms in R */
bool has_r_terms(const QuadTermQR *quad_qr, int col);

/* Check if variable k has any quadratic terms in QR decomposition (Q or R matrix).
 * Returns true if column k has non-zero entries in Q matrix or R matrix.
 */
bool has_quadratic_terms_qr(const QuadTermQR *quad_qr, int k);

/* Check if variable k has ONLY diagonal entry in Q (no Q off-diagonal, no R entries).
 * This is useful for safe variable substitution in Doubleton Equality:
 * - Q[k][k] > 0 (diagonal entry exists)
 * - No other Q[k][j] or Q[j][k] entries (j != k)
 * - No R[:,k] entries (column k of R is empty)
 * 
 * Returns true if variable is "pure Q diagonal" and can be safely substituted
 * without creating cross-terms in the objective.
 */
bool has_only_q_diag(const QuadTermQR *quad_qr, int k);

/* Compute bounds on (P·x)_k for dual fix in QP.
 * Returns true if bounds can be computed, false if infinite bounds make it impossible.
 */
/* Compute bounds on (P·x)_k for dual fix in QP.
 * Returns true if bounds can be computed, false if infinite bounds make it impossible.
 */
struct Bound;  /* Forward declaration */
bool compute_Px_bounds(const QuadTermQR *qr, int k,
                       const struct Bound *bounds, size_t n_cols,
                       double *min_Px, double *max_Px);

/* Clean small entries (noise) from R matrix using a threshold.
 * Entries with |value| < tol are treated as zeros and removed.
 * Typical tolerance: 1e-13 to 1e-15
 */
void clean_zero_entries_qr(QuadTermQR *qr, double tol);

/* Extract rows of R that have only a single non-zero entry and move them to Q diagonal.
 * This reduces the rank of R by extracting "diagonal-like" terms.
 * Returns number of rows extracted.
 */
size_t extract_single_nnz_rows_from_R(QuadTermQR *qr);

/* Detect and merge collinear rows in R matrix.
 * 
 * If rows r_i and r_j are collinear (r_j = alpha * r_i), we can merge them
 * to reduce the rank of R. The merged row becomes:
 *   r_new = sqrt(alpha^2 + 1) * r_i
 * 
 * Returns number of rows merged.
 */
size_t merge_collinear_rows_in_R(QuadTermQR *qr, double tol);

/* Forward declaration for Problem struct */
struct Problem;

/* Detect and fix isolated quadratic variables to their optimal values.
 * 
 * A variable is isolated if:
 * 1. It has a diagonal entry in Q (Q[i][i] > 0)
 * 2. It has no entries in R
 * 3. It has no entries in constraint matrix A
 * 
 * For such variable x_i, the problem is:
 *   min 0.5 * Q[i][i] * x_i^2 + c_i * x_i
 *   s.t. lb_i <= x_i <= ub_i
 * 
 * Optimal solution: x_i* = -c_i / Q[i][i] (clipped to bounds)
 * 
 * Returns REDUCED if any variables were fixed, UNCHANGED otherwise.
 * May return INFEASIBLE if the computed optimal value conflicts with bounds.
 */
PresolveStatus fix_isolated_quadratic_vars_qp(struct Problem *prob);

#endif // CORE_PROBLEM_H
