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

#include <math.h>  /* For isnan, NAN */

#include "Simple_dual_fix.h"
#include "Bounds.h"
#include "Constraints.h"
#include "CoreTransformations.h"
#include "Debugger.h"
#include "Locks.h"
#include "Problem.h"
#include "SimpleReductions.h"
#include "State.h"
#include "Workspace.h"
#include "iVec.h"

// a column that becomes fixed because of dual fix should have zk = ck - ak^T y
static inline PresolveStatus _simple_dual_fix(Constraints *constraints, double ck,
                                              double lb, double ub, Lock lock,
                                              ColTag col_tag, int col,
                                              iVec *cols_to_inf, Objective *obj)
{
    assert(!HAS_TAG(col_tag, C_TAG_INACTIVE));

    // --------------------------------------------------
    // see if we can fix the variable to its lower bound
    // --------------------------------------------------
    if (ck > 0 && lock.down == 0)
    {
        if (HAS_TAG(col_tag, C_TAG_LB_INF))
        {
            return UNBNDORINFEAS;
        }

        assert(!IS_ABS_INF(lb));
        fix_col(constraints, col, lb, ck, obj);
        return UNCHANGED;
    }

    // ---------------------------------------------------
    // see if we can fix the variable to its upper bound
    // ---------------------------------------------------
    if (ck < 0 && lock.up == 0)
    {
        if (ck != 0 && HAS_TAG(col_tag, C_TAG_UB_INF))
        {
            return UNBNDORINFEAS;
        }

        assert(!IS_ABS_INF(ub));
        fix_col(constraints, col, ub, ck, obj);
        return UNCHANGED;
    }

    // ------------------------------------------------------------------------
    // If the objective coefficient is 0 and there are no downlocks,
    // then the variable may be set to -infinity if it is unbounded from below,
    // making the constraints it appears in redundant. If the variable is
    // bounded from below we can set the variable to any value. Our design
    // choice is to set the variable to its lower bound.
    // A similar comment applies to the case when there are no uplocks.
    // ------------------------------------------------------------------------
    if (ck == 0)
    {
        if (lock.down == 0)
        {
            if (HAS_TAG(col_tag, C_TAG_LB_INF))
            {
                iVec_append(cols_to_inf, -col);
            }
            else
            {
                assert(!IS_ABS_INF(lb));
                fix_col(constraints, col, lb, ck, obj);
            }
            return UNCHANGED;
        }

        if (lock.up == 0)
        {
            if (HAS_TAG(col_tag, C_TAG_UB_INF))
            {
                iVec_append(cols_to_inf, col);
            }
            else
            {
                assert(!IS_ABS_INF(ub));
                fix_col(constraints, col, ub, ck, obj);
            }
            return UNCHANGED;
        }
    }

    return UNCHANGED;
}

/* Forward declarations - these are in Problem.h */
extern bool has_quadratic_terms_qr(const QuadTermQR *quad_qr, int k);
extern bool compute_Px_bounds(const QuadTermQR *qr, int k, 
                               const Bound *bounds, size_t n_cols,
                               double *min_Px, double *max_Px);

/* Helper function for QP dual fix with quadratic terms.
 * Returns the effective coefficient for dual fix decision.
 * For LP: returns c_k
 * For QP with known (P·x)_k bounds: returns c_k + adjusted_Px
 * For QP with unbounded (P·x)_k: returns NaN (cannot determine)
 */
static inline double get_effective_coeff_for_dual_fix(const QuadTermQR *qr, int k,
                                                       double c_k, 
                                                       const Bound *bounds, 
                                                       size_t n_cols,
                                                       bool *can_fix_to_lb,
                                                       bool *can_fix_to_ub)
{
    *can_fix_to_lb = false;
    *can_fix_to_ub = false;
    
    if (qr == NULL || !qr->has_quad || !has_quadratic_terms_qr(qr, k))
    {
        /* Linear variable: use standard LP dual fix */
        *can_fix_to_lb = (c_k > 0);
        *can_fix_to_ub = (c_k < 0);
        return c_k;
    }
    
    /* QP variable with quadratic terms: compute bounds on (P·x)_k */
    double min_Px, max_Px;
    if (!compute_Px_bounds(qr, k, bounds, n_cols, &min_Px, &max_Px))
    {
        /* Cannot determine bounds (infinite bounds in involved variables) */
        return NAN;
    }
    
    /* For dual fix to lb: need c_k + (P·x)_k > 0 for all feasible x
     * The minimum of c_k + (P·x)_k is c_k + min_Px
     */
    *can_fix_to_lb = (c_k + min_Px > 0);
    
    /* For dual fix to ub: need c_k + (P·x)_k < 0 for all feasible x
     * The maximum of c_k + (P·x)_k is c_k + max_Px
     */
    *can_fix_to_ub = (c_k + max_Px < 0);
    
    return c_k;  /* Return original c_k for additional checks */
}

PresolveStatus simple_dual_fix(Problem *prob)
{
    assert(prob->constraints->state->empty_cols->len == 0);
    DEBUG(verify_problem_up_to_date(prob->constraints));

    Constraints *constraints = prob->constraints;
    const double *c = prob->obj->c;
    const Bound *bounds = constraints->bounds;
    const ColTag *col_tags = constraints->col_tags;
    const Lock *locks = constraints->state->col_locks;
    size_t n_cols = constraints->n;
    iVec *cols_to_inf = constraints->state->work->int_vec;
    iVec_clear_no_resize(cols_to_inf);
    PresolveStatus status = UNCHANGED;

    for (size_t k = 0; k < n_cols; k++)
    {
        // if it is too slow to loop through all columns, we can keep a list
        // of candidate columns to process
        if (HAS_TAG(col_tags[k], C_TAG_INACTIVE) ||
            (locks[k].up > 0 && locks[k].down > 0))
        {
            continue;
        }

        /* For QP: compute effective coefficient including quadratic term contributions */
        bool can_fix_to_lb, can_fix_to_ub;
        double effective_c = get_effective_coeff_for_dual_fix(
            prob->obj->quad_qr, (int) k, c[k], bounds, n_cols, 
            &can_fix_to_lb, &can_fix_to_ub);
        
        /* If we cannot determine bounds (NaN), skip this variable */
        if (isnan(effective_c))
        {
            continue;
        }
        
        /* For LP or QP with determined bounds: use the can_fix_to_lb/ub flags */
        
        // --------------------------------------------------
        // see if we can fix the variable to its lower bound
        // --------------------------------------------------
        if (can_fix_to_lb && locks[k].down == 0)
        {
            if (HAS_TAG(col_tags[k], C_TAG_LB_INF))
            {
                return UNBNDORINFEAS;
            }

            assert(!IS_ABS_INF(bounds[k].lb));
            fix_col(constraints, (int) k, bounds[k].lb, c[k], prob->obj);
            continue;  /* Handled, move to next variable */
        }

        // ---------------------------------------------------
        // see if we can fix the variable to its upper bound
        // ---------------------------------------------------
        if (can_fix_to_ub && locks[k].up == 0)
        {
            if (HAS_TAG(col_tags[k], C_TAG_UB_INF))
            {
                return UNBNDORINFEAS;
            }

            assert(!IS_ABS_INF(bounds[k].ub));
            fix_col(constraints, (int) k, bounds[k].ub, c[k], prob->obj);
            continue;  /* Handled, move to next variable */
        }

        // ------------------------------------------------------------------------
        // If the objective coefficient is 0 and there are no downlocks,
        // then the variable may be set to -infinity if it is unbounded from below,
        // making the constraints it appears in redundant. If the variable is
        // bounded from below we can set the variable to any value. Our design
        // choice is to set the variable to its lower bound.
        // A similar comment applies to the case when there are no uplocks.
        // ------------------------------------------------------------------------
        if (effective_c == 0)
        {
            if (locks[k].down == 0)
            {
                if (HAS_TAG(col_tags[k], C_TAG_LB_INF))
                {
                    iVec_append(cols_to_inf, -(int) k);
                }
                else
                {
                    assert(!IS_ABS_INF(bounds[k].lb));
                    fix_col(constraints, (int) k, bounds[k].lb, c[k], prob->obj);
                }
                continue;
            }

            if (locks[k].up == 0)
            {
                if (HAS_TAG(col_tags[k], C_TAG_UB_INF))
                {
                    iVec_append(cols_to_inf, (int) k);
                }
                else
                {
                    assert(!IS_ABS_INF(bounds[k].ub));
                    fix_col(constraints, (int) k, bounds[k].ub, c[k], prob->obj);
                }
                continue;
            }
        }
    }

    // now fix the columns that can be fixed to inf
    // (processing the inf columns after processing the other columns
    // simplifies the dual postsolve a lot)
    for (size_t i = 0; i < cols_to_inf->len; ++i)
    {
        int col = cols_to_inf->data[i];

        if (col < 0)
        {
            fix_col_to_negative_inf(constraints, -col);
        }
        else if (col > 0)
        {
            fix_col_to_positive_inf(constraints, col);
        }
        else
        {
            assert(col == 0);
            assert(locks[col].up == 0 || locks[col].down == 0);
            if (locks[col].down == 0)
            {
                fix_col_to_negative_inf(constraints, col);
            }
            else if (locks[col].up == 0)
            {
                fix_col_to_positive_inf(constraints, col);
            }
        }
    }

    // after the simple dual fix we might have inactive rows, and fixed
    // columns
    delete_inactive_rows(prob->constraints);
    delete_fixed_cols_from_problem(prob);
    delete_inactive_cols_from_A_and_AT(prob->constraints);

    return status;
}
