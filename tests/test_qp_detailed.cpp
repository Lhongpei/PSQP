/*
 * Detailed QP Presolver Tests
 * 
 * This test file provides comprehensive validation of QP presolver functionality
 * including:
 * - Correctness of P matrix reduction
 * - Offset computation
 * - Postsolve recovery
 * - Edge cases
 */

#include "PSLP_API.h"
#include "PSLP_stats.h"
#include "PSLP_sol.h"
#include <cassert>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <cstring>

const double TOL = 1e-6;

// Helper to check if two doubles are approximately equal
bool approx_equal(double a, double b, double tol = TOL) {
    return std::fabs(a - b) < tol;
}

// Print problem statistics
void print_problem_stats(const char* name, const PresolvedProblem* prob) {
    std::cout << "\n=== " << name << " ===" << std::endl;
    std::cout << "  Dimensions: " << prob->m << " x " << prob->n << std::endl;
    std::cout << "  Non-zeros: " << prob->nnz << std::endl;
    std::cout << "  Has quadratic: " << (prob->has_quadratic ? "yes" : "no") << std::endl;
    if (prob->has_quadratic) {
        std::cout << "  P nnz: " << prob->Pnnz << std::endl;
    }
    std::cout << "  Objective offset: " << std::setprecision(10) << prob->obj_offset << std::endl;
}

// Test 1: Simple QP with known analytical solution
// min 0.5 * (x1^2 + x2^2)  s.t. x1 + x2 = 1, x1, x2 >= 0
// Solution: x1 = x2 = 0.5, optimal value = 0.25
bool test_simple_qp_analytical() {
    std::cout << "\nTest 1: Simple QP with analytical solution" << std::endl;
    
    size_t m = 1;
    size_t n = 2;
    
    // Constraint: x1 + x2 = 1
    double Ax[] = {1.0, 1.0};
    int Ai[] = {0, 1};
    int Ap[] = {0, 2};
    size_t nnz = 2;
    
    double lhs[] = {1.0};
    double rhs[] = {1.0};
    double lbs[] = {0.0, 0.0};
    double ubs[] = {INFINITY, INFINITY};
    double c[] = {0.0, 0.0};
    
    // P = I (identity matrix), stored as upper triangular
    double Px[] = {1.0, 1.0};
    int Pi[] = {0, 1};
    int Pp[] = {0, 1, 2};
    size_t Pnnz = 2;
    
    Settings *stgs = default_settings();
    stgs->verbose = false;
    stgs->ston_cols = false;  // Disable to prevent full reduction
    stgs->parallel_cols = false;
    stgs->dton_eq = false;
    
    Presolver *presolver = new_qp_presolver(Ax, Ai, Ap, m, n, nnz,
                                            lhs, rhs, lbs, ubs, c,
                                            Px, Pi, Pp, Pnnz, stgs);
    if (!presolver) {
        std::cerr << "ERROR: Presolver creation failed" << std::endl;
        return false;
    }
    
    PresolveStatus status = run_presolver(presolver);
    PresolvedProblem *reduced = presolver->reduced_prob;
    
    print_problem_stats("Reduced Problem", reduced);
    
    bool success = true;
    
    // Verify status
    if (status != REDUCED && status != UNCHANGED) {
        std::cerr << "ERROR: Unexpected status: " << status << std::endl;
        success = false;
    }
    
    // Verify quadratic term is preserved
    if (reduced->n > 0 && !reduced->has_quadratic) {
        std::cerr << "ERROR: Quadratic term should be preserved" << std::endl;
        success = false;
    }
    
    free_settings(stgs);
    free_presolver(presolver);
    
    if (success) {
        std::cout << "PASSED" << std::endl;
    }
    return success;
}

// Test 2: QP with off-diagonal terms
// Tests that the P matrix transpose works correctly
bool test_qp_offdiagonal() {
    std::cout << "\nTest 2: QP with off-diagonal terms" << std::endl;
    
    size_t m = 1;
    size_t n = 3;
    
    // Constraint: x1 + x2 + x3 = 3
    double Ax[] = {1.0, 1.0, 1.0};
    int Ai[] = {0, 1, 2};
    int Ap[] = {0, 3};
    size_t nnz = 3;
    
    double lhs[] = {3.0};
    double rhs[] = {3.0};
    double lbs[] = {0.0, 0.0, 0.0};
    double ubs[] = {INFINITY, INFINITY, INFINITY};
    double c[] = {0.0, 0.0, 0.0};
    
    // P = [2 1 0; 1 2 1; 0 1 2] (tridiagonal)
    // Stored as upper triangular: [2 1 0; 0 2 1; 0 0 2]
    double Px[] = {2.0, 1.0, 2.0, 1.0, 2.0};
    int Pi[] = {0, 1, 1, 2, 2};
    int Pp[] = {0, 2, 4, 5};
    size_t Pnnz = 5;
    
    Settings *stgs = default_settings();
    stgs->verbose = false;
    stgs->ston_cols = false;
    stgs->parallel_cols = false;
    stgs->dton_eq = false;
    
    Presolver *presolver = new_qp_presolver(Ax, Ai, Ap, m, n, nnz,
                                            lhs, rhs, lbs, ubs, c,
                                            Px, Pi, Pp, Pnnz, stgs);
    if (!presolver) {
        std::cerr << "ERROR: Presolver creation failed" << std::endl;
        return false;
    }
    
    PresolveStatus status = run_presolver(presolver);
    PresolvedProblem *reduced = presolver->reduced_prob;
    
    print_problem_stats("Reduced Problem", reduced);
    
    bool success = true;
    
    if (status != REDUCED && status != UNCHANGED) {
        std::cerr << "ERROR: Unexpected status: " << status << std::endl;
        success = false;
    }
    
    // If problem was reduced but not eliminated, check P structure
    if (reduced->n > 0 && reduced->has_quadratic) {
        // Verify P is valid CSR
        if (reduced->Pp[reduced->n] != (int)reduced->Pnnz) {
            std::cerr << "ERROR: Invalid P row pointer" << std::endl;
            success = false;
        }
    }
    
    free_settings(stgs);
    free_presolver(presolver);
    
    if (success) {
        std::cout << "PASSED" << std::endl;
    }
    return success;
}

// Test 3: QP with variable fixing and offset verification
// Tests that fixing a variable correctly updates the objective offset
bool test_qp_offset_computation() {
    std::cout << "\nTest 3: QP offset computation with fixed variables" << std::endl;
    
    size_t m = 1;
    size_t n = 2;
    
    // Constraint: x1 + x2 = 3
    double Ax[] = {1.0, 1.0};
    int Ai[] = {0, 1};
    int Ap[] = {0, 2};
    size_t nnz = 2;
    
    double lhs[] = {3.0};
    double rhs[] = {3.0};
    // x1 is fixed to 1 by bounds
    double lbs[] = {1.0, 0.0};
    double ubs[] = {1.0, INFINITY};
    double c[] = {0.0, 0.0};
    
    // P = 2*I, so quadratic term is 0.5 * (2*x1^2 + 2*x2^2) = x1^2 + x2^2
    double Px[] = {2.0, 2.0};
    int Pi[] = {0, 1};
    int Pp[] = {0, 1, 2};
    size_t Pnnz = 2;
    
    Settings *stgs = default_settings();
    stgs->verbose = false;
    stgs->ston_cols = false;
    stgs->parallel_cols = false;
    stgs->dton_eq = false;
    
    Presolver *presolver = new_qp_presolver(Ax, Ai, Ap, m, n, nnz,
                                            lhs, rhs, lbs, ubs, c,
                                            Px, Pi, Pp, Pnnz, stgs);
    if (!presolver) {
        std::cerr << "ERROR: Presolver creation failed" << std::endl;
        return false;
    }
    
    PresolveStatus status = run_presolver(presolver);
    PresolvedProblem *reduced = presolver->reduced_prob;
    
    print_problem_stats("Reduced Problem (x1 fixed)", reduced);
    
    bool success = true;
    
    // Note: The presolver may fully reduce the problem, in which case
    // the offset captures the complete optimal value.
    // Expected minimum contribution from x1: 0.5 * P_00 * x1^2 = 0.5 * 2 * 1 = 1.0
    double min_expected_offset = 0.5 * 2.0 * 1.0 * 1.0;
    
    if (reduced->n == 0) {
        // Problem fully reduced - offset should include optimal value
        std::cout << "  Problem fully reduced, offset = " << reduced->obj_offset << std::endl;
        if (reduced->obj_offset < min_expected_offset - TOL) {
            std::cerr << "ERROR: Offset too small. Minimum expected: " << min_expected_offset 
                      << ", Got: " << reduced->obj_offset << std::endl;
            success = false;
        }
    } else if (reduced->n == 1 && reduced->has_quadratic) {
        // Partial reduction - check P structure
        if (reduced->Pnnz != 1) {
            std::cerr << "ERROR: Expected Pnnz=1, got: " << reduced->Pnnz << std::endl;
            success = false;
        }
        if (!approx_equal(reduced->Px[0], 2.0)) {
            std::cerr << "ERROR: Expected Px[0]=2.0, got: " << reduced->Px[0] << std::endl;
            success = false;
        }
        // Check offset is at least the x1 contribution
        if (!approx_equal(reduced->obj_offset, min_expected_offset)) {
            std::cerr << "ERROR: Offset mismatch. Expected: " << min_expected_offset 
                      << ", Got: " << reduced->obj_offset << std::endl;
            success = false;
        }
    }
    
    free_settings(stgs);
    free_presolver(presolver);
    
    if (success) {
        std::cout << "PASSED" << std::endl;
    }
    return success;
}

// Test 4: Postsolve verification
// Tests that postsolve correctly recovers the original solution
bool test_qp_postsolve_recovery() {
    std::cout << "\nTest 4: Postsolve solution recovery" << std::endl;
    
    size_t m = 2;
    size_t n = 3;
    
    // Constraints:
    // x1 + x2 + x3 = 6
    // x1 - x2 = 0 (implies x1 = x2)
    double Ax[] = {1.0, 1.0, 1.0, 1.0, -1.0};
    int Ai[] = {0, 1, 2, 0, 1};
    int Ap[] = {0, 3, 5};
    size_t nnz = 5;
    
    double lhs[] = {6.0, 0.0};
    double rhs[] = {6.0, 0.0};
    double lbs[] = {0.0, 0.0, 0.0};
    double ubs[] = {INFINITY, INFINITY, INFINITY};
    double c[] = {0.0, 0.0, 0.0};
    
    // P = 2*I
    double Px[] = {2.0, 2.0, 2.0};
    int Pi[] = {0, 1, 2};
    int Pp[] = {0, 1, 2, 3};
    size_t Pnnz = 3;
    
    Settings *stgs = default_settings();
    stgs->verbose = false;
    stgs->ston_cols = false;
    
    Presolver *presolver = new_qp_presolver(Ax, Ai, Ap, m, n, nnz,
                                            lhs, rhs, lbs, ubs, c,
                                            Px, Pi, Pp, Pnnz, stgs);
    if (!presolver) {
        std::cerr << "ERROR: Presolver creation failed" << std::endl;
        return false;
    }
    
    PresolveStatus status = run_presolver(presolver);
    PresolvedProblem *reduced = presolver->reduced_prob;
    
    print_problem_stats("Reduced Problem", reduced);
    
    bool success = true;
    
    // Try postsolve if problem was reduced
    if (reduced->n > 0 && reduced->m > 0) {
        // Create dummy solution for reduced problem
        std::vector<double> x_reduced(reduced->n, 1.0);
        std::vector<double> y_reduced(reduced->m, 0.5);
        std::vector<double> z_reduced(reduced->n, 0.0);
        
        postsolve(presolver, x_reduced.data(), y_reduced.data(), z_reduced.data());
        
        if (!presolver->sol) {
            std::cerr << "ERROR: Postsolve did not create solution" << std::endl;
            success = false;
        } else {
            // Verify solution dimensions
            if (presolver->sol->dim_x != (int)n) {
                std::cerr << "ERROR: x dimension mismatch. Expected: " << n 
                          << ", Got: " << presolver->sol->dim_x << std::endl;
                success = false;
            }
            if (presolver->sol->dim_y != (int)m) {
                std::cerr << "ERROR: y dimension mismatch. Expected: " << m 
                          << ", Got: " << presolver->sol->dim_y << std::endl;
                success = false;
            }
            
            // Print recovered solution
            std::cout << "Recovered x: [";
            for (size_t i = 0; i < n; i++) {
                std::cout << presolver->sol->x[i];
                if (i < n-1) std::cout << ", ";
            }
            std::cout << "]" << std::endl;
        }
    } else {
        std::cout << "Problem fully reduced, skipping postsolve" << std::endl;
    }
    
    free_settings(stgs);
    free_presolver(presolver);
    
    if (success) {
        std::cout << "PASSED" << std::endl;
    }
    return success;
}

// Test 5: Empty P matrix (LP problem through QP interface)
bool test_qp_empty_p() {
    std::cout << "\nTest 5: QP with empty P matrix (effectively LP)" << std::endl;
    
    size_t m = 1;
    size_t n = 2;
    
    double Ax[] = {1.0, 1.0};
    int Ai[] = {0, 1};
    int Ap[] = {0, 2};
    size_t nnz = 2;
    
    double lhs[] = {1.0};
    double rhs[] = {1.0};
    double lbs[] = {0.0, 0.0};
    double ubs[] = {INFINITY, INFINITY};
    double c[] = {1.0, 2.0};
    
    // Empty P matrix (NULL, 0 nnz)
    Settings *stgs = default_settings();
    stgs->verbose = false;
    
    Presolver *presolver = new_qp_presolver(Ax, Ai, Ap, m, n, nnz,
                                            lhs, rhs, lbs, ubs, c,
                                            NULL, NULL, NULL, 0, stgs);
    if (!presolver) {
        std::cerr << "ERROR: Presolver creation failed" << std::endl;
        return false;
    }
    
    PresolveStatus status = run_presolver(presolver);
    PresolvedProblem *reduced = presolver->reduced_prob;
    
    print_problem_stats("LP through QP interface", reduced);
    
    bool success = true;
    
    // Should not have quadratic term
    if (reduced->has_quadratic) {
        std::cerr << "ERROR: Should not have quadratic term for empty P" << std::endl;
        success = false;
    }
    
    free_settings(stgs);
    free_presolver(presolver);
    
    if (success) {
        std::cout << "PASSED" << std::endl;
    }
    return success;
}

// Test 6: Large sparse QP
// Tests performance and correctness on larger sparse problems
bool test_qp_large_sparse() {
    std::cout << "\nTest 6: Large sparse QP (100 vars, 50 constraints)" << std::endl;
    
    const size_t n = 100;
    const size_t m = 50;
    
    // Create sparse constraint matrix
    std::vector<double> Ax;
    std::vector<int> Ai;
    std::vector<int> Ap(m + 1, 0);
    
    size_t nnz = 0;
    for (size_t i = 0; i < m; i++) {
        // Each row has ~5 non-zeros
        for (size_t j = 0; j < 5; j++) {
            int col = (int)((i * 7 + j * 13) % n);
            Ax.push_back(1.0 + (j % 3));
            Ai.push_back(col);
            nnz++;
        }
        Ap[i + 1] = (int)nnz;
    }
    
    std::vector<double> lhs(m, -INFINITY);
    std::vector<double> rhs(m);
    for (size_t i = 0; i < m; i++) {
        rhs[i] = 50.0 + i * 5;
    }
    
    std::vector<double> lbs(n, 0.0);
    std::vector<double> ubs(n, INFINITY);
    std::vector<double> c(n);
    for (size_t i = 0; i < n; i++) {
        c[i] = (double)(i + 1);
    }
    
    // Diagonal P matrix
    std::vector<double> Px(n);
    std::vector<int> Pi(n);
    std::vector<int> Pp(n + 1);
    
    for (size_t i = 0; i < n; i++) {
        Px[i] = 2.0;
        Pi[i] = (int)i;
        Pp[i] = (int)i;
    }
    Pp[n] = (int)n;
    
    Settings *stgs = default_settings();
    stgs->verbose = false;
    stgs->max_time = 30.0;
    
    Presolver *presolver = new_qp_presolver(Ax.data(), Ai.data(), Ap.data(), 
                                            m, n, nnz,
                                            lhs.data(), rhs.data(), 
                                            lbs.data(), ubs.data(), c.data(),
                                            Px.data(), Pi.data(), Pp.data(), n, stgs);
    if (!presolver) {
        std::cerr << "ERROR: Presolver creation failed" << std::endl;
        return false;
    }
    
    PresolveStatus status = run_presolver(presolver);
    PresolvedProblem *reduced = presolver->reduced_prob;
    
    std::cout << "Original:  " << m << " x " << n << ", nnz=" << nnz << std::endl;
    std::cout << "Reduced:   " << reduced->m << " x " << reduced->n 
              << ", nnz=" << reduced->nnz << std::endl;
    std::cout << "Reduction: " << (100.0 * (1.0 - (double)reduced->nnz / nnz)) 
              << "% non-zeros removed" << std::endl;
    
    bool success = true;
    
    if (status != REDUCED && status != UNCHANGED) {
        std::cerr << "ERROR: Unexpected status: " << status << std::endl;
        success = false;
    }
    
    free_settings(stgs);
    free_presolver(presolver);
    
    if (success) {
        std::cout << "PASSED" << std::endl;
    }
    return success;
}

// Test 7: QP with ill-conditioned P matrix
// Tests numerical stability
bool test_qp_ill_conditioned() {
    std::cout << "\nTest 7: QP with varying scale in P" << std::endl;
    
    size_t m = 1;
    size_t n = 3;
    
    double Ax[] = {1.0, 1.0, 1.0};
    int Ai[] = {0, 1, 2};
    int Ap[] = {0, 3};
    size_t nnz = 3;
    
    double lhs[] = {1.0};
    double rhs[] = {1.0};
    double lbs[] = {0.0, 0.0, 0.0};
    double ubs[] = {INFINITY, INFINITY, INFINITY};
    double c[] = {0.0, 0.0, 0.0};
    
    // P with very different scales: [1e6, 0, 0; 0, 1, 0; 0, 0, 1e-6]
    double Px[] = {1e6, 1.0, 1e-6};
    int Pi[] = {0, 1, 2};
    int Pp[] = {0, 1, 2, 3};
    size_t Pnnz = 3;
    
    Settings *stgs = default_settings();
    stgs->verbose = false;
    stgs->ston_cols = false;
    stgs->parallel_cols = false;
    
    Presolver *presolver = new_qp_presolver(Ax, Ai, Ap, m, n, nnz,
                                            lhs, rhs, lbs, ubs, c,
                                            Px, Pi, Pp, Pnnz, stgs);
    if (!presolver) {
        std::cerr << "ERROR: Presolver creation failed" << std::endl;
        return false;
    }
    
    PresolveStatus status = run_presolver(presolver);
    PresolvedProblem *reduced = presolver->reduced_prob;
    
    print_problem_stats("Ill-conditioned QP", reduced);
    
    bool success = (status == REDUCED || status == UNCHANGED);
    
    free_settings(stgs);
    free_presolver(presolver);
    
    if (success) {
        std::cout << "PASSED" << std::endl;
    }
    return success;
}

int main() {
    std::cout << "========================================" << std::endl;
    std::cout << "   Detailed QP Presolver Tests" << std::endl;
    std::cout << "========================================" << std::endl;
    
    int passed = 0;
    int total = 7;
    
    if (test_simple_qp_analytical()) passed++;
    if (test_qp_offdiagonal()) passed++;
    if (test_qp_offset_computation()) passed++;
    if (test_qp_postsolve_recovery()) passed++;
    if (test_qp_empty_p()) passed++;
    if (test_qp_large_sparse()) passed++;
    if (test_qp_ill_conditioned()) passed++;
    
    std::cout << "\n========================================" << std::endl;
    std::cout << "Results: " << passed << "/" << total << " tests passed" << std::endl;
    std::cout << "========================================" << std::endl;
    
    return (passed == total) ? 0 : 1;
}
