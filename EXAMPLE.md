# PSQP Quick Start Example

This example shows how to use PSQP to presolve a simple quadratic program (QP).

## Problem Formulation

Consider the following QP:

```
minimize    0.5 * x^T P x + c^T x
subject to  lhs <= A x <= rhs
            lbs <= x <= ubs
```

where:
- `P` is a sparse symmetric positive semi-definite matrix (only upper triangular part stored)
- `A` is the constraint matrix in CSR format
- `c`, `lhs`, `rhs`, `lbs`, `ubs` are vectors

## Example Code

```c
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "PSQP_API.h"

int main() {
    // Problem dimensions
    size_t n = 3;  // variables
    size_t m = 2;  // constraints
    
    // Constraint matrix A (2x3, CSR format)
    // A = [1.0  1.0  0.0]
    //     [0.0  1.0  1.0]
    double Ax[] = {1.0, 1.0, 1.0, 1.0};
    int Ai[] = {0, 1, 1, 2};
    int Ap[] = {0, 2, 4};
    size_t nnz = 4;
    
    // Constraint bounds: 1 <= A*x <= inf
    double lhs[] = {1.0, 1.0};
    double rhs[] = {INFINITY, INFINITY};
    
    // Variable bounds: 0 <= x <= inf
    double lbs[] = {0.0, 0.0, 0.0};
    double ubs[] = {INFINITY, INFINITY, INFINITY};
    
    // Linear objective: c^T x
    double c[] = {1.0, 1.0, 1.0};
    
    // Quadratic term P (diagonal matrix, upper triangular CSR)
    // P = diag(2.0, 2.0, 2.0)
    double Px[] = {2.0, 2.0, 2.0};
    int Pi[] = {0, 1, 2};
    int Pp[] = {0, 1, 2, 3};
    size_t Pnnz = 3;
    
    // Create default settings
    Settings *settings = default_settings();
    settings->verbose = true;
    
    // Initialize presolver for QP
    Presolver *presolver = new_qp_presolver(
        Ax, Ai, Ap, m, n, nnz,
        lhs, rhs, lbs, ubs, c,
        Px, Pi, Pp, Pnnz,
        settings
    );
    
    if (!presolver) {
        fprintf(stderr, "Failed to create presolver\n");
        return 1;
    }
    
    // Run presolve
    PresolveStatus status = run_presolver(presolver);
    
    printf("Presolve status: %d\n", status);
    printf("Original:  %zu rows, %zu columns, %zu nonzeros\n", m, n, nnz);
    printf("Reduced:   %zu rows, %zu columns, %zu nonzeros\n",
           presolver->reduced_prob->m,
           presolver->reduced_prob->n,
           presolver->reduced_prob->nnz);
    
    // Allocate solution arrays (size of original problem)
    double *x = calloc(n, sizeof(double));
    double *y = calloc(m, sizeof(double));
    double *z = calloc(n, sizeof(double));
    
    // Set reduced solution (would come from solver)
    // For this example, we just set zeros
    
    // Postsolve to recover original solution
    postsolve(presolver, x, y, z);
    
    // Print recovered solution
    printf("\nRecovered primal solution:\n");
    for (size_t i = 0; i < n; i++) {
        printf("  x[%zu] = %g\n", i, x[i]);
    }
    
    // Cleanup
    free(x); free(y); free(z);
    free_presolver(presolver);
    free_settings(settings);
    
    return 0;
}
```

## Compilation

```bash
# Assuming PSQP is installed
cc -o example example.c -lPSQP -lm

# Or with local build
cc -o example example.c -I/path/to/PSQP/include/PSQP -L/path/to/PSQP/build -lPSQP -lm -Wl,-rpath,/path/to/PSQP/build
```

## Expected Output

```
Running presolver (PSQP x.x.x)...
  status          : REDUCED (or UNCHANGED)
  presolve time   : 0.000xxx sec
  reduced problem : X rows, Y columns, Z nonzeros
Original:  2 rows, 3 columns, 4 nonzeros
Reduced:   X rows, Y columns, Z nonzeros

Recovered primal solution:
  x[0] = ...
  x[1] = ...
  x[2] = ...
```

## API Summary

### Key Functions

| Function | Description |
|----------|-------------|
| `default_settings()` | Create default presolve settings |
| `new_qp_presolver()` | Create QP presolver with P matrix |
| `new_qp_presolver_qr()` | Create QP presolver with Q+R·Rᵀ decomposition |
| `run_presolver()` | Execute presolve routines |
| `postsolve()` | Recover original solution from reduced solution |
| `free_presolver()` | Free presolver memory |
| `free_settings()` | Free settings memory |

### Data Structures

- **`PresolvedProblem`**: Contains reduced problem data (A, P, bounds, etc.)
- **`Presolver`**: Main presolver object with internal state
- **`Settings`**: Configuration options (toggles for different reductions)

## More Examples

For more comprehensive examples, see:

- [`tests/qp_example.c`](tests/qp_example.c) - QP examples with standard P matrix
  - Portfolio optimization
  - Fixed variables handling
  - Postsolve demonstration
  - LP through QP interface

- [`tests/qr_example.c`](tests/qr_example.c) - QR decomposition examples
  - Factor model portfolio optimization (P = D + F·Fᵀ)
  - Low-rank plus diagonal structure
  - QR matrix shrinking during presolve

## See Also

- [README.md](README.md) - Full documentation
