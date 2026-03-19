# PSQP — A Lightweight C Presolver for Quadratic Programs

[![License: Apache 2.0](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![Language: C](https://img.shields.io/badge/Language-C-blue.svg)](https://en.wikipedia.org/wiki/C_(programming_language))

**PSQP** is a fast, dependency-free presolver for quadratic programs (QPs) of the form

$$
\begin{array}{ll}
\text{minimize} & \frac{1}{2} x^T (Q + RR^T) x + c^T x \\
\text{subject to} & \underline{b} \leq Ax \leq \overline{b} \\
& \underline{x} \leq x \leq \overline{x}.
\end{array}
$$

where $Q$ is a sparse symmetric matrix and $R$ is a sparse $n \times k$ matrix with $k \ll n$. This $P = Q + RR^T$ structure enables memory-efficient handling of large-scale QPs.

**PSQP** extends [PSLP](https://github.com/dance858/PSLP) with full QP support. It is written in C99 with no external dependencies.

---

## Installation

**PSQP** is built using CMake:

```bash
mkdir build && cd build
cmake ..
make
sudo make install
```

For testing:
```bash
cmake .. -DPSQP_BUILD_TESTING=ON
make all_tests
./all_tests
```

---

## API

The public C API is defined in `PSLP_API.h`. The API consists of three main operations:

1. **Initialization** — performed using `new_qp_presolver_qr()`, which creates and initializes internal data structures.

2. **Presolve** — performed using `run_presolver()`, which executes presolve routines to reduce the QP. The reduced problem is available in `presolver->reduced_prob`.

3. **Postsolve** — performed using `postsolve()`, which recovers a primal-dual solution to the original problem from a solution to the reduced problem.


---

## Citation

```bibtex
@software{Cederberg2025,
  author = {Cederberg, Daniel and Boyd, Stephen},
  title = {PSLP — A Lightweight C Presolver for Linear Programs},
  year = {2025},
  url = {https://github.com/dance858/PSLP},
}


@software{PSQP2026,
  author = {Hongpei Li},
  title = {PSQP — A Lightweight C Presolver for Quadratic Programs},
  year = {2026},
  url = {https://github.com/Lhongpei/PSQP}
}
```

---

## License

Licensed under the Apache License, Version 2.0.  
Original PSLP code Copyright 2025 Daniel Cederberg.  
QP extensions Copyright 2026 Hongpei Li.
