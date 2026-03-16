# Quadratic Programming Solver

This repository implements a **Quadratic Programming (QP) solver** using a transformation to a **Linear Complementarity Problem (LCP)** and solving it with a **Revised Simplex method**.

The solver handles quadratic optimization problems of the form:

\[
\max \; \frac{1}{2} x^T Q x + c^T x
\]

subject to

\[
Ax \le b
\]

\[
x \ge 0
\]

where:

- \(Q \in \mathbb{R}^{n \times n}\) is a symmetric matrix  
- \(c \in \mathbb{R}^{n}\) is the linear coefficient vector  
- \(A \in \mathbb{R}^{m \times n}\) is the constraint matrix  
- \(b \in \mathbb{R}^{m}\) is the constraint vector  

The solver checks whether **Q is negative semidefinite** to ensure the problem is **convex (CMP)**.

---

# Files in the Repository

- `quadratic_programming.py`  
  Main program that converts the quadratic program into a form solvable by simplex.

- `Simplex_revised.py`  
  Implementation of the **Revised Simplex algorithm** with complementary pivot rules.

- `input.txt`  
  Input file containing the quadratic program.

---

# How to Run

Make sure the following files are in the **same directory**:

- `quadratic_programming.py`
- `Simplex_revised.py`
- `input.txt`

Then run:

```bash
python quadratic_programming.py
```

---

# Input File Format

`input.txt`

```
# n
# Q[0][0] Q[0][1] ... Q[0][n-1]
# ...
# Q[n-1][0] ... Q[n-1][n-1]

# c[0] c[1] ... c[n-1]

# m

# A[0][0] ... A[0][n-1] b[0]
# ...
# A[m-1][0] ... A[m-1][n-1] b[m-1]
```

### Description

- `n` → number of decision variables  
- `Q` → quadratic coefficient matrix  
- `c` → linear coefficient vector  
- `m` → number of constraints  
- `A` → constraint matrix  
- `b` → right-hand side vector  

---

# Output

The solver prints:

- the **optimal value of the objective function**
- the **solution vector**

The **first `n` elements of the output correspond to the decision variables**:

\[
x = (x_1, x_2, ..., x_n)
\]

---

# Example

Example quadratic program:

\[
\max \frac{1}{2} x^T Q x + c^T x
\]

subject to

\[
Ax \le b
\]

The solver reads the problem from `input.txt`, runs the simplex-based solver, and outputs the optimal solution.

---

# Requirements

- Python 3
- NumPy

Install NumPy if needed:

```bash
pip install numpy
```

