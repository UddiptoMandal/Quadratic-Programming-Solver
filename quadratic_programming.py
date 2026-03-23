
import numpy as np
from Simplex_revised import revised_simplex

# input.txt format:
# n
# Q[0][0] Q[0][1] ... Q[0][n-1]
# ...
# Q[n-1][0] ... Q[n-1][n-1]
# c[0] c[1] ... c[n-1]
# m
# A[0][0] ... A[0][n-1] b[0]
# ...
# A[m-1][0] ... A[m-1][n-1] b[m-1]

with open("input.txt", "r") as f:
    data = f.read().split()
ptr = 0
def get_int():
    global ptr
    val = int(data[ptr]); ptr += 1
    return val

def get_float():
    global ptr
    val = float(data[ptr]); ptr += 1
    return val


n = get_int()

Q = [[get_float() for _ in range(n)] for _ in range(n)]
Q = np.array(Q, dtype=float)
c = [get_float() for _ in range(n)]
c = np.array(c, dtype=float)

# max 1/2 xTQx + cTx
# st Ax <= b
m = get_int()
A = [[0.0 for _ in range(n)] for _ in range(m)]
b = [0.0 for _ in range(m)]
print("Input A and b: ")
for i in range(m):
    for j in range(n):
        A[i][j] = get_float()
    b[i] = get_float()
A = np.array(A, dtype=float)
b = np.array(b, dtype=float)
 
assert np.allclose(Q.T, Q), "Q must be symmetric"
eigenvalues = np.linalg.eigvals(Q)
if np.all(eigenvalues <= 0):
    print("Q is negative semidefinite, the problem is CMP.")
else:
    print("Q is not negative semidefinite, the problem is not CMP.")
    exit(1)
# equations: 
# Q = nxn
# c = nx1
# A = mxn
# b = mx1
# variables: x = nx1, s = mx1, u = mx1
# 1. Qx + c = AT u
# 2. Ax + s = b
# 3. s >=0, u >= 0, x >= 0
# 4. ut s = 0
# x, then s, then u, then a, then w
# max -M(a1+ a2..an + w1 + w2..wm) 
# Qx + c - AT u + [a1..an] <= 0
# Ax + s - b + [w1..wm] <= 0
# a1, a2..an >= 0, w1, w2..wm >= 0
# ut s = 0
num_var = n + m + m + n + m
num_eq = n + m
num_complementary = m
M = 1e6
c_lp = [0 for _ in range(n+m+m)] + [-M for _ in range(n+m)]
b_lp = [0 for _ in range(num_eq)]
A_lp = np.zeros((num_eq, num_var))
for _ in range(n):
    A_lp[_, 0:n] = -1*Q[_, :]
    A_lp[_, n: n+m] = 0
    A_lp[_, n+m: n+m+m] = A[:, _]
    A_lp[_, n+m+m+_] = 1
    A_lp[_, n+m+m+n: ] = 0
    b_lp[_] = c[_]



for _ in range(m):
    A_lp[n+_, 0:n] = A[_, :]
    A_lp[n+_, n+_] = 1
    A_lp[n+_, n+m: n+m+m] = 0
    A_lp[n+_, n+m+m: n+m+m+n] = 0
    A_lp[n+_, n+m+m+n+_] = 1
    b_lp[n+_] = b[_]
rh = 0

def Gauss_jordan_converter(A_lp, b_lp, rh, c_lp, n, m):
    # rows: n+m
    # cols: n+m+m+n+m
    # basic: n+m+m to n+m+m+n+m
    # A from n+m+m to n+m+m+n+m is identity matrix
    for i in range(n+m):
        c_lp += M * A_lp[i]
        rh += M*b_lp[i]
    return A_lp, b_lp, rh, c_lp



ind1_lp = [n+i for i in range(m)]
ind2_lp = [n+m+i for i in range(m)]
print("num_var:", num_var, "num_eq:", num_eq)
print("c_lp", c_lp)
print("A:", A_lp)
print("b:", b_lp)

# here c, A, rh , b makes a simplex table not in gauss jordan form
A_lp, b_lp, rh, c_lp = Gauss_jordan_converter(A_lp, b_lp, rh, c_lp, n, m)
print("complementary indices", ind1_lp, ind2_lp)
print("num_var:", num_var-num_eq, "num_eq:", num_eq)
print("c_lp", c_lp)
print("A:", A_lp)
print("b:", b_lp)
print(num_complementary)
print("rh:", rh)

A_lp = A_lp.tolist()
solver = revised_simplex(num_var-num_eq, num_eq, c_lp, A_lp, b_lp, num_complementary, ind1_lp, ind2_lp, rh)
x, rh = solver.solve()
output_objective = 0.5 * x[:n].T @ Q @ x[:n] + c @ x[:n]
print(f"Val of objective function: {output_objective}, x: {x}")






    