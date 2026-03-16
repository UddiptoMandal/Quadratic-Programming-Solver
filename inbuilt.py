# save as solve_qp_scipy.py
import numpy as np
from scipy.optimize import minimize

Q = np.array([[-2.0, 0.0],[0.0, -2.0]])
c = np.array([4.0, 6.0])

A = np.array([[1.0,1.0],
              [1.0,0.0],
              [0.0,1.0]])
b = np.array([4.0, 3.0, 3.0])

def obj_to_min(x):
    # minimize negative of original objective
    return - (0.5 * x.T @ Q @ x + c @ x)

cons = [{'type': 'ineq', 'fun': lambda x, Ai=A[i], bi=b[i]: bi - Ai.dot(x)} for i in range(A.shape[0])]
bounds = [(0, None) for _ in range(Q.shape[0])]

x0 = np.clip(np.linalg.lstsq(A, b, rcond=None)[0], 0, None)  # a reasonable starting point
res = minimize(obj_to_min, x0, method='SLSQP', bounds=bounds, constraints=cons,
               options={'ftol':1e-9, 'maxiter':1000})

print("success:", res.success)
print("x:", res.x)
print("objective value:", 0.5 * res.x.T @ Q @ res.x + c @ res.x)