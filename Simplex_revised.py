import sys
import numpy  as np

# input: 
# n
# m
# c1 c2 ... cn
# a11 a12 ... a1n b1
# ...
# am1 am2 ... amn bm
# k
# ind1 ind2
# ...

class revised_simplex:
    def __init__(self, n, m, c, A, b, k, ind1, ind2, rh):
        self.n = n
        self.m = m
        self.c_new = -1*c
        self.A_new = A
        self.b = b
        self.rh = rh
        self.k = k
        self.ind1 = ind1
        self.ind2 = ind2
        for _ in range(self.k):
            assert self.ind1[_] != self.ind2[_], "Indices must be different"

        self.basis = [self.n + i for i in range(self.m)]

        self.ind1_basic = [0 for _ in range(self.k)]
        self.ind2_basic = [0 for _ in range(self.k)]


        self.c_new = np.array(self.c_new, dtype=float)
        self.A_new = np.array(self.A_new, dtype=float)
        self.b = np.array(self.b, dtype=float)  

    def is_in_basis(self, var_index):
        return var_index in self.basis


    def enterring_basic(self, index):
        for counter in range(self.k):
            if self.ind1[counter] == index and self.is_in_basis(self.ind2[counter]):      
                return False

            if self.ind2[counter] == index and self.is_in_basis(self.ind1[counter]):
                return False
        return True
    def optimality(self):
        for i in range(len(self.c_new)):
            if self.c_new[i] < 0 and self.enterring_basic(i):
                return False
        return True
    

    def solve(self):

        iteration = 0
        while(True):
            if self.optimality():
                print("Optimal solution found")
                break
            min_val = 0
            index = -1
            for i in range(0, len(self.c_new)):
                if self.c_new[i] < min_val and self.enterring_basic(i):
                    min_val = self.c_new[i]
                    index = i
            if index == -1:
                print("Problem is infeasible or stuck in a loop")
                break
            
            minratio = sys.maxsize
            index2 = -1

            for _ in range(self.m):
                if self.A_new[_][index] > 0:

                    ratio = self.b[_] / self.A_new[_][index]
                    if ratio < minratio:
                        minratio = ratio
                        index2 = _
            if index2 == -1:
                print("Unbounded solution")
                break
            leaving_var = self.basis[index2]
            print(f"Iteration: {iteration}, Entering variable: {index}, Leaving variable: {leaving_var}")
            self.basis[index2] = index

            pivot = self.A_new[index2][index]

            self.A_new[index2] = self.A_new[index2] / pivot
            self.b[index2] = self.b[index2] / pivot
            factor = self.c_new[index]

            self.c_new = self.c_new - factor * self.A_new[index2]
            self.rh = self.rh - factor * self.b[index2]
            for i in range(self.m):
                if self.A_new[i][index] != 0 and i!= index2:
                    c1 = self.A_new[i][index]
                    self.A_new[i] -= c1 * self.A_new[index2]
                    self.b[i] -= c1 * self.b[index2]
            iteration += 1
            print(f"After iteration {iteration}, basis: {self.basis}, rh: {self.rh}, c_new: {self.c_new}")
            print(f"A_new: {self.A_new}, b: {self.b}")

        x = np.zeros(self.n + self.m)
        for _ in range(self.m):
            x[self.basis[_]] = self.b[_]
        print(f"Final solution: {x}, objective value: {self.rh}")
        print(f"Final basis: {self.basis}")
        print(f"Final c_new: {self.c_new}")
        print(f"Final A_new: {self.A_new}, b: {self.b}")


        return x, self.rh
    
    def print_solution(self, x):
        print("Variable values:")
        for i in range(self.n):
            print(f"x[{i}] = {x[i]:.6f}")
        print(f"Objective value : {self.rh:.6f}")
        print(f"Reduced costs: {self.c_new}")

