# CFD Homework 20171842 김희성
# Programming Language/ver : Python 3.9
# Due: 08. December, 2021
# Example 7.1 Tri-Diagonal Matrix Algorithm(TDMA) for example 4.3

# Imports
import numpy as np
import math
import pandas as pd
import matplotlib.pyplot as plt

# Introduction
print('CFD Homework 20171842 김희성')
print('Programming Language/ver : Python 3.9')
print('Due: 08. December, 2021')
print('Example 7.1 Tri-Diagonal Matrix Algorithm(TDMA) for example 4.3')
print('-----------------------------------------------------------------------')

# Constants
L = 1           # meter
q = 0           # Watt per meter cubic
T_B = 100       # deg Celsius
T_ambient = 20  # deg Celsius
n_square = 25   # Watt

# Input Variables
n = int(input('Type the # of CVs: '))    # Number of Control Volume
delta_x = L / n

# Initialize Variables
x = np.zeros(n)
sol = np.zeros(n)
T = np.zeros(n)
alpha = np.zeros(n)
beta = np.zeros(n)
S_P = np.zeros(n)
D = np.zeros(n)
A = np.zeros(n)
C = np.zeros(n)
Cprime = np.zeros(n)
error = np.zeros(n)

# Solve
for j in range(n):
    if j == 0:                # Exterior Node (Node 1)
        x[j] = delta_x / 2
        sol[j] = math.cosh(n_square ** 0.5 * (L - x[j])) * (T_B - T_ambient) / math.cosh(
            n_square ** 0.5 * L) + T_ambient
        beta[j] = 0
        alpha[j] = 1 / delta_x
        C[j] = n_square * delta_x * T_ambient + 2 * T_B / delta_x
        S_P[j] = -n_square * delta_x - 2 / delta_x
        D[j] = beta[j] + alpha[j] - S_P[j]
        A[j] = alpha[j] / D[j]
        Cprime[j] = C[j] / D[j]
    elif j < n-1:             # Interior Nodes
        x[j] = x[0] + j * delta_x
        sol[j] = math.cosh(n_square ** 0.5 * (L - x[j])) * (T_B - T_ambient) / math.cosh(
            n_square ** 0.5 * L) + T_ambient
        beta[j] = 1 / delta_x
        alpha[j] = 1 / delta_x
        C[j] = n_square * delta_x * T_ambient
        S_P[j] = -n_square * delta_x
        D[j] = beta[j] + alpha[j] - S_P[j]
        A[j] = alpha[j] / (D[j] - A[j-1] * beta[j])
        Cprime[j] = (beta[j] * Cprime[j-1] + C[j]) / (D[j] - A[j-1] * beta[j])
    elif j == n-1:             # Exterior Node (Last Node)
        x[j] = x[0] + j * delta_x
        sol[j] = math.cosh(n_square ** 0.5 * (L - x[j])) * (T_B - T_ambient) / math.cosh(
            n_square ** 0.5 * L) + T_ambient
        beta[j] = 1 / delta_x
        alpha[j] = 0
        C[j] = n_square * delta_x * T_ambient
        S_P[j] = -n_square * delta_x
        D[j] = beta[j] + alpha[j] - S_P[j]
        A[j] = alpha[j] / (D[j] - A[j-1] * beta[j])
        Cprime[j] = (beta[j] * Cprime[j-1] + C[j]) / (D[j] - A[j-1] * beta[j])
    else:                       # if something went wrong...
        break

# Solution
for k in range(n-1, -1, -1):
    if k == n-1:
        T[n - 1] = Cprime[n - 1]
    else:
        T[k] = A[k] * T[k + 1] + Cprime[k]
    error[k] = abs(T[k] - sol[k]) / sol[k] * 100

# Drawing Table
nodes = [_ + 1 for _ in range(n)]
data1 = {'Node': nodes, 'beta': beta, 'D': D, 'alpha': alpha, 'C': C,
         'A': A, 'Cprime': Cprime}
data2 = {'Node': nodes, 'x(m)': x, 'TDMA solution': T,
         'Exact solution': sol, 'Percentage error': error}
Table1 = pd.DataFrame(data1)
Table2 = pd.DataFrame(data2)
print('[Table 7.1]')
print(Table1)
print('-----------------------------------------------------------------------')
print('[Comparison with the TDMA solution]')
print(Table2)

# Drawing Graph
x_min = 0
x_max = L
y_min = 20
y_max = 100
x_sol = np.linspace(x_min, x_max, 1000)
y_sol = np.zeros(1000)
for space in range(len(y_sol)):
    y_sol[space] = math.cosh(n_square ** 0.5 * (L - x_sol[space])) \
                   * (T_B - T_ambient) / math.cosh(n_square ** 0.5 * L) + T_ambient

for m in range(n):
    plt.plot(x[m], T[m], "b*")
    # style : :(dot),-(line),--(dash line),^(tri),*(bold)
    # color : R, G, B, C, M, Y, K
plt.plot(0, T_B, "m^")
plt.plot(L, T_ambient, "m^")
plt.plot(x_sol, y_sol, "g-", label="analytic solution")
plt.title("Tri-diagonal matrix algorithm for example 4.3 with %d nodes" % n)
plt.xlabel("Distance x (m)")
plt.ylabel("Temperature (Celsius)")
plt.legend()
plt.axis([x_min, x_max, y_min, y_max])
plt.show()
