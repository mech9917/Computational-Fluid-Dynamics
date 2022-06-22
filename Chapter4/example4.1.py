# CFD Homework 20171842 김희성
# Programming Language/ver : Python 3.9
# Due: 17. November, 2021
# Example 4.1 Source-free heat conduction in an insulated rod

print('CFD Homework 20171842 김희성')
print('Programming Language/ver : Python 3.9')
print('Due: 17. November, 2021')
print('Example 4.1 Source-free heat conduction in an insulated rod')
print('-----------------------------------------------------------------------')

# Imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Constants
k = 1000        # Watt per meter Kelvin
A = 0.01        # meter square
L = 0.5         # meter
T_A = 100       # deg Celsius
T_B = 500       # deg Celsius

# Input Variables
n = int(input('Type the # of CVs: '))    # Number of Control Volume
delta_x = L / n

# Initialize Variables
x = np.zeros(n)
sol = np.zeros(n)
a_W = np.zeros(n)
a_E = np.zeros(n)
S_u = np.empty((0, 1), int)
S_u_arr = np.zeros(n)
S_P = np.zeros(n)
a_P = np.zeros(n)
M = np.empty((0, n), int)
error = np.zeros(n)

for i in range(n):
    if i == 0:                # Exterior Node (Node 1)
        x[i] = delta_x / 2
        sol[i] = (T_B - T_A) / L * x[i] + T_A
        a_W[i] = 0
        a_E[i] = k * A / delta_x
        S_u = np.vstack((S_u, 2 * k * A * T_A / delta_x))
        S_u_arr[i] = 2 * k * A * T_A / delta_x
        S_P[i] = -2 * k * A / delta_x
        a_P[i] = a_W[i] + a_E[i] - S_P[i]
        m = np.zeros(n)
        m[i] = a_P[i]
        m[i + 1] = -a_E[i]
        M = np.vstack((M, m))
    elif i < n-1:             # Interior Nodes
        x[i] = x[0] + i * delta_x
        sol[i] = (T_B - T_A) / L * x[i] + T_A
        a_W[i] = k * A / delta_x
        a_E[i] = k * A / delta_x
        S_u = np.vstack((S_u, 0))
        S_u_arr[i] = 0
        S_P[i] = 0
        a_P[i] = a_W[i] + a_E[i] - S_P[i]
        m = np.zeros(n)
        m[i-1] = -a_W[i]
        m[i] = a_P[i]
        m[i + 1] = -a_E[i]
        M = np.vstack((M, m))
    elif i == n-1:             # Exterior Node (Last Node)
        x[i] = x[0] + i * delta_x
        sol[i] = (T_B - T_A) / L * x[i] + T_A
        a_W[i] = k * A / delta_x
        a_E[i] = 0
        S_u = np.vstack((S_u, 2 * k * A * T_B / delta_x))
        S_u_arr[i] = 2 * k * A * T_A / delta_x
        S_P[i] = -2 * k * A / delta_x
        a_P[i] = a_W[i] + a_E[i] - S_P[i]
        m = np.zeros(n)
        m[i - 1] = -a_W[i]
        m[i] = a_P[i]
        M = np.vstack((M, m))
    else:                       # if something went wrong...
        break

# Solution
T = np.linalg.solve(M, S_u)
T_arr = np.zeros(n)
for j in range(n):
    T_arr[j] = T[j][0]
    error[j] = abs(T_arr[j] - sol[j]) / sol[j] * 100

# Drawing Table
nodes = [_ + 1 for _ in range(n)]
data1 = {'Node': nodes, 'a_W': a_W, 'a_E': a_E, 'S_u': S_u_arr, 'S_P': S_P, 'a_P=a_W+a_E-S_P': a_P}
data2 = {'Node': nodes, 'x(m)': x, 'Finite volume solution': T_arr, 'Exact solution': sol, 'Percentage error': error}
Table1 = pd.DataFrame(data1)
Table2 = pd.DataFrame(data2)
print('[Table 4.1]')
print(Table1)
print('-----------------------------------------------------------------------')
print('[Comparison with the analytical solution]')
print(Table2)

# Drawing Graph
x_min = 0
x_max = L
y_min = 0
y_max = 600
x_sol = np.linspace(x_min, x_max, 1000)
y_sol = (T_B - T_A) / L * x_sol + T_A

for m in range(n):
    plt.plot(x[m], T[m][0], "b*")
    # style : :(dot),-(line),--(dash line),^(tri),*(bold)
    # color : R, G, B, C, M, Y, K
plt.plot(0, T_A, "m^")
plt.plot(L, T_B, "m^")
plt.plot(x_sol, y_sol, "g-", label = "analytic solution")
plt.title("1-D Source-Free Heat Conduction")
plt.xlabel("Distance x (m)")
plt.ylabel("Temperature (Celsius)")
plt.legend()
plt.axis([x_min, x_max, y_min, y_max])
plt.show()
