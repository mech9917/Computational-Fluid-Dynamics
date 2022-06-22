# CFD Homework 20171842 김희성
# Programming Language/ver : Python 3.9
# Due: 24. November, 2021
# Example 5.1 1-D convection-diffusion problem with phi_0=1, phi_L=0

# Imports
import numpy as np
import math
import pandas as pd
import matplotlib.pyplot as plt

# Introduction
print('CFD Homework 20171842 김희성')
print('Programming Language/ver : Python 3.9')
print('Due: 24. November, 2021')
print('Example 5.1 1-D convection-diffusion problem with phi_0=1, phi_L=0')
print('-----------------------------------------------------------------------')

# Constants
rho = 1         # kilogram per meter cubic
u = 0.1         # meter per seconds
F = rho * u
L = 1           # meter
Gamma = 0.1
phi_0 = 1       # deg Celsius
phi_L = 0       # deg Celsius

# Input Variables
n = int(input('Type the # of CVs: '))    # Number of Control Volume
delta_x = L / n
D = Gamma / delta_x

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
        sol[i] = (phi_L - phi_0) * (math.exp(rho * u * x[i] / Gamma) - 1) / (math.exp(rho * u * L / Gamma) - 1) + phi_0
        a_W[i] = 0
        a_E[i] = D - F / 2
        S_u = np.vstack((S_u, (2 * D + F) * phi_0))
        S_u_arr[i] = (2 * D + F) * phi_0
        S_P[i] = -(2 * D + F)
        a_P[i] = a_W[i] + a_E[i] - S_P[i]
        m = np.zeros(n)
        m[i] = a_P[i]
        m[i + 1] = -a_E[i]
        M = np.vstack((M, m))
    elif i < n-1:             # Interior Nodes
        x[i] = x[0] + i * delta_x
        sol[i] = (phi_L - phi_0) * (math.exp(rho * u * x[i] / Gamma) - 1) / (math.exp(rho * u * L / Gamma) - 1) + phi_0
        a_W[i] = D + F / 2
        a_E[i] = D - F / 2
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
        sol[i] = (phi_L - phi_0) * (math.exp(rho * u * x[i] / Gamma) - 1) / (math.exp(rho * u * L / Gamma) - 1) + phi_0
        a_W[i] = D + F / 2
        a_E[i] = 0
        S_u = np.vstack((S_u, (2 * D - F) * phi_L))
        S_u_arr[i] = (2 * D - F) * phi_L
        S_P[i] = -(2 * D - F)
        a_P[i] = a_W[i] + a_E[i] - S_P[i]
        m = np.zeros(n)
        m[i - 1] = -a_W[i]
        m[i] = a_P[i]
        M = np.vstack((M, m))
    else:                       # if something went wrong...
        break

# Solution
phi = np.linalg.solve(M, S_u)
phi_arr = np.zeros(n)
for j in range(n):
    phi_arr[j] = phi[j][0]
    error[j] = abs(phi_arr[j] - sol[j]) / sol[j] * 100

# Drawing Table
nodes = [_ + 1 for _ in range(n)]
data1 = {'Node': nodes, 'a_W': a_W, 'a_E': a_E, 'S_u': S_u_arr, 'S_P': S_P, 'a_P=a_W+a_E-S_P': a_P}
data2 = {'Node': nodes, 'x(m)': x, 'Finite volume solution': phi_arr, 'Analytical solution': sol, 'Percentage error': error}
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
y_max = 1
x_sol = np.linspace(x_min, x_max, 1000)
y_sol = np.zeros(1000)
for space in range(len(y_sol)):
    y_sol[space] = (phi_L - phi_0) * (math.exp(rho * u * x_sol[space] / Gamma) - 1) / (math.exp(rho * u * L / Gamma) - 1) + phi_0


for m in range(n):
    plt.plot(x[m], phi[m][0], "b*")
    # style : :(dot),-(line),--(dash line),^(tri),*(bold)
    # color : R, G, B, C, M, Y, K
plt.plot(0, phi_0, "m^")
plt.plot(L, phi_L, "m^")
plt.plot(x_sol, y_sol, "g-", label = "Exact Solution")
plt.title("1-D convection-diffusion")
plt.xlabel("Distance (m)")
plt.ylabel("phi")
plt.legend()
plt.axis([x_min, x_max, y_min, y_max])
plt.show()
