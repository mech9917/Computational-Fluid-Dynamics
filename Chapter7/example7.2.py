# CFD Homework 20171842 김희성
# Programming Language/ver : Python 3.9
# Due: 08. December, 2021
# Example 7.2 Two-Dimensional Diffusion Problem

# Imports
import numpy as np
import pandas as pd

# Introduction
print('CFD Homework 20171842 김희성')
print('Programming Language/ver : Python 3.9')
print('Due: 08. December, 2021')
print('Example 7.2 Two-Dimensional Diffusion Problem')
print('A two-dimensional line-by-line application of TDMA. A two-dimensional plate of thickness 1cm is shown.')
print('The thermal conductivity of the plate material is k = 1000W/m/K.')
print('The west boundary receives a steady heat flux of 500kW/m^2 and the south and east boundaries are insulated.')
print('If the north boundary is maintained at a temperature of 100 deg Celsius,')
print('use uer defined uniform grid to calculate the steady state temperature distribution at nodes 1,2,3,4,...etc.')
print('------------------------------------------------------------------------------------------------------------')

# Constants
L_x = 0.3       # meter
L_y = 0.4       # meter
t = 0.01        # meter, thickness of plate
q_W = 500000    # Watt per meter square
T_N = 100       # deg Celsius
k = 1000        # Watt per meter Kelvin

# Input Variables
n_x = int(input('Type the # of nodes in x-direction: '))
n_y = int(input('Type the # of nodes in y-direction: '))
delta_x = L_x / n_x
delta_y = L_y / n_y
A_we = delta_y * t
A_sn = delta_x * t

# Initialize Variables
x = np.zeros(shape=(n_x, n_y))
y = np.zeros(shape=(n_x, n_y))
x_arr = np.array([])
y_arr = np.array([])
pointnum = [_ + 1 for _ in range(n_x * n_y)]
a_W = np.zeros(shape=(n_x, n_y))
a_W_arr = np.array([])
a_E = np.zeros(shape=(n_x, n_y))
a_E_arr = np.array([])
a_S = np.zeros(shape=(n_x, n_y))
a_S_arr = np.array([])
a_N = np.zeros(shape=(n_x, n_y))
a_N_arr = np.array([])
b_W = np.zeros(shape=(n_x, n_y))
b_E = np.zeros(shape=(n_x, n_y))
b_S = np.zeros(shape=(n_x, n_y))
b_N = np.zeros(shape=(n_x, n_y))
S_u = np.zeros(shape=(n_x, n_y))
S_u_arr = np.array([])
S_P = np.zeros(shape=(n_x, n_y))
a_P = np.zeros(shape=(n_x, n_y))
a_P_arr = np.array([])
M = np.zeros(shape=(n_x * n_y, n_x * n_y))
T_FVM = np.zeros(n_x * n_y)
A = np.zeros(shape=(n_x, n_y))
C = np.zeros(shape=(n_x, n_y))
Cprime = np.zeros(shape=(n_x, n_y))
T_TDMA_arr = np.array([])
T_TDMA_old = np.zeros(shape=(n_x, n_y+1))
T_TDMA_new = np.zeros(shape=(n_x, n_y+1))
iteration = 0
itererror_TDMA = 0

# Solve with Direct FVM Matrix Method
for i in range(n_x):
    for j in range(n_y):
        x[i, j] = delta_x / 2 + i * delta_x
        y[i, j] = delta_y / 2 + j * delta_y
        a_W[i, j] = k * A_we / delta_x
        a_E[i, j] = k * A_we / delta_x
        a_S[i, j] = k * A_sn / delta_y
        a_N[i, j] = k * A_sn / delta_y
        if i == 0:
            a_W[i, j] = 0
            b_W[i, j] = q_W * A_we
        elif i == n_x - 1:
            a_E[i, j] = 0
        if j == 0:
            a_S[i, j] = 0
        elif j == n_y - 1:
            a_N[i, j] = 0
            b_N[i, j] = 2 * k * A_sn / delta_y * T_N
            S_P[i, j] = -2 * k * A_sn / delta_y
        a_P[i, j] = a_W[i, j] + a_E[i, j] + a_S[i, j] + a_N[i, j] - S_P[i, j]
        S_u[i, j] = b_W[i, j] + b_E[i, j] + b_S[i, j] + b_N[i, j]
    x_arr = np.hstack((x_arr, x[i]))
    y_arr = np.hstack((y_arr, y[i]))
    a_W_arr = np.hstack((a_W_arr, a_W[i]))
    a_E_arr = np.hstack((a_E_arr, a_E[i]))
    a_S_arr = np.hstack((a_S_arr, a_S[i]))
    a_N_arr = np.hstack((a_N_arr, a_N[i]))
    a_P_arr = np.hstack((a_P_arr, a_P[i]))
    S_u_arr = np.hstack((S_u_arr, S_u[i]))

# Generate M Matrix
for i in range(n_x):
    for j in range(n_y):
        M[i * n_y + j, i * n_y + j] = a_P[i, j]
        if j != 0:
            M[i * n_y + j - 1, i * n_y + j] = -a_S[i, j]
        if j != n_y - 1:
            M[i * n_y + j + 1, i * n_y + j] = -a_N[i, j]
        if i != 0:
            M[(i - 1) * n_y + j, i * n_y + j] = -a_W[i, j]
        if i != n_x - 1:
            M[(i + 1) * n_y + j, i * n_y + j] = -a_E[i, j]

# Solve Matrix Algebra
T_FVM = np.linalg.solve(M, S_u_arr)
print('[FVM Tensor Matrix]')
print(M)

# Drawing Table
data1 = {'Point': pointnum, 'a_N': a_N_arr, 'a_S': a_S_arr, 'a_W': a_W_arr,
         'a_E': a_E_arr, 'a_P': a_P_arr, 'S_u': S_u_arr}
data2 = {'Point': pointnum, 'x': x_arr, 'y': y_arr, 'T': T_FVM}
Table1 = pd.DataFrame(data1)
Table2 = pd.DataFrame(data2)
print('[Table 7.3]')
print(Table1)
print('[Direct FVM Matrix Solution]')
print(Table2)


# Second Method
# Solve with TDMA Iterative Method
D = a_P
alpha = a_N
beta = a_S

# Solve until Numerical Solution Converges
while True:
    T_TDMA_arr = np.array([])
    x_arr = np.array([])
    y_arr = np.array([])
    itererror_TDMA = 0
    iteration += 1
    for i in range(n_x):
        for j in range(n_y):
            x[i, j] = delta_x / 2 + i * delta_x
            y[i, j] = delta_y / 2 + j * delta_y
            T_TDMA_old[i, j] = T_TDMA_new[i, j]
            if i == 0:
                C[i, j] = a_E[i, j] * T_TDMA_old[i+1, j] + S_u[i, j]
            elif i < n_x-1:
                C[i, j] = a_W[i, j] * T_TDMA_new[i-1, j] + a_E[i, j] * T_TDMA_old[i+1, j] + S_u[i, j]
            elif i == n_x-1:
                C[i, j] = a_W[i, j] * T_TDMA_new[i-1, j] + S_u[i, j]
            if j == 0:
                A[i, j] = alpha[i, j] / D[i, j]
                Cprime[i, j] = C[i, j] / D[i, j]
            else:
                A[i, j] = alpha[i, j] / (D[i, j] - beta[i, j] * A[i, j-1])
                Cprime[i, j] = (beta[i, j] * Cprime[i, j-1] + C[i, j]) / (D[i, j] - beta[i, j] * A[i, j-1])
        for j in range(n_y-1, -1, -1):
            T_TDMA_new[i, j] = A[i, j] * T_TDMA_new[i, j + 1] + Cprime[i, j]
            itererror_TDMA = itererror_TDMA + abs(T_TDMA_new[i, j] - T_TDMA_old[i, j])
        x_arr = np.hstack((x_arr, x[i]))
        y_arr = np.hstack((y_arr, y[i]))
        T_TDMA_arr = np.hstack((T_TDMA_arr, T_TDMA_new[i]))
        T_TDMA_arr = np.delete(T_TDMA_arr, (i + 1) * n_y)
    if iteration == 1:
        data = {'Point': pointnum, 'x': x_arr, 'y': y_arr, 'T': T_TDMA_arr}
        Table = pd.DataFrame(data)
        print('[TDMA Solution]')
        print('Iteration %d' % iteration)
        print(Table)
    # Define Conditions of Solution Convergence
    elif iteration % 2 != 0 and itererror_TDMA < 0.1:
        # Drawing Table
        data = {'Point': pointnum, 'x': x_arr, 'y': y_arr, 'T': T_TDMA_arr}
        Table = pd.DataFrame(data)
        print('Iteration %d' % iteration)
        print(Table)
        print('Converged solution after %d iterations' % iteration)
        break
    print('Iteration %d not converged' % iteration)
end = input('Press any key to end: ')
