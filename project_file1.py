# -*- coding: utf-8 -*-
"""
Created on Sat Dec 11 14:19:57 2021

@author: K Raghavendiran
"""

import numpy as np
import math 

def INPUT(n):
    E = np.zeros((n,1))  
    for i in range(0, n):
        temp = float(input())
        E[i] += temp
    return E

#Engineering constants at plane
print("Enter E1,E2,G12,v12 in GPa:")
E = INPUT(4)


# Compliance matrix calculation
S = np.zeros((3,3))

S[0][0] += 1/E[0]
S[1][1] += 1/E[1]
S[2][2] += 1/E[2]
S[1][0] -= E[3]/E[0]
S[0][1] -= E[3]/E[0]



# Stiffness matrix calculation
Q = np.linalg.inv(S)


# Function for Transformation matrix calculation
def Theta(theta):
    m = math.cos((theta*math.pi)/180)
    n = math.sin((theta*math.pi)/180)
    T = np.zeros((3,3))
    T[0][0] = m*m
    T[0][1] = n*n
    T[0][2] = 2*m*n
    T[1][0] = n*n
    T[1][1] = T[0][0]
    T[1][2] = -2*m*n
    T[2][0] = -m*n
    T[2][1] = m*n
    T[2][2] = (m*m) - (n*n)
    return T


print("Enter multiply orientation (fibre orientation angle) one by one:")
theta = INPUT(4)


# Calculation of Transformation matrix for each orientation angle  
T1 = Theta(theta[0])
inv_T1 = np.linalg.inv(T1)
T2 = Theta(theta[1])
inv_T2 = np.linalg.inv(T2)
T3 = Theta(theta[2])
inv_T3 = np.linalg.inv(T3)
T4 = Theta(theta[3])
inv_T4 = np.linalg.inv(T4)


#Function for Calculation of stiffness matrix
def stiffness(T,Q,inv_T):
    Q_xy = np.dot(np.dot(inv_T,Q),T)
    return Q_xy


#Calculation of Stiffness matrix 
Q_xy_1 = stiffness(T1,Q,inv_T1)

Q_xy_2 = stiffness(T2,Q,inv_T2)

Q_xy_3 = stiffness(T3,Q,inv_T3)
for i in range(3):
    k = Q_xy_3[i][2]
    Q_xy_3[i][2] = (k/2)
    
Q_xy_4 = stiffness(T4,Q,inv_T4)
for i in range(3):
    k = Q_xy_4[i][2]
    Q_xy_4[i][2] = (k/2)

print("Enter ply thickness ratios:")
ratio = INPUT(4)
h = float(input("Enter total thickness of the laminate:"))


# Total Stiffness matrix
K = np.zeros((3,3)) 
K += np.multiply(Q_xy_1,(h*ratio[0])) + np.multiply(Q_xy_2,(h*ratio[1])) + np.multiply(Q_xy_3,(h*ratio[2])) + np.multiply(Q_xy_4,(h*ratio[3]))

def eng_const_lamina(S):
    E_xx = 1/S[0][0]
    print("E_xx:",E_xx)
    E_yy = 1/S[1][1]
    print("E_yy:",E_yy)
    G_xy = 1/S[2][2]
    print("G_xy:",G_xy)
    v_xy = -S[1][0]/S[0][0]
    print("v_xy:",v_xy)
    n_sx = S[0][0]/S[2][2]
    print("n_sx:",n_sx)
    n_ys = S[2][1]/S[2][2]
    print('n_ys:',n_ys)

A = np.linalg.inv(K)

E_x = 1/(h*A[0][0])
print("E_x:",E_x)
E_y = 1/(h*A[1][1])
print("E_y:",E_y)
G_xy = 1/(h*A[2][2])
print("G_xy:",G_xy)
v_xy = -A[1][0]/A[0][0]
print("v_xy:",v_xy)
v_yx = -A[0][1]/A[1][1]
print("v_yx:",v_yx)
n_sx = A[0][2]/A[2][2]
print("n_sx:",n_sx)
n_xs = A[2][0]/A[0][0]
print("n_xs:",n_xs)
n_ys = A[2][1]/A[1][1]
print('n_ys:',n_ys)
n_sy = A[1][2]/A[2][2]
print('n_sy:',n_sy)
