# -*- coding: utf-8 -*-
"""
Created on Mon Dec  6 23:00:31 2021

@author: K Raghavendiran
"""
#import modules 
import numpy as np
import math 


def INPUT(n):
    E = np.zeros((n,1))  
    
    for i in range(0, n):
        temp = float(input())
        E[i] += temp
    return E


print("Enter E11,E22,E33 in GPa:")
E = INPUT(3)
print("ENTER G23,G31,G12 in GPa:")
G = INPUT(3)
print("ENTER v12,v13,v23:")
v = INPUT(3)

'''
E = [200,10,8.5]
G = [12.5,6.5,8.5]
v = [0.3,0,0]
'''



# Compliance matrix
S = np.zeros((6,6))

for i in range(3):
    if E[i]!= 0:
        S[i][i] += 1/E[i]
    else:
        S[i][i] += 0
    
for i in range(3):
    if G[i]!=0:
        S[i+3][i+3] += 1/G[i]
    else:
        S[i+3][i+3] += 0
    
S[0][1] -= v[0]/E[0]
S[0][2] -= v[1]/E[0]
S[1][2] -= v[2]/E[1]

S[1][0] += S[0][1]
S[2][0] += S[0][2]
S[2][1] += S[1][2]



# Stiffness matrix
Q = np.linalg.inv(S)
            


'''
# Finding Strain
print("Enter the stresses:")
stress = INPUT(6)

strain = np.dot(S,stress)  

#Finding stress
print("Enter the strain:")
strain = INPUT(6)
stress = np.dot(Q,strain) 
      
'''  
stress = np.array([15,10,0,0,0,23])
strain = np.dot(S,stress)



#Transformation Matrix
theta = int(input("enter Fibre orientation in Degrees:"))
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
inv_T = np.linalg.inv(T)


sigma1 = np.zeros((3,1))
sigma1[0] += stress[0]
sigma1[1] += stress[1]
sigma1[2] += stress[5]

strain1 = np.zeros((3,1))
strain1[0] += strain[0]
strain1[1] += strain[1]
strain1[2] += strain[5]

sigma_xy = np.dot(inv_T,sigma1)
print("Stress in reference axes:",sigma_xy)
strain_xy =np.dot(inv_T,strain1)
print("Strain in reference axes:",strain_xy)


# 3*3 Stiffness matrix
Q_1 = np.zeros((3,3))
for i in range(2):
    for j in range(2):
        Q_1[i][j] += Q[i][j]
        
Q_1[2][2] += Q[5][5]

Q_xy = np.dot(np.dot(inv_T,Q_1),T)
  
S_1 = np.zeros((3,3))
for i in range(2):
    for j in range(2):
        S_1[i][j] += S[i][j]
        
S_1[2][2] += S[5][5]

S_xy = np.dot(np.dot(inv_T,S_1),T)

print("Engineering constants after transformations:")
E_xx = 1/S_xy[0][0]
print("E_xx:",E_xx)
E_yy = 1/S_xy[1][1]
print("E_yy:",E_yy)
G_xy = 1/S_xy[2][2]
print("G_xy:",G_xy)
v_xy = -S_xy[1][0]/S_xy[0][0]
print("v_xy:",v_xy)
n_sx = S_xy[0][0]/S_xy[2][2]
print("n_sx:",n_sx)
n_ys = S_xy[2][1]/S_xy[2][2]
print('n_ys:',n_ys)
