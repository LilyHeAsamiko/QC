# -*- coding: utf-8 -*-
"""
Created on Thu Sep 24 00:25:27 2020

@author: LilyHeAsamiko
"""
import numpy as np
from scipy.stats import *
from random import *
import matplotlib.pyplot as plt
U1 = np.ones((8,8))
U2 = U1
R1 = np.zeros((8,8))
#R2 = R1
R2 = 1j*R1
CNOT = R1
Err =[]
Theta = []
#Theta1 = []
#Theta2 = []
#Theta3 = []
#Theta4 = []
Phi = []
Psi = []
errs = 64
phis = 0
psis = 0
theta1s = 0
theta2s = 0
theta3s = 0
theta4s = 0
theta1 = 0.0
theta2 = 0.0
theta3 = 0.0
theta4 = 0.0
for L in range(200):
    for theta1 in np.linspace(0,2*np.pi,10):
        for theta2 in np.linspace(0,2*np.pi,10):
            for theta3 in np.linspace(0,2*np.pi,10):
                for theta4 in np.linspace(0,2*np.pi,10):               
 #                   z = [0 0 1 1]
                    z = np.array([0, 0, 1, 1, 0, 0, 1, 1])
                    R1[0,0] = np.exp(-theta1*1j/2)
                    R1[1,1] = np.exp(theta1*1j/2)
                    R1[2,2] = np.exp(-theta2*1j/2)
                    R1[3,3] = np.exp(theta2*1j/2)
                    R1[4,4] = np.exp(-theta3*1j/2)
                    R1[5,5] = np.exp(theta3*1j/2)
                    R1[6,6] = np.exp(-theta4*1j/2)
                    R1[7,7] = np.exp(theta4*1j/2)
                    U1 = U1*R1
                    CNOT[0,0] = 1
                    CNOT[1,1] = 1
                    CNOT[2,3] = 1
                    CNOT[3,2] = 1
                    CNOT[4,4] = 1
                    CNOT[5,5] = 1
                    CNOT[6,7] = 1
                    CNOT[7,6] = 1                   
                    U1 = U1*CNOT
  #                  R2[0,0] = np.cos(theta1/2)
  #                  R2[1,1] = np.cos(theta1/2)
  #                  R2[0,1] = -np.sin(theta1/2)
  #                  R2[1,0] = np.sin(theta1/2)                
  #                  R2[2,2] = np.cos(theta2/2)
  #                  R2[3,3] = np.cos(theta2/2)
  #                  R2[2,3] = -np.sin(theta2/2)
                    #R2[3,2] = np.sin(theta2/2)  
  #                  R2[3,2] = np.sin(theta2/2)  
  #                  R2[4,4] = np.cos(theta3/2)
  #                  R2[5,5] = np.cos(theta3/2)
  #                  R2[4,5] = -np.sin(theta3/2)
  #                 R2[5,4] = np.sin(theta3/2)
  #                  R2[6,6] = np.cos(theta4/2)
  #                  R2[7,7] = np.cos(theta4/2)
  #                  R2[6,7] = -np.sin(theta4/2)
  #                  R2[7,6] = np.sin(theta4/2)
                    # For Bonus Questions: substitute Rx with Ry:
                    R2[0,0] = np.cos(theta1/2)
                    R2[1,1] = np.cos(theta1/2)
                    R2[0,1] = 1j*np.sin(theta1/2)
                    R2[1,0] = 1j*np.sin(theta1/2)                
                    R2[2,2] = np.cos(theta2/2)
                    R2[3,3] = np.cos(theta2/2)
                    R2[2,3] = 1j*np.sin(theta2/2)
                    #R2[3,2] = np.sin(theta2/2)  
                    R2[3,2] = 1j*np.sin(theta2/2)  
                    R2[4,4] = np.cos(theta3/2)
                    R2[5,5] = np.cos(theta3/2)
                    R2[4,5] = 1j*np.sin(theta3/2)
                    R2[5,4] = 1j*np.sin(theta3/2)
                    R2[6,6] = np.cos(theta4/2)
                    R2[7,7] = np.cos(theta4/2)
                    R2[6,7] = 1j*np.sin(theta4/2)
                    R2[7,6] = 1j*np.sin(theta4/2)
                    U2 = U2*R2
                    phi = (U1*z.T + U2*z.T)/2*z.T
                   # psi = (binom.pmf(k=range(64),n=10,p=0.5)).reshape(8,8)*z.T
#```````1      
                    psi = np.array(choices([0,1],[1000,1000],k=64)).reshape(8,8)*z.T
                    errs0 = phi - psi
                    if sum(sum(abs(errs0))) <= errs:
                        errs = sum(sum(abs(errs0)))
                        phis = phi
                        psis = psi
                        theta1s = theta1
                        theta2s = theta2
                        theta3s = theta3
                        theta4s = theta4
    Theta.append([theta1s,theta2s,theta3s,theta4s])
#    Theta1.append(theta1s)
#    Theta2.append(theta2s)
#    Theta3.append(theta3s)
#    Theta4.append(theta4s)                    
    Phi.append(phis)
    Psi.append(psis)    
    Err.append(errs)
                    
plt.figure()
plt.plot(Err)
plt.title('sum of error in 200 epochs(bonus)')  

plt.figure()
plt.plot(np.array(Theta)/np.pi*180)
plt.legend(['theta1','theta2','theta3','theta4'])
plt.title('Best Theta(digree) for each epoch in 200 epochs(bonus)')                                     
                    
plt.figure()
plt.plot(np.array(Theta))
plt.legend(['theta1','theta2','theta3','theta4'])
plt.title('Best Theta for each epoch in 200 epochs(bonus)')                                     
         

