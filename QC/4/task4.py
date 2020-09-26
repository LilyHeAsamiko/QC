# -*- coding: utf-8 -*-
"""
Created on Thu Sep 24 14:31:58 2020

@author: LilyHeAsamiko
"""
import numpy as np
from scipy.stats import *
from random import *
import matplotlib.pyplot as plt

#solve the H = H1 + H2 + H3 +H4= a*XX + b*YY + c*ZZ + d*II
#try [-1/2, -1, 0, 1, 1/2].
param = np.array([-1/2, -1, 0, 1, 1/2])
H = np.zeros((4,4))
thetas = param
H[0,0] = 1
H[1,2] = -1
H[2,1] = -1
H[3,3] = 1
H0 = H
Theta = []
#thetas = []
HH = []
II = np.eye(4)
for i in range(len(param)):
    a = param[i]
    b = a
    theta0 = np.arcsin(-(1j+0.0001)/2/(a+0.0001)
    
#    thetas[i] = 1j*theta0
#    EXP = np.exp(-1j*theta0)
#     cos(theta0) = 
    c = [param[np.argmin(abs(1/(2* np.exp(1j*theta0)+np.exp(-1j*theta0))-param))]]
    d = param[np.argmin(abs(1+1j*np.sin(theta0)/np.cos(theta0)-1/(2+exp(-1j*theta0))-param))]
    XX = np.array([[np.cos(theta0), 0, 0, -1j*np.sin(theta0)],[0, np.cos(theta0), -1j*np.sin(theta0), 0],[0, -1j*np.sin(theta0), np.cos(theta0), 0],[-1j*np.sin(theta0), 0, 0, np.cos(theta0)]])
    YY = np.array([[np.cos(theta0), 0, 0, 1j*np.sin(theta0)],[0, np.cos(theta0), -1j*np.sin(1j*theta0), 0],[0, -1j*np.sin(theta0), np.cos(theta0), 0],[-1j*np.sin(theta0), 0, 0, np.cos(theta0)]])
#    ZZ = [[np.exp(1j*theta/2), 0, 0, 0],[0, np.exp(-1j*theta/2),0, 0],[0, 0, np.exp(-1j*theta/2), 0],[0, 0, 0, np.exp(1j*theta/2)]]
    ZZ = np.array([[np.exp(1j*theta0), 0, 0, 0],[0, np.exp(-1j*theta0),0, 0],[0, 0, np.exp(-1j*theta0), 0],[0, 0, 0, np.exp(1j*theta0)]])
    Ht = a*XX + b*YY + c*ZZ + d*II
    
    if np.mean(Ht - H - 0.5*H)>0:
        H0 = Ht
        HH.append(H0)
        Theta.append(theta0)
        

