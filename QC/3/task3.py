# -*- coding: utf-8 -*-
"""
Created on Thu Sep 24 13:17:43 2020

@author: LilyHeAsamiko
"""
import numpy as np

def Rz(theta):
    RZ = np.array([[np.exp(-1j*theta/2), 0],[0,np.exp(-1j*theta/2)]])
    return RZ

def Rx(theta):
    RX = np.array([[np.cos(theta/2),  -np.sin(theta/2)],[np.cos(theta/2), np.sin(theta/2)]])
    return RX

def Ry(theta):
    RY = np.array([[np.cos(theta/2),  1j*np.sin(theta/2)],[np.cos(theta/2), 1j*np.sin(theta/2)]])
    return RY

def H(theta):
    H = np.array([[1, 1],[1,-1]])/np.sqrt(2)
    return H

def I(theta):
    if np.size(theta)>1:
        I = np.eye(2*len(theta))
    else:
        I = np.eye(2)
    return I

def X(theta):
    X = np.array([[0,1],[1,0]])
    return X

def Y(theta):
    Y = np.array([[0,-1j],[1j,0]])
    return Y

def Z(theta):
    Z = np.array([[1,0],[0,-1]])
    return Z

def CX(theta):
    CX = np.array([1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1])
    return CX

def CZ(theta):
    CZ = np.array([1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,-1])
    return CZ

#test 
I(np.pi/2) + np.array(Rz(2*np.pi),int)

X(np.pi/2) +Rx(np.pi)
#X2 = [[0, 0, 0, 1],[0 ,0, 1, 0],[0, 1, 0, 0],[1, 0, 0, 0]]
