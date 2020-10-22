# -*- coding: utf-8 -*-
"""
Created on Wed Oct 21 23:16:12 2020

@author: LilyHeAsamiko
"""
import numpy as np
import matplotlib.pyplot as plt
# Reflection coefficient for weak probe fields
Cl = 192
Cr = 257
Ck = 78
Cj = 29
w0 = 7.5*2*np.pi

for n in np.linspace(0.1,1,10):
    J = Cj*w0/4/Cr
    argtao = []
    for w in np.linspace(2*np.pi*6.5,2*np.pi*8.5,11):
        k = J/n
        gamma = k/10
        delta = w0-w
        tao = -1+k*(gamma/2-1j*delta)/(J**2-(1j*k/2+delta)*(1j*gamma/2+delta))
        argtao.append(np.arctan(np.imag(tao)/np.real(tao)))  
        
    plt.figure()
    plt.plot(np.linspace(6.5,8.5,11),argtao)
    plt.title('w/2*pi when J/k = {}'.format(n))
               
# Amplifier gain
# assume bin = np.cos(w/2*pi*t)-1j*np.sin(w/2*pi*t)
k = J/0.7
b = np.array([np.cos(w/2/np.pi)-1j*np.sin(w/2*np.pi),np.cos(w/2/np.pi)+1j*np.sin(w/2/np.pi),0,0]).T    
alpha0 = J**2 -(1j*k/2+delta)*(1j*gamma/2+delta)
Lj = 1.2/30
h = 6.62607015*10**(-34)
delta0 = 4150
U0 = -80*2*np.pi #(-80=-Ec/h/M**2 kHz)
count = 1
gS = []
gI = []

for w in np.linspace(2*np.pi*6.5,2*np.pi*8.5,22):
    k = J/0.7
    gamma = k/10
    delta = w0-w
    if count == 1:
        b = np.array([np.cos(w/2/np.pi)-1j*np.sin(w/2*np.pi),np.cos(w/2/np.pi)+1j*np.sin(w/2/np.pi),0,0]).T    
        alpha0 = J**2 -(1j*k/2+delta)*(1j*gamma/2+delta)
        alphaL = np.sqrt(k)*alpha0
        alphaR = -1j*J*alphaL
        S = np.array([[-1j*(-delta-1j*k/2), 0, -1j*J, 0],[0, 1j*(-delta-1j*k/2), 0, 1j*J],[-1j*J, 0, -1j*(-delta-1j*gamma/2), 0],[0, 1j*J, 0, 1j*(-delta-1j*gamma/2)]])
        G = np.linalg.inv(-1j*np.eye(4,4)*delta-S)
    else:
        UL = (1j*(delta+1j*k/2)*alphaL - 1j*J*alphaR + 1j*np.sqrt(k)*alphaL)/1j/alphaL**3
        UR = (1j*(delta+1j*gamma/2)*alphaR - 1j*J*alphaL)/1j/alphaR**3
        deltaLext = -delta - 1j*k/2 + 2*UL*alphaL**2
        deltaRext = -delta - 1j*gamma/2 + 2*UR*alphaR**2        
        S = np.array([[-1j*deltaLext, -1j*UL*alphaL**2, -1j*J, 0],[-1j*UL*alphaL**2, 1j*deltaLext, 0, 1j*J],[-1j*J, 0, -1j*deltaRext, -1j*UR*alphaR**2],[0, 1j*J, -1j*UR*alphaR**2, 1j*deltaRext]])
        G = np.linalg.inv(-1j*np.eye(4,4)*delta-S)
    count += 1
    gS.append(k*G[1,1]-1) 
    gI.append(k*G[1,2])     

wL = (2*np.pi*6.5-w0)/k
wR = (2*np.pi*8.5-w0)/k
plt.figure()
plt.plot(np.linspace(wL,wR,22),abs(np.array(gI)))
plt.title('|gI|^2 when J/k = {}'.format(delta/k))

plt.figure()
plt.plot(np.linspace(wL,wR,22),np.array(gS-min(gS)))
plt.title('|gs|^2 when J/k = {}'.format(delta/k))

plt.figure()
plt.plot(np.linspace(wL,wR,22),np.array(gS-min(gS))/abs(np.array(gI)))
plt.title('GS/I when J/k = {}'.format(delta/k))
