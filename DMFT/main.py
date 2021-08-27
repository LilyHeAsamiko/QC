'''

                            Online Python Compiler.
                Code, Compile, Run and Debug python program online.
Write your code in this editor and press "Run" button to execute it.

'''

#print("Hello World")
import numpy as np


def RXL(theta,XL):
    RXL = np.exp(-1j*theta*XL/2)
    return RXL
    
def RZL(theta,ZL):
    RZL = np.exp(-1j*theta*ZL/2)
    return RZL

PX = np.array([[0,1,0,0],[1,0,0,0],[0,0,0,1],[0,0,1,0]])
PY = np.array([[0,-1j,0,0],[1j,0,0,0],[0,0,0,-1j],[0,0,1j,0]])
PZ = np.array([[1,0,0,0],[0,-1,0,0],[0,0,1,0],[0,0,0,-1]])
V12 = np.array([[1,0,0,0],[0,-1,0,0],[0,0,1,0],[0,0,0,-1]])#ZL 
V23 = np.array([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]]) 


nm = 10**(-9) #m
h = 4.136/10**(-15) #eV*s
hbar = h/np.pi/2
Lambda = 514*nm
m = 1/2
ms = 1000
# DFS qubit
zero = np.array([[0,0,0,0],[0,-1,0,0],[0,0,0,0],[0,0,0,-1]])/2
one = (np.array([[2,0,0,0],[0,2,0,0],[0,0,2,0],[0,0,0,2]])-np.array([[0,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]))/2/np.sqrt(3)

E_R = hbar**2*(2*np.pi)**2/2/m/Lambda**2
J0 = 0.033*E_R
#initialization
J2 = 0
J1 = J0
J3 = J0
J4 = 0
J = np.array([J1,J2,J3,J4])
U = 75*J0
H0 = -J0**2/U*(PX+PY+PZ)
t0 = np.pi*hbar*U/8/min(J+0.00001)**2
state = np.array([[1j,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1j]])/2
state = np.array([[np.exp(-1j*np.pi/4)*1j,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,-1j*np.exp(-1j*np.pi/4*1j)]])/2
state = np.array([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]])/2
print('initialization need'+str(t0*ms)+'ms')

#RZ
J2 = 0
J3 = 0
J = np.array([J1,J2,J3,J4])
ZL = -V12
HZ = -J**2/U*ZL
#print(HZ)

# theta = np.pi example
theta = np.pi
RZL = RZL(theta,ZL)
tZ = theta*U*hbar/4/min(J+0.00001)**2 #8.7 ms
print('Rotation on Z axis need '+str(tZ*ms)+'ms')
#special case: interaction free
tZ1 = theta*hbar/2/min(J+0.00001) #0.23 ms
print('Under interaction free case, Rotation on Z axis need'+str(tZ1*ms)+'ms')
TZ1 = np.linspace(0,tZ1,100)
#print(TZ1)
print(HZ)
print(np.max(HZ,0))
fZ = abs(np.exp(-1j*sum(np.max(HZ,0))*TZ1/hbar))**2

import matplotlib.pyplot as plt
plt.figure()
plt.plot(fZ)

#RX
J1= J0
J2 = np.sqrt(2)*J0
J = np.array([J1,J2,J3,J4])
XL = (V12+2*V23)/np.sqrt(3)
HX = -2*np.sqrt(3)*J**2/U*XL*4
# theta = -np.pi example
tX1 = -theta*U*hbar/np.sqrt(3)/4/min(J+0.00001)**2 #5 ms
TX1 = np.linspae(0,tx1,100)
fX = np.exp(-1j*np.max(HX,0)*TX1/hb*hbar)
U = 100*J

print('Rotation on X axis need'+str(tx1*ms)+'ms')



