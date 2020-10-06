# -*- coding: utf-8 -*-
"""
Created on Sat Oct  3 01:44:37 2020

@author: LilyHeAsamiko
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy as sci
from scipy import *
import scipy.stats as stats

def PX = X:
    PX = [[0,1],[1,0]];
    return PX
    
def PY = Y:
    PY = [[0,-1j],[1j,0]];
    return PY
    
def PZ = Z:
    PZ = [[1,0], [0,-1]];
    return PZ

def PXX = XX:
    PXX = [[1,0,0,-1j],[0,1,-1j,0],[0,-1j,1,0],[-1j,0,0,1]];
    return PXX
    
def PYY = YY:
    PYY = [[1,0,0,1j],[0,1,-1j,0],[0,-1j,1,0],[1j,0,0,1]];
    return PYY
    
def PZZ = ZZ:
    PZZ = [[1,0,0,0], [0,-1,0,0],[0,0,-1,0],[0,0,0,1]];
    return PZZ

def C_A_I_qr = operatorSC(PX,PY,PZ,A,I,q,r):
    C_A_I_qr = 1j/2
    C_A_I_qrR = 1j/2
    for M in range(q+1,A):  #M = range(q+1,A)
        if r == 'Re':
            C_A_I_qrR *= np.array(PZ)*(np.array(PY).T*np.array(PX)-np.array(PX)*np.array(PY))
            return C_A_I_qrR
        else:
            C_A_I_qr *= np.array(PZ)*(np.array(PY).T*np.array(PX)+np.array(PX)*np.array(PY))
            return C_A_I_qr

def C_AB_IJ_qr = operatorDC(PXX,PYY,PZZ,A,B,I,J,q,r):
    PX = np.array(PX)
    PY = np.array(PY)
    PZ = np.array(PX)    
    C_AB_IJ_qr = np.ones((4,4))*1j/2;
    C_AB_IJ_qrR = np.ones((4,4))*1j/2;
    for P in range(B+1,A):
        for Q in range(J+1,I): #range(J+1,I)
#            if r == 'Re':
                C_AB_IJ_qrR = C_AB_IJ_qrR*np.array(PZZ)*np.array(PZZ)*np.array(PYY)*np.array(PYY)*np.array(PYY)*np.array(PXX)+np.array(PXX)*np.array(PYY)*np.array(PXX)*np.array(PXX)+np.array(PYY)*np.array(PYY)*np.array(PXX)*np.array(PYY)+np.array(PYY)*np.array(PXX)*np.array(PXX)*np.array(PXX)-np.array(PYY).T*np.array(PXX)*np.array(PYY)*np.array(PYY)-np.array(PXX).T*np.array(PXX)*np.array(PXX)*np.array(PYY)-np.array(PXX).T*np.array(PYY)*np.array(PYY)*np.array(PYY)-np.array(PXX).T*np.array(PXX)*np.array(PYY)*np.array(PXX);
#                return C_AB_IJ_qrR
#            elif r == 'Im':
                C_AB_IJ_qr = C_AB_IJ_qr*np.array(PZZ)*np.array(PZZ)*np.array(PYY)*np.array(PYY)*np.array(PYY)*np.array(PYY)+np.array(PXX)*np.array(PXX)*np.array(PXX)*np.array(PXX)+np.array(PXX)*np.array(PYY)*np.array(PXX)*np.array(PYY)+np.array(PYY)*np.array(PXX)*np.array(PXX)*np.array(PXX)+np.array(PXX)*np.array(PYY)*np.array(PYY)*np.array(PXX)-np.array(PYY).T*np.array(PXX)*np.array(PXX)*np.array(PYY)-np.array(PYY).T*np.array(PYY)*np.array(PXX)*np.array(PXX)-np.array(PXX).T*np.array(PYY)*np.array(PYY)*np.array(PYY);
 #               return C_AB_IJ_qr

#intergral
def H = operatorH(sd):
    if sd == 's':
        Hs = [[0,0],[1,0]]
        return Hs 
    elif sd =='d':
        Hd = [[0,0,0,0],[1,0,0,0],[0,0,np.sqrt(2),0],[0,0,0,np.sqrt(3)]];   
        return Hd 

def G = operatorG(sd):
    if sd =='s':
        Gs = [[0,1],[0,0]];
        return Gs 
    elif sd =='d':
        Gd = [[0,1,0,0],[0,0,np.sqrt(2),0],[0,0,0,np.sqrt(3)],[0,0,0,0]];
        return Gd

# compare phase 0 and pi
#p3, table I

L = 2
d = 0.75 #(pA)
k = np.pi/2.5/d
#R1 = [[-[0,2,1],[2,0,1],-[1，3，2],[3，1，2]],[[0,2,1],[2,0,1],[1,3,2],[3,1,2]],[[-0,2,1],[0,2,1],[-1,2,3],[1,2,3]],[[0,2,1],[2,0,1],[1,3,2],[3,1,2]]]
#R2 = [[-[0,2,3],[0,0,0], [0,1，2]，[0,1,3]],[-[3,2,0],[2,3,0],[3,2,1],[3,2,1]]]
#Alternstively, operator projection in U instead
R = 1
As = 0
Ad = 0
HR = []
HT = []
M0 = []
M1 = []
M2 = []
M3 = []
H = []
Res = []
U = []
M = []
#A = 16
#B = 9
Magnitude = []
Phase = []
Intercept = []
STATS = []
#I= 16
#J = 9
for n in range(1,101):
#    for q = 1: 8  #0-7
    for theta in [-np.pi/2,0, np.pi/2, np.pi]:
        As = As + 1/np.sqrt(L)*sum(np.exp(1j*q*k*R)*C_A_I_qr+np.reshape(np.random.randn(4),(2,2))) #operatorSC(PX,PY,PZ,A,I,q,'Im')
        Asconj = -(As + 1/np.sqrt(L)*sum(np.exp(1j*q*k*R)*C_A_I_qrR+np.reshape(np.random.randn(4),(2,2)))) #operatorSC(PX,PY,PZ,A,I,q,'Re')
        Ad = Ad + 1/np.sqrt(L)*sum(np.exp(1j*q*k*R)*C_AB_IJ_qr+np.reshape(np.random.randn(16),(4,4))) #operatorDC(PX,PY,PZ,A,I,q,'Im')
        Adconj = -Ad - 1/np.sqrt(L)*sum(np.exp(1j*q*k*R)*C_AB_IJ_qrR-np.reshape(np.random.randn(16),(4,4))) #operatorDC(PX,PY,PZ,A,I,q,'Re')
        HR.append(sum(sum(np.reshape(Hs,(2,2))*As*Asconj)) + 0.5*sum(sum(np.reshape(Gs,(2,2))*As*Asconj)))
        HT.append(sum(sum(np.reshape(Hd,(4,4))*Ad*Adconj)) + 0.5*sum(sum(np.reshape(Gd,(4,4))*Ad*Adconj)))

        M.append(np.reshape(Hs,(2,2))*As*Asconj+ 0.5*np.reshape(Gs,(2,2))*As*Asconj)                
    # toatal
    U.append(np.cos(theta/2)*np.eye(np.shape(M[n-1:n+3])[0])-1j*np.sin(theta/2)) 
    Etemp = np.exp(-np.array(U[n-1])*np.reshape(M[n-1:n+3],(4,4))*np.conj(np.array(U)))
#    m0 = Etemp/k*(-0.5*1j*R1)*np.array(Hd)*Ad*Adconj 
#    m1 = Etemp/k*(-0.25*1j*R2)*np.array(Hd)*Ad*Adconj 
#    m2 = Etemp/k*(-0.5*R1)*np.array(Hd)*Ad*Adconj 
#    m3 = Etemp/k*(-0.25*R2)*np.array(Hd)*Ad*Adconj 
    m0 = M[0]
    m1 = M[1]
    m2 = M[2]
    m3 = M[3]
    #HR = np.array(HR)
    Etemp[Etemp == 0] = 0.000000000001    
    sigma = np.log(Etemp)
    HT.append(sum(np.reshape(Hd,(4,4))*Ad*Adconj) + 0.5*sum(np.reshape(Gd,(4,4))*Ad*Adconj))
    #HT = np.array(HT)
    Hext = np.exp(-sigma)*np.array(HT[-1])*np.exp(sigma)        
#   Hext = exp(sum(-2.*i.*PY.*PZ.*PY+2.*i.*PY.*PZ.*PX,[1,2])+sum(i.*PX.*PX.*PX.*PX+i.*PY.*PY.*PY+randn(4).*PY-i.*PX.*PX.*PY.*PY,'all'))   
  
    M0 =  np.reshape(Etemp[0,0,:],(2,2))/k*(-0.25*1j)*np.array(Hs)*As*Asconj 
    M1 =  np.reshape(Etemp[0,1,:],(2,2))/k*(0.25*1j)*np.array(Hs)*As*Asconj
    M2 =  np.reshape(Etemp[0,2,:],(2,2))/k*(-0.5)*np.array(Hs)*As*Asconj
    M3 =  np.reshape(Etemp[0,3,:],(2,2))/k*(0.5)*np.array(Hs)*As*Asconj
    # 
    # 0, pi, -pi/2, pi/2, optimization with magnitude A, phase B, intercept C
    A = 0.5*np.sqrt((M1-M3)**2+(M2-M0))**2
    B = np.arctan(2*np.linspace((M1-M3),(M2-M0),np.size(A))) 
    C = (M1+M3)/2
    A0 = 0.5*np.sqrt(((m1-m3)**2+(m2-m0))**2)
    B0 = np.arctan(2*np.linspace((m1-m3),(m2-m0),np.size(A)))
    C0 = (m1+m3)/2 
    res0 = A0 - B0 + C0
    res = abs(A-A0) - abs(B-B0) + (C-C0) 
    if sum(sum(sum(abs(res)-abs(res0)))>0):
        Res.append(res)
        Magnitude.append(A)
        Phase.append(B)
        Intercept.append(C)
#        [h,p,ci,stats] = 
#        [t1,t2]=sci.stats.ttest_ind(res,res0,variance = False)
        [t1,t2]=stats.ttest_ind(res,res0,equal_var = False)
        STATS.append([t1,t2])

plt.figure()
plt.plot(mean(abs(np.array(Res))[np.argmax(np.array(Res)>0)],2)/mean(mean(abs(np.array(Res))[np.argmax(np.array(Res)>0)],2)).flatten())
plt.title('Residule')
plt.figure()
plt.plot(abs(np.array(Magnitude))[np.array(Magnitude)>0]/mean(abs(np.array(Magnitude))[np.array(Magnitude)>0]).flatten())
plt.title('Magnitude')
plt.figure()
plt.plot(abs(np.array(Phase))[np.array(Phase)>0]/mean(abs(np.array(Phase))[np.array(Phase)>0]).flatten())
plt.title('Phase')
plt.figure()
plt.plot(abs(np.array(Intercept))[np.array(Intercept)>0]/mean(abs(np.array(Intercept))[np.array(Intercept)>0]).flatten())
plt.title('Intercept')
#A = H 

#convergency
Mag_mean = []
Mag_std = []
Mag_conv = abs(np.array(Magnitude))[np.array(Magnitude)>0]/mean(abs(np.array(Magnitude))[np.array(Magnitude)>0]).flatten()
for i in range(len(Mag_conv)):
    Mag_mean.append(np.mean(Mag_conv[0:i]))
    Mag_std.append(np.std(Mag_conv[0:i]))
plt.figure()
plt.errorbar(x = range(len(Mag_conv)),y = Mag_mean, yerr = Mag_std)
plt.title('Magnitude convergence')

Phase_mean = []
Phase_std = []
Phase_conv = abs(np.array(Phase))[np.array(Phase)>0]/mean(abs(np.array(Phase))[np.array(Phase)>0]).flatten()
for i in range(len(Phase_conv)):
    Phase_mean.append(np.mean(Phase_conv[0:i]))
    Phase_std.append(np.std(Phase_conv[0:i]))
plt.figure()
plt.errorbar(x = range(len(Phase_conv)), y = Phase_mean, yerr = Phase_std)
plt.title('Phase convergence')

