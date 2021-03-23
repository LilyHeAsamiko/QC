# -*- coding: utf-8 -*-
"""
Created on Mon Mar 22 07:37:05 2021

@author: LilyHeAsamiko
"""
import numpy as np
import matplotlib.pyplot as plt
#with the parameters in TPN: 
k5 = [6.15,7.75,10.25,11.25,12.65,13.25]
NI = 4096

gGABAE = 0.0006681
gGABAI = 0.0005120
gAMPAE = 0.0001905
gAMPAI = 0.0001460
gNMDAE = np.array(k5,dtype = float)/NI*0.001
gNMDAI = 60/NI*0.001
#in TNN: 
gGABAE = 0.0006681*1.5 
gGABAI = 0.000512*1.5
gAMPAE = 0.0001905*4 
gAMPAI = 0.0001460*4 
gNMDAE = 60/NI*0.001 
gNMDAI = np.array(k5,dtype = float)/NI*0.001
k5 = [6.15,7.75,10.25,11.25,12.65,13.25] 
NI = 4096
#assuming staring from inhibitory population, the synaptic reversal potentials of excitatory being 0. According to the experiment data of synaptic conductance g between pyramidal and interneuron population[]
VGABA = 0
Iaff = -0.5

dt = 1
Vm = -60
Cm = 0.5
gL = 0.025
VL = -70
sigma = 11.25
J = 3.62
W = np.ones((NI,1))
dtheta = 360/NI
#W[range(int(NI/2))] = -J +np.exp(-(np.linspace(0,180,NI/2))**2/2/sigma)
W[range(int(NI/2)),0] = J +(J+J/NI)*np.exp(-(np.linspace(0,180,NI/2))**2/2/sigma)
W[range(int(NI/2),int(NI)),0] = 360/dtheta - W[range(int(NI/2)),0]
I = np.ones((700,len(k5)))
Iaff = np.ones((700,len(k5)))
IGABA = np.ones((700,len(k5)))
V = Vm*np.ones((700,len(k5)))
t = 0
T2 = 10
I1 = -0.3
i = 0
for K5 in k5:
    for v in np.linspace(-90, 0, 700):
        I[t,i] = gGABAI*(Vm-VGABA)+gAMPAI*Vm+gNMDAI[i]*Vm/(1+1*np.exp(-0.062/3.57))
        # in exitatory population:
        Iaff[t,i] += 1/100
        I1 += 0.7/0.7*1/100
        IGABA[t,i] = I1
        dVm = (-gL*(Vm-VL)-IGABA[t,i])/Cm
        Vm += dVm 
        V[int(t),i] = Vm
        Vm_thr = -50
        R = np.log(I1/(Vm_thr/VL))*np.ones((NI,1))
        R[np.isnan(R)] = 0
        if I1 >0:
            kI2 = (-IGABA[t,i]+sum(W*R)-np.log(I1/(Vm_thr/VL)))/T2
        else:
            kI2 = (-IGABA[t,i]+sum(W*R))/T2
        Gsc = -(I1 + T2*kI2-sum(W*R))/R
        IGABA[int(t),i] = sum(-T2*kI2+sum(W*R)-Gsc*R+Iaff[t,i])
        if Iaff[int(t),i] > 0.4:
            Iaff[int(t),i] -= 1/100
            IGABA[int(t),i] = sum(-T2*kI2+sum(W*R)-Gsc*R+Iaff[t,i])
            if Iaff[int(t),i]<-0.1:
                Iaff[int(t),i] += 100
                IGABA[int(t),i] = sum(-T2*kI2+sum(W*R)-Gsc*R+Iaff[t,i])
        # R = np.log(I[t+dt,i]/(Vm_thr/VL))
        # R[np.isnan(R)] = 0
        VGABA =  Vm-I[int(t),i]/gGABAI
        t += dt
        if t >= 700:
            break
    i += 1

Vv = V
Iv = IGABA    
plt.figure(),   
plt.plot(Vv,Iv)
plt.plot(Vv[139],Iv[139],'ro')
plt.plot(Vv[456],Iv[456],'ro')
plt.plot(Vv[742],Iv[742],'ro')
#plt.axhline(y = 0,'b-')
plt.axhline(y = 0,xmin=-90,xmax=0,ls='--',color='black')
plt.axvline(x = Vv[139], ls='--',color='black')
plt.axvline(x = Vv[456], ls='--',color='black')
plt.axvline(x = Vv[742], ls='--',color='black')
plt.title('I v.s. V')
plt.xlabel('Membrane potential(mV)')
plt.ylabel('Current(uA/cm^2)')
