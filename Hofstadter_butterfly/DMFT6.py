# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
#%matplotlib inline
#%matplotlib auto

import qiskit
#qiskit.__version__
from qiskit import IBMQ
MY_API_TOKEN = 'a93830f80226030329fc4e2e4d78c06bdf1942ce349fcf8f5c8021cfe8bd5abb01e4205fbd7b9c34f0b26bd335de7f1bcb9a9187a2238388106d16c6672abea2'
provider = IBMQ.enable_account(MY_API_TOKEN)

from qiskit.compiler import transpile, assemble
from qiskit.tools.jupyter import *
from qiskit.visualization import *

import matplotlib.pyplot as plt
import numpy as np
import math
import os
import sys
import io
import requests
import urllib
#store pi amplitudes
#given the drive and target indices, and the option to either start with the drive qubit in the ground or excited state, returns a list of experiments for observing the oscillations.
from IPython import display
import time
import pandas as pd

# importing Qiskit
import qiskit
from qiskit import IBMQ, Aer
from qiskit import QuantumRegister, QuantumCircuit, execute

# import basic plot tools
from qiskit.visualization import plot_histogram

from random import *
from qiskit.visualization.bloch import Bloch
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis

import qiskit.pulse as pulse
import qiskit.pulse.pulse_lib as pulse_lib
from qiskit.pulse.pulse_lib import Gaussian, GaussianSquare
from qiskit.compiler import assemble
from qiskit.ignis.characterization.calibrations import rabi_schedules, RabiFitter
#from qiskit.pulse.commands import SamplePulse
from qiskit.pulse import *
from qiskit.tools.monitor import job_monitor
# function for constructing duffing models
from qiskit.providers.aer.pulse import duffing_system_model

#We will experimentally find a π-pulse for each qubit using the following procedure: 
#- A fixed pulse shape is set - in this case it will be a Gaussian pulse. 
#- A sequence of experiments is run, each consisting of a Gaussian pulse on the qubit, followed by a measurement, with each experiment in the sequence having a subsequently larger amplitude for the Gaussian pulse. 
#- The measurement data is fit, and the pulse amplitude that completely flips the qubit is found (i.e. the π-pulse amplitude).
import warnings
warnings.filterwarnings('ignore')
from qiskit.tools.jupyter import *
get_ipython().run_line_magic('matplotlib', 'inline')

provider = IBMQ.get_provider(hub='ibm-q', group='open', project='main')
backend = provider.get_backend('ibmq_armonk')
backend_config = backend.configuration()
assert backend_config.open_pulse, "Backend doesn't support Pulse"
#from qiskit.providers.aer import QasmSimulator
from qiskit.providers.aer import PulseSimulator

from qiskit import QuantumCircuit, execute, Aer
from qiskit.visualization import plot_histogram
import qiskit.providers.aer.noise as noise

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
from itertools import combinations


def CZ(phi):
    CZ = [[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,np.exp(1j*phi)]]
    return CZ

def alpha(j, c1u_D, c1d_D, GS):#XX YY
    a = j*(c1u_D+c1d_D)*GS[j-1][j-1]               
    return a

def Slm(theta, alpha, l, m, sigma1x, sigma1y, sigma1z):#XX YY
    A = np.array(np.zeros((4,4)))
#    assert(alpha in ['x','y','z']) 
    if alpha == 'x':
        S = sigma1x[0][1]
    elif alpha == 'y':
        S = sigma1y[0][1]
    elif alpha == 'z':
        S = sigma1z[1][1]
    for i in range(4):
        #if i in [l,m]:
        if i==l or i == m: 
            A[i][i] = np.exp(1j*S*theta/2)
        else:
            A[i][i] = 1
    if alpha == 'y':
        A = A.transpose()            
    return A

def Ra(theta, alpha, l):
    if alpha == 'x':
        if l == 1:
            Ra = np.exp(-1j*theta/2*sigma1x)  
        elif l == 2:
            Ra = np.exp(-1j*theta/2*sigma2x)
        elif l == 3:
            Ra = np.exp(-1j*theta/2*sigma3x)
        elif l == 4:
            Ra = np.exp(-1j*theta/2*sigma4x)
    if alpha == 'y':
        if l == 1:
            Ra = np.exp(-1j*theta/2*sigma1y)  
        elif l == 2:
            Ra = np.exp(-1j*theta/2*sigma2y)
        elif l == 3:
            Ra = np.exp(-1j*theta/2*sigma3y)
        elif l == 4:
            Ra = np.exp(-1j*theta/2*sigma4y)
    if alpha == 'z':
        if l == 1:
            Ra = np.exp(-1j*theta/2*sigma1z)  
        elif l == 2:
            Ra = np.exp(-1j*theta/2*sigma2z)
        elif l == 3:
            Ra = np.exp(-1j*theta/2*sigma3z)
        elif l == 4:
            Ra = np.exp(-1j*theta/2*sigma4z)              
    return Ra        

def omega(G0, G):
#    E_GS = np.fft(G)**2
    E_GS = (np.fft.fft(G*np.array([[1,0],[1,0]]))[1][1])**2
    E = 1/G0 -1/G
    wj =  E - E_GS
    return wj
   
def U3(theta, phi, Lambda):
    if theta == 0 and phi ==0:
        U.name = 'U1'
        U.theta = 0
        U.phi = 0
        U.Lambda = 0
        U.matrix = [[1,0],[0,np.exp(1j*Lambda)]]
    if theta == np.pi/2:
        U.name = 'U2'
        U.theta = np.pi/2
        U.phi = phi
        U.Lambda = Lambda
        U.matrix = [[1,-np.exp(1j*Lambda)],[np.exp(1j*phi),np.exp(1j*(phi+Lambda))]]/np.sqrt(2)
    else:
        U.name = 'U3'
        U.theta = np.pi/2
        U.phi = phi
        U.Lambda = Lambda
        U.matrix = [[np.cos(theta/2), -np.exp(1j*Lambda)*np.sin(theta/2)],[np.exp(1j*phi)*np.sin(theta/2), np.exp(1j*Lambda+1j*phi)*np.cos(theta/2)]]
    return U

# def GS(gates, n):                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     (gates,n):
#     """Creates a 4-qubit QFT circuit with n qubit operators"""
#     for ni in range(2^n):
#         assert(gates[ni].name in ['U1', 'U2', 'U3', 'Rx', 'Ry', 'Rz'])
#         theta = gates[ni].theta
#         phi = gates[ni].phi
#         Lambda = gates[ni].Lambda        
#         circuit[ni].u3(theta, phi, Lambda)
#         theta = gates[ni+1].theta
#         phi = gates[ni+1].phi
#         Lambda = gates[ni+1].Lambda 
#         circuit[ni+1].u3(theta, phi, Lambda)
#         circuit.cx([ni],[ni+1])  
#         theta = gates[ni+2].theta
#         phi = gates[ni+2].phi
#         Lambda = gates[ni+2].Lambda        
#         circuit[ni+2].u3(theta, phi, Lambda)
#         theta = gates[ni+3].theta
#         phi = gates[ni+3].phi
#         Lambda = gates[ni+3].Lambda 
#         circuit[ni+3].u3(theta, phi, Lambda)        
#     circuit = sub_circuit.to_instruction()
#     return circuit

def MyPuzzle(gates, n, GS):
    """Small circuit solution for puzzle"""
    # Do circuit.o
    qr = QuantumRegister(2^n+1)
    qc = QuantumCircuit(qr)
    gs = GS(gates,n)
    qc.append(GS, [qr[1],qr[2],qr[3],qr[4]])
    qc.draw()
    # Oracle
    qc.h(1)
    qc.x([1,3])
    qc.h(3)
    # Diffuser
    qc.h(range(3))
    qc.z(3)
    qc.cx(1,0)
    qc.mct([1,0,2],3)
    qc.cx(1,2)
#    qc.h(range(3))
    qc.z(3)
    return qc

def qft(n):
    """Creates an n-qubit QFT circuit"""
    circuit = QuantumCircuit(4)
    return circuit

def calculate_M(measured_int, t, n):
    """For Processing Output of Quantum Counting"""
    # Calculate Theta
    theta = (measured_int/(2**t))*math.pi*2
    print("Theta = %.5f" % theta)
    # Calculate No. of Solutions
    N = 2**n
    M = N * (math.sin(theta/2)**2)
    print("No. of Solutions = %.1f" % (N-M))
    # Calculate Upper Error Bound
    m = t - 1 #Will be less than this (out of scope) 
    err = (math.sqrt(2*M*N) + N/(2**(m-1)))*(2**(-m))
    print("Error < %.2f" % err)

def get_job_data(job, average):
    """Retrieve data from a job that has already run.
    Args:
        job (Job): The job whose data you want.
        average (bool): If True, gets the data assuming data is an average.
                        If False, gets the data assuming it is for single shots.
    Return:
        list: List containing job result data. 
    """
    job_results = job.result(timeout=120) # timeout parameter set to 120 s
    result_data = []
    for i in range(len(job_results.results)):
        if average: # get avg data
            result_data.append(job_results.get_memory(i)[qubit]*scale_factor) 
        else: # get single data
            result_data.append(job_results.get_memory(i)[:, qubit]*scale_factor)  
    return result_data

def get_closest_multiple_of_16(num):
    """Compute the nearest multiple of 16. Needed because pulse enabled devices require 
    durations which are multiples of 16 samples.
    """
    return (int(num) - (int(num)%16))

def exp2target(y):
    return np.ones((np.size(y)))/(1+np.exp(y))

def KL(p,q):#p:used distribution for data, q: designed distribution for model
    return -sum(p*np.log(q/p))

def get_noise(prob_1,prob_2,qc):
    
    # Error probabilities
#    prob_1 = 0.001  # 1-qubit gate
#    prob_2 = 0.01   # 2-qubit gate
    
    # Depolarizing quantum errors
    error_1 = noise.depolarizing_error(prob_1, 1)
    error_2 = noise.depolarizing_error(prob_2, 2)
    
    # Add errors to noise model
    noise_model = noise.NoiseModel()
    noise_model.add_all_qubit_quantum_error(error_1, ['u1', 'u2', 'u3'])
    noise_model.add_all_qubit_quantum_error(error_2, ['cx'])
    
    # Get basis gates from noise model
    basis_gates = noise_model.basis_gates
    
    # Make a circuit
    # Perform a noise simulation
    new_result = execute(qc, Aer.get_backend('qasm_simulator'),
                         basis_gates=basis_gates,
                         noise_model=noise_model).result()
    if np.empty_like(new_result)==0:
        new_counts = new_result.get_counts(0)
        plot_histogram(new_counts)        
        return [noise_model, new_counts]
    else:
        return noise_model

# with depth D = 2 and with periodic boundary conditions,there are only two pos-
# sible causal cones: a 4-qubit cone enclosing Nb = 3 blocks
# time step tau = 0.1

def get_raw_results(code,noise_model):#=None
    circuits = code.get_circuit_list()
#    raw_results = {}
#    for log in range(2):
#        job = execute( circuits[log], Aer.get_backend('qasm_simulator'), noise_model=noise_model)
#        raw_results[str(log)] = job.result().get_counts(str(log))
    table_results = {}
    for log in range(2):
        job = execute( circuits[log], Aer.get_backend('qasm_simulator'), noise_model=noise_model, shots=10000 )
        table_results[str(log)] = job.result().get_counts(str(log))

    P = lookuptable_decoding(raw_results,table_results)
    print('P =',P)
    return [P,table_results]#raw_results

import numpy as np    
#GS = GS(Gates, n)
# assume 4 qubits intialized as 0
sigma1x = np.array([[0, 1],[1, 0]], dtype = int) #|0><1| + |1><0|
sigma2x = sigma1x
sigma3x = sigma1x
sigma4x = sigma1x

sigma1y = np.array([[0, -1],[1,0]], dtype = int)*[[0, 1j],[1j, 0]] #-i|0><1| + i|1><0|
sigma2y = sigma1y
sigma3y = sigma1y
sigma4y = sigma1y

sigma1z = np.array([[1, 0],[0, -1]], dtype = int)  #|0><0| - |1><1|
sigma2z = sigma1z
sigma3z = sigma1z
sigma4z = sigma1z

I = np.array([[1, 0],[0, 1]], dtype = int) #|0><1| + |1><0|

# star z = infinite
# quantum:
# ground state: gamma1 and gamm2(half-filled Hubbard model) Ja = Jb = const
# steady state： Ja = Jb = 0 


c1d_D = (sigma1x-1j*sigma1y)/2 #= sigma1n
c2d_D = sigma1z*(sigma2x-1j*sigma2y)/2 #= sigma1z*sigma2n 
c1u_D = sigma1z*sigma2z*(sigma3x-1j*sigma3y)/2#= sigma1z*sigma2z*sigma3n
c2u_D = sigma1z*sigma2z*sigma3z*(sigma4x-1j*sigma4y)/2#= sigma1z*sigma2z*sigma3n
n1d = (I-sigma1z)/2
n2d = (I-sigma2z)/2
n1u = (I-sigma3z)/2
n2u = (I-sigma4z)/2
c1d = np.conj((sigma1x-1j*sigma1y)/2) #= sigma1n
c2d = np.conj(sigma1z*(sigma2x-1j*sigma2y)/2) #= sigma1z*sigma2n 
c1u = np.conj(sigma1z*sigma2z*(sigma3x-1j*sigma3y)/2)#= sigma1z*sigma2z*sigma3n
c2u = np.conj(sigma1z*sigma2z*sigma3z*(sigma4x-1j*sigma4y)/2)#= sigma1z*sigma2z*sigma3n
#generally c1sgt = np.exp(1j*tau*HSAIM)*c1sigma*np.exp(-1j*tau*HSAIM)

Sp = c1u_D*c1d + c2u_D*c2d
Sn = n1u+n2u-(n1d+n2d)
Sz = Sn
etap = c1u_D*c1d_D - c2u_D*c2d_D
etan = n1u+n2u+n1d+n2d-2
etaz = etan


#Fist fix U and mu and give a initial guess for epsc and V
#XY gate 
#epsc = 0 #half filled
epsc = 10**5
N = 24 #total 1D chain' Total time steps
ts = 10 # [10, 100, 0.04, 100, 1.1]
tau = 6/ts #step_width
w = 2*np.pi/N
U =5*ts #XYgate #8*ts 
V = ts
mu = U/2 #half filled
dw = V**2/(w-epsc)
rho0 = np.sqrt(4*ts**2-epsc**2)/(2*np.pi*ts**2)  

Lambda = 4*10^(-5) # [4*10^(-5), 4*10^(-6), 4*10^(-3), 4*10^(-6), 4*10^(-4)]    


Sw_XYO_t = []
Aw_XYO_t = []
A_XYO_t = []
B_XYO_t = []
HSIAM1_XYO_t = []
HSIAM2_XYO_t = []
HSIAM_XYO_t = []
Hlas_t = []

HSIAM1_XYO_t_R = []
HSIAM2_XYO_t_R = []                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
HSIAM_XYO_t_R = []
U_XYO_t_R = []
Hlas_t_R = []

c1sgtu = c1u
c1sgtd = c1d
Epsc = []
#first 
Sigx12 = Slm(np.pi/2, 'x', 1, 2, sigma1x, sigma1y, sigma1z)
Sigy12 = Slm(np.pi/2, 'y', 1, 2, sigma1x, sigma1y, sigma1z)
Sigz13 = Slm(np.pi/2, 'z', 1, 3, sigma1x, sigma1y, sigma1z)

GS = -1j*V/2*(Sigx12*Sigy12)/2/dw #4*4
#a1 = alpha(1, c1u_D, c1d_D, GS)[1][0] 
#a2 = alpha(2, c1u_D, c1d_D, GS)[1][0]
#initialize
G0_i = 1/(w + mu -dw) 
Sw_latt = 1
Gw = (4*ts**2-2*epsc**2)/(2*np.pi*ts**2*np.sqrt(4*ts**2-ts**2))*(w + mu - epsc -Sw_latt)-epsc*rho0/(w+mu-Sw_latt-1)
G1_i = 1*Gw
G2_i = 2*Gw
w1 = omega(G0_i, G1_i)
w2 = omega(G0_i, G2_i)
Sw = w +mu - dw - 1/(a1*(1/(w-w1)+1/(w+w1))+a2*(1/(w-w2)+1/(w+w2)))
Aw = rho0*(w + mu - Sw)


# update
G1_i = 1/(w1 + mu -dw)
G2_i = 1/(w2 + mu -dw)
Sw = w +mu - dw - 1/(a1*(1/(w-w1)+1/(w+w1))+a2*(1/(w-w2)+1/(w+w2)))
Aw = rho0*(w + mu - Sw)
w1 = omega(G0_i, G1_i)
w2 = omega(G0_i, G2_i)
iGimp_tau = 2*(a1*np.cos(w1*tau)+a2*np.cos(w2*tau))*tau/N

#HSIAM = U/4*(sigma1z*sigma3z-sigma1z-sigma3z)+mu/2*(sigma1z+sigma3z)-epsc/2*(sigma2z+sigma4z)+V/2*(sigma1x*sigma2x+sigma1y*sigma2y+sigma3x, *sigma4x+sigma3y*sigma4y)
#consider the elctrons with charge e moving on a lattice
# in an external magnetic filed 
e = 2.718281828
A = a1*a2
alpha = 1/2/np.pi
B = 2*np.pi*alpha/A/e
M = 1 #assume mass energy to be 1
Hlas_temp = 1

for t in np.linspace(0.0001,tau,N+1):
    V = t
    GS = -1j*V/2*(Sigx12*Sigy12)/2/dw #4*4
#    a1 = alpha(1, c1u_D, c1d_D, GS)[1][0] 
#    a2 = alpha(2, c1u_D, c1d_D, GS)[1][0]
#    U1 = Ra(-np.pi/2, 'x', 1)*Ra(-np.pi/2, 'y', 2)*Ra(-np.pi/2, 'x', 3)*Ra(-np.pi/2, 'y', 4)*np.exp(-1j*V/2*Sigz13)[0][1]*np.exp(-1j*V/2*Sigz13)[0][1]*Ra(-np.pi/2, 'x', 1)*Ra(-np.pi/2, 'x', 2)*Ra(-np.pi/2, 'x', 3)*Ra(-np.pi/2, 'x', 4)
    a1 = (c1u_D+c1d_D)*GS[0][0]   
    a2 = 2*(c1u_D+c1d_D)*GS[1][1] 
    U1 = np.exp(-1j*V/2*(sigma1x*sigma2x+sigma1y*sigma2y)*t)*np.exp(-1j*V/2*(sigma3x*sigma4x+sigma3y*sigma4y)*t)
    U2 = np.exp(-1j*U/4*Sigz13*t)[1][1]*np.exp(-1j*(U/4-mu/2)*sigma1z*t)[1][1]*np.exp(-1j*(U/4-mu/2)*sigma3z*t)[1][1]*np.exp(1j*epsc/2*sigma2z*t)[1][1]*np.exp(1j*epsc/2*sigma4z*t)[1][1]
    U_XYO_t = np.array([U1[0][0], U2])
#    U_XYO_t_R.append(np.array([np.real(U1[1][1]), np.real(U2)]))
    Hsiamu = U_XYO_t*n1u-mu*sum(0.5*(I-sigma1z)+0.5*(I-sigma3z))+epsc*sum(c2u_D*c2u)+V*sum(c1u_D*c2u)
    Hsiamd = U_XYO_t*n1d-mu*sum(0.5*(I-sigma1z)+0.5*(I-sigma3z))+epsc*sum(c2d_D*c2d)+V*sum(c1d_D*c2d)
    n = n1d+n1u
    Gimp_t = []
    nimp = 0
    for ts in np.linspace(0.0001,6*t,N+1):
#    ts = 6/t
        rho0 = np.sqrt(abs(4*ts**2-epsc**2))/(2*np.pi*ts**2)  
        nimp += 2*w*np.sqrt(abs(4*ts**2-epsc**2))/(2*np.pi*(6*t)**2)  
        c1sgtu = np.exp(1j*ts*Hsiamu)*c1sgtu*np.exp(-1j*ts*Hsiamu)
        c1sgtd = np.exp(1j*ts*Hsiamd)*c1sgtd*np.exp(-1j*ts*Hsiamd)
        Gimp_t.append(-1j*(c1sgtu*c1u_D)-1j*(c1u_D*c1sgtu)-1j*(c1sgtd*c1d_D)-1j*(c1d_D*c1sgtd))
    Gimp_wi = np.fft.fft(Gimp_t)
    S = Gimp_wi[1][1]
    epsc_temp = abs(nimp-n[1][1])
    if epsc_temp >10**5:
        epsc_temp = 10**5-1
    if epsc_temp < epsc:
        epsc = epsc_temp
        dw = V**2/(w-epsc)
        S = w + mu - dw- 1/(w+mu-dw-S)
        print(S)
        Epsc.append(epsc)
#    if sum(U_XYO_t>0)<3:
#        Gimp_wi = a1*(1/(w+1j*etan-w1)+1/(w+1j*etan-w1))+a2*(1/(w+1j*etan-w2)+1/(w+1j*etan-w2))
#    else:
#        Gimp_wi = a1*(1/(w+1j*etap-w1)+1/(w+1j*etap-w1))+a2*(1/(w+1j*etan-w2)+1/(w+1j*etap-w2))
#    print(Gimp_wi)
#    print(w +mu - dw - 1/Gimp_wi[1][1])
#    print(w + mu - Sw_XY_t[-1])
#    print(dw + 1/Gimp_wi[1][1])
#    print(-np.log(U_t)/tau/1j)
        Sw_XYO_t.append(S[0])
         #   Aw_XY_t.append(rho0*(w + mu - Sw_XY_t[-1]))
        Aw_XYO_t.append(rho0*(dw + 1/S[0]))
        HSIAM1_XYO_t.append(-np.log(U1)/tau/1j)
        HSIAM2_XYO_t.append(-np.log(U2)/tau/1j)
        HSIAM_XYO_t.append(-np.log(U_XYO_t)/tau/1j) 
        HSIAM1_XYO_t_R.append(np.real(-np.log(U1)/tau/1j))
        HSIAM2_XYO_t_R.append(np.real(-np.log(U2)/tau/1j))
        HSIAM_XYO_t_R.append(np.real(-np.log(U_XYO_t)/tau/1j))  
    #consider the elctrons with charge e moving on a lattice
    # in an external magnetic filed inhomogeeous electric field
        x = np.linspace(0, 10, 11)
        y = x
        m = 1
        nn = 2
        lambdax = 2*x/nn
        lambday = 4*y/m
        Lambda = lambdax
        kx = 2*np.pi/lambdax
        ky = 2*np.pi/lambday
        Ex  = kx**2//M
        Ey  = ky**2//M
        Ex[np.isnan(Ex)] = 0.0001
        E = Ex
        dE = E[2]-E[1]
        delta = mu*dE*Lambda/4
        Er = Ex
        vx = np.sqrt(4*Er*V)
        q = kx+ky
        if min(abs(delta -vx[0])) <max(0.5*vx):
            A_XYO_t.append(a1*a2)
            Btemp = 2*np.pi*alpha/A_XYO_t[-1]/e
            Btemp[np.isnan(Btemp)] = 0
            B_XYO_t.append(Btemp)
            delta2 = ky-kx
            gamma = np.sqrt((q**2-delta2**2)*(4*kx**2-q**2-4*kx*delta2+delta2**2))
            gamma[np.isnan(gamma)] = 0.00001
            cospy = gamma/2/q/(kx-delta2)
            cospx = gamma/2/q/kx
            q = kx*np.sin(np.arccos(cospx)) + kx*np.sin(np.arccos(cospy))  
            q[[np.isnan(q)]] = 0.00001
            qc = np.sqrt(4*kx*(kx-delta2)+delta2)
            qc[np.isnan(qc)] = 0.00001
            print('q')
            if min(q) < max(qc):
                print('qc')
                w = omega(G0_i, G1_i)
                J = w*gamma**2/2
                Hlas1 = J*np.exp(2*np.pi*1j*alpha*m)
                Hlas1[np.isnan(Hlas1)] = 0.00001
                Hlas_temp = np.mean(Hlas1)*(c2d_D+c2u_D)*(c1d+c1u)-np.mean(delta)*(c1d_D+c1u_D)*(c1d+c1u)-np.mean(delta)*(c2d_D+c2u_D)*(c2d+c2u)
                Hlas_temp[np.isnan(Hlas_temp)] = 0.0001
                Hlas_t.append(Hlas_temp)
                Hlas_t_R.append(np.real(Hlas_t)) 
                
HSIAM1_XYO_t = np.reshape(HSIAM1_XYO_t, 4, 1)
Hlas_t = np.reshape(Hlas_t, 4, 1)
HSIAM1_XYO_t_M = np.zeros((np.shape(HSIAM1_XYO_t)))
HSIAM2_XYO_t_M = np.zeros((np.shape(HSIAM1_XYO_t)))
HSIAM_XYO_t_M = np.zeros((np.shape(HSIAM1_XYO_t)))
Hlas_t_M = np.zeros((np.shape(Hlas_t)))
HSIAM1_XYO_t_R_M = np.zeros((np.shape(HSIAM1_XYO_t)))
HSIAM2_XYO_t_R_M = np.zeros((np.shape(HSIAM1_XYO_t)))
HSIAM_XYO_t_R_M = np.zeros((np.shape(HSIAM1_XYO_t)))
Hlas_t_R_M = np.zeros((np.shape(Hlas_t)))
#U_CZO_t_R_M = np.zeros((np.shape(U_CZO_t_R)))
HSIAM1_XYO_t_std = np.ones((np.shape(HSIAM1_XYO_t)))
HSIAM2_XYO_t_std = np.ones((np.shape(HSIAM1_XYO_t)))
HSIAM_XYO_t_std = np.ones((np.shape(HSIAM1_XYO_t)))
Hlas_t_std = np.ones((np.shape(Hlas_t)))
HSIAM1_XYO_t_R_std = np.ones((np.shape(HSIAM1_XYO_t)))
HSIAM2_XYO_t_R_std = np.ones((np.shape(HSIAM1_XYO_t)))
HSIAM_XYO_t_R_std = np.ones((np.shape(HSIAM1_XYO_t)))
Hlas_t_R_std = np.ones((np.shape(Hlas_t)))
#U_CZO_t_R_std = np.ones((np.shape(U_CZO_t_R)))

#initiate
HSIAM1_XYO_t_M[0] = HSIAM1_XYO_t[0]
HSIAM2_XYO_t_M[0] = HSIAM2_XYO_t[0]
HSIAM1_XYO_t_R_M[0] = np.real(HSIAM1_XYO_t_M[0])
HSIAM2_XYO_t_R_M[0] = np.real(HSIAM2_XYO_t_M[0])
Hlas_t_M[0] = Hlas_t[0]
Hlas_t_R_M[0] = np.real(Hlas_t_M[0])

HSIAM1_XYO_t_std[0] = HSIAM1_XYO_t[0]
HSIAM2_XYO_t_std[0] = HSIAM2_XYO_t[0]
HSIAM_XYO_t_std[0] = HSIAM_XYO_t[0][0]
Hlas_t_std[0] = Hlas_t[0]
HSIAM1_XYO_t_R_std[0] = np.real(HSIAM1_XYO_t_std[0])
HSIAM2_XYO_t_R_std[0] = np.real(HSIAM2_XYO_t_std[0])
HSIAM_XYO_t_R_std[0] = np.real(HSIAM_XYO_t_std[0])
Hlas_t_R_std[0] = np.real(Hlas_t_std[0])

for tst in range(1,4):
    HSIAM2_XYO_t_M[tst] = abs(HSIAM2_XYO_t_M[0] +np.random.normal(0,1))
    HSIAM2_XYO_t_R_M[tst] = np.real(HSIAM2_XYO_t_M[tst])

HSIAM_XYO_t_M[0] = HSIAM_XYO_t[0][0]
HSIAM_XYO_t_M[1] = abs(0.5*(HSIAM_XYO_t[0][0]+HSIAM_XYO_t[0][1]))
HSIAM_XYO_t_R_M[0] = np.real(HSIAM_XYO_t_M[0])
HSIAM_XYO_t_R_M[1] = np.real(HSIAM_XYO_t_M[0])

for tst in range(2,4):
    HSIAM_XYO_t_M[tst] = abs(HSIAM_XYO_t_M[0] +np.random.normal(0,1))
    HSIAM_XYO_t_R_M[tst] = np.real(HSIAM_XYO_t_M[tst])

for tst in range(1,4):
    HSIAM1_XYO_t_M[tst] = np.mean(HSIAM1_XYO_t_M[0:tst])
    HSIAM1_XYO_t_R_M[tst] = np.real(HSIAM1_XYO_t_R_M[tst])
    HSIAM2_XYO_t_M[tst] = abs(HSIAM2_XYO_t[0]+np.random.normal(0,1))
    HSIAM2_XYO_t_R_M[tst] = np.real(HSIAM2_XYO_t_M[tst])
    Hlas_t_M[tst] = np.mean(Hlas_t[0:tst])
    Hlas_t_R_M[tst] = np.real(Hlas_t_M[tst])

for tst in range(1,4):
    HSIAM1_XYO_t_std[tst] = np.std(HSIAM1_XYO_t_M[0:tst])
    HSIAM2_XYO_t_std[tst] = np.std(HSIAM2_XYO_t_M[0:tst])
    HSIAM_XYO_t_std[tst] = np.std(HSIAM_XYO_t_M[0:tst])
    Hlas_t_std[tst] = np.std(Hlas_t[0:tst])
    HSIAM1_XYO_t_R_std[tst] = np.real(HSIAM1_XYO_t_std[tst])
    HSIAM2_XYO_t_R_std[tst] = np.real(HSIAM2_XYO_t_std[tst])
    HSIAM_XYO_t_R_std[tst] = np.real(HSIAM_XYO_t_std[tst])
    Hlas_t_R_std[tst] = np.real(Hlas_t[tst])

HSIAM1_XYO_t_std[HSIAM1_XYO_t_std==0] = 10**(-5)
HSIAM2_XYO_t_std[HSIAM2_XYO_t_std==0] = 10**(-5)
HSIAM_XYO_t_std[HSIAM_XYO_t_std==0] = 10**(-5)
HSIAM1_XYO_t_R_std[HSIAM1_XYO_t_R_std==0] = 10**(-5)
HSIAM2_XYO_t_R_std[HSIAM2_XYO_t_R_std==0] = 10**(-5)
HSIAM_XYO_t_R_std[HSIAM_XYO_t_R_std==0] = 10**(-5)
#U_CZO_t_R_std[U_CZO_t_R_std==0] = 10**(-5)
Hlas_t_std[Hlas_t_std==0] = 10**(-5)
Hlas_t_R_std[Hlas_t_R_std==0] = 10**(-5)
    
HSIAM1_XYO_t_CV = np.array(HSIAM1_XYO_t_M, dtype = float)/np.array(HSIAM1_XYO_t_R_std, dtype = float)
HSIAM2_XYO_t_CV = np.array(HSIAM2_XYO_t_M, dtype = float)/np.array(HSIAM2_XYO_t_R_std, dtype = float)
HSIAM_XYO_t_CV = np.array(HSIAM_XYO_t_M, dtype = float)/np.array(HSIAM_XYO_t_R_std, dtype = float)
Hlas_t_CV = np.array(Hlas_t_M, dtype = float)/np.array(Hlas_t_std, dtype = float)
HSIAM1_XYO_R_CV = np.array(HSIAM1_XYO_t_M, dtype = float)/np.array(HSIAM1_XYO_t_R_std, dtype = float)
HSIAM2_XYO_R_CV = np.array(HSIAM2_XYO_t_M, dtype = float)/np.array(HSIAM2_XYO_t_R_std, dtype = float)
HSIAM_XYO_R_CV = np.array(HSIAM_XYO_t_M, dtype = float)/np.array(HSIAM_XYO_t_R_std, dtype = float)
#U_CZO_t_R_CV = np.array(U_CZO_t_R_M, dtype = float)/np.array(U_CZO_t_R_std, dtype = float)
Hlas_t_R_CV = np.array(np.real(Hlas_t_R_M), dtype = float)/np.array(np.real(Hlas_t_R_std), dtype = float)

#HSIAM_tau_CV = np.array(HSIAM_tau_M, dtype = float)/np.array(HSIAM_tau_std, dtype = float)
plt.figure()
plt.plot(HSIAM1_XYO_t_CV)
plt.title('HSIAM1_XYO_odd_t_CV')
plt.figure()
plt.plot(HSIAM2_XYO_t_CV)
plt.title('HSIAM2_XYO_odd_t_CV')
plt.figure()
plt.plot(HSIAM_XYO_t_CV)
plt.title('HSIAM_XYO_odd_t_CV')


plt.figure()
plt.plot(HSIAM1_XYO_t_M)
plt.title('HSIAM1_XYO_odd_t_M')
plt.figure()
plt.plot(HSIAM2_XYO_t_M)
plt.title('HSIAM2_XYO_odd_t_M')
plt.figure()
plt.plot(HSIAM_XYO_t_M)
plt.title('HSIAM_XYO_odd_t_M')

plt.figure()
plt.plot(HSIAM1_XYO_R_CV)
plt.title('HSIAM1_XYO_odd_R_CV')
plt.figure()
plt.plot(HSIAM2_XYO_R_CV)
plt.title('HSIAM2_XYO_odd_R_CV')
plt.figure()
plt.plot(HSIAM_XYO_R_CV)
plt.title('HSIAM_XYO_odd_R_CV')

plt.figure()
plt.plot(HSIAM1_XYO_t_R_M)
plt.title('HSIAM1_XYO_odd_R_M')
plt.figure()
plt.plot(HSIAM2_XYO_t_R_M)
plt.title('HSIAM2_XYO_odd_R_M')
plt.figure()
plt.plot(HSIAM_XYO_t_R_M)
plt.title('HSIAM_XYO_odd_R_M')

plt.figure()
plt.plot(Hlas_t_M)
plt.title('Hlas_t_M')
plt.figure()
plt.plot(Hlas_t_R_M)
plt.title('Hlas_t_R_M')
plt.figure()
plt.plot(Hlas_t_CV)
plt.title('Hlas_t_CV')
plt.figure()
plt.plot(Hlas_t_R_CV)
plt.title('Hlas_t_R_CV')

A_XYO_t = np.reshape(A_XYO_t,4,1)
B_XYO_t = np.reshape(B_XYO_t,4,1)
plt.figure()
plt.plot(A_XYO_t)
plt.title('A_XYO_t')
plt.figure()
plt.plot(B_XYO_t)
plt.title('B_XYO_t')

plt.figure()
plt.errorbar(x = range(4), y = HSIAM1_XYO_R_CV, yerr = HSIAM1_XYO_t_R_std/10)
plt.title('HSIAM1_XYO_odd_R_CV_err')
plt.figure()
plt.errorbar(x = range(4), y = HSIAM2_XYO_R_CV, yerr = HSIAM2_XYO_t_R_std*10**10)
plt.title('HSIAM2_XYO_odd_R_CV_err')
plt.figure()
plt.errorbar(x = range(4), y = Hlas_t_R_CV, yerr = Hlas_t_R_std*10**10)
plt.title('Hlas_t_R_CV_err')
plt.figure()
plt.errorbar(x = range(4), y = HSIAM_XYO_R_CV, yerr= HSIAM_XYO_t_R_std*10**10)
plt.title('HSIAM_XYO_odd_R_CV_err')

plt.figure()
plt.contour(np.reshape(Hlas_t_R,[1,4])) 
plt.ylabel('U_XYO_odd_t_R')
plt.xlabel('tau')
plt.title('Countour of U_XYO_odd_t_R versus tau')
plt.figure()
plt.contour(np.reshape(HSIAM_XYO_t_R,[1,2])) 
plt.ylabel('HSIAM_XY_odd_t_R')
plt.xlabel('tau')
plt.title('Countour of HSIAM_XY_odd_t_R versus tau')

plt.figure()
plt.polar(Sw_XYO_t)
plt.ylabel('Sw_XYO_t')
plt.xlabel('tau')
plt.title('PolarPlot of S_XY_odd_t versus tau')

plt.figure()
plt.polar(Aw_XYO_t)
plt.ylabel('Aw_O_t')
plt.xlabel('tau')
plt.title('PolarPlot of a XY_odd_t versus tau')

plt.figure()
plt.polar(A_XYO_t)
plt.ylabel('A_XYO_t')
plt.xlabel('tau')
plt.title('PolarPlot of A_XY_odd_t versus tau')

plt.figure()
plt.polar(B_XYO_t)
plt.ylabel('B_O_t')
plt.xlabel('tau')
plt.title('PolarPlot of B XY_odd_t versus tau')

plt.figure()
plt.plot(Hlas_t)
plt.ylabel('Hlas_t')
plt.xlabel('tau')
plt.title('Hlas versus tau')

