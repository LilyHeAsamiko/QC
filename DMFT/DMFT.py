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

def GS(gates, n):                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     (gates,n):
    """Creates a 4-qubit QFT circuit with n qubit operators"""
    for ni in range(2^n):
        assert(gates[ni].name in ['U1', 'U2', 'U3', 'Rx', 'Ry', 'Rz'])
        theta = gates[ni].theta
        phi = gates[ni].phi
        Lambda = gates[ni].Lambda        
        circuit[ni].u3(theta, phi, Lambda)
        theta = gates[ni+1].theta
        phi = gates[ni+1].phi
        Lambda = gates[ni+1].Lambda 
        circuit[ni+1].u3(theta, phi, Lambda)
        circuit.cx([ni],[ni+1])  
        theta = gates[ni+2].theta
        phi = gates[ni+2].phi
        Lambda = gates[ni+2].Lambda        
        circuit[ni+2].u3(theta, phi, Lambda)
        theta = gates[ni+3].theta
        phi = gates[ni+3].phi
        Lambda = gates[ni+3].Lambda 
        circuit[ni+3].u3(theta, phi, Lambda)        
    circuit = sub_circuit.to_instruction()
    return circuit

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

#draw puzzle
# Create QuantumCircuit
t = 4   # no. of counting qubits
n = 4   # no. of searching qubits
qc = QuantumCircuit(n+t, t) # Circuit with n+t qubits and t classical bits


# Initialise all qubits to |+>
for qubit in range(t+n):
    qc.h(qubit)
    
# Create controlled-Grover
gi = MyPuzzle().to_gate()
cgi = gi.control()
cgi.label = "ControlGroverIteration"

# If Include registered 0-3 bits
qc.append(gi, range(t))

# Compute Energy
# Measure counting qubits
qc.measure(range(t), range(t))

# Display the circuit
qc.draw()

# Error probabilities
prob_1 = 0.001  # 1-qubit gate
prob_2 = 0.01   # 2-qubit gate

#     ┌───┐┌────────────┐┌─┐         
#q_0: ┤ H ├┤0           ├┤M├─────────
#     ├───┤│            │└╥┘┌─┐      
#q_1: ┤ H ├┤1           ├─╫─┤M├──────
#     ├───┤│  circuit11 │ ║ └╥┘┌─┐   
#q_2: ┤ H ├┤2           ├─╫──╫─┤M├───
#     ├───┤│            │ ║  ║ └╥┘┌─┐
#q_3: ┤ H ├┤3           ├─╫──╫──╫─┤M├
#     ├───┤└────────────┘ ║  ║  ║ └╥┘
#q_4: ┤ H ├───────────────╫──╫──╫──╫─
#     ├───┤               ║  ║  ║  ║ 
#q_5: ┤ H ├───────────────╫──╫──╫──╫─
#     ├───┤               ║  ║  ║  ║ 
#q_6: ┤ H ├───────────────╫──╫──╫──╫─
#     ├───┤               ║  ║  ║  ║ 
#q_7: ┤ H ├───────────────╫──╫──╫──╫─
#     └───┘               ║  ║  ║  ║ 
#c: 4/════════════════════╩══╩══╩══╩═
#                         0  1  2  3 

# Execute and see results
emulator = Aer.get_backend('qasm_simulator')
job1 = execute(qc, emulator, shots=1024)
hist1 = job1.result().get_counts()
plot_histogram(hist1)

output = get_noise(prob_1, prob_2, qc)

# Measure Control bits only
qc = QuantumCircuit(n+t, t) # Circuit with n+t qubits and t classical bits

# Initialise all qubits to |+>
for qubit in range(t+n):
    qc.h(qubit)
    
# Create controlled-Grover
gi = MyPuzzle().to_gate()
cgi = gi.control()
cgi.label = "ControlGroverIteration"

# Begin controlled Grover iterations
iterations = 1
for qubit in range(t):
    for i in range(iterations):
        qc.append(cgi, [qubit] + [*range(t, n+t)])
    iterations *= 2

# Display the circuit
qc.draw()

#     ┌───┐┌────────────┐┌─┐┌───┐               ┌─────────────────────────┐                           »
#q_0: ┤ H ├┤0           ├┤M├┤ H ├───────────────┤0                        ├───────────────────────────»
#     ├───┤│            │└╥┘└┬─┬┘┌───┐          │                         │┌─────────────────────────┐»
#q_1: ┤ H ├┤1           ├─╫──┤M├─┤ H ├──────────┤                         ├┤0                        ├»
#     ├───┤│  circuit11 │ ║  └╥┘ └┬─┬┘┌───┐     │                         ││                         │»
#q_2: ┤ H ├┤2           ├─╫───╫───┤M├─┤ H ├─────┤                         ├┤                         ├»
#     ├───┤│            │ ║   ║   └╥┘ └┬─┬┘┌───┐│                         ││                         │»
#q_3: ┤ H ├┤3           ├─╫───╫────╫───┤M├─┤ H ├┤                         ├┤                         ├»
#     ├───┤└───┬───┬────┘ ║   ║    ║   └╥┘ └───┘│  ControlGroverIteration ││                         │»
#q_4: ┤ H ├────┤ H ├──────╫───╫────╫────╫───────┤1                        ├┤1 ControlGroverIteration ├»
#     ├───┤    ├───┤      ║   ║    ║    ║       │                         ││                         │»
#q_5: ┤ H ├────┤ H ├──────╫───╫────╫────╫───────┤2                        ├┤2                        ├»
#     ├───┤    ├───┤      ║   ║    ║    ║       │                         ││                         │»
#q_6: ┤ H ├────┤ H ├──────╫───╫────╫────╫───────┤3                        ├┤3                        ├»
#     ├───┤    ├───┤      ║   ║    ║    ║       │                         ││                         │»
#q_7: ┤ H ├────┤ H ├──────╫───╫────╫────╫───────┤4                        ├┤4                        ├»
#     └───┘    └───┘      ║   ║    ║    ║       └─────────────────────────┘└─────────────────────────┘»
#c: 4/════════════════════╩═══╩════╩════╩═════════════════════════════════════════════════════════════»
#                         0   1    2    3                                                             »
#«                                                                                                                 »
# Execute and see results after qft
emulator = Aer.get_backend('qasm_simulator')
job2 = execute(qc, emulator, shots=1024)
hist2 = job2.result().get_counts()
plot_histogram(hist2)
output = get_noise(prob_1, prob_2, qc)

# Measure Control bits only ffor complete process
qc = QuantumCircuit(n+t, t) # Circuit with n+t qubits and t classical bits

# Initialise all qubits to |+>
for qubit in range(t+n):
    qc.h(qubit)
    
# Create controlled-Grover
gi = MyPuzzle().to_gate()
cgi = gi.control()
cgi.label = "ControlGroverIteration"

# Begin controlled Grover iterations
iterations = 1
for qubit in range(t):
    for i in range(iterations):
        qc.append(cgi, [qubit] + [*range(t, n+t)])
    iterations *= 2

qft_dagger = qft(4).to_gate().inverse()
qft_dagger.label = "QFT†"
    
# Do inverse QFT on counting qubits
qc.append(qft_dagger, range(t))

# Measure counting qubits
qc.measure(range(t), range(t))

# Display the circuit
qc.draw()

n = 2
# According to Rotoselect, apply the gates all as RY gates
class gates:
    def __init__(self):
        self.name = 'RUXY'
        self.theta = 0
        self.phi = 0
        self.Lambda = 4 #default
        self.matrix = [[np.cos(0), -np.sin(0)],[np.sin(0), np.cos(0)]]
        return self

Gates = []
for ni in range(2^n):
    temp = gates()
    temp.name = 'RY'
    temp.theta = 0
    temp.phi = 0
    temp.Lambda = 4 
    temp.matrix = [[np.cos(temp.theta/2), -np.sin(temp.theta/2)],[np.sin(temp.theta/2), np.cos(temp.theta/2)]]       
    Gates.append(temp)
    
GS = GS(Gates, n)
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

N = 100 #total 1D chain'
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


Sp = c1u_D*c1d + c2u_D*c2d
Sn = n1u+n2u-(n1d+n2d)
Sz = Sn
etap = c1u_D*c1d_D - c2u_D*c2d_D
etan = n1u+n2u+n1d+n2d-2
etaz = etan

t = 10 # [10, 100, 0.04, 100, 1.1]
Lambda = 4*10^(-5) # [4*10^(-5), 4*10^(-6), 4*10^(-3), 4*10^(-6), 4*10^(-4)]    
t = 100
Lambda = 4*10^(-6)
ts = t
#in this test, we assume U initialized as 1 
U = 1

dLambda = Lambda
if U > 0:
    U == 5#S>0 spin
    O1p = Sp/2
    O2p = Sp/2
    O1n = Sn/2
    O2n = Sn/2    
    N1 = N/2
    N2 = N/2
    dN = N/2
else:
    U == -5
    O1p = etap/2
    O2p = etap/2
    O1n = etan/2
    O2n = etan/2  
    Lambda_hat = []
    for N1 in np.linspace(N/2, 3*N/2, N):
        N2 = N - N1
        Lambda_temp = (N1-N2)/(N1+N2)
        if Lambda_temp - Lambda <= dLambda:
            Lambda_hat.append(Lambda_temp)
            dLambda = Lambda_temp - Lambda 
    if N1 == N2:
        dN = N1
    else:
        dN = abs(N1-N2) #here only O1p, O1n, O2p, O2n 

C = (O1p*O2n+O1n*O2p)/dN #N is the number of pairs of sites satidfying |a-b|=d
V =ts
tau = 6/ts
Z = 1
HSIAM = U/4*(sigma1z*sigma3z)+V*(sigma1x*sigma2x+sigma1y*sigma2y+sigma3x*sigma4x+sigma3y*sigma4y)/2    
hbar = 1 #h/2m
Uhat = np.exp(-1j*HSIAM*tau/hbar)
Uhat_D = np.exp(HSIAM*tau/hbar)
iGimpR = np.real(sigma1x*Uhat_D*sigma1x*Uhat)

#second
#classically
w1 = 2*np.pi/N1
w2 = 2*np.pi/N2
w = 2*np.pi/N
mu = U/2
ec = 0
iGimpR_tau = []
Z_w = []
Vnew_tau =[]
HSIAM_tau = []
Ts = []
Tau = [] #

for alpha in np.linspace(0,1,N):
    iGimpR_hat = alpha*np.cos(w1*tau)+(1-alpha)*np.cos(w2*tau)
    iGimpR_tau.append(iGimpR_hat)
    if ts is not t:
        V = Vnew
    if 1/(V**4*(alpha/w1**4+(1-alpha)/w2**4)) is not Z:
        if V**4*(alpha/w1**4+(1-alpha)/w2**4) < V*10:
            Z = 1/(V**4*(alpha/w1**4+(1-alpha)/w2**4))
        else:
            Z = Z*10 + np.random.normal(0,1) 
        Z_w.append(Z)
        ts = t*np.sqrt(Z)
        if abs(Vnew-V) < V*10: 
            Vnew = np.sqrt(Z)*ts
        else:
            Vnew = V*10 + np.random.normal(0,1) 
        if Vnew >  V*10**5:
            Vnew = t*10**5
        Vnew_tau.append(Vnew)
        dw = Vnew**2/(w-ec) #ignore first
        HSIAM = Uhat*(sigma1z*sigma3z-sigma1z-sigma3z)/4+mu*(sigma1z+sigma3z)/2-ec*(sigma2z+sigma4z)/2+Vnew*(sigma1x*sigma2x+sigma1y*sigma2y+sigma3x*sigma4x+sigma3y*sigma4y)/2
        HSIAM_tau.append(HSIAM) 
        Alpha1 = c1d_D 
    ts += 1/N
    tau = 6/ts
    Tau.append(tau)
    
rhoiGimpR, piGimpR = stats.spearmanr(np.array(iGimpR_tau).reshape(100,1), np.array(iGimpR_tau_M).reshape(100,1))
rhoZ, pZ = stats.spearmanr(np.array(Z_w).reshape(100,1), np.array(Z_w_M).reshape(100,1))
rhoVnew, pVnew = stats.spearmanr(np.array(Vnew_tau).reshape(100,1), np.array(Vnew_tau_M).reshape(100,1))
    

#result visualisation
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
from itertools import combinations
iGimpR_tau_M = np.zeros((np.shape(iGimpR_tau)))
Z_w_M = np.zeros((np.shape(Z_w)))
Vnew_tau_M = np.zeros((np.shape(Vnew_tau)))
HSIAM_tau_M = np.zeros((np.shape(HSIAM_tau)))
iGimpR_tau_std = np.ones((np.shape(iGimpR_tau)))
Z_w_std = np.ones((np.shape(Z_w)))
Vnew_tau_std = np.ones((np.shape(Vnew_tau)))
HSIAM_tau_std = np.ones((np.shape(HSIAM_tau)))

for tst in range(1,len(Tau)):
    iGimpR_tau_M[tst] = np.mean(iGimpR_tau[0:tst])
    Z_w_M[tst] = np.mean(Z_w[0:tst])
    Vnew_tau_M[tst] = np.mean(Vnew_tau[0:tst])
#    HSIAM_tau_M[tst] = np.mean(HSIAM_tau[0:tst])
    iGimpR_tau_std[tst] = np.std(iGimpR_tau[0:tst])
    Z_w_std[tst] = np.std(Z_w[0:tst])
    Vnew_tau_std[tst] = np.std(Vnew_tau[0:tst])
#    HSIAM_tau_std[tst] = np.std(HSIAM_tau[0:tst])
iGimpR_tau_CV = np.array(iGimpR_tau_M+, dtype = float)/np.array(iGimpR_tau_std, dtype = float)
Z_w_CV = np.array(Z_w_M, dtype = float)/np.array(Z_w_std, dtype = float)
Vnew_tau_CV = np.array(Vnew_tau_M, dtype = float)/np.array(Vnew_tau_std, dtype = float)
#HSIAM_tau_CV = np.array(HSIAM_tau_M, dtype = float)/np.array(HSIAM_tau_std, dtype = float)
plt.figure()
plt.plot(iGimpR_tau_CV)
plt.title('iGimpR_tau_CV')
plt.figure()
plt.plot(Z_w_CV)
plt.title('Z_w_CV')
plt.figure()
plt.plot(Vnew_tau_CV)
plt.title('Vnew_tau_CV')

plt.figure()
plt.plot(iGimpR_tau_M)
plt.title('iGimpR_tau_M')
plt.figure()
plt.plot(Z_w_M)
plt.title('Z_w_M')
plt.figure()
plt.plot(Vnew_tau_M)
plt.title('Vnew_tau_M')

plt.figure()
plt.errorbar(x = Tau, y = iGimpR_tau_M, yerr = iGimpR_tau_std)
plt.title('iGimpR_tau_err')
plt.figure()
plt.errorbar(x = Tau, y = Z_w_M, yerr = Z_w_std)
plt.title('Z_w_err')
plt.figure()
plt.errorbar(x = Tau, y = Vnew_tau_M, yerr= Vnew_tau_std)
plt.title('Vnew_tau_err')
