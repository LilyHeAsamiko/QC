# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 06:25:53 2020

@author: LilyHeAsamiko
"""

#%matplotlib inline
#The model is specified in terms of the following parameters:
#Each Duffing oscillator is specified by a frequency ν, anharmonicity α, and drive strength r, which result in the Hamiltonian terms:

#2πνa†a+παa†a(a†a−1)+2πr(a+a†)×D(t),
#where D(t) is the signal on the drive channel for the qubit, and a† and a are, respectively, the creation and annihilation operators for the qubit. Note that the drive strength r sets the scaling of the control term, with D(t) assumed to be a complex and unitless number satisfying |D(t)|≤1. - A coupling between a pair of oscillators (l,k) is specified by the coupling strength J, resulting in an exchange coupling term:

#2πJ(ala†k+a†lak),
#where the subscript denotes which qubit the operators act on. - Additionally, for numerical simulation, it is necessary to specify a cutoff dimension; the Duffing oscillator model is infinite dimensional, and computer simulation requires restriction of the operators to a finite dimensional subspace.


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
from qiskit import QuantumCircuit, execute

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

def MyPuzzle():
    """Small circuit solution for puzzle"""
    # Do circuit
    qc = QuantumCircuit(4)
    # Oracle
    qc.h(3)
    qc.ccx(0,1,2)
    qc.h(2)
    qc.x([1,3])
    qc.h(2)
    qc.mct([0,1,3],2)
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
    def swap_registers(circuit, n):
        for qubit in range(n//2):
            circuit.swap(qubit, n-qubit-1)
        return circuit
    def qft_rotations(circuit, n):
        """Performs qft on the first n qubits in circuit (without swaps)"""
        if n == 0:
            return circuit
        n -= 1
        circuit.h(n)
        for qubit in range(n):
            circuit.cu1(np.pi/2**(n-qubit), qubit, n)
        qft_rotations(circuit, n)
    
    qft_rotations(circuit, n)
    swap_registers(circuit, n)
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

#      ┌───┐┌─────────────────────────┐                                                                                 »
# q_0: ┤ H ├┤0                        ├─────────────────────────────────────────────────────────────────────────────────»
#      ├───┤│                         │┌─────────────────────────┐┌─────────────────────────┐                           »
# q_1: ┤ H ├┤                         ├┤0                        ├┤0                        ├───────────────────────────»
#      ├───┤│                         ││                         ││                         │┌─────────────────────────┐»
# q_2: ┤ H ├┤                         ├┤                         ├┤                         ├┤0                        ├»
#      ├───┤│                         ││                         ││                         ││                         │»
# q_3: ┤ H ├┤                         ├┤                         ├┤                         ├┤                         ├»
#      ├───┤│  ControlGroverIteration ││                         ││                         ││                         │»
# q_4: ┤ H ├┤1                        ├┤1 ControlGroverIteration ├┤1 ControlGroverIteration ├┤1                        ├»
#      ├───┤│                         ││                         ││                         ││  ControlGroverIteration │»
# q_5: ┤ H ├┤2                        ├┤2                        ├┤2                        ├┤2                        ├»
#      ├───┤│                         ││                         ││                         ││                         │»
# q_6: ┤ H ├┤3                        ├┤3                        ├┤3                        ├┤3                        ├»
#      ├───┤│                         ││                         ││                         ││                         │»
# q_7: ┤ H ├┤4                        ├┤4                        ├┤4                        ├┤4                        ├»
#      └───┘└─────────────────────────┘└─────────────────────────┘└─────────────────────────┘└─────────────────────────┘»
# c: 4/═════════════════════════════════════════════════════════════════════════════════════════════════════════════════»
#                                                                                                                       »
# «                                                                                                                 »
# «q_0: ────────────────────────────────────────────────────────────────────────────────────────────────────────────»
# «                                                                                                                 »
# «q_1: ────────────────────────────────────────────────────────────────────────────────────────────────────────────»
# «     ┌─────────────────────────┐┌─────────────────────────┐┌─────────────────────────┐                           »
# «q_2: ┤0                        ├┤0                        ├┤0                        ├───────────────────────────»
# «     │                         ││                         ││                         │┌─────────────────────────┐»
# «q_3: ┤                         ├┤                         ├┤                         ├┤0                        ├»
# «     │                         ││                         ││                         ││                         │»
# «q_4: ┤1                        ├┤1                        ├┤1                        ├┤1                        ├»
# «     │  ControlGroverIteration ││  ControlGroverIteration ││  ControlGroverIteration ││                         │»
# «q_5: ┤2                        ├┤2                        ├┤2                        ├┤2 ControlGroverIteration ├»
# «     │                         ││                         ││                         ││                         │»
# «q_6: ┤3                        ├┤3                        ├┤3                        ├┤3                        ├»
# «     │                         ││                         ││                         ││                         │»
# «q_7: ┤4                        ├┤4                        ├┤4                        ├┤4                        ├»
# «     └─────────────────────────┘└─────────────────────────┘└─────────────────────────┘└─────────────────────────┘»
# «c: 4/════════════════════════════════════════════════════════════════════════════════════════════════════════════»
# «                                                                                                                 »
# «                                                                                                                 »
# «q_0: ────────────────────────────────────────────────────────────────────────────────────────────────────────────»
# «                                                                                                                 »
# «q_1: ────────────────────────────────────────────────────────────────────────────────────────────────────────────»
# «                                                                                                                 »
# «q_2: ────────────────────────────────────────────────────────────────────────────────────────────────────────────»
# «     ┌─────────────────────────┐┌─────────────────────────┐┌─────────────────────────┐┌─────────────────────────┐»
# «q_3: ┤0                        ├┤0                        ├┤0                        ├┤0                        ├»
# «     │                         ││                         ││                         ││                         │»
# «q_4: ┤1                        ├┤1                        ├┤1                        ├┤1                        ├»
# «     │                         ││                         ││                         ││                         │»
# «q_5: ┤2 ControlGroverIteration ├┤2 ControlGroverIteration ├┤2 ControlGroverIteration ├┤2 ControlGroverIteration ├»
# «     │                         ││                         ││                         ││                         │»
# «q_6: ┤3                        ├┤3                        ├┤3                        ├┤3                        ├»
# «     │                         ││                         ││                         ││                         │»
# «q_7: ┤4                        ├┤4                        ├┤4                        ├┤4                        ├»
# «     └─────────────────────────┘└─────────────────────────┘└─────────────────────────┘└─────────────────────────┘»
# «c: 4/════════════════════════════════════════════════════════════════════════════════════════════════════════════»
# «                                                                                                                 »
# «                                                                                      ┌───────┐┌─┐         
# «q_0: ─────────────────────────────────────────────────────────────────────────────────┤0      ├┤M├─────────
# «                                                                                      │       │└╥┘┌─┐      
# «q_1: ─────────────────────────────────────────────────────────────────────────────────┤1      ├─╫─┤M├──────
# «                                                                                      │  QFT† │ ║ └╥┘┌─┐   
# «q_2: ─────────────────────────────────────────────────────────────────────────────────┤2      ├─╫──╫─┤M├───
# «     ┌─────────────────────────┐┌─────────────────────────┐┌─────────────────────────┐│       │ ║  ║ └╥┘┌─┐
# «q_3: ┤0                        ├┤0                        ├┤0                        ├┤3      ├─╫──╫──╫─┤M├
# «     │                         ││                         ││                         │└───────┘ ║  ║  ║ └╥┘
# «q_4: ┤1                        ├┤1                        ├┤1                        ├──────────╫──╫──╫──╫─
# «     │                         ││                         ││                         │          ║  ║  ║  ║ 
# «q_5: ┤2 ControlGroverIteration ├┤2 ControlGroverIteration ├┤2 ControlGroverIteration ├──────────╫──╫──╫──╫─
# «     │                         ││                         ││                         │          ║  ║  ║  ║ 
# «q_6: ┤3                        ├┤3                        ├┤3                        ├──────────╫──╫──╫──╫─
# «     │                         ││                         ││                         │          ║  ║  ║  ║ 
# «q_7: ┤4                        ├┤4                        ├┤4                        ├──────────╫──╫──╫──╫─
# «     └─────────────────────────┘└─────────────────────────┘└─────────────────────────┘          ║  ║  ║  ║ 
# «c: 4/═══════════════════════════════════════════════════════════════════════════════════════════╩══╩══╩══╩═
# «                                                                                                0  1  2  3 

## Execute and see results
emulator = Aer.get_backend('qasm_simulator')
job0 = execute(qc, emulator, shots=2048 )
hist0 = job0.result().get_counts()
plot_histogram(hist0)
output = get_noise(prob_1, prob_2, qc)

#Put into EEG
P0 = hist0['0000']/2048
P1 = hist0['1000']/2048
P2 = 1-P0-P1

P0_noisy = output[1]['0000']/2048
P1_noisy = output[1]['1000']/2048
P2_noisy = 1-P0_noisy-P1_noisy

#     ┌───┐┌───┐          
#q_0: ┤ H ├┤ X ├──■───────
#     ├───┤└─┬─┘  │       
#q_1: ┤ H ├──■────■────■──
#     ├───┤       │  ┌─┴─┐
#q_2: ┤ H ├───────■──┤ X ├
#     └───┘     ┌─┴─┐├───┤
#q_3: ──────────┤ X ├┤ Z ├
#               └───┘└───┘

# each signal as two labels marked by two different physicians.
# consider VQE on each signal contains two qubits and the 0 is 0, 1 is 1. The reuslt is the higher possibility one.(Each signal is independent.)

#Number of qubits
# nQ = 3
# #There are 3 qubits: Q0,Q1,Q2.
# #Number of seeds (random sequences)
# nseeds = 5
# #Number of Cliffords in the sequence (start, stop, steps)
# nCliffs = np.arange(1,200,20)
# #2Q RB on Q0,Q2 and 1Q RB on Q1
# rb_pattern = [[0,1],[2]]
# #Do three times as many 1Q Cliffords
# length_multiplier = [1,3]

#Generate the RB sequences
# rb_opts = {}
# rb_opts['length_vector'] = nCliffs
# rb_opts['nseeds'] = nseeds
# rb_opts['rb_pattern'] = rb_pattern
# rb_opts['length_multiplier'] = length_multiplier
# rb_circs, xdata = rb.randomized_benchmarking_seq(**rb_opts)

# #The H is an arbitary state (with a global phase)
# #The U is an identity (with a global phase)
# backend = qiskit.Aer.get_backend('qasm_simulator')
# basis_gates = ['u1','u2','u3','cx'] # use U,CX 
#print(np.around(job.result().get_unitary(), 3))    

