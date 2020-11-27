# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 07:24:45 2020

@author: LilyHeAsamiko
"""

#%matplotlib inline
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

def grover_iteration_5_16():
    """Small circuit with 5/16 solutions"""
    # Do circuit
    qc = QuantumCircuit(4)
    # Oracle
    qc.h([2,3])
    qc.ccx(0,1,2)
    qc.h(2)
    qc.x(2)
    qc.ccx(0,2,3)
    qc.x(2)
    qc.h(3)
    qc.x([1,3])
    qc.h(2)
    qc.mct([0,1,3],2)
    qc.x([1,3])
    qc.h(2)
    # Diffuser
    qc.h(range(3))
    qc.x(range(3))
    qc.z(3)
    qc.mct([0,1,2],3)
    qc.x(range(3))
    qc.h(range(3))
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




# Create QuantumCircuit
t = 4   # no. of counting qubits
n = 4   # no. of searching qubits
qc = QuantumCircuit(n+t, t) # Circuit with n+t qubits and t classical bits

# Initialise all qubits to |+>
for qubit in range(t+n):
    qc.h(qubit)
    
# Create controlled-Grover
gi = grover_iteration_5_16().to_gate()
cgi = gi.control()
cgi.label = "ControlGroverIteration"


# If Include registered 0-3 bits
qc.append(gi, range(t))

# Compute Energy
# Measure counting qubits
qc.measure(range(t), range(t))

# Display the circuit
qc.draw()
#
#                                                                                                                3 »
#     ┌───┐┌─────────────┐┌─┐         
#q_0: ┤ H ├┤0            ├┤M├─────────
#     ├───┤│             │└╥┘┌─┐      
#q_1: ┤ H ├┤1            ├─╫─┤M├──────
#     ├───┤│  circuit561 │ ║ └╥┘┌─┐   
#q_2: ┤ H ├┤2            ├─╫──╫─┤M├───
#     ├───┤│             │ ║  ║ └╥┘┌─┐
#q_3: ┤ H ├┤3            ├─╫──╫──╫─┤M├
#     ├───┤└─────────────┘ ║  ║  ║ └╥┘
#q_4: ┤ H ├────────────────╫──╫──╫──╫─
#     ├───┤                ║  ║  ║  ║ 
#q_5: ┤ H ├────────────────╫──╫──╫──╫─
#     ├───┤                ║  ║  ║  ║ 
#q_6: ┤ H ├────────────────╫──╫──╫──╫─
#     ├───┤                ║  ║  ║  ║ 
#q_7: ┤ H ├────────────────╫──╫──╫──╫─
#     └───┘                ║  ║  ║  ║ 
#c: 4/═════════════════════╩══╩══╩══╩═
#                          0  1  2  3 
# Execute and see results
emulator = Aer.get_backend('qasm_simulator')
job1 = execute(qc, emulator, shots=1024)
hist1 = job1.result().get_counts()
plot_histogram(hist1)


# Measure Control bits only
qc = QuantumCircuit(n+t, t) # Circuit with n+t qubits and t classical bits

# Initialise all qubits to |+>
for qubit in range(t+n):
    qc.h(qubit)
    
# Create controlled-Grover
gi = grover_iteration_5_16().to_gate()
cgi = gi.control()
cgi.label = "ControlGroverIteration"

# Begin controlled Grover iterations
iterations = 1
for qubit in range(t):
    for i in range(iterations):
        qc.append(cgi, [qubit] + [*range(t, n+t)])
    iterations *= 2

#«q_0: ───────────────────────────────────────────────────────────────────────────────────────────────────────────────
#«                                                                                                                    
#«q_1: ───────────────────────────────────────────────────────────────────────────────────────────────────────────────
#«                                                                                                                    
#«q_2: ───────────────────────────────────────────────────────────────────────────────────────────────────────────────
#«     ┌─────────────────────────┐┌─────────────────────────┐┌─────────────────────────┐┌─────────────────────────┐┌─┐
#«q_3: ┤0                        ├┤0                        ├┤0                        ├┤0                        ├┤M├
#«     │                         ││                         ││                         ││                         │└╥┘
#«q_4: ┤1                        ├┤1                        ├┤1                        ├┤1                        ├─╫─
#«     │                         ││                         ││                         ││                         │ ║ 
#«q_5: ┤2 ControlGroverIteration ├┤2 ControlGroverIteration ├┤2 ControlGroverIteration ├┤2 ControlGroverIteration ├─╫─
#«     │                         ││                         ││                         ││                         │ ║ 
#«q_6: ┤3                        ├┤3                        ├┤3                        ├┤3                        ├─╫─
#«     │                         ││                         ││                         ││                         │ ║ 
#«q_7: ┤4                        ├┤4                        ├┤4                        ├┤4                        ├─╫─
#«     └─────────────────────────┘└─────────────────────────┘└─────────────────────────┘└─────────────────────────┘ ║ 
#«c: 4/═════════════════════════════════════════════════════════════════════════════════════════════════════════════╩═
#«   

qft_dagger = qft(4).to_gate().inverse()
qft_dagger.label = "QFT†"
    
# Do inverse QFT on counting qubits
qc.append(qft_dagger, range(t))

# Measure counting qubits
qc.measure(range(t), range(t))

# Display the circuit
qc.draw()

#┌───────┐┌─┐         
#«q_0: ─────────────────────────────────────────────────────────────────────────────────┤0      ├┤M├─────────
#«                                                                                      │       │└╥┘┌─┐      
#«q_1: ─────────────────────────────────────────────────────────────────────────────────┤1      ├─╫─┤M├──────
#«                                                                                      │  QFT† │ ║ └╥┘┌─┐   
#«q_2: ─────────────────────────────────────────────────────────────────────────────────┤2      ├─╫──╫─┤M├───
#«     ┌─────────────────────────┐┌─────────────────────────┐┌─────────────────────────┐│       │ ║  ║ └╥┘┌─┐
#«q_3: ┤0                        ├┤0                        ├┤0                        ├┤3      ├─╫──╫──╫─┤M├
#«     │                         ││                         ││                         │└───────┘ ║  ║  ║ └╥┘
#«q_4: ┤1                        ├┤1                        ├┤1                        ├──────────╫──╫──╫──╫─
#«     │                         ││                         ││                         │          ║  ║  ║  ║ 
#«q_5: ┤2 ControlGroverIteration ├┤2 ControlGroverIteration ├┤2 ControlGroverIteration ├──────────╫──╫──╫──╫─
#«     │                         ││                         ││                         │          ║  ║  ║  ║ 
#«q_6: ┤3                        ├┤3                        ├┤3                        ├──────────╫──╫──╫──╫─
#«     │                         ││                         ││                         │          ║  ║  ║  ║ 
#«q_7: ┤4                        ├┤4                        ├┤4                        ├──────────╫──╫──╫──╫─
#«     └─────────────────────────┘└─────────────────────────┘└─────────────────────────┘          ║  ║  ║  ║ 
#«c: 4/═══════════════════════════════════════════════════════════════════════════════════════════╩══╩══╩══╩═
#«                                                                                                0  1  2  3 
## Execute and see results
emulator = Aer.get_backend('qasm_simulator')
job0 = execute(qc, emulator, shots=2048 )
hist0 = job0.result().get_counts()
plot_histogram(hist0)
#Put into EEG
P0 = hist0[0]
P1 = hist0[3]





#train and test


#readin data
url1 = 'http://homepages.cae.wisc.edu/~ece539/data/eeg/nic23a1.txt'
nic23a1 = urllib.request.urlopen(url1)
#s1=requests.get(url1).content
#nic23a1=pd.read_csv(io.StringIO(s1.decode('utf-8')))
#nic23a1 = pd.read_csv('D:/PhD in Oxford,Ethz,KI,others/OxfordPh/QC/ML(hybrid and momentum-space UCC)/nic23a1.txt')
url2 = 'http://homepages.cae.wisc.edu/~ece539/data/eeg/nic23a3.txt'
nic23a3 = urllib.request.urlopen(url2)
#another obervations
url21 = 'http://homepages.cae.wisc.edu/~ece539/data/eeg/nic8a1.txt'
nic8a1 = urllib.request.urlopen(url21)
url22 = 'http://homepages.cae.wisc.edu/~ece539/data/eeg/nic8a3.txt'
nic8a3 = urllib.request.urlopen(url22)

Dims = 29
Labels = 8
# Do not consider the effect of time first(Spacially only)
Width = int(Dims/Labels)+1 # particle number (batch size)
tx1 = []
tx2 = []
#print(np.shape(nic23a1.readlines()))
for line1 in nic23a1.readlines():
    tx1.append(line1.split())
for line2 in nic23a3.readlines():
    tx2.append(line2.split())
tx1 = np.array(tx1)
tx2 = np.array(tx2) #same data as tx1 but with different labels
rows,cols = np.shape(tx1) #cols = Dims + Labels
print(rows)
print(cols)
phi1 = np.ones((rows,Labels))
phi2 = np.ones((rows,Labels))
dataset = [] #896*(29+8)， consider about single person's data first(phi1 only in PSO, phi2 = 0)
label1 = []
label2 = []
for i in range(rows):
    dataset.append(tx1[i][range(Dims)])
    label1.append(tx1[i][range(Dims,Dims+Labels)])
    label2.append(tx2[i][range(Dims,Dims+Labels)])
label1 = np.array(label1,dtype = float)
label2 = np.array(label2,dtype = float)
Tl = 100
T = int(rows/Tl)+1

sim_result = []
correct_result = []
train = 0.75
#trainidx = sample(range(np.shape(tx1)[0]),k = int(0.75*np.shape(tx1)[0]))
#labeltrain = label1[trainidx][:]
#testidx = sample(range(np.shape(tx1)[0]),k = int(0.25*np.shape(tx1)[0]))
#labeltest = label1[testidx][:]
Tr_Acc = np.zeros((int(np.shape(dataset)[0]*train),1))
Tr_TP = np.zeros((int(np.shape(dataset)[0]*train),1))
Tr_TN = np.zeros((int(np.shape(dataset)[0]*train),1))
Tr_FP = np.zeros((int(np.shape(dataset)[0]*train),1))
Tr_FN = np.zeros((int(np.shape(dataset)[0]*train),1))
Tr_F = np.zeros((int(np.shape(dataset)[0]*train),1))

Te_Acc = np.zeros(((np.shape(dataset)[0]-np.shape(Tr_Acc)[0]),1))
Te_TP = np.zeros(((np.shape(dataset)[0]-np.shape(Tr_Acc)[0]),1))
Te_TN = np.zeros(((np.shape(dataset)[0]-np.shape(Tr_Acc)[0]),1))
Te_FP = np.zeros(((np.shape(dataset)[0]-np.shape(Tr_Acc)[0]),1))
Te_FN = np.zeros(((np.shape(dataset)[0]-np.shape(Tr_Acc)[0]),1))
Te_F = np.zeros(((np.shape(dataset)[0]-np.shape(Tr_Acc)[0]),1))

W = np.ones((np.shape(label1)[0],np.shape(label1)[1]))
g = 2*np.log(np.sqrt(2))
#test_result = []

# cutoff dimension
dim_oscillators = 3
# frequencies for transmon drift terms
# Number of oscillators in the model is determined from len(oscillator_freqs)
oscillator_freqs = [5.0e9, 5.2e9] #harmonic term 
anharm_freqs = [-0.33e9, -0.33e9] #anharmonic term

# drive strengths
drive_strengths = [0.02e9, 0.02e9]

# specify coupling as a dictionary (qubits 0 and 1 are coupled with a coefficient 0.002e9)
coupling_dict = {(0,1): 0.002e9}

#sample duration for pulse instructions in accordance
dt = backend_config.dt #1e-9
dt
backend_defaults = backend.defaults()

# unit conversion factors -> all backend properties returned in SI (Hz, sec, etc)
GHz = 1.0e9 # Gigahertz
MHz = 1.0e6 # Megahertz
us = 1.0e-6 # Microseconds
ns = 1.0e-9 # Nanoseconds

qubit = 0 # qubit we will analyze
default_qubit_freq = backend_defaults.qubit_freq_est[qubit] # Default qubit frequency in Hz. 
#print(f"Qubit {qubit} has an estimated frequency of {default_qubit_freq/ GHz} GHz.")

# scale data (specific to each device)
scale_factor = 1e-14

# number of shots for our experiments
NUM_SHOTS = 1024

### Collect the necessary channels
drive_chan = pulse.DriveChannel(qubit)
meas_chan = pulse.MeasureChannel(qubit)
acq_chan = pulse.AcquireChannel(qubit)

# Drive pulse parameters (us = microseconds)
drive_sigma_us = 0.075                     # This determines the actual width of the gaussian
drive_samples_us = drive_sigma_us*8        # This is a truncating parameter, because gaussians don't have 
                                           # a natural finite length
drive_sigma = get_closest_multiple_of_16(drive_sigma_us * us /dt)       # The width of the gaussian in units of dt
drive_samples = get_closest_multiple_of_16(drive_samples_us * us /dt) 
# The truncating parameter in units of dt
# Find out which measurement map index is needed for this qubit
meas_map_idx = None
for i, measure_group in enumerate(backend_config.meas_map):
    if qubit in measure_group:
        meas_map_idx = i
        break
assert meas_map_idx is not None, f"Couldn't find qubit {qubit} in the meas_map!"

# Get default measurement pulse from instruction schedule map
inst_sched_map = backend_defaults.instruction_schedule_map
measure = inst_sched_map.get('measure', qubits=backend_config.meas_map[meas_map_idx])
# We will sweep 40 MHz around the estimated frequency, with 75 frequencies
#num_freqs = 75
#ground_sweep_freqs = default_qubit_freq + np.linspace(-20*MHz, 20*MHz, num_freqs)
#ground_freq_sweep_program = create_ground_freq_sweep_program(ground_sweep_freqs, drive_power=0.3)

# create the Duffing oscillator system model(returns a PulseSystemModel object, which is a general object for storing model information required for simulation with the PulseSimulator.)
two_qubit_model = duffing_system_model(dim_oscillators=dim_oscillators,
                                       oscillator_freqs=oscillator_freqs,
                                       anharm_freqs=anharm_freqs,
                                       drive_strengths=drive_strengths,
                                       coupling_dict=coupling_dict,
                                       dt=dt)


# <span id="papermill-error-cell" style="color:red; font-family:Helvetica Neue, Helvetica, Arial, sans-serif; font-size:2em;">Execution using papermill encountered an exception here and stopped:</span>

#calibrate pi pulse on each qubit using Ihnis
#4.1 Constructing the schedules
# list of qubits to be used throughout the notebook

qubits = [0, 1]
# Construct a measurement schedule and add it to an InstructionScheduleMap
meas_amp = 0.025
meas_samples = 1200
meas_sigma = 4
meas_width = 1150
meas_pulse = GaussianSquare(duration=meas_samples, amp=meas_amp,
                            sigma=meas_sigma, width=meas_width)

acq_sched = pulse.Acquire(meas_samples, pulse.AcquireChannel(0), pulse.MemorySlot(0))
acq_sched += pulse.Acquire(meas_samples, pulse.AcquireChannel(1), pulse.MemorySlot(1))

measure_sched = pulse.Play(meas_pulse, pulse.MeasureChannel(0)) | pulse.Play(meas_pulse, pulse.MeasureChannel(1)) | acq_sched

inst_map = pulse.InstructionScheduleMap()
inst_map.add('measure', qubits, measure_sched)

#Rabi schedules
#recall: Rabii oscillation
# The magnetic moment is thus {\displaystyle {\boldsymbol {\mu }}={\frac {\hbar }{2}}\gamma {\boldsymbol {\sigma }}}{\boldsymbol {\mu }}={\frac {\hbar }{2}}\gamma {\boldsymbol {\sigma }}. 
# The Hamiltonian of this system is then given by {H} =-{{\mu }}\cdot{B} =-{\frac {\hbar }{2}}\omega _{0}\sigma _{z}-{\frac {\hbar }{2}}\omega _{1}(\sigma _{x}\cos \omega t-\sigma _{y}\sin \omega t)}\mathbf {H} =-{\boldsymbol {\mu }}\cdot \mathbf {B} =-{\frac {\hbar }{2}}\omega _{0}\sigma _{z}-{\frac {\hbar }{2}}\omega _{1}(\sigma _{x}\cos \omega t-\sigma _{y}\sin \omega t) where {\displaystyle \omega _{0}=\gamma B_{0}}\omega _{0}=\gamma B_{0} and {\displaystyle \omega _{1}=\gamma B_{1}}\omega _{1}=\gamma B_{1}
# Now, let the qubit be in state {\displaystyle |0\rangle }{\displaystyle |0\rangle } at time {\displaystyle t=0}t=0. Then, at time {\displaystyle t}t, the probability of it being found in state {\displaystyle |1\rangle }|1\rangle  is given by {\displaystyle P_{0\to 1}(t)=\left({\frac {\omega _{1}}{\Omega }}\right)^{2}\sin ^{2}\left({\frac {\Omega t}{2}}\right)}{\displaystyle P_{0\to 1}(t)=\left({\frac {\omega _{1}}{\Omega }}\right)^{2}\sin ^{2}\left({\frac {\Omega t}{2}}\right)} where {\displaystyle \Omega ={\sqrt {(\omega -\omega _{0})^{2}+\omega _{1}^{2}}}}\Omega ={\sqrt {(\omega -\omega _{0})^{2}+\omega _{1}^{2}}}
# the qubit oscillates between the {\displaystyle |0\rangle }|0\rang and {\displaystyle |1\rangle }|1\rangle  states. 
# The maximum amplitude for oscillation is achieved at {\displaystyle \omega =\omega _{0}}\omega =\omega _{0}, which is the condition for resonance. 
# At resonance, the transition probability is given by {\displaystyle P_{0\to 1}(t)=\sin ^{2}\left({\frac {\omega _{1}t}{2}}\right)}{\displaystyle P_{0\to 1}(t)=\sin ^{2}\left({\frac {\omega _{1}t}{2}}\right)}

# pi pulse: 
# To go from state {\displaystyle |0\rangle }|0\rang to state {\displaystyle |1\rangle }|1\rangle  it is sufficient to adjust the time {\displaystyle t}t during which the rotating field acts such that {\displaystyle {\frac {\omega _{1}t}{2}}={\frac {\pi }{2}}}{\frac {\omega _{1}t}{2}}={\frac {\pi }{2}} or {\displaystyle t={\frac {\pi }{\omega _{1}}}}t={\frac {\pi }{\omega _{1}}}
#  If a time intermediate between 0 and {\displaystyle {\frac {\pi }{\omega _{1}}}}{\frac {\pi }{\omega _{1}}} is chosen, we obtain a superposition of {\displaystyle |0\rangle }|0\rang and {\displaystyle |1\rangle }|1\rangle . 
# Approximation with pi/2 pulse: 
# In particular for {\displaystyle t={\frac {\pi }{2\omega _{1}}}}t={\frac {\pi }{2\omega _{1}}}, we have a {\displaystyle {\frac {\pi }{2}}}{\frac {\pi }{2}} pulse, which acts as: {\displaystyle |0\rangle \to {\frac {|0\rangle +i|1\rangle }{\sqrt {2}}}}{\displaystyle |0\rangle \to {\frac {|0\rangle +i|1\rangle }{\sqrt {2}}}}.
# The equations are essentially identical in the case of a two level atom in the field of a laser when the generally well satisfied rotating wave approximation is made. Then {\displaystyle \hbar \omega _{0}}\hbar \omega _{0} is the energy difference between the two atomic levels, {\displaystyle \omega }\omega  is the frequency of laser wave and Rabi frequency {\displaystyle \omega _{1}}\omega _{1} is proportional to the product of the transition electric dipole moment of atom {\displaystyle {\vec {d}}}{\vec {d}} and electric field {\displaystyle {\vec {E}}}{\vec {E}} of the laser wave that is {\displaystyle \omega _{1}\propto \hbar \ {\vec {d}}\cdot {\vec {E}}}{\displaystyle \omega _{1}\propto \hbar \ {\vec {d}}\cdot {\vec {E}}}.


#experiments

meas_amps = np.linspace(0, 0.9, 48)
meas_sigma = 16
meas_duration = 128

meas_channels = [pulse.DriveChannel(0), pulse.DriveChannel(1)]

rabi_experiments, rabi_amps = rabi_schedules(amp_list=meas_amps,
                                                 qubits=qubits,
                                                 pulse_width=meas_duration,
                                                 pulse_sigma=meas_sigma,
                                                 drives=meas_channels,
                                                 inst_map=inst_map,
                                                 meas_map=[[0, 1]])
rabi_experiments[10].draw()

#ground_freq_sweep_job = backend.run(ground_freq_sweep_program)
#print(ground_freq_sweep_job.job_id())
#job_monitor(ground_freq_sweep_job)
# Get the job data (average)
#ground_freq_sweep_data = get_job_data(ground_freq_sweep_job, average=True)
#To simulate the Rabi experiments, assemble the Schedule list into a qobj. When assembling, pass the PulseSimulator as the backend.#To simulate the Rabi experiments, assemble the Schedule list into a qobj. When assembling, pass the PulseSimulator as the backend.

# instantiate the pulse simulator
backend_sim = PulseSimulator()

# compute frequencies from the Hamiltonian
qubit_lo_freq = two_qubit_model.hamiltonian.get_qubit_lo_from_drift()

rabi_qobj = assemble(rabi_experiments,
                     backend=backend,
                     qubit_lo_freq=qubit_lo_freq,
                     meas_level=1,
                     meas_return='avg',
                     shots=256)

# run the simulation
rabi_result = backend_sim.run(rabi_qobj, two_qubit_model).result()

rabifit = RabiFitter(rabi_result, rabi_amps, qubits, fit_p0 = [0.5,0.5,0.6,1.5])

plt.figure(figsize=(15, 10))
q_offset = 0
multiplier = 0.5
for qubit in qubits:
    ax = plt.subplot(2, 2, qubit + 1)
    #Xvmin, Xvmax = ax.xaxis.get_data_interval()
    #Yvmin, Yvmax = ax.yaxis.get_data_interval()
    #print(Xvmin, Xvmax,Yvmin, Yvmax)
    Xvmin = multiplier * np.floor(0.1 / multiplier)
    Xvmax = multiplier * np.ceil(0.5 / multiplier)
    ax.set_xlim([Xvmin, Xvmax])
    rabifit.plot(qubit, ax=ax)
    print('Pi Amp: %f'%rabifit.pi_amplitude(qubit))
plt.show()

count = 0
p = 0.5
for l in range(np.shape(label1)[0]):
    pg = []
    F = 0
    for i in range(np.shape(label1)[1]):
        if label1[l][i] == 1 and label2[l][i] == 1: #label1[l][i] == 1 & label2[l][i] == 1
            sim_result.append(1)
            flag = 1
        elif label1[l][i] == 0 and label2[l][i] == 0:
            sim_result.append(0)
            flag = 0
        else:
            t = 4   # no. of counting qubits
            n = 4   # no. of searching qubits
            # Measure Control bits only
            qc = QuantumCircuit(n+t, t) # Circuit with n+t qubits and t classical bits
            
            # Initialise all qubits to |+>
            for qubit in range(t+n):
                qc.h(qubit)
                
            # Create controlled-Grover
            gi = grover_iteration_5_16().to_gate()
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
            #simulate
            emulator = Aer.get_backend('qasm_simulator')
            job0 = execute(qc, emulator, shots=256)
            hist0 = job0.result().get_counts()
#            plot_histogram(hist0)
            #Put into EEG
            P0 = hist0['0101']/256
            P1 = hist0['1011']/256            
            p = P0*0+P1*1
            flag = round(p)
            
            eta = 0.1
            tol = 0.1
            #run simulation
            # consider single person data first, no pgd, use p instead
            if i < np.shape(label1)[1]:
                xtemp = np.array(dataset[l],dtype = float)
                x = xtemp[int(i*Width):int(i*Width+Width)]
                Wtemp = np.repeat(W[l][i],Width)
            else:
                xtemp = np.array(dataset[l],dtype = float)
                x = xtemp[int(-Width-1):-1]
                Wtemp = np.repeat(W[l][i],Width)                       
            for steps in range(50):
                dw = x*(1-2*label1[l][i])
                dw *= eta/len(x)
                Wtemp = Wtemp - dw
                if np.mean(abs((np.ones((np.shape(x)))+0.1)/(1+0.1-np.exp(Wtemp*x)) - label1[l][i])) < tol:
                    W[l][i] = np.mean(Wtemp)                    
                    break
            if np.mean(abs(np.ones((np.shape(x)))+0.1)/(1+0.1-np.exp(Wtemp*x)))>0.5:
                pred = 1
            else:
                pred = 0  
            sim_result.append(pred) 
        print(np.size(sim_result))
        if np.mod(l,100)<75:
            if sim_result[-1] == 1 and flag == 1:
                Tr_TP[l//100*75+np.mod(l,100)] += 1 
                Tr_Acc[l//100*75+np.mod(l,100)] += 1
            elif sim_result[-1] == 0 and flag == 1:
                Tr_TN[l//100*75+np.mod(l,100)] += 1 
            elif sim_result[-1] == 1 and flag == 0:
                Tr_FP[l//100*75+np.mod(l,100)] += 1 
            elif sim_result[-1] == 0 and flag == 0:
                Tr_FN[l//100*75+np.mod(l,100)] += 1       
                Tr_Acc[l//100*75+np.mod(l,100)] += 1
        else:
            if sim_result[-1] == 1 and flag == 1:
                Te_TP[l//100*25+np.mod(l,100)-75] += 1 
                Te_Acc[l//100*25+np.mod(l,100)-75] += 1
            elif sim_result[-1] == 0 and flag == 1:
                Te_TN[l//100*25+np.mod(l,100)-75] += 1 
            elif sim_result[-1] == 1 and flag == 0:
                Te_FP[l//100*25+np.mod(l,100)-75] += 1 
            elif sim_result[-1] == 0 and flag == 0:
                Te_FN[l//100*25+np.mod(l,100)-75] += 1       
                Te_Acc[l//100*25+np.mod(l,100)-75] += 1
        pg.append(sim_result[-1])
    pg = np.array(pg)
    qubits = [0, 1]

# N =  int(Dims*dim_oscillators/2)

# meas_amp = Dims/128*0.025/2
# meas_samples = int(rows/16*1200/10)
# meas_sigma = int(Dims/16*4)
# meas_width = int(rows/128*1150/2)
    meas_sigma = 16
    meas_duration = 128
    
    meas_channels = [pulse.DriveChannel(0), pulse.DriveChannel(1)]
    
    rabi_experiments, rabi_amps = rabi_schedules(amp_list=pg,
                                                 qubits=qubits,
                                                 pulse_width=meas_duration,
                                                 pulse_sigma=meas_sigma,
                                                 drives=meas_channels,
                                                 inst_map=inst_map,
                                                 meas_map=[[0, 1]])
    
    meas_pulse = GaussianSquare(duration=meas_samples, amp=meas_amp,
                                sigma=meas_sigma, width=meas_width)

    acq_sched = pulse.Acquire(meas_samples, pulse.AcquireChannel(0), pulse.MemorySlot(0))
    acq_sched += pulse.Acquire(meas_samples, pulse.AcquireChannel(1), pulse.MemorySlot(1))
    
    measure_sched = pulse.Play(meas_pulse, pulse.MeasureChannel(0)) | pulse.Play(meas_pulse, pulse.MeasureChannel(1)) | acq_sched
    
    inst_map = pulse.InstructionScheduleMap()
    inst_map.add('measure', qubits, measure_sched)
 
    backend_sim = PulseSimulator()

    # compute frequencies from the Hamiltonian
    qubit_lo_freq = two_qubit_model.hamiltonian.get_qubit_lo_from_drift()
    
    rabi_qobj = assemble(rabi_experiments,
                         backend=backend,
                         qubit_lo_freq=qubit_lo_freq,
                         meas_level=1,
                         meas_return='avg',
                         shots=256)
    
    # run the simulation
    rabi_result = backend_sim.run(rabi_qobj, two_qubit_model).result()
    
    rabifit = RabiFitter(rabi_result, rabi_amps, qubits, fit_p0 = [1.5,1.5,1.6,2.5])

    plt.figure(figsize=(15, 10))
    q_offset = 0
    multiplier = 0.5
    for qubit in qubits:
#        ax = plt.subplot(2, 2, qubit + 1)
        #Xvmin, Xvmax = ax.xaxis.get_data_interval()
        #Yvmin, Yvmax = ax.yaxis.get_data_interval()
        #print(Xvmin, Xvmax,Yvmin, Yvmax)
#        Xvmin = multiplier * np.floor(0.1 / multiplier)
#        Xvmax = multiplier * np.ceil(0.5 / multiplier)
#        ax.set_xlim([Xvmin, Xvmax])
#        rabifit.plot(qubit, ax=ax)
        print('Pi Amp: %f'%rabifit.pi_amplitude(qubit))
#    plt.show()
    rabi_pred = (rabifit.ydata['0'][0]['mean']+rabifit.ydata['0'][1]['mean'])/2
    
    #Benchmark maximizing F-score(PSO)
    count += 1
    if np.mod(count,Tl) < Tl*train:
#            pg = sim_result[-9:-1] 
        precision = Tr_TP[l//100*75+np.mod(l,100)]/(Tr_TP[l//100*75+np.mod(l,100)]+Tr_FP[l//100*75+np.mod(l,100)])
        recall = Tr_TP[l//100*75+np.mod(l,100)]/(Tr_TP[l//100*75+np.mod(l,100)]+Tr_FN[l//100*75+np.mod(l,100)])
        if 2*precision*recall/(precision+recall) >= F:
            F = 2*precision*recall/(precision+recall)
            Tr_F[l//100*75+np.mod(l,100)] = F
        u = np.random.sample(1)
#           u = sample(list(np.linspace(0,1,100)),k=1)  
        L = (1/g)*abs(np.array(rabi_pred) - p) 
        if p > 0.5:
            rabi_pred = p - L*np.log(1/u) 
        else:
            rabi_pred = p + L*np.log(1/u)
    else: 
        precision = Te_TP[l//100*25+np.mod(l,100)-75]/(Te_TP[l//100*25+np.mod(l,100)-75]+Te_FP[l//100*25+np.mod(l,100)-75])
        recall = Te_TP[l//100*25+np.mod(l,100)-75]/(Te_TP[l//100*25+np.mod(l,100)-75]+Te_FN[l//100*25+np.mod(l,100)-75])
#            recall = Te_TP[l//100*25+np.mod(l,100)]/(Te_TP[l//100*25+np.mod(l,100)]+Te_FN[l//100*25+np.mod(l,100)])
        if 2*precision*recall/(precision+recall) >= F:               
            F = 2*precision*recall/(precision+recall)
            Te_F[l//100*25+np.mod(l,100)-75] = F
        if p > 0.5:
            rabi_pred = p - L*np.log(1/u) 
        else:
            rabi_pred = p + L*np.log(1/u)
        
    correct_result.append(np.array((pg-label1[l])*rabi_pred,dtype = int))

correct_result = np.array(correct_result)    
sim_result = np.reshape(np.array(sim_result),(np.shape(correct_result)))

Tr_Acc = Tr_Acc/8
Tr_Acc_mean = np.ones((np.shape(Tr_Acc)))
Tr_Acc_std = np.zeros((np.shape(Tr_Acc)))
for i in range(1,len(Tr_Acc)):
    Tr_Acc_mean[i] = np.mean(Tr_Acc[0:i])
    Tr_Acc_std[i] = np.std(Tr_Acc[0:i])
#show the training result
plt.figure()
plt.errorbar(range(len(Tr_Acc_mean)),Tr_Acc_mean[:], Tr_Acc_std[:])
plt.title('training accuracy')

plt.figure()
plt.boxplot(Tr_TP.T)
plt.title('traning TP')

plt.figure()
plt.boxplot(Tr_TN.T)
plt.title('traning TN')

plt.figure()
plt.boxplot(Tr_FP.T)
plt.title('traning FP')

plt.figure()
plt.boxplot(Tr_FN.T)
plt.title('traning FN')

#show the 200 epochs of test data
Te_Acc = Te_Acc[range(200)]/8
Te_Acc_mean = np.ones((np.shape(Te_Acc)))
Te_Acc_std = np.zeros((np.shape(Te_Acc)))
for i in range(1,len(Te_Acc)):
    Te_Acc_mean[i] = np.mean(Te_Acc[0:i])
    Te_Acc_std[i] = np.std(Te_Acc[0:i])

plt.figure()
plt.errorbar(range(len(Te_Acc_mean)),Te_Acc_mean[:], Te_Acc_std[:])
plt.title('testing accuracy')

plt.figure()
plt.boxplot(Te_TP[range(200)].T)
plt.title('testing TP')

plt.figure()
plt.boxplot(Te_TN[range(200)].T)
plt.title('testing TN')

plt.figure()
plt.boxplot(Te_FP[range(200)].T)
plt.title('testing FP')

plt.figure()
plt.boxplot(Te_FN[range(200)].T)
plt.title('testing FN')


#Time_related_PSO aggregated Results
Time_Agg_Acc = np.ones((np.shape(correct_result)[0],1))
Time_Agg_TP = np.zeros((np.shape(correct_result)[0],1))
Time_Agg_TN = np.zeros((np.shape(correct_result)[0],1))
Time_Agg_FP = np.zeros((np.shape(correct_result)[0],1))
Time_Agg_FN = np.zeros((np.shape(correct_result)[0],1))

for l in range(np.shape(correct_result)[0]):
    Time_Agg_TP[l] = sum((correct_result[l]==1)*(label1[l] == 1))
    Time_Agg_TN[l] = sum((correct_result[l]==1)*(label1[l] == 0))
    Time_Agg_FP[l] = sum((correct_result[l]==0)*(label1[l] == 1))
    Time_Agg_FN[l] = sum((correct_result[l]==0)*(label1[l] == 0))
    Time_Agg_Acc[l] = Time_Agg_TP[l]+Time_Agg_FN[l]

Time_Agg_Acc = Time_Agg_Acc[range(200)]/8
Time_Agg_Acc_mean = np.ones((np.shape(Time_Agg_Acc)))
Time_Agg_Acc_std = np.zeros((np.shape(Time_Agg_Acc)))
for i in range(1,len(Time_Agg_Acc)):
    Time_Agg_Acc_mean[i] = np.mean(Time_Agg_Acc[0:i])
    Time_Agg_Acc_std[i] = np.std(Time_Agg_Acc[0:i])

plt.figure()
plt.errorbar(range(len(Time_Agg_Acc_mean)),Time_Agg_Acc_mean[:], Time_Agg_Acc_std[:])
plt.title('Aggregated accuracy')

plt.figure()
plt.boxplot(Time_Agg_TP.T)
plt.title('Aggregated TP')

plt.figure()
plt.boxplot(Time_Agg_TN.T)
plt.title('Aggregated TN')

plt.figure()
plt.boxplot(Time_Agg_FP.T)
plt.title('Aggregated FP')

plt.figure()
plt.boxplot(Time_Agg_FN.T)
plt.title('Aggregated FN')
