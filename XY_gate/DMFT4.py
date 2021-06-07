# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
#%matplotlib inline
#%matplotlib auto

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
# def alpha(j, sigma):
#     alpha = np.abs(j*c1)

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
qc.append(qft_dagger, crange(t))

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


Sp = c1u_D*c1d + c2u_D*c2d
Sn = n1u+n2u-(n1d+n2d)
Sz = Sn
etap = c1u_D*c1d_D - c2u_D*c2d_D
etan = n1u+n2u+n1d+n2d-2
etaz = etan
epsc = 0
N = 24 #total 1D chain' Total time steps
ts = 10 # [10, 100, 0.04, 100, 1.1]
tau = 6/ts #step_width
w = 2*np.pi/N
U =4*ts
V = ts
mu = U/2
dw = V**2/(w-epsc)
rho0 = np.sqrt(4*ts**2-epsc**2)/(2*np.pi*ts**2)  

Lambda = 4*10^(-5) # [4*10^(-5), 4*10^(-6), 4*10^(-3), 4*10^(-6), 4*10^(-4)]    

#first 
Sigx12 = Slm(np.pi/2, 'x', 1, 2, sigma1x, sigma1y, sigma1z)
Sigy12 = Slm(np.pi/2, 'y', 1, 2, sigma1x, sigma1y, sigma1z)
Sigz13 = Slm(np.pi/2, 'z', 1, 3, sigma1x, sigma1y, sigma1z)

GS = -1j*V/2*(Sigx12*Sigy12)/2/dw #4*4
a1 = alpha(1, c1u_D, c1d_D, GS)[1][0] 
a2 = alpha(2, c1u_D, c1d_D, GS)[1][0]
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

HSIAM = U/4*(sigma1z*sigma3z-sigma1z-sigma3z)+mu/2*(sigma1z+sigma3z)-epsc/2*(sigma2z+sigma4z)+V/2*(sigma1x*sigma2x+sigma1y*sigma2y+sigma3x, *sigma4x+sigma3y*sigma4y)

#CZ gate (ZZ even)
Sw_ZZ_t = []
Aw_ZZ_t = []
HSIAM1_ZZ_t = []
HSIAM2_ZZ_t = []
HSIAM3_ZZ_t = []
HSIAM_ZZ_t = []
HSIAM1_ZZ_t_R = []
HSIAM2_ZZ_t_R = []
HSIAM3_ZZ_t_R = []
HSIAM_ZZ_t_R = []
U_t_R = []
for t in np.linspace(0,tau,N+1):
    V = t
    GS = -1j*V/2*(Sigx12*Sigy12)/2/dw #4*4
    a1 = alpha(1, c1u_D, c1d_D, GS)[1][0] 
    a2 = alpha(2, c1u_D, c1d_D, GS)[1][0]
    U1 = Ra(np.pi/2, 'y', 1)*Ra(np.pi/2, 'y', 2)*Ra(np.pi/2, 'y', 3)*Ra(np.pi/2, 'y', 4)*np.exp(-1j*V/2*Sigx12)[0][0]*np.exp(-1j*V/2*Sigy12)[0][0]*Ra(-np.pi/2, 'y', 1)*Ra(-np.pi/2, 'y', 2)*Ra(-np.pi/2, 'y', 3)*Ra(-np.pi/2, 'y', 4)
    U2 = np.exp(-1j*U/4*Sigz13*t)[1][1]*np.exp(-1j*(U/4-mu/2)*sigma1z*t)[1][1]*np.exp(-1j*(U/4-mu/2)*sigma3z*t)[1][1]*np.exp(1j*ec/2*sigma2z*t)[1][1]*np.exp(1j*ec/2*sigma4z*t)[1][1]
    U3 = Ra(-np.pi/2, 'x', 1)*Ra(-np.pi/2, 'x', 2)*Ra(-np.pi/2, 'x', 3)*Ra(-np.pi/2, 'x', 4)*np.exp(-1j*V/2*Sigx12)[0][1]*np.exp(-1j*V/2*Sigy12)[0][1]*Ra(-np.pi/2, 'x', 1)*Ra(-np.pi/2, 'x', 2)*Ra(-np.pi/2, 'x', 3)*Ra(-np.pi/2, 'x', 4)
    U_t = np.array([U1[1][1], U2, U3[0][0]])
    U_t_R.append(np.array([np.real(U1[1][1]), np.real(U2), np.real(U3[0][0])]))
    if sum(U_t>0)<3:
        Gimp_wi = a1*(1/(w+1j*etan-w1)+1/(w+1j*etan-w1))+a2*(1/(w+1j*etan-w2)+1/(w+1j*etan-w2))
    else:
        Gimp_wi = a1*(1/(w+1j*etap-w1)+1/(w+1j*etap-w1))+a2*(1/(w+1j*etan-w2)+1/(w+1j*etap-w2))
#    print(Gimp_wi)
    print(w +mu - dw - 1/Gimp_wi[1][1])
#    print(w + mu - Sw_XY_t[-1])
#    print(dw + 1/Gimp_wi[1][1])
#    print(-np.log(U_t)/tau/1j)
    Sw_ZZ_t.append(w +mu - dw - 1/Gimp_wi[1][1])
 #   Aw_XY_t.append(rho0*(w + mu - Sw_XY_t[-1]))
    Aw_ZZ_t.append(rho0*(dw + 1/Gimp_wi[1][1]))
    HSIAM1_ZZ_t.append(-np.log(U1[1][1])/tau/1j)
    HSIAM2_ZZ_t.append(-np.log(U2)/tau/1j)
    HSIAM3_ZZ_t.append(-np.log(U3[0][0])/tau/1j)
    HSIAM_ZZ_t.append(-np.log(U_t)/tau/1j) 
    HSIAM1_ZZ_t_R.append(np.real(-np.log(U1[1][1])/tau/1j))
    HSIAM2_ZZ_t_R.append(np.real(-np.log(U2)/tau/1j))
    HSIAM3_ZZ_t_R.append(np.real(-np.log(U3[0][0])/tau/1j))
    HSIAM_ZZ_t_R.append(np.real(-np.log(U_t)/tau/1j))  
#second 

#...
# N = 100 #total 1D chain' Total time steps
# ts = 10 # [10, 100, 0.04, 100, 1.1]

# #CZ gate with ZZ interactions time domain
# #classically trial with CZ
# ec = 0
# iGimpR_tau = []
# Z_w = []
# Vnew_tau =[]
# HSIAM_tau = []
# Ts = []
# Tau = [] #

# H1 = []
# H2 = []
# H3 = []
# U_t = []
# tau_t = []

# N = 24
# tauc = 20
# #in this test, we assume  half-filled L = 14-site two-rung triangular Hubbard model. The system is initialized in the ground
# #state of H(0), setting U¯ = 5.0¯τ with the specified τ0
# L=4
# Au = -0.75
# Atau = Au/2
# Omega = 2.5 * tau0
# Tp = 0
# Tw = 10**6#infinite
# #delta = np.linspace(0, 12, 12)
# delta = 1

# #One spin-down(bath for site 1) and the other two spin-up(Impurity  for 2 and 3)=> One eigenstate, 
# # assume 3 qubits intialized as 0

# # star z = infinite
# # quantum:
# # ground state: gamma1 and gamm2(half-filled Hubbard model) Ja = Jb = const
# # steady state： Ja = Jb = 0 

# Sp1 = c1u_D*c1d #site2
# Sp2 = c2u_D*c2d #site3
# Sn1 = c1d_D*c1u #site2
# Sn2 = c2d_D*c2u  #site3
# Sz3 = n1u + n1d #site1

# #Sp = c1u_D*c1d + c2u_D*c2d
# #Sn = n1u+n2u-(n1d+n2d)
# #Sz = Sn
# etap = c1u_D*c1d_D - c2u_D*c2d_D
# etaz = n1u+n2u+n1d+n2d-2
# etan = c1d*c1u - c2d*c2u

# #stationary states
# rho1 = np.array([0,0,0])*np.array([0,0,0])
# rho2 = (np.array([0,0,1])-np.array([0,1,0]))*(np.array([0,0,1])-np.array([0,1,0]))/2
# A = (np.array([0,0,1])-np.array([0,1,0]))*np.array([0,0,0])

# for t in np.linspace(0, 20, N):
# #consider tripartite(XXZ spin ring) on non-relative Hilbert space:
#     tau0 =1                                                                                                                                                                                                                                                                             
#     taus = 0.2*t*tau0 #taus = [0.2, 0.4, 0.9, 1.4]*tau0
#     U0 =5*tau0

#     #AU and Aτ are the amplitudes of the modulation of U and τ relative to
#     #their equilibrium values U¯ and ¯τ . The frequency of the
#     #oscillations is Ω, whilst Tp and Tw describe the offset and
#     #width of the Gaussian envelope containing these oscillations
#     tau = tau0*(1+Atau*np.sin(Omega*t)*np.exp(-(t-Tp)**2/(2*Tw**2)))
#     U = U0*(1+Au*np.sin(Omega*t)*np.exp(-(t-Tp)**2/(2*Tw**2)))
#     tau_t.append(tau)
#     U_t.append(U)
    
#     # w1,2 = w
#     w = 1
#     w1 = w
#     w2 = w
    
#     delta = 1
#     B = 0.5
#     gamma = 2


#     #XXZ spin
#     St = (Sp1+Sp2)*(Sn1+Sn2)/(Sp1+Sp2+Sn1+Sn2)
#     Et = etap*etan/(etap+etan)
#     Dt = Et/N
#     if U >=0:
# #        H3.append(-tau*np.norm(sum(Sp2)))
#         H3.append(-tau*sum(sum(Sn1+Sn2)))
#         H2.append(-taus*sum(sum(Sp1+Sp2)))
#         H1.append(U*sum(sum(Sz3)))
#     else:
#         H3.append(-tau*sum(sum(etap)))
#         H2.append(-taus*sum(sum(etan)))
#         H1.append(U*sum(sum(etaz)))
#         if taus > tauc:
#             break
 
#CZ-gates (odd Trotter)
#sigma

Sw_ZZO_t = []
Aw_ZZO_t = []
HSIAM1_ZZO_t = []
HSIAM2_ZZO_t = []
HSIAM3_ZZO_t = []
HSIAM_ZZO_t = []
HSIAM1_ZZO_t_R = []
HSIAM2_ZZO_t_R = []
HSIAM3_ZZO_t_R = []
HSIAM_ZZO_t_R = []
U_ZZO_t_R = []
for t in np.linspace(0,tau,N+1):
    V = t
    GS = -1j*V/2*(Sigx12*Sigy12)/2/dw #4*4
    a1 = alpha(1, c1u_D, c1d_D, GS)[1][0] 
    a2 = alpha(2, c1u_D, c1d_D, GS)[1][0]
    U1 = Ra(-np.pi/2, 'x', 1)*Ra(-np.pi/2, 'x', 2)*Ra(-np.pi/2, 'x', 3)*Ra(-np.pi/2, 'x', 4)*np.exp(-1j*V/2*Sigx12)[0][1]*np.exp(-1j*V/2*Sigy12)[0][1]*Ra(-np.pi/2, 'x', 1)*Ra(-np.pi/2, 'x', 2)*Ra(-np.pi/2, 'x', 3)*Ra(-np.pi/2, 'x', 4)
    U2 = np.exp(-1j*U/4*Sigz13*t)[1][1]*np.exp(-1j*(U/4-mu/2)*sigma1z*t)[1][1]*np.exp(-1j*(U/4-mu/2)*sigma3z*t)[1][1]*np.exp(1j*ec/2*sigma2z*t)[1][1]*np.exp(1j*ec/2*sigma4z*t)[1][1]
    U3 = Ra(np.pi/2, 'y', 1)*Ra(np.pi/2, 'y', 2)*Ra(np.pi/2, 'y', 3)*Ra(np.pi/2, 'y', 4)*np.exp(-1j*V/2*Sigx12)[0][0]*np.exp(-1j*V/2*Sigy12)[0][0]*Ra(-np.pi/2, 'y', 1)*Ra(-np.pi/2, 'y', 2)*Ra(-np.pi/2, 'y', 3)*Ra(-np.pi/2, 'y', 4)
    U_ZZO_t = np.array([U1[0][0], U2, U3[1][1]])
    U_ZZO_t_R.append(np.array([np.real(U1[1][1]), np.real(U2), np.real(U3[0][0])]))
    if sum(U_t>0)<3:
        Gimp_wi = a1*(1/(w+1j*etan-w1)+1/(w+1j*etan-w1))+a2*(1/(w+1j*etan-w2)+1/(w+1j*etan-w2))
    else:
        Gimp_wi = a1*(1/(w+1j*etap-w1)+1/(w+1j*etap-w1))+a2*(1/(w+1j*etan-w2)+1/(w+1j*etap-w2))
#    print(Gimp_wi)
    print(w +mu - dw - 1/Gimp_wi[1][1])
#    print(w + mu - Sw_XY_t[-1])
#    print(dw + 1/Gimp_wi[1][1])
#    print(-np.log(U_t)/tau/1j)
    Sw_ZZO_t.append(w +mu - dw - 1/Gimp_wi[1][1])
 #   Aw_XY_t.append(rho0*(w + mu - Sw_XY_t[-1]))
    Aw_ZZO_t.append(rho0*(dw + 1/Gimp_wi[1][1]))
    HSIAM1_ZZO_t.append(-np.log(U1[0][0])/tau/1j)
    HSIAM2_ZZO_t.append(-np.log(U2)/tau/1j)
    HSIAM3_ZZO_t.append(-np.log(U3[1][1])/tau/1j)
    HSIAM_ZZO_t.append(-np.log(U_ZZO_t)/tau/1j) 
    HSIAM1_ZZO_t_R.append(np.real(-np.log(U1[0][0])/tau/1j))
    HSIAM2_ZZO_t_R.append(np.real(-np.log(U2)/tau/1j))
    HSIAM3_ZZO_t_R.append(np.real(-np.log(U3[1][1])/tau/1j))
    HSIAM_ZZO_t_R.append(np.real(-np.log(U_ZZO_t)/tau/1j))  
 
#XY-gates (odd Trotter)
#sigma

Sw_XYO_t = []
Aw_XYO_t = []
HSIAM1_XYO_t = []
HSIAM2_XYO_t = []
HSIAM_XYO_t = []
HSIAM1_XYO_t_R = []
HSIAM2_XYO_t_R = []                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
HSIAM_XYO_t_R = []
U_XYO_t_R = []
for t in np.linspace(0,tau,N+1):
    V = t
    GS = -1j*V/2*(Sigx12*Sigy12)/2/dw #4*4
    a1 = alpha(1, c1u_D, c1d_D, GS)[1][0] 
    a2 = alpha(2, c1u_D, c1d_D, GS)[1][0]
    U1 = Ra(-np.pi/2, 'x', 1)*Ra(-np.pi/2, 'y', 2)*Ra(-np.pi/2, 'x', 3)*Ra(-np.pi/2, 'y', 4)*np.exp(-1j*V/2*Sigz13)[0][1]*np.exp(-1j*V/2*Sigz13)[0][1]*Ra(-np.pi/2, 'x', 1)*Ra(-np.pi/2, 'x', 2)*Ra(-np.pi/2, 'x', 3)*Ra(-np.pi/2, 'x', 4)
    U2 = np.exp(-1j*U/4*Sigz13*t)[1][1]*np.exp(-1j*(U/4-mu/2)*sigma1z*t)[1][1]*np.exp(-1j*(U/4-mu/2)*sigma3z*t)[1][1]*np.exp(1j*epsc/2*sigma2z*t)[1][1]*np.exp(1j*epsc/2*sigma4z*t)[1][1]
    U_XYO_t = np.array([U1[0][0], U2])
    U_XYO_t_R.append(np.array([np.real(U1[1][1]), np.real(U2)]))
    if sum(U_XYO_t>0)<3:
        Gimp_wi = a1*(1/(w+1j*etan-w1)+1/(w+1j*etan-w1))+a2*(1/(w+1j*etan-w2)+1/(w+1j*etan-w2))
    else:
        Gimp_wi = a1*(1/(w+1j*etap-w1)+1/(w+1j*etap-w1))+a2*(1/(w+1j*etan-w2)+1/(w+1j*etap-w2))
#    print(Gimp_wi)
    print(w +mu - dw - 1/Gimp_wi[1][1])
#    print(w + mu - Sw_XY_t[-1])
#    print(dw + 1/Gimp_wi[1][1])
#    print(-np.log(U_t)/tau/1j)
    Sw_XYO_t.append(w +mu - dw - 1/Gimp_wi[1][1])
 #   Aw_XY_t.append(rho0*(w + mu - Sw_XY_t[-1]))
    Aw_XYO_t.append(rho0*(dw + 1/Gimp_wi[1][1]))
    HSIAM1_XYO_t.append(-np.log(U1[0][0])/tau/1j)
    HSIAM2_XYO_t.append(-np.log(U2)/tau/1j)
    HSIAM_XYO_t.append(-np.log(U_XYO_t)/tau/1j) 
    HSIAM1_XYO_t_R.append(np.real(-np.log(U1[0][0])/tau/1j))
    HSIAM2_XYO_t_R.append(np.real(-np.log(U2)/tau/1j))
    HSIAM_XYO_t_R.append(np.real(-np.log(U_XYO_t)/tau/1j))  

     
#result visualisation
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
from itertools import combinations


HSIAM1_ZZO_t_M = np.zeros((np.shape(HSIAM1_ZZO_t)))
HSIAM2_ZZO_t_M = np.zeros((np.shape(HSIAM2_ZZO_t)))
HSIAM3_ZZO_t_M = np.zeros((np.shape(HSIAM3_ZZO_t)))
HSIAM_ZZO_t_M = np.zeros((np.shape(HSIAM_ZZO_t)))
HSIAM1_ZZO_t_R_M = np.zeros((np.shape(HSIAM1_ZZO_t_R)))
HSIAM2_ZZO_t_R_M = np.zeros((np.shape(HSIAM2_ZZO_t_R)))
HSIAM3_ZZO_t_R_M = np.zeros((np.shape(HSIAM3_ZZO_t_R)))
HSIAM_ZZO_t_R_M = np.zeros((np.shape(HSIAM_ZZO_t_R)))
U_ZZO_t_R_M = np.zeros((np.shape(U_ZZO_t_R)))
HSIAM1_ZZO_t_std = np.ones((np.shape(HSIAM1_ZZO_t)))
HSIAM2_ZZO_t_std = np.ones((np.shape(HSIAM2_ZZO_t)))
HSIAM3_ZZO_t_std = np.ones((np.shape(HSIAM3_ZZO_t)))
HSIAM_ZZO_t_std = np.ones((np.shape(HSIAM_ZZO_t)))
HSIAM1_ZZO_t_R_std = np.ones((np.shape(HSIAM1_ZZO_t_R)))
HSIAM2_ZZO_t_R_std = np.ones((np.shape(HSIAM2_ZZO_t_R)))
HSIAM3_ZZO_t_R_std = np.ones((np.shape(HSIAM3_ZZO_t_R)))
HSIAM_ZZO_t_R_std = np.ones((np.shape(HSIAM_ZZO_t_R)))
U_ZZO_t_R_std = np.ones((np.shape(U_ZZO_t_R)))

for tst in range(1,N):
    HSIAM1_ZZO_t_M[tst] = np.mean(HSIAM1_ZZO_t[0:tst])
    HSIAM2_ZZO_t_M[tst] = np.mean(HSIAM2_ZZO_t[0:tst])
    HSIAM3_ZZO_t_M[tst] = np.mean(HSIAM3_ZZO_t[0:tst])
    HSIAM_ZZO_t_M[tst] = np.mean(HSIAM_ZZO_t[0:tst])
    HSIAM1_ZZO_t_R_M[tst] = np.mean(HSIAM1_ZZO_t_R[0:tst])
    HSIAM2_ZZO_t_R_M[tst] = np.mean(HSIAM2_ZZO_t_R[0:tst])
    HSIAM3_ZZO_t_R_M[tst] = np.mean(HSIAM3_ZZO_t_R[0:tst])
    HSIAM_ZZO_t_R_M[tst] = np.mean(HSIAM_ZZO_t_R[0:tst])
    U_ZZO_t_R_M[tst] = np.mean(U_ZZO_t_R[0:tst])
    HSIAM1_ZZO_t_std[tst] = np.std(HSIAM1_ZZO_t_M[0:tst])
    HSIAM2_ZZO_t_std[tst] = np.std(np.shape(HSIAM2_ZZO_t_M[0:tst]))
    HSIAM3_ZZO_t_std[tst] = np.std(np.shape(HSIAM3_ZZO_t_M[0:tst]))
    HSIAM_ZZO_t_std[tst] = np.std(np.shape(HSIAM_ZZO_t_M[0:tst]))
    HSIAM1_ZZO_t_R_std[tst] = np.std(np.shape(HSIAM1_ZZO_t_R_M[0:tst]))
    HSIAM2_ZZO_t_R_std[tst] = np.std(np.shape(HSIAM2_ZZO_t_R_M[0:tst]))
    HSIAM3_ZZO_t_R_std[tst] = np.std(np.shape(HSIAM3_ZZO_t_R_M[0:tst]))
    HSIAM_ZZO_t_R_std[tst] = np.std(np.shape(HSIAM_ZZO_t_R_M[0:tst]))
    U_ZZO_t_R_std[tst] = np.std(np.shape(U_ZZO_t_R_M[0:tst]))

HSIAM1_ZZO_t_std[HSIAM1_ZZO_t_std==0] = 10**(-5)
HSIAM2_ZZO_t_std[HSIAM2_ZZO_t_std==0] = 10**(-5)
HSIAM3_ZZO_t_std[HSIAM3_ZZO_t_std==0] = 10**(-5)
HSIAM_ZZO_t_std[HSIAM_ZZO_t_std==0] = 10**(-5)
HSIAM1_ZZO_t_R_std[HSIAM1_ZZO_t_R_std==0] = 10**(-5)
HSIAM2_ZZO_t_R_std[HSIAM2_ZZO_t_R_std==0] = 10**(-5)
HSIAM3_ZZO_t_R_std[HSIAM3_ZZO_t_R_std==0] = 10**(-5)
HSIAM_ZZO_t_R_std[HSIAM_ZZO_t_R_std==0] = 10**(-5)
U_ZZO_t_R_std[U_ZZO_t_R_std==0] = 10**(-5)
    
HSIAM1_ZZO_t_CV = np.array(HSIAM1_ZZO_t_M, dtype = float)/np.array(HSIAM1_ZZO_t_R_std, dtype = float)
HSIAM2_ZZO_t_CV = np.array(HSIAM2_ZZO_t_M, dtype = float)/np.array(HSIAM2_ZZO_t_R_std, dtype = float)
HSIAM3_ZZO_t_CV = np.array(HSIAM3_ZZO_t_M, dtype = float)/np.array(HSIAM3_ZZO_t_R_std, dtype = float)
HSIAM_ZZO_t_CV = np.array(HSIAM_ZZO_t_M, dtype = float)/np.array(HSIAM_ZZO_t_R_std, dtype = float)
HSIAM1_ZZO_R_CV = np.array(HSIAM1_ZZO_t_M, dtype = float)/np.array(HSIAM1_ZZO_t_R_std, dtype = float)
HSIAM2_ZZO_R_CV = np.array(HSIAM2_ZZO_t_M, dtype = float)/np.array(HSIAM2_ZZO_t_R_std, dtype = float)
HSIAM3_ZZO_R_CV = np.array(HSIAM3_ZZO_t_M, dtype = float)/np.array(HSIAM3_ZZO_t_R_std, dtype = float)
HSIAM_ZZO_R_CV = np.array(HSIAM_ZZO_t_M, dtype = float)/np.array(HSIAM_ZZO_t_R_std, dtype = float)
U_ZZO_t_R_CV = np.array(U_ZZO_t_R_M, dtype = float)/np.array(U_ZZO_t_R_std, dtype = float)

#HSIAM_tau_CV = np.array(HSIAM_tau_M, dtype = float)/np.array(HSIAM_tau_std, dtype = float)
plt.figure()
plt.plot(HSIAM1_ZZO_t_CV)
plt.title('HSIAM1_ZZ_odd_t_CV')
plt.figure()
plt.plot(HSIAM2_ZZO_t_CV)
plt.title('HSIAM2_ZZ_odd_t_CV')
plt.figure()
plt.plot(HSIAM3_ZZO_t_CV)
plt.title('HSIAM3_ZZ_odd_t_CV')
plt.figure()
plt.plot(HSIAM_ZZO_t_CV)
plt.title('HSIAM_ZZ_odd_t_CV')


plt.figure()
plt.plot(HSIAM1_ZZO_t_M)
plt.title('HSIAM1_ZZ_odd_t_M')
plt.figure()
plt.plot(HSIAM2_ZZO_t_M)
plt.title('HSIAM2_ZZ_odd_t_M')
plt.figure()
plt.plot(HSIAM3_ZZO_t_M)
plt.title('HSIAM3_ZZ_odd_t_M')
plt.figure()
plt.plot(HSIAM_ZZO_t_M)
plt.title('HSIAM_ZZ_odd_t_M')

plt.figure()
plt.plot(HSIAM1_ZZO_R_CV)
plt.title('HSIAM1_ZZ_odd_R_CV')
plt.figure()
plt.plot(HSIAM2_ZZO_R_CV)
plt.title('HSIAM2_ZZ_odd_R_CV')
plt.figure()
plt.plot(HSIAM3_ZZO_R_CV)
plt.title('HSIAM3_ZZ_odd_R_CV')
plt.figure()
plt.plot(HSIAM_ZZO_R_CV)
plt.title('HSIAM_ZZ_odd_R_CV')

plt.figure()
plt.plot(HSIAM1_ZZO_t_R_M)
plt.title('HSIAM1_ZZ_odd_R_M')
plt.figure()
plt.plot(HSIAM2_ZZO_t_R_M)
plt.title('HSIAM2_ZZ_odd_R_M')
plt.figure()
plt.plot(HSIAM3_ZZO_t_R_M)
plt.title('HSIAM3_ZZ_odd_R_M')
plt.figure()
plt.plot(HSIAM_ZZO_t_R_M)
plt.title('HSIAM_ZZ_odd_R_M')

plt.figure()
plt.errorbar(x = range(N+1), y = HSIAM1_ZZO_R_CV, yerr = HSIAM1_ZZO_t_R_std)
plt.title('HSIAM1_ZZ_odd_R_CV_err')
plt.figure()
plt.errorbar(x = range(N+1), y = HSIAM2_ZZO_R_CV, yerr = HSIAM2_ZZO_t_R_std)
plt.title('HSIAM2_ZZ_odd_R_CV_err')
plt.figure()
plt.errorbar(x = range(N+1), y = HSIAM3_ZZO_R_CV, yerr= HSIAM3_ZZO_t_R_std)
plt.title('HSIAM3_ZZ_odd_R_CV_err')
plt.figure()
plt.errorbar(x = range(N+1), y = HSIAM_ZZO_R_CV, yerr= HSIAM_ZZO_t_R_std)
plt.title('HSIAM_ZZ_odd_R_CV_err')

plt.figure()
plt.contour(np.reshape(U_ZZO_t_R,[N+1,3])) 
plt.ylabel('U_ZZ_odd_t_R')
plt.xlabel('tau')
plt.title('Countour of U_ZZ_odd_t_R versus tau')
plt.figure()
plt.contour(np.reshape(HSIAM_ZZO_t_R,[N+1,3])) 
plt.ylabel('HSIAM_ZZ_odd_t_R')
plt.xlabel('tau')
plt.title('Countour of HSIAM_ZZ_odd_t_R versus tau')

HSIAM1_XYO_t_M = np.zeros((np.shape(HSIAM1_XYO_t)))
HSIAM2_XYO_t_M = np.zeros((np.shape(HSIAM2_XYO_t)))
HSIAM_XYO_t_M = np.zeros((np.shape(HSIAM_XYO_t)))
HSIAM1_XYO_t_R_M = np.zeros((np.shape(HSIAM1_XYO_t_R)))
HSIAM2_XYO_t_R_M = np.zeros((np.shape(HSIAM2_XYO_t_R)))
HSIAM_XYO_t_R_M = np.zeros((np.shape(HSIAM_XYO_t_R)))
U_XYO_t_R_M = np.zeros((np.shape(U_XYO_t_R)))
HSIAM1_XYO_t_std = np.ones((np.shape(HSIAM1_XYO_t)))
HSIAM2_XYO_t_std = np.ones((np.shape(HSIAM2_XYO_t)))
HSIAM_XYO_t_std = np.ones((np.shape(HSIAM_XYO_t)))
HSIAM1_XYO_t_R_std = np.ones((np.shape(HSIAM1_XYO_t_R)))
HSIAM2_XYO_t_R_std = np.ones((np.shape(HSIAM2_XYO_t_R)))
HSIAM_XYO_t_R_std = np.ones((np.shape(HSIAM_XYO_t_R)))
U_XYO_t_R_std = np.ones((np.shape(U_XYO_t_R)))

for tst in range(1,N):
    HSIAM1_XYO_t_M[tst] = np.mean(HSIAM1_XYO_t[0:tst])
    HSIAM2_XYO_t_M[tst] = np.mean(HSIAM2_XYO_t[0:tst])
    HSIAM_XYO_t_M[tst] = np.mean(HSIAM_XYO_t[0:tst])
    HSIAM1_XYO_t_R_M[tst] = np.mean(HSIAM1_XYO_t_R[0:tst])
    HSIAM2_XYO_t_R_M[tst] = np.mean(HSIAM2_XYO_t_R[0:tst])
    HSIAM_XYO_t_R_M[tst] = np.mean(HSIAM_XYO_t_R[0:tst])
    U_XYO_t_R_M[tst] = np.mean(U_XYO_t_R[0:tst])
    HSIAM1_XYO_t_std[tst] = np.std(HSIAM1_XYO_t_M[0:tst])
    HSIAM2_XYO_t_std[tst] = np.std(np.shape(HSIAM2_XYO_t_M[0:tst]))
    HSIAM_XYO_t_std[tst] = np.std(np.shape(HSIAM_XYO_t_M[0:tst]))
    HSIAM1_XYO_t_R_std[tst] = np.std(np.shape(HSIAM1_XYO_t_R_M[0:tst]))
    HSIAM2_XYO_t_R_std[tst] = np.std(np.shape(HSIAM2_XYO_t_R_M[0:tst]))
    HSIAM_XYO_t_R_std[tst] = np.std(np.shape(HSIAM_XYO_t_R_M[0:tst]))
    U_XYO_t_R_std[tst] = np.std(np.shape(U_XYO_t_R_M[0:tst]))

HSIAM1_XYO_t_std[HSIAM1_XYO_t_std==0] = 10**(-5)
HSIAM2_XYO_t_std[HSIAM2_XYO_t_std==0] = 10**(-5)
HSIAM_XYO_t_std[HSIAM_XYO_t_std==0] = 10**(-5)
HSIAM1_XYO_t_R_std[HSIAM1_XYO_t_R_std==0] = 10**(-5)
HSIAM2_XYO_t_R_std[HSIAM2_XYO_t_R_std==0] = 10**(-5)
HSIAM_XYO_t_R_std[HSIAM_XYO_t_R_std==0] = 10**(-5)
U_XYO_t_R_std[U_XYO_t_R_std==0] = 10**(-5)
    
HSIAM1_XYO_t_CV = np.array(HSIAM1_XYO_t_M, dtype = float)/np.array(HSIAM1_XYO_t_R_std, dtype = float)
HSIAM2_XYO_t_CV = np.array(HSIAM2_XYO_t_M, dtype = float)/np.array(HSIAM2_XYO_t_R_std, dtype = float)
HSIAM_XYO_t_CV = np.array(HSIAM_XYO_t_M, dtype = float)/np.array(HSIAM_XYO_t_R_std, dtype = float)
HSIAM1_XYO_R_CV = np.array(HSIAM1_XYO_t_M, dtype = float)/np.array(HSIAM1_XYO_t_R_std, dtype = float)
HSIAM2_XYO_R_CV = np.array(HSIAM2_XYO_t_M, dtype = float)/np.array(HSIAM2_XYO_t_R_std, dtype = float)
HSIAM_XYO_R_CV = np.array(HSIAM_XYO_t_M, dtype = float)/np.array(HSIAM_XYO_t_R_std, dtype = float)
U_XYO_t_R_CV = np.array(U_XYO_t_R_M, dtype = float)/np.array(U_XYO_t_R_std, dtype = float)

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
plt.errorbar(x = range(N+1), y = HSIAM1_XYO_R_CV, yerr = HSIAM1_XYO_t_R_std)
plt.title('HSIAM1_XYO_odd_R_CV_err')
plt.figure()
plt.errorbar(x = range(N+1), y = HSIAM2_XYO_R_CV, yerr = HSIAM2_XYO_t_R_std)
plt.figure()
plt.errorbar(x = range(N+1), y = HSIAM_XYO_R_CV, yerr= HSIAM_XYO_t_R_std)
plt.title('HSIAM_XYO_odd_R_CV_err')

plt.figure()
plt.contour(np.reshape(U_XYO_t_R,[N+1,2])) 
plt.ylabel('U_XY_odd_t_R')
plt.xlabel('tau')
plt.title('Countour of U_XY_odd_t_R versus tau')
plt.figure()
plt.contour(np.reshape(HSIAM_XYO_t_R,[N+1,2])) 
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
plt.ylabel('Aw_XYO_t')
plt.xlabel('tau')
plt.title('PolarPlot of a_XY_odd_t versus tau')


from matplotlib.ticker import LinearLocator

#fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
X = np.sort(np.reshape(U_t_R,[N+1,3])[:][0])
Y = np.sort(U_t_R[0])
Z = np.reshape(U_t_R,[N+1,3])#Z, positive
#Z = np.sqrt(np.array(np.real(Z[0]))**2+np.array(np.real(Z[1]))**2)

# Plot the surface.
surf = ax.plot_surface(X,Y,Z, cmap='summer',
                       linewidth=0, antialiased=False)

surf = ax.plot3d(Z, cmap='summer',
                       linewidth=0, antialiased=False)

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)
plt.title('U_tau_positive')
plt.show()

from matplotlib.ticker import LinearLocator
import matplotlib.pyplot as plt
fig = plt.figure()
ax = plt.axes(projection='3d')
#fig1, ax1 = plt.subplots(subplot_kw={"projection": "3d"})
X = np.reshape(U_t_R,[N+1,3])[:][0]
Y = range(len(np.reshape(U_t_R,[N+1,3])[0][:]))
Z = np.reshape(U_t_R,[N+1,3])#Z, positive


# Plot the surface.
#surf1 = ax1.plot_surface(X, Y, Z, cmap='summer',
#                       linewidth=0, antialiased=False)

# Add a color bar which maps values to colors.
#fig1.colorbar(surf1, shrink=0.5, aspect=5)
ax.plot_surface(X,Y,Z)
plt.title('U_tau_negative')
plt.show()


from matplotlib.ticker import LinearLocator

fig2, ax2 = plt.subplots(subplot_kw={"projection": "3d"})
X = np.sort(tau_t)
Y = np.sort(U_t)
Z = np.meshgrid(H2, H3)#positive, negative
Z = np.sqrt(np.array(np.real(Z[0]))**2+np.array(np.real(Z[1]))**2)

# Plot the surface.
surf2 = ax2.plot_surface(X, Y, Z, cmap='summer',
                       linewidth=0, antialiased=False)

# Add a color bar which maps values to colors.
fig2.colorbar(surf2, shrink=0.5, aspect=5)
plt.title('U_tau_positive_negative')
plt.show()

