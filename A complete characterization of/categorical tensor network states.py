# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np
import matplotlib.pyplot as plt

def AND(a,b):
    A000 = 0
    A010 = 0
    A100 = 0
    A111 = 0
    for n in range(len(a)):
        if [a[n],b[n]] == [0, 0]:
            AND = 0
            A000 += 1
        elif [a[n],b[n]] == [0, 1]:
            AND = 0
            A010 += 1
        elif [a[n],b[n]] == [1, 0]:
            AND = 0
            A100 += 1
        elif [a[n], b[n]] == [1, 1]:
            AND = 1
            A111 += 1
    return AND

def OR(a,b):
    A001 = 0
    A011 = 0
    A101 = 0
    A111 = 0
    for n in range(len(a)):
        if [a[n],b[n]] == [0, 0]:
            OR = 0
            A001 += 1
        elif [a[n],b[n]] == [0, 1]:
            OR = 1
            A011 += 1
        elif [a[n],b[n]] == [1, 0]:
            OR = 1
            A101 += 1
        elif [a[n],b[n]] == [1, 1]:
            OR = 1
            A111 += 1
    return OR,[A001, A011, A101, A111]
        
def XOR(a,b):#note that XOR is used in negation of a boolean x = 1[+] x
    A000 = 0
    A011 = 0
    A101 = 0
    A110 = 0
    for n in range(len(a)):
        if [a[n],b[n]] == [0, 0]:
            XOR = 0
            A000 += 1
        elif [a[n],b[n]] == [0, 1]:
            OR = 1
            A011 += 1
        elif [a[n],b[n]] == [1, 0]:
            XOR = 1
            A101 += 1
        elif [a[n], b[n]] == [1, 1]:
            XOR = 0
            A110 += 1
    return XOR,[A000, A011, A101, A110]
  
def XNOR(a,b):
    A001 = 0
    A010 = 0
    A100 = 0
    A111 = 0
    for n in range(len(a)):
        if [a[n],b[n]] == [0, 0]:
            XNOR = 1
            A001 += 1
        elif [a[n],b[n]] == [0, 1]:
            XNOR = 0
            A010 += 1
        elif [a[n],b[n]] == [1, 0]:
            XNOR = 0
            A100 += 1
        elif [a[n], b[n]] == [1, 1]:
            XNOR = 1
            A111 += 1
    return XNOR,[A001, A010, A100, A111]

def fWB(a,b,c):
    for n in range(len(a)):
        if [a[n],b[n],c[n]] == [0, 0, 0]:
            fWB = 0
        elif [a[n],b[n],c[n]] == [0, 0, 1]:
            fWB = 1
        elif [a[n],b[n],c[n]] == [0, 1, 0]:
            fWB = 1
        elif [a[n],b[n],c[n]] == [0, 1, 1]:
            fWB = 0
        elif [a[n],b[n],c[n]] == [1, 1, 0]:
            fWB = 0
        elif [a[n],b[n],c[n]] == [1, 1, 1]:
            fWB = 0
        elif [a[n],b[n],c[n]] == [1, 0, 0]:
            fWB = 1        
        elif [a[n],b[n],c[n]] == [1, 0, 1]:
            fWB = 0
    return fWB

def fGHZB(a,b,c):
    for n in range(len(a)):
        if [a[n],b[n],c[n]] == [0, 0, 0]:
            fGHZB = 1
        elif [a[n],b[n],c[n]] == [0, 0, 1]:
            fGHZB = 0
        elif [a[n],b[n],c[n]] == [0, 1, 0]:
            fGHZB = 1
        elif [a[n],b[n],c[n]] == [0, 1, 1]:
            fGHZB = 0
        elif [a[n],b[n],c[n]] == [1, 1, 0]:
            fGHZB = 1
        elif [a[n],b[n],c[n]] == [1, 1, 1]:
            fGHZB = 1
        elif [a[n],b[n],c[n]] == [1, 0, 0]:
            fGHZB = 1        
        elif [a[n],b[n],c[n]] == [1, 0, 1]:
            fGHZB = 0
    return fGHZB 
        
def fWS(a,b,c):
    fWS = 1-a-b+a*b-c+a*c+b*c
    return fWS

def fGHZS(a,b,c):
    fGHZS = a+b+c-2*a*b*2*a*c*2*b*c+3*a*b*c
    return fGHZS

def CTNRS(a,b):
    if np.size(a) == 1:
        ak = 1
        if [a, b] in [[0,0], [1,1]]: 
            CTNRS = 1 
        else:
            CTNRS = 0
    else:
        ak = len(a)
        for n in range(len(a)):
            if [a[n], b[n]] in [[0,0], [1,1]]: 
                CTNRS = XOR(a,b)[0]
                ak -= 1 
            else:
                CTNRS = 0
    return ak*CTNRS

def normalize(a):
    if np.linalg.norm(a)/len(a) >=0.5:
        b = 1
    else:
        b = 0
    return b

# test with [X1, Z2, Z3, Z4; Z1, X2, I3, Z4; Z1, I2, X3, Z4;Z1, Z2, Z3, X4];
NEIBORM= [[1, 1, 1, 1],[1, 1, 0, 1],[1, 0, 1, 1],[1, 1, 1, 1]]
AvN1= [[1, -1, -1, -1],[-1, 0, 1, -1],[-1, -1, -1, 1]]
AvN2= [[1, -1, -1, -1],[-1, 1, 0, -1],[-1, -1, -1, 1]]
AvN3= [[1, -1, -1, -1],[0.5, 0.5, -1, 0],[0.5, -1, 0.5, 0]]
AvN4= [[-1, -1, -1, 1],[0, 0.5, -1, 0.5],[0, -1, 0.5, 0.5]]
CTNRS(normalize(AvN1[0]), normalize(AvN1[1])) == normalize(AvN1[2])
CTNRS(normalize(AvN2[0]), normalize(AvN2[1])) == normalize(AvN2[2])
CTNRS(normalize(AvN3[0]), normalize(AvN3[1])) == normalize(AvN3[2])
CTNRS(normalize(AvN4[0]), normalize(AvN4[1])) == normalize(AvN4[2])

AVN30 = np.repeat(normalize(AvN3[0]),100)
AVN20 = np.repeat(normalize(AvN2[0]),100)
AvN31 = np.repeat(AvN3[1],100).reshape(100,4)
AvN21 = np.repeat(AvN2[1],100).reshape(100,4)
AvN32 = np.repeat(AvN3[2],100).reshape(100,4)
AvN22 = np.repeat(AvN2[2],100).reshape(100,4)
for steps in range(100):   
    AvN31[steps,0] = np.ceil(np.random.normal(0.5,0.1))
    AvN31[steps,1] = np.ceil(np.random.normal(0.5,0.1))
    AvN32[steps,0] = AvN31[steps,0]
    AvN32[steps,2] = AvN31[steps,1]
    AvN21[steps,0] = np.ceil(np.random.normal(0.5,0.1))
    AvN21[steps,1] = np.ceil(np.random.normal(0.5,0.1))
    AvN22[steps,0] = AvN21[steps,0]
    AvN22[steps,2] = AvN21[steps,1]
 
#temp3 = np.repeat(normalize(AvN3[0]),100)
#for s in range(len(temp3)):
#    temp3[s] = np.ceil(np.random.normal(0.5,0.1))
N3 = CTNRS(AVN30, AvN31[:,0]) == normalize(AvN32)
N2 = CTNRS(AVN20, AvN21[:,0]) == normalize(AvN22)
     
plt.figure()
plt.plot(temp3)  
plt.ylabel('AvN31')
plt.xlabel('0')

     
plt.figure()
plt.plot(temp2)  
plt.ylabel('AvN21')
plt.xlabel('0')

from scipy import stats
t23, o23 = stats.spearmanr(np.array([temp2, temp3]).flatten(),np.array([origin2, origin3]).flatten())
rho23, p23 = stats.spearmanr(np.array([temp2, origin2]).flatten(),np.array([temp3, origin3]).flatten())

#t23=0.027409372236958444
#o23=0.7000375367247231
#rho23=0.07272247479895004
#p23=0.3061374773622768

#test for preprocessed signal on band Alpha and band Beta
restemp23 = stats.linregress(temp2.flatten(), temp3.flatten())
#LinregressResult(slope=0.12499999999999944, intercept=0.7500000000000004, rvalue=0.11706544723733367, pvalue=0.2460774319063414, stderr=0.10712030787245214)resori23 = stats.linregress(origin2.flatten(), origin3.flatten())
resorigin23 = stats.linregress(origin2.flatten(), origin3.flatten())
#LinregressResult(slope=0.030303030303030346, intercept=0.8333333333333333, rvalue=0.02837950236056591, pvalue=0.7792600307394673, stderr=0.10781850024443329)resTb = stats.linregress(ASMTb.flatten(), ASMb.flatten())

plt.plot(temp2.flatten(), temp3.flatten(), 'o', label='original data')
plt.plot(temp2.flatten(), restemp23.intercept + restemp23.slope*temp3.flatten(), 'r', label='fitted line')
plt.legend()
plt.title('dot 2 and dot 3 errorcorrected')
plt.show()

plt.plot(origin2.flatten(), origin3.flatten(), 'o', label='original data')
plt.plot(origin2.flatten(), resorigin23.intercept + restemp23.slope*origin3.flatten(), 'r', label='fitted line')
plt.legend()
plt.title('dot 2 and dot 3 origin')
plt.show()



#Calculate 95% confidence interval on slope and intercept:
from scipy.stats import t
tinv = lambda p, df: abs(t.ppf(p/2, df))
ts = tinv(0.05, len(ASMFRa.flatten())-2)
print(f"slope (95%): {resFRa.slope:.6f} +/- {ts*resFRa.stderr:.6f}")
#slope (95%): 1.243373 +/- 0.020287
print(f"intercept (95%): {resFRa.intercept:.6f}"
      f" +/- {ts*resFRa.stderr:.6f}")
#intercept (95%): -1851.389757 +/- 0.010143

tinv = lambda p, df: abs(t.ppf(p/2, df))
ts = tinv(0.05, len(ASMTa.flatten())-2)
print(f"slope (95%): {resTa.slope:.6f} +/- {ts*resTa.stderr:.6f}")
#slope (95%): 2.467063 +/- 0.122289
print(f"intercept (95%): {resTa.intercept:.6f}"
      f" +/- {ts*resTsa.tderr:.6f}")
#intercept (95%): -4371.021167 +/- 0.061145

tinv = lambda p, df: abs(t.ppf(p/2, df))
ts = tinv(0.05, len(ASMFRb.flatten())-2)
print(f"slope (95%): {resFRb.slope:.6f} +/- {ts*resFRb.stderr:.6f}")
#slope (95%): slope (95%): 1.222878 +/- 0.066106
print(f"intercept (95%): {resFRb.intercept:.6f}"
      f" +/- {ts*resFRb.stderr:.6f}")
#intercept (95%): -10985.633548 +/- 0.033053

tinv = lambda p, df: abs(t.ppf(p/2, df))
ts = tinv(0.05, len(ASMTb.flatten())-2)
print(f"slope (95%): {resTb.slope:.6f} +/- {ts*resTb.stderr:.6f}")
#slope (95%): 1.188525 +/- 0.055916
print(f"intercept (95%): {resTb.intercept:.6f}"
      f" +/- {ts*resTb.stderr:.6f}")
#intercept (95%): -16654.500426 +/- 0.055916