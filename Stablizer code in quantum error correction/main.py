'''

                            Online Python Compiler.
                Code, Compile, Run and Debug python program online.
Write your code in this editor and press "Run" button to execute it.

'''

def KO(pa):#Kraus Operator 
    if pa == 'x': KO = [[0,1],[1,0]] 
    if pa == 'y': KO = [[0,-1j],[1j,0]] 
    if pa == 'z': KO = [[1,0],[0,-1]] 
    KO =np.array(KO) 
    return KO
    
print('\n Stablizer Code')
print('--------------')
kB = 1
ka = 0.02
dT = 1/0.2 #dt/T
alpha = ka*np.exp(-dT)/(1+np.exp(-dT))
print('default alpha describes finite tamperature: ')
print(alpha)
beta = ka/(1+np.exp(-dT))
print('default beta describes finite tamperature: ')
print(beta)t = 10**(-5)#default alpha describes finite tamperature: #0.00013385701848569713
#default beta describes finite tamperature: #0.019866142981514304
Pu = beta/(alpha+beta+0.0001)
Pd = alpha/(alpha+beta+0.0001)
krsx = KO('x')
krsy = KO('y')
krsz = KO('z')
ii = np.array([[1,0],[0,1]]) #
gamma = np.linspace(0.1,1,100)
N = np.linspace(1,51,50)
i = 0F = []F0 = []Fdec = []Fdep = []Fdepl = []FT = []
if len(N)>1: 
while N[i]<51: 
Ni = N[i] 
pse = 0.99*Ni 
gamma = pse*alpha 
zs = (beta - alpha)/(alpha + beta+0.0001)
rou0 = 0.5*(ii + zs*krsz + zs*krsy + zs*krsx) 
rou = 0.5*(ii + zs*krsz + np.exp((alpha+beta+2*gamma)*t/2)*(zs*krsy + zs*krsx)+np.exp(-(alpha+beta+0.0001)*t)*(krsz-zs)*krsz) #assume for single-qubit channel with dimension db=2 db = 2 
Sol =0.5*gamma* sum(sum(krsz*rou*krsz-rou))#superoperator
lambda1 = 0.25*(1+2*np.exp((-alpha+beta+2*gamma)*t/2)+np.exp((-alpha+beta)*t)) lambda2 = 0.25*(1-np.exp((-alpha+beta)*t)) lambda3 = 0.25*(1-np.exp((-alpha+beta)*t)) lambda4 = 0.25*(1-2*np.exp((-alpha+beta+2*gamma)*t/2)+np.exp((-alpha+beta)*t)) mu = krsz/4*(1-np.exp((-alpha+beta)*t)) x = np.sqrt(4*mu**2+(lambda1-lambda4)**2) 
gamma1 = 0.5*(x-2*mu-lambda1+lambda4)*np.sqrt((lambda1+lambda3-x)/(4*mu**2-(lambda1-lambda4)*(x-lambda1+lambda4))) 
gamma2 = 0.5*(x+2*mu-lambda1+lambda4)*np.sqrt((lambda1+lambda3-x)/(4*mu**2-(lambda1-lambda4)*(x-lambda1+lambda4)+0.0001)) 
gamma3 = 0.5*(x+2*mu+lambda1-lambda4)*np.sqrt((lambda1+lambda3+x)/(4*mu**2+(lambda1-lambda4)*(x+lambda1-lambda4)+0.0001)) 
gamma4 = 0.5*(x-2*mu+lambda1-lambda4)*np.sqrt((lambda1+lambda3+x)/(4*mu**2+(lambda1-lambda4)*(x+lambda1-lambda4))) 
E1 = np.array([[gamma1,0],[0, gamma2]]) 
E2 = np.sqrt(Pu)*np.array([[0,np.sqrt(1-np.exp(-(alpha+beta)*t))],[0,0]]) 
E3 = [[gamma3,0],[0, gamma4]] 
E4 = np.sqrt(Pd)*np.array([[0,0],[np.sqrt(1-np.exp(-(alpha+beta)*t)),0]]) #F = (abs(E1)**2+abs(E2)**2+abs(E3)**2+abs(E4)**2+d)/db/(db+1) F.append(np.min((2+(lambda1-lambda4+x)**2*(lambda1+lambda4+x)/(4*mu**2+(lambda1-lambda4)*(lambda1-lambda4+x))+(lambda1-lambda4+x)**2*abs(lambda1+lambda4-x)/abs(4*mu**2-(lambda1-lambda4)*(lambda1-lambda4+x)))/6)) 
print('Fidelity describes general solution finite tamperature: ') 
print(F) 
# special cases 
# mu = 0 
E1 = np.sqrt(lambda1)*ii 
E2 = np.sqrt(lambda2)*krsx 
E3 = np.sqrt(lambda3)*krsy 
E4 = np.sqrt(lambda4)*krsz 
#F = (abs(E1)**2+abs(E2)**2+abs(E3)**2+abs(E4)**2+d)/db/(db+1) F0.append((4*lambda1+2)/6) print('Fidelity describes mu=0 finite tamperature: ') 
print(F0)
Rou.append(0.5*(ii + zs*krsz + np.exp((alpha+beta+2*gamma)*t/2)*(zs*krsy + zs*krsx)+np.exp(-(alpha+beta+0.0001)*t)*(krsz-zs)*krsz)
#pure decay noise(T = 0) 
alpha = 0 
beta = ka 
gamma = 0 
p1 = np.exp(-ka*t) 
E1 = np.array([[1,0],[0, np.sqrt(p1)]]) 
E2 = np.array([[0,np.sqrt(1-p1)],[0,0]]) 
E3 = np.array([[0,0],[0,0]]) 
E4 = np.array([[0,0],[0,0]]) Fdec.append((np.linalg.norm(E1)**2+np.linalg.norm(E2)**2+np.linalg.norm(E3)**2+np.linalg.norm(E4)**2+db)/db/(db+1)) 
print('Fidelity describes pure decay noise finite tamperature: ') 
print(Fdec) # pure dephasing noise(T -> infinite) 
p2 = np.sqrt(0.5*(1+np.exp(-ka*t))) 
alpha = 0 
beta = 0 
gamma = ka 
E1 = np.sqrt(p2)*ii 
E2 = np.sqrt(1-p2)*krsz 
E3 = [[0,0],[0,0]] 
E4 = [[0,0],[0,0]] Fdep.append(2/3+np.exp(-ka*t)/3) print('Fidelity describes pure dephasing finite tamperature: ') 
print(Fdep) 
# depolarizing noise(T -> infinite) 
alpha = 0.5*ka 
beta = 0.5*ka 
gamma = 0.5*ka 
p3 = 0.25*(1-np.exp(-ka*t)) 
E1 = np.sqrt(1-3*p3)*ii 
E2 = np.sqrt(p3)*krsx 
E3 = np.sqrt(p3)*krsy 
E4 = np.sqrt(p3)*krsz Fdepl.append(1/2+np.exp(-ka*t)/2) print('Fidelity describes depolarizing finite tamperature: ') 
print(Fdepl) 
# finite-T noise(alpha+beta = ka) 
alpha = 0.5*ka 
beta = 0.5*ka 
gamma = 0 
E1 = np.sqrt(Pu)*np.array([[1,0],[0,np.sqrt(np.exp(-ka*t))]]) E2 = np.sqrt(Pu)*np.array([[0,1-np.sqrt(np.exp(-ka*t))],[0,0]]) E3 = np.sqrt(Pd)*np.array([[np.sqrt(np.exp(-ka*t)),0],[0,1]]) E4 = np.sqrt(Pd)*np.array([[0,0],[1-np.sqrt(np.exp(-ka*t)),0]]) FT.append(0.5+np.exp(-ka*t/2)+np.exp(-ka*t)/6) print('Fidelity describes one special finite tamperature: ') print(FT) i = i+1 plt.figure() plt.plot(F-max(F)+0.99)plt.title('general fidelity v.s. N') plt.figure() plt.plot(F0-max(F0)+0.99)plt.title('general original fidelity v.s. N') plt.figure() plt.plot(Fdec)plt.title('Pure Decay noise fidelity v.s. N') plt.figure() plt.plot(Fdep)plt.title('Pure Dephasing noise fidelity v.s. N') plt.figure() plt.plot(Fdepl)plt.title('Pure Deplorization fidelity v.s. N') import pandas as pddf = np.array([F-max(F)+0.99,F0-max(F0)+0.99,Fdec, Fdep, Fdepl])#df.columns = ['F','F0','Fdecay','Fdephasing','Fdepolarizaing']#df.head()import seaborn as snsplt.figure() #ax = sns.boxplot(x='F', y='value', data=df, color='#99c2a2')plt.boxplot(np.transpose(df))plt.show() import scipy.stats as stats# stats f_oneway functions takes the groups as input and returns ANOVA F and p valuefvalue, pvalue = stats.f_oneway(df[0,:],df[1,:],df[2,:],df[3,:],df[4,:])print(fvalue, pvalue)# 4568774956086.299 0.0 # # get ANOVA table as R like output# import statsmodels.api as sm# from statsmodels.formula.api import ols# model = ols('C(1) ~ C(2):C(5)', data=df).fit()# from statsmodels.formula.api import ols# anova_table = sm.stats.anova_lm(model, typ=2)# anova_table
