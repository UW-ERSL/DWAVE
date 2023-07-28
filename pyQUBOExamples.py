#from dwave.system import DWaveSampler, EmbeddingComposite
from dimod.reference.samplers import ExactSolver, SimulatedAnnealingSampler
from pyqubo import Binary, Array

#%% Use pyQubo to solve a simple quadratic 
a,b,c = Binary("a"), Binary("b"),Binary("c")
H  = -0.1*a - 0.5*b -0.2*c +0.6*a*b +0.2*a*c + 0.03*b*c
model = H.compile()
Q = model.to_bqm()

sampler = SimulatedAnnealingSampler()
results = sampler.sample(Q, num_reads=100) 

print(results.first)
#%% Solving a 2x2 linear system using a total of 8 qubits via QUBO
A = [[3,1],[-1,2]]
b = [-1,12]

nDimensions = 2
q11,q12,q13,q14 = Binary("q11"), Binary("q12"),Binary("q13"),Binary("q14")
q21,q22,q23,q24 = Binary("q21"), Binary("q22"),Binary("q23"),Binary("q24")

x = nDimensions*[None]
x[0] = q11 + 2*q12 - q13 - 2*q14
x[1] = q21 + 2*q22 - q23 - 2*q24

A = [[3,1],[-1,2]]
b = [-1,5]
H = 0
for  i in range(nDimensions):
    Ax = sum(A[i][j]*x[j] for j in range(nDimensions))
    H = H + (Ax-b[i])**2

model = H.compile()
bqm = model.to_bqm()

sampler = ExactSolver()
results = sampler.sample(bqm) 

s = results.first.sample

xSol = nDimensions*[None]
xSol[0] = s["q11"] + 2*s["q12"] - s["q13"] - 2*s["q14"]
xSol[1] = s["q21"] + 2*s["q22"] - s["q23"] - 2*s["q24"]

print("Computed Solution:", xSol)

residual = 0
for  i in range(nDimensions):
    Ax = sum(A[i][j]*xSol[j] for j in range(nDimensions))
    residual = residual + (Ax-b[i])**2

print("Residual:",residual) #should be zero if answer is correct

#%% Random matrix of higher order
m = 3 #precision
print("Precision:",m)

A = [[2,-1,0],[-1,2,-1],[0,-1,2]]
b = [2,4,2]
nDimensions = len(A)
print("Total qubits:",2*nDimensions*m)
qPlus = nDimensions*[None]
qMinus = nDimensions*[None]
for i in range(nDimensions):
    qPlus[i] = Array.create("qPlus[" + str(i)+"]",shape = m,vartype = "BINARY")
    qMinus[i] = Array.create("qMinus["+ str(i)+"]",shape = m,vartype = "BINARY")

x = nDimensions*[None]
for i in range(nDimensions):
    x[i] = sum( ((2**k)*qPlus[i][k] - (2**k)*qMinus[i][k]) for k in range(m))
  
H = 0
for  i in range(nDimensions):
    Ax = sum(A[i][j]*x[j] for j in range(nDimensions))
    H = H + (Ax-b[i])**2

model = H.compile()
bqm = model.to_bqm()

if (0):
    sampler = SimulatedAnnealingSampler()
    results = sampler.sample(bqm, num_reads=100) 
else:
    sampler = ExactSolver()
    results = sampler.sample(bqm) 

xSol = nDimensions*[None]
for i in range(nDimensions):
    xSol[i] = sum(((2**k)*results.first.sample["qPlus[" + str(i)+"][" + str(k) + "]"]- \
        (2**k)*results.first.sample["qMinus[" + str(i)+"][" + str(k) + "]"]) for k in range(m))
        
print("Computed Solution:", xSol)
residual = 0
for  i in range(nDimensions):
    Ax = sum(A[i][j]*xSol[j] for j in range(nDimensions))
    residual = residual + (Ax-b[i])**2

print("Residual:",residual)  
    