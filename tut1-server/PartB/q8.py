import numpy as np 
from matplotlib import pyplot as plt 
from math import factorial as fact 

F1 = lambda j1, j2, j3, m1, m2, m3 : np.sqrt(fact(j1+j2-j3)*fact(j3+j1-j2)*fact(j2+j3-j1)*(2*j3+1)/fact(j1+j2+j3+1))
F2 = lambda j1, j2, j3, m1, m2, m3 : np.sqrt(fact(j3+m3)*fact(j3-m3)*fact(j2+m2)*fact(j2-m2)*fact(j1+m1)*fact(j1-m1))
def F3(j1, j2, j3, m1, m2, m3):
    smax = min(j1-m1, j2+m2, j1+j2-j3)
    smin = abs(min(min(j3-j2+m1, j3-j1-m2), 0))
    ans = 0
    for s in np.arange(smin, smax+1):
        ans += (-1)**s/(fact(j1-m1-s)*fact(j2+m2-s)*fact(j3-j2+m1+s)*fact(j3-j1-m2+s)*fact(j1+j2-j3-s)*fact(s))
    return ans

CG = lambda j1, j2, j3, m1, m2, m3 : int(not (m1+m2-m3))*F1(j1, j2, j3, m1, m2, m3)*F2(j1, j2, j3, m1, m2, m3)*F3(j1, j2, j3, m1, m2, m3)

j1 = -.5 ; j2 = .5 
m1 = np.arange(-j1, j1+1) ; m2 = np.arange(-j2, j2+1)
j3 = np.arange(abs(j1-j2), j1+j2+1)

for M1 in m1 : 
    for M2 in m2 :
        for J in j3:
            for M in np.arange(-J, J+1) :
                cg = CG(j1, j2, J, M1, M2, M)
                print("C.G. Coefficient for j1 = ", j1, " j2 = ", j2,  " m1 = ", M1, " m2 = ", M2, " J = ", J, " M = ", M, " : ", cg)
