import numpy as np 
from matplotlib import pyplot as plt 
from math import factorial as fact 

F1 = lambda j1, j2, j3, m1, m2, m3 : np.sqrt(fact(j1+j2-j3)*fact(j3+j1-j2)*fact(j2+j3-j1)*(2*j3+1)/fact(j1+j2+j3+1))
F2 = lambda j1, j2, j3, m1, m2, m3 : np.sqrt(fact(j3+m3)*fact(j3-m3)*fact(j2+m2)*fact(j2-m2)*fact(j1+m1)*fact(j1-m1)*1.0)
def F3(j1, j2, j3, m1, m2, m3):
    smax = min(j1-m1, j2+m2, j1+j2-j3)
    smin = abs(min(min(j3-j2+m1, j3-j1-m2), 0))
    ans = 0
    for s in np.arange(smin, smax+1):
        ans += (-1)**s/(fact(j1-m1-s)*fact(j2+m2-s)*fact(j3-j2+m1+s)*fact(j3-j1-m2+s)*fact(j1+j2-j3-s)*fact(s))
    return ans

CG = lambda j1, j2, j3, m1, m2, m3 : F1(j1, j2, j3, m1, m2, m3)*F2(j1, j2, j3, m1, m2, m3)*F3(j1, j2, j3, m1, m2, m3)

BE2 = lambda I, Qt : CG(I, 2, I-2, 2.5, 0, 2.5)**2*5*Qt**2/(16*np.pi)

print(BE2(8.5, 224))
print(BE2(10.5, 202))