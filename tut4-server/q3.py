import numpy as np 
from math import factorial

fact = lambda x : factorial(int(x))

F1 = lambda j1, j2, j3, m1, m2, m3 : np.sqrt(fact(j1+j2-j3)*fact(j3+j1-j2)*fact(j2+j3-j1)*(2*j3+1)/fact(j1+j2+j3+1))
F2 = lambda j1, j2, j3, m1, m2, m3 : np.sqrt(fact(j3+m3)*fact(j3-m3)*fact(j2+m2)*fact(j2-m2)*fact(j1+m1)*fact(j1-m1))
def F3(j1, j2, j3, m1, m2, m3):
    smax = min(j1-m1, j2+m2, j1+j2-j3)
    smin = abs(min(min(j3-j2+m1, j3-j1-m2), 0))
    ans = 0
    for s in np.arange(smin, smax+1):
        ans += (-1)**s/(fact(j1-m1-s)*fact(j2+m2-s)*fact(j3-j2+m1+s)*fact(j3-j1-m2+s)*fact(j1+j2-j3-s)*fact(s))
    return ans

CG = lambda j1, j2, j3, m1, m2, m3 : int(not (m1+m2-m3))*F1(j1, j2, j3, m1, m2, m3)*F2(j1, j2, j3, m1, m2, m3)*F3(j1, j2, j3, m1, m2, m3) if (abs(m1)<=j1 and abs(m2)<=j2 and abs(m3)<=j3) else 0

def tbme(oa, ob, oc, od, J, T):
    At=1
    na, nb, nc, nd = oa[0], ob[0], oc[0], od[0]
    la, lb, lc, ld = oa[1], ob[1], oc[1], od[1]
    ja, jb, jc, jd = oa[2], ob[2], oc[2], od[2]
    if J not in np.intersect1d(np.arange(abs(ja-jb), ja+jb), np.arange(abs(jc-jd), jc+jd)) : raise Exception("Check J") 
    delta_ab = 1 if oa==ob else 0 
    delta_cd = 1 if oc==od else 0 
    if (delta_ab or delta_cd) and (J+T)%2==0 : raise Exception("Check T")
    term1 = (-1)**(na+nb+nc+nd) * At / (2*(2*J+1))
    term2 = np.sqrt((2*ja+1)*(2*jb+1)*(2*jc+1)*(2*jd+1) / ((1+delta_ab)*(1+delta_cd)) )
    term3 = (-1)**(jb+jd+lb+ld) * CG(jb, ja, J, -0.5, 0.5, 0) * CG(jd, jc, J, -0.5, 0.5, 0) * (1 - (-1)**(la+lb+J+T))
    term4 = CG(jb, ja, J, 0.5, 0.5, 1) * CG(jd, jc, J, 0.5, 0.5, 1) * (1 + (-1)**T)
    return term1 * term2 * (term3 - term4)

s1 = [0, 0, 1/2]
s2 = [0, 1, 3/2] 
s3 = [0, 1, 1/2]
print(tbme(s1, s1, s2, s2, 0, 1))