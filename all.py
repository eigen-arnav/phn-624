import numpy as np 
from matplotlib import pyplot as plt 

C = -0.15 ; D = -0.0225
orig_energy = lambda N : (N+3/2) 
split_energy = lambda N, l, j : orig_energy(N) + D*l*(l+1) + (C*l/2 if j== l+1/2 else -C/2*(l+1))

for N in range(5):
    E0 = orig_energy(N)
    plt.axhline(E0, 0, 0.4, color='blue')
    for l in range(N, -1, -2):
        for j in [abs(l-1/2), l+1/2]:
            E = split_energy(N, l, j) 
            plt.plot([0.4, 0.6], [E0, E], color='violet', linestyle="dashed")
            plt.axhline(E, 0.6, 1.0, color='red')
plt.xlim([0,1])
plt.show()

""" Question 2 """


import numpy as np 
from matplotlib import pyplot as plt 

E = lambda N, nz, delosc : N + 1.5 -delosc/3*(3*nz - N)

color_list = ['red', 'blue', 'green', 'orange', 'violet']
delosc_arr = np.linspace(-1,1,100)
for N in range(5):
    for nz in range(N+1):
        energies = E(N, nz, delosc_arr)
        plt.plot(delosc_arr, energies, color=color_list[N])
plt.show()

A = 80
delosc_arr = np.linspace(-1,1,100)
hw0 = 41*A**(-1/3)*((1+(2/3)*delosc_arr)**2*(1-(4/3)*delosc_arr))**(-1/6)
color_list = ['red', 'blue', 'green', 'orange', 'violet']
for N in range(5):
    for nz in range(N+1):
        energies = hw0*E(N, nz, delosc_arr)
        plt.plot(delosc_arr, energies, color=color_list[N])
plt.show()

""" Question 3 """

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

""" Question 4 """

from math import factorial
import numpy as np 
from matplotlib import pyplot as plt 

fact = lambda x : factorial(int(x))
sqrt = lambda x : np.sqrt(abs(x))*np.sign(x)

F1 = lambda j1, j2, j3, m1, m2, m3 : sqrt(fact(j1+j2-j3)*fact(j3+j1-j2)*fact(j2+j3-j1)*(2*j3+1)/fact(j1+j2+j3+1))
F2 = lambda j1, j2, j3, m1, m2, m3 : sqrt(fact(j3+m3)*fact(j3-m3)*fact(j2+m2)*fact(j2-m2)*fact(j1+m1)*fact(j1-m1))
def F3(j1, j2, j3, m1, m2, m3):
    smax = min(j1-m1, j2+m2, j1+j2-j3)
    smin = abs(min(min(j3-j2+m1, j3-j1-m2), 0))
    ans = 0
    for s in np.arange(smin, smax+1):
        ans += (-1)**s/(fact(j1-m1-s)*fact(j2+m2-s)*fact(j3-j2+m1+s)*fact(j3-j1-m2+s)*fact(j1+j2-j3-s)*fact(s))
    return ans

def CG(j1, j2, j3, m1, m2, m3) : 
    if not m1+m2==m3 : return 0 
    if abs(m1) > j1 or abs(m2) > j2 or abs(m3) > j3 : return 0 
    if j3 > j1+j2 or j3< abs(j1-j2) : return 0 
    else : return F1(j1, j2, j3, m1, m2, m3)*F2(j1, j2, j3, m1, m2, m3)*F3(j1, j2, j3, m1, m2, m3)

delta = lambda x, y : int(x==y)

A = 80
kappa = 0.05 
mu = [0, 0, 0, 0.35, 0.625, 0.63, 0.448, 0.434]
fdelta = lambda delta: ((1+(2/3)*delta)**2*(1-(4/3)*delta))**(-1/6) 
hw0 = 1
hw00_func = lambda delta : hw0/fdelta(delta)
C_func = lambda hw : -2*kappa*hw
D_func = lambda C : [C*m/2 for m in mu]

#hw00 = 41*A**(-1/3)
#hw0 = lambda delta: hw00*fdelta(delta) 
#C = -2*kappa*hw00 
#D = [C*m/2 for m in mu]

def shell_basis(N, omega):
    if omega-0.5>N : raise Exception("Omega exceeds N+0.5")
    basis_set = [] 
    for l in range(N, -1, -2):
        for lam in range(-l, l+1):
            sigma = omega - lam 
            if sigma==.5 or sigma ==-.5 : 
                basis_set.append([l, lam, sigma]) 
    return basis_set

def nilsson_hamiltonian(N, omega, deltaf):
    hw00 = hw00_func(deltaf) 
    C = C_func(hw00) 
    D = D_func(C)
    basis = shell_basis(N, omega) 
    dim = len(basis)
    H = np.zeros((dim, dim))
    for i, a in enumerate(basis): 
        for j, b in enumerate(basis): 
            la, lb = a[0], b[0] 
            lama, lamb = a[1], b[1] 
            siga, sigb = a[2], b[2]
            if i==j: H[i, j] += (hw0*(N+1.5) + D[N]*la*(la+1))
            H[i, j]+= C*(.5*sqrt((la-lamb)*(la+lamb+1))*delta(lama, lamb+1)*delta(siga, sigb-1)+\
                .5*sqrt((la+lamb)*(la-lamb+1))*delta(lama, lamb-1)*delta(siga, sigb+1)+\
                    lama*siga*delta(lama, lamb)*delta(siga, sigb))*delta(la, lb)
            ex_r2 = sqrt((N-lb+2)*(N+lb+1))*delta(la, lb-2) + sqrt((N-lb)*(N+lb+3))*delta(la, lb+2) + (N+1.5)*delta(la, lb)
            ex_Y = sqrt((2*lb+1)/(2*la+1))*CG(lb, 2, la, lamb, 0, lama)*CG(lb, 2, la, 0, 0, 0) 
            H[i, j] += -deltaf*hw0*(2/3)*ex_r2*ex_Y*delta(lama, lamb)*delta(siga,sigb) 
    return H 

def cranked_hamiltonian(N, omega, freq):
    deltaf = 0
    hw00 = hw00_func(deltaf) 
    C = C_func(hw00) 
    D = D_func(C)
    basis = shell_basis(N, omega) 
    dim = len(basis)
    H = np.zeros((dim, dim))
    for i, a in enumerate(basis): 
        for j, b in enumerate(basis): 
            la, lb = a[0], b[0] 
            lama, lamb = a[1], b[1] 
            siga, sigb = a[2], b[2]
            if i==j: H[i, j] += hw0*(N+1.5) + D[N]*la*(la+1) - hw0*freq*omega
            H[i, j]+= C*(.5*sqrt((la-lamb)*(la+lamb+1))*delta(lama, lamb+1)*delta(siga, sigb-1)+\
                .5*sqrt((la+lamb)*(la-lamb+1))*delta(lama, lamb-1)*delta(siga, sigb+1)+\
                    lama*siga*delta(lama, lamb)*delta(siga, sigb))*delta(la, lb)
            ex_r2 = sqrt((N-lb+2)*(N+lb+1))*delta(la, lb-2) + sqrt((N-lb)*(N+lb+3))*delta(la, lb+2) + (N+1.5)*delta(la, lb)
            ex_Y = sqrt((2*lb+1)/(2*la+1))*CG(lb, 2, la, lamb, 0, lama)*CG(lb, 2, la, 0, 0, 0) 
            H[i, j] += -deltaf*hw0*(2/3)*ex_r2*ex_Y*delta(lama, lamb)*delta(siga,sigb) 
    return H 

def plot_energy(N, omega, deltafs, type="nilsson"):

    Es = np.zeros((deltafs.shape[0], len(shell_basis(N, omega)))) 
    for i, deltaf in enumerate(deltafs):
        if type=="cranked" : H = cranked_hamiltonian(N, omega, deltaf)
        else : H = nilsson_hamiltonian(N, omega, deltaf)
        E = np.linalg.eigvalsh(H) 
        E = sorted(E)
        Es[i, :] = E 
    plt.plot(deltafs, Es) 

deltafs = np.linspace(-0.3, 0.3, 100)
for N in range(3) : 
    for o in range(N+1) : 
        omega = o + .5 
        plot_energy(N, omega, deltafs)
plt.show()

deltafs = np.linspace(0, 0.2, 100)
for N in range(3) : 
    for o in range(N+1) : 
        omega = o + .5 
        plot_energy(N, omega, deltafs, "cranked")
plt.show()

""" Question 5 """

N_max = 7
Np = 40 
temps = np.arange(0.5, 5.5, .5) 

def energies(N_max, deltaf=0.2):
    Es = []
    for N in range(N_max+1): 
        for o in range(N+1):
            omega = o+0.5 
            H = hamiltonian(N, omega, deltaf) 
            E = np.linalg.eigvalsh(H)
            Es.extend(E) 
    return np.array(sorted(Es))

fd = lambda E, lam, kT : 1/(1+np.exp((E-lam)/kT))
def conservation(lam, Es, kT, Np): 
    return Np/2 - np.sum(fd(Es, lam, kT))
Ef = lambda Es, kT, Np : newton(conservation, x0 = 4, args=(Es, kT, Np))

for kT in temps : 
    Es = energies(N_max)
    lam = Ef(Es, kT, Np) 
    ni = fd(Es, lam, kT) 
    plt.plot(Es, ni)
plt.show()

for kT in temps : 
    Es = energies(N_max)
    lam = Ef(Es, kT, Np) 
    ni = fd(Es, lam, kT) 
    si = -ni*np.log(ni) - (1-ni)*np.log(1-ni)
    plt.plot(Es, si)
plt.show()

F = []
Ts = np.linspace(1,10,100)
for T in Ts: 
    Es = energies(N_max)
    lam = Ef(Es, kT, Np) 
    ni = fd(Es, lam, kT) 
    si = -ni*np.log(ni) - (1-ni)*np.log(1-ni)
    S = np.sum(si) 
    ET = np.sum(Es*ni)
    F.append(ET-T*S)
plt.plot(Ts, F)
plt.show()


""" Question 6"""