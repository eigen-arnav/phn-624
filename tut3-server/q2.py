import numpy as np 
from matplotlib import pyplot as plt 
from tqdm import tqdm 

def runge_kutta(f, t0, tf, v0, h=0.01):
    if type(v0)==list : v0 = np.array(v0)
    N = int((tf-t0)/h)
    v = np.zeros((N, v0.shape[0]))
    v[0, :] = v0 
    for i in range(N-1): 
        k1 = f(t0 + i*h, v[i, :])
        k2 = f(t0 + (i+1/2)*h, v[i, :] + h*k1/2)
        k3 = f(t0 + (i+1/2)*h, v[i, :] + h*k2/2)
        k4 = f(t0 + (i+1)*h, v[i, :] + h*k3)
        v[i+1, :] = v[i, :] + h*(k1+2*(k2+k3)+k4)/6
    return v 

def zero_bs(f, lb, ub, precision = 1e-4):
    low = lb ; high = ub ; mid = .5*(low+high) ; mid_prev=0
    while abs(f(mid)-f(mid_prev))>precision:
        mid_prev = mid
        if f(low)*f(mid)>0 : low = mid 
        else : high = mid 
        mid = .5*(low+high)
    return mid

def root_finder(f, lb, ub, precision=1e-4, delta=0.1):
    roots = []
    start = lb ; stop = lb 
    print("Finding Roots")
    while(stop<=ub): 
        stop += delta
        if(f(start)*f(stop)<0):
            roots.append(zero_bs(f, start, stop, precision))
            start = stop 
        print(str((start-lb)/(ub-lb)*100) + " percent done")
    return roots

m = 938.272 # MeV/c**2 
hbarc  = 197.326931 # MeV-fm 
# m = 1 ; hbarc=1

class inf_well:
    def __init__(self, a=10, h=0.01):
        self.Vo = np.inf
        self.E = 0
        self.a = a
        self.h = h
        self.mesh = np.arange(-self.a, self.a, h)
    def V(self, x):
        return 0 if x>=-self.a and x<self.a else self.Vo
    def set_E(self, E): 
        self.E=E
    def rk_return(self, x, v):
        return np.array([v[1], (2*m/hbarc**2)*(self.V(x) - self.E)*v[0]])
    def energy_state(self, E):
        self.set_E(E) 
        return runge_kutta(self.rk_return, -self.a, self.a, [0,0.001], self.h)[:, 0]
    def wf_end(self, E):
        return self.energy_state(E)[-1]
    def eigen_energies(self, Estart, Eend, precision=1e-4, delta=0.1): 
        return root_finder(self.wf_end, Estart, Eend, precision, delta)
    def eigen_info(self, Estart, Eend, precision=1e-4, delta=0.1):
        Es = self.eigen_energies(Estart, Eend, precision, delta)
        wfs = np.zeros((self.mesh.shape[0], len(Es)))
        for i, E in enumerate(Es) : 
            wfs[:, i] = self.energy_state(E)
        return Es, wfs   
    
well_obj = inf_well() 
En, Psi = well_obj.eigen_info(0, 20) 
legend_info = [f"n = {i} ; E = {round(En[i], 3)}" for i in range(len(En))]
plt.plot(well_obj.mesh, Psi)
plt.legend(legend_info)
plt.show()


