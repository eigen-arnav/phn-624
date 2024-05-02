import numpy as np 

def simpson(f, a, b, n):
    if n%2==1 : raise Exception("n must be even")
    h = (b-a)/n
    int = f(a)+f(b)
    for i in  range(1, n, 2):
        int+=4*f(a+h*i)
    for i in range(2, n-1, 2):
        int+=2*f(a+h*i)
    return h*int/3

A = 238
Z = 92 
r1 = 1.07*A**(1/3)*1e-15
r2 = 8.0*A**(1/3)*1e-15
hbar = 6.626*1e-34/(2*np.pi)
m = 6.644*1e-27 
E = 5*1.6*1e-13
V_coulomb = lambda r : 9*1e9*2*Z*(1.6*1e-19)**2/r # input in m and output in J
f_integral = lambda r : np.sqrt(2*m*(V_coulomb(r) - E))/hbar
integral = simpson(f_integral, r1, r2, 10000)
T = np.exp(-2*integral) 
print("Transition Probability : ", T)