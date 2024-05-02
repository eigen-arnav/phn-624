import numpy as np 
from matplotlib import pyplot as plt

def hermite(n, x):
    if n==0 : return np.ones((x.shape[0],))
    if n==1 : return 2*x 
    return 2*x*hermite(n-1, x) + 2*(n-1)*hermite(n-2,x) 
    
def qho(n, x):
    return (1/(np.sqrt(2**n*np.math.factorial(n))))*np.pi**(-1/4)*np.exp(-x**2/2)*hermite(n, x)

x = np.arange(-2,2.001, 0.001)
n = 3
y = qho(n,x)
print(y)

plt.plot(x, y) 

plt.show()