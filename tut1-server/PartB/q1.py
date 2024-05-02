import numpy as np 
from matplotlib import pyplot as plt

def hermite(n, x):
    if n==0 : return np.ones((x.shape[0],))
    if n==1 : return 2*x 
    return 2*x*hermite(n-1, x) + 2*(n-1)*hermite(n-2,x) 

x = np.arange(-2,2.001, 0.001)
n = 3
y = hermite(n,x)
print(y)

plt.plot(x, y) 

plt.show()