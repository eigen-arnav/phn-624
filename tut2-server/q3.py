import numpy as np 
from matplotlib import pyplot as plt 

av = 15.75 ; aS = 17.80 ; ac = 0.71 ; aA = 23.70 ; ap = 33.50
En = lambda A, Z : (av*A - aS*A**(2/3) - aA*(A-2*Z)**2/A - ac*Z**2/A**(1/3) + ap* (1-A%2)*(-1)**Z/(A**(3/4)))/A
Ens = []
for A in range(1, 251): 
    Ens.append(En(A, A//2))
plt.title("Binding Energy Per Nucleon")
plt.xlabel("Mass No. - A")
plt.ylabel("En (MeV)")
plt.plot(range(1, 251), Ens, '-+')
plt.show()