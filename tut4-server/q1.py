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