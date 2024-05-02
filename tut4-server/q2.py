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
