import numpy as np 
from matplotlib import pyplot as plt 

def laguerre(k, x):
    if k==0 : return np.ones((x.shape[0], ))
    if k==1 : return 1 - x 
    if k>1 : return ((2*k - 1 - x)*laguerre(k-1, x) - (k-1)*laguerre(k-2, x))/k
    else : return np.exp(x)*laguerre(-k-1, -x)

def ass_lag(n, l, x) : 
    if n < 0 or l < 0 : raise Exception("!!!Error!!! Please recheck values of n and l ")
    if n==0 : return  np.ones((x.shape[0], ))
    if n==1 : return 1 + l - x 
    else : return ((2*n - 1 + l - x)*ass_lag(n-1, l, x) - (n-1+l)*ass_lag(n-2, l, x))/n

x = np.arange(-2, 10, 0.005)

legend_list = [] # List for defining legends in the plot
for n in range(5):
    y = laguerre(n,x)
    legend_list.append("n = "+str(n))
    plt.plot(x,y)
plt.xlabel("x")
plt.ylabel("laguerre_n(x)")
plt.title("Laguerre Polynomials")
plt.legend(legend_list)
plt.grid()
plt.show()

legend_list = [] # List for defining legends in the plot
for n in range(4):
    for l in range(n+1):
        y = ass_lag(n,l,x)
        legend_list.append("n = "+str(n)+", l = "+str(l))
        plt.plot(x,y)
plt.ylim(-5, 10)
plt.xlabel("x")
plt.ylabel("assosciated_laguerre_n(x)")
plt.title("Assosciated Laguerre Polynomials")
plt.legend(legend_list)
plt.grid()
plt.show()