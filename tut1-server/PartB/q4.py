import numpy as np 
from matplotlib import pyplot as plt
from scipy import special

def legendre(n, x):
    if n==0 : return np.ones((x.shape[0],))
    if n==1 : return x 
    return (2-1/n)*x*legendre(n-1, x) - (1-1/n)*legendre(n-2,x) 
def ass_leg(l, m, x): 
    if l<0 or abs(m)>l : raise Exception("Error!!! Please check the values of l=" + l + "and m=" + m)
    if m==l : return (-1)**l*(1-x**2)**(l/2)*special.factorial2(2*l-1) 
    if m==l-1 : return x*(2*l-1)*ass_leg(l-1, l-1, x) 
    if m>=0 : return ((2*l-1)*x*ass_leg(l-1, m, x) - (l+m-1)*ass_leg(l-2, m, x))/(l+m)
    if m<0 : return (-1)**m*np.math.factorial(l+m)*ass_leg(l, -m, x)/np.math.factorial(l-m) 

x = np.arange(-1,1.001, 0.001)

legend_list = [] # List for defining legends in the plot
for n in range(5):
    y = legendre(n,x)
    legend_list.append("n = "+str(n))
    plt.plot(x,y)
plt.xlabel("x")
plt.ylabel("legendre_n(x)")
plt.title("Legendre Polynomials")
plt.legend(legend_list)
plt.grid()
plt.show()

legend_list = [] # List for defining legends in the plot
for l in range(1, 3) : 
  for m in range(-l, l+1):
      y = ass_leg(l, m, x)
      legend_list.append("l = "+str(l)+", m = "+str(m))
      plt.plot(x,y)
plt.ylim(-5, 10)
plt.xlabel("x")
plt.ylabel("assosciated_legendre_n(x)")
plt.title("Assosciated Legendre Polynomials")
plt.legend(legend_list)
plt.grid()
plt.show()
