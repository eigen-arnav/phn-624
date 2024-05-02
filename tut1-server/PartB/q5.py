import numpy as np 
from matplotlib import pyplot as plt
from scipy import special

def ass_leg(l, m, x): 
    if l<0 or abs(m)>l : raise Exception("Error!!! Please check the values of l=" + l + "and m=" + m)
    if m==l : return (-1)**l*(1-x**2)**(l/2)*special.factorial2(2*l-1) 
    if m==l-1 : return x*(2*l-1)*ass_leg(l-1, l-1, x) 
    if m>=0 : return ((2*l-1)*x*ass_leg(l-1, m, x) - (l+m-1)*ass_leg(l-2, m, x))/(l+m)
    if m<0 : return (-1)**m*np.math.factorial(l+m)*ass_leg(l, -m, x)/np.math.factorial(l-m) 
def spherical_harmonics(l, m, theta, phi):
    Plm = ass_leg(l,m,np.cos(theta))
    phi_part = np.exp(1j*m*phi)
    mat = np.outer(Plm, phi_part)
    Ylm = mat*(-1)**m*np.sqrt(((2*l+1)*(np.math.factorial(l-m)))/(4*np.pi*np.math.factorial(l+m)))
    return Ylm

theta = np.linspace(0, np.pi, 1000) 
phi = np.linspace(-np.pi, np.pi, 2000)

fig = plt.figure()
fig.suptitle("Real and Imaginary part of Spherical Harmonics")

T, P = np.meshgrid(theta, phi, indexing='ij')
R = spherical_harmonics(3, 1, theta, phi)
R_r, R_i = np.real(R), np.imag(R)

ax1 = fig.add_subplot(121, projection='3d')
ax1.set_title("Real Part")
X = np.sin(T)*np.cos(P)*np.abs(R_r)
Y = np.sin(T)*np.sin(P)*np.abs(R_r)
Z_r = np.cos(T)*np.abs(R_r)
ax1.plot_surface(X, Y, Z_r)
ax1.set_xlabel("x")
ax1.set_ylabel("y")
ax1.set_zlabel("Real(Ylm)")

ax2 = fig.add_subplot(122, projection='3d')
ax2.set_title("Imaginary Part")
X = np.sin(T)*np.cos(P)*np.abs(R_i)
Y = np.sin(T)*np.sin(P)*np.abs(R_i)
Z_i = np.cos(T)*np.abs(R_i)
ax2.plot_surface(X, Y, Z_i)
ax2.set_xlabel("x")
ax2.set_ylabel("y")
ax2.set_zlabel("Imag(Ylm)")

plt.show()
