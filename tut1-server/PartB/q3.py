import numpy as np 

def hermite(n, x):
    if n==0 : return 1
    if n==1 : return 2*x 
    return 2*x*hermite(n-1, x) - 2*(n-1)*hermite(n-2,x) 
    
def qho(n, x):
    return (1/(np.sqrt(2**n*np.math.factorial(n)*np.sqrt(np.pi))))*np.exp(-x**2/2)*hermite(n, x) 

def simpson(f, a, b, n):
    if n%2!=0 : raise Exception("Please enter an even number of integral points")
    h = (b-a)/n
    ans = f(a) + f(b) 
    for i in range(1, n):
        if i%2==0 : ans+=2*f(a+h*i)
        else : ans+=4*f(a+h*i)
    ans*=h/3
    return ans 

def phidashdash(n, x):
    if n==0 : return (x**2-1)*qho(0, x) 
    if n==1 : return -(8)**(1/2)*x*qho(0, x) + (x**2-1)*qho(1,x)
    return 2*(n**2-n)**(1/2)*qho(n-2, x) -(8*n)**(1/2)*x*qho(n-1, x) + (x**2-1)*qho(n, x) 

def KE_element(m, n): 
    f = lambda x : qho(m, x)*phidashdash(n, x)
    return -.5*simpson(f, -1000, 1000, 10000)
def PE_element(m, n): 
    f = lambda x : qho(m,x)*qho(n,x)*x**2
    return .5*simpson(f, -1000, 1000, 10000)

f = lambda x : qho(2, x)*qho(2, x)
print(simpson(f, -1000, 1000, 10000))

N = 4
KE = np.zeros((N, N))
PE = np.zeros((N, N))
for m in range(N):
    for n in range(N):
        KE[m ,n] = KE_element(m, n)
        PE[m ,n] = PE_element(m, n)
print(KE); print(PE); print(KE+PE)