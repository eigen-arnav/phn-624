def factorial(n): 
    ans = 1
    for i in range(1, n+1):
        ans*= i 
    return ans 

def double_factorial(n):
    ans=1
    for i in range(n, 0, -2):
        ans*=i
    return ans

print(factorial(7))
print(double_factorial(7))