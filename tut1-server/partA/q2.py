def sum_square(n):
    ans = 0
    for i in range(1, n+1):
        ans+=i**2
    return ans 

print(sum_square(5))
