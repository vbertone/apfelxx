import numpy as np
from scipy.special import gamma, factorial, gegenbauer
import math

# Anomalous dimension
def V(n):
    return 3./2. - 1. / ( n + 1 ) - 1. / ( n + 2 ) - 2 * sum([1. / k for k in range(1, n + 1)])

# Gegenbauer coefficients
def a(n, k, alpha):
    l = n - 2 * k
    return (-1)**k * 2**l * gamma(l + k + alpha) / gamma(alpha) / factorial(k) / factorial(l)

# Expected value of the integral
def Iexp(n, z):
    return V(n) * gegenbauer(n, 3/2)(z)

# Computed integral
def Icom(n, z):
    summation = 0
    for k in range(math.floor(n/2)+1):
        l = n - 2 * k
        summation += a(n, k, 3/2) * ( sum([ ( 1. / ( l + 2 ) + 1. / ( l + 1 ) - 2. / ( l - j ) ) * ( 1 + (-1)**(l-j) ) * z**j / 2. for j in range(l)]) + ( 1. / ( l + 2 ) + 1. / ( l + 1 ) + 2 * sum([ 1. / ( j + 1 ) for j in range(l)]) ) * z**l )
    return 3 * gegenbauer(n, 3/2)(z) / 2 - summation

# Computed integral
def Icom2(n, z):
    coeffs = []
    for h in range(math.floor(n/2)+1):
        terms = [3. / 2, - 1. / ( n - 2 * h + 2 ), - 1. / ( n - 2 * h + 1 )]
        terms.extend([ - 2. / j for j in range(1,n-2*h+1)])
        terms.extend([ - a(n, h-j, 3/2) / a(n, h, 3/2) * ( 1. / ( n - 2 * h + 2 * j + 2 ) + 1. / ( n - 2 * h + 2 * j + 1 ) - 1. / j ) for j in range(1,h+1)])
        coeffs.append(terms)
        #print(h, terms)
    for h in range(math.floor(n/2)+1):
        print(sum(coeffs[h]), V(n))
    return sum([a(n, h, 3/2) * z**(n-2*h) * sum(coeffs[h]) for h in range(math.floor(n/2)+1)])

# Check
n = 12
z = 0.7
print(Iexp(n, z), Icom(n, z), Icom2(n, z))
