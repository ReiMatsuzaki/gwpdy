#from sympy import *
from scipy.integrate import quad
from numpy import sqrt, exp, pi

#x = Symbol('x')

gA = 0.5
rA = 0.0
pA = 0.0
gB = gA
rB = 0.3
pB = 0.4

nA = (2*gA/pi)**0.25
nB = (2*gB/pi)**0.25

GA = lambda x: nA*exp(-gA*(x-rA)**2 + 1j*pA*(x-rA))
GB = lambda x: nB*exp(-gB*(x-rB)**2 + 1j*pB*(x-rB))
real_GG = lambda x: (GA(x).conjugate()*GB(x)).real
imag_GG = lambda x: (GA(x).conjugate()*GB(x)).imag
x0 = 10
print(quad(lambda x: abs(GA(x))**2, -x0, x0))
print(quad(real_GG, -x0, x0, epsabs=10.0**(-15)))
print(quad(imag_GG, -x0, x0, epsabs=10.0**(-15)))
