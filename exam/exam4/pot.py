import numpy as np
from numpy import exp, sqrt
import pandas as pd

from naewdy2.math import mat2csv

A = 0.01
B = 1.6
C = 0.005
D = 1.0

xs = np.linspace(-10, 10, 200)
nx = len(xs)
dx = xs[1]-xs[0]

v11 = np.where(xs>0, A*(1-exp(-B*xs)), -A*(1-exp(B*xs)))
v12 = C*exp(-D*xs*xs)
v21 = v12
v22 = -v11

de = sqrt((v11+v22)**2 - 4*(v11*v22-v12*v21))/2
e1 = (v22+v11)/2 - de
e2 = (v22+v11)/2 + de

us = [[v12,    v21],
      [e1-v11, e2-v11]]
n0 = 1/sqrt(us[0][0]**2 + us[1][0]**2)
n1 = 1/sqrt(us[0][1]**2 + us[1][1]**2)
uus = np.array([[ n0*us[0][0], n1*us[0][1] ],
                [ n0*us[1][0], n1*us[1][1]] ])
dx_uus = (np.roll(uus, -1, axis=2) - np.roll(uus, +1, axis=2)) / (2*dx)

i = 0
j = 1
x12 = uus[0,i,:]*dx_uus[0,j,:] + uus[1,i,:]*dx_uus[1,j,:]

HeIJ = np.zeros((nx, 2, 2))
XkIJ = np.zeros((nx, 1, 2, 2))
HeIJ[:,0,0] = e1
HeIJ[:,1,1] = e2
XkIJ[:,0,0,1] = x12
XkIJ[:,0,1,0] = -x12

pd.DataFrame({"val":xs}).to_csv("xs.csv", index=None)
mat2csv(HeIJ, "heij.csv")
mat2csv(XkIJ, "xkij.csv")


