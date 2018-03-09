import numpy as np
from numpy import exp
import pandas as pd

from naewdy2.math import mat2csv

A = 0.01
B = 1.6
C = 0.005
D = 1.0

xs = np.linspace(-10, 10, 200)
nx = len(xs)

HeIJ = np.zeros((nx, 2, 2))
XkIJ = np.zeros((nx, 1, 2, 2))

HeIJ[:,0,0] = np.where(xs>0, A*(1-exp(-B*xs)), -A*(1-exp(B*xs)))
HeIJ[:,0,1] = C*exp(-D*xs*xs)
HeIJ[:,1,0] = +HeIJ[:,0,1]
HeIJ[:,1,1] = -HeIJ[:,0,0]

pd.DataFrame({"val":xs}).to_csv("xs.csv", index=None)
mat2csv(HeIJ, "heij.csv")
mat2csv(XkIJ, "xkij.csv")


