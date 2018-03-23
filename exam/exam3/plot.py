import os
join = os.path.join
exists = os.path.exists
expanduser = os.path.expanduser
import sys
import matplotlib as mpl
mpl.use('PDF')
import matplotlib.pylab as plt
import numpy as np
tr = np.transpose
import pandas as pd
import json

from naewdy2.math import csv2mat
from gwpdy.dydata import correct_csvs


j = json.load(open("out/common.json", "r"))
nt = j["nt"]
dt = j["dt"]
ts = np.arange(nt)*dt

r1s = []
prob1 = []
prob2 = []
for it in range(nt):
    dir_it = join("out", str(it))
    if(not exists(dir_it)):
        break
    r = csv2mat(join(dir_it, "r.csv"))
    r1s.append(r[0,0])
    c = csv2mat(join(dir_it, "c.csv"))
    prob1.append(abs(c[0,0])**2)
    prob2.append(abs(c[0,1])**2)

n = len(r1s)    
plt.plot(ts[:n], r1s[:])
if(not exists("fig")):
    os.makedirs("fig")
plt.savefig("fig/r.pdf")
plt.close()

prob1 = np.array(prob1)
prob2 = np.array(prob2)
plt.plot(ts, prob1, label="1")
plt.plot(ts, prob2, label="2")
plt.plot(ts, prob1+prob2, label="all")
plt.savefig("fig/prob.pdf")
plt.close()

"""
xs = csv2mat("xs.csv")
heij = csv2mat("heij.csv")
xkij = csv2mat("xkij.csv")
plt.plot(xs, heij[:,0,0], label="E1")
plt.plot(xs, heij[:,1,1], label="E2")
plt.plot(xs[1:-2], 0.02*xkij[1:-2,0,0,1], label="x12")
plt.legend()
plt.savefig("fig/pot.pdf")
plt.close()
"""

