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

rs = correct_csvs("r", nt)
cs = correct_csvs("c", nt)

plt.plot(ts, rs[:,0])
if(not exists("fig")):
    os.makedirs("fig")
plt.savefig("fig/r.pdf")
plt.close()

plt.plot(ts, abs(cs[:,0])**2)
plt.savefig("fig/prob.pdf")
plt.close()

xs = csv2mat("xs.csv")
heij = csv2mat("heij.csv")
xkij = csv2mat("xkij.csv")
plt.plot(xs, heij[:,0,0], label="E1")
plt.plot(xs, heij[:,1,1], label="E2")
plt.plot(xs[1:-2], 0.02*xkij[1:-2,0,0,1], label="x12")
plt.legend()
plt.savefig("fig/pot.pdf")
plt.close()

