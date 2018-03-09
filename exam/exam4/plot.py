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

rs = []
ps = []
prob1 = []
for it in range(nt):
    dir_it = join("out", str(it))
    if(not exists(dir_it)):
        break
    r = csv2mat(join(dir_it, "r.csv"))
    p = csv2mat(join(dir_it, "p.csv"))
    c = csv2mat(join(dir_it, "c.csv"))
    cc = csv2mat(join(dir_it, "cc.csv"))
    npath = len(cc)    

    if(npath==1):
        rs.append([r[0,0], 0.0])
        ps.append([p[0,0], 0.0])
    else:
        rs.append([r[0,0], r[1,0]])
        ps.append([p[0,0], p[1,0]])
        
    p1 = abs(np.sum(cc[:]*c[:,0]))**2
    prob1.append(p1)

rs = np.array(rs)
n = len(rs[:,0])
plt.plot(ts[:n], rs[:,0], label="1")
plt.plot(ts[:n], rs[:,1], label="2")
if(not exists("fig")):
    os.makedirs("fig")
plt.savefig("fig/r.pdf")
plt.close()

ps = np.array(ps)
plt.plot(ts[:n], ps[:,0], label="1")
plt.plot(ts[:n], ps[:,1], label="2")
plt.ylim(14,16)
plt.savefig("fig/p.pdf")
plt.close()

plt.plot(ts, prob1)
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

