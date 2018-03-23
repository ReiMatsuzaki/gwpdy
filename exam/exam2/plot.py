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


plt.plot(ts, abs(cs[:,0])**2, label="1")
plt.plot(ts, abs(cs[:,1])**2, label="2")
plt.plot(ts, abs(cs[:,1])**2+abs(cs[:,0])**2, label="all")
plt.ylim(-0.05, 1.05)
plt.savefig("fig/prob.pdf")
plt.close()

