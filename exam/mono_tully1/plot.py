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

j = json.load(open("out/common.json", "r"))
nt = j["nt"]
dt = j["dt"]
ts = np.arange(nt)*dt

rs = []
for it in range(nt):
    t = it*dt
    out_it = join("out", str(it))
    rs.append( pd.read_csv(join(out_it, "gwp_r.csv"))["val"] )

plt.plot(ts, rs)
if(not exists("fig")):
    os.makedirs("fig")
plt.savefig("fig/r.pdf")

