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

def correct_csvs(name, nt, dir_out="out"):
    """
    """
    data = []
    for it in range(nt):
        dir_out_it = join(dir_out, str(it))
        df = pd.read_csv(join(dir_out_it, name+".csv"))
        if("val" in df.columns):
            data.append(df["val"])
        elif("re" in df.columns):
            data.append(np.array(df["re"])+1j*np.array(df["im"]))
    return np.array(data)
