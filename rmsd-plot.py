#! /usr/bin

import matplotlib as mpl
mpl.use('Agg')

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_style('whitegrid')

df = pd.read_csv('rmsd.xvg', sep = '\s+', skiprows = 18, header = None, engine = 'python')
plt.title('RMSD for C-alpha')
plt.xlabel('Time (ns)')
plt.ylabel('RMSD (nm)')
plt.plot(df[0],df[1])
#plt.legend()
plt.savefig("rmsd-xvg.png")
plt.show()
