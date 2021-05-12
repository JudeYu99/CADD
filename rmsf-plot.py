#! /usr/bin

import matplotlib as mpl
mpl.use('Agg')

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_style('whitegrid')

df = pd.read_csv('rmsf.xvg', sep = '\s+', skiprows = 17, header = None, engine = 'python')
df['ResID'] = df.index + 1
plt.title('RMSF for C-alpha')
plt.xlabel('Residue')
plt.ylabel('Fluctuation (nm)')
plt.plot(df['ResID'], df[1])
#plt.legend()
plt.savefig("rmsf-xvg.png")
plt.show()

