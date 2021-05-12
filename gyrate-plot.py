#! /usr/bin

import matplotlib as mpl
mpl.use('Agg')

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_style('whitegrid')

df = pd.read_csv('gyrate.xvg', sep = '\s+', skiprows = 27, header = None, engine = 'python')
plt.title('Radius of Gyration')
plt.xlabel('Time (ps)')
plt.ylabel('Rg (nm)')
plt.plot(df[0],df[1])
#plt.legend()
plt.savefig("gyrate-xvg.png")
plt.show()
