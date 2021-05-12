#! /usr/bin

import matplotlib as mpl
mpl.use('Agg')

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_style('whitegrid')

df = pd.read_csv('temperature.xvg', sep = '\s+', skiprows = 24, header = None, engine = 'python')
plt.title('Temperature Changes in NVT Equilibration')
plt.xlabel('Time (ps)')
plt.ylabel('Temperature (K)')
plt.plot(df[0],df[1])
#plt.legend()
plt.savefig("temperature-xvg.png")
plt.show()
