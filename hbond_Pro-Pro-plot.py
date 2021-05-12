#! /usr/bin

import matplotlib as mpl
mpl.use('Agg')

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_style('whitegrid')

df = pd.read_csv('hbnum_Pro-Pro.xvg', sep = '\s+', skiprows = 25, header = None, engine = 'python')
plt.title('Hydrogen Bonds within Protein')
plt.xlabel('Time (ps)')
plt.ylabel('Number')
plt.plot(df[0],df[1])
#plt.legend()
plt.savefig("hbnum_Pro-Pro-xvg.png")
plt.show()
