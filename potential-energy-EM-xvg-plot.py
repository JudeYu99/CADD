#! /usr/bin

import matplotlib as mpl
mpl.use('Agg')

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_style('whitegrid')

df = pd.read_csv('potential-energy-EM.xvg', sep = '\s+', skiprows = 24, header = None, engine = 'python')
plt.title('GROMACS Energies')
plt.xlabel('Step')
plt.ylabel('Potential Energy (kJ/mol)')
plt.plot(df[0],df[1], label = 'Potential')
plt.legend()
plt.savefig("potential-energy-EM-xvg.png")
plt.show()
