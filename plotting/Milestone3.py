import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as scc
import tabulate as tab
import latextable as lt
import texttable as txttab

# Relevant paths to get data
path = r'/home/antonabr/AST5220/data'
savefig_path = r'/home/antonabr/AST5220/figures/Milestone3'

# Changing standard color cycle to preferred cycle  
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
colors.insert(6, colors[1]) ; colors[1:5] = colors[2:5]
plt.rcParams['axes.prop_cycle'] = plt.cycler(color=colors) 

# Get data
perturbations_k001_data = np.loadtxt(fr'{path}' + r'/perturbations_k0.01.txt')

# Assign data
x, Theta0, Theta1, Theta2, Phi, Psi, Source_T, Source_T5, Source_T50, Source_T500= perturbations_k001_data.T

plt.plot(x, Theta0)
plt.show()