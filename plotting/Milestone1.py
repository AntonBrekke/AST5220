import numpy as np
import matplotlib.pyplot as plt 
import scipy.constants as scc

# If I want to plot in Unix or Windows
Unix_path = r'/home/antonabr/AST5220/data'
Windows_path = r'C:\Users\anton\OneDrive\Skrivebord\Python\AST5220\data'
path = Unix_path

supernovadata = np.loadtxt(fr'{path}' + r'/supernovadata.txt')
results_supernovafitting = np.loadtxt(fr'{path}' + r'/results_supernovafitting.txt')
data = np.loadtxt(fr'{path}' + r'/cosmology.txt')

# Constants
c = scc.c

x, eta_of_x, detadx_of_x, t_of_x, Hp_of_x, dHpdx_of_x, OmegaB, OmegaCDM, \
OmegaLambda, OmegaR, OmegaNu, OmegaK, luminosity_distance_of_x = data.T

x_start = -14
x_axis_index = np.where(x >= x_start)

fig = plt.figure()
ax = fig.add_subplot()

ax.plot(x[x_axis_index], (eta_of_x*Hp_of_x / c)[x_axis_index])
# ax.plot(x, OmegaR + OmegaNu)
# ax.plot(x, OmegaB + OmegaCDM)
# ax.plot(x, OmegaLambda)
# plt.savefig('Hp_plot.png')
plt.show()