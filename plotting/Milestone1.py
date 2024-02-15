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

x, eta_of_x, detadx_of_x, t_of_x, Hp_of_x, dHpdx_of_x, ddHpddx_of_x, OmegaB, OmegaCDM, \
OmegaLambda, OmegaR, OmegaNu, OmegaK, luminosity_distance_of_x = data.T

x_start = -12
x_axis_index = np.where(x >= x_start)

# fig = plt.figure()
# ax = fig.add_subplot()

plt.plot(x, dHpdx_of_x/Hp_of_x)
plt.show()
plt.plot(x, ddHpddx_of_x/Hp_of_x)
plt.show()
plt.plot(x[x_axis_index], (eta_of_x*Hp_of_x/c)[x_axis_index], label=r'$\eta(x)\mathcal{H}(x)/c$')
plt.show()
plt.plot(x[x_axis_index], Hp_of_x[x_axis_index], label=r'$\mathcal{H}(x)$')
plt.yscale('log')
plt.show()
plt.plot(x[x_axis_index], t_of_x[x_axis_index], label='t(x)')
plt.yscale('log')
plt.show()
plt.plot(x[x_axis_index], (eta_of_x/c)[x_axis_index], label=r'$\eta(x)$')
plt.yscale('log')
plt.show()
plt.plot(x, OmegaR + OmegaNu)
plt.plot(x, OmegaB + OmegaCDM)
plt.plot(x, OmegaLambda)
plt.show()