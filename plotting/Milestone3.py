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
perturbations_k01_data = np.loadtxt(fr'{path}' + r'/perturbations_k0.1.txt')
perturbations_k001_data = np.loadtxt(fr'{path}' + r'/perturbations_k0.01.txt')
perturbations_k0001_data = np.loadtxt(fr'{path}' + r'/perturbations_k0.001.txt')

# Assign data
x_k01, Theta0_k01, Theta1_k01, Theta2_k01, Phi_k01, Psi_k01, Source_T_k01, v_cdm_k01, v_b_k01 = perturbations_k01_data.T
x_k001, Theta0_k001, Theta1_k001, Theta2_k001, Phi_k001, Psi_k001, Source_T_k001, v_cdm_k001, v_b_k001 = perturbations_k001_data.T
x_k0001, Theta0_k0001, Theta1_k0001, Theta2_k0001, Phi_k0001, Psi_k0001, Source_T_k0001, v_cdm_k0001, v_b_k0001 = perturbations_k0001_data.T

plot_index = np.where(x_k01 <= 0)
plot_start = -15
plot_end = x_k01[plot_index][-1] + 0.005

plt.plot(x_k01, Phi_k01)
plt.plot(x_k001, Phi_k001)
plt.plot(x_k0001, Phi_k0001)
plt.xlim(plot_start, plot_end)
plt.show()

plt.plot(x_k01, Phi_k01 + Psi_k01)
plt.plot(x_k001, Phi_k001 + Psi_k001)
plt.plot(x_k0001, Phi_k0001 + Psi_k0001)
plt.xlim(plot_start, plot_end)
plt.show()

plt.plot(x_k01, Theta0_k01 + Psi_k01)
plt.plot(x_k001, Theta0_k001 + Psi_k001)
plt.plot(x_k0001, Theta0_k0001 + Psi_k0001)
plt.xlim(plot_start, plot_end)
plt.show()
