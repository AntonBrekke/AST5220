import numpy as np
import matplotlib.pyplot as plt

# Changing standard color cycle to preferred cycle  
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
colors.insert(6, colors[1]) ; colors[1:5] = colors[2:5]
plt.rcParams['axes.prop_cycle'] = plt.cycler(color=colors) 

# Relevant paths to get data
path = r'/home/antonabr/AST5220/data'
savefig_path = r'/home/antonabr/AST5220/figures/Milestone1'

recombination_data = np.loadtxt(fr'{path}' + r'/recombination.txt')

print(recombination_data.shape)

x, Xe_of_x, ne_of_x, tau_of_x, dtaudx_of_x, ddtauddx_of_x, g_tilde_of_x, dgdx_tilde_of_x, ddgddx_tilde_of_x = recombination_data.T

plt.plot(x, tau_of_x)
plt.plot(x, -dtaudx_of_x)
plt.plot(x, ddtauddx_of_x)
plt.show()