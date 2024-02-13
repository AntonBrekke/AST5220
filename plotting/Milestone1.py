import numpy as np
import matplotlib
# matplotlib.use('TkAgg')
import matplotlib.pyplot as plt 


supernovadata = np.loadtxt('/home/antonabr/AST5220/data/supernovadata.txt')
results_supernovafitting = np.loadtxt('/home/antonabr/AST5220/data/results_supernovafitting.txt')

data = np.loadtxt('/home/antonabr/AST5220/data/cosmology.txt')

x = data[:, 0]
eta_of_x = data[:, 1]
t_of_x = data[:, 2]
Hp_of_x = data[:, 3]

Hp_plot = plt.plot(x, Hp_of_x)
plt.savefig('Hp_plot.png')
plt.show()