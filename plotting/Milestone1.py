import numpy as np
import matplotlib
# matplotlib.use('TkAgg')
import matplotlib.pyplot as plt 

# Unix
# supernovadata = np.loadtxt('/home/antonabr/AST5220/data/supernovadata.txt')
# results_supernovafitting = np.loadtxt('/home/antonabr/AST5220/data/results_supernovafitting.txt')

# Windows
supernovadata = open(r"Users\anton\OneDrive\Skrivebord\Python\AST5220\data\supernovadata.txt")
results_supernovafitting = np.loadtxt(r'C:/Users\anton\OneDrive\Skrivebord\Python\AST5220\data\results_supernovafitting.txt')

data = np.loadtxt('/home/antonabr/AST5220/data/cosmology.txt')

x = data[:, 0]
eta_of_x = data[:, 1]
detadx_of_x = data[:, 2]
t_of_x = data[:, 3]
Hp_of_x = data[:, 4]

fig = plt.figure()
ax = fig.add_subplot()

Hp_plot = ax.plot(x, detadx_of_x*Hp_of_x)
# ax.set_yscale('log')
plt.savefig('Hp_plot.png')
plt.show()