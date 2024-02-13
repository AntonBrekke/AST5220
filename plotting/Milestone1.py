import numpy as np
import matplotlib
# matplotlib.use('TkAgg')
import matplotlib.pyplot as plt 

Unix_path = r'/home/antonabr/AST5220/data'
Windows_path = r'C:\Users\anton\OneDrive\Skrivebord\Python\AST5220\data'
# Unix
supernovadata = np.loadtxt(fr'{Unix_path}' + r'\supernovadata.txt')
results_supernovafitting = np.loadtxt(fr'{Unix_path}' + r'\results_supernovafitting.txt')
data = np.loadtxt(fr'{Unix_path}' + r'\cosmology.txt')

# Windows
supernovadata = np.loadtxt(fr'{Windows_path}' + r'\supernovadata.txt')
results_supernovafitting = np.loadtxt(fr'{Windows_path}' + r'\results_supernovafitting.txt')
data = np.loadtxt(fr'{Windows_path}' + r'\cosmology.txt')

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