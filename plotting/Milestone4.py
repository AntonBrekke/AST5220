import numpy as np
import matplotlib.pyplot as plt 

# Relevant paths to get data
path = r'/home/antonabr/AST5220/data'
savefig_path = r'/home/antonabr/AST5220/figures/Milestone3'

# Changing standard color cycle to preferred cycle  
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
colors.insert(6, colors[1]) ; colors[1:5] = colors[2:5]
plt.rcParams['axes.prop_cycle'] = plt.cycler(color=colors) 

Mpc = 3.08567758e22

photon = np.loadtxt(fr'{path}' + '/cells.txt')
matter = np.loadtxt(fr'{path}' + '/Pk.txt')
transfer_func = np.loadtxt(fr'{path}' + '/transfer_function.txt')
bessel_func = np.loadtxt(fr'{path}' + '/transfer_function.txt')

l, C_l = photon.T
k, P_k = matter.T

keta0, Theta = transfer_func[:, 0], transfer_func[:, 1:]

z, bessel = bessel_func[:, 0], bessel_func[:, 1:]

plt.plot(z, bessel)
plt.xlim(0, 1000)
plt.ylim(-0.02, 0.04)
plt.show()

plt.plot(keta0, Theta)
plt.xlim(0, 1500)
plt.ylim(-0.005, 0.015)
plt.show()

# plt.plot(k, abs(Theta)**2/k)
# plt.xlim(0, 1500)
# plt.show()

plt.plot(l, l*(l+1)/(2*np.pi)*C_l)
plt.show()

plt.plot(k, P_k)
plt.show()



