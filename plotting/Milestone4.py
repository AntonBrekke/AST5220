import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.colors as mc
import healpy as hp

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
low_l_data = np.loadtxt(fr'{path}' + '/low_l_TT_data.txt')
high_l_data = np.loadtxt(fr'{path}' + '/high_l_TT_data.txt')

ell, C_ell = photon.T
k, P_k = matter.T
keta0, Theta = transfer_func[:, 0], transfer_func[:, 1:]
z, bessel = bessel_func[:, 0], bessel_func[:, 1:]
ell_TT_low, C_ell_TT_low, err_up_low, err_down_low = low_l_data.T 
ell_TT_high, C_ell_TT_high, err_up_high, err_down_high, bestfit_high = high_l_data.T 

plt.plot(z, bessel, lw=0.8)
plt.xlim(0, 1000)
plt.ylim(-0.02, 0.04)
plt.show()

plt.plot(keta0, Theta, lw=0.8)
plt.xlim(0, 1500)
plt.ylim(-0.005, 0.015)
plt.show()

plt.plot(keta0, 6*(6+1)*abs(Theta[:, 0])**2/k, lw=0.8)
plt.plot(keta0, 100*(100+1)*abs(Theta[:, 1])**2/k, lw=0.8)
plt.plot(keta0, 200*(200+1)*abs(Theta[:, 2])**2/k, lw=0.8)
plt.plot(keta0, 500*(500+1)*abs(Theta[:, 3])**2/k, lw=0.8)
plt.xlim(0, 1500)
plt.show()

plt.plot(ell, C_ell)
plt.errorbar(ell_TT_low, C_ell_TT_low, yerr=[err_down_low, err_up_low], ls='', marker='o', ms=1.5, ecolor='tab:red', color='k', capsize=2)
plt.errorbar(ell_TT_high, C_ell_TT_high, yerr=[err_down_high, err_up_high], ls='', marker='o', ms=1.5, ecolor='tab:red', color='k', capsize=2)
plt.xscale('log')
plt.yscale('log')
plt.show()

plt.plot(k, P_k)
plt.xscale('log')
plt.yscale('log')
plt.show()

nside = int(1028)
lmax = int(np.max(ell))

# Make my own CMAP for the CMB 
color_factor = 2

get_cmap = plt.get_cmap('RdYlBu_r')
RdYlBu_r = get_cmap(np.linspace(0, 1, C_ell.size))**color_factor
CMB_cmap = mc.ListedColormap(RdYlBu_r)

# Generate Gaussian random fields
alm = hp.synalm(C_ell/(ell*(ell+1)), lmax=lmax)

# Generate the map
cmb_map = hp.alm2map(alm, nside=nside, lmax=lmax)

# Plot the map
hp.mollview(cmb_map, cmap=CMB_cmap, title='CMB Map')
plt.show()



