import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.colors as mc
import healpy as hp

# Relevant paths to get data
path = r'/home/antonabr/AST5220/data'
savefig_path = r'/home/antonabr/AST5220/figures/Milestone4'

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

def plot_bessel(show=True):
    plt.plot(z, bessel, lw=0.8)
    plt.xlim(0, 1000)
    plt.ylim(-0.02, 0.04)
    if show is True: plt.show()

    plt.plot(keta0, Theta, lw=0.8)
    plt.xlim(0, 1500)
    plt.ylim(-0.005, 0.015)
    if show is True: plt.show()

def plot_integrand(show=True):
    plt.plot(keta0, 6*(6+1)*abs(Theta[:, 0])**2/k, lw=0.8)
    plt.plot(keta0, 100*(100+1)*abs(Theta[:, 1])**2/k, lw=0.8)
    plt.plot(keta0, 200*(200+1)*abs(Theta[:, 2])**2/k, lw=0.8)
    plt.plot(keta0, 500*(500+1)*abs(Theta[:, 3])**2/k, lw=0.8)
    plt.xlim(0, 1500)
    if show is True: plt.show()

def plot_power_spectrum(show=True):
    plt.plot(ell, C_ell)
    plt.errorbar(ell_TT_low, C_ell_TT_low, yerr=[err_down_low, err_up_low], ls='', marker='o', ms=1.5, ecolor='tab:red', color='k', capsize=2)
    plt.errorbar(ell_TT_high, C_ell_TT_high, yerr=[err_down_high, err_up_high], ls='', marker='o', ms=1.5, ecolor='tab:red', color='k', capsize=2)
    plt.xscale('log')
    plt.yscale('log')
    if show is True: plt.show()

    plt.plot(k, P_k)
    plt.xscale('log')
    plt.yscale('log')
    if show is True: plt.show()

def plot_CMB(show=True):
    # Normalize spectrum back
    C_ell_normal = C_ell/(ell*(ell+1))*(2*np.pi)/(1e6*2.7255)**2

    # Make my own CMAP for the CMB 
    color_factor = 2
    get_cmap = plt.get_cmap('RdYlBu_r')
    colors = get_cmap(np.linspace(0, 1, 256))       # 256 not too important, just because print(get_cmap.N) -> 256
    # RdYlBu_r gets too yellow when amplifying colors. Modify map
    cmap = np.array([*colors[:30:15, :3],                           # Steal from RdYlBu_r
                    [0, 0.8, 1], [1, 1, 0.9], [1, 0.7, 0],        # Light blue, Yellow, Orange
                     *colors[256-30::15, :3]])**color_factor       # Steal from RdYlBu_r

    CMB_cmap = mc.LinearSegmentedColormap.from_list('CMB_cmap', cmap, N=300)

    # Generate Gaussian random fields
    nside = int(2**10)       # Resolution of image    
    print("Approximate resolution at NSIDE {} is {:.2} deg".format(nside, hp.nside2resol(nside, arcmin=True) / 60))     # stole this from https://healpy.readthedocs.io/en/latest/tutorial.html
    lmax = int(np.max(ell))
    seed = 15012001         # Set seed to get same plot each time.
    np.random.seed(seed)
    alm = hp.synalm(C_ell_normal, lmax=lmax)

    # Generate and plot the map
    cmb_map = hp.alm2map(alm, nside=nside, lmax=lmax)

    fig = plt.figure()
    fig.suptitle('CMB map', fontsize=20, y=0.95)
    mollview = hp.mollview(cmb_map, fig=fig, cmap=CMB_cmap, title='', cbar=False, format='%.3e', return_projected_map=True)

    # Get max and min - dummy proof (there is -inf at least)
    mollview_min = np.nanmin(mollview.data[np.logical_and(mollview.data!=np.inf, mollview.data!=-np.inf)])
    mollview_max = np.nanmax(mollview.data[np.logical_and(mollview.data!=np.inf, mollview.data!=-np.inf)])

    # Make my OWN cbar since the one from mollview is small
    ax = plt.gca()      # Get cbar from mollview
    image = ax.get_images()[0]

    format_func = lambda x, pos: f'{x:.3e}'
    cbar_L = 0.7
    cbar_ax = fig.add_axes([0.5 - cbar_L/2, 0.1, cbar_L, 0.04]) 
    cbar = fig.colorbar(image, cax=cbar_ax, ticks=[mollview_min, mollview_max], format=format_func, orientation='horizontal')
    cbar.ax.tick_params(labelsize=14)
    plt.savefig(savefig_path + r'/CMB.pdf')
    if show is True: plt.show()

# plot_bessel(show=False)
# plot_integrand(show=False)
# plot_power_spectrum(show=False)
plot_CMB(show=True)

