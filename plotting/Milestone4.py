import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.colors as mc
import healpy as hp

# Relevant paths to get data
path = r'/home/antonabr/AST5220/data'
savefig_path = r'/home/antonabr/AST5220/figures/Milestone4'

# Changing standard color cycle to preferred cycle  
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
colors.insert(5, colors[1]) ; colors[1:5] = colors[2:5]
plt.rcParams['axes.prop_cycle'] = plt.cycler(color=colors) 

Mpc = 3.08567758e22

# Load data
photon = np.loadtxt(fr'{path}' + '/cells.txt')
matter = np.loadtxt(fr'{path}' + '/Pk.txt')
transfer_func = np.loadtxt(fr'{path}' + '/transfer_function.txt')
bessel_func = np.loadtxt(fr'{path}' + '/transfer_function.txt')
low_l_data = np.loadtxt(fr'{path}' + '/low_l_TT_data.txt')
high_l_data = np.loadtxt(fr'{path}' + '/high_l_TT_data.txt')
matter_spectrum_SDSS = np.loadtxt(fr'{path}' + '/matter_spectrum_SDSS.txt')
matter_spectrum_WMAP = np.loadtxt(fr'{path}' + '/matter_spectrum_WMAP.txt')
los_integrand_max = np.loadtxt(fr'{path}' + '/los_integrand_max.txt')

# Assign data
ell, C_ell = photon.T
k, P_k = matter.T
keta0, Theta = transfer_func[:, 0], transfer_func[:, 1:]
z, bessel = bessel_func[:, 0], bessel_func[:, 1:]
ell_TT_low, C_ell_TT_low, err_up_low, err_down_low = low_l_data.T 
ell_TT_high, C_ell_TT_high, err_up_high, err_down_high, bestfit_high = high_l_data.T
k_SDSS, P_k_SDSS, error_SDSS = matter_spectrum_SDSS.T  
k_WMAP, P_k_WMAP, error_WMAP = matter_spectrum_WMAP.T
x_los_max, I_los_max = los_integrand_max[:,0], los_integrand_max[:, 1:], 

test_ells = [6, 100, 200, 500, 1000]

def plot_bessel(show=True):

    fig = plt.figure()
    ax = fig.add_subplot()

    ax.set_title(r'Bessel functions $j_\ell(z)$', fontsize=16)
    ax.plot(z, bessel, lw=0.8)
    ax.set_xlim(0, 1000)
    ax.set_ylim(-0.02, 0.04)
    ax.set_xlabel('z', fontsize=16)
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.grid(True)
    leg = ax.legend([fr'$\ell$={i}' for i in test_ells], prop={'size': 14}, frameon=False)
    for legobj in leg.legend_handles:
        legobj.set_linewidth(2)
    fig.tight_layout()
    plt.savefig(savefig_path + r'/Bessel_functions.pdf')
    if show is True: plt.show()

    fig = plt.figure()
    ax = fig.add_subplot()

    ax.set_title(r'Transfer functions $\Theta_\ell(k, \eta_0)$', fontsize=16)
    ax.plot(keta0, Theta, lw=0.8)
    ax.set_xlim(0, 1500)
    ax.set_ylim(-0.005, 0.015)
    ax.set_xlabel(r'$k\eta_0$', fontsize=16)
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.grid(True)
    leg = ax.legend([fr'$\ell$={i}' for i in test_ells], prop={'size': 14}, frameon=False)
    for legobj in leg.legend_handles:
        legobj.set_linewidth(2)
    fig.tight_layout()
    plt.savefig(savefig_path + r'/transfer_function.pdf')
    if show is True: plt.show()

def plot_integrand(show=True):

    fig = plt.figure()
    ax = fig.add_subplot()

    ax.set_title(r'Integrand $|\Theta_\ell|^2/k\cdot\ell(\ell+1)$', fontsize=16)
    ax.plot(keta0, abs(Theta[:, 0])**2/k * test_ells[0]*(test_ells[0] + 1), lw=0.8)
    ax.plot(keta0, abs(Theta[:, 1])**2/k * test_ells[1]*(test_ells[1] + 1), lw=0.8)
    ax.plot(keta0, abs(Theta[:, 2])**2/k * test_ells[2]*(test_ells[2] + 1), lw=0.8)
    ax.plot(keta0, abs(Theta[:, 3])**2/k * test_ells[3]*(test_ells[3] + 1), lw=0.8)
    ax.plot(keta0, abs(Theta[:, 4])**2/k * test_ells[4]*(test_ells[4] + 1), lw=0.8)
    ax.set_xlim(0, 1500)
    ax.set_xlabel(r'$k\eta_0$', fontsize=16)
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.grid(True)
    leg = ax.legend([fr'$\ell$={i}' for i in test_ells], prop={'size': 14}, frameon=False)
    for legobj in leg.legend_handles:
        legobj.set_linewidth(2)
    fig.tight_layout()
    plt.savefig(savefig_path + r'/integrand.pdf')
    if show is True: plt.show()

def plot_line_of_sight_integrand(show=True):
    
    fig = plt.figure()
    ax = fig.add_subplot()
    ax.set_title(r'Line of sight integrand $\tilde S(k_{\rm max},x)j_{\ell}[k_{\rm max}(\eta_0-\eta(x))]$', fontsize=16)

    ax.plot(x_los_max, I_los_max, lw=0.8)
    ax.set_xlabel(r'x=ln(a)', fontsize=16)
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.set_xlim(-8, 0)
    ax.set_ylim(-5e-5, 5e-5)
    ax.grid(True)
    leg = ax.legend([fr'$\ell$={i}' for i in test_ells], prop={'size': 14}, frameon=False)
    for legobj in leg.legend_handles:
        legobj.set_linewidth(2)
    fig.tight_layout()
    plt.savefig(savefig_path + r'/line_of_sight_integrand.pdf')
    if show is True: plt.show()

def plot_power_spectrum(show=True):

    fig = plt.figure()
    ax = fig.add_subplot()

    cosmic_variance = np.sqrt(2/(2*ell + 1))*C_ell

    ax.set_title(r'CMB Power spectrum $C_\ell\cdot\frac{\ell(\ell+1)}{2\pi}(10^6T_{\rm CMB0})^2$', fontsize=16)
    ax.plot(ell, C_ell, lw=2)
    ax.fill_between(x=ell, y1=C_ell-cosmic_variance, y2=C_ell+cosmic_variance, label=r'$\sqrt{\text{Var}(C_\ell)}$', color='limegreen', edgecolor='darkgreen', alpha=0.3)
    # ax.plot(ell**(1.018), C_ell * 1.13*np.exp(-0.05*(ell/200)**1.5), lw=2)        # Cheating values
    ax.errorbar(ell_TT_low, C_ell_TT_low, yerr=[err_down_low, err_up_low], label='Planck spectrum', ls='', marker='o', ms=3, ecolor='tab:red', color='k', capsize=2, lw=1.5)
    ax.errorbar(ell_TT_high, C_ell_TT_high, yerr=[err_down_high, err_up_high], ls='', marker='o', ms=3, ecolor='tab:red', color='k', capsize=2, lw=1.5)
    ax.set_xscale('log')
    ax.set_xlabel(r'Multipole $\ell$', fontsize=16)
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.grid(True)
    ax.legend(prop={'size':14}, frameon=False, loc='upper left')
    fig.tight_layout()
    plt.savefig(savefig_path + r'/power_spectrum.pdf')
    if show is True: plt.show()

    fig = plt.figure()
    ax = fig.add_subplot()
    
    k_eq = 3.3634935539305665e-25 * Mpc/(0.67)          # k_eq = a_eq*H(a_eq) / c * (Mpc/h). a_eq*H(a_eq) printed in Milestone1.py, h=0.67.

    ax.set_title(r'Matter Power spectrum $P(k)\;(\text{Mpc}/h)^3$', fontsize=16)
    ax.plot(k, P_k ,lw=2)
    ax.axvline(k_eq, color='k', ls='--', lw=2, label=r'$k_{\rm eq}$')
    ax.errorbar(k_SDSS, P_k_SDSS, yerr=error_SDSS, label='SDSS Galaxies (DR7 LRG)', ls='', marker='o', ms=3, ecolor='tab:red', color='k', capsize=2, lw=1.5)
    ax.errorbar(k_WMAP, P_k_WMAP, yerr=abs(P_k_WMAP - error_WMAP), label='WMAP+ACT', ls='', marker='o', ms=3, ecolor='tab:orange', color='k', capsize=2, lw=1.5)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(2e-3, k[-1])
    ax.set_xlabel(r'$k\;(h/\rm Mpc)$', fontsize=16)
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.grid(True)
    ax.legend(prop={'size':14}, frameon=False, loc='lower left')
    fig.tight_layout()
    plt.savefig(savefig_path + r'/matter_spectrum.pdf')
    if show is True: plt.show()

def plot_CMB(show=True):
    # Normalize spectrum back
    C_ell_normal = C_ell/(ell*(ell+1))*(2*np.pi)

    # Make my own CMAP for the CMB 
    color_factor = 2      # Scales amplitude of total colormap
    get_cmap = plt.get_cmap('RdYlBu_r')
    colors = get_cmap(np.linspace(0, 1, 256))       # 256 not too important, just because print(get_cmap.N) -> 256
    # RdYlBu_r gets too yellow when amplifying colors. Modify map
    # cmap = np.array([[0, 0, 0.4], *colors[10:40:15, :3],                           # Steal from RdYlBu_r
    #                 [0, 0.8, 1], [1, 1, 0.9], [1, 0.7, 0],         # Light blue, Yellow, Orange
    #                  *(colors)[256-30::15, :3], [0,0,0]])**color_factor       # Steal from RdYlBu_r

    cmap = np.array([*colors[10:40:15, :3],                           # Steal from RdYlBu_r
                    [0, 0.8, 1], [1, 1, 0.9], [1, 0.7, 0],         # Light blue, Yellow, Orange
                     *(colors**(2))[256-30::15, :3]])**color_factor       # Steal from RdYlBu_r

    CMB_cmap = mc.LinearSegmentedColormap.from_list('CMB_cmap', cmap, N=300)
    # CMB_cmap = 'jet'

    # Generate Gaussian random fields
    nside = int(2**10)       # Resolution of image    
    print("Approximate resolution at NSIDE {} is {:.2} deg".format(nside, hp.nside2resol(nside, arcmin=True) / 60))     # stole this from https://healpy.readthedocs.io/en/latest/tutorial.html
    lmax = int(np.max(ell))
    seed = 5220         # Set seed to get same plot each time.
    np.random.seed(seed)
    cmb_map = hp.synfast(C_ell_normal, nside=nside, lmax=lmax, mmax=lmax)
    """
    If you are analyzing a map from an instrument with a specific beam width, 
    you can correct the power spectrum by the smoothing factor caused by that 
    beam and obtain a better approximation 
    of the power spectrum of the original sky.
    https://www.zonca.dev/posts/2021-04-27-correct-beam-healpy
    """
    deg = 9/60     # arcmin
    cmb_map = hp.smoothing(cmb_map, sigma=deg*np.pi/180)

    # Generate and plot the map
    # cmb_map = hp.alm2map(alm, nside=nside, lmax=lmax, mmax=lmax)

    fig = plt.figure()
    fig.suptitle('CMB map', fontsize=20, y=0.95)
    mollview = hp.mollview(cmb_map, fig=fig, cmap=CMB_cmap, title='', cbar=False, format='%.3e', return_projected_map=True)

    # Get max and min - dummy proof (there is -inf at least)
    mollview_min = np.nanmin(mollview.data[np.logical_and(mollview.data!=np.inf, mollview.data!=-np.inf)])
    mollview_max = np.nanmax(mollview.data[np.logical_and(mollview.data!=np.inf, mollview.data!=-np.inf)])

    # Make my OWN cbar since the one from mollview is small
    ax = plt.gca()      # Get cbar from mollview
    image = ax.get_images()[0]

    format_func = lambda x, pos: f'{x:.2f}'
    cbar_L = 0.8
    cbar_ax = fig.add_axes([0.5 - cbar_L/2, 0.1, cbar_L, 0.04]) 
    cbar = fig.colorbar(image, cax=cbar_ax, ticks=[mollview_min, mollview_max], format=format_func, orientation='horizontal')
    cbar.ax.tick_params(labelsize=14)
    cbar.ax.set_xlabel(r'$[\mu K]$', rotation=0, fontsize=16, horizontalalignment='center', verticalalignment='baseline')
    plt.savefig(savefig_path + r'/CMB.pdf')
    if show is True: plt.show()

# plot_bessel(show=True)
plot_integrand(show=True)
# plot_line_of_sight_integrand(show=True)
# plot_power_spectrum(show=True)
# plot_CMB(show=True)

