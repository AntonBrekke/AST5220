import numpy as np
import matplotlib.pyplot as plt 
import seaborn as sns 
import scipy.constants as scc

# Changing standardcolor cycle to preferred cycle  
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
colors.insert(6, colors[1]) ; colors[1:5] = colors[2:5]
plt.rcParams['axes.prop_cycle'] = plt.cycler(color=colors) 

# If I want to plot in Unix or Windows
Unix_path = r'/home/antonabr/AST5220/data'
Windows_path = r'C:\Users\anton\OneDrive\Skrivebord\Python\AST5220\data'
path = Unix_path

# Get data
supernovadata = np.loadtxt(fr'{path}' + r'/supernovadata.txt')
results_supernovafitting = np.loadtxt(fr'{path}' + r'/results_supernovafitting.txt', skiprows=200)
cosmology = np.loadtxt(fr'{path}' + r'/cosmology.txt')

# Assign data
x, eta_of_x, detadx_of_x, t_of_x, Hp_of_x, dHpdx_of_x, ddHpddx_of_x, OmegaB, OmegaCDM, \
OmegaLambda, OmegaR, OmegaNu, OmegaK, luminosity_distance_of_x, z_of_x = cosmology.T

chi2, h, OmegaM, OmegaK = results_supernovafitting.T

# (z, Gpc, Gpc), z is redshift 
z, luminosity_distance_of_z, luminosity_error = supernovadata.T

# Constants
c = scc.c       # m/s
parsec = scc.parsec     # 1 parsec in m

# Units of convertion. In General, we plot H(B) / B, so just insert B -> A for given A. 
# Going from Gyr to s 
Gyr_to_seconds = 365*24*60*60*1e9             
# Going from 1/s to 100km/(Mpc*s)     
pr_second_to_km_pr_second_pr_Mparsec = 100e3 / (parsec*1e6)    


# Make small class to not copy too much code. Not meant to be as general as possible
class make_plot:
    def __init__(self, x_start=None, x_end=None, title=''):
        self.x_start = x_start
        self.x_end = x_end
        self.title = title
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot() 
        self.ax.set_title(self.title, fontsize=16)

    def plot(self, x, data, label='', **kwargs):
        if (self.x_start) is None: self.x_start = x[0]
        if (self.x_end) is None: self.x_end = x[-1]
        self.x_axis_index = np.logical_and(x >= self.x_start, x <= self.x_end)
        self.ax.plot(x[self.x_axis_index], (data)[self.x_axis_index], label=label, **kwargs)
        if label != '':
            self.ax.legend(prop={'size': 12})
        # Call show when ready 

    def hist(self, data, bins, label='', density=False):
        self.counts, bins = np.histogram(data, bins=bins)
        self.n, self.bins, self.patches = self.ax.hist(bins[:-1], bins, weights=self.counts, density=density)
        cm = plt.get_cmap('Spectral_r')
        cmap = cm((self.n - np.min(self.n))/(np.max(self.n) - np.min(self.n)))
        for color, patch in zip(cmap, self.patches):
            plt.setp(patch, 'facecolor', color)
        if label != '':
            self.ax.legend(prop={'size': 12})
        # Call show when ready

    def format_plot(self, xlabel='', ylabel='', xscale='linear', yscale='linear'):
        self.ax.set_xlabel(xlabel, fontsize=16)
        self.ax.set_ylabel(ylabel, fontsize=16)
        self.ax.set_xscale(xscale)
        self.ax.set_yscale(yscale)
        self.ax.set_xlim(self.x_start, self.x_end)
        self.fig.tight_layout()
# End of class 

a = np.exp(x)       # Scaling factor instead of x may be more intuitive

# Call all plots 
def plot_demonstrate_code():
    label_dHpdx = r'$\frac{1}{\mathcal{H}(x)}\frac{d\mathcal{H}(x)}{dx}\;$'
    label_ddHpddx = r'$\frac{1}{\mathcal{H}(x)}\frac{d^2\mathcal{H}(x)}{dx^2}\;$'
    title = r'Evolution of derivatives of $\mathcal{H}(x)\equiv aH(x)$'
    plot_Hp_derivatives = make_plot(title=title)
    plot_Hp_derivatives.plot(x, dHpdx_of_x/Hp_of_x, label=label_dHpdx)
    plot_Hp_derivatives.plot(x, ddHpddx_of_x/Hp_of_x, label=label_ddHpddx)

    plot_Hp_derivatives.format_plot('x=ln(a)')
    plt.show()

    title = r'$\frac{\eta(x)\mathcal{H}(x)}{c}\;$'
    plot_etaHp_pr_c = make_plot(-14, 0, title=title)
    plot_etaHp_pr_c.plot(x, eta_of_x*Hp_of_x/c)
    plot_etaHp_pr_c.format_plot('x=ln(a)')
    plt.show()


def plot_conformal_Hubble():
    title = r'$\mathcal{H}(x)\;\left(\frac{100km/s}{Mpc}\right)$'
    plot_Hp = make_plot(-12, 0, title=title)
    plot_Hp.plot(x, Hp_of_x / pr_second_to_km_pr_second_pr_Mparsec)
    plot_Hp.format_plot('x=ln(a)', yscale='log')
    plt.show()


def plot_conformal_time_pr_c():
    title = r'$\frac{\eta(x)}{c}\;\left(\frac{Mpc}{c}\right)$'
    plot_eta_pr_c = make_plot(-12, 0, title=title)
    plot_eta_pr_c.plot(x, eta_of_x/(c*parsec))
    plot_eta_pr_c.format_plot('x=ln(a)', yscale='log')
    plt.show()


def plot_time():
    plot_t = make_plot(-12, 0, title=r'$t(x)\;(Gyr)$')
    plot_t.plot(x, t_of_x / Gyr_to_seconds)
    plot_t.format_plot('x=ln(a)', yscale='log')
    plt.show()


def plot_densities():
    plot_densities = make_plot(x[0], x[-1], title=r'$\Omega_i$')
    plot_densities.plot(x, OmegaR + OmegaNu, label=r'$\Omega_R + \Omega_\nu$')
    plot_densities.plot(x, OmegaB + OmegaCDM, label=r'$\Omega_B + \Omega_{CDM}$')
    plot_densities.plot(x, OmegaLambda, label=r'$\Omega_\Lambda$')
    plot_densities.format_plot('x=ln(a)')
    plt.show()


def plot_luminosity_distance_of_z():
    fig = plt.figure()
    ax = fig.add_subplot()
    ax.set_title(r'Luminosity distance $d_L / z$', fontsize=16)
    ax.errorbar(z, luminosity_distance_of_z / z, ls='', marker='o', ms=3, yerr=luminosity_error/z, ecolor='k', capsize=0, color="k", label=r'Supernovadata')
    ax.plot(z_of_x, luminosity_distance_of_x / (z_of_x*parsec*1e9), color='r', label='Theoretical')
    ax.set_xlim(1.1*z[0], 1.1*z[-1])
    ax.set_ylim(4, 8)
    ax.set_xscale('log')
    ax.set_xlabel('Redshift z', fontsize=16)
    ax.set_ylabel(r'$d_L$/z (Mpc)', fontsize=16)

    ax.legend(prop={'size': 12})
    fig.tight_layout()
    plt.show()


def plot_supernovadata_MCMC_fits():
    chi2_min = np.min(chi2)

    chi2_1sigma = np.where((chi2 - chi2_min) < 3.53)    # 1 sigma
    chi2_2sigma = np.where((chi2 - chi2_min) < 8.02)    # 2 sigma

    OmegaM_selected_1sigma = OmegaM[chi2_1sigma]
    OmegaM_selected_2sigma = OmegaM[chi2_2sigma]
    OmegaLambda_selected_1sigma = (1 - (OmegaM + OmegaK))[chi2_1sigma]
    OmegaLambda_selected_2sigma = (1 - (OmegaM + OmegaK))[chi2_2sigma]

    line = -1*np.linspace(0, 1, len(OmegaLambda_selected_2sigma)) + 1
    plt.scatter(OmegaM_selected_2sigma, OmegaLambda_selected_2sigma, label=r'$2\sigma$')
    plt.scatter(OmegaM_selected_1sigma, OmegaLambda_selected_1sigma, label=r'$1\sigma$')
    plt.plot((0,1), (1,0), 'k', ls='--', label='Flat Universe')

    plt.xlabel(r'$\Omega_{M}$', fontsize=16)
    plt.ylabel(r'$\Omega_{\Lambda}$', fontsize=16)
    plt.legend(prop={'size':12})
    plt.show()


def plot_posterior_PDF_Hubble_param():
    # H0 = h / Constants.H0_over_h, see BackgroundCosmology.cpp. PDF will be same.
    posterior_H0_pdf = make_plot(title=r'Posterior PDF for $H_0$')
    posterior_H0_pdf.hist(h, bins=75, density=True)
    bins = posterior_H0_pdf.bins
    # print(np.sum(np.diff(posterior_H0_pdf.bins) * posterior_H0_pdf.n))   # Testing if prop. dist sum to 1 

    sigma = np.std(h)
    mu = np.mean(h)
    gaussian = 1/(sigma*np.sqrt(2*np.pi))*np.exp(-(bins-mu)**2/(2*sigma**2))
    posterior_H0_pdf.plot(bins, gaussian, color='k', lw=2.5)
    posterior_H0_pdf.format_plot(xlabel=r'$H_0$')

    plt.show()

def plot_posterior_PDF_OmegaLambda():
    OmegaLambda = 1 - (OmegaM + OmegaK)
    posterior_OmegaLambda_pdf = make_plot(title=r'Posterior PDF for $\Omega_{\Lambda}$')
    posterior_OmegaLambda_pdf.hist(OmegaLambda, bins=75, density=True)
    bins = posterior_OmegaLambda_pdf.bins
    # print(np.sum(np.diff(posterior_OmegaLambda_pdf.bins) * posterior_OmegaLambda_pdf.n))   # Testing if prop. dist sum to 1 


    sigma = np.std(OmegaLambda)
    mu = np.mean(OmegaLambda)
    gaussian = 1/(sigma*np.sqrt(2*np.pi))*np.exp(-(bins-mu)**2/(2*sigma**2))
    posterior_OmegaLambda_pdf.plot(bins, gaussian, color='k', lw=2.5)
    posterior_OmegaLambda_pdf.format_plot(xlabel=r'$\Omega_{\Lambda}$')

    plt.show()

# Control unit for plotting 
# plot_demonstrate_code()
# plot_conformal_Hubble()
# plot_conformal_time_pr_c()
# plot_time()
plot_densities()
plot_luminosity_distance_of_z()
plot_supernovadata_MCMC_fits()
plot_posterior_PDF_Hubble_param()
plot_posterior_PDF_OmegaLambda()

