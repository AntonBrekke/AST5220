import numpy as np
import matplotlib.pyplot as plt 
import scipy.constants as scc

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
OmegaLambda, OmegaR, OmegaNu, OmegaK, luminosity_distance_of_x = cosmology.T

chi2, h, OmegaM, OmegaK = results_supernovafitting.T

# (z, Gpc, Gpc), z is redshift 
z, luminosity_distance_of_z, Error = supernovadata.T

# Constants
c = scc.c       # m/s
parsec = scc.parsec     # 1 parsec in m

# Units of convertion. In General, we plot H(B) / B, so just insert B -> A for given A. 
# Going from Gyr to s 
Gyr_to_seconds = 365*24*60*60*1e9             
# Going from 1/s to 100km/(Mpc*s)     
pr_second_to_km_pr_second_pr_Mparsec = 100e3 / (parsec*1e6)        


# Make small class to not copy too much code 
class make_plot:
    def __init__(self, x_start=0, x_end=1, title=''):
        self.x_start = x_start
        self.x_end = x_end
        self.title = title 

    def plot(self, x, data, use_prev_fig=False, label='', **kwargs):
        if use_prev_fig == False:
            self.fig = plt.figure()
            self.ax = self.fig.add_subplot()
            self.ax.set_title(self.title, fontsize=16)
        self.x_axis_index = np.logical_and(x >= self.x_start, x <= self.x_end)
        self.ax.plot(x[self.x_axis_index], (data)[self.x_axis_index], label=label, **kwargs)
        if label != '':
            self.ax.legend()
        # Call show when ready 

    def hist(self, data, bins, label='', density=False, use_prev_fig=False):
        self.counts, self.bins = np.histogram(data, bins=bins)
        if use_prev_fig == False:
            self.fig = plt.figure()
            self.ax = self.fig.add_subplot()
            self.ax.set_title(self.title, fontsize=16)
        self.ax.hist(self.bins[:-1], self.bins, weights=self.counts, density=density)
        if label != '':
            self.ax.legend()
        # Call show when ready

    def format_plot(self, xlabel='', ylabel='', scale='linear'):
        self.ax.set_xlabel(xlabel)
        self.ax.set_ylabel(ylabel)
        self.ax.set_yscale(scale)
        self.ax.set_xlim(self.x_start, self.x_end)
        self.fig.tight_layout()
# End of class 

# Call all plots 
def plot_demonstrate_code():
    plot_dHpdx_pr_Hp = make_plot(-12, 0)
    plot_dHpdx_pr_Hp.plot(x, dHpdx_of_x/Hp_of_x*pr_second_to_km_pr_second_pr_Mparsec)
    plot_dHpdx_pr_Hp.format_plot()
    plt.show()

    plot_ddHpddx_pr_Hp = make_plot(-12, 0)
    plot_ddHpddx_pr_Hp.plot(x, ddHpddx_of_x/Hp_of_x*pr_second_to_km_pr_second_pr_Mparsec)
    plot_ddHpddx_pr_Hp.format_plot()
    plt.show()

    plot_etaHp_pr_c = make_plot(-12, 0, title=r'$\frac{\eta(x)\mathcal{H}}{c}$')
    plot_etaHp_pr_c.plot(x, eta_of_x*Hp_of_x/c)
    plot_etaHp_pr_c.format_plot()
    plt.show()


def plot_conformal_Hubble():
    plot_Hp = make_plot(-12, 0, title=r'$\mathcal{H}(x)\;\left(\frac{100km/s}{Mpc}\right)$')
    plot_Hp.plot(x, Hp_of_x / pr_second_to_km_pr_second_pr_Mparsec)
    plot_Hp.format_plot(scale='log')
    plt.show()


def plot_conformal_time_pr_c():
    plot_eta_pr_c = make_plot(-12, 0, title=r'$\frac{\eta(x)}{c}\;\left(\frac{Mpc}{c}\right)$')
    plot_eta_pr_c.plot(x, eta_of_x/(c*parsec))
    plot_eta_pr_c.format_plot(scale='log')
    plt.show()


def plot_time():
    plot_t = make_plot(-12, 0, title=r'$t(x)\;(Gyr)$')
    plot_t.plot(x, t_of_x / Gyr_to_seconds)
    plot_t.format_plot(scale='log')
    plt.show()


def plot_densities():
    plot_densities = make_plot(x[0], x[-1], title=r'$\Omega_i$')
    plot_densities.plot(x, OmegaR + OmegaNu, label=r'$\Omega_R + \Omega_\nu$')
    plot_densities.plot(x, OmegaB + OmegaCDM, label=r'$\Omega_B + \Omega_{CDM}$', use_prev_fig=True)
    plot_densities.plot(x, OmegaLambda, label=r'$\Omega_\Lambda$', use_prev_fig=True)
    plot_densities.format_plot()
    plt.show()


def plot_luminosity_distance_of_z():
    fig = plt.figure()
    ax = fig.add_subplot()
    ax.set_title('')
    ax.errorbar(z, luminosity_distance_of_z, yerr=Error, ecolor='r', capsize=4, color="k")

    fig.tight_layout()
    plt.show()


def plot_supernovadata_MCMC_fits():
    chi2_min = np.min(chi2)
    print(chi2_min)

    chi2_1sigma = np.where((chi2 - chi2_min) < 3.53)    # 1 sigma
    chi2_2sigma = np.where((chi2 - chi2_min) < 8.02)    # 2 sigma

    OmegaM_selected_1sigma = OmegaM[chi2_1sigma]
    OmegaM_selected_2sigma = OmegaM[chi2_2sigma]
    OmegaLambda_selected_1sigma = (1 - (OmegaM + OmegaK))[chi2_1sigma]
    OmegaLambda_selected_2sigma = (1 - (OmegaM + OmegaK))[chi2_2sigma]

    plt.scatter(OmegaM_selected_2sigma, OmegaLambda_selected_2sigma, label=r'$2\sigma$')
    plt.scatter(OmegaM_selected_1sigma, OmegaLambda_selected_1sigma, label=r'$1\sigma$')
    plt.plot((0,1), (1,0), 'k', ls='--', label='Flat Universe')
    plt.legend()
    plt.show()


def plot_posterior_PDF_Hubble_param():
    # H0 = h / Constants.H0_over_h, see BackgroundCosmology.cpp. PDF will be same.
    posterior_H0_pdf = make_plot()
    posterior_H0_pdf.hist(h, bins=75, density=True)
    bins = posterior_H0_pdf.bins

    sigma = np.std(h)
    mu = np.mean(h)
    gaussian = 1/(sigma*np.sqrt(2*np.pi))*np.exp(-(bins-mu)**2/(2*sigma**2))
    posterior_H0_pdf.plot(bins, gaussian, color='r', lw=2, use_prev_fig=True)

    plt.show()

# Control unit for plotting 
plot_demonstrate_code()
plot_conformal_Hubble()
plot_conformal_time_pr_c()
plot_time()
plot_densities()
plot_luminosity_distance_of_z()
plot_supernovadata_MCMC_fits()
plot_posterior_PDF_Hubble_param()