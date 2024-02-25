import numpy as np
import matplotlib.pyplot as plt 
import scipy.constants as scc
import tabulate as tab
import latextable as lt
import texttable as txt

# Changing standardcolor cycle to preferred cycle  
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
colors.insert(6, colors[1]) ; colors[1:5] = colors[2:5]
plt.rcParams['axes.prop_cycle'] = plt.cycler(color=colors) 

# Relevant paths to get data
path = r'/home/antonabr/AST5220/data'
savefig_path = r'/home/antonabr/AST5220/figures/Milestone1'

# Get data
supernovadata = np.loadtxt(fr'{path}' + r'/supernovadata.txt')
results_supernovafitting = np.loadtxt(fr'{path}' + r'/results_supernovafitting.txt', skiprows=200)
cosmology = np.loadtxt(fr'{path}' + r'/cosmology.txt')

# Assign data
x, eta_of_x, detadx_of_x, t_of_x, Hp_of_x, dHpdx_of_x, ddHpddx_of_x, OmegaB, OmegaCDM, \
OmegaLambda, OmegaR, OmegaNu, OmegaK, luminosity_distance_of_x, z_of_x = cosmology.T

chi2, h, OmegaM_sn, OmegaK_sn = results_supernovafitting.T

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

    def format_plot(self, xlabel='', ylabel='', xscale='linear', yscale='linear', **scalekwargs):
        self.ax.set_xlabel(xlabel, fontsize=16)
        self.ax.set_ylabel(ylabel, fontsize=16)
        self.ax.set_xscale(xscale, **scalekwargs)
        self.ax.set_yscale(yscale, **scalekwargs)
        self.ax.set_xlim(self.x_start, self.x_end)
        self.fig.tight_layout()
# End of class 


# Function which finds index of data for given abs(F - value) < 0. Made this when Newtons method failed me 
def find_index(data, condition, start=None, end=None, eps=0.1, step=1e-5):
    if start is None: start = np.min(data)
    if end is None: end = np.max(data)
    domain = (data > start) * (data < end)
    condition = condition[domain]
    index_satisfied_condition = np.where(condition < eps)[0]
    N_satisfied_condition = len(index_satisfied_condition)

    while N_satisfied_condition > 1:
        eps -= step 
        index_satisfied_condition = np.where(condition < eps)[0]
        N_satisfied_condition = len(index_satisfied_condition)
    if N_satisfied_condition < 1: raise Exception("Could not find single value.")
    
    index = np.where(data[domain][index_satisfied_condition] == data)[0]
    return index[0]


# Approximated functions in different domination eras 

def etaHp_R_pr_c(x):
    # Radiation dominated approx
    return np.ones(len(x))

def etaHp_M_pr_c(x):
    # Matter dominated approx
    return np.exp(x_RelM_dom - 0.5*x) + 2*(1 - np.exp(0.5*(x_RelM_dom - x))) 

def etaHp_L_pr_c(x):
    # Dark energy dominated approx
    return np.exp(x + x_RelM_dom) + np.exp(x - x_MLambda_dom) + 2*(np.exp(x + 0.5*x_MLambda_dom) - np.exp(x + 0.5*x_RelM_dom)) - 1


# Define plotting functions: 

def plot_demonstrate_code():
    label_dHpdx = r'$\frac{1}{\mathcal{H}(x)}\frac{d\mathcal{H}(x)}{dx}\;$'
    label_ddHpddx = r'$\frac{1}{\mathcal{H}(x)}\frac{d^2\mathcal{H}(x)}{dx^2}\;$'
    title = r'Evolution of derivatives of $\mathcal{H}(x)\equiv aH(x)$'
    plot_Hp_derivatives = make_plot(title=title)
    plot_Hp_derivatives.plot(x, dHpdx_of_x/Hp_of_x, label=label_dHpdx)
    plot_Hp_derivatives.plot(x, ddHpddx_of_x/Hp_of_x, label=label_ddHpddx)
    plot_Hp_derivatives.plot(x, np.ones(len(x)), color='k', ls=':')
    plot_Hp_derivatives.plot(x, -np.ones(len(x)), color='k', ls=':')
    plot_Hp_derivatives.plot(x, -1/2*np.ones(len(x)), color='k', ls=':')
    plot_Hp_derivatives.plot(x, 1/4*np.ones(len(x)), color='k', ls=':')

    plot_Hp_derivatives.format_plot('x=ln(a)')
    plot_Hp_derivatives.fig.tight_layout()
    plt.savefig(savefig_path + r'/Hp_derivatives.pdf')
    plt.show()

    # Define plotting domains for dominations 
    domain_Rel_dom = np.where(x < x_RelM_dom) 
    domain_M_dom = np.logical_and(x > x_RelM_dom, x < x_MLambda_dom) 
    domain_DE_dom = np.where(x > x_MLambda_dom) 

    title = r'$\frac{\eta(x)\mathcal{H}(x)}{c}\;$'
    plot_etaHp_pr_c = make_plot(title=title)
    plot_etaHp_pr_c.plot(x, eta_of_x*Hp_of_x/c, color='k', label='Numerical', lw=2, ls='--')
    plot_etaHp_pr_c.plot(x[domain_Rel_dom], etaHp_R_pr_c(x[domain_Rel_dom]), color='tab:red', label='Rad. dom. approx.')
    plot_etaHp_pr_c.plot(x[domain_M_dom], etaHp_M_pr_c(x[domain_M_dom]), color='tab:orange', label='Matter dom. approx.')
    plot_etaHp_pr_c.plot(x[domain_DE_dom], etaHp_L_pr_c(x[domain_DE_dom]), color='tab:green', label='DE dom. approx.')

    plot_etaHp_pr_c.format_plot('x=ln(a)', yscale='log')
    plot_etaHp_pr_c.fig.tight_layout()
    plt.savefig(savefig_path + r'/etaHp_pr_c.pdf')
    plt.show()

def plot_conformal_Hubble():
    title = r'Evolution of conformal Hubble factor $\mathcal{H}(x)\;\left(\frac{100km/s}{Mpc}\right)$'
    plot_Hp = make_plot(-12, 0, title=title)
    plot_Hp.plot(x, Hp_of_x / pr_second_to_km_pr_second_pr_Mparsec)
    plot_Hp.format_plot('x=ln(a)', yscale='log')

    plot_Hp.fig.tight_layout()
    plt.savefig(savefig_path + r'/Hp.pdf')
    plt.show()

def plot_conformal_time_pr_c():
    title = r'Evolution of conformal time $\frac{\eta(x)}{c}\;\left(\frac{Mpc}{c}\right)$'
    plot_eta_pr_c = make_plot(title=title)
    plot_eta_pr_c.plot(x, eta_of_x/(1e6*parsec*c))
    plot_eta_pr_c.format_plot('x=ln(a)', yscale='log')

    plot_eta_pr_c.fig.tight_layout()
    plt.savefig(savefig_path + r'/eta_pr_c.pdf')
    plt.show()

def plot_time():
    plot_t = make_plot(title=r'Evolution of cosmic time $t(a)\;(Gyr)$')
    plot_t.plot(x, t_of_x / Gyr_to_seconds)
    plot_t.format_plot('x', yscale='log')

    plot_t.fig.tight_layout()
    plt.savefig(savefig_path + r'/cosmic_time.pdf')
    plt.show()

def plot_densities():
    plot_densities = make_plot(x[0], x[-1], title=r'$\Omega_i$')
    plot_densities.plot(x, OmegaRel, label=r'$\Omega_{rel}=\Omega_R + \Omega_\nu$')
    plot_densities.plot(x, OmegaM, label=r'$\Omega_{matter}=\Omega_B + \Omega_{CDM}$')
    plot_densities.plot(x, OmegaLambda, label=r'$\Omega_\Lambda$')

    # Find density parameter equalities 
    Rel_M_equality = abs(OmegaRel - OmegaM)
    M_Lambda_equality = abs(OmegaLambda - OmegaM)

    index_RelM_eq = find_index(x, Rel_M_equality, start=-10, end=-5)
    index_MLambda_eq = find_index(x, M_Lambda_equality, start=-1, end=1)
    x_RelM_eq = x[index_RelM_eq]
    x_MLambda_eq = x[index_MLambda_eq]

    plot_densities.ax.vlines(x_RelM_eq, 0, 1, color='k', ls='--')
    plot_densities.ax.vlines(x_MLambda_eq, 0, 1, color='k', ls='--')

    plot_densities.format_plot('x=ln(a)')

    plot_densities.fig.tight_layout()
    plt.savefig(savefig_path + r'/densities.pdf')
    plt.show()

def plot_luminosity_distance_of_z():
    fig = plt.figure()
    ax = fig.add_subplot()
    ax.set_title(r'Luminosity distance $d_L / z$ compared to supernovadata', fontsize=16)
    ax.errorbar(z, luminosity_distance_of_z / z, ls='', marker='o', ms=3, yerr=luminosity_error/z, ecolor='k', capsize=0, color="k", label=r'Supernovadata')
    ax.plot(z_of_x, luminosity_distance_of_x / (z_of_x*parsec*1e9), color='r', label='Theoretical')
    ax.set_xlim(1.1*z[0], 1.1*z[-1])
    ax.set_ylim(4, 8)
    ax.set_xscale('log')
    ax.set_xlabel('Redshift z', fontsize=16)
    ax.set_ylabel(r'$d_L$/z (Mpc)', fontsize=16)

    ax.legend(prop={'size': 12})
    fig.tight_layout()
    plt.savefig(savefig_path + r'/luminosity_distance.pdf')
    plt.show()

def plot_supernovadata_MCMC_fits():
    fig = plt.figure()
    ax = fig.add_subplot()

    chi2_min = np.min(chi2)

    chi2_1sigma = np.where((chi2 - chi2_min) < 3.53)    # 1 sigma
    chi2_2sigma = np.where((chi2 - chi2_min) < 8.02)    # 2 sigma

    OmegaM_selected_1sigma = OmegaM_sn[chi2_1sigma]
    OmegaM_selected_2sigma = OmegaM_sn[chi2_2sigma]
    OmegaLambda_selected_1sigma = (1 - (OmegaM_sn + OmegaK_sn))[chi2_1sigma]
    OmegaLambda_selected_2sigma = (1 - (OmegaM_sn + OmegaK_sn))[chi2_2sigma]

    line = -1*np.linspace(0, 1, len(OmegaLambda_selected_2sigma)) + 1
    ax.scatter(OmegaM_selected_2sigma, OmegaLambda_selected_2sigma, label=r'$2\sigma$')
    ax.scatter(OmegaM_selected_1sigma, OmegaLambda_selected_1sigma, label=r'$1\sigma$')
    ax.plot((0,1), (1,0), 'k', ls='--', label='Flat Universe')

    ax.set_xlabel(r'$\Omega_{M}$', fontsize=16)
    ax.set_ylabel(r'$\Omega_{\Lambda}$', fontsize=16)
    ax.legend(prop={'size':12})
    fig.tight_layout()
    plt.savefig(savefig_path + r'/supernovadata_MCMC_fits.pdf')
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

    posterior_H0_pdf.fig.tight_layout()
    plt.savefig(savefig_path + r'/posterior_PDF_Hubble_param.pdf')
    plt.show()

def plot_posterior_PDF_OmegaLambda():
    OmegaLambda_sn = 1 - (OmegaM_sn + OmegaK_sn)
    posterior_OmegaLambda_pdf = make_plot(title=r'Posterior PDF for $\Omega_{\Lambda}$')
    posterior_OmegaLambda_pdf.hist(OmegaLambda_sn, bins=75, density=True)
    bins = posterior_OmegaLambda_pdf.bins
    # print(np.sum(np.diff(posterior_OmegaLambda_pdf.bins) * posterior_OmegaLambda_pdf.n))   # Testing if prop. dist sum to 1 

    sigma = np.std(OmegaLambda_sn)
    mu = np.mean(OmegaLambda_sn)
    gaussian = 1/(sigma*np.sqrt(2*np.pi))*np.exp(-(bins-mu)**2/(2*sigma**2))
    posterior_OmegaLambda_pdf.plot(bins, gaussian, color='k', lw=2.5)
    posterior_OmegaLambda_pdf.format_plot(xlabel=r'$\Omega_{\Lambda}$')

    posterior_OmegaLambda_pdf.fig.tight_layout()
    plt.savefig(savefig_path + r'/posterior_PDF_OmegaLambda.pdf')
    plt.show()

def make_table(latex=False):
    # Find density parameters equalities 
    Rel_M_equality = abs(OmegaRel - OmegaM)
    M_Lambda_equality = abs(OmegaLambda - OmegaM)

    index_RelM_eq = find_index(x, Rel_M_equality, start=-10, end=-5)
    index_MLambda_eq = find_index(x, M_Lambda_equality, start=-1, end=1)
    x_RelM_eq = x[index_RelM_eq]
    z_RelM_eq = z_of_x[index_RelM_eq]
    t_RelM_eq = (t_of_x/Gyr_to_seconds)[index_RelM_eq]

    x_MLambda_eq = x[index_MLambda_eq]
    z_MLambda_eq = z_of_x[index_MLambda_eq]
    t_MLambda_eq = (t_of_x/Gyr_to_seconds)[index_MLambda_eq]

    # Find age of universe today (x=0)
    age_index = find_index(t_of_x/Gyr_to_seconds, abs(x))
    x_age = x[age_index]
    z_age = z_of_x[age_index]
    t_age = (t_of_x/Gyr_to_seconds)[age_index]

    # Universe accelerating (when a_double_dot is positive)
    a_double_dot = np.exp(-x)*dHpdx_of_x*Hp_of_x
    acc_index = np.min(np.where(a_double_dot >= 0))
    # acc_index = find_index(x, abs(a_double_dot/np.min(a_double_dot)), step=1e-10)
    x_acc  = x[acc_index]
    z_acc = z_of_x[acc_index]
    t_acc = (t_of_x/Gyr_to_seconds)[acc_index]

    # Make tables
    tab_print = [['', 'x', 'z', 't (Gyr)'],
                 ['Rad. Mat. eq', f'{x_RelM_eq:.3f}', f'{z_RelM_eq:.3f}', f'{t_RelM_eq:.3e}'],
                 ['Mat. DE. eq', f'{x_MLambda_eq:.3f}', f'{z_MLambda_eq:.3f}', f'{t_MLambda_eq:.3f}'],
                 ['Uni. today', f'{x_age:.3f}', f'{z_age:.3f}', f'{t_age:.3f}'],
                 ['Uni. acc', f'{x_acc:.3f}', f'{z_acc:.3f}', f'{t_acc:.3f}']]

    tab_data = tab.tabulate(tab_print, tablefmt="simple_grid")
    print(tab_data)
    print(f'Conformal time today: {(eta_of_x/(c*Gyr_to_seconds))[age_index]:.3f} (Gyr)\n')
    
    if latex is True:
        # Latex table
        print('Beautiful LateX-table, just for you! (double check small values)')
        print('================================================================\n')
        tab_data_latex = txt.Texttable()
        tab_data_latex.set_cols_align(["l", "c", "c", "c"])
        tab_data_latex.add_rows(tab_print)
        print(lt.draw_latex(tab_data_latex))



# Scaling factor instead of x may be more intuitive
a = np.exp(x)       

# # Define shorthand for these sums 
OmegaRel = OmegaR + OmegaNu
OmegaM = OmegaB + OmegaCDM

# Find density parameters domination. Defined to be dominating if 0.75% is left. 
percent_for_domination = 0.75
Rel_dom = abs(OmegaRel - percent_for_domination)
M_dom = abs(OmegaM - percent_for_domination)

index_Rel_dom = find_index(x, Rel_dom, start=-10, end=-5)
index_M_dom = find_index(x, M_dom, start=-1, end=1)

# These are needed for plotting
x_RelM_dom = x[index_Rel_dom]
x_MLambda_dom = x[index_M_dom]

# Control unit for plotting 
plot_demonstrate_code()
plot_conformal_Hubble()
plot_conformal_time_pr_c()
plot_time()
plot_densities()
plot_luminosity_distance_of_z()
plot_supernovadata_MCMC_fits()
plot_posterior_PDF_Hubble_param()
plot_posterior_PDF_OmegaLambda()

# make_table(latex=True)

