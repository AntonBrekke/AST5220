import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt 
import scipy.constants as scc
import tabulate as tab
import latextable as lt
import texttable as txttab
import pynverse as pyinv

# Changing standard color cycle to preferred cycle  
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

    def plot(self, x, data, label='', loc='best', **kwargs):
        if (self.x_start) is None: self.x_start = x[0]
        if (self.x_end) is None: self.x_end = x[-1]
        self.x_axis_index = np.logical_and(x >= self.x_start, x <= self.x_end)
        self.ax.plot(x[self.x_axis_index], (data)[self.x_axis_index], label=label, **kwargs)
        self.ax.tick_params(axis='both', which='major', labelsize=14)
        if label != '':
            self.ax.legend(prop={'size': 14}, loc=loc, frameon=False)
        # Call show when ready 

    def hist(self, data, bins, label='', loc='best', density=False):
        self.counts, bins = np.histogram(data, bins=bins)
        self.n, self.bins, self.patches = self.ax.hist(bins[:-1], bins, weights=self.counts, density=density)
        cm = plt.get_cmap('Spectral_r')
        cmap = cm((self.n - np.min(self.n))/(np.max(self.n) - np.min(self.n)))
        for color, patch in zip(cmap, self.patches):
            plt.setp(patch, 'facecolor', color)
        self.ax.tick_params(axis='both', which='major', labelsize=14)
        # # Normalizer 
        norm = mpl.colors.Normalize(vmin=np.min(self.counts), vmax=np.max(self.counts)) 
        # creating ScalarMappable 
        sm = plt.cm.ScalarMappable(cmap=cm, norm=norm) 
        sm.set_array([])
        cbar_ax = self.fig.add_axes([0.64, 0.8, 0.3, 0.03]) 
        cbar = plt.colorbar(sm, cax=cbar_ax, ticks=np.linspace(np.min(self.counts), np.max(self.counts), 6), format=lambda x, pos: f'{x:.0f}', cmap='Spectral_r', orientation='horizontal') 
        cbar.ax.tick_params(labelsize=12)
        if label != '':
            self.ax.legend(prop={'size': 14}, loc=loc, frameon=False)
        # Call show when ready

    def format_plot(self, xlabel='', ylabel='', xscale='linear', yscale='linear', grid=False, tight=True, **scalekwargs):
        self.ax.set_xlabel(xlabel, fontsize=16)
        self.ax.set_ylabel(ylabel, fontsize=16)
        self.ax.set_xscale(xscale, **scalekwargs)
        self.ax.set_yscale(yscale, **scalekwargs)
        self.ax.set_xlim(self.x_start, self.x_end)
        self.ax.grid(grid)
        if tight is True: self.fig.tight_layout()
# End of class 


# Function which finds index of data
def find_index(data, condition, start=None, end=None, eps=0.1, step=1e-5):
    """
    data: array you want index for. 
    condition: abs(F - value). 
    start, end: limit domain of search.
    eps: treshold for acceptance. 
    step: make treshold smaller each iteration.
    return -> index for data when F = value  
    """
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
# End of index-function


# Approximated functions in different domination eras:

def etaHp_R_pr_c(x):
    # Radiation dominated approx
    return np.ones(len(x))

def etaHp_M_pr_c(x):
    # Matter dominated approx
    return np.exp(x_M_dom - 0.5*x) + 2*(1 - np.exp(0.5*(x_M_dom - x))) 

def etaHp_DE_pr_c(x):
    # Dark energy dominated approx
    return np.exp(x + x_M_dom) + np.exp(x - x_Lambda_dom) + 2*(np.exp(x + 0.5*x_Lambda_dom) - np.exp(x + 0.5*x_M_dom)) - 1


# Define plotting functions: 

def plot_demonstrate_code():
    # Define plotting domains for dominations 
    domain_Rel_dom = np.where(x < x_M_dom) 
    domain_M_dom = np.logical_and(x > x_M_dom, x < x_Lambda_dom) 
    domain_DE_dom = np.where(x > x_Lambda_dom) 

    label_dHpdx = r'$\frac{1}{\mathcal{H}(x)}\frac{d\mathcal{H}(x)}{dx}\;$'
    label_ddHpddx = r'$\frac{1}{\mathcal{H}(x)}\frac{d^2\mathcal{H}(x)}{dx^2}\;$'
    title = r'Evolution of derivatives of $\mathcal{H}(x)\equiv aH(x)$'
    plot_Hp_derivatives = make_plot(title=title)
    # Hp'
    plot_Hp_derivatives.plot(x, dHpdx_of_x/Hp_of_x, label=label_dHpdx, lw=2)
    plot_Hp_derivatives.plot(x[domain_Rel_dom], -np.ones(len(x[domain_Rel_dom])), color='tab:red', lw=2)
    plot_Hp_derivatives.plot(x[domain_M_dom], -1/2*np.ones(len(x[domain_M_dom])), color='tab:orange', lw=2)
    plot_Hp_derivatives.plot(x[domain_DE_dom], np.ones(len(x[domain_DE_dom])), color='mediumorchid', lw=2)
    # Hp''
    plot_Hp_derivatives.plot(x, ddHpddx_of_x/Hp_of_x, label=label_ddHpddx, lw=2)
    plot_Hp_derivatives.plot(x[domain_Rel_dom], np.ones(len(x[domain_Rel_dom])), color='tab:red', label='Rad. dom. approx.', lw=2)
    plot_Hp_derivatives.plot(x[domain_M_dom], 1/4*np.ones(len(x[domain_M_dom])), color='tab:orange', label='Matter dom. approx.', lw=2)
    plot_Hp_derivatives.plot(x[domain_DE_dom], np.ones(len(x[domain_DE_dom])), color='mediumorchid', label='DE dom. approx.', lw=2)
    plot_Hp_derivatives.ax.axvline(0, color='k', ls='--')

    plot_Hp_derivatives.format_plot('x=ln(a)', grid=True)
    plot_Hp_derivatives.fig.tight_layout()
    plt.savefig(savefig_path + r'/Hp_derivatives.pdf')
    plt.show()

    # Skip some points for plotting 
    N_points_skip = int(0.045*len(x))
    x_Rel_relevant = x[domain_Rel_dom][::N_points_skip]
    x_M_relevant = x[domain_M_dom][::N_points_skip]
    x_DE_relevant = x[domain_DE_dom][::N_points_skip]

    # Want nice spacing between points. Since matter dominated curve a lot, linear spacing makes it ugly. Distribute points according to function 
    y_M_relevant = np.linspace(etaHp_M_pr_c(x_M_relevant[0]), etaHp_M_pr_c(x_M_relevant[-1]), len(x_M_relevant))
    x_M_relevant = pyinv.inversefunc(etaHp_M_pr_c, y_values=y_M_relevant)

    title = r'$\frac{\eta(x)\mathcal{H}(x)}{c}\;$'
    plot_etaHp_pr_c = make_plot(title=title)
    plot_etaHp_pr_c.plot(x, eta_of_x*Hp_of_x/c, color='k', label='Numerical', lw=2)
    plot_etaHp_pr_c.plot(x_Rel_relevant, etaHp_R_pr_c(x_Rel_relevant), color='tab:red', label='Rad. dom. approx.', ls='--', lw=2, marker='X')
    plot_etaHp_pr_c.plot(x_M_relevant, etaHp_M_pr_c(x_M_relevant), color='tab:orange', label='Matter dom. approx.', ls='--', lw=2, marker='d')
    plot_etaHp_pr_c.plot(x_DE_relevant, etaHp_DE_pr_c(x_DE_relevant), color='tab:green', label='DE dom. approx.', ls='--', lw=2, marker='P')
    plot_etaHp_pr_c.ax.axvline(x_acc, color='k', ls='--')
    plot_etaHp_pr_c.ax.axhline((eta_of_x*Hp_of_x/c)[acc_index], color='k', ls='--')

    plot_etaHp_pr_c.format_plot('x=ln(a)', yscale='log', grid=True)
    plot_etaHp_pr_c.fig.tight_layout()
    plt.savefig(savefig_path + r'/etaHp_pr_c.pdf')
    plt.show()

def plot_conformal_Hubble():
    title = r'Evolution of conformal Hubble factor $\mathcal{H}(x)\;\left(\frac{100km/s}{Mpc}\right)$'
    plot_Hp = make_plot(-12, 1, title=title)
    plot_Hp.plot(x, Hp_of_x / pr_second_to_km_pr_second_pr_Mparsec, lw=2)
    plot_Hp.format_plot('x=ln(a)', yscale='log', grid=True)
    plot_Hp.ax.axvline(x_acc, color='k', ls='--')
    plot_Hp.ax.axhline((Hp_of_x / pr_second_to_km_pr_second_pr_Mparsec)[acc_index], color='k', ls='--')

    plot_Hp.fig.tight_layout()
    plt.savefig(savefig_path + r'/Hp.pdf')
    plt.show()

def plot_conformal_time():
    title = r'Evolution of conformal time $\frac{\eta(x)}{c}\;\left(\frac{Mpc}{c}\right)$'
    plot_eta_pr_c = make_plot(title=title)
    plot_eta_pr_c.plot(x, eta_of_x/(1e6*parsec*c), lw=2)
    plot_eta_pr_c.format_plot('x=ln(a)', yscale='log', grid=True)
    plot_eta_pr_c.ax.axvline(0, color='k', ls='--')
    plot_eta_pr_c.ax.axhline((eta_of_x/(1e6*parsec*c))[today_index], color='k', ls='--')

    plot_eta_pr_c.fig.tight_layout()
    plt.savefig(savefig_path + r'/eta_pr_c.pdf')
    plt.show()

def plot_time():
    plot_t = make_plot(title=r'Evolution of cosmic time $t(a)\;(Gyr)$')
    plot_t.plot(x, t_of_x / Gyr_to_seconds, lw=2)
    plot_t.format_plot('x=ln(a)', yscale='log', grid=True)
    plot_t.ax.axvline(0, color='k', ls='--')
    plot_t.ax.axhline((t_of_x / Gyr_to_seconds)[today_index], color='k', ls='--')

    plot_t.fig.tight_layout()
    plt.savefig(savefig_path + r'/cosmic_time.pdf')
    plt.show()

def plot_densities():
    plot_densities = make_plot(title=r'$\Omega_i$')
    plot_densities.plot(x, OmegaRel, label=r'$\Omega_{\text{R}}$', lw=2)
    plot_densities.plot(x, OmegaM, label=r'$\Omega_{\text{M}}$', lw=2)
    plot_densities.plot(x, OmegaLambda, label=r'$\Omega_\Lambda$', lw=2, loc='center left')

    # Find density parameter equalities 
    Rel_M_equality = abs(OmegaRel - OmegaM)
    M_Lambda_equality = abs(OmegaLambda - OmegaM)

    index_RelM_eq = find_index(x, Rel_M_equality, start=-10, end=-5)
    index_MLambda_eq = find_index(x, M_Lambda_equality, start=-1, end=1)
    x_RelM_eq = x[index_RelM_eq]
    x_MLambda_eq = x[index_MLambda_eq]

    plot_densities.ax.axvline(x_RelM_eq, color='k', ls='--')
    plot_densities.ax.axvline(x_MLambda_eq, color='k', ls='--')

    plot_densities.format_plot('x=ln(a)', grid=True)

    plot_densities.fig.tight_layout()
    plt.savefig(savefig_path + r'/densities.pdf')
    plt.show()

def plot_luminosity_distance_of_z():
    fig = plt.figure()
    ax = fig.add_subplot()

    ax.set_title(r'Luminosity distance $d_L / z$', fontsize=16)
    ax.errorbar(z, luminosity_distance_of_z / z, ls='', marker='o', ms=3, yerr=luminosity_error/z, ecolor='k', capsize=0, color="k", label=r'Supernovadata')
    ax.plot(z_of_x, luminosity_distance_of_x / (z_of_x*parsec*1e9), color='r', label='Theoretical', lw=2)

    ax.set_xlim(1.1*z[0], 1.1*z[-1])
    ax.set_ylim(4, 8)
    ax.set_xscale('log')
    ax.set_xlabel('Redshift z', fontsize=16)
    ax.set_ylabel(r'$d_L$/z (Gpc)', fontsize=16)
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.grid(True)

    ax.legend(prop={'size': 14}, frameon=False)
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
    # ax.scatter(OmegaM_selected_2sigma, OmegaLambda_selected_2sigma, c=chi2_2sigma, cmap='winter', label=r'$2\sigma$', s=10, alpha=0.5, rasterized=True)
    sigma1_plot = ax.scatter(OmegaM_selected_1sigma, OmegaLambda_selected_1sigma, c=chi2[chi2_1sigma], cmap='RdYlBu_r', label=r'$1\sigma$', s=10, alpha=0.9, rasterized=True)
    ax.plot((0,1), (1,0), 'k', ls='--', label='Flat Universe')

    cbar_ax = fig.add_axes([0.53, 0.25, 0.4, 0.03]) 
    cbar = plt.colorbar(sigma1_plot, cax=cbar_ax, orientation='horizontal') 
    cbar.ax.set_title(r'$\chi^2$', fontsize=16)
    cbar.ax.tick_params(labelsize=12)

    ax.set_xlim(0, 0.5)
    ax.set_ylim(0.3, 1)
    ax.set_xlabel(r'$\Omega_{M}$', fontsize=16)
    ax.set_ylabel(r'$\Omega_{\Lambda}$', fontsize=16)
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.set_title('Supernovadata MCMC fits', fontsize=16)

    ax.legend(prop={'size':14}, loc='center left', bbox_to_anchor=(0, 0.7), frameon=False)
    fig.tight_layout()
    plt.savefig(savefig_path + r'/supernovadata_MCMC_fits.pdf')
    plt.show()

def plot_posterior_PDF_Hubble_param():
    # H0 = h / Constants.H0_over_h, see BackgroundCosmology.cpp. PDF will be same.
    H0 = h * H0_over_h
    posterior_H0_pdf = make_plot(title=r'Posterior PDF for $H_0$')
    posterior_H0_pdf.hist(H0, bins=75, density=True)
    bins = posterior_H0_pdf.bins
    n = posterior_H0_pdf.n
    # print(np.sum(np.diff(posterior_H0_pdf.bins) * posterior_H0_pdf.n))   # Testing if prop. dist sum to 1 

    sigma = np.std(H0)
    mu = np.mean(H0)
    gaussian = 1/(sigma*np.sqrt(2*np.pi))*np.exp(-(bins-mu)**2/(2*sigma**2))
    posterior_H0_pdf.plot(bins, gaussian, color='k', lw=2.5)
    posterior_H0_pdf.ax.axvline(H0_best_fit, color='k', ls='--', lw=2, label=r'Best fit value of $H_0$')
    posterior_H0_pdf.ax.axvline(mu, color='k', ls='-', lw=2)
    posterior_H0_pdf.ax.legend(prop={'size':14}, frameon=False)
    # print(mu, H0_best_fit)

    posterior_H0_pdf.format_plot(xlabel=r'$H_0$', tight=True)
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

    # Make tables
    tab_print = [['Event', '$x$', '$z$', '$t$ (Gyr)'],
                 ['Radiation-Matter equality', f'{x_RelM_eq:.3f}', f'{z_RelM_eq:.3f}', f'{t_RelM_eq:.3e}'],
                 ['Matter-Dark Energy equality', f'{x_MLambda_eq:.3f}', f'{z_MLambda_eq:.3f}', f'{t_MLambda_eq:.3f}'],
                 ['Universe accelerate', f'{x_acc:.3f}', f'{z_acc:.3f}', f'{t_acc:.3f}'],
                 ['Today', f'{x_age:.3f}', f'{z_age:.3f}', f'{t_age:.3f}']]

    tab_data = tab.tabulate(tab_print, tablefmt="simple_grid")
    print(tab_data)
    print(f'Conformal time today: {(eta_of_x/(c*Gyr_to_seconds))[age_index]:.3f} (Gyr)\n')
    
    if latex is True:
        # Latex table
        print('Beautiful LateX-table, just for you! (double check small values)')
        print('================================================================\n')
        tab_data_latex = txttab.Texttable()
        tab_data_latex.set_cols_align(["l", "c", "c", "c"])
        tab_data_latex.add_rows(tab_print)
        print(lt.draw_latex(tab_data_latex))
# End of plotting functions 


# Scaling factor instead of x may be more intuitive
a = np.exp(x)    
h_best_fit = 0.70189
H0_over_h = 100
H0_best_fit = H0_over_h * h_best_fit
today_index = find_index(x, abs(x))

# # Define shorthand for these sums 
OmegaRel = OmegaR + OmegaNu
OmegaM = OmegaB + OmegaCDM

# Find density parameters domination. Defined to be dominating if 50% is quantity. 
percent_for_domination = 0.5
M_dom = abs(OmegaM - percent_for_domination)
Lambda_dom = abs(OmegaLambda - percent_for_domination)

index_M_dom = find_index(x, M_dom, start=-10, end=-5)
index_Lambda_dom = find_index(x, Lambda_dom, start=-1, end=1)

# These are needed for plotting
x_M_dom = x[index_M_dom]
x_Lambda_dom = x[index_Lambda_dom]

# Universe accelerating (when a_double_dot is positive)
a_double_dot = np.exp(-x)*dHpdx_of_x*Hp_of_x
acc_index = np.min(np.where(a_double_dot >= 0))
# acc_index = find_index(x, abs(a_double_dot/np.min(a_double_dot)), step=1e-10)
x_acc  = x[acc_index]
z_acc = z_of_x[acc_index]
t_acc = (t_of_x/Gyr_to_seconds)[acc_index]

print(OmegaLambda[acc_index])
print(OmegaR[acc_index])
print(OmegaM[acc_index])

# Control unit for plotting 
# plot_demonstrate_code()

# plot_conformal_Hubble()
# plot_conformal_time()
# plot_time()
# plot_densities()
# plot_luminosity_distance_of_z()
# plot_supernovadata_MCMC_fits()
# plot_posterior_PDF_Hubble_param()
# plot_posterior_PDF_OmegaLambda()

# make_table(latex=False)

