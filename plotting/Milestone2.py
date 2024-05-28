import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as scc
import tabulate as tab
import latextable as lt
import texttable as txttab

# Relevant paths to get data
path = r'/home/antonabr/AST5220/data'
savefig_path = r'/home/antonabr/AST5220/figures/Milestone2'

# Changing standard color cycle to preferred cycle  
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
colors.insert(6, colors[1]) ; colors[1:5] = colors[2:5]
plt.rcParams['axes.prop_cycle'] = plt.cycler(color=colors) 

# Get data
recombination_data = np.loadtxt(fr'{path}' + r'/recombination.txt')

# Assign data 
x, t_of_x, z_of_x, Xe_of_x, XeSaha_of_x, ne_of_x, tau_of_x, dtaudx_of_x, \
ddtauddx_of_x, g_tilde_of_x, dgdx_tilde_of_x, ddgddx_tilde_of_x, T_of_x, s_of_x = recombination_data.T

# Practical constants
c = scc.c       # m/s
parsec = scc.parsec     # 1 parsec in m
Myr_to_seconds = 365*24*60*60*1e6
m_to_Mpc = 1 / (1e6*parsec)             

# Less scuffed version of last index-finder (Milestone 1)
def find_index(data, x=None, value=0, x_start=None, x_end=None):
    if x is None:
        index = np.where(np.min(abs(data - value)) == abs(data - value))[0][0]
        return index
    else: # Specify what data is function of - used if data = value at multiple x-values
        if x_start is None: x_start = np.min(x)
        if x_end is None: x_end = np.max(x)
        domain = (x >= x_start) * (x <= x_end)
        data_domain = data[domain]
        index = np.where(np.min(abs(data_domain - value)) == abs(data_domain - value))[0][0]
        index = len(x[(x <= x_start)]) + index      # Translate to index on original domain
        return index
# End of index-function


def plot_fractional_electron_density(show=True):
    # Plot fractional electron density Xe
    fig = plt.figure()
    ax = fig.add_subplot()
    ax.set_title(r'Evolution of fractional electron density $X_e$', fontsize=16)

    # Plotting 
    ax.plot(x, Xe_of_x, label=r'$X_e$', lw=2)
    ax.plot(x, XeSaha_of_x, label=r'$X_e$ (Saha)', lw=2, ls='--')
    ax.axvline(x_recombination, color='k', label='Recombination', ls='--')
    ax.axhline(Xe_freeze_out, color='purple', label='Freeze-out', ls='-.')

    # Formatting 
    ax.set_xlim(-8, -4)
    ax.set_ylim(1e-4, 1.25)
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.set_xlabel('x=ln(a)', fontsize=16)
    ax.set_yscale('log')
    ax.grid(True)
    ax.legend(prop={'size': 14}, loc='best', frameon=False)
    fig.tight_layout()

    # Show 
    plt.savefig(savefig_path + r'/frac_electron_density.pdf')
    if show is True: plt.show()


def plot_optical_depth(show=True):
    # Plot optical depth tau
    fig = plt.figure()
    ax = fig.add_subplot()
    ax.set_title(r'Evolution of optical depth $\tau$', fontsize=16)

    # Plotting 
    ax.plot(x, tau_of_x, label=r'$\tau(x)$', lw=2)
    ax.plot(x, -dtaudx_of_x, label=r'$-\tau^\prime(x)$', lw=2)
    ax.plot(x, ddtauddx_of_x, label=r'$\tau^{\prime\prime}(x)$', lw=2)
    ax.axvline(x_last_scatter, color='k', label='Last scattering', ls='--')

    # Formatting
    ax.set_xlim(-10, -4)
    ax.set_ylim(1e-5, 1e5)
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.set_xlabel('x=ln(a)', fontsize=16)
    ax.set_yscale('log')
    ax.grid(True)
    ax.legend(prop={'size': 14}, loc='best', frameon=False)
    fig.tight_layout()

    # Show 
    plt.savefig(savefig_path + r'/optical_depth.pdf')
    if show is True: plt.show()

def plot_visibility_function(show=True):
    # Plot visibility function g
    fig = plt.figure()
    ax = fig.add_subplot()
    ax.set_title(r'Evolution of visibility function $\tilde g$', fontsize=16)

    # Normalizing each graph so they can be plottet together
    g_tilde_max = np.max(abs(g_tilde_of_x))
    dgdx_tilde_max = np.max(abs(dgdx_tilde_of_x))
    ddgddx_tilde_max = np.max(abs(ddgddx_tilde_of_x))

    # Plotting 
    ax.plot(x, g_tilde_of_x, label=r'$\tilde g(x)$', lw=2)
    ax.plot(x, dgdx_tilde_of_x * g_tilde_max / dgdx_tilde_max, label=r'$\tilde g^\prime(x) \cdot |\tilde g|_{\text{max}} / |\tilde g^\prime|_{\text{max}}$', lw=2)
    ax.plot(x, ddgddx_tilde_of_x * g_tilde_max / ddgddx_tilde_max, label=r'$\tilde g^{\prime\prime}(x) \cdot |\tilde g|_{\text{max}} / |\tilde g^{\prime\prime}|_{\text{max}}$', lw=2)
    ax.axvline(x_last_scatter, color='k', label='Last scattering', ls='--')

    # Formatting 
    ax.set_xlim(-8, -6)
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.set_xlabel('x=ln(a)', fontsize=16)
    ax.grid(True)
    ax.legend(prop={'size': 14}, loc='best', frameon=False)
    fig.tight_layout()

    # Show
    plt.savefig(savefig_path + r'/visibility_function.pdf')
    if show is True: plt.show()
# End of plotting functions 

def find_values(latex=False):
    # Definition of recombination
    Xe_recombination = 0.1
    recombination_index = find_index(Xe_of_x, value=Xe_recombination)
    # print(Xe_of_x[recombination_index])     # Check

    x_recombination = x[recombination_index]
    t_recombination = (t_of_x / Myr_to_seconds)[recombination_index]
    z_recombination = z_of_x[recombination_index]
    s_recombination = s_of_x[recombination_index] * m_to_Mpc

    # Definition of last scattering 
    last_scatter_index = find_index(g_tilde_of_x, value=np.max(g_tilde_of_x))

    x_last_scatter = x[last_scatter_index]
    t_last_scatter = (t_of_x / Myr_to_seconds)[last_scatter_index]
    z_last_scatter = z_of_x[last_scatter_index]
    s_last_scatter = s_of_x[last_scatter_index] * m_to_Mpc

    # Definition of freeze-out
    freeze_out_index = find_index(x, value=0)
    Xe_freeze_out = Xe_of_x[freeze_out_index]
    print(tau_of_x[last_scatter_index])
    print(f'Photon temperature recombination: {T_of_x[recombination_index]}K')
    print(f'Photon energy recombination: {scc.physical_constants["Boltzmann constant in eV/K"][0]*T_of_x[recombination_index]}eV')
    # print(x[freeze_out_index])      # Check

    # Make tables
    p = 5
    tab_print = [['Event', '$x$', '$z$', '$t$ (Myr)', '$r_s$ (Mpc)'],
                 ['Recombination', f'{x_recombination:.{p}f}', f'{z_recombination:.{p}f}', f'{t_recombination:.{p}f}', f'{s_recombination:.{p}f}'],
                 ['Last scatter', f'{x_last_scatter:.{p}f}', f'{z_last_scatter:.{p}f}', f'{t_last_scatter:.{p}f}', f'{s_last_scatter:.{p}f}']]

    tab_data = tab.tabulate(tab_print, tablefmt="simple_grid")
    print(tab_data)
    print(f'Freeze-out abundance free electrons: {Xe_freeze_out:.{p}e}\n')
    
    if latex is True:
        # Latex table
        print('Beautiful LateX-table, just for you! (double check small values)')
        print('================================================================\n')
        tab_data_latex = txttab.Texttable()
        tab_data_latex.set_cols_align(["l", "c", "c", "c", "c"])
        tab_data_latex.add_rows(tab_print)
        print(lt.draw_latex(tab_data_latex))

    return [x_recombination, x_last_scatter, Xe_freeze_out]
# End of find values 


x_recombination, x_last_scatter, Xe_freeze_out = find_values(latex=False)

if __name__ == "__main__":
    # Control unit for calling functions
    plot_fractional_electron_density(show=True)
    plot_optical_depth(show=True)
    plot_visibility_function(show=True)
