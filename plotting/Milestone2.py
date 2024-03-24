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
ddtauddx_of_x, g_tilde_of_x, dgdx_tilde_of_x, ddgddx_tilde_of_x, s_of_x = recombination_data.T

# Practical constants
c = scc.c       # m/s
parsec = scc.parsec     # 1 parsec in m
Myr_to_seconds = 365*24*60*60*1e6
m_to_Mpc = 1 / (1e6*parsec)             

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

def plot_fractional_electron_density():
    # Plot fractional electron density Xe
    fig = plt.figure()
    ax = fig.add_subplot()
    ax.set_title(r'Evolution of fractional electron density $X_e$', fontsize=16)

    # Plotting 
    ax.plot(x, Xe_of_x, label=r'$X_e$', lw=2)
    ax.plot(x, XeSaha_of_x*(XeSaha_of_x > np.min(Xe_of_x)), label=r'$X_e$ (Saha)', lw=2, ls='--')

    # Formatting 
    ax.set_xlim(-8, -4)
    ax.set_ylim(1.7e-4, 2)
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.set_xlabel('x=ln(a)', fontsize=16)
    ax.set_yscale('log')
    ax.grid(True)
    ax.legend(prop={'size': 14}, loc='best', frameon=False)
    fig.tight_layout()

    # Show 
    plt.savefig(savefig_path + r'/frac_electron_density.pdf')
    plt.show()


def plot_optical_depth():
    # Plot optical depth tau
    fig = plt.figure()
    ax = fig.add_subplot()
    ax.set_title(r'Evolution of optical depth $\tau$', fontsize=16)

    # Plotting 
    ax.plot(x, tau_of_x, label=r'$\tau(x)$', lw=2)
    ax.plot(x, -dtaudx_of_x, label=r'$-\tau^\prime(x)$', lw=2)
    ax.plot(x, ddtauddx_of_x, label=r'$\tau^{\prime\prime}(x)$', lw=2)

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
    plt.show()

def plot_visibility_function():
    # Plot visibility function g
    fig = plt.figure()
    ax = fig.add_subplot()
    ax.set_title(r'Evolution of visibility function $\tilde g$', fontsize=16)

    # Normalizing each graph so they can be plottet together
    g_tilde_max = np.max(abs(g_tilde_of_x))
    dgdx_tilde_max = np.max(abs(dgdx_tilde_of_x))
    ddgddx_tilde_max = np.max(abs(ddgddx_tilde_of_x))

    # Plotting 
    ax.plot(x, g_tilde_of_x / g_tilde_max, label=r'$\tilde g(x) / |\tilde g|_{\text{max}}$', lw=2)
    ax.plot(x, dgdx_tilde_of_x / dgdx_tilde_max, label=r'$\tilde g^\prime(x) / |\tilde g^\prime|_{\text{max}}$', lw=2)
    ax.plot(x, ddgddx_tilde_of_x / ddgddx_tilde_max, label=r'$\tilde g^{\prime\prime}(x) / |\tilde g^{\prime\prime}|_{\text{max}}$', lw=2)

    # Formatting 
    ax.set_xlim(-8, -6)
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.set_xlabel('x=ln(a)', fontsize=16)
    ax.grid(True)
    ax.legend(prop={'size': 14}, loc='best', frameon=False)
    fig.tight_layout()

    # Show
    plt.savefig(savefig_path + r'/visibility_function.pdf')
    plt.show()
# End of plotting functions 

def find_values(latex=False):
    # Definition of recombination
    Xe_recombination = 0.1
    recombination_index = find_index(x, abs(Xe_of_x - Xe_recombination))

    x_recombination = x[recombination_index]
    t_recombination = (t_of_x / Myr_to_seconds)[recombination_index]
    z_recombination = z_of_x[recombination_index]
    s_recombination = s_of_x[recombination_index] * m_to_Mpc

    # Definition of last scattering 
    last_scatter_index = np.where(g_tilde_of_x == np.max(g_tilde_of_x))[0][0]

    x_last_scatter = x[last_scatter_index]
    t_last_scatter = (t_of_x / Myr_to_seconds)[last_scatter_index]
    z_last_scatter = z_of_x[last_scatter_index]
    s_last_scatter = s_of_x[last_scatter_index] * m_to_Mpc

    # Definition of freeze-out
    freeze_out_index = find_index(Xe_of_x, abs(x))
    Xe_freeze_out = Xe_of_x[freeze_out_index]


    # Make tables
    tab_print = [['Event', '$x$', '$z$', '$t$ (Myr)', '$r_s$ (Mpc)'],
                 ['Recombination', f'{x_recombination:.4f}', f'{z_recombination:.4f}', f'{t_recombination:.4f}', f'{s_recombination:.4f}'],
                 ['Last scatter', f'{x_last_scatter:.4f}', f'{z_last_scatter:.4f}', f'{t_last_scatter:.4f}', f'{s_last_scatter:.4f}']]

    tab_data = tab.tabulate(tab_print, tablefmt="simple_grid")
    print(tab_data)
    print(f'Freeze-out abundance free electrons: {Xe_freeze_out:.3e}\n')
    
    if latex is True:
        # Latex table
        print('Beautiful LateX-table, just for you! (double check small values)')
        print('================================================================\n')
        tab_data_latex = txttab.Texttable()
        tab_data_latex.set_cols_align(["l", "c", "c", "c"])
        tab_data_latex.add_rows(tab_print)
        print(lt.draw_latex(tab_data_latex))
# End of find values 

# Control unit for calling functions
plot_fractional_electron_density()
plot_optical_depth()
plot_visibility_function()

find_values()