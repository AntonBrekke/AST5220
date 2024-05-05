import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mpcl
import scipy.constants as scc
import tabulate as tab
import latextable as lt
import texttable as txttab

# Relevant paths to get data
path = r'/home/antonabr/AST5220/data'
savefig_path = r'/home/antonabr/AST5220/figures/Milestone3'

# Changing standard color cycle to preferred cycle  
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
colors.insert(6, colors[1]) ; colors[1:5] = colors[2:5]
plt.rcParams['axes.prop_cycle'] = plt.cycler(color=colors) 

Mpc = 3.08567758e22

def find_horizon_crossing(eta, k):
    k /= Mpc
    return np.where(np.min(abs(eta*k - 1)) == abs(eta*k - 1))[0][0]

# Function that plots one of the given datas (way overcomplicated)
def plot_data(data_list=['Theta_0'], factor=[1], kvals=[0.1, 0.01, 0.001], title='', ylabel='', ylim=None, xscale='linear', yscale='linear', loc='best', show=True, print_hc=False):
    if type(data_list) != list: data_list = [data_list] 
    if type(factor) != list: factor = [factor] 
    if len(factor) < len(data_list): factor = len(data_list)*factor
    k_list = kvals
    styles = ['-', '--', ':', '-.']
    N_k = len(k_list)
    N_styles = len(styles)

    fig = plt.figure()
    ax = fig.add_subplot()

    ax.set_title(title, fontsize=16)

    # For each data to be plotted in same plot, plot each k
    for data_index, data in enumerate(data_list):
        linestyle = styles[data_index % N_styles]
        for k_index, ki in enumerate(k_list):
            # Get data
            perturbations_data = np.loadtxt(fr'{path}' + f'/perturbations_k{ki}.txt')
            x, Theta0, Theta1, Theta2, Phi, Psi, Source_T, delta_cdm, delta_b, v_cdm, v_b, eta_of_x, T5, T50, T500 = perturbations_data.T

            data_dict = {'x': x, 'Theta_0': Theta0, 'Theta_1': Theta1, 'Theta_2': Theta2,
                        'Phi': Phi, 'Psi': Psi, 'Source_T': Source_T, 'delta_CDM': abs(delta_cdm),
                        'delta_b': abs(delta_b), 'v_CDM': abs(v_cdm), 'v_b': abs(v_b), 'eta': eta_of_x, 
                        'T5': T5, 'T50': T50, 'T500': T500, 'Phi + Psi': Phi + Psi}

            horizon_crossing_index = find_horizon_crossing(eta_of_x, ki)
            color = colors[k_index % N_k]
            if data_index == 0: 
                label = f'k={ki}/Mpc'
                ax.axvline(x[horizon_crossing_index], color=0.75*np.array(mpcl.to_rgb(color)), lw=2, linestyle='-.')
            else: 
                label = None 
            ax.plot(x, factor[data_index]*data_dict[data], color=color, label=label, lw=2, linestyle=linestyle)
            if print_hc:
                print(f'k={ki}: x_hc={x[horizon_crossing_index]}')
        # Checklist for which types should not have '\' in $$.
        checklist = ['v_CDM', 'v_b', 'T5', 'T50', 'T500']
        # Tell which data is which linestyle 
        if len(data_list) > 1: 
            string =  r"_{\rm ".join(data.split("_")) + "}"*('_' in data)
            if data in checklist: label =  f'{factor[data_index]}'*(factor[data_index] != 1) + rf'${string}$'
            else: label = f'{factor[data_index]}'*(factor[data_index] != 1) + rf'$\{string}$'
            ax.plot([], [], color='dimgray', label=label, linestyle=linestyle)

    ax.axvspan(-7.25, -6.5, color='grey', label='Recombination', alpha=0.3)
    ax.set_xlim(-14, 0)
    if ylim is not None: ax.set_ylim(ylim)
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.set_xlabel('x=ln(a)', fontsize=16)
    ax.set_ylabel(ylabel, fontsize=16)
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    ax.grid(True)

    ax.legend(prop={'size': 14}, loc=loc, frameon=False)
    fig.tight_layout()
    plt.savefig(savefig_path + rf'/{"_".join(data_list)}.pdf')
    if show is True: plt.show()

def plot_density_perturbations(show=True):
    # Plot density perturbations
    plot_data('Theta_0', factor=4, title=r'Density perturbation $\delta_{\gamma}=4\Theta_0$', loc='lower right', show=show, print_hc=True)
    plot_data(['delta_CDM', 'delta_b'], title=r'Density perturbations $\delta_{\rm CDM}$, $\delta_{\rm b}$', yscale='log', loc=(0.6, 0), show=show)
    # plot_data('delta_CDM', title=r'Density perturbation $\delta_{\rm CDM}$', yscale='log', show=show)
    # plot_data('delta_b', title=r'Density perturbation $\delta_{\rm b}$', yscale='log', show=show)

def plot_velocity_perturbations(show=True):
    # Plot velocity perturbations
    plot_data('Theta_1',  factor=-3, title=r'Velocity perturbation $v_{\gamma}=-3\Theta_1$', show=show)
    plot_data(['v_CDM', 'v_b'], title=r'Velocity perturbations $v_{\rm CDM}$, $v_{\rm b}$', yscale='log', show=show)
    # plot_data(['Theta_0', 'Theta_1'], factor=[4, -3], show=show)
    # plot_data('v_CDM', title=r'Velocity perturbation $v_{\rm CDM}$', yscale='log', show=show)
    # plot_data('v_b', title=r'Velocity perturbation $v_{\rm b}$', yscale='log', show=show)

def plot_photon_quadropole(show=True):
    # Plot photon quadropole
    plot_data('Theta_2', title=r'Photon quadropole $\Theta_2$', loc='upper left', show=show)

def plot_potentials(show=True):
    # Plot potentials 
    plot_data('Phi', title=r'Gravitational potential $\Phi$', loc=(0.6,0.2), show=show)
    plot_data('Phi + Psi', title=r'Sum of potentials: $\Phi + \Psi$', loc='upper left', show=show)

def compare_baryon_photon(show=True):
    k_list = [0.1, 0.01, 0.001]
    ylims_1 = [[-3, 4], [0.5, 4], [0.75, 3]]
    ylims_2 = [[-1.5, 1.5], [-0.75, 1], [-0.2, 0.4]]
    loc1 = ['upper right', 'upper left', 'upper left']
    loc2 = ['lower right', 'upper left', 'upper left']
    for i, ki in enumerate(k_list):
        perturbations_data = np.loadtxt(fr'{path}' + f'/perturbations_k{ki}.txt')
        x, Theta0, Theta1, Theta2, Phi, Psi, Source_T, delta_cdm, delta_b, v_cdm, v_b, eta_of_x, T5, T50, T500 = perturbations_data.T

        horizon_crossing_index = find_horizon_crossing(eta_of_x, ki)

        fig = plt.figure(figsize=(6.4, 7))
        ax1 = fig.add_subplot(211)
        ax2 = fig.add_subplot(212)
        ax1.set_title(r'Density perturbations $\delta_\gamma,\delta_{\rm b}$', fontsize=16)
        ax2.set_title(r'Velocity perturbations $v_\gamma,v_{\rm b}$', fontsize=16)

        label = f'k={ki}/Mpc'

        # Densities
        ax1.plot(x, 4*Theta0, label=r'$\delta_\gamma$(' + f'k={ki}/Mpc)', lw=2, linestyle='-')
        ax1.plot(x, delta_b, color='tab:red', label=r'$\delta_{\rm b}$(' + f'k={ki}/Mpc)', lw=2, linestyle='--')
        # Velocities
        ax2.plot(x, -3*Theta1, label=r'$v_\gamma$(' + f'k={ki}/Mpc)', lw=2, linestyle='-')
        ax2.plot(x, v_b, color='tab:red', label=r'$v_{\rm b}$(' + f'k={ki}/Mpc)', lw=2, linestyle='--')

        ax1.axvline(x[horizon_crossing_index], color='k', lw=2, linestyle='-.')
        ax2.axvline(x[horizon_crossing_index], color='k', lw=2, linestyle='-.')

        ax1.axvspan(-7.25, -6.5, color='grey', label='Recombination', alpha=0.3)
        ax2.axvspan(-7.25, -6.5, color='grey', label='Recombination', alpha=0.3)
        ax1.set_xlim(-14, 0)
        ax1.set_ylim(ylims_1[i])
        ax2.set_xlim(-14, 0)
        ax2.set_ylim(ylims_2[i])

        ax1.tick_params(axis='both', which='major', labelsize=14)
        ax2.tick_params(axis='both', which='major', labelsize=14)
        ax2.set_xlabel('x=ln(a)', fontsize=16)

        ax1.legend(prop={'size': 14}, frameon=False, loc=loc1[i])
        ax2.legend(prop={'size': 14}, frameon=False, loc=loc2[i])
        ax1.grid(True)
        ax2.grid(True)
        fig.tight_layout()
        plt.savefig(savefig_path + f'/compare_photon_baryon_k={ki}.pdf')
        if show is True: plt.show()

# Control center for running plots
plot_density_perturbations(show=True)
plot_velocity_perturbations(show=True)
plot_photon_quadropole(show=True)
plot_potentials(show=True)
compare_baryon_photon(show=True)