import numpy as np
import matplotlib.pyplot as plt
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

def plot_data(datatype='Theta0', factor=1, title='', ylabel='', xscale='linear', yscale='linear', show=True):
    k_list = [0.1, 0.01, 0.001]

    fig = plt.figure()
    ax = fig.add_subplot()

    ax.set_title(title, fontsize=16)

    for ki in k_list:
        # Get data
        perturbations_data = np.loadtxt(fr'{path}' + f'/perturbations_k{ki}.txt')
        x, Theta0, Theta1, Theta2, Phi, Psi, Source_T, delta_cdm, delta_b, v_cdm, v_b, T5, T50, T500 = perturbations_data.T
        data_dict = {'x': x, 'Theta0': Theta0, 'Theta1': Theta1, 'Theta2': Theta2,
                    'Phi': Phi, 'Psi': Psi, 'Source_T': Source_T, 'delta_cdm': delta_cdm,
                    'delta_b': delta_b, 'v_cdm': v_cdm, 'v_b': v_b, 'T5': T5, 
                    'T50': T50, 'T500': T500, 'Phi + Psi': Phi + Psi}

        ax.plot(x, factor*data_dict[datatype], label=f'k={ki}', lw=2)

    ax.set_xlim(-15, 0)
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.set_xlabel('x=ln(a)', fontsize=16)
    ax.set_ylabel(ylabel, fontsize=16)
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    ax.grid(True)

    ax.legend(prop={'size': 14}, loc='best', frameon=False)
    fig.tight_layout()
    plt.savefig(savefig_path + rf'/{datatype}.pdf')
    if show is True: plt.show()

def plot_density_perturbations(show=True):
    # Plot density perturbations
    plot_data('Theta0', factor=4, title=r'Density perturbation $\delta_{\gamma}=4\Theta_0$', show=show)
    plot_data('delta_cdm', title=r'Density perturbation $\delta_{\rm CDM}$', yscale='log', show=show)
    plot_data('delta_b', title=r'Density perturbation $\delta_{\rm b}$', yscale='log', show=show)

def plot_velocity_perturbations(show=True):
    # Plot velocity perturbations
    plot_data('Theta1', factor=-3, title=r'Velocity perturbation $v_{\gamma}=-3\Theta_1$', show=show)
    plot_data('v_cdm', title=r'Velocity perturbation $v_{\rm CDM}$', yscale='log', show=show)
    plot_data('v_b', title=r'Velocity perturbation $v_{\rm b}$', yscale='log', show=show)

def plot_photon_quadropole(show=True):
    # Plot photon quadropole
    plot_data('Theta2', title=r'Photon quadropole $\Theta_2$', show=show)

def plot_potentials(show=True):
    # Plot potentials 
    plot_data('Phi', title=r'Potential $\Phi$', show=show)
    plot_data('Phi + Psi', title=r'Sum: $\Phi + \Psi$', show=show)

plot_density_perturbations(show=False)
plot_velocity_perturbations(show=False)
plot_photon_quadropole(show=False)
plot_potentials(show=False)