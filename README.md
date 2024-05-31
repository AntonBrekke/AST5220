# Cosmology II, AST5220
## Abstract 
> The Cosmic Microwave Background (CMB) radiation is a relic from the early universe, providing a snapshot of the cosmos approximately 380,000 years after the Big Bang. This faint glow, permeating the entire sky, holds valuable information about the composition and evolution of the Universe. \\
Through the $\Lambda\rm CDM$-model with initial values from \cite{Planck_data}, we develop an Einstein-Boltzmann solver to make a prediction of the CMB and matter power spectrum by considering linear perturbation theory on a flat FLRW universe as our background cosmology. Comparing our model to data from the supernova \cite{Supernovadata} we find that demanding a flat Universe forces us to include dark energy $\Omega_{\rm DE}\sim0.66-0.76$. We also predict the acceleration of the Universe to happen when $\Omega_{\Lambda}\approx \frac{1}{3}$. We also found that during recombination, the event of last scattering happened at time $t=0.37801\,\rm Myr$, as well as a freeze-out abundance $X_e(x=0)=2.02\cdot10^{-4}$ today. Finally, we obtain results for the power spectrums, giving us an esitimated theoretical upper bound $\lambda\approx 486\,\rm Mpc=1.62\,\rm Gly$ for the size of cosmological structures. Comparing the power spectrums to data from \cite{Planck_data}, \cite{Galaxy-survey} and \cite{WMAP} we also find inaccuracies in our predicted model due to ignoring neutrinos, polarization and matter heavier than hydrogen. 

The report can be found [here](https://github.com/AntonBrekke/AST5220/tree/main/report).

This project was done during the course AST5220 at ITA University of Oslo, by the help of templates found at [this GitHub](https://github.com/HAWinther/AST5220-Cosmology/tree/master). The project is mainly done in C++, with expceptions to plotting/analyzing numerical results. During the project, we develop an Einstein-Boltzmann solver, and essentially solve for the CMB power spectrum. 

### Website
All relevant information about the project and the different milestones can be found on this [website](https://cmb.wintherscoming.no/).

### Running the codes
Compile files by typing "make cmb" in the terminal while being in the same folder as the Makefile. 
Run file in main folder by typing "./cmb" in terminal. To run the project, you would have to install Gnu Scientific Library (GSL) (see how below). 

#### Milestone 1 
Calcualte the Background Cosmology.
To run Milestone1.py, make sure to change/comment out the path of "path" and "savefig_path" in the code.

#### Milestone 2
Calculate the Recombination History for our universe. 
To run Milestone2.py, make sure to change/comment out the path of "path" and "savefig_path" in the code.

#### Milestone 3
Calculate Perturbations in the metric of the Universe.
To run Milestone3.py, make sure to change/comment out the path of "path" and "savefig_path" in the code.

#### Milestone 4
Calculate Powerspectrum of the CMB.
To run Milestone4.py, make sure to change/comment out the path of "path" and "savefig_path" in the code.

## How to install GSL

See [this](https://solarianprogrammer.com/) for how to install it on a Windows machine. On Linux or a Mac you can either use a package manager or install it directly as follows:

- Go the the home directory:

cd $HOME

- Make a local folder to hold libraries:

mkdir local

- Enter this directory:

cd local

- Download the code (if you don't have wget you need to get the file to this dir by other means):

wget ftp://ftp.gnu.org/gnu/gsl/gsl-2.6.tar.gz

- Untar the code:

tar -xvf gsl-2.6.tar.gz

- You should now have the gsl-2.6 folder. Enter it:

cd gsl-2.6

- Run the configure script:

./configure --prefix=$HOME/local

- Compile and install it:

make ; make install

- In the CMB code Makefile change the include and lib paths to point to the library:

INC  = -I$(HOME)/local/include
LIBS = -L$(HOME)/local/lib -lgsl -lgslcblas

- If this fails with "libgsl.so not found" then run the command:

export LD\_LIBRARY\_PATH="$LD\_LIBRARY\_PATH:$HOME/local/lib"

and try to run ./cmb again and it should work. To avoid having
to run this command every time you open a new terminal open
the $HOME/.bashrc file and add this line to the end of the file
and it will load everytime you open a new window.
