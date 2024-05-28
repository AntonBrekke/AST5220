# Cosmology II, AST5220
## This repo contains the project for the AST5220 course. 

Compile files by running "make cmb" in the same folder as the Makefile. 
Run file in main folder by "./cmb". To run the project, you would have to install Gnu Scientific Library (GSL). 

### Milestone 1 
Calcualte the Background Cosmology.
To run Milestone1.py, make sure to change/comment out the path of "path" and "savefig_path".

### Milestone 2
Calculate the Recombination History for our universe. 
To run Milestone2.py, make sure to change/comment out the path of "path" and "savefig_path".

### Milestone 3
Calculate perturbations in the metric of the Universe.
To run Milestone3.py, make sure to change/comment out the path of "path" and "savefig_path".

### Milestone 4
Calculate Powerspectrum of the CMB.
To run Milestone4.py, make sure to change/comment out the path of "path" and "savefig_path".


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
