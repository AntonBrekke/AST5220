# Hans A. Winther (2020) (hans.a.winther@gmail.com)

SHELL := /bin/bash

# Set compiler (use =c++17 if you have this availiable)
CC = g++ -std=c++11 

# Paths to GSL library
INC  = -I ./include -I$(HOME)/local/include
LIBS = -L$(HOME)/local/lib -lgsl -lgslcblas


#=======================================================
# Options
#=======================================================
OPTIONS = 

# Add bounds checking
OPTIONS += -D_GLIBCXX_DEBUG

# Show warnings if atempting to evaluate a spline out of bounds
OPTIONS += -D_SPLINE_WARNINGS_ON

# Show info about the solution as we integrate
# OPTIONS = -D_FIDUCIAL_VERBOSE_ODE_SOLVER_TRUE

# Add OpenMP parallelization
# OPTIONS += -D_USEOPEMP
# CC += -fopenmp

#=======================================================

C = -O3 -g $(OPTIONS)

#=======================================================

VPATH=src/
TARGETS := cmb
all: $(TARGETS)

# OBJECT FILES
OBJS = Main.o Utils.o BackgroundCosmology.o RecombinationHistory.o Perturbations.o PowerSpectrum.o Spline.o ODESolver.o

# DEPENDENCIES
Main.o                  : include/BackgroundCosmology.h include/RecombinationHistory.h include/Perturbations.h include/PowerSpectrum.h
Spline.o                : include/Spline.h
ODESolver.o             : include/ODESolver.h
Utils.o                 : include/Utils.h include/Spline.h include/ODESolver.h
BackgroundCosmology.o   : include/BackgroundCosmology.h include/Utils.h include/Spline.h include/ODESolver.h
RecombinationHistory.o  : include/RecombinationHistory.h include/BackgroundCosmology.h
Perturbations.o         : include/Perturbations.h include/BackgroundCosmology.h include/RecombinationHistory.h
PowerSpectrum.o         : include/PowerSpectrum.h include/BackgroundCosmology.h include/RecombinationHistory.h include/Perturbations.h
Examples.o              : include/Utils.h include/Spline.h include/ODESolver.h

examples: Examples.o Utils.o Spline.o ODESolver.o
	${CC} -o $@ $^ $C $(INC) $(LIBS)

cmb: $(OBJS)
	${CC} -o $@ $^ $C $(INC) $(LIBS)

%.o: %.cpp
	${CC}  -c -o $@ $< $C $(INC) 

clean:
	rm -rf $(TARGETS) *.o 
