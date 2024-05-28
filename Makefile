# Makefile for project in AST5220

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
OPTIONS += -D_USEOPEMP
CC += -fopenmp

#=======================================================

C = -O3 -g $(OPTIONS)

#=======================================================

# Check both src and include for dependencies
VPATH=src/ : include/  
TARGETS := cmb
all: $(TARGETS)

# Create obj directory if it doesn't exist
$(shell mkdir -p obj)

# OBJECT FILES
OBJS = obj/Main.o obj/Utils.o obj/BackgroundCosmology.o obj/RecombinationHistory.o obj/Perturbations.o obj/PowerSpectrum.o obj/Spline.o obj/ODESolver.o

# DEPENDENCIES
obj/Main.o                  : BackgroundCosmology.h RecombinationHistory.h Perturbations.h PowerSpectrum.h
obj/Spline.o                : Spline.h
obj/ODESolver.o             : ODESolver.h
obj/Utils.o                 : Utils.h Spline.h ODESolver.h
obj/BackgroundCosmology.o   : BackgroundCosmology.h Utils.h Spline.h ODESolver.h
obj/RecombinationHistory.o  : RecombinationHistory.h BackgroundCosmology.h
obj/Perturbations.o         : Perturbations.h BackgroundCosmology.h RecombinationHistory.h
obj/PowerSpectrum.o         : PowerSpectrum.h BackgroundCosmology.h RecombinationHistory.h Perturbations.h
obj/Examples.o              : Utils.h Spline.h ODESolver.h

examples: obj/Examples.o obj/Utils.o obj/Spline.o obj/ODESolver.o
	${CC} -o $@ $^ $C $(INC) $(LIBS)

cmb: $(OBJS)
	${CC} -o $@ $^ $C $(INC) $(LIBS)

obj/%.o: %.cpp
	${CC}  -c -o $@ $< $C $(INC) 

clean:
	rm -rf $(TARGETS) obj/*.o 
