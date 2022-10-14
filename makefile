EXEC = newkry
HOST = gcc/11-openmp


# For linux systems, it is assumed that the environment
# variable LD_LIBRARY_PATH contains the locations to libfmm3d.so
# and libsolvers3d.so, for Macosx, these .so files also need to be
# copied over /usr/local/lib


FMMBIE_INSTALL_DIR= $(PREFIX)
FMM_INSTALL_DIR   = $(PREFIX_FMM)
LFMMLINKLIB       = -lfmm3d
LLINKLIB          = -lfmm3dbie
LBLAS             =  flexiblas
FMM_INSTALL_DIR   =${HOME}/lib
FMMBIE_INSTALL_DIR=${HOME}/lib

FC    =  f2py3 
FFLAGS= -finit-real=zero -fPIC -O3 -funroll-loops -march=broadwell -fopenmp -std=legacy 
FEND  = -l$(LBLAS) -L${FMMBIE_INSTALL_DIR} $(LLINKLIB)

.PHONY: all clean 


OBJECTS =  constants.o
SRC     =  constants.f90  oocyterout.f90 axisymm.f90 normalcal.f90  matroutine.f90 chebcoeff.f90 bingham.f90 stokgmres.f gradgmres.f90 activenewton.f90
		

# use only the file part of the filename, then manually specify
# the build location


%.o : %.f90  
	$(FC) -c $(FFLAGS) $< -o $@ $(FEND)

all: $(SRC)
	$(FC) $(FEND) -c -m $(EXEC) $(SRC) 

clean:
	rm -f $(OBJECTS)

list: $(SOURCES)
	$(warning Requires:  $^)


