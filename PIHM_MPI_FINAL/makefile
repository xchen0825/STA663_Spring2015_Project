# -----------------------------------------------------------------
# Version: 2.0
# Date: "Nov, 2007" 
# -----------------------------------------------------------------
# Programmer: Mukesh Kumar (muk139@psu.edu)@ PSU
# -----------------------------------------------------------------
# Makefile for PIHM 
#
# cvode/pihm/Makefile.  
# -----------------------------------------------------------------

SHELL = /bin/sh

srcdir       = .
builddir     = .
top_builddir = ../../
top_builddir = ../../
prefix       = /dscrhome/xc46/SUNDIALS/sundials 
exec_prefix  = ${prefix}
includedir   = /dscrhome/xc46/SUNDIALS/sundials/include
libdir       = /dscrhome/xc46/SUNDIALS/sundials/lib


CPP         = /usr/bin/cc -E
CPPFLAGS    =
CC          = /usr/bin/cc 
CFLAGS      = -g
LDFLAGS     =
LIBS        = -lm
MPICC       = /opt/apps/intelcsl//impi/5.0.2.044/intel64/bin/mpicc
MPI_INC_DIR = /opt/apps/intelcsl//impi/5.0.2.044/intel64/bin/../include
MPI_LIB_DIR = /opt/apps/intelcsl//impi/5.0.2.044/intel64/bin/../lib
MPI_LIBS    =
MPI_FLAGS   =

# -DSUB_SURF_RIV is used while solving all the components of the model (surface, subsurface and river states)
# -DDIFFUSION is used for diffusion wave approximation of St. Venant's equation. Otherwise kinematic wave assumption will be used
 MACRO    = -DSUB_SURF_RIV -DDIFFUSION -DLAYER3
# MACRO    = -DSUB_SURF_RIV -DDIFFUSION -DLAYER2
# MACRO    = -DSUB_SURF_RIV -DKINEMATIC -DLAYER2
# -DSURF_RIV is used to run the model without sub-surface component. Generally used for calibration purposes
# MACRO    = -DSURF_RIV -DDIFFUSION -DLAYER2
# -DNO_UNSAT is used to run the model without unsaturated zone component. Overlandflow is set to zero 0 initially. DY for overland flow, and unsat zone is set to 0. DY for groundwater = ElePrep. Our goal is to use -DSUB_SURF_RIV -DDIFFUSION -DLAYER2. Long simulations to steady state can be performed. Note: DY for overland should not be set to zero if we want to track where lakes and rivers will be created. Exfiltration will have to be incorporated then
# MACRO    = -DSUB_SURF_RIV -DDIFFUSION -DLAYER2 -DNO_UNSAT 
SRC    = pihm.c f.c read_alloc.c initialize.c is_sm_et.c update.c		 
#SRC    = pihm.c f.c read_alloc.c initialize.c is_sm_et.c update.c print.c  
 

COMPILER_PREFIX = 
LINKER_PREFIX   = 

SUNDIALS_INC_DIR = $(includedir)
SUNDIALS_LIB_DIR = $(libdir)
SUNDIALS_LIBS    = -lsundials_cvode -lsundials_shared -lsundials_nvecparallel

all:
	@(echo)
	@(echo '       make pihm     - make pihm        ')
	@(echo '       make clean    - remove all executable files')
	@(echo)

pihm:
	@echo '...Compiling PIHM (IntelMPI) ...'
	@$(MPICC) $(CFLAGS) -I$(SUNDIALS_INC_DIR) -I$(MPI_INC_DIR) -L$(SUNDIALS_LIB_DIR) -o $(builddir)/pihm $(SRC) $(SUNDIALS_LIBS) $(LIBS) $(MACRO) -L$(MPI_LIB_DIR) $(MPI_LIBS) $(LIBS)

clean:
	@rm -f *.o
	@rm -f pihm

