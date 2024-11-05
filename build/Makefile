#Makefile Make_model
############################## Change Log ##################################
# 5.2.5
#
############################################################################
# Define path and compilation


.SUFFIXES: .C .f90 .F90 .c .o

version        = 5.0.1
PREFIX         = /home/lufla/bin
HAVE_HDF5      = @HAVE_HDF5@
HDF5LIB        = 

CHEM	       = RELACS_TUV
AER 	       = SIMPLE
COLOR        = @COLOR@
MATRIXLEV	   = 1

CCOMP	       = gcc
CPCOMP	       = /opt/app/bin/mpicc
FCOMP	       = gfortran
FPCOMP	       = /opt/app/bin/mpif90

CMACH=       PC_LINUX1
PAR_DEFS=    -DRAMS_MPI

HAVE_ZLIB      = no
ZLIB        = 

WGRIBDIR = 
NETCDFFDIR = /opt/app/
NETCDFCDIR = /opt/app/
HDF5DIR = 
MPIDIR = @MPIDIR@

ENABLEJULES = true

MOD_EXT   = mod

include utils.mk

ifeq ($(ZLIB),)
ZLIB_LIB = 
else
ZLIB_LIB =  -L$(ZLIB) -lz
endif

ARCHIVE=ar rsu

BINDIR = $(PREFIX)/bin

PAR_INCS=./
PAR_LIBS=./

include paths.mk
include objects_model.mk

ifeq ($(WGRIBDIR),)
WGRIB_LIB  =
else
WGRIB_LIB  =-L$(WGRIBDIR) -lwgrib2_api -lwgrib2 
endif

ifeq ($(MPIDIR),)
MPI_LIB  =
else
MPI_LIB_LIB  =-L$(MPIDIR) -lmpi -lmpifort
endif


ifeq ($(HDF5DIR),)
HDF5_LIB = 
HDF5_INC =
else
HDF5_LIB = -L$(HDF5DIR)/lib -lhdf5 -lhdf5_hl -lhdf5hl_fortran -lhdf5_fortran
HDF5_INC = -I$(HDF5DIR)/include
endif

ifeq ($(NETCDFFDIR),)
NETCDFF_LIB = 
NETCDFC_LIB = 
else
NETCDFF_LIB = -L$(NETCDFFDIR)/lib -lnetcdff 
NETCDFC_LIB = -L$(NETCDFCDIR)/lib -lnetcdf
endif

ifeq ($(ENABLEJULES), true)
JULES_DEFS0 =JULES SCMA BL_DIAG_HACK MPI_DUMMY
else
JULES_DEFS0 =
JULESID    =
JULES_DEFS0 = SCMA BL_DIAG_HACK GNU_FORTRAN MPI_DUMMY
LIBJULES   =
JULES_LIB  =
JULES_INC  =
endif

#including ifs cumulus scheme
ifeq ($(IFS), ifson)
LIBIFS=ifs.a
else
LIBIFS=
endif

ifndef FPCOMP
$(error You need to set the variable FPCOMP either in configure command line)
endif

ifeq ($(notdir $(FPCOMP)),mpif90)
FC=$(FPCOMP) #-f90=$(FCOMP)
FC_LOADER=$(FC)
    # For intel compiler debugging
	#FC_FLAGS=-O0 -mcmodel=large -check bounds -traceback -fpe0 -check uninit -init=snan,arrays
	#CCFLAGS=-O0
endif

ifeq ($(notdir $(FPCOMP)),ftn)
JULES_DEFS = $(JULES_DEFS0) CRAY_FORTRAN
FC=$(FPCOMP)
FC_LOADER=$(FC)
endif

ifeq ($(notdir $(FPCOMP)),h5pfc)
JULES_DEFS = $(JULES_DEFS0)
FC=$(FPCOMP)
FC_LOADER=$(FC)
endif


ifeq ($(notdir $(FCOMP)),pgf90)
JULES_DEFS = $(JULES_DEFS0)
FC_FLAGS=-fast -O3 -fastsse -Munroll=n:4 #-Mipa=fast,inline,safe  -Msmart
CCFLAGS=-O2
FPP=-Mpreprocess
endif

ifeq ($(notdir $(FCOMP)),ifort)
JULES_DEFS = $(JULES_DEFS0) INTEL_FORTRAN
FC_FLAGS=-O3
CCFLAGS=-O3
    # For intel compiler debugging
	#FC_FLAGS=-O0 -mcmodel=large -check bounds -traceback -fpe0 -check uninit -init=snan,arrays
	#CCFLAGS=-O0
FPP=-fpp
FPP_VARSI = INTEL
endif

ifeq ($(notdir $(FPCOMP)),mpifort) 
JULES_DEFS = $(JULES_DEFS0) INTEL_FORTRAN
FC=$(FPCOMP) #-traceback
CC=$(CPCOMP)

FC_LOADER=$(FC) #-traceback
CC_LOADER=$(CC)
endif

ifeq ($(notdir $(FPCOMP)),mpiifort) 
JULES_DEFS = $(JULES_DEFS0) INTEL_FORTRAN
FC=$(FPCOMP) #-traceback
CC=$(CPCOMP)
FC_LOADER=$(FC) #-traceback
CC_LOADER=$(CC)
endif

ifeq ($(notdir $(FCOMP)),gfortran)
JULES_DEFS = $(JULES_DEFS0) GNU_FORTRAN
#FC_FLAGS=-O3 -ffree-form -ffree-line-length-none -fbacktrace -g -fcheck=bounds -fcheck=all -finit-real=nan # -Wall -ffpe-trap=invalid,zero
FC_FLAGS=-ffree-form -ffree-line-length-none -fbacktrace -g -fno-range-check -fcheck=all -fcheck=bounds

ifeq ($(GCC_VER_GTE10),)
     FC_FLAGS +=-fallow-argument-mismatch
endif

OT_FLAGS=-O3
#CCFLAGS=-O2
CCFLAGS=-O0 -fpermissive -std=gnu89
FPP=-cpp
FC=$(FPCOMP) $(OT_FLAGS)
endif

#FCFLAGS=$(FC_FLAGS) -Wl,--no-as-needed -ldl
FCFLAGS=$(FC_FLAGS) -ldl


LOADERFLAGS=$(FCFLAGS)

#For jules
FPP_INC=-I
FPP_FDEF=-D
INCLUDE_DIR=$(JULES_DIR)/preprocess/include/

SFIRE_INC = -I$(SFIRE)

#########################################################
## Build pre-processor variable definitions for jules  ##
#########################################################
ifeq ($(NETCDFFDIR),)
FPP_VARS0=$(FPP_VARSI) $(AER) $(JULES_DEFS)	
else
FPP_VARS0=$(FPP_VARSI) $(AER) $(JULES_DEFS) cdf
endif

ifeq ($(WGRIBDIR),)
FPP_VARS1=$(FPP_VARS0)
else
FPP_VARS1=$(FPP_VARS0) GRIB2
endif

FPP_VARS=$(FPP_VARS1)
ifeq ($(HDF5DIR),)
FPP_VARS=$(FPP_VARS1)
else
FPP_VARS=$(FPP_VARS1) HDF5
endif


FPP_DEFS=$(foreach FPP_VAR,$(FPP_VARS),$(FPP_FDEF)$(FPP_VAR))
FPP_INC_PATH=$(FPP_INC)$(INCLUDE_DIR)

ifeq ($(NETCDFFDIR),)
ifeq ($(WGRIBDIR),)
INCLUDES  = -I$(UTILS_INCS) $(JULES_INC) $(HDF5_INC) $(FPP_INC_PATH) $(SFIRE_INC)
else
INCLUDES  = -I$(UTILS_INCS) $(JULES_INC) $(HDF5_INC) $(FPP_INC_PATH) -I$(WGRIBDIR)/include $(SFIRE_INC)
endif
else
ifeq ($(WGRIBDIR),)
INCLUDES  = -I$(UTILS_INCS) $(JULES_INC) $(HDF5_INC) $(FPP_INC_PATH) -I$(NETCDFFDIR)/include $(SFIRE_INC)
else
INCLUDES  = -I$(UTILS_INCS) $(JULES_INC) $(HDF5_INC) $(FPP_INC_PATH) -I$(WGRIBDIR)/include -I$(NETCDFFDIR)/include $(SFIRE_INC)
endif
endif


## Define fancy Make  ###
ifndef V
	VERBOSE_FC = @printf '%b ' $(BOLD_BLUE)$(FC)$(NO_COLOUR) \
	                           $(BOLD)$(FCFLAGS)$(NO_COLOUR) \
														 $(BOLD_GREEN)$(INCLUDES)$(NO_COLOUR) \
														 $(BOLD_YELLOW)$(FPP_DEFS)$(NO_COLOUR) \
														 $(NO_COLOUR)$?$(NO_COLOUR)    \
	                           $(CYAN)"=="$(NO_COLOUR)       \
	                           $(BOLD_RED) $@$ $(NO_COLOUR);   \
	                           printf '\n';

endif

# Compiler commands.

EXTRAFLAGSF = $(INCLUDES) $(FCFLAGS) $(FPP_DEFS)
#EXTRAFLAGSF = $(JULES_DEFS) $(INCLUDES) $(FCFLAGS) $(FPP)
F_COMMAND = $(VERBOSE_FC) $(FC) -c
F_COMMAND_MACRO = $(F_COMMAND) -D$(AER) -D$(MATRIXLEV) -D$(COLOR)
F_COMMAND_LIGHT = $(FPCOMP) -O0 -c

# Define objects.
include objects_model.mk

# Define targets.

all: $(EXE)

$(EXE): libjules $(LIBUTILS) $(LIBIFS) $(LIBMODEL) $(MAINOBJ) $(WGRIBDIR) $(NETCDFFDIR)
	@echo ""
	#$(FC_LOADER) -o $(EXE) rammain.o $(LOADERFLAGS) $(LIBMODEL) $(LIBUTILS) $(LIBIFS) $(JULES_LIB) -L$(PAR_LIBS) $(HDF5_LIB)
	$(FC_LOADER) -o $(EXE) rammain.o $(LOADERFLAGS) $(LIBMODEL) $(LIBUTILS) $(LIBIFS) -L$(PAR_LIBS) $(WGRIB_LIB) $(NETCDFF_LIB) $(NETCDFC_LIB) $(HDF5_LIB) $(MPI_LIB)
	@echo ""
	@echo Finished building === $(EXE)
	@echo ""

$(MAINOBJ):  $(MAIN)
	./getGitInfo.sh
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@rm -f $(<F:.f90=.f90)

$(LIBMODEL): $(OBJ_MODEL)
	$(ARCHIVE) $(LIBMODEL) $(OBJ_MODEL) $(JULES_OBJ_SFCLYR) $(JULES_OBJ_MEM)

$(LIBUTILS):
	$(MAKE) -f Make_utils

$(LIBIFS):
	$(MAKE) -f Make_ifs

libjules:
ifeq ($(ENABLEJULES), true)
	#cd ../src/jules/LIB && $(MAKE)
endif
ifeq ($(notdir $(FCOMP)),gfortran)

endif

filter:
	cd ./ && $(MAKE) -f Make_filter

install-filter:
	\cp -f -R $(RAMS_ROOT)/build/filter $(BINDIR)/filter

pre-brams:
	cd ./ && $(MAKE) -f Make_pre

install-pre-brams:
	\cp -f -R $(RAMS_ROOT)/build/pre $(BINDIR)/pre

clean:
	$(MAKE) -f Make_utils  clean
	$(MAKE) -f Make_ifs clean
	rm -f $(OBJ_MODEL) $(MOD_MODEL) $(MAINOBJ) $(JULES_OBJ_SFCLYR) $(JULES_OBJ_MEM) $(EXE) $(LIBMODEL) *.oo

install:
	install -D $(EXE) $(BINDIR)/$(notdir $(EXE))
	\cp -f -R $(RAMS_ROOT)/jules.in $(BINDIR)/jules.in
	\cp -f $(RAMS_ROOT)/variables.csv $(BINDIR)/variables.csv
	\cp -f -R $(RAMS_ROOT)/namelists $(BINDIR)/
	\cp -f -R $(RAMS_ROOT)/namelists/tables/* $(BINDIR)/tables/
	\cp -f $(RAMS_ROOT)/scripts/* $(BINDIR)/

install-strip:
	install -D -s $(EXE) $(BINDIR)/$(notdir $(EXE))
	\cp -f $(RAMS_ROOT)/jules.in $(BINDIR)/jules.in
#	\cp -rf $(RAMS_ROOT)/tables $(BINDIR)/tables
	\cp -f $(RAMS_ROOT)/RAMSIN* $(BINDIR)
	\cp -f $(RAMS_ROOT)/variables.csv $(BINDIR)/variables.csv

uninstall:
	rm -f $(BINDIR)/$(notdir $(EXE))
	rm -f $(BINDIR)/jules.in
#	rm -rf $(BINDIR)/tables
	rm -f $(BINDIR)/RAMSIN*
	rm -f $(BINDIR)/variables.csv

# Define dependencies
ifeq ($(NETCDFFDIR),)
jules_netcdf_dummy.o: $(JULES_UNCDF)/jules_netcdf_dummy.f90
	@cp -f  $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src	
endif
include depend_model.mk

.PHONY: all libjules install install-strip uninstall clean filter
