# BRAMS - The New Brazilian developments on the Regional Atmospheric Modeling System

BRAMS is a numerical modeling system designed for regional scale atmospheric forecasting and research, with focus on atmospheric chemistry, air quality and biogeochemical cycles.

BRAMS is based on the Regional Atmospheric Modeling System (RAMS) originally developed at CSU/USA. BRAMS software is under a free license (CC-GPL). Currently, it is developed and maintained by CPTEC/INPE, USP and others institutions in Brazil and abroad. BRAMS/RAMS are multipurpose numerical weather prediction models that were designed to simulate atmospheric circulations spanning from hemispheric scales down to large eddy simulations (LES) of the planetary boundary layer. Several new functionalities and modifications were included to improve the numerical representation of key physical processes over tropical and subtropical regions (Freitas et al., 2005, 2009), a complete in-line module for atmospheric chemistry and aerosol processes (Longo et al., 2013), as well as a state-of-the-art model for simulation of atmosphere-surface exchange of water, energy, momentum, carbon and others biogeochemical tracers (Moreira et al., 2013).

The current version of BRAMS contains the following additional physical parameterizations:

    *Radiation:  CARMA (Toon et al., 1988) and RRTMG (Iacono et al., 2008) schemes for long- and short-wave, including aerosols effects and coupled with microphysics and convection schemes
    *Microphysics: a double moment from RAMS CSU version, Thompson double moment and aerosol aware (Thompson and Eidhammer, 2014).
    *Convection schemes:  Souza (1999) for shallow convection, Grell and Deveny (2002) for deep convection and Grell and Freitas (2014) scale and aerosol aware for deep and shallow convection including convective transport and wet removal of tracers.
    *Turbulence parameterizations: Nakanishi & Nino (2004) TKE based formulation, Taylor’s theory based formulation (Campos Velho, 1998).
    *Surface interaction and carbon cycle: TEB (Town Energy Budget) scheme to simulate urban areas (Freitas et al., 2007) and the Joint UK Land Environment Simulator (JULES) model (Moreira et al. 2013).
    *In-line emission, deposition, transport and chemistry model (CCATT) for atmospheric composition and air pollution studies (Longo et al., 2013).

Besides, BRAMS includes the following features:

    *Highly accurate and monotonic scheme (Walcek 2000, Freitas et al. 2012) for advection of scalars.
    *Complete, mass conservative formulation for the Exner function prognostic equation (Medvigy et al. 2005).
    *Computational Scalability: high efficiency up to ∼10000 cores with the MPI approach.
    *Digital filter for model initialization.
    *Model output at GrADS format during the model run time.
    *Coupling with the STILT Lagrangian Particle Dispersion Modelling.
    *Output to drive air parcels trajectory model calculation.

And ongoing work includes:

    *Fire spread model.
    *MATRIX (Bauer et al., 2008) aerosol model.
    *Data assimilation with GSI / 3D-VAR technique.
    *Time integration schemes (Runge-Kutta 3rd order, Adams-Bashforth-Moulton 3rd order, leapfrog with Robert-Asselin-Williams filter 2nd order).
    *Advection schemes (from 2nd to 6th order) with positivity and monotonicity constraints and WENO (Weighted Essentially Non-Oscillatory) formulation.
    *Code reformulation for fully compliant with the current FORTRAN standards.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

*BRAMS has been tested with the compilers: INTEL® compilers, PGI® compilers and GNU Fortran compiler (GPL). Follow the instructions of each site to install the compilers.
*BRAMS works only in parallel mode. One can run the model using a single processor/core, but must to do it using MPI with the MPIRUN command. We recommend download and install the last version of MPICH stable release. The last working version is 3.1.4. You may have problems with new versions, such as version 3.2. Also take care to choose the correct version to your OS.
After download the MPICH, please, proceed the installation:

```
Uncompress: ~\> tar -zxvf mpich-3.1.4.tar.gz
goto mpich directory: ~\> cd mpich-3.1.4
Configure mpich’s makefile: ~\> ./configure -disable-fast CFLAGS=-O2 FFLAGS=-O2 CXXFLAGS=-O2 FCFLAGS=-O2 -prefix=/opt/mpich3 CC=icc FC=ifort F77=ifort
make and install: ~\> make; ~\> sudo make install

```


### Making & Installing

After cloning, compile the model:

```
go to build directory: ~\>cd BRAMS/build/
configure the model’s makefile: ~\> ./configure -program-prefix=BRAMS -prefix=/home/xuser -enable-jules -with-chem=RELACS_TUV -with-aer=SIMPLE -with-fpcomp=/opt/mpich3/bin/mpif90 -with-cpcomp=/opt/mpich3/bin/mpicc -with-fcomp=ifort -with-ccomp=icc
make and make install: ~\>make; sudo make install
```

## Running the tests

Before run the model, You must get the input data. There are two packages for initial tests:

a) Small meteorological case for laptops and desktops (722MB)

* [](ftp://ftp.cptec.inpe.br/brams/BRAMS/data/meteo-only.tgz) - Meteorological small case

b) Small chemical case (using RELACS_TUV) for laptops and desktops (945 MB)

* [](ftp://ftp.cptec.inpe.br/brams/BRAMS/data/meteo-chem.tgz) - Meteorological/Chemical small case

Both cases are just for testing and learning processes. To get data for any different runs, visit the Input Data page.

```
    After downloading a test case, uncompress it using the command “tar”: ~> tar -zxvf meteo-only.tgz
    Goto the test directory: ~>cd meteo-only
    Make a tmp directory (necessary to Jules): ~>mkdir ./tmp
    Export the tmp directory (especially if you are running NOT locally !): ~>export TMPDIR=./tmp
    Take care about stack size of your machine. Make it at least 65536 or unlimited : ~>ulimit -s 65536
    Run the model using mpirun installed on step 3: ~>/opt/mpich3/bin/mpirun -np 4 ./brams-5.2.5

```
