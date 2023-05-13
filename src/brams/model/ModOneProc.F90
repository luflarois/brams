module ModOneProc
	!# Initialization whenever all processes compute
	!#
	!# @note
	!# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
	!#
	!# **Brief**: Initialization whenever all processes compute; master also does io_ 
	!#
	!# **Documentation**: <http://twixar.me/kW2T>
	!#
	!# **Author**: Jairo Panetta **&#9993;**<mailto:jairo.panetta@gmail.com>
	!#
	!# **Date**: 2020-04-16
	!# @endnote
	!#
	!# @changes
	!#**Changelogs:**
	!# 
	!# @endchanges
	!# @bug
	!# **Open Bugs:**
	!#
	!# @endbug
	!#
	!# @todo
	!# **Todo list:**
	!#
	!# @endtodo
	!#
	!# @warning
	!# Now is under CC-GPL License, please see
	!# &copy; <https://creativecommons.org/licenses/GPL/2.0/legalcode.pt>
	!# @endwarning
	!#
	!#--- ----------------------------------------------------------------------------------------
	
  use ModVarfFile, only: &
       VarfReadStoreOwnChunk

  use module_dry_dep, only: &
       dep_init        ! Subroutine

  use memSoilMoisture, only : &
       SOIL_MOIST ! INTENT(IN)

  use mem_gaspart, only   : &
       gaspart_g ! intent(inout)

  use mem_teb, only       : &
       teb_g     ! intent(inout)

  use mem_teb_common, only: &
       tebc_g    ! intent(inout)

  use teb_vars_const, only: &
       iteb, &      ! intent(in)
       StoreNamelistFileAtTeb_vars_const

  use var_tables, only: &
       num_var,         &
       vtab_r

  use mem_micro, only: &
       micro_g

  use ModTimestep    , only: timestep
  use ModTimestep_RK , only: timestep_rk
  use ModTimestep_ABM, only: timestep_abm
  use ISO_FORTRAN_ENV, only: OUTPUT_UNIT
  use io_params, only: &
       ndvitime2, &
       ndviflg, &
       iupdndvi, &
       ssttime2, &
       isstflg, &
       iupdsst, &
       frqprt, &
       timstr, &
       iopunt,   & ! intent(out)
       srctime2, &
       ipastin,  &
       iinput,  &
       ioutput, &
       frqlite, &
       avgtim,  &
       initfld, &
       StoreNamelistFileAtIo_params

  use Isan_coms, only: &
      StoreNamelistFileAtIsan_coms, &
      ICFILETYPE

  use mem_cuparm, only: &
       ncufl, &
       cu_times, &
       NNQPARM, &
       if_cuinv, &
       cuparm_g, &
       StoreNamelistFileAtMem_cuparm

  use Mem_globrad, only: &
       master_read_carma_data, &
       StoreNamelistFileAtMem_globrad

  use Mem_grell_param, only: StoreNamelistFileAtMem_grell_param

  use Mem_leaf, only: &
       StoreNamelistFileAtMem_leaf, &
       isfcl, &
       leaf_g

  use Mem_oda, only: &
       if_oda, &
       StoreNamelistFileAtMem_oda

  use mem_radiate, only: &
       iswrtyp,          &
       ilwrtyp, &
       StoreNamelistFileAtMem_radiate

  use SoilMoisture, only: &
       soilMoistureInit, &
       StoreNamelistFileAtSoilMoisture

  use Mem_turb, only: &
       AKMIN,         &
       if_urban_canopy, &
       StoreNamelistFileAtMem_turb

  use Mem_varinit, only: &
       vtime2,&
       htime2,&
       condtime2,&
       nud_type, &
       nud_cond, &
       StoreNamelistFileAtMem_varinit

  use Micphys, only: &
       level, &
       StoreNamelistFileAtMicphys, &
       mcphys_type, &
       icloud,idriz,ipris,idust,imd1flg,imd2flg

  use Shcu_vars_const, only: StoreNamelistFileAtShcu_vars_const

  use Mem_scalar, only: &
       scalar_g, &
       StoreNamelistFileAtMem_scalar

  use Teb_spm_start, only: &
       TEB_SPM, &
       StoreNamelistFileAtTeb_spm_start

  use mem_emiss, only: &
       ichemi,   &                    ! intent(inout)
       isource,  &                    ! intent(inout)
       ichemi_in, &                   ! intent(inout)
       StoreNamelistFileAtMem_emiss

  use Domain_decomp, only: StoreNamelistFileAtDomain_decomp
  use ModPostProcess, only: AllPostTypes, CreatePostProcess, &
       PostProcess, DestroyPostProcess
#ifdef cdf
  use MOdPostGridNetCDF, only: netCDFFirstTime
#endif
  !--(DMK-CCATT-INI)------------------------------------------------------------------
  USE monotonic_adv, ONLY: StoreNamelistFileAtRadvc_mnt
  !  USE newComm, ONLY: findAndFillGhostZone

  use ccatt_start, only: &
       StoreNamelistFileAtCCatt_start, &
       ccatt

  use Extras, only: &
       extra2d, &
       na_extra2d , &
       StoreNamelistFileAtExtras

  use mem_chem1, only: &
       chem1_src_g,                  & ! (INOUT) read_sourcemaps()
       chemistry,                    & ! (IN) initOneProc(), (IN) read_sourcemaps()
       recycle_tracers,              & ! (IN) initOneProc() (conflito com 'use mem_scalar')
       nsrc,                         & ! (IN) read_sourcemaps()
       nvert_src=>chem1_src_z_dim_g, & ! (IN) read_sourcemaps()
       bburn, antro, bioge,  geoge,  & ! (IN) read_sourcemaps()
       src_name,                     & ! (IN) read_sourcemaps()
       ntimes_src,                   & ! (IN) read_sourcemaps()
       diur_cycle, &                   ! (IN) read_sourcemaps()
       StoreNamelistFileAtMem_chem1, &
       chemistry, &
       chem_assim

  use mem_aer1, only: aerosol, StoreNamelistFileAtMem_aer1
  use mem_chem1aq, only: StoreNamelistFileAtMem_chem1aq
  use mem_plume_chem1, only: StoreNamelistFileAtMem_plumeChem1
  use mem_volc_chem1, only: StoreNamelistFileAtMem_volcChem1
  use chem_sources, only: &
       emiss_cycle_time,  &
       srcmapfn, &
       read_sourcemaps, &
       StoreNamelistFileAtChemSources

  use mem_stilt, only: StoreNamelistFileAtMem_stilt
  !lfr - tuv
  use tuvParameter, only: wstart,wstop,nwint,listFiles
  use ModTuv, only: initTuv,sw,ns,slabel,ks,kw,wc,wbioStart,wBioEnd !subroutine
  use ModTuvDriver, only: InitTuvDriver

  use wind_Farm, only: StoreNamelistFileAtWindFarm, &
			output_windFarms, &
			windfarm

  !lfr - tuv
  !--(DMK-CCATT-OLD)-----------------------------------------------------
  !  use Catt_start, only: StoreNamelistFileAtCatt_start
  !--(DMK-CCATT-FIM)------------------------------------------------------------------

  use cuparm_grell3, only: &
       init_weights, &
       StoreNamelistFileAtCup_grell3

  use ModNamelistFile, only : &
       namelistFile, &
       CreateNamelistFile, &
       DestroyNamelistFile, &
       GetNamelistFileName, &
       ReadNamelistFile, &
       BroadcastNamelistFile, &
       DumpNamelistFile, &
       TimeUnitsToSeconds

  use io_params, only: afilout

  ! For specific optimization depending the type of machine
  use machine_arq, only: &
       machine ! INTENT(IN)

  use ParLib, only: &
       parf_barrier    ! DEBUG

  use ModTimeStamp, only: SynchronizedTimeStamp ! DEBUG

  use advect_kit, only :   &
       advect_first_alloc, &  ! Subroutine
       prepare_inv            ! Subroutine

  use grid_dims, only : &
       nzpmax, &          ! (IN) read_sourcemaps()
       maxgrds            ! INTENT(IN)

  use io_params, only : & ! Include by Alvaro L.Fazenda
       createIoData,    & ! Subroutine
       destroyIoData,   & ! Subroutine
       IOUTPUT,         & !INTENT(IN)
       IPOS,            & !INTENT(IN)
       prtcputime,      & ! intent(in)
       avgtim,          & !INTENT(IN)
       frqanl,          & !INTENT(IN)
       frqlite,         & !INTENT(IN)
       frqmean,         & !INTENT(IN)
       frqhis,          & !INTENT(IN)
       frqboth            !INTENT(IN)

  use grid_dims, only: &
       maxgrds,        &
       maxmach,        &
       BRAMS_version


  use mem_grid, only:  &
       iyear1,         & ! (IN) read_sourcemaps()
       imonth1,        & ! (IN) read_sourcemaps()
       idate1,         & ! (IN) read_sourcemaps()
       itime1,         & ! (IN) read_sourcemaps()
       dzt,            & ! (IN) read_sourcemaps()
       zm,             & ! (IN) read_sourcemaps()
       zt,             & ! (IN) read_sourcemaps()
       get_akmin2d, &
       akminvar,         & ! intent(inout)
       iversion,      & ! intent(out)
       initial,       &
       nzs,           &
       ngridsh,       &
       StoreNamelistFileAtMem_grid, &
       timeunit,       &
       timmax,         &
       npatch,         & ! INTENT(IN)
       createMemGrid,  & !Subroutine
       destroyMemGrid, & !Subroutine
       allocAkmin2d,   & !Subroutine
       istp,           &
       grid_g,         &
       ngbegun,        & ! INTENT(OUT)
       nnxp,           &
       nnyp,           &
       nzg,            & ! INTENT(IN)
       f_thermo_e,     & ! INTENT(OUT)
       f_thermo_w,     & ! INTENT(OUT)
       f_thermo_s,     & ! INTENT(OUT)
       f_thermo_n,     & ! INTENT(OUT)
       dtlt,           & ! INTENT(IN)
       isched,         & ! INTENT(IN) - a ser (OUT) na chamada a modsched
       isstp,          & ! INTENT(OUT)
       dtlongn,        & ! INTENT(IN) - a ser Modifyed in DTSET
       iflag,          & ! INTENT(IN) - ALF - Modifyed in DTSET
       ideltat,        & ! INTENT(IN)
       runtype,        &
       dyncore_flag,   &
       expnme,         &
       ngrids,         &
       ngrid,          & ! INTENT(OUT)
       nzp,            &
       nxp,            &
       nyp,            &
       nnzp,           &
       isched,         &
       maxsched,       &
       maxschent,      &
       nxtnest,        &
       nndtrat,        &
       nsubs,          &
       cflxy,          & ! intent(out)
       cflz,           & ! intent(out)
       time,           &
       timmax,	       &
       begtime,        &
       order_h

  use ModParallelEnvironment, only : &
       ParallelEnvironment,    & ! type
       CreateParallelEnvironment, & ! Subroutine
       MsgDump, & ! Subroutine
       DestroyParallelEnvironment   ! Subroutine

  use node_mod, only:  &
       mynum, &
       mchnum, &
       master_num, &
       mxp,  &
       myp,  &
       mzp,  &
       nodemzp, &
       nodemxp, &
       nodemyp, &
       nodei0, &
       nodej0, &
       StoreNamelistFileAtNode_mod, &
       StoreDomainDecompAtNode_mod, &
       StoreParallelEnvironmentAtNode_mod, &
       alloc_bounds,   & ! Subroutine
       dealloc_bounds, & ! Subroutine
       ibcon,          & ! INTENT(IN)
       ia, iz, izu,    & ! INTENT(IN)
       ja, jz, jzv,    & ! INTENT(IN)
       i0, j0,         & ! INTENT(IN)
       nodemxp,        & ! INTENT(IN)
       nodemyp,        & ! INTENT(IN)
       load_bal,       & ! INTENT(IN)
       nodei0,         & ! INTENT(IN)
       nodej0,         & ! INTENT(IN)

!--(DMK-CCATT-INI)-----------------------------------------------------------
       nodeia, &
       nodeiz, &
       nodeja, &
       nodejz, &
!--(DMK-CCATT-FIM)-----------------------------------------------------------

       nmachs,         & ! INTENT(IN)
       mynum,          & ! INTENT(IN)
       mxp,            &
       myp,            &
       mzp,            &
       mchnum,         &
       master_num,     &

!--(DMK-CCATT-INI)-----------------------------------------------------------
       ixb, ixe, iyb, iye, & !To advect_mnt
!--(DMK-CCATT-FIM)-----------------------------------------------------------

       ProcessOrder      ! procedure

  use dtset, only: &
       SetDt,      &
       dtset_new,  &
       dtSet_firstTime, &
       sscourn, &
       maxCflPercent, &
       dump_courn    ! subroutine

  use ReadBcst, only: &
       MaxCFLOverall

  use mem_scratch, only : &
       scratch, &
       destroyVctr ! subroutine

  use ref_sounding, only : &
       createRefSounding,  & ! subroutine
       destroyRefSounding, &! subroutine
       StoreNamelistFileAtRef_sounding

  use mem_scratch, only: &
       createvctr, &
       destroyvctr

  use digitalFilter, only: &
       StoreNamelistFileAtdigitalFilter,&
       applyDF,          &
       timeWindowDF,     &
       frqanlDF, iposDF

  !--(DMK-CCATT-INI)------------------------------------------------------------------
  USE AdvectData, ONLY: InitAdvect
  USE modComm, ONLY: initExtraComm

  USE monotonic_adv, ONLY: advmnt, &
       GhostZoneLength

  ! OBS: MODULOS NECESSARIOS PARA LEITURA DE EMISSAO
  !-----------------------------------------------------------------------------------
  use module_dry_dep, only: dep_init        ! Subroutine

  use chem_sources, only:   read_sourcemaps ! Subroutine

  use chem1_list, only:           &
       dvj,                       & ! (IN) dep_init()
       chem_nspecies=>nspecies,   & ! (IN) dep_init(),(IN) read_sourcemaps()
       spc_chem_alloc=>spc_alloc, & ! (IN) read_sourcemaps()
       src,                       & ! (IN) read_sourcemaps()
       off,                       & ! (IN) read_sourcemaps()
       spc_chem_name =>spc_name,  & ! (IN) read_sourcemaps()
       on,                        & ! (IN) read_sourcemaps()
       chemical_mechanism,        & ! (IN) read_sourcemaps()
       emiss_ajust,               & ! (IN) read_sourcemaps()
       transport,                 & ! (IN) read_sourcemaps()
       co,                        &
       co2,                       &
       PhotojMethod,              &
       chemical_mechanism



  use aer1_list, only:           &
       aer_nspecies=>nspecies,   & ! (IN) read_sourcemaps()
       spc_aer_alloc=>spc_alloc, & ! (IN) read_sourcemaps()
       spc_aer_name =>aer_name,  & ! (IN) read_sourcemaps()
       nmodes,                   & ! (IN) read_sourcemaps()
       mass_bin_dist,            &
       aerosol_mechanism,        &
       emiss_ajust_aer

!- for matrix
!       ,v_ash,                     &! (IN) read_sourcemaps()
!       urban, nucle, accum,      & ! (IN) read_sourcemaps()
!       ,aer_bburn => bburn  ! (IN) read_sourcemaps()

  use mem_aer1, only:  &
       aer1_g  ! (INOUT) read_sourcemaps()

  use mem_plume_chem1, only: &
       plume_g,              & ! (INOUT) read_sourcemaps()
       plume_mean_g,         & ! (INOUT) read_sourcemaps()
       plumerise,            & ! (IN) read_sourcemaps()
       nveg_agreg,           & ! (IN) read_sourcemaps()
       tropical_forest,      & ! (IN) read_sourcemaps()
       boreal_forest,        & ! (IN) read_sourcemaps()
       savannah,             & ! (IN) read_sourcemaps()
       grassland,            & ! (IN) read_sourcemaps()
       plume_fre_g             ! (IN) read_sourcemaps()

  use mem_stilt, only: &
       stilt_g,        & ! (OUT) initOneProc()
       iexev             ! (IN) initOneProc()

  use mem_volc_chem1, only: &
       volc_mean_g,         &  ! (INOUT) read_sourcemaps()
       volcanoes               ! (IN) read_sourcemaps()

  use mem_basic, only: &
       basic_g

  use grid_dims, only: &
       nzpmax            ! (IN) read_sourcemaps()
  !--(DMK-CCATT-FIM)------------------------------------------------------------------

  use ModGridTree, only: &
       GridTree, &
       CreateGridTree, &
       DumpGridTree, &
       GridTreeRoot, &
       NextOnGridTree, &
       DestroyGridTree

  use ModGrid, only: &
       Grid, &
       DumpGrid, &
       InsertMessagePassingAtOneGrid


  use meteogram, only:              &
      InitMeteogram,                &
      ProcessLocalMeteogram,        &
      StoreNamelistFileAtmeteogram, &
      applyMeteogram,               &
      meteogramFreq,                &
      meteogramMap

  use rams_microphysics_2m, only :  &
        jnmbinit_2M=>jnmbinit       &
       ,initqin_2M  =>initqin  	    &
       ,initqin2_2M =>initqin2 	    &
       ,initqin3_2M =>initqin3 	    &
       ,initqin4_2M =>initqin4 	    &
       ,initqin5_2M =>initqin5 	    &
       ,micro_master_2M =>micro_master

   use mem_carma, only: &
       read_aotMap, &
       StoreNamelistFileAtmem_carma, &
       carma_aotMap

   use dump

   use dam, only: &
       damModule, &
       initDams, &
       acumPrecipInDam, &
       outputDamPrecip, &
       StoreNamelistFileAtDams

  use aerClimMod, only: &
      gradsRead

  use initMicThompson, only: &
    readDataFriendly, &
    adJustFriendlyForMonth

  use ModEvaluation, only: &
    StoreNamelistFileAtEvaluate, &
    RMSE_average, &
    evaluate, &
    timeCount
 
 use modIau,only: &
    StoreNamelistFileAtIAU, &
    initComIau, &
    applyIAU

 use modTimeLineFRN, only: &
    readSites, &
    createSitesFile

 use mem_sfire, only: &
    StoreNamelistFileAtMem_sfire


  implicit none
  private
  public :: OneProc

contains

  ! ****************************************************************************
  !  OneProc: Initialization whenever all processes compute; master also does io
  !

  subroutine OneProc(nmachs_in, mchnum_in, master_num_in)

    !MB: for testing only
    use ModTimestep_RK, only: flag_mb_adv_test
    use dump, only: &
      dumpMessage !dump function

    !If compiler INTEL sometimes is necessary unconmment the line bellow:
    !USE IFPORT
    !

    include "i8.h"
    include "files.h"
    include "tsNames.h"
    include "mpif.h"
    include "constants.f90"

    ! Arguments:
    integer, intent(in) :: nmachs_in           ! number of processes (0 iff sequential run)
    integer, intent(in) :: master_num_in        ! this process rank (0:nmachs_in-1); 0 on sequential runs
    integer, intent(in) :: mchnum_in            ! this process rank (0:nmachs_in-1); 0 on sequential runs
    ! Local variables:
    integer :: isendflg, isendlite, isendmean, isendboth, nt, npass, &
               icm, ifm, nfeed
    integer :: isendbackflg ! ALF - For local processing
    integer :: isendiv, isendsst, isendndvi
    integer :: isendsrc

    logical :: instFlag, liteFlag, meanFlag, histFlag, posFlag
    integer :: ng
    integer :: i_xyz
    integer(kind=i8) :: l_xyz
    character(len=12) :: c0
    character(len=12) :: c1
    integer :: i,nndtflg
    integer :: j,k
    real :: w2, w3, w4, wtime_start, wtime_end  ! wall time !w1
    real, external :: walltime    ! wall time
    real :: t1, t2, t6, totcpu       ! cpu time
    real, allocatable :: dxtmax_local(:)
    integer :: ierr
    integer :: nn2, nn3
    logical :: AKMIN_ALLOC
    character(len=f_name_length) :: namelistFileName ! namelist file name
    character(len=*), parameter :: h="**(OneProc)**"
    character(len=*), parameter :: header="**(OneProc)**"
    character(len=*), parameter :: version="5.4"

    type(parallelEnvironment), pointer :: oneParallelEnvironment => null()
    type(namelistFile), pointer :: oneNamelistFile => null()
    type(GridTree), pointer :: AllGrids => null()
    type(GridTree), pointer :: OneGridTreeNode => null()
    type(Grid), pointer :: OneGrid => null()
    type(AllPostTypes), pointer :: oneAllPostTypes => null()

    logical, parameter :: dumpLocal=.false.

    logical :: dirExist
    character(len=255) :: tmpdir

!XXXsrf    integer :: iau_phase

    dtSet_firstTime=.true.

    !Open the file to log the run
    ierr=openLogFile()

    open(unit=66,file='jules.log',status='replace',action='write')
    write(unit=66,fmt='(A)') 'BRAMS - JULES - LOG'
    close(unit=66)
    open(unit=66,file='brams.log',status='replace',action='write')
    write(unit=66,fmt='(A)') 'BRAMS - LOG'
    close(unit=66)

    !MB: only for testing:
    !flag_mb_adv_test = .FALSE.

    ! store MPI rank, size and master at mpi/node_mod.f90

    call CreateParallelEnvironment(nmachs_in, mchnum_in, master_num_in, &
         MPI_COMM_WORLD, oneParallelEnvironment)
    call StoreParallelEnvironmentAtNode_mod(oneParallelEnvironment)

    ! wall time at the beginning of execution

    wtime_start = walltime()

    ! create namelistFile object
    call CreateNamelistFile(oneNamelistFile)

    ! master gets namelist file name
    ! master reads file
    ! master convert namelist time units into seconds

    if (oneParallelEnvironment%mchnum==oneParallelEnvironment%master_num) then
       call GetNamelistFileName(oneNamelistFile)
       call ReadNamelistFile(oneNamelistFile)
       call TimeUnitsToSeconds(oneNamelistFile)
       call DumpNamelistFile(oneNamelistFile,oneParallelEnvironment%nmachs,oneParallelEnvironment%mchnum &
                         ,oneParallelEnvironment%master_num)
    end if



    if(trim(oneNamelistFile%runtype)=='MAKESFC' .and. oneParallelEnvironment%nmachs>1)&
       iErrNumber=dumpMessage(c_tty,c_yes,header,c_modelVersion,c_fatal,'For run type '// &
                  trim(oneNamelistFile%runtype)//' the number of processes must be 1!')

    ! Broadcast namelist
    call BroadcastNamelistFile(oneNamelistFile, oneParallelEnvironment)

    call StoreNamelistFileAtIo_Params(oneNamelistFile)
    call StoreNamelistFileAtIsan_coms(oneNamelistFile)
    call StoreNamelistFileAtMem_cuparm(oneNamelistFile)
    call StoreNamelistFileAtCup_grell3(oneNameListFile)
    call StoreNamelistFileAtMem_globrad(oneNamelistFile)
    call StoreNamelistFileAtMem_grell_param(oneNamelistFile)
    call StoreNamelistFileAtMem_grid(oneNamelistFile)
    call StoreNamelistFileAtMem_leaf(oneNamelistFile)
    call StoreNamelistFileAtMem_oda(oneNamelistFile)
    call StoreNamelistFileAtMem_radiate(oneNamelistFile)
    call StoreNamelistFileAtMem_sfire(oneNamelistFile)
    call StoreNamelistFileAtSoilMoisture(oneNamelistFile)
    call StoreNamelistFileAtMem_turb(oneNamelistFile)
    call StoreNamelistFileAtMem_varinit(oneNamelistFile)
    call StoreNamelistFileAtMicphys(oneNamelistFile)
    call StoreNamelistFileAtNode_mod(oneNamelistFile)
    call StoreNamelistFileAtRef_sounding(oneNamelistFile)
    call StoreNamelistFileAtShcu_vars_const(oneNamelistFile)

    !--(DMK-CCATT-INI)------------------------------------------------------------------
    !  call StoreNamelistFileAtCatt_start(oneNamelistFile)
    !--(DMK-CCATT-FIM)------------------------------------------------------------------

    call StoreNamelistFileAtMem_scalar(oneNamelistFile)

    !--(DMK-CCATT-INI)------------------------------------------------------------------
    call StoreNamelistFileAtExtras(oneNamelistFile)
    call StoreNamelistFileAtCCatt_start(oneNamelistFile)
    call StoreNamelistFileAtMem_chem1(oneNamelistFile)
    call StoreNamelistFileAtMem_aer1(oneNamelistFile)
    call StoreNamelistFileAtMem_chem1aq(oneNamelistFile)
    call StoreNamelistFileAtChemSources(oneNamelistFile)
    call StoreNamelistFileAtMem_plumeChem1(oneNamelistFile)
    call StoreNamelistFileAtMem_volcChem1(oneNamelistFile)
    call StoreNamelistFileAtMem_stilt(oneNamelistFile)
    !--(DMK-CCATT-FIM)------------------------------------------------------------------

    call StoreNamelistFileAtTeb_spm_start(oneNamelistFile)
    call StoreNamelistFileAtMem_emiss(oneNamelistFile)
    call StoreNamelistFileAtTeb_vars_const(oneNamelistFile)
    call StoreNamelistFileAtDomain_decomp(oneNamelistFile)

    !--(DMK-CCATT-INI)-----------------------------------------------------------
    CALL StoreNamelistFileAtradvc_mnt(oneNamelistFile)
    !--(DMK-CCATT-FIM)-----------------------------------------------------------

    call StoreNamelistFileAtdigitalFilter(oneNamelistFile)

    call StoreNamelistFileAtmeteogram(oneNamelistFile)
    call StoreNamelistFileAtmem_carma(oneNamelistFile)
    call StoreNamelistFileAtWindFarm(oneNamelistFile)
    call StoreNamelistFileAtDams(oneNamelistFile)
    call StoreNamelistFileAtEvaluate(oneNamelistFile)

    call StoreNamelistFileAtIAU(oneNamelistFile)
    ! build and dump all grids
    call CreateGridTree(oneNamelistFile, oneParallelEnvironment, AllGrids)

    ! Allocating dxtmax_local
    allocate(dxtmax_local(ngrids), STAT=ierr)
    !if (ierr/=0)  call fatal_error("ERROR allocating dxtmax_local (oneproc)",header,version)
    if (ierr/=0) iErrNumber=dumpMessage(c_tty,c_yes,header,modelVersion,c_fatal &
        ,'ERROR allocating dxtmax_local (oneproc)')

    ! Allocating "mem_grid" data
    call createMemGrid(ngrids, nnxp, nnyp, nnzp)

    ! Allocating bounds
    call alloc_bounds(ngrids, nmachs)

    call StoreDomainDecompAtNode_mod(AllGrids)

    ! Allocating IO Data
    call createIoData(ngrids)

    ! Allocating Sounding Data
    call createRefSounding(ngrids, nnzp)

    ! create process numbering (basic initialization of module node_mod @mpi/node_mod.f90)
    call ProcessOrder()

    ! Various option settings that should normally not be changed
    ! (previously namelist parameters)
    !**(JP)** is this really necessary?

    call eng_params()

    ! First check of options, mainly for numbers of grid points

    call opspec1(nmachs, mchnum, master_num)

    ! Basic grid coordinate setup for statically allocated data structures:
    ! number of grid points, deltas, coordinate and nesting coefficients,
    ! all stored at mem_grid

    call createVctr(ngrids, nnxp, nnyp, nnzp)
    call GridSetup(1)
    call destroyVctr()

    ! Additional checks, mainly for nesting

    call opspec2()

    ! Parallel domain decomposition:
    ! compute position of all sub-domains on the global grid (global indices ixb, ixe, iyb, iye at node_mod);
    ! include ghost zone (global indices nxbeg, nxend, nybeg, nyend at node_mod);
    ! verify if any sub-domain boundary is also a global domain boundary (ibcon at node_mod);
    ! compute first and last local index of each sub-domain without ghost zone (nodeia, nodeiz, nodeja, nodejz at node_mod);
    !
    ! The execution of decomp_node was antecipated (replacing the execution of node_decomp during initialization by the master)
    ! since memory has to be allocated for any process' grid before initialization.
    ! In the master-slave execution case, memory was allocated by the master for the full grid (by rams_mem_alloc(0)) prior to compute
    ! domain decomposition. After computing domain decomposition, sub-domain size and position were sent to each slave,
    ! that allocates memory based on the received data and then initializes.
    !
    ! Procedure decomp_node does not fully replace node_decomp, since parts of node_decomp
    ! (par_node_paths and following code) required memory allocation (not done yet).
    ! Remaining code was moved to procedure NodePathsBuffAlloc.
    call decomp_node(1)

    ! master dumps domain decomposition at stdout and at selected file

    if (mchnum==master_num) then
       call domain_decomposition_dump(OUTPUT_UNIT)
    end if

    ! Check sfc,sst,ndvi files; remake if needed
    call MakeSfcfiles()

    
    ! Behave accordingly to run typ
!================================================================================================
!================================================================================================

    if (runtype(1:7)=='MAKESFC') then

       ! done on a MAKESFC run
       if (mchnum==master_num) then
          write(OUTPUT_UNIT,"(a)") ' MAKESFC run complete'
       end if

!================================================================================================
!================================================================================================

    else if (runtype(1:9)=='MAKEVFILE') then
#ifndef netcdf 
      if(ICFILETYPE==2 .or. ICFILETYPE==3) &
          iErrNumber=dumpMessage(c_tty,c_yes,header,version,c_fatal," To use NetCDF (GEOS/ECMWF,etc) the code must be compiled with NetCDF!!! ")
      if(IPOS==3) &
          iErrNumber=dumpMessage(c_tty,c_yes,header,version,c_fatal," To use NetCDF output the code must be compiled with NetCDF!!! ")
#endif
   time  = 0.
   ! on a "MAKEVFILE" run, call ISAN, then exit.
      !if (mchnum==master_num) then
        call chem_isan_driver(namelistFileName)
        if(ccatt==1 .and. chem_assim==1 .and. chemistry >= 0)then
             iErrNumber=dumpMessage(c_tty,c_yes,header,version,c_notice," CHEM_ISAN complete ")
	      else
             iErrNumber=dumpMessage(c_tty,c_yes,header,version,c_notice," ISAN complete ")
	      endif
      !endif
!================================================================================================
!================================================================================================

    else  ! for RUNTYPE = INITIAL/HISTORY (endif @ line ~1536)
    
     
#ifdef cdf
      !To open netCDF output (if used) just once
      netCDFFirstTime=.true.
#endif

      ! If we got here, we are doing an actual simulation (INITIAL)
      ! Initialize micro arrays. May need to change some settings which affect memory.

      ! initiate tuv tables
      if (trim(PhotojMethod) == 'FAST-TUV' .and. chemistry >= 0) then
          CALL initTuv(wstart,wstop,nwint,listFiles,mchnum,chemical_mechanism)
       !	  TODO - to review code below for bio vars, probably just debug stuff
       !          write(77,*) wbioStart,wBioEnd
       !          do i = wbioStart,wBioEnd
       !            write (77,fmt='(I3.3,1X,A)') i,slabel(i)
       !            do j=1,kw
       !              write (77,fmt='(I3.3,F10.4,1X,F10.4)') j,sw(i,j),wc(j)
       !            end do
       !          end do
          CALL InitTuvDriver()
      endif

      if(mcphys_type==0) then
            CALL jnmbinit()

       elseif(mcphys_type==1) then
            CALL jnmbinit_2M()

       endif

       ! Allocate memory for this process sub-domain only
       !**(JP)** This should allocate memory for all modules (to be certified!!!)
       call rams_mem_alloc(2)
       
       if(ioutput == 5)then
       	call setInitial4Vtable(1)
	     ioutput = 2
       end if

       ! Allocate AKMIN2D if necessary
       AKMIN_ALLOC = .false.
       do ifm=1,ngrids
          if (AKMIN(ifm)<0.) AKMIN_ALLOC = .true.
       end do
       if (AKMIN_ALLOC) then
          call allocAkmin2d(ngrids, nodemxp(mynum,1:ngrids), nodemyp(mynum,1:ngrids))
       endif

       ! Communication paths, sizes and buffers
       !
       ! Procedure NodePathsBuffAlloc finishes replacing node_decomp, since memory was already allocated
       ! (does the part from call to par_node_paths to the end). It computes values of node_mod module.
       !

       call NodePathsBuffAlloc()

       ! build message passing
       !print *,'LFR-DBG inside oneproc 4: ',nodemxp(mynum,1); call flush(6)
       OneGridTreeNode => GridTreeRoot(AllGrids)
       do while (associated(OneGridTreeNode))
          call InsertMessagePassingAtOneGrid(OneGridTreeNode%curr)
          OneGridTreeNode => NextOnGridTree(OneGridTreeNode)
       end do

       if (dumpLocal) then
          OneGridTreeNode => GridTreeRoot(AllGrids)
          OneGrid => OneGridTreeNode%curr
          call DumpGrid(OneGrid)
       end if


!======XXXXXsrf
!    IAU_phase = 0
!    do while (IAU_phase <= applyIAU)
!     IAU_phase = IAU_phase + 1
!======XXXXXsrf
       !**(JP)** Code brought from old "rams_node" initialization. Its execution
       ! was antecipated to allow message passing during initialization, required
       ! for binary reproducibility (see subroutine FilDn0uv)
       !call dumpVarAllLatLonk(leaf_g(1)%patch_area,'patch_area',846,0,0,1,nodemxp(mynum,1),1,nodemyp(mynum,1),1,npatch,0.0,600.0,h) !
       call InitFields(1)

       ! initialization driver
       call initOneProc(AllGrids, namelistFileName)

       ! Compute Courant numbers cflxy and cflz, get maximum over all processes  and dump
       do ifm=1,ngrids
          call newgrid(ifm)
          call cfl(mzp, mxp, myp, 0, 0)
       enddo
       call MaxCFLOverall(cflxy, cflz)

       ! Initialize dtlongn, nndtrat, and nnacoust, and compute the timestep
       ! schedule for all grid operations.
       call SetDt(mynum, nndtflg, dxtmax_local)
       if (mchnum==master_num) then
          call dump_dtset(nndtflg)
       end if


       call modsched(isched, maxsched, maxschent, ngrids, nxtnest, &
            nndtrat, nsubs)
       if (mchnum==master_num) then
          call dump_modsched(isched, maxsched, maxschent, ngrids, &
               nxtnest, nndtrat, nsubs)
          call dump_courn()
       end if

       !=========================================================================

       !**(JP)** old subroutine rams_node starts

       ! finish initialization

       isendflg     = 0
       isendlite    = 0
       isendmean    = 0
       isendboth    = 0
       isendbackflg = 0 
       isendiv      = 0 
       isendsst     = 0 
       isendndvi    = 0 
       isendsrc     = 0 

       ! correctness of local value of x*y*z

       i_xyz = maxval(nodemxp(mynum,1:ngrids))* &
               maxval(nodemyp(mynum,1:ngrids))*maxval(nnzp(1:ngrids))
       l_xyz = maxval(nodemxp(mynum,1:ngrids))* &
               maxval(nodemyp(mynum,1:ngrids))*maxval(nnzp(1:ngrids))
       if (i_xyz/=l_xyz) then
          print *, "NODE:", mynum
          print *, "Checking number of processors and total of points"
          print *, "Integer Total(x*y*z) =", i_xyz, &
               " - Long Total(x*y*z) =", l_xyz
          !call fatal_error("Integer - Long not equal to Zero! Use more process",header,version)
          iErrNumber=dumpMessage(c_tty,c_yes,header,modelVersion,c_fatal, &
                    'Integer - Long not equal to Zero! Use more process')
       end if
      !print *,'lfr-dbg:', nodemxp(mynum,1:ngrids),nodemyp(mynum,1:ngrids),1+2*order_h,dyncore_flag; call flush(6)
      if(dyncore_flag==2 .and. nodemxp(mynum,1)<1+2*order_h) then
        write (*,fmt='(A)') '******** WARNING *******'
        write (*,fmt='(A,I4.4,A,I3,A,I6)') 'Number of columns to be process in X is ',nodemxp(mynum,1),' that is less than ',&
        1+2*order_h,' for processor #',mynum
      endif
      if(dyncore_flag==2 .and. nodemyp(mynum,1)<1+2*order_h) then
        write (*,fmt='(A)') '******** WARNING *******'
        write (*,fmt='(A,I4.4,A,I3,A,I6)') 'Number of columns to be process in Y is ',nodemxp(mynum,1),' that is less than ',& 
        1+2*order_h,' for processor #',mynum
      endif



       ! Calculating DXTMAX
       nt = 0
       do ifm=1,ngrids
          nn2 = nodemxp(mynum,ifm)
          nn3 = nodemyp(mynum,ifm)
          dxtmax_local(ifm) = max(grid_g(ifm)%dxt(1,1), grid_g(ifm)%dxt(nn2,1), &
               grid_g(ifm)%dxt(nn2,nn3), grid_g(ifm)%dxt(1,nn3))
       enddo

       !tst LFR
       Call initExtraComm(nmachs,mynum,GhostZoneLength,nnxp,nnyp,nnzp,ixb,ixe,iyb,iye,master_num, &
                   nodei0,nodej0,nodemxp,nodemyp,nodemzp)
       !CALL InitComm(ngrids,nmachs,mynum,GhostZoneLength,nnxp,nnyp,nnzp,ixb,ixe,iyb,iye)
       !end tst
       
       IF(advmnt>=1) then
       !- monotonic advection
          CALL InitAdvect(ngrids,nmachs,mynum,GhostZoneLength,nnxp,nnyp,nnzp,ixb,ixe,iyb,iye)
       ENDIF

       if(damModule==1) then
         call initDams(nodemxp(mynum,1),nodemyp(mynum,1) &
                   ,grid_g(1)%glat,grid_g(1)%glon,mchnum,master_num)
       endif

      if(mcphys_type==3) call readDataFriendly()

      if (aerosol==-1 .and. .not. (CCATT==1 .and. chemistry >= 1)) &
       call gradsRead('./tables/aerClim/','aerosols.gra',grid_g(1)%glat,grid_g(1)%glon)

       IF(machine == 1) then
         ! Allocate and initialize data for new advection scheme if used
         ! Actually it is only used if machine=1, defining NEC-SC system
         call advect_first_alloc(ngrids, nnzp(1:ngrids), &
                        nodemxp(mynum,1:ngrids), nodemyp(mynum,1:ngrids))
         call prepare_inv(ngrids)
       ENDIF

       ! Checking if the actual node have to run thermo on the boundaries

       f_thermo_e = .false.
       f_thermo_w = .false.
       f_thermo_n = .false.
       f_thermo_s = .false.
       do ng=1,ngrids
          call newgrid(ng)
          call node_index()
          if (i0==0)         f_thermo_e(ng) = .true.   ! east  bondary
          if ((mxp+i0)==nxp) f_thermo_w(ng) = .true.   ! west  boundary
          if (j0==0)         f_thermo_s(ng) = .true.   ! south boundary
          if ((myp+j0)==nyp) f_thermo_n(ng) = .true.   ! north boundary
       end do

       if (IPOS/=0) then
          ! create post processing
          oneAllPostTypes => null()
          call CreatePostProcess(oneNamelistFile, oneAllPostTypes)
       endif
       
       !LFR: Fazendo a inicialização dos sites para séries temporais
       !print *,mchnum,master_num,'Chamando....'
       !if(mchnum == master_num) then
         !print *,'Chamando readsites ::::::::'
       if(IPOS == 10 .or. IPOS == 11) call readSites(mchnum,master_num,oneNamelistFile)
       !endif


       select case (IPOS)

      !!$     case (1)
      !!$        ! post process initial state of the atmosphere in HDF5
      !!$        call PostProcessHDF5(oneNamelistFile, oneAllPostTypes)

       case (2,3,10,11)
          ! post process initial state of the atmosphere in Grads
          call PostProcess(AllGrids, oneNamelistFile, oneAllPostTypes)
       end select

       ! wall time at the end of initialization

       call parf_barrier(777)
       if (mchnum==master_num) then
          w2=walltime()
          write(c1,"(f12.2)") w2-wtime_start
          write(*,fmt='(A)') c_empty
#ifdef color
  iErrNumber=dumpMessage(c_tty,c_yes,'','',c_notice," === Finish"//&
     " initialization; Wall(sec) = "//c_purple//trim(adjustl(c1))//&
     c_peach//" [s]"//c_noColor)
#else
  iErrNumber=dumpMessage(c_tty,c_no,'','',c_notice," === Finish"//&
     " initialization; Wall(sec) = "//trim(adjustl(c1))//&
     " [s]")
#endif
          write(*,fmt='(A)') c_empty
          !write(OUTPUT_UNIT,"(/,a,/)") " === Finish initialization; Wall(sec)="// &
          !     trim(adjustl(c1))
       end if


       !--Digital Filter ------------------------------------------------------
       IF(applyDF)then
          frqanldf=frqanl
          frqanl= 0.
          iposDF=ipos
          ipos=0
       ENDIF
       ! force first grid on grid tree

       OneGridTreeNode => GridTreeRoot(AllGrids)
       OneGrid => OneGridTreeNode%curr
       

       IF(applyMeteogram) call InitMeteogram(oneGrid%meteoPolygons, oneGrid%id, trim(meteogramMap))

       ! loop over time, advancing all grids one long timestep forward


       !MB>>
       !if ( flag_mb_adv_test ) then
       !  !basic_g(ngrid)%uc(:,:,:) = 267.0  ! for horizontal advection test
       ! basic_g(ngrid)%uc(:,:,:) =   0.0
       ! basic_g(ngrid)%vc(:,:,:) =   0.0
       !  basic_g(ngrid)%wc(:,:,:) =   1.0   ! for vertical advection test
       !  basic_g(ngrid)%pc(:,:,:) = 0.0
       !  basic_g(ngrid)%thc(:,:,:) = 300.0
       !       do k=1, mzp
       !         do i=1, mxp
       !           do j=1, myp
       !                    ! test horizontal scalar advection:
       !                    ! basic_g(ngrid)%thc(k,i,j) = 300.0 + max( 0.0,  1.0 - abs( (i-10)/5.0 )  )
       !		            ! test vertical scalar advection:
       !	            basic_g(ngrid)%thc(k,i,j) = 300.0 + max( 0.0,  1.0 - abs( zt(k)-800.0)/500.0 )
       !             enddo
       !         enddo
       !     enddo
       !
       !end if
       !MB<<

       !Uncoment to calculate execution time and set noInstrumentation = false in ModTimestamp.f90
       ! call SynchronizedTimeStamp(TS_INIT) ! timestamp before: initialization

       ! loop over time, advancing all grids one long timestep forward

!---> here is the 'big' loop for the model time integration

       istp = 0

       do while (time<timmax)

          !print*,"time=====", time
       
          istp    = istp+1
          totcpu  = 0
          w3      = walltime()
          call timing(1,t1)
          nt      = nt + 1
          begtime = time


          ! if input time, get new fields from master

          if (isendbackflg==1) then
             if (isendiv==1) then

                !print*,"===> going to VarfReadStoreOwnChunk",time
                call VarfReadStoreOwnChunk(AllGrids, 2, nud_type)
                
		
             endif
             if (isendsst==1) then
                do ifm=1,ngrids
                   call SstReadStoreOwnChunk(3,ifm,ierr)
                enddo
             endif
             if (isendndvi==1) then
                do ifm=1,ngrids
                   call NdviReadStoreOwnChunk(3,ifm,ierr)
                enddo
             endif

             if (isendsrc==1) then
                do ifm = 1, ngrids
                   call newgrid(ifm)
                   !if(mynum==1) print*,'->1 chem/aer sources maps reading: ',time/3600.,' grid=',ifm

                   call read_sourcemaps(ifm,nnzp(ifm),nodemxp(mynum,ifm),nodemyp(mynum,ifm), &
                        ia,iz,ja,jz, time,iyear1,imonth1,idate1,itime1,ngrids,timmax,        &
                        chem_nspecies,spc_chem_alloc,src,off,nsrc,nvert_src,chem1_src_g, &
                        bburn,antro, bioge,  geoge,spc_chem_name,on,chemical_mechanism,  &
                        emiss_ajust,co,aer_nspecies,spc_aer_alloc,spc_aer_name,          &
                        src_name,chemistry,ntimes_src,aer1_g,nmodes,aerosol,plumerise,   &
                        nveg_agreg,plume_mean_g,nzpmax,dzt,grid_g(ifm)%rtgt,grid_g(ifm)%topt, &
                        transport,plume_g,tropical_forest,boreal_forest,savannah,         &
                        grassland,diur_cycle,volcanoes,volc_mean_g,basic_g(ifm)%dn0,zt,zm,&
                        mchnum, master_num,mass_bin_dist,CO2,ISFCL,aerosol_mechanism,     &
			plume_fre_g,emiss_ajust_aer)

                enddo
                !print*,'----------------------------------------------'
             endif

          end if

          ! If Generate only Pos-Proc read History file
          !if (IPOS/=0 .and. RUNTYPE=='POS') then
          !!$           call readFieldsHis(vtab_r, ngrids, WillGather, maxNFields)
          !endif

          !Uncoment to calculate execution time and set noInstrumentation = false in ModTimestamp.f90
          !   call SynchronizedTimeStamp(TS_INPUT) ! timestamp In

          ! compute output flags for this iteration and input flags for
          ! next iteration

          call comm_time(isendflg, isendlite, isendmean, isendboth, &
               isendbackflg, isendiv, isendsst, isendndvi, isendsrc)

          if(applyDF) then
             if( abs(time+dtlongn(1) - timeWindowDF)<0.0001 .or. & ! timemf>=vtime2 (=timewindowDF)
                  abs(time            - timeWindowDF)<0.0001    ) then
                isendbackflg=0
                frqanl=frqanldf
                ipos  =iposdf
             endif
          endif

          ! examine Courant numbers to verify if model needs to be stopped
          ! or (if ideltat < 0), to update dtlongn, nndtrat,
          ! nnacoust, sspct and isched.

          call dtset_new(mynum, nndtflg, dxtmax_local)

          ! CFL exausted: model will stop

          if (iflag>0) then
             isendflg  = 1
             isendlite = 1
             isendmean = 1
             isendboth = 1
          end if

          ! new long deltat required by CFL; re-schedule grids

          if (nndtflg>0) then
             call modsched(isched, maxsched, maxschent, ngrids, nxtnest, &
                  nndtrat, nsubs)
          end if

          ! loop through all grids and advance them in time by a dtlong

          if (RUNTYPE/='POS') then

             do npass=1,nsubs

                ! make scheduled grid the current grid

                isstp = isched(npass, 3)
                ngrid = isched(npass, 1)
                call newgrid(ngrid)
                call node_index()

                ! advance current grid forward by corresponding deltat

                time = begtime + (isched(npass,5)-1)*dtlt
                if(mcphys_type==3) call adJustFriendlyForMonth(time)

                !MB: callin the timestep routine
                !srf print*,"going to timestep ", time
                if ( ( dyncore_flag == 0 ) .or. ( dyncore_flag == 1 ) ) then
                  ! Leapfrog/forward-time based scheme
                  call timestep(OneGrid,oneNamelistFile)
                else if ( dyncore_flag == 2 ) then
                  ! Runge-Kutta based scheme
                  call timestep_rk(OneGrid,oneNamelistFile)
                else if ( dyncore_flag == 3 ) then
                  ! ABM3 based scheme
                  call timestep_abm(OneGrid,oneNamelistFile)
                else
                  !call fatal_error("ERROR in subroutine OneProc: value of dyncore_flag is unknown",header,version)
                  iErrNumber=dumpMessage(c_tty,c_yes,header,modelVersion,c_fatal &
                         ,'ERROR in subroutine OneProc: value of dyncore_flag is unknown')
                end if

                ngbegun(ngrid) = 1

                ! whenever required, send the current grid domain points to
                ! the nested grid to interpolate a nested grid's boundaries

                if (isched(npass,2)/=0) then

                   !**(JP)** not converted yet

                   call fatal_error(h//" multiple grids not converted yet")
                   icm = ngrid
                   do ifm=1,ngrids
                      if (nxtnest(ifm)==icm) then
                         ngrid = ifm            ! JP: is this necessary?
                         isstp = isched(npass,3)! JP: is this necessary (maybe incorrect)!
                         call newgrid(ifm)
                         call node_sendnbc(ifm, icm)
                         call node_getnbc(ifm, icm)
                      end if
                   end do
                end if

                ! whenever required, send feedback fields from current grids to
                ! coarser grids

                if (isched(npass,4)/=0) then

                   !**(JP)** not converted yet

                   !call fatal_error(h//" multiple grids not converted yet")
                   iErrNumber=dumpMessage(c_tty,c_yes,header,modelVersion,c_fatal, &
                     "multiple grids not converted yet")
                   ngrid = isched(npass,1)
                   do nfeed=1,isched(npass, 4)
                      call newgrid(ngrid)
                      call node_sendfeed(ngrid)
                      call node_getfeed(nxtnest(ngrid), ngrid)
                      ngrid = nxtnest(ngrid)
                   end do
                end if
             end do

          end if

          ! done advancing all grids one dtlong
          ! update main time variable by a long timestep.

          time = begtime + dtlongn(1)
	  
	  ! average analysis variables over time and
          ! update cfl numbers

          do ngrid=1,ngrids
             call newgrid(ngrid)
             if ((avgtim/=0.) .and. (frqmean/=0. .or. frqboth/=0.))  &
                  call anlavg(mzp,mxp,myp,nzg)
             call cfl(mzp, mxp, myp, nodei0(mynum,ngrid), nodej0(mynum,ngrid))
!print *,'LFR->call CFL',mynum,dtlt
          end do

          ! get max CFL from all processes to probe numerical stability
          ! and to recalculate each grid's deltat if necessary

          if (ideltat<0 .or. ideltat==3) then
             call reduce_max_cfl_and_broadcast(master_num, mynum, &
                  ngrids, cflxy, cflz)
          else
             call reduce_max_cfl_to_master(master_num, mynum, &
                  ngrids, cflxy, cflz)
          end if

          !		Uncoment to calculate execution time and set noInstrumentation = false in ModTimestamp.f90
          !          call SynchronizedTimeStamp(TS_INTEG) ! timestamp CFL

          ! output, if required

          instFlag = (frqanl>0.0 .and. IOUTPUT/=0)
          if (instFlag) then
             instFlag = &
                  mod(time,frqanl)<dtlongn(1)    .or.  &
                  time>=timmax - 0.01*dtlongn(1) .or.  &
                  iflag==1
          end if

          liteFlag = (frqlite>0.0)
          if (liteFlag) then
             liteFlag = &
                  mod(time,frqlite)<dtlongn(1)   .or.  &
                  time>=timmax - 0.01*dtlongn(1) .or.  &
                  iflag==1
          end if

          meanFlag = (frqmean>0.0)
          if (meanFlag) then
             meanFlag = &
                  avgtim>0.0                              .and. &
                  mod(time-avgtim/2., frqmean)<dtlongn(1) .and. &
                  time>=avgtim
             meanFlag = meanFlag .or. &
                  avgtim<0.0     .and. &
                  mod(time,frqmean)<dtlongn(1)
          end if

          histFlag = (frqhis>0.0 .and. IOUTPUT/=0)
          if (histFlag) then
             histFlag = &
                  mod(time,frqhis)<dtlongn(1)    .or.  &
                  time>=timmax - 0.01*dtlongn(1) .or.  &
                  iflag==1
          end if

          posFlag = (frqanl>0.0)
          if (posFlag) then
             posFlag = &
                  mod(time,frqanl)<dtlongn(1)    .or.  &
                  time>=timmax - 0.01*dtlongn(1) .or.  &
                  iflag==1
          end if

          call OutputFields(histFlag, instFlag, liteFlag, meanFlag )

          ! Uncoment to calculate execution time and set noInstrumentation = false in ModTimestamp.f90
          ! call SynchronizedTimeStamp(TS_OUTPUT)

          ! post-process file

!srf----------------------------------------
!IF(applyIAU == 1) posFlag = .FALSE. 
!IF(applyIAU == 2 .and. time > 10800.) then 
!-- this is for IAU =2 and model evaluation, for normal operation comment these lines.
!-- Model starts at 21UTC, after 3 hours the 00UTC forecast is saved.
!-- then  the forecasts files are saved only each 6 hours.

!  posFlag = (mod(time-10800.,10800.)<dtlongn(1)) .and. (mod(time-10800.,21600.)<dtlongn(1)) .or. &
!             time>=timmax - 0.01*dtlongn(1)  
!if(mynum==1)print*,"modoneproc time IAU=",time,posFlag,applyIAU
!ENDIF 
!srf----------------------------------------


          if (posFlag) then
             if(damModule==1) then
                call acumPrecipInDam(nodemxp(mynum,1),nodemyp(mynum,1) &
                ,ia,iz,ja,jz,mcphys_type &
                ,cuparm_g(1)%aconpr &
                ,micro_g(1)%accpr &
                ,micro_g(1)%accpp &
                ,micro_g(1)%accps &
                ,micro_g(1)%accpa &
                ,micro_g(1)%accpg &
                ,micro_g(1)%accph)

                call outputDamPrecip(time,dtlongn(1),timmax,mchnum,master_num)
              endif

             select case (IPOS)

            !!$              call PostProcessHDF5(oneNamelistFile, oneAllPostTypes)

             case (2,3,10,11)
                call PostProcess(AllGrids, oneNamelistFile, oneAllPostTypes)
             end select

          end if

         !Uncoment to calculate execution time and set noInstrumentation = false in ModTimestamp.f90
         ! call SynchronizedTimeStamp(TS_POST)

         !!$ METEOGRAMA
	 if(applyMeteogram)then
		if(mod(time,meteogramFreq) .lt. dtlongn(1) .or.  &
                   time .ge. timmax - 0.01*dtlongn(1)      .or. &
		   iflag .eq. 1)call ProcessLocalMeteogram(oneGrid%meteoPolygons)
	  end if

	 IF(windfarm==1) THEN
		call output_windFarms()
	 END IF

         ! execution time info and dump
          call timing(2,t6)
          w4     = walltime()
          totcpu = totcpu + t6 - t1

     !#ifdef color
    if (mchnum==master_num) then
      if(evaluate/=1) then
        if(maxCflPercent*100>90.0) then 
          write(OUTPUT_UNIT, "(a,i6,a,f9.1,a,f8.3,a,f7.2,a,f7.2,a,f7.2,a)", advance="no") &
            c_lightAqua//" Timestep #"//c_pink, istp,&
            c_noColor//";"//c_lightAqua//" Sim Time"//c_purple, time,&
            c_peach//" [s]"//c_noColor//";"//c_lightAqua//" Wall Time"//c_purple &
            , max(0.001, w4-w3), c_peach//" [s]"//c_noColor//";"//c_lightAqua//" DT"//c_purple &
            ,dtlongn(1),c_peach//" [s]"//c_noColor//";"//c_lightAqua//" sscourn"//c_purple &
            ,sscourn(1),c_noColor//";"//c_lightAqua//" MaxCFL"//c_purple//c_blink &
            ,maxCflPercent*100,c_peach//" [%]"//c_noColor
        else
          write(OUTPUT_UNIT, "(a,i6,a,f9.1,a,f8.3,a,f7.2,a,f7.2,a,f7.2,a)", advance="no") &
            c_lightAqua//" Timestep #"//c_pink, istp,&
            c_noColor//";"//c_lightAqua//" Sim Time"//c_purple, time,&
            c_peach//" [s]"//c_noColor//";"//c_lightAqua//" Wall Time"//c_purple &
            , max(0.001, w4-w3), c_peach//" [s]"//c_noColor//";"//c_lightAqua//" DT"//c_purple &
            ,dtlongn(1),c_peach//" [s]"//c_noColor//";"//c_lightAqua//" sscourn"//c_purple &
            ,sscourn(1),c_noColor//";"//c_lightAqua//" MaxCFL"//c_purple &
            ,maxCflPercent*100,c_peach//" [%]"//c_noColor
        endif
      else
        write(OUTPUT_UNIT, "(a,i6,a,f9.1,a,f11.3,a,f11.3,a,f11.3,a,f11.3,a)", advance="no") &
            c_lightAqua//" Timestep #"//c_pink, istp,&
            c_noColor//";"//c_lightAqua//" Sim Time"//c_purple, time,&
            c_peach//" [s]"//c_noColor//";"//c_lightAqua//" Wall Time"//c_purple &
            , max(0.001, w4-w3), c_peach//" [s]"//c_noColor//";"//c_lightAqua//" DT"//c_purple &
            ,dtlongn(1),c_peach//" [s]"//c_noColor//";"//c_lightAqua//" sscourn"//c_purple &
            ,sscourn(1),c_noColor//";"//c_lightAqua//" rmse_a"//c_purple &
            ,RMSE_average(timeCount),c_noColor
      endif
    endif
!#else
!    if (mchnum==master_num) then
!      write(OUTPUT_UNIT, "(a,i6,a,f9.1,a,f11.3,a,f11.3,a,f11.3,a)", advance="no") &
!         " Timestep #", istp                  &
!        ,"; Sim Time", time                   &
!        ," [s]; Wall Time", max(0.001, w4-w3) &
!        ," [s]; DT ",dtlongn(1)                &
!        ," [s]; sscourn",sscourn(1)
!    endif
!#endif


          ! if selected, gather cpu times and print

          if (prtcputime==0) then
             if (mchnum==master_num) then
                write(OUTPUT_UNIT,"(1x)")
             end if
          else
             call gather_cpu_time_master_print(master_num, mynum, nmachs, t6-t1)
          end if

       end do   !---- do while (time<timmax)


       ! DEBUG-ALF - Barreira artificial
       call parf_barrier(99995)

!XXXXXsrf     end do   !---- do while (IAU_phase =< applyIAU)

    end if  !---- runtype == INITIAL or HISTORY
!================================================================================================
!================================================================================================

!==== simulation ended =======

    wtime_end = walltime()
    if(IPOS == 10 .or. IPOS == 11) iErrNumber = createSitesFile(oneNamelistFile)
    if (mchnum==master_num) then

      write(c0,"(f10.1)") wtime_end - wtime_start
      write(*,fmt='(A)') c_empty
#ifdef color
      iErrNumber=dumpMessage(c_tty,c_yes,'','',c_notice," === Time "//&
           "integration ends; Total run time = "//c_purple//trim(adjustl(c0))//&
           c_peach//" [s]"//c_noColor)
#else
      iErrNumber=dumpMessage(c_tty,c_no,'','',c_notice," === Time "//&
            "integration ends; Total run time = "//trim(adjustl(c0))//" [s]")
#endif
      write(*,fmt='(A)') c_empty
    end if






    ! Deallocating dynamic arrays

    if (runtype(1:7)/='MAKESFC' .and. runtype(1:9)/='MAKEVFILE') &
         call destroyVctr()

    call dealloc_bounds()

    ! Deallocating dxtmax_local
    deallocate(dxtmax_local, STAT=ierr)
    if (ierr/=0) &!call fatal_error("ERROR deallocating dxtmax_local (oneproc)")
      iErrNumber=dumpMessage(c_tty,c_yes,header,modelVersion,c_fatal, &
      "error deallocating dxtmax_local (oneproc)")
    ! Deallocating "mem_grid" data
    call destroyMemGrid()

    ! Deallocating IO Data
    call destroyIoData()

    ! Deallocating Sounding Data
    call destroyRefSounding()

    ! Deallocating Post process info
    if (associated(oneAllPostTypes)) then
       call DestroyPostProcess(oneNamelistFile, oneAllPostTypes)
    end if

    call DestroyGridTree(AllGrids)
    call DestroyNamelistFile(oneNamelistFile)
    call DestroyParallelEnvironment(oneParallelEnvironment)

!--(inspxe)-------------------------------------------------------
!    call rams_mem_dealloc()
!--(inspxe)-------------------------------------------------------

  end subroutine OneProc




  subroutine initOneProc (AllGrids, name_name)
    use dump, only: &
      dumpMessage
    use mem_grid, only: &
       grid_g, &
       oneGlobalGridData
    use mem_aer1, only: dumpAer

    include "constants.f90"
    type(GridTree), pointer :: AllGrids
    character(len=*), intent(in) :: name_name
    character(len=*), parameter :: h="**(initOneProc)**"
    character(len=*), parameter :: header="**(initOneProc)**"

    !---------------------------------------------------------------------
    !     *** This routine is the driver of all initialization packages
    !---------------------------------------------------------------------

    character(len=8) :: rest
    integer :: ifm,icm,ihm,ngr,nv,ierr
    logical :: histFlag, instFlag, meanFlag, liteFlag
    !- for changing initial soil moisture
    LOGICAL, PARAMETER ::  change_soilm=.true.



    ! Set version number for common blocks.

    iversion = 2

    ! Unit numbers for some I/O files.

    iopunt=6


    call opspec3()

    if (runtype(1:7) == 'INITIAL') then

       time=0.
       ngbegun(1:ngrids) = 0

       !----------------------------------------------------------------------
       !                 Initial startup
       !----------------------------------------------------------------------

       ! Read surface, topo, sst, and ndvi files for all grids. Store required
       ! portions of the file at appropriate memory locations. These are the sub-domain
       ! parts.
       ! All the files were checked earlier, so they must be correct.

       do ifm = 1,ngrids

          !--(DMK-LFR NEC-SX6)----------------------------------------------
          call TRSFFieldAndOwnChunk(ifm)
          !--(DMK-LFR NEC-SX6)----------------------------------------------

       enddo
       do ifm = 1,ngrids
          call SfcReadStoreOwnChunk(ifm)
       enddo
       !     Define grid topography, transform, latitude-longitude,
       !        and map factor arrays.

       call newgrid(1)
       call GridSetup(2)

       ! read SST files
       do ifm = 1,ngrids
          call SstReadStoreOwnChunk(1,ifm,ierr)
       enddo

       ! read NDVI files
       do ifm = 1,ngrids
          call NdviReadStoreOwnChunk(1,ifm,ierr)
       enddo

       ! Initialize snowcover arrays
       do ifm = 1,ngrids
          call snowinit(nodemxp(mynum,ifm),nodemyp(mynum,ifm)  &
               ,leaf_g(ifm)%snow_mass(1,1),leaf_g(ifm)%snow_depth(1,1))
       enddo

       ! TEB_SPM
       if (TEB_SPM==1) then
          ! read FUSO (Local Time) files
          do ifm = 1,ngrids
             call FusoReadStoreOwnChunk(ifm)
          enddo
       endif

       ! The following things will be done for INITIAL = 1 or 3...
       if ( initial == 1 .or. initial == 3) then

          !**(JP)** not worked yet
          !call fatal_error(h//"**(JP)** initial==1 or initial==3 was not worked yet")
          iErrNumber=dumpMessage(c_tty,c_yes,header,c_modelVersion,c_fatal, &
          "**(JP)** initial==1 or initial==3 was not worked yet")
          ! If horizontally homogeneous initialization,
          !    subroutine INITHH loops through all grids and initializes
          !    those for which nxtnest = 0.

          if(initial == 1) then
             print*,'Horizontally-homogeneous-INITIAL start of grid- 1'
             call inithh()
          endif

          !If "history" initialization, call INITHIS.
          !This will define initial fields and reference state on grid 1 from
          !history file. Other grids will be interpolated as in a INITIAL=1 start

          if (initial == 3) then
             print*,'History-INITIAL start of grid- 1'
             call inithis()
          endif

          !  On all fine grids, initialize the surface layer characteristics,
          !  the 1-D reference state arrays, the 3-D reference state arrays,
          !  and the prognostic atmospheric fields by interpolation.

          call fmrefs1d(2,ngrids)

          do ifm = 2,ngrids
             icm = nxtnest(ifm)
             if (icm  >=  1) then
                call fmrefs3d(ifm)

                call prgintrp(nnzp(icm),nnxp(icm),nnyp(icm)  &
                     ,nnzp(icm),nnxp(icm),nnyp(icm),0,0,ifm,1,mynum)

                call fmdn0(ifm)

                print*,'Initial interpolation of grid-',ifm
             endif
          enddo


          !--(DMK-CCATT-INI)-----------------------------------------------------
       elseif(initial == 2 .or. initial == 4) then
          !--(DMK-CCATT-FIM)-----------------------------------------------------

          ! If "variable initialization", do it all here

          !--(DMK-CCATT-INI)-----------------------------------------------------
          call VarfReadStoreOwnChunk(AllGrids, 0, initial)
          !--(DMK-CCATT-FIM)-----------------------------------------------------
       endif

       !     Initialize past time level velocity and perturbation Exner function
       !     on all grids.

       do ifm=1,ngrids
          call newgrid(ifm)

          call FieldInit(1)

          call negadj1(mzp,mxp,myp)

          call thermo(mzp,mxp,myp,1,mxp,1,myp,'THRM_ONLY')

          if(ilwrtyp==6 .or. iswrtyp==6 ) THEN
            if (level  ==  3) &
	      call effective_radius(mzp,mxp,myp &
                  ,micro_g(ifm)%rei             &
                  ,micro_g(ifm)%rel             )
          endif
          if (mcphys_type == 0) then
	    if (level  ==  3) then
               call initqin(mzp,mxp,myp        &
                  ,micro_g(ifm)%q2      &
                  ,micro_g(ifm)%q6      &
                  ,micro_g(ifm)%q7      &
                  ,basic_g(ifm)%pi0     &
                  ,basic_g(ifm)%pp      &
                  ,basic_g(ifm)%theta   &
                  ,basic_g(ifm)%dn0     &
                  ,micro_g(ifm)%cccnp   &
                  ,micro_g(ifm)%cifnp   )
            endif

	  elseif(mcphys_type == 1) then

	   if (level  ==  3) then
	     call initqin_2M(mzp,mxp,myp        &
            ,micro_g(ifm)%q2   &
	    ,micro_g(ifm)%q6      &
            ,micro_g(ifm)%q7      &
            ,basic_g(ifm)%pi0     &
            ,basic_g(ifm)%pp      &
            ,basic_g(ifm)%theta   &
            ,basic_g(ifm)%dn0     )


            if(icloud >= 5) call initqin2_2M(mzp,mxp,myp        &
            ,micro_g(ifm)%cccnp   &
            ,micro_g(ifm)%cccmp   &
            ,basic_g(ifm)%dn0   )

            if(idriz  >= 5) call initqin3_2M(mzp,mxp,myp        &
            ,micro_g(ifm)%gccnp   &
            ,micro_g(ifm)%gccmp   &
            ,basic_g(ifm)%dn0   )

            if(ipris  >= 5) call initqin4_2M(mzp,mxp,myp        &
            ,micro_g(ifm)%cifnp   &
            ,basic_g(ifm)%dn0   )

            if(idust > 0 .or. imd1flg > 0 .or. imd2flg > 0)  &
             call initqin5_2M(mzp,mxp,myp    &
            ,micro_g(ifm)%md1np   &
            ,micro_g(ifm)%md2np )
          endif
	 endif
	 !-- initialization of current theta_il field
	 !-- only for RK time integration
	 if(DYNCORE_FLAG==2) then
	    basic_g(ifm)%thc(:,:,:)=basic_g(ifm)%thp(:,:,:)
	 endif
       enddo

       ! If initializing some fields from previous runs...

       ! Do recycle procedure for normal RAMS recycle or only to do
       ! tracers assimilation


       !--(DMK-CCATT-INI)------------------------------------------------------------------
       if ((ipastin == 1) .or. (CCATT==1 .and. RECYCLE_TRACERS==1)) then

!!$        !**(JP)** not worked yet
!!$        call fatal_error(h//"**(JP)** ipastin==1 or (CATT==1 .and. RECYCLE_TRACERS==1) was not worked yet")
          call recycle()

       endif

!call dumpAer('Aer_pos1')

       ! Fill land surface data for all grids that have no standard input files

       ! ALF - For use with SiB

       !DSM if (isfcl <= 2) then
       if (isfcl <= 2  .or. isfcl == 5) then
          call sfcdata
       elseif (isfcl == 3) then

          !**(JP)** not worked yet
          !call fatal_error(h//"**(JP)** isfcl==3 was not worked yet")
          iErrNumber=dumpMessage(c_tty,c_yes,header,c_modelVersion,c_fatal, &
          "**(JP)**isfcl==3 was not worked yet")
       endif

       ! Initialize various LEAF variables.

       if (ipastin == 0) then
          call GeonestNoFile(1,ngrids)
       end if


       !- change initial soil moisture if desired
       if(change_soilm) then

        do ifm=1,ngrids
         call change_soil_moisture_init(nnzp(ifm),nodemxp(mynum,ifm),nodemyp(mynum,ifm)    &
             ,nzg,nzs,npatch,ifm	 &
             ,basic_g(ifm)%theta (1,1,1) &
             ,basic_g(ifm)%pi0   (1,1,1) &
             ,basic_g(ifm)%pp    (1,1,1) &
             ,leaf_g(ifm)%soil_water	   (1,1,1,1)  &
             ,leaf_g(ifm)%soil_energy      (1,1,1,1)  &
             ,leaf_g(ifm)%soil_text	   (1,1,1,1)  &
             ,leaf_g(ifm)%sfcwater_mass    (1,1,1,1)  &
             ,leaf_g(ifm)%sfcwater_energy  (1,1,1,1)  &
             ,leaf_g(ifm)%sfcwater_depth   (1,1,1,1)  &
             ,grid_g(ifm)%glat   (1,1)       &
             ,grid_g(ifm)%glon   (1,1)	     &
             ,grid_g(ifm)%lpw    (1,1)	     &
             ,leaf_g(ifm)%leaf_class(1,1,1)  )
        enddo
       endif



       ! TEB
       if (TEB_SPM==1) then

          !**(JP)** not worked yet
!!$        call fatal_error(h//"**(JP)** TEB_SPM==1 was not worked yet")
          ! Initial values for Common use TEB vars
          ! For now, it is only being used to zero out this four common
          ! use variables.
          ! This variables will be used for other purposes later.
          ! Edmilson D. Freitas 07/07/2006

          do ifm=1,ngrids
             call TEBC_INIT(nodemxp(mynum,ifm), nodemyp(mynum,ifm), npatch, &
                  leaf_g(ifm)%G_URBAN(1,1,1),             &
                  tebc_g(ifm)%EMIS_TOWN(1,1),             &
                  tebc_g(ifm)%ALB_TOWN(1,1),              &
                  tebc_g(ifm)%TS_TOWN(1,1)                )
          enddo

          if (iteb==1) then
             do ifm=1,ngrids
                call TEB_INIT(                      &
                     nnzp(ifm),                     &
                     nodemxp(mynum,ifm),            &
                     nodemyp(mynum,ifm),            &
                     npatch,                        &
                     leaf_g(ifm)%leaf_class(1,1,1), &
                     basic_g(ifm)%theta    (1,1,1), &
                     basic_g(ifm)%rv       (1,1,1), &
                     basic_g(ifm)%pi0      (1,1,1), &
                     basic_g(ifm)%pp       (1,1,1), &
                     teb_g(ifm)%T_ROOF     (1,1,1), &
                     teb_g(ifm)%T_ROAD     (1,1,1), &
                     teb_g(ifm)%T_WALL     (1,1,1), &
                     teb_g(ifm)%TI_BLD       (1,1), &
                     teb_g(ifm)%TI_ROAD      (1,1), &
                     teb_g(ifm)%T_CANYON     (1,1), &
                     teb_g(ifm)%R_CANYON     (1,1), &
                     teb_g(ifm)%TS_ROOF      (1,1), &
                     teb_g(ifm)%TS_ROAD      (1,1), &
                     teb_g(ifm)%TS_WALL      (1,1), &
                     teb_g(ifm)%H_TRAFFIC    (1,1), &
                     teb_g(ifm)%LE_TRAFFIC   (1,1), &
                     teb_g(ifm)%H_INDUSTRY   (1,1), &
                     teb_g(ifm)%LE_INDUSTRY  (1,1), &
                     teb_g(ifm)%WS_ROOF      (1,1), &
                     teb_g(ifm)%WS_ROAD      (1,1), &
                     tebc_g(ifm)%EMIS_TOWN   (1,1), &
                     tebc_g(ifm)%ALB_TOWN    (1,1), &
                     tebc_g(ifm)%TS_TOWN     (1,1), &
                     leaf_g(ifm)%G_URBAN   (1,1,1)  )

             enddo
          endif

          ! Initialize gases and particulate matter

          if (isource==1) then
             do ifm=1,ngrids
                call init_conc1(1, ifm,           &
                     nnzp(ifm),                   &
                     nodemxp(mynum,ifm),          &
                     nodemyp(mynum,ifm),          &
                     npatch,                      &
                     leaf_g(ifm)%G_URBAN (1,1,1), &
                     gaspart_g(ifm)%pno  (1,1,1), &
                     gaspart_g(ifm)%pno2 (1,1,1), &
                     gaspart_g(ifm)%ppm25(1,1,1), &
                     gaspart_g(ifm)%pco  (1,1,1), &
                     gaspart_g(ifm)%pvoc (1,1,1), &
                     gaspart_g(ifm)%pso2 (1,1,1), &
                     gaspart_g(ifm)%pso4 (1,1,1), &
                     gaspart_g(ifm)%paer (1,1,1), &
                     zt                           )

                if (ichemi==1) then  !calling more added scalars for chemistry
                   if (ichemi_in==1) then !reading init.values from previous run
                      call init_conc_prev()
                   else
                      call init_conc2(1, ifm,           &
                           nnzp(ifm),                   &
                           nodemxp(mynum,ifm),          &
                           nodemyp(mynum,ifm),          &
                           npatch,                      &
                           leaf_g(ifm)%G_URBAN(1,1,1),  &
                           gaspart_g(ifm)%po3(1,1,1),   &
                           gaspart_g(ifm)%prhco(1,1,1), &
                           gaspart_g(ifm)%pho2(1,1,1),  &
                           gaspart_g(ifm)%po3p(1,1,1),  &
                           gaspart_g(ifm)%po1d(1,1,1),  &
                           gaspart_g(ifm)%pho(1,1,1),   &
                           gaspart_g(ifm)%proo(1,1,1),  &
                           zt                           )
                   endif
                endif
             enddo
          endif

       endif

       !--(DMK-CCATT-INI)-----------------------------------------------------
       ! CCATT - CHEMISTRY

       !-srf  Initialize the true air density
       if (iexev == 2) then
          do ifm=1,ngrids
             call newgrid(ifm)
             stilt_g(ifm)%dnp(:,:,:)= basic_g(ifm)%dn0(:,:,:)
          enddo
       endif

       if (ccatt == 1 .and. chemistry >= 0 .and. aerosol > 0) then

          !srf- initialize dry-dep constants
          call dep_init(chem_nspecies,dvj) ! (DMK) deposicao seca (cod. limpo)

          !-srf: initialize mixing ratios (only if chem assim is off)
          !-srf: and sources
          do ifm=1,ngrids
             call newgrid(ifm)

             call initial_condition(ifm,nzp,nxp,nyp)
             !if(mynum==1) print*,'->2 chem/aer sources maps reading: ',time/3600.,' grid=',ifm
             call read_sourcemaps(ifm,nnzp(ifm),nodemxp(mynum,ifm),nodemyp(mynum,ifm), &
             	  ia,iz,ja,jz, time,iyear1,imonth1,idate1,itime1,ngrids,timmax,		   &
             	  chem_nspecies,spc_chem_alloc,src,off,nsrc,nvert_src,chem1_src_g, &
             	  bburn,antro, bioge,  geoge,spc_chem_name,on,chemical_mechanism,  &
             	  emiss_ajust,co,aer_nspecies,spc_aer_alloc,spc_aer_name,	   &
             	  src_name,chemistry,ntimes_src,aer1_g,nmodes,aerosol,plumerise,   &
             	  nveg_agreg,plume_mean_g,nzpmax,dzt,grid_g(ifm)%rtgt,grid_g(ifm)%topt, &
             	  transport,plume_g,tropical_forest,boreal_forest,savannah,	    &
             	  grassland,diur_cycle,volcanoes,volc_mean_g,basic_g(ifm)%dn0,zt,zm,&
             	  mchnum, master_num,mass_bin_dist,CO2,ISFCL,aerosol_mechanism     ,&
		  plume_fre_g,emiss_ajust_aer)

             call aer_background(ifm,nnzp(ifm),nodemxp(mynum,ifm),nodemyp(mynum,ifm),&
                  1,nodemxp(mynum,ifm),1,nodemyp(mynum,ifm))
          enddo
       end if
!call dumpAer('Aer_pos2')
       ! Read Radiation Parameters if CARMA or RRTMG Radiation is selected
       if (ilwrtyp==4 .or. iswrtyp==4 .or. ilwrtyp==6 .or. iswrtyp==6 ) then
         CALL master_read_carma_data()
         CALL read_aotMap()
       endif

       ! AKMIN variable:
       do ifm=1,ngrids
          !srf- initialize akmin = f(x,y)
          if (AKMIN(ifm)<0.) then
             call newgrid(ifm)
             call get_akmin2d(ifm, nodemxp(mynum,ifm), nodemyp(mynum,ifm), &
                  akminvar(ifm)%akmin2d, mynum, nodei0, nodej0,akmin(ifm),grid_g(ifm)%topt)
          endif
       enddo
       !--(DMK-CCATT-FIM)-----------------------------------------------------


       !srf-g3d: training for G3d
       do ifm=1,ngrids
          if(nnqparm(ifm) == 3 .or. nnqparm(ifm) == 5) then ! and training==1
             call newgrid(ifm)
             call init_weights(ifm,nodemxp(mynum,ifm),nodemyp(mynum,ifm),nnqparm(ifm))
          endif
       enddo
       !srf-g3d

       if (initial == 3) then

          !**(JP)** not worked yet
          !call fatal_error(h//"**(JP)** sfcinit_hstart was not worked yet")
          iErrNumber=dumpMessage(c_tty,c_yes,header,c_modelVersion,c_fatal, &
          "**(JP)** sfcinit_hstart was not worked yet")
          call sfcinit_hstart()

       end if



    else if(runtype(1:7) == 'HISTORY') then

       !**(JP)** not worked yet
       !call fatal_error(h//"**(JP)** runtype==HISTORY was not worked yet")

       !                  History file start

       call history_start(name_name)

       ! Checking latter if possible to change "grid_setup" by "gridSetup"
       call gridSetup(1)

       ! Check surface,topo,sst,ndvi files. Remake if necessary.
       call MakeSfcFiles()

       ! Read surface and topo files for any added grids
!       do ifm = ngridsh+1,ngrids
!       do ifm = ngridsh+1,ngrids
       do ifm = 1,ngrids
          call SfcReadStoreOwnChunk(ifm)
       enddo
!       do ifm = ngridsh+1,ngrids
       do ifm = 1,ngrids
          !--(DMK-LFR NEC-SX6)----------------------------------------------
          !        call TopReadStoreFullFieldAndOwnChunk(ifm)
          call TRSFFieldAndOwnChunk(ifm)
          !--(DMK-LFR NEC-SX6)----------------------------------------------

       enddo

       ! Checking latter if possible to change "grid_setup" by "gridSetup"
       call gridSetup(2)

       ! Read in sst and ndvi files for all grids
       do ifm = 1,ngrids
          call SstReadStoreOwnChunk(1, ifm, ierr)
       enddo

       do ifm = 1,ngrids
          call NdviReadStoreOwnChunk(1, ifm, ierr)
       enddo

       ! TEB

       if (TEB_SPM==1) then
          ! Read FUSO (Local Time) files for any added grids
          do ifm = ngridsh+1,ngrids
             call FusoReadStoreOwnChunk(ifm)
          enddo
       endif

       do ifm = 1,ngrids
          icm = nxtnest(ifm)
          if (icm  ==  0) then
             call newgrid(ifm)
              call refs3d (mzp,mxp,myp  &
                  ,basic_g(ifm)%pi0  (1,1,1),basic_g(ifm)%dn0  (1,1,1)  &
                  ,basic_g(ifm)%dn0u (1,1,1),basic_g(ifm)%dn0v (1,1,1)  &
                  ,basic_g(ifm)%th0  (1,1,1),grid_g(ifm)%topt  (1,1)    &
                  ,grid_g(ifm)%rtgt  (1,1)  )
          endif
       enddo

       do ifm = 1,min(ngrids,ngridsh)
          icm = nxtnest(ifm)
          if (icm  >  0) call fmrefs3d(ifm)
          call negadj1(mzp,mxp,myp)
       enddo

       ! ALF - For use with SiB

       !DSM if (isfcl <= 2) then
       if (isfcl <= 2  .or. isfcl == 5) then
          call sfcdata
       endif

       ! Heterogenous Soil Moisture Init.

       if ((SOIL_MOIST == 'h').or.(SOIL_MOIST == 'H').or.  &
            (SOIL_MOIST == 'a').or.(SOIL_MOIST == 'A')) then

          do ifm = 1,min(ngrids,ngridsh)
             call newgrid(ifm)
             call soilMoistureInit(nnzp(ifm), nodemxp(mynum,ifm),         &
                  nodemyp(mynum,ifm), nzg, nzs, npatch, ifm,              &
                  basic_g(ifm)%theta, basic_g(ifm)%pi0, basic_g(ifm)%pp,  &
                  leaf_g(ifm)%soil_water, leaf_g(ifm)%soil_energy,        &
                  leaf_g(ifm)%soil_text,                                  &
                  grid_g(ifm)%glat, grid_g(ifm)%glon, grid_g(ifm)%lpw     &
                 ,leaf_g(ifm)%seatp, leaf_g(ifm)%seatf                    )

          enddo

       endif


       !     If any grids are being added for this run, initialize their
       !     surface layer variables, 1-D reference state variables, and
       !     prognostic atmospheric and soil model fields.

       if (ngrids  >  ngridsh) then
          print*,' +-------------------------------------'
          print*,'            !      New grids will be added.       '
          print*,'            !'
          print*,'            ! ',ngridsh,' grid(s) on history file.'
          print*,'            ! ',ngrids, ' grids to be run.        '
          print*,' +-------------------------------------'
          call fmrefs1d(ngridsh+1,ngrids)
          do ifm = ngridsh+1,ngrids
             icm = nxtnest(ifm)
             if (icm  ==  0) then
                !call fatal_error(h//"Attempted to add "//&
                !     "a hemispheric grid on a history restart; "//&
                !     "this cannot be done.")
                iErrNumber=dumpMessage(c_tty,c_yes,header,c_modelVersion,c_fatal, &
                     "Attempted to add "//&
                     "a hemispheric grid on a history restart; "//&
                     "this cannot be done.")
             endif
             call fmrefs3d(ifm)
             call prgintrp(nnzp(icm),nnxp(icm),nnyp(icm)  &
                  ,nnzp(icm),nnxp(icm),nnyp(icm),0,0,ifm,1,mynum)
             print*,'History start interpolation of added grid-',ngrid

             call fmdn0(ifm)
             call newgrid(ifm)

             call FieldInit(0)

             call negadj1(nzp,nxp,nyp)

             call thermo(nzp,nxp,nyp,1,nxp,1,nyp,'THRM_ONLY')

            if (mcphys_type == 0) then
	     if (level  ==  3) then
               call initqin(mzp,mxp,myp        &
                  ,micro_g(ifm)%q2      &
                  ,micro_g(ifm)%q6      &
                  ,micro_g(ifm)%q7      &
                  ,basic_g(ifm)%pi0     &
                  ,basic_g(ifm)%pp      &
                  ,basic_g(ifm)%theta   &
                  ,basic_g(ifm)%dn0     &
                  ,micro_g(ifm)%cccnp   &
                  ,micro_g(ifm)%cifnp   )
             endif

	    elseif(mcphys_type == 1) then

	    if (level  ==  3) then
	     call initqin_2M(mzp,mxp,myp        &
            ,micro_g(ifm)%q2      &
            ,micro_g(ifm)%q6      &
            ,micro_g(ifm)%q7      &
            ,basic_g(ifm)%pi0     &
            ,basic_g(ifm)%pp      &
            ,basic_g(ifm)%theta   &
            ,basic_g(ifm)%dn0     )

            if(icloud >= 5) call initqin2_2M(mzp,mxp,myp        &
            ,micro_g(ifm)%cccnp   &
            ,micro_g(ifm)%cccmp   &
            ,basic_g(ifm)%dn0   )

            if(idriz  >= 5) call initqin3_2M(mzp,mxp,myp        &
            ,micro_g(ifm)%gccnp   &
            ,micro_g(ifm)%gccmp   &
            ,basic_g(ifm)%dn0   )

            if(ipris  >= 5) call initqin4_2M(mzp,mxp,myp        &
            ,micro_g(ifm)%cifnp   &
            ,basic_g(ifm)%dn0   )

            if(idust > 0 .or. imd1flg > 0 .or. imd2flg > 0)  &
             call initqin5_2M(mzp,mxp,myp    &
            ,micro_g(ifm)%md1np   &
            ,micro_g(ifm)%md2np )
          endif
	 endif

             ! Heterogenous Soil Moisture Init.
             if ((SOIL_MOIST == 'h').or.(SOIL_MOIST == 'H').or.  &
                  (SOIL_MOIST == 'a').or.(SOIL_MOIST == 'A')) then
                call soilMoistureInit(nnzp(ifm), nodemxp(mynum,ifm),         &
                     nodemyp(mynum,ifm), nzg, nzs, npatch, ifm,              &
                     basic_g(ifm)%theta, basic_g(ifm)%pi0, basic_g(ifm)%pp,  &
                     leaf_g(ifm)%soil_water, leaf_g(ifm)%soil_energy,        &
                     leaf_g(ifm)%soil_text,                                  &
                     grid_g(ifm)%glat, grid_g(ifm)%glon, grid_g(ifm)%lpw     &
                    ,leaf_g(ifm)%seatp, leaf_g(ifm)%seatf                    )

             endif


          enddo


          !Fill land surface data for all grids that have no standard input files

          call GeonestNofile(ngridsh+1,ngrids)

       elseif (ngrids  <  ngridsh) then
          print*,' +-------------------------------------'
          print*,'            !      Grids will be subtracted.       '
          print*,'            !'
          print*,'            ! ',NGRIDSH,' grid(s) on history file.'
          print*,'            ! ',NGRIDS, ' grids to be run.        '
          print*,' +-------------------------------------'
       endif

       ! Read Radiation Parameters if CARMA or RRTMG Radiation is selected
       if (ilwrtyp==4 .or. iswrtyp==4 .or. ilwrtyp==6 .or. iswrtyp==6 ) then
               CALL master_read_carma_data()
               CALL read_aotMap()
       endif

       !--(DMK-CCATT-INI)-----------------------------------------------------
       if (ccatt == 1 .and. chemistry >= 0) then

          !srf- initialize dry-dep constants
          call dep_init(chem_nspecies,dvj) ! (DMK) deposicao seca (cod. limpo)

          !-srf: initialize mixing ratios (only if chem assim is off)
          !-srf: and sources
          do ifm=1,ngrids
             call newgrid(ifm)
             !if(mynum==1) print*,'->3 chem/aer sources maps reading: ',time/3600.,' grid=',ifm
             call read_sourcemaps(ifm,nnzp(ifm),nodemxp(mynum,ifm),nodemyp(mynum,ifm), &
             	  ia,iz,ja,jz,  						   &
		  time,iyear1,imonth1,idate1,itime1,ngrids,timmax,		   &
             	  chem_nspecies,spc_chem_alloc,src,off,nsrc,nvert_src,chem1_src_g, &
             	  bburn,antro, bioge,  geoge,spc_chem_name,on,chemical_mechanism,  &
             	  emiss_ajust,co,aer_nspecies,spc_aer_alloc,spc_aer_name,	   &
             	  src_name,chemistry,ntimes_src,aer1_g,nmodes,aerosol,plumerise,   &
             	  nveg_agreg,plume_mean_g,nzpmax,dzt,grid_g(ifm)%rtgt,grid_g(ifm)%topt, &
             	  transport,plume_g,tropical_forest,boreal_forest,savannah,	    &
             	  grassland,diur_cycle,volcanoes,volc_mean_g,basic_g(ifm)%dn0,zt,zm,&
             	  mchnum, master_num,mass_bin_dist,CO2,ISFCL,aerosol_mechanism,&
		  plume_fre_g,emiss_ajust_aer)

	   enddo
       end if

    else
       !call fatal_error(h//" wrong runtype")
       iErrNumber=dumpMessage(c_tty,c_yes,header,c_modelVersion,c_fatal, &
            "wrong runtype")
    endif

    ! micro initialization

    call newgrid(1)

    if     (mcphys_type == 0) then
     call micro_master()
    elseif (mcphys_type == 1) then
     call micro_master_2M()
    endif

    !       Fill latitude-longitude, map factor, and Coriolis arrays.
    do ifm = 1,ngrids
       call newgrid(ifm)
       call fcorio(mxp,myp           &
            ,basic_g(ifm)%fcoru (1,1)  &
            ,basic_g(ifm)%fcorv (1,1)  &
            ,grid_g(ifm)%glat   (1,1)  )
    enddo


    !  If we are doing one-way nesting or varfile nudging, inventory,
    !     prepare history/varfile files
    !     and fill past/future nudging arrays for start of simulation

!call dumpAer('Aer_pos3')

    if ( nud_type == 1 ) then

       !**(JP)** not worked yet
       !call fatal_error(h//"**(JP)** nud_type==1 was not worked yet")
       iErrNumber=dumpMessage(c_tty,c_yes,header,c_modelVersion,c_fatal, &
            "**(JP)** nud_type==1 was not worked yet")
       call nud_read(1)

       !--(DMK-CCATT-INI)-----------------------------------------------------
    elseif(nud_type == 2 .or. nud_type == 4) then
       !--(DMK-CCATT-FIM)-----------------------------------------------------

       !--(DMK-CCATT-INI)-----------------------------------------------------
       call VarfReadStoreOwnChunk(AllGrids, 1, nud_type)
       !--(DMK-CCATT-FIM)-----------------------------------------------------

    endif

    ! Do same if doing condensate nudging

    if ( nud_cond == 1 ) then

       !**(JP)** not worked yet
       !call fatal_error(h//"**(JP)** nud_cond==1 was not worked yet")
       iErrNumber=dumpMessage(c_tty,c_yes,header,c_modelVersion,c_fatal, &
            "**(JP)** nud_cond==1 was not worked yet")
       call cond_read(1)

    end if

    ! Process and read observations for ODA - observational data assimilation

    if (if_oda == 1) then

       !**(JP)** not worked yet
       !call fatal_error(h//"**(JP)** if_oda==1 was not worked yet")
       iErrNumber=dumpMessage(c_tty,c_yes,header,c_modelVersion,c_fatal, &
            "**(JP)** if_oda==1 was not worked yet")
       call oda_read()

    end if


    ! Read cumulus heating fields

    if (if_cuinv == 1) then

       !**(JP)** not worked yet
       !call fatal_error(h//"**(JP)** if_cuinv==1 was not worked yet")
       iErrNumber=dumpMessage(c_tty,c_yes,header,c_modelVersion,c_fatal, &
            "**(JP)** if_cuinv==1 was not worked yet")
       call cu_read(1)

    end if

    ! Initialize urban canopy drag coefficients

    if (if_urban_canopy == 1) then

       !**(JP)** not worked yet
       !call fatal_error(h//"**(JP)** if_urban_canopy==1 was not worked yet")
       iErrNumber=dumpMessage(c_tty,c_yes,header,c_modelVersion,c_fatal, &
            "**(JP)** if_urban_canopy==1 was not worked yet")
       call urb_drag_init()

    end if

    ! one process prints locations of all grids

    if (mchnum == master_num) then
       call gridloc_prt()
    end if

    ! Save initial fields on history and analysis files

    iinput=ioutput
    if (runtype  ==  'HISTORY') then
       rest = 'yes'
    else
       rest  ='no'
    end if

    ! produce analysis and history output
    !**(JP)** lite and mean output options not converted

!!$  histFlag=.true.; instFlag=.true.; liteFlag=frqlite > 0.; meanFlag=frqmean >0.
!!$  histFlag=.true.; instFlag=.true.; liteFlag=.false.; meanFlag=.false.
    if (IOUTPUT/=0) then
       histFlag=.true.; instFlag=.true.
    else
       histFlag=.false.; instFlag=.false.
    endif
    if (FRQLITE/=0.) then
       liteFlag=.true.
    else
       liteFlag=.false.
    endif
    if (AVGTIM/=0.) then
       meanFlag=.true.
    else
       meanFlag=.false.
    endif
    !srf
    !--(DMK-CCATT-INI)--------------------------------------------------------
    iF(applyDF)then
       !--(DMK-CCATT-FIM)--------------------------------------------------------
       histFlag=.false.; instFlag=.false.
       liteFlag=.false.
       meanFlag=.false.
    endif
    !--(DMK-CCATT-INI)--------------------------------------------------------
    if( applyDF .and. time <0.00001) instFlag=.true.
    !--(DMK-CCATT-FIM)--------------------------------------------------------
    !srf

    call OutputFields(histFlag, instFlag, liteFlag, meanFlag)

    ! Save initial fields into the averaged arrays

    if(avgtim /= 0.)then

       !**(JP)** not converted yet
       !call fatal_error(h//" avgtim /= 0 was not converted yet")
       iErrNumber=dumpMessage(c_tty,c_yes,header,c_modelVersion,c_fatal, &
            "**(JP)** avgtim/=0 was not worked yet")

       do ngr=1,ngrids
          do nv=1,num_var(ngr)
             if(vtab_r(nv,ngr)%imean == 1) then
                call atob_long(vtab_r(nv,ngr)%npts, vtab_r(nv,ngr)%var_p, &
                     vtab_r(nv,ngr)%var_m)
             endif
          enddo
       enddo
    endif

    ! Print the header information and initial fields
    if (mchnum == master_num) then
       call PrtOpt()
    end if

    !**(JP)** ver o que fazer com esses prints
!call dumpAer('Aer_pos4')
!!$  ngrid=1
!!$  call prtopt(6)
!!$  if (initfld  ==  1) then
!!$     do ifm = 1,ngrids
!!$        call newgrid(ifm)
!!$        call prtout()
!!$     enddo
!!$  endif

!!$  call opspec3()
  end subroutine initOneProc







  subroutine comm_time(isendflg, isendlite, isendmean, isendboth, &
       isendbackflg, isendiv, isendsst, isendndvi, isendsrc)

!!$    use mem_varinit
!!$    use mem_cuparm
!!$    use io_params
!!$    use mem_grid


    integer, intent(out) :: isendflg, isendlite, isendmean, isendboth, &
         isendbackflg, isendiv, isendsst, isendndvi

    !--(DMK-CCATT-INI)-----------------------------------------------------
    integer, intent(out) :: isendsrc
    !--(DMK-CCATT-FIM)-----------------------------------------------------

    real :: timemf
    integer :: ifm
    real :: frqqueim  ! CATT

    !         ISENDFLG designates whether nodes should send back
    !            stuff things it normally doesn't have to
    !            at the end of timestep for history/analysis write,
    !            load balancing, etc.

    !         isendflg  = the usual RAMS stuff
    !         isendlite = the "lite" variables
    !         isendmean = the "mean" variasbles
    !         isendboth = Both the "mean" and "lite" variables

    !            Determines whether nodes send stuff back at the END of the
    !            timestep!!!!!

    timemf = time + dtlongn(1)

    isendflg     = 0
    isendlite    = 0
    isendmean    = 0
    isendboth    = 0
    isendbackflg = 0 ! ALF
    isendiv      = 0 ! ALF
    isendsst     = 0 ! ALF
    isendndvi    = 0 ! ALF

    !--(DMK-CCATT-INI)-----------------------------------------------------
    isendsrc     = 0
    !--(DMK-CCATT-FIM)-----------------------------------------------------


    !--(DMK-CCATT-INI)---------------------------------------------------------
    if (ccatt == 1) then
       if (chemistry >= 0 &
            .and. srcmapfn(1:len_trim(srcmapfn)) /= 'NONE' &
            .and. srcmapfn(1:len_trim(srcmapfn)) /= 'none') then

          !time para leitura dos mapas de queimadas
          !frqqueim=24.*3600.
          !inicializando as fontes as 00UTC
          !print *,'LFR-DBG: ',timemf,itime1,srctime2,dtlongn(1)
          !if ( mod(timemf + 0.01*itime1*3600.,srctime2) .lt. dtlongn(1) .and. & !LFR !!!!!!!!!!
          if ( mod(timemf + 0.01*itime1*3600.,86400.0) .lt. dtlongn(1) .and. &
               timemf.ge.0.*3600.) then
             !       if (timemf>=srctime2 .and. timemf<timmax) then
             isendbackflg = 1
             write(*,*) "**** Setting isendsrc to 0. The source map will be the last one! (must be reorganized!!!) ****"
             isendsrc = 0!1
          endif

       endif
    endif
    !--(DMK-CCATT-OLD)---------------------------------------------------------
    !  if (CATT == 1) then
    !     !----srf-
    !     !time para leitura dos mapas de queimadas
    !     frqqueim=24.*3600.
    !     !inicializando as fontes as 00UTC
    !     if ( mod(timemf-0.*3600.,frqqueim) .lt. dtlongn(1) .and. &
    !          timemf.ge.0.*3600.) then
    !        isendflg = 1
    !        isendbackflg = 1 ! ALF
    !        !        print*,'QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ'
    !        !        print*,'BURN MAP READING: ',isendflg,timemf/3600.
    !        !        print*,'QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ'
    !     endif
    !     !----srf-
    !  endif
    !--(DMK-CCATT-FIM)---------------------------------------------------------

    if(frqlite > 0.) then
       if (mod(timemf,frqlite) < dtlongn(1)) isendlite=1
    endif

    if (frqmean > 0.) then
       if(avgtim > 0.)then
          if(mod(timemf-avgtim/2.,frqmean) < dtlongn(1) .and.  &
               timemf >= avgtim) isendmean=1
       elseif(avgtim < 0.)then
          if(mod(timemf,frqmean) < dtlongn(1)) isendmean=1
       endif
    endif
    if(runtype(1:7) == 'INITIAL'.and.timemf < dtlongn(1)) isendmean=0
    if(runtype(1:7) == 'HISTORY'.and.timemf <= timstr) isendmean=0

    if (frqboth > 0.) then
       if(avgtim > 0.)then
          if(mod(timemf-avgtim/2.,frqboth) < dtlongn(1) .and.  &
               timemf >= avgtim) isendboth=1
       elseif(avgtim < 0.)then
          if(mod(timemf,frqboth) < dtlongn(1)) isendboth=1
       endif
    endif
    if(runtype(1:7) == 'INITIAL'.and.timemf < dtlongn(1))isendboth=0
    if(runtype(1:7) == 'HISTORY'.and.timemf <= timstr)isendboth=0

    if (ioutput  /=  0) then
       if ( mod(timemf,frqanl)  <  dtlongn(1) .or.  &
            mod(timemf,frqhis)  <  dtlongn(1)) then
          isendflg = 1
          !return
       endif
    endif

    if( timemf  >=  timmax - .01*dtlongn(1) ) then
       isendflg = 1
       !return
    endif

    if  ( (nud_type == 2 .or. nud_type == 4 ) .and. timemf  >=  vtime2  &
         .and. timemf  <  timmax) then
!!$     isendflg = 1 ! Not used if sending just the necessary input data
       isendbackflg = 1 ! ALF
       isendiv = 1 ! ALF
       return
    endif

    if  (nud_type == 1 .and. timemf  >=  htime2  &
         .and. timemf  <  timmax) then
       isendflg = 1
       isendbackflg = 1 ! ALF
       !return
    endif

    if  ( nud_cond == 1 .and. timemf  >=  condtime2  &
         .and. timemf  <  timmax) then
       isendflg = 1
       isendbackflg = 1 ! ALF
       !return
    endif

    if (mod(timemf,frqprt)  <  dtlongn(1) ) then
       isendflg = 1
       !return
    endif

    if (iupdsst  ==  1 ) then
       do ifm = 1,ngrids
          if (isstflg(ifm)  ==  1) then
             if (timemf  >=  ssttime2(ifm) .and.  &
                  timemf  <  timmax) then
                !isendflg = 1
                isendbackflg = 1 ! ALF
                isendsst     = 1 ! ALF
                !return
             endif
          endif
       enddo
    endif

    if (iupdndvi  ==  1 ) then
       do ifm = 1,ngrids
          if (ndviflg(ifm)==1) then
             if (timemf>=ndvitime2(ifm) .and. timemf<timmax) then
                !isendflg = 1
                isendbackflg = 1 ! ALF
                isendndvi    = 1 ! ALF
                !return
             endif
          endif
       enddo
    endif

    if (if_cuinv  ==  1 ) then
       do ifm = 1,ngrids
          if (timemf  >=  cu_times(ncufl+1) .and.  &
               timemf  <  timmax) then
             isendflg = 1
             isendbackflg = 1 ! ALF
             !return
          endif
       enddo
    endif

    return
  end subroutine comm_time




end module ModOneProc

!subroutine lfr_debug(header,version,message,value)
!  implicit NONE
!  character(len=*), intent(in) :: header,version, message
!  integer, intent(in) :: value
!
!  print *,'LFR-DBG '//header//' % '//version//' % '//message//' : ',value
!  call flush(6)
!
!end subroutine lfr_debug
