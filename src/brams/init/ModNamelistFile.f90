!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

module ModNamelistFile

  use grid_dims
  use dump

  private
  public :: namelistFile
  public :: CreateNamelistFile
  public :: DestroyNamelistFile
  public :: GetNamelistFileName
  public :: ReadNamelistFile
  public :: BroadcastNamelistFile
  public :: DumpNamelistFile
  public :: TimeUnitsToSeconds


  include "files.h"
  integer, parameter :: maxodagrids=10
  integer, parameter :: maxcugrids=10
  integer, parameter :: ncat_dummy=8 ! for micro_2m
  integer, parameter :: maxsndg=200
  integer, parameter :: maxisn=100
  integer, parameter :: maxagrds=10

!--(DMK-CCATT-INI)-----------------------------------------------------------
  integer, parameter :: nsrc=4   !  number_sources
!--(DMK-CCATT-FIM)-----------------------------------------------------------

  type namelistFile
     character(len=f_name_length) :: fileName

     ! namelist /MODEL_GRID/

     character(len=64) :: expnme
     character(len=16) :: runtype
     character(len=1)  :: timeunit
     real    :: timmax
     integer :: load_bal
     integer :: imonth1
     integer :: idate1
     integer :: iyear1
     integer :: itime1
     integer :: ngrids
     integer :: nnxp(maxgrds)
     integer :: nnyp(maxgrds)
     integer :: nnzp(maxgrds)
     integer :: nzg
     integer :: nzs
     integer :: nxtnest(maxgrds)
     character(len=f_name_length) :: domain_fname
     integer :: if_adap
     integer :: ihtran
     real :: deltax
     real :: deltay
     real :: deltaz
     real :: dzrat
     real :: dzmax
     integer :: fixLevels
     real :: zz(nzpmax)
     real :: dtlong
     integer :: nacoust
     integer :: ideltat
     integer :: nstratx(maxgrds)
     integer :: nstraty(maxgrds)
     integer :: nndtrat(maxgrds)
     integer :: nestz1
     integer :: nstratz1(nzpmax)
     integer :: nestz2
     integer :: nstratz2(nzpmax)
     real :: polelat
     real :: polelon
     real :: centlat(maxgrds)
     real :: centlon(maxgrds)
     integer :: ninest(maxgrds)
     integer :: njnest(maxgrds)
     integer :: nknest(maxgrds)
     integer :: nnsttop(maxgrds)
     integer :: nnstbot(maxgrds)
     real :: gridu(maxgrds)
     real :: gridv(maxgrds)

     ! namelist /CCATT_INFO/

     integer :: ccatt
     integer :: chemistry
     character(len=20) :: split_method
     real    :: chem_timestep
     integer :: chemistry_aq
     !integer :: aerosol
     integer :: chem_assim
     character(len=256) :: srcmapfn
     integer :: recycle_tracers
     character(len=32) :: def_proc_src
     integer :: diur_cycle(nsrc)
     integer :: na_extra2d
     integer :: na_extra3d
     integer :: plumerise
     real :: prfrq
     integer :: volcanoes
     !Matrix
     REAL :: aer_timestep
     INTEGER :: aerosol
     integer :: aer_assim
     integer :: mech


!--(DMK-CCATT-OLD)-----------------------------------------------------------
!     ! namelist /CATT_INFO/
!
     integer :: catt
     character(len=256) :: firemapfn
!    integer :: recycle_tracers
!    integer :: plumerise
     integer :: define_proc
!    real :: prfrq
!--(DMK-CCATT-FIM)-----------------------------------------------------------

     ! namelist /TEB_SPM_INFO/                                              &

     integer :: teb_spm
     character(len=f_name_length) :: fusfiles
     integer :: ifusflg(maxgrds)
     character(len=f_name_length) :: ifusfn(maxgrds)
     integer             :: ichemi
     integer             :: ichemi_in
     character (len=f_name_length) :: chemdata_in
     integer             :: isource
     character(len=3)    :: weekdayin
     real            :: rushh1
     real            :: rushh2
     real            :: daylight
     real            :: efsat
     real            :: efsun
     real            :: eindno
     real            :: eindno2
     real            :: eindpm
     real            :: eindco
     real            :: eindso2
     real            :: eindvoc
     real            :: eveino
     real            :: eveino2
     real            :: eveipm
     real            :: eveico
     real            :: eveiso2
     real            :: eveivoc
     integer         :: iteb
     real            :: tminbld
     integer         :: nteb
     real            :: hc_roof(maxsteb)
     real            :: tc_roof(maxsteb)
     real            :: d_roof(maxsteb)
     real            :: hc_road(maxsteb)
     real            :: tc_road(maxsteb)
     real            :: d_road(maxsteb)
     real            :: hc_wall(maxsteb)
     real            :: tc_wall(maxsteb)
     real            :: d_wall(maxsteb)
     integer         :: nurbtype
     integer         :: ileafcod(maxubtp)
     real            :: z0_town(maxubtp)
     real            :: bld(maxubtp)
     real            :: bld_height(maxubtp)
     real            :: bld_hl_ratio(maxubtp)
     real            :: aroof(maxubtp)
     real            :: eroof(maxubtp)
     real            :: aroad(maxubtp)
     real            :: eroad(maxubtp)
     real            :: awall(maxubtp)
     real            :: ewall(maxubtp)
     real            :: htraf(maxubtp)
     real            :: hindu(maxubtp)
     real            :: pletraf(maxubtp)
     real            :: pleindu(maxubtp)


     !namelist /MODEL_FILE_INFO/                                           &

     integer :: initial
     integer :: nud_type
     character(len=f_name_length)       :: varfpfx
     real :: vwait1
     real :: vwaittot
     character(len=f_name_length) :: nud_hfile
     integer :: nudlat
     real :: timeWindowIAU
     real :: ramp
     real :: tnudlat
     real :: tnudcent
     real :: tnudtop
     real :: znudtop
     real :: wt_nudge_grid(maxgrds)
     real :: wt_nudge_uv
     real :: wt_nudge_th
     real :: wt_nudge_pi
     real :: wt_nudge_rt
     integer :: applyIAU
     character(len=f_name_length) :: fileNameIAU
     integer :: nud_cond
     character(len=f_name_length) :: cond_hfile
     real :: tcond_beg
     real :: tcond_end
     real :: t_nudge_rc
     real :: wt_nudgec_grid(maxgrds)
     integer            :: if_oda
     character(len=128) :: oda_upaprefix
     character(len=128) :: oda_sfcprefix
     real               :: frqoda
     real               :: todabeg
     real               :: todaend
     real               :: tnudoda
     real               :: wt_oda_grid(maxodagrids)
     real               :: wt_oda_uv
     real               :: wt_oda_th
     real               :: wt_oda_pi
     real               :: wt_oda_rt
     real               :: roda_sfce(maxodagrids)
     real               :: roda_sfc0(maxodagrids)
     real               :: roda_upae(maxodagrids)
     real               :: roda_upa0(maxodagrids)
     real               :: roda_hgt(maxodagrids)
     real               :: roda_zfact(maxodagrids)
     real               :: oda_sfc_til
     real               :: oda_sfc_tel
     real               :: oda_upa_til
     real               :: oda_upa_tel
     integer :: if_cuinv
     character(len=128) :: cu_prefix
     real :: tnudcu
     real :: wt_cu_grid(maxcugrids)
     real :: tcu_beg
     real :: tcu_end
     real :: cu_tel
     real :: cu_til
     real :: timstr
     character(len=f_name_length) :: hfilin
     integer :: ipastin
     character(len=f_name_length) :: pastfn
     integer :: ioutput
     character(len=f_name_length) :: hfilout
     character(len=f_name_length) :: afilout
     integer :: iclobber
     integer :: ihistdel
     real :: frqhis
     real :: frqanl
     real :: frqlite
     integer :: ipos
     character(len=20)  :: xlite
     character(len=20)  :: ylite
     character(len=20)  :: zlite
     integer :: nlite_vars
     character(len=32)  :: lite_vars(maxlite)
     real :: avgtim
     real :: frqmean
     real :: frqboth
     integer :: kwrite
     real             :: frqprt
     integer            :: initfld
     integer                :: prtcputime
     character(len=f_name_length) :: topfiles
     character(len=f_name_length) :: sfcfiles
     character(len=f_name_length) :: sstfpfx
     character(len=f_name_length) :: ndvifpfx
     integer :: itoptflg(maxgrds)
     integer :: isstflg(maxgrds)
     integer :: ivegtflg(maxgrds)
     integer :: isoilflg(maxgrds)
     integer :: ndviflg(maxgrds)
     integer :: nofilflg(maxgrds)
     integer            :: iupdndvi
     integer            :: iupdsst
     character(len=f_name_length) :: itoptfn(maxgrds)
     character(len=f_name_length) :: isstfn(maxgrds)
     character(len=f_name_length) :: ivegtfn(maxgrds)
     character(len=f_name_length) :: isoilfn(maxgrds)
     character(len=f_name_length) :: ndvifn(maxgrds)
     integer :: itopsflg(maxgrds)
     real  :: toptenh(maxgrds)
     real  :: toptwvl(maxgrds)
     integer :: iz0flg(maxgrds)
     real  :: z0max(maxgrds)
     real  :: z0fact
     integer            :: mkcoltab
     character(len=f_name_length) :: coltabfn
     character(len=f_name_length) :: mapaotfile
     character(len=f_name_length) :: julesin

     ! namelist /MODEL_OPTIONS/ &
     integer :: advmnt
     integer :: ghostzonelength
     integer :: naddsc
     integer :: icorflg
     integer :: dyncore_flag
     integer :: pd_or_mnt_constraint
     integer :: order_h
     integer :: order_v

!--(DMK-CCATT-INI)-----------------------------------------------------------
     integer :: iexev
     integer :: imassflx
     integer :: vveldamp
!--(DMK-CCATT-FIM)-----------------------------------------------------------
     integer :: ibnd
     integer :: jbnd
     real :: cphas
     integer :: lsflg
     integer :: nfpt
     real :: distim
     integer :: sfire
     integer :: iswrtyp
     integer :: ilwrtyp
     character(LEN=f_name_length) :: raddatfn
     real    :: radfrq
     real    :: radtun
     integer :: lonrad
     integer :: nnqparm(maxgrds)
     character(len=2) :: closure_type
     integer :: nnshcu(maxgrds)
     real :: confrq
     real    :: shcufrq
     real :: wcldbs

     integer :: g3d_spread
     integer :: g3d_smoothh
     integer :: g3d_smoothv

     integer :: npatch
     integer :: nvegpat
     integer :: isfcl
     integer :: isfcl_ocean

     integer :: nvgcon
     real    :: pctlcon
     integer :: nslcon
     real    :: drtcon
     real    :: zrough
     real    :: albedo
     real    :: seatmp
     real    :: dthcon
     character(len=1)   :: soil_moist
     character(len=1)   :: soil_moist_fail
     character(len=f_name_length) :: usdata_in
     character(len=f_name_length) :: usmodel_in
     real    :: slz(nzgmax)
     real    :: slmstr(nzgmax)
     real    :: stgoff(nzgmax)
     integer :: if_urban_canopy
     integer :: idiffk(maxgrds)
     integer :: ihorgrad
     real    :: csx(maxgrds)
     real    :: csz(maxgrds)
     real    :: xkhkm(maxgrds)
     real    :: zkhkm(maxgrds)
     real    :: akmin(maxgrds)
     integer            :: mcphys_type
     integer            :: idriz
     integer            :: iccnlev
     integer            :: irime
     integer            :: iplaws

     integer            :: level
     integer            :: icloud
     integer            :: irain
     integer            :: ipris
     integer            :: isnow
     integer            :: iaggr
     integer            :: igraup
     integer            :: ihail
     real               :: cparm
     real               :: rparm
     real               :: pparm
     real               :: sparm
     real               :: aparm
     real               :: gparm
     real               :: hparm
     real               :: dparm
     real               :: cnparm
     real               :: gnparm
     real               :: epsil

     real               :: gnu(ncat_dummy)
     integer 		:: windfarm
     character(len=f_name_length) :: wfFile
     integer :: damModule
     real :: frqPrecip
     character(len=256) :: damOutPrefix
     integer :: evaluate
     character(len=256) :: evaluatePrefix


     !namelist /MODEL_SOUND/ &

     integer             :: ipsflg
     integer             :: itsflg
     integer             :: irtsflg
     integer             :: iusflg
     real                :: hs(maxsndg)
     real                :: ps(maxsndg)
     real                :: ts(maxsndg)
     real                :: rts(maxsndg)
     real                :: us(maxsndg)
     real                :: vs(maxsndg)

     !namelist /MODEL_PRINT/ &

     integer            :: nplt
     character(len=16)  :: iplfld(50)
     integer            :: ixsctn(50)
     integer            :: isbval(50)

     !namelist /ISAN_CONTROL/ &

     integer :: iszstage
     integer :: ivrstage
     integer :: isan_inc
     character(len=8)   :: guess1st
     integer :: i1st_flg
     integer :: iupa_flg
     integer :: isfc_flg
     character(len=256) :: iapr
     character(len=256) :: iarawi
     character(len=256) :: iasrfce
     character(len=256) :: varpfx
     integer :: ioflgisz
     integer :: ioflgvar

     !namelist /ISAN_ISENTROPIC/ &

     integer                    :: nisn
     integer                    :: levth(maxisn)
     integer                    :: nigrids
     real                       :: topsigz
     real                       :: hybbot
     real                       :: hybtop
     real                       :: sfcinf
     real                       :: sigzwt
     integer                    :: nfeedvar
     integer                    :: maxsta
     integer                    :: maxsfc
     integer                    :: notsta
     character(len=8)           :: notid(50)
     integer                    :: iobswin
     real                       :: stasep
     integer                    :: igridfl
     real                       :: gridwt(maxagrds)
     real                       :: gobsep
     real                       :: gobrad
     real                       :: wvlnth(maxagrds)
     real                       :: respon(maxagrds)
     real                       :: swvlnth(maxagrds)
     !
     integer                    :: icFileType
     character(len=256)         :: icPrefix 
     character(len=32)          :: wind_u_varname  
     character(len=32)          :: wind_v_varname 
     character(len=32)          :: temperature_varname
     character(len=32)          :: geo_varname
     character(len=32)          :: ur_varname    
     real                       :: initial_latitude
     real                       :: final_latitude
     real                       :: initial_longitude
     real                       :: final_longitude
     integer                    :: z_max_level
     real                       :: dlimit(5)
     real                       :: ulimit(5)
     real                       :: scale_factor(5)
     character(len=256)         :: icGradsPrefix
     integer                    :: ccGradsWrite


     !namelist POST

     integer                       :: nvp
     character(len=f_name_length)  :: vp(200)
     character(len=f_name_length)  :: gprefix
     character(len=f_name_length)  :: csvFile
     character(len=10)             :: anl2gra
     character(len=10)             :: proj
     character(len=3)              :: mean_type
     real                          :: lati(maxgrds)
     real                          :: loni(maxgrds)
     real                          :: latf(maxgrds)
     real                          :: lonf(maxgrds)
     integer                       :: zlevmax(maxgrds)
     integer                       :: ipresslev
     integer                       :: inplevs
     integer                       :: iplevs(nzpmax)
     character(len=10)             :: mechanism
     character(len=20)             :: ascii_data
     real                          :: site_lat
     real                          :: site_lon

     !namelist digital filter
     logical :: applyDigitalFilter
     real    :: digitalFilterTimeWindow

     !namelist meteogram
     logical                      :: applyMeteogram
     real                         :: meteogramFreq
     character(len=f_name_length) :: meteogramMap
     character(len=f_name_length) :: meteogramDir
  end type namelistFile

contains





  subroutine CreateNamelistFile(oneNamelistFile)
    implicit none
    type(namelistFile), pointer :: oneNamelistFile
    if (associated(oneNamelistFile)) then
       deallocate(oneNamelistFile)
    end if
    allocate(oneNamelistFile)
  end subroutine CreateNamelistFile




  subroutine DestroyNamelistFile(oneNamelistFile)
    implicit none
    type(namelistFile), pointer :: oneNamelistFile
    if (associated(oneNamelistFile)) then
       deallocate(oneNamelistFile)
    end if
    nullify(oneNamelistFile)
  end subroutine DestroyNamelistFile



  subroutine GetNamelistFileName(oneNamelistFile)
    implicit none
    type(namelistFile), pointer :: oneNamelistFile

    integer :: nargs
    integer :: iarg
    integer :: lenArg
    integer :: status
    logical :: flagName ! true iff arg="-f"; next arg is file name
    character(len=f_name_length) :: arg

    oneNamelistFile%fileName="RAMSIN" ! default namelist

    ! search command line for "-f " <namelist file name>
    ! return default if not found
    nargs = command_argument_count()
    if (nargs >= 0) then
       flagName = .false.
       do iarg = 0, nargs
          call get_command_argument(iarg, arg, lenArg, status)
          !print *, iarg,arg(iarg); call flush(6)
          if (status == 0) then
             if (flagName) then
                oneNamelistFile%fileName = arg(1:lenArg)
                exit
             else
                flagName = arg(1:lenArg) == "-f"
             end if
          end if
       end do
    end if
  end subroutine GetNamelistFileName


  ! ReadNamelistFile:
  !    open, reads and close namelist file
  !    implements defaults for namelist variables
  !    check input options consistency


  subroutine ReadNamelistFile(oneNamelistFile)

    implicit none
    type(namelistFile), pointer :: oneNamelistFile

    include "files.h"



    integer :: i                        ! loop count
    integer :: iunit,iunit2             ! io unit number
    integer, parameter :: firstUnit=20  ! lowest io unit number available
    integer, parameter :: lastUnit=99   ! highest io unit number available
    logical :: op                       ! io unit number opened or not
    logical :: ex                       ! namelist file exists?
    integer :: err                      ! return code on iostat
    character(len=10) :: c0             ! scratch
    character(len=*), parameter :: h="**(ReadNamelistFile)**"  ! program unit name
    character(len=*), parameter :: header=h
    character(len=*), parameter :: version="5.4"

    ! namelist /MODEL_ADV_RAMSIN/
    character(len=f_name_length) :: advanced_ramsin

     namelist /MODEL_ADV_RAMSIN/                                               &
         advanced_ramsin


    ! namelist /MODEL_GRID/

    character(len=64) :: expnme
    character(len=16) :: runtype
    character(len=1)  :: timeunit
    real    :: timmax
    integer :: load_bal
    integer :: imonth1
    integer :: idate1
    integer :: iyear1
    integer :: itime1
    integer :: ngrids
    integer :: nnxp(maxgrds)
    integer :: nnyp(maxgrds)
    integer :: nnzp(maxgrds)
    integer :: nzg
    integer :: nzs
    integer :: nxtnest(maxgrds)
    character(len=f_name_length) :: domain_fname
    integer :: if_adap
    integer :: ihtran
    real :: deltax
    real :: deltay
    real :: deltaz
    real :: dzrat
    real :: dzmax
    integer :: fixLevels
    real :: zz(nzpmax)
    real :: dtlong
    integer :: nacoust
    integer :: ideltat
    integer :: nstratx(maxgrds)
    integer :: nstraty(maxgrds)
    integer :: nndtrat(maxgrds)
    integer :: nestz1
    integer :: nstratz1(nzpmax)
    integer :: nestz2
    integer :: nstratz2(nzpmax)
    real :: polelat
    real :: polelon
    real :: centlat(maxgrds)
    real :: centlon(maxgrds)
    integer :: ninest(maxgrds)
    integer :: njnest(maxgrds)
    integer :: nknest(maxgrds)
    integer :: nnsttop(maxgrds)
    integer :: nnstbot(maxgrds)
    real :: gridu(maxgrds)
    real :: gridv(maxgrds)

    namelist /MODEL_GRIDS/                                               &
         expnme, runtype, timeunit, timmax, imonth1, idate1,   &
         iyear1, itime1, nnxp, nnyp, nnzp, nzg, nzs, &
         deltax, deltay, deltaz, dzrat, dzmax, fixLevels, &
         zz,dtlong, polelat, polelon, centlat, centlon

    namelist /MODEL_GRIDS2/                                               &
         ngrids, load_bal, nxtnest,domain_fname,if_adap, ihtran, &
         nacoust, ideltat, nstratx, nstraty, nndtrat, nestz1,    &
         nstratz1, nestz2, nstratz2, &
         ninest, njnest, nknest, nnsttop, nnstbot, gridu, gridv
!--(DMK-CCATT-FIM)-----------------------------------------------------------


!    ! namelist /CATT_INFO/
!
!    integer :: catt
!    character(len=256) :: firemapfn
!    integer :: recycle_tracers
!    integer :: plumerise
!    integer :: define_proc
!    real :: prfrq

!--(DMK-CCATT-INI)-----------------------------------------------------------
    ! namelist /CCATT_INFO/

    integer :: ccatt
    integer :: chemistry
    character(len=20) :: split_method
    real    :: chem_timestep
    integer :: chemistry_aq
    !integer :: aerosol
    integer :: chem_assim
    character(len=256) :: srcmapfn
    integer :: recycle_tracers
    character(len=32) :: def_proc_src
    integer :: diur_cycle(nsrc)
    integer :: na_extra2d
    integer :: na_extra3d
    integer :: plumerise
    real :: prfrq
    integer :: volcanoes
    !Matrix
    INTEGER :: aerosol
    REAL :: aer_timestep
    integer :: aer_assim
    integer :: mech
    integer :: windfarm
    character(len=f_name_length) :: wfFile
    integer :: damModule
    real :: frqPrecip
    character(len=256) :: damOutPrefix
     integer :: evaluate
     character(len=256) :: evaluatePrefix

!--(DMK-CCATT-FIM)-----------------------------------------------------------

!--(DMK-CCATT-INI)-----------------------------------------------------------
    namelist /CCATT_INFO/                                                &
         ccatt, chemistry, chem_timestep, chem_assim, srcmapfn, plumerise, &
         aerosol,aer_timestep, aer_assim

    namelist /CCATT_INFO2/                                                &
         split_method, chemistry_aq,  &
         chem_assim, recycle_tracers, def_proc_src, diur_cycle,&
         na_extra2d, na_extra3d, plumerise, prfrq, volcanoes, mech
!--(DMK-CCATT-OLD)-----------------------------------------------------------
!    namelist /CATT_INFO/                                                 &
!         catt,                                                           &
!         firemapfn, recycle_tracers,                                     &
!         plumerise, define_proc, prfrq
!--(DMK-CCATT-FIM)-----------------------------------------------------------

    ! namelist /TEB_SPM_INFO/

    integer :: teb_spm
    character(len=f_name_length) :: fusfiles
    integer :: ifusflg(maxgrds)
    character(len=f_name_length) :: ifusfn(maxgrds)
    integer             :: ichemi
    integer             :: ichemi_in
    character (len=f_name_length) :: chemdata_in
    integer             :: isource
    character(len=3)    :: weekdayin
    real            :: rushh1
    real            :: rushh2
    real            :: daylight
    real                :: efsat
    real                :: efsun
    real                :: eindno
    real                :: eindno2
    real                :: eindpm
    real                :: eindco
    real                :: eindso2
    real                :: eindvoc
    real                :: eveino
    real                :: eveino2
    real                :: eveipm
    real                :: eveico
    real                :: eveiso2
    real                :: eveivoc
    integer         :: iteb
    real            :: tminbld
    integer         :: nteb
    real            :: hc_roof(maxsteb)
    real            :: tc_roof(maxsteb)
    real            :: d_roof(maxsteb)
    real            :: hc_road(maxsteb)
    real            :: tc_road(maxsteb)
    real            :: d_road(maxsteb)
    real            :: hc_wall(maxsteb)
    real            :: tc_wall(maxsteb)
    real            :: d_wall(maxsteb)
    integer         :: nurbtype
    integer         :: ileafcod(maxubtp)
    real            :: z0_town(maxubtp)
    real            :: bld(maxubtp)
    real            :: bld_height(maxubtp)
    real            :: bld_hl_ratio(maxubtp)
    real            :: aroof(maxubtp)
    real            :: eroof(maxubtp)
    real            :: aroad(maxubtp)
    real            :: eroad(maxubtp)
    real            :: awall(maxubtp)
    real            :: ewall(maxubtp)
    real            :: htraf(maxubtp)
    real            :: hindu(maxubtp)
    real            :: pletraf(maxubtp)
    real            :: pleindu(maxubtp)

    namelist /TEB_SPM_INFO/                                              &
         teb_spm,                                                        &
         fusfiles, ifusflg, ifusfn,                                      &
         ichemi, ichemi_in, chemdata_in, isource, weekdayin, rushh1,     &
         rushh2, daylight, efsat, efsun, eindno, eindno2, eindpm,        &
         eindco, eindso2, eindvoc, eveino, eveino2, eveipm, eveico,      &
         eveiso2, eveivoc, iteb, tminbld, nteb, hc_roof, tc_roof,        &
         d_roof, hc_road, tc_road, d_road, hc_wall, tc_wall, d_wall,     &
         nurbtype, ileafcod, z0_town, bld, bld_height, bld_hl_ratio,     &
         aroof, eroof, aroad, eroad, awall, ewall, htraf, hindu,         &
         pletraf, pleindu


    !namelist /MODEL_FILE_INFO/

    integer :: initial
    integer :: nud_type
    character(len=f_name_length)       :: varfpfx
    real :: vwait1
    real :: vwaittot
    character(len=f_name_length) :: nud_hfile
    integer :: nudlat
    real :: timeWindowIAU
    real :: ramp
    real :: tnudlat
    real :: tnudcent
    real :: tnudtop
    real :: znudtop
    real :: wt_nudge_grid(maxgrds)
    real :: wt_nudge_uv
    real :: wt_nudge_th
    real :: wt_nudge_pi
    real :: wt_nudge_rt
    
    integer :: applyIAU
    character(len=f_name_length) :: fileNameIAU

    integer :: nud_cond
    character(len=f_name_length) :: cond_hfile
    real :: tcond_beg
    real :: tcond_end
    real :: t_nudge_rc
    real :: wt_nudgec_grid(maxgrds)
    integer            :: if_oda
    character(len=128) :: oda_upaprefix
    character(len=128) :: oda_sfcprefix
    real               :: frqoda
    real               :: todabeg
    real               :: todaend
    real               :: tnudoda
    real               :: wt_oda_grid(maxodagrids)
    real               :: wt_oda_uv
    real               :: wt_oda_th
    real               :: wt_oda_pi
    real               :: wt_oda_rt
    real               :: roda_sfce(maxodagrids)
    real               :: roda_sfc0(maxodagrids)
    real               :: roda_upae(maxodagrids)
    real               :: roda_upa0(maxodagrids)
    real               :: roda_hgt(maxodagrids)
    real               :: roda_zfact(maxodagrids)
    real               :: oda_sfc_til
    real               :: oda_sfc_tel
    real               :: oda_upa_til
    real               :: oda_upa_tel
    integer :: if_cuinv
    character(len=128) :: cu_prefix
    real :: tnudcu
    real :: wt_cu_grid(maxcugrids)
    real :: tcu_beg
    real :: tcu_end
    real :: cu_tel
    real :: cu_til
    real :: timstr
    character(len=f_name_length) :: hfilin
    integer :: ipastin
    character(len=f_name_length) :: pastfn
    integer :: ioutput
    character(len=f_name_length) :: hfilout
    character(len=f_name_length) :: afilout
    integer :: iclobber
    integer :: ihistdel
    real :: frqhis
    real :: frqanl
    real :: frqlite
    integer :: ipos
    character(len=20)  :: xlite
    character(len=20)  :: ylite
    character(len=20)  :: zlite
    integer :: nlite_vars
    character(len=32)  :: lite_vars(maxlite)
    real :: avgtim
    real :: frqmean
    real :: frqboth
    integer :: kwrite
    real             :: frqprt
    integer            :: initfld
    integer                :: prtcputime
    character(len=f_name_length) :: topfiles
    character(len=f_name_length) :: sfcfiles
    character(len=f_name_length) :: sstfpfx
    character(len=f_name_length) :: ndvifpfx
    integer :: itoptflg(maxgrds)
    integer :: isstflg(maxgrds)
    integer :: ivegtflg(maxgrds)
    integer :: isoilflg(maxgrds)
    integer :: ndviflg(maxgrds)
    integer :: nofilflg(maxgrds)
    integer            :: iupdndvi
    integer            :: iupdsst
    character(len=f_name_length) :: itoptfn(maxgrds)
    character(len=f_name_length) :: isstfn(maxgrds)
    character(len=f_name_length) :: ivegtfn(maxgrds)
    character(len=f_name_length) :: isoilfn(maxgrds)
    character(len=f_name_length) :: ndvifn(maxgrds)
    integer :: itopsflg(maxgrds)
    real  :: toptenh(maxgrds)
    real  :: toptwvl(maxgrds)
    integer :: iz0flg(maxgrds)
    real  :: z0max(maxgrds)
    real  :: z0fact
    integer            :: mkcoltab
    character(len=f_name_length) :: coltabfn
    character(len=f_name_length) :: mapaotfile
    character(len=f_name_length) :: julesin


    namelist /MODEL_FILE_INFO/                                           &
         initial, varfpfx, nudlat,tnudlat, tnudcent, tnudtop, znudtop,   &
         ioutput, hfilout, afilout, frqhis, frqanl, ipos, topfiles,      &
         sfcfiles, sstfpfx, ndvifpfx, itoptfn, isstfn, ivegtfn, isoilfn, &
         ndvifn

    namelist /MODEL_FILE_INFO2/                                           &
         nud_type, varfpfx, vwait1, vwaittot, nud_hfile, &
         timeWindowIAU,ramp, wt_nudge_grid, wt_nudge_uv,&
         wt_nudge_th, wt_nudge_pi, wt_nudge_rt, applyIAU,fileNameIAU,    &
         nud_cond, cond_hfile,    &
         tcond_beg, tcond_end, t_nudge_rc, wt_nudgec_grid, if_oda,       &
         oda_upaprefix,oda_sfcprefix, frqoda, todabeg, todaend, tnudoda, &
         wt_oda_grid, wt_oda_uv, wt_oda_th, wt_oda_pi, wt_oda_rt,        &
         roda_sfce, roda_sfc0, roda_upae,roda_upa0, roda_hgt,            &
         roda_zfact, oda_sfc_til, oda_sfc_tel, oda_upa_til, oda_upa_tel, &
         if_cuinv, cu_prefix, tnudcu, wt_cu_grid, tcu_beg, tcu_end,      &
         cu_tel, cu_til, timstr, hfilin, ipastin, pastfn, iclobber,      &
         ihistdel, frqlite, xlite, ylite, zlite, nlite_vars, lite_vars,  &
         avgtim, frqmean, frqboth, kwrite, frqprt, initfld, prtcputime,  &
         itoptflg, isstflg, ivegtflg,isoilflg, ndviflg, nofilflg,        &
         iupdndvi, iupdsst, itopsflg, toptenh, toptwvl, iz0flg,   &
         z0max, z0fact, mkcoltab, coltabfn, mapaotfile, julesin

    ! namelist /MODEL_OPTIONS/

    integer :: naddsc
    integer :: icorflg
    integer :: dyncore_flag
    integer :: pd_or_mnt_constraint
    integer :: order_h
    integer :: order_v

!--(DMK-CCATT-INI)-----------------------------------------------------------
    integer :: iexev
    integer :: imassflx
    integer :: vveldamp
!--(DMK-CCATT-FIM)-----------------------------------------------------------
    integer :: ibnd
    integer :: jbnd
    real :: cphas
    integer :: lsflg
    integer :: nfpt
    real :: distim
    integer :: sfire
    integer :: iswrtyp
    integer :: ilwrtyp
    character(LEN=f_name_length) :: raddatfn
    real    :: radfrq
    real    :: radtun
    integer :: lonrad
    integer :: nnqparm(maxgrds)
    character (len=2) :: closure_type
    integer :: nnshcu(maxgrds)
    real :: confrq
    real    :: shcufrq
    real :: wcldbs
    integer :: g3d_spread
    integer :: g3d_smoothh
    integer :: g3d_smoothv
    integer :: npatch
    integer :: nvegpat
    integer :: isfcl
    integer :: isfcl_ocean
    integer :: nvgcon
    real    :: pctlcon
    integer :: nslcon
    real    :: drtcon
    real    :: zrough
    real    :: albedo
    real    :: seatmp
    real    :: dthcon
    character(len=1)   :: soil_moist
    character(len=1)   :: soil_moist_fail
    character(len=f_name_length) :: usdata_in
    character(len=f_name_length) :: usmodel_in
    real    :: slz(nzgmax)
    real    :: slmstr(nzgmax)
    real    :: stgoff(nzgmax)
    integer :: if_urban_canopy
    integer :: idiffk(maxgrds)
    integer :: ihorgrad
    real    :: csx(maxgrds)
    real    :: csz(maxgrds)
    real    :: xkhkm(maxgrds)
    real    :: zkhkm(maxgrds)
    real    :: akmin(maxgrds)
    integer            :: mcphys_type
    integer            :: idriz
    integer	       :: iccnlev
    integer	       :: irime
    integer	       :: iplaws
    integer            :: level
    integer            :: icloud
    integer            :: irain
    integer            :: ipris
    integer            :: isnow
    integer            :: iaggr
    integer            :: igraup
    integer            :: ihail
    real               :: cparm
    real               :: rparm
    real               :: pparm
    real               :: sparm
    real               :: aparm
    real               :: gparm
    real               :: hparm
    real               :: dparm
    real	       :: cnparm
    real	       :: gnparm
    real	       :: epsil
    real               :: gnu(ncat_dummy)
    INTEGER            :: advmnt
    INTEGER            :: GhostZoneLength

    namelist /MODEL_OPTIONS/ &
         sfire, iswrtyp, ilwrtyp,radfrq, nnqparm, closure_type,       &
         nnshcu, confrq, shcufrq,  isfcl, isfcl_ocean, soil_moist_fail, &
         usdata_in, usmodel_in, mcphys_type, level


    namelist /MODEL_OPTIONS2/ &
         dyncore_flag, pd_or_mnt_constraint, order_h, order_v ,advmnt,  &
         GhostZoneLength,naddsc, icorflg,iexev, imassflx, vveldamp,     &
         ibnd, jbnd, cphas, lsflg, nfpt, distim,       &
         raddatfn, radtun, lonrad,wcldbs, g3d_spread, g3d_smoothh,      &
         g3d_smoothv, npatch, nvegpat, nvgcon,      &
         pctlcon, nslcon, drtcon, zrough, albedo, seatmp, dthcon,       &
         soil_moist, slz,       &
         slmstr, stgoff, if_urban_canopy, idiffk, ihorgrad, csx, csz,   &
         xkhkm, zkhkm, akmin, idriz, icloud, irain, &
         ipris, isnow, iaggr, igraup, ihail,irime, iplaws,iccnlev,      &
         cparm, rparm, pparm, sparm, aparm, gparm, hparm,  dparm,cnparm,&
         epsil,gnparm,gnu,windfarm,wfFile,damModule,frqPrecip,          &
         damOutPrefix,evaluate,evaluatePrefix


    !namelist /MODEL_SOUND/

    integer             :: ipsflg
    integer             :: itsflg
    integer             :: irtsflg
    integer             :: iusflg
    real                :: hs(maxsndg)
    real                :: ps(maxsndg)
    real                :: ts(maxsndg)
    real                :: rts(maxsndg)
    real                :: us(maxsndg)
    real                :: vs(maxsndg)

    namelist /MODEL_SOUND/ &
         ipsflg, itsflg, irtsflg, iusflg, hs, ps, ts, rts, us, vs

    !namelist /MODEL_PRINT/

    integer            :: nplt
    character(len=16)  :: iplfld(50)
    integer            :: ixsctn(50)
    integer            :: isbval(50)

    namelist /MODEL_PRINT/ &
         nplt, iplfld, ixsctn, isbval

    !namelist /ISAN_CONTROL/

    integer :: iszstage
    integer :: ivrstage
    integer :: isan_inc
    character(len=8)   :: guess1st
    integer :: i1st_flg
    integer :: iupa_flg
    integer :: isfc_flg
    character(len=256) :: iapr
    character(len=256) :: iarawi
    character(len=256) :: iasrfce
    character(len=256) :: varpfx
    integer :: ioflgisz
    integer :: ioflgvar

    namelist /ISAN_CONTROL/ &
          isan_inc,iapr,varpfx

    namelist /ISAN_CONTROL2/ &
         iszstage, ivrstage, guess1st, i1st_flg, iupa_flg,       &
         isfc_flg, iarawi, iasrfce, ioflgisz, ioflgvar


    !namelist /ISAN_ISENTROPIC/

    integer                    :: nisn
    integer                    :: levth(maxisn)
    integer                    :: nigrids
    real                       :: topsigz
    real                       :: hybbot
    real                       :: hybtop
    real                       :: sfcinf
    real                       :: sigzwt
    integer                    :: nfeedvar
    integer                    :: maxsta
    integer                    :: maxsfc
    integer                    :: notsta
    character(len=8)           :: notid(50)
    integer                    :: iobswin
    real                       :: stasep
    integer                    :: igridfl
    real                       :: gridwt(maxagrds)
    real                       :: gobsep
    real                       :: gobrad
    real                       :: wvlnth(maxagrds)
    real                       :: respon(maxagrds)
    real                       :: swvlnth(maxagrds)
    !
    integer                    :: icFileType
    character(len=256)         :: icPrefix
    character(len=32)          :: wind_u_varname  
    character(len=32)          :: wind_v_varname 
    character(len=32)          :: temperature_varname
    character(len=32)          :: geo_varname
    character(len=32)          :: ur_varname    
    real                       :: initial_latitude
    real                       :: final_latitude
    real                       :: initial_longitude
    real                       :: final_longitude
    integer                    :: z_max_level
    real                       :: dlimit(5)
    real                       :: ulimit(5)
    real                       :: scale_factor(5)
    character(len=256)         :: icGradsPrefix
    integer                    :: ccGradsWrite
    include "constants.f90"


    namelist /ISAN_ISENTROPIC/ &
         icFileType, icPrefix, &
         wind_u_varname, wind_v_varname, temperature_varname, geo_varname, &
         ur_varname, initial_latitude, final_latitude, initial_longitude, &
         final_longitude, z_max_level, scale_factor

    namelist /ISAN_ISENTROPIC2/ &
         nisn, levth, nigrids, topsigz, hybbot, hybtop, sfcinf, sigzwt,    &
         nfeedvar, maxsta, maxsfc, notsta, notid, iobswin, stasep, igridfl,&
         gridwt, gobsep, gobrad, wvlnth, swvlnth, respon, dlimit, ulimit,  &
         ccGradsWrite, icGradsPrefix

    !namelist POST

    integer                       :: nvp
    character(len=f_name_length)  :: vp(200)
    character(len=f_name_length)  :: gprefix
    character(len=f_name_length)  :: csvFile
    character(len=10)             :: anl2gra
    character(len=10)             :: proj
    character(len=3)              :: mean_type
    real                          :: lati(maxgrds)
    real                          :: loni(maxgrds)
    real                          :: latf(maxgrds)
    real                          :: lonf(maxgrds)
    integer                       :: zlevmax(maxgrds)
    integer                       :: ipresslev
    integer                       :: inplevs
    integer                       :: iplevs(nzpmax)
    character(len=10)             :: mechanism
    character(len=20)             :: ascii_data
    real                          :: site_lat
    real                          :: site_lon
    namelist /POST/ &
         nvp, vp, gprefix, csvFile, anl2gra, proj, mean_type, lati, loni, &
         latf, lonf, zlevmax, ipresslev, inplevs, iplevs, &
         mechanism, ascii_data, site_lat, site_lon

	!namelist digital filter
	logical :: applyDigitalFilter
	real	:: digitalFilterTimeWindow
     namelist /DIGITALFILTER/ &
	       applyDigitalFilter, digitalFilterTimeWindow


     logical                      :: applyMeteogram
     real                         :: meteogramFreq
     character(len=f_name_length) :: meteogramMap
     character(len=f_name_length) :: meteogramDir

     namelist /METEOGRAM/ &
     applyMeteogram,      &
     meteogramFreq,       &
     meteogramMap,        &
     meteogramDir

!--(DMK-CCATT-INI)-----------------------------------------------------------
    ! CCATT_INFO

    advanced_ramsin = './RAMSIN_ADVANCED'

    ccatt               = 0
    chemistry           = 0
    split_method        = ''
    chem_timestep       = 0.
    chemistry_aq        = 0
!    aerosol             = 0
    chem_assim          = 0
    srcmapfn            = ''
    recycle_tracers     = 0
    def_proc_src        = 'STOP'
    diur_cycle(1:nsrc)  = 1
    na_extra2d          = 0
    na_extra3d          = 0
    plumerise           = 0
    prfrq               = 3600. ! Initial Value for PlumeRise Frequency - CCATT
    volcanoes           = 0
!Matrix
    aerosol             = 0
    aer_timestep        = 0.
    aer_assim           = 0
    mech                = 8
!--(DMK-CCATT-OLD)-----------------------------------------------------------
!    ! CATT_INFO
!    catt                = 0
!    firemapfn           = ''
!    recycle_tracers     = 0
!    plumerise           = 0
!    define_proc         = 0
!    prfrq               = 3600. ! Initial Value for PlumeRise Frequency - CATT
!--(DMK-CCATT-FIM)-----------------------------------------------------------

    ! ISAN_CONTROL
    iszstage            = 1
    ivrstage            = 1
    isan_inc            = 0600
    guess1st	      = 'PRESS'
    i1st_flg	      = 1
    iupa_flg	      = 3
    isfc_flg	      = 3
    iapr		      = './dprep/dp' ! 2
    iarawi	      = ''
    iasrfce	      = ''
    varpfx	      = './ivar/iv-brams4' ! 2
    ioflgisz 	      = 0
    ioflgvar 	      = 1

    ! ISAN_ISENTROPIC
    nisn                = 43
    levth               = 800
    levth(1:nisn)       = (/280,282,284,286,288,290,292,294,296,298,300,303,306,309,&
         312,315,318,321,324,327,330,335,340,345,350,355,360,380,400,420,440, &
         460,480,500,520,540,570,600,630,670,700,750,800/)
    nigrids             = 1
    topsigz	      = 20000.
    hybbot	      = 4000.
    hybtop	      = 6000.
    sfcinf	      = 1000.
    sigzwt	      = 1.
    nfeedvar            = 1.
    maxsta 	      = 150
    maxsfc 	      = 1000
    notsta 	      = 0
    notid               = ''
    iobswin	      = 1800
    stasep	      = .1
    igridfl	      = 3
    gridwt	      = .01
    gobsep	      = 5.
    gobrad	      = 5.
    wvlnth	      = 1000
    swvlnth	      = 500.
    respon	      = .90
    icFileType         = 0
    icPrefix   = '../GFS/gfs.t00z.pgrb2.0p25.f'
    wind_u_varname  = "UGRD"
    wind_v_varname  = "VGRD"
    temperature_varname = "TMP"
    geo_varname = "HGT"
    ur_varname =  "RH"  
    initial_latitude = -70.
    final_latitude = 29.
    initial_longitude = 250.
    final_longitude = 358.
    z_max_level = 26
    dlimit = (/500.,-500.,   0.,  -1000.,   0./)
    ulimit = (/500., 500., 500., 500000., 100./)
    scale_factor = (/1.0,1.0,1.0,1.0,1.0/)
    ccGradsWrite = 0
    icGradsPrefix = './icGrads'

    ! MODEL_FILE_INFO
    initial	      = 2 ! 2
    nud_type	      = 2 ! 2
    varfpfx	      = varpfx ! will be rewrited on the namelist read, and after too
    vwait1	      = 0.
    vwaittot	      = 0.
    nud_hfile	      = '' ! 2
    nudlat	      = 5 ! 2
    timeWindowIAU     = 0.0
    ramp              = 0.0
    tnudlat	      = 1800. ! 2
    tnudcent	      = 0. ! 2
    tnudtop	      = 10800. ! 2
    znudtop	      = 16000. ! 2
    wt_nudge_grid       = 1. ! 2
    wt_nudge_grid(1:4)  = (/1.,0.8,0.7,0.5/) ! 2
    wt_nudge_uv	      = 1.
    wt_nudge_th	      = 1.
    wt_nudge_pi	      = 1.
    wt_nudge_rt	      = 1.
    applyIAU          = 0
    fileNameIAU       = 'IAUFILE'
 
 
    nud_cond	      = 0
    cond_hfile	      = ''
    tcond_beg	      = 0.
    tcond_end	      = 21600.
    t_nudge_rc	      = 3600.
    wt_nudgec_grid      = 0.5
    wt_nudgec_grid(1:4) = (/1.,0.8,0.7,0.5/)
    if_oda	      = 0
    oda_upaprefix       = ''
    oda_sfcprefix       = ''
    frqoda 	      = 300.
    todabeg	      = 0.
    todaend	      = 9999999.
    tnudoda	      = 900.
    wt_oda_grid         = 1.
    wt_oda_grid(1:4)    = (/1.,0.8,0.7,0.5/)
    wt_oda_uv	      = 1.
    wt_oda_th	      = 1.
    wt_oda_pi	      = 1.
    wt_oda_rt	      = 1.
    roda_sfce	      = 0.
    roda_sfce(1:4)      = (/50000.,100.,100.,100./)
    roda_sfc0           = 0.
    roda_sfc0(1:4)      = (/100000.,100000.,100000.,100000./)
    roda_upae           = 0.
    roda_upae(1:4)      = (/100000.,200.,200.,200./)
    roda_upa0           = 0.
    roda_upa0(1:4)      = (/200000.,2000.,2000.,2000./)
    roda_hgt            = 0.
    roda_hgt(1:4)       = (/3000.,3000.,3000.,3000./)
    roda_zfact          = 0.
    roda_zfact(1:4)     = (/100.,100.,100.,100./)
    oda_sfc_til 	      = 21600.
    oda_sfc_tel 	      = 900.
    oda_upa_til 	      = 43200.
    oda_upa_tel 	      = 21600.
    oda_upa_tel 	      = 21600.
    if_cuinv            = 0
    cu_prefix           = ''
    tnudcu              = 900.
    wt_cu_grid          = 1.
    wt_cu_grid(1:4)     = (/1.,1.,0.5,0.5/)
    tcu_beg	      = 0.
    tcu_end	      = 7200.
    cu_tel 	      = 3600.
    cu_til 	      = 21600.
    timstr 	      = 0
    hfilin 	      = ''
    ipastin	      = 0 ! 2
    pastfn 	      = '' ! 2
    ioutput	      = 2 ! 2
    hfilout	      = './H/H-brams4' ! 2
    afilout	      = './A/A-brams4' ! 2
    iclobber 	      = 0 ! 2
    ihistdel 	      = 0 ! 2
    ipos              = 0
    frqhis	      = 21600. ! 2
    frqanl	      = 10800. ! 2
    frqlite             = 0.
    xlite 	      = '/0:0/'
    ylite 	      = '/0:0/'
    zlite 	      = '/0:0/'
    nlite_vars	      = 4
    lite_vars 	      = ''
    lite_vars(1)        = 'UP'
    lite_vars(2)        = 'VP'
    lite_vars(3)        = 'WP'
    lite_vars(4)        = 'THETA'
    avgtim 	      = 0.
    frqmean	      = 0.
    frqboth	      = 0.
    kwrite 	      = 0
    frqprt 	      = 21600.
    initfld	      = 1
    topfiles	      = './data/toph-'
    sfcfiles	      = './data/sfc-'
    sstfpfx 	      = './data/sst-'
    ndvifpfx	      = './data/ndvi-'
    itoptflg	      = 2 ! 2
    itoptflg(1:4)       = (/2,2,2,2/) ! 2
    isstflg             = 2 ! 2
    isstflg(1:4)        = (/2,2,2,2/) ! 2
    ivegtflg            = 2 ! 2
    ivegtflg(1:4)       = (/2,2,2,2/) ! 2
    isoilflg            = 2 ! 2
    isoilflg(1:4)       = (/2,2,2,2/) ! 2
    ndviflg             = 2 ! 2
    ndviflg(1:4)        = (/2,2,2,2/) ! 2
    nofilflg            = 2
    nofilflg(1:4)       = (/2,2,2,2/)
    iupdndvi            = 1
    iupdsst	      = 1
    itoptfn	      = ''
    itoptfn(1)          = './topo10km/H'
    itoptfn(2:4)        = (/'./topo/EL','./topo/EL','./topo/EL'/)
    isstfn              = ''
    isstfn(1:4)         = (/'./sst/S','./sst/S','./sst/S','./sst/S'/)
    ivegtfn             = ''
    ivegtfn(1:4)        = (/'./soil-fao/FAO','./soil-fao/FAO','./soil-fao/FAO','./soil-fao/FAO'/)
    isoilfn             = ''
    isoilfn(1:4)        = (/'./oge_brams4/OGE','./oge_brams4/OGE','./oge_brams4/OGE','./oge_brams4/OGE'/)
    ndvifn              = ''
    ndvifn(1:4)         = (/'./ndvi-modis/N','./ndvi-modis/N','./ndvi-modis/N','./ndvi-modis/N'/)
    itopsflg            = 0
    itopsflg(1:4)       = (/0,0,0,0/)
    toptenh             = 1
    toptenh(1:4)        = (/1,1,1,1/)
    toptwvl             = 1
    toptwvl(1:4)        = (/1,1,1,1/)
    iz0flg              = 0
    iz0flg(1:4)         = (/0,0,0,0/)
    z0max               = 5.
    z0max(1:4)          = (/5.,5.,5.,5./)
    z0fact              = 0.005
    mkcoltab            = 0 ! duvida!
    coltabfn            = './micro/ct2.0' ! duvida!
    mapaotfile          = './tables2/rad_carma/infMapAOT.vfm'
    julesin             = './julesin'
    prtcputime          = 0

    ! MODEL_GRIDS
    expnme              = 'BRAMS 41' ! 2
    timeunit 	      = 'h' ! 2
    load_bal 	      = 0
    ngrids 	      = 1 ! 2
    nzg    	      = 9 ! 2
    nzs    	      = 32
    nxtnest	      = 1 ! 2
    nxtnest(1:4)        = (/0,1,2,3/) ! 2
    if_adap	      = 0 ! 2
    ihtran 	      = 1
    deltaz 	      = 250 ! 2
    dzrat  	      = 1.2 ! 2
    dzmax  	      = 1000.
    fixLevels        = 0
    zz     	      = 19700.
    zz(1:41)            = (/0.0,20.0,46.0,80.0,120.0,165.0,220.0,290.0,380.0,480.0,590.0, &
         720.0,870.0,1030.0,1200.0,1380.0,1595.0,1850.0,2120.0,2410.0,2715.0,  &
         3030.0,3400.0,3840.0,4380.0,5020.0,5800.0,6730.0,7700.0,8700.0,9700.0,&
         10700.,11700.,12700.,13700.,14700.,15700.,16700.,17700.,18700.,19700./)
    dtlong  	      = 30. ! 2
    nacoust 	      = 3
    ideltat 	      = 1 ! 2
    nstratx 	      = 3
    nstratx(1:4)        = (/1,3,3,3/)
    nstraty             = 3
    nstraty(1:4)        = (/1,3,3,3/)
    nndtrat             = 3
    nndtrat(1:4)        = (/1,3,3,3/)
    nestz1              = 0
    nstratz1            = 1
    nstratz1(1:4)       = (/2,2,2,1/)
    nestz2              = 0
    nstratz2            = 1
    nstratz2(1:4)       = (/3,3,2,1/)
    ninest              = 3 ! 2
    ninest(1:4)         = (/1,1,2,3/) ! 2
    njnest              = 3 ! 2
    njnest(1:4)         = (/1,1,2,3/) ! 2
    nknest              = 1
    nknest(1:4)         = (/1,1,1,1/)
    nnsttop             = 1
    nnsttop(1:4)        = (/1,1,1,1/)
    nnstbot             = 1
    nnstbot(1:4)        = (/1,1,1,1/)
    gridu               = 0.
    gridu(1:4)          = (/0.,0.,0.,0./)
    gridv               = 0.
    gridv(1:4)          = (/0.,0.,0.,0./)

!--(DMK-CCATT-INI)-----------------------------------------------------------
    dyncore_flag        = 0 !default value
    pd_or_mnt_constraint =0
    order_h = 3
    order_v = 3
    advmnt              = 0
    ghostzonelength     = 1
    windfarm = 0
    wfFile = ''
    damModule = 0
    frqPrecip = 3600.
    damOutPrefix = './damOutput'
    evaluate = 0
    evaluatePrefix = './statistic-'
!--(DMK-CCATT-FIM)-----------------------------------------------------------

    domain_fname        = ''

    ! MODEL_OPTIONS
    naddsc  	      = 0
    icorflg 	      = 1
!--(DMK-CCATT-INI)-----------------------------------------------------------
    iexev             = 1
    imassflx          = 0
    vveldamp          = 0
!--(DMK-CCATT-FIM)-----------------------------------------------------------
    ibnd    	      = 1
    jbnd    	      = 1
    cphas   	      = 20.
    lsflg   	      = 0
    nfpt    	      = 0
    distim  	      = 400.
    sfire            = 0
    iswrtyp 	      = 1
    ilwrtyp 	      = 1
    raddatfn	      = "./carma/rad_param.data" ! 2
    radfrq  	      = 900.
    radtun  	      = 1. 
    lonrad  	      = 1
    nnqparm 	      = 2
    closure_type      = "EN" ! 2
    nnshcu            = 1 ! 2
    nnshcu(1:4)       = (/1,1,1,1/) !2
    confrq            = 600. ! 2
    shcufrq           = 600. ! 2
    wcldbs            = .0005
    g3d_spread	      = 0
    g3d_smoothh	      = 0
    g3d_smoothv	      = 0
    npatch            = 2 ! 2
    nvegpat           = 1
    isfcl             = 1 ! 2
    isfcl_ocean       = -999 !-srf: do not change this.
    nvgcon            = 6 ! 2
    pctlcon           = 1.
    nslcon 	      = 6 ! 2
    zrough 	      = .05
    albedo 	      = .2
    seatmp 	      = 298. ! 2
    dthcon 	      = 0.
    drtcon 	      = 0.
    soil_moist          = "i" ! 2
    soil_moist_fail     = "l" ! 2
    usdata_in  	      = "./soil-moisture/GL_SM.GMNR." ! 2
    usmodel_in 	      = "./data/SM-" ! 2
    slz        	      = -0.1 ! 2
    slz(1:9)   	      = (/-2.0,-1.75,-1.50,-1.25,-1.00,-0.75,-0.50,-0.25,-0.1/)
    slmstr     	      = 0.28 ! 2
    slmstr(1:9)	      = (/0.40,0.37,0.35,0.33,0.32,0.31,0.30,0.29,0.28/) ! 2
    stgoff     	      = 0.0
    stgoff(1:9)	      = (/0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/)
    if_urban_canopy     = 0
    idiffk              = 1 ! 2
    idiffk(1:4)         = (/1,1,1,1/) ! 2
    ihorgrad            = 2
    csx                 = .2
    csx(1:4) 	      = (/.2,.2,.2,.2/)
    csz      	      = .2
    csz(1:4) 	      = (/.35,.35,.35,.35/)
    xkhkm    	      = 3.
    xkhkm(1:4) 	      = (/3.,3.,3.,3./)
    zkhkm      	      = 3.
    zkhkm(1:4) 	      = (/3.,3.,3.,3./)
    akmin      	      = 1.
    akmin(1:4) 	      = (/1.,1.,1.,1./)
    mcphys_type       = 0
    level  	      = 3 ! 2
    icloud 	      = 4 ! 2
    idriz 	      = 1 ! 2
    irime             = 0
    iplaws            = 0
    iccnlev           = 0
    irain  	      = 2 ! 2
    ipris  	      = 5
    isnow  	      = 2 ! 2
    iaggr  	      = 2 ! 2
    igraup 	      = 2 ! 2
    ihail  	      = 2 ! 2
    cparm  	      = .1e9
    rparm  	      = 1e-3
    pparm  	      = 0.
    sparm  	      = 1e-3
    aparm  	      = 1e-3
    gparm  	      = 1e-3
    hparm  	      = 3e-3
    dparm             = 1.e-5
    cnparm            = 0.04e-4
    gnparm            = 3.00e-4
    epsil             = 0.1
    gnu    	      = 2.
    gnu(1:8)	      = (/2.,2.,2.,2.,2.,2.,2.,2./)
    windfarm          = 0
    wfFile           = ''
    damModule        = 0
    frqPrecip = 3600.
    damOutPrefix = './damOutput'
    evaluate = 0
    evaluatePrefix = 'statistic-'

    ! MODEL_PRINT
    nplt   	      = 0
    iplfld 	      = ""
    ixsctn 	      = 3
    ixsctn(1:4)         = (/3,3,3,3/)
    isbval              = 2
    isbval(1:4)         = (/2,2,2,2/)

    ! MODEL_SOUND
    ipsflg  	      = 1 ! 2
    itsflg  	      = 0 ! 2
    irtsflg 	      = 3 ! 2
    iusflg  	      = 0 ! 2
    hs      	      = 0. ! 2
    ps(1:11)	      = (/1010.,1000.,2000.,3000.,4000.,6000.,8000.,11000.,15000.,20000.,25000./) ! 2
    ts(1:11)	      = (/25.,18.5,12.,4.5,-11.,-24.,-37.,-56.5,-56.5,-56.5,-56.5/) ! 2
    rts(1:11) 	      = (/70.,70.,70.,70.,20.,20.,20.,20.,10.,10.,10./) ! 2
    us(1:11)  	      = (/10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10./) ! 2
    vs(1:11)  	      = (/0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0./) ! 2

    ! TEB
    teb_spm  	      = 0
    fusfiles 	      = ''
    ifusflg  	      = 0
    ifusflg(1:4)        = (/1,1,1,1/)
    ifusfn 	      = ''
    ifusfn(1:4)	      = (/'./fusos/fuso','./fusos/fuso','./fusos/fuso','./fusos/fuso'/) ! 2
    ichemi     	      = 0    !Photochemical module activation - 1=on, 0=off
    ichemi_in  	      = 1    !Use initial values from previous run (1=yes,0=no) ! 2
    chemdata_in	      = ''
    isource    	      = 1    !Emission module activation - 1=on, 0=off ! 2
    weekdayin  	      = 'SUN'  !Initial weeakday of the simulation
    rushh1     	      = 7.81  !Morning Rush Hour (Local Time in Hours)
    rushh2     	      = 16.0  !Afternoon/Evening Rush Hour (Local Time)
    daylight   	      = 0.    !Daylight saving time (horario de verao)
    ! Emission factor (fraction of weekdays) for Saturdays and Sundays
    ! They are used in the emission module and TEB. - EDF
    efsat     	      = 0.8
    efsun     	      = 0.5
    ! Input GMT difference time variable (To define local time)
    !Industrial emissions (kg/s/m2)
    eindno     	      = 2.6636227e-10
    eindno2    	      = 2.9595805e-11
    eindpm     	      = 4.3421278e-10
    eindco     	      = 8.1599860e-10
    eindso2    	      = 3.6149164e-10
    eindvoc    	      = 2.5367833e-10
    !Veicular emissions (kg/day/m2)
    eveino     	      = 4.3196708e-04
    eveino2    	      = 6.8566209e-05
    eveipm     	      = 6.2648396e-06
    eveico     	      = 7.5433785e-03
    eveiso2    	      = 4.0730592e-05
    eveivoc    	      = 1.1892237e-03
    !----- Urban canopy parameterization using TEB (Masson, 2000)-------------
    iteb      	      = 0     !1=on, 0=off
    tminbld   	      = 12.   !Minimum internal building temperature (degrees Celsius)
    nteb      	      = 3     !Number of roof,road and wall layers used in TEB, Max.3
    ! ROOF layers properties
    ! heat capacity
    hc_roof             = 0.
    hc_roof(1:3)        = (/2110000.,280000.,290000./)
    ! thermal conductivity
    tc_roof             = 0.
    tc_roof(1:3)        = (/0.41,0.05,0.03/)
    ! depth
    d_roof              = 0.
    d_roof(1:3)         = (/0.05,0.4,0.05/)
    ! ROAD layers properties
    ! heat capacity
    hc_road             = 0.
    hc_road(1:3)        = (/1240000.,1280000.,1280000./)
    ! thermal conductivity 1.01
    tc_road             = 1.0103
    ! depth
    d_road              = 0.
    d_road(1:3)         = (/0.05,0.1,1.0/)
    ! WALL layers properties
    ! heat capacity J/m3/K 10e6
    hc_wall             = 1000000.
    ! thermal conductivity 0.81 W/m/K
    tc_wall             = 0.81
    ! depth
    d_wall              = 0.
    d_wall(1:3)         = (/0.02,0.125,0.02/)
    nurbtype            = 2	  !Number of urban types (maximum of 3)

    !Leaf class code to identify each urban type
    ileafcod            = 19
    ileafcod(1:2)       = (21,19)
    !Urban type properties
    !Urban type roughness length 5 e 1
    z0_town             = 0.0
    z0_town(1:2)        = (/3.0,0.5/)
    !Fraction occupied by buildings in the grid cell
    bld      	      = 0.0
    bld(1:2) 	      = (/0.5,0.7/)
    !Building Height
    bld_height          = 0
    bld_height(1:2)     = (/50.,5.0/)
    !Vertical/Horizontal rate 3 e 0.5
    bld_hl_ratio        = 2.4
    bld_hl_ratio(1:2)   = (/4.4,2.4/)
    !Roof albedo
    aroof               = 0.15
    !Roof emissivitiy
    eroof               = 0.9
    !Road albedo
    aroad               = 0.1
    !Road emissivity 90% masson
    eroad               = 0.9
    !Wall albedo
    awall               = 0.25
    !Wall emissivity
    ewall               = 0.85
    !Maximum value of sensible heat
    htraf      	      = 0.
    htraf(1:2) 	      = (/90.0,60.0/)
    !Maximum value of sensible heat
    hindu      	      = 14.
    hindu(1:2) 	      = (/10.0,14.0/)
    !released by Industry (W/m2)
    !Maximum value of latent heat
    pletraf             = 5.
    pletraf(1:2)        = (/10.0,5.0/)
    !released by Traffic (W/m2)
    !Maximum value of latent heat
    pleindu             = 50.
    pleindu(1:2)        = (/30.0,50.0/)

    !PPL defaults for necessary variables
    ! (values from the first-time users)
    ! ATENTION: those variables below should be declared on the simplest RAMSIN
    runtype 	      = "initial"
    timmax  	      = 24
    imonth1 	      = 01
    idate1  	      = 25
    iyear1  	      = 2005
    itime1  	      = 0000
    nnxp    	      = 35
    nnyp    	      = 34
    nnzp    	      = 32
    deltax  	      = 112000.
    deltay  	      = 112000.
    polelat 	      = -23.
    polelon 	      = -52.5
    centlat 	      = -22.
    centlon             = -56.

    !namelist POST
    mechanism = " "
    proj="YES"
    anl2gra="ONE"
    mean_type="VMP"
    lati(:)=-90
    latf(:)=+90
    loni(:)=-180
    lonf(:)=+180
    ascii_data="NO"

    ! PPL end defaults

    ! select unused i/o unit

    do iunit = firstUnit, lastUnit
       inquire(iunit,opened=op)
       if (.not. op) exit
    end do

    if (iunit > lastUnit) then
       call fatal_error(h//" all i/o units in use")
    end if

    ! if namelist file exists, open, read each section and close

    inquire(file=trim(oneNamelistFile%fileName), exist=ex)
    if (.not. ex) then
      print *," namelist file "//trim(oneNamelistFile%fileName)//&
           " does not exist"; call flush(6)
       call fatal_error(" namelist file "//trim(oneNamelistFile%fileName)//&
            " does not exist")

    end if

    open(iunit, file=trim(oneNamelistFile%fileName), status="old", action="read",&
         iostat=err)
    if (err /= 0) then
       write(c0,"(i10)") err
       call fatal_error(h//" open namelist file "//trim(oneNamelistFile%fileName)//&
            " returned iostat="//trim(adjustl(c0)))
    end if

!LFR - INicio da leitura do RAMSIN em partes. Leitura do nome do RAMSIN avançado
    read(iunit, iostat=err, NML=MODEL_ADV_RAMSIN)
    if (err /= 0) then
      write(*,"(a)") h//"**(ERROR)** reading section MODEL_ADV_RAMSIN "//&
            &"of namelist file "//trim(oneNamelistFile%fileName)
      call fatal_error(h//" reading namelist")
    endif 
    inquire(file=trim(advanced_ramsin), exist=ex)
    if (.not. ex) then
      print *," namelist file "//trim(advanced_ramsin)//&
           " does not exist"; call flush(6)
       call fatal_error(" namelist file "//trim(advanced_ramsin)//&
            " does not exist")
    end if
    do iunit2 = firstUnit, lastUnit
       inquire(iunit2,opened=op)
       if (.not. op) exit
    end do
    if (iunit2 > lastUnit) then
       call fatal_error(h//" all i/o units in use")
    end if    
    open(iunit2, file=trim(advanced_ramsin), status="old", action="read",&
         iostat=err)
    if (err /= 0) then
       write(c0,"(i10)") err
       call fatal_error(h//" open namelist file "//trim(advanced_ramsin)//&
            " returned iostat="//trim(adjustl(c0)))
    end if    
!
    read (iunit, iostat=err, NML=MODEL_GRIDS)
    if (err /= 0) then
       write(*,"(a)") h//"**(ERROR)** reading section MODEL_GRIDS "//&
            &"of namelist file "//trim(oneNamelistFile%fileName)
       write(*,"(a)") h//" compare values read with file contents:"
       write(*,*) "expnme=",trim(expnme)
       write(*,*) "runtype=",trim(runtype)
       write(*,*) "timeunit=",trim(timeunit)
       write(*,*) "timmax=",timmax
       write(*,*) "load_bal=",load_bal
       write(*,*) "imonth1=",imonth1
       write(*,*) "idate1=",idate1
       write(*,*) "iyear1=",iyear1
       write(*,*) "itime1=",itime1
       write(*,*) "nnxp=",nnxp
       write(*,*) "nnyp=",nnyp
       write(*,*) "nnzp=",nnzp
       write(*,*) "nzg=",nzg
       write(*,*) "nzs=",nzs
       write(*,*) "deltax=",deltax
       write(*,*) "deltay=",deltay
       write(*,*) "deltaz=",deltaz
       write(*,*) "dzrat=",dzrat
       write(*,*) "dzmax=",dzmax
       write(*,*) "fixLevels=",fixLevels
       write(*,*) "zz=",zz
       write(*,*) "dtlong=",dtlong
       write(*,*) "polelat=",polelat
       write(*,*) "polelon=",polelon
       write(*,*) "centlat=",centlat
       write(*,*) "centlon=",centlon

       call fatal_error(h//" reading namelist")
    else
       ! namelist /MODEL_GRID/
       oneNamelistFile%expnme=expnme
       oneNamelistFile%runtype=runtype
       oneNamelistFile%timeunit=timeunit
       oneNamelistFile%timmax=timmax
       oneNamelistFile%imonth1=imonth1
       oneNamelistFile%idate1=idate1
       oneNamelistFile%iyear1=iyear1
       oneNamelistFile%itime1=itime1
       oneNamelistFile%ngrids=ngrids
       oneNamelistFile%nnxp=nnxp
       oneNamelistFile%nnyp=nnyp
       oneNamelistFile%nnzp=nnzp
       oneNamelistFile%nzg=nzg
       oneNamelistFile%nzs=nzs
       oneNamelistFile%deltax=deltax
       oneNamelistFile%deltay=deltay
       oneNamelistFile%deltaz=deltaz
       oneNamelistFile%dzrat=dzrat
       oneNamelistFile%dzmax=dzmax
       oneNamelistFile%fixLevels=fixLevels
       oneNamelistFile%zz=zz
       oneNamelistFile%dtlong=dtlong
       oneNamelistFile%polelat=polelat
       oneNamelistFile%polelon=polelon
       oneNamelistFile%centlat=centlat
       oneNamelistFile%centlon=centlon
    end if


    read (iunit2, iostat=err, NML=MODEL_GRIDS2)
    if (err /= 0) then
       write(*,"(a)") h//"**(ERROR)** reading section MODEL_GRIDS "//&
            &"of namelist file "//trim(advanced_ramsin)
       write(*,"(a)") h//" compare values read with file contents:"
       write(*,*) "ngrids=",ngrids
       write(*,*) "nxtnest=",nxtnest
       write(*,*) "DOMAIN_FNAME=", DOMAIN_FNAME
       write(*,*) "if_adap=",if_adap
       write(*,*) "ihtran=",ihtran
       write(*,*) "nacoust=",nacoust
       write(*,*) "ideltat=",ideltat
       write(*,*) "nstratx=",nstratx
       write(*,*) "nstraty=",nstraty
       write(*,*) "nndtrat=",nndtrat
       write(*,*) "nestz1=",nestz1
       write(*,*) "nstratz1=",nstratz1
       write(*,*) "nestz2=",nestz2
       write(*,*) "nstratz2=",nstratz2
       write(*,*) "ninest=",ninest
       write(*,*) "njnest=",njnest
       write(*,*) "nknest=",nknest
       write(*,*) "nnsttop=",nnsttop
       write(*,*) "nnstbot=",nnstbot
       write(*,*) "gridu=",gridu
       write(*,*) "gridv=",gridv

       call fatal_error(h//" reading namelist")
    else
       ! namelist /MODEL_GRID/
       oneNamelistFile%ngrids=ngrids
       oneNamelistFile%nxtnest=nxtnest
       oneNamelistFile%domain_fname=domain_fname
       oneNamelistFile%if_adap=if_adap
       oneNamelistFile%ihtran=ihtran
       oneNamelistFile%nacoust=nacoust
       oneNamelistFile%ideltat=ideltat
       oneNamelistFile%nstratx=nstratx
       oneNamelistFile%nstraty=nstraty
       oneNamelistFile%nndtrat=nndtrat
       oneNamelistFile%nestz1=nestz1
       oneNamelistFile%nstratz1=nstratz1
       oneNamelistFile%nestz2=nestz2
       oneNamelistFile%nstratz2=nstratz2
       oneNamelistFile%ninest=ninest
       oneNamelistFile%njnest=njnest
       oneNamelistFile%nknest=nknest
       oneNamelistFile%nnsttop=nnsttop
       oneNamelistFile%nnstbot=nnstbot
       oneNamelistFile%gridu=gridu
       oneNamelistFile%gridv=gridv
    end if



!--(DMK-CCATT-INI)-----------------------------------------------------------
    read (iunit, iostat=err, NML=CCATT_INFO)
    if (err /= 0) then
       write(*,"(a)") h//"**(ERROR)** reading section CCATT_INFO "//&
            &"of namelist file "//trim(oneNamelistFile%fileName)
       write(*,"(a)") h//" compare values read with file contents:"
       print *, "CCATT=", CCATT
       write(*,*) "CHEMISTRY=", chemistry
       write(*,*) "CHEM_TIMESTEP=", chem_timestep
       write(*,*) "CHEM_ASSIM=", chem_assim
       write(*,*) "SRCMAPFN=", srcmapfn
       write(*,*) "AEROSOL=", aerosol
       write(*,*) "AER_TIMESTEP=", aer_timestep
       write(*,*) "AER_ASSIM=", aer_assim

       call fatal_error(h//" reading namelist")
    else
       oneNamelistFile%ccatt=ccatt
       oneNamelistFile%chemistry=chemistry
       oneNamelistFile%chem_timestep=chem_timestep
       oneNamelistFile%chem_assim=chem_assim
       oneNamelistFile%srcmapfn=srcmapfn
       oneNamelistFile%aerosol=aerosol
       oneNamelistFile%aer_timestep=aer_timestep
       oneNamelistFile%aer_assim=aer_assim

    end if


!--(DMK-CCATT-INI)-----------------------------------------------------------
    read (iunit2, iostat=err, NML=CCATT_INFO2)
    if (err /= 0) then
       write(*,"(a)") h//"**(ERROR)** reading section CCATT_INFO "//&
            &"of namelist file "//trim(advanced_ramsin)
       write(*,"(a)") h//" compare values read with file contents:"
       write(*,*) "SPLIT_METHOD=", trim(split_method)
       write(*,*) "CHEMISTRY_AQ=", chemistry_aq
       write(*,*) "RECYCLE_TRACERS=", RECYCLE_TRACERS
       write(*,*) "DEF_PROC_SRC=", trim(def_proc_src)
       write(*,*) "DIUR_CYCLE=", diur_cycle
       write(*,*) "NA_EXTRA2D=", na_extra2d
       write(*,*) "NA_EXTRA3D=", na_extra3d
       write(*,*) "PLUMERISE=", PLUMERISE
       write(*,*) "PRFRQ=", PRFRQ
       write(*,*) "VOLCANOES=", volcanoes
       write(*,*) "MECH=", mech

       call fatal_error(h//" reading namelist")
    else
       oneNamelistFile%split_method=split_method
       oneNamelistFile%chemistry_aq=chemistry_aq
       oneNamelistFile%recycle_tracers=recycle_tracers
       oneNamelistFile%def_proc_src=def_proc_src
       oneNamelistFile%diur_cycle=diur_cycle
       oneNamelistFile%na_extra2d=na_extra2d
       oneNamelistFile%na_extra3d=na_extra3d
       oneNamelistFile%plumerise=plumerise
       oneNamelistFile%prfrq=prfrq
       oneNamelistFile%volcanoes=volcanoes
       oneNamelistFile%mech=mech

    end if

    read (iunit2, iostat=err, NML=TEB_SPM_INFO)
    if (err /= 0) then
       write(*,"(a)") h//"**(ERROR)** reading section TEB_SPM_INFO "//&
            &"of namelist file "//trim(advanced_ramsin)
       write(*,"(a)") h//" compare values read with file contents:"
       print *, "TEB_SPM=", TEB_SPM
       print *, "ifusflg=", ifusflg
       print *, "ifusfn=", ifusfn
       print *, "fusfiles=", trim(fusfiles)
       print *, "ICHEMI=", ICHEMI
       print *, "ICHEMI_IN=", ICHEMI_IN
       print *, "CHEMDATA_IN=", CHEMDATA_IN
       print *, "ISOURCE=", ISOURCE
       print *, "WEEKDAYIN=", trim(WEEKDAYIN)
       print *, "RUSHH1=", RUSHH1
       print *, "RUSHH2=", RUSHH2
       print *, "DAYLIGHT=", DAYLIGHT
       print *, "EFSAT=", EFSAT
       print *, "EFSUN=", EFSUN
       print *, "EINDNO=", EINDNO
       print *, "EINDNO2=", EINDNO2
       print *, "EINDPM=", EINDPM
       print *, "EINDCO=", EINDCO
       print *, "EINDSO2=", EINDSO2
       print *, "EINDVOC=", EINDVOC
       print *, "EVEINO=", EVEINO
       print *, "EVEINO2=", EVEINO2
       print *, "EVEIPM=", EVEIPM
       print *, "EVEICO=", EVEICO
       print *, "EVEISO2=", EVEISO2
       print *, "EVEIVOC=", EVEIVOC
       print *, "ITEB=", ITEB
       print *, "TMINBLD=", TMINBLD
       print *, "NTEB=", NTEB
       print *, "HC_ROOF=", HC_ROOF
       print *, "TC_ROOF=", TC_ROOF
       print *, "D_ROOF=", D_ROOF
       print *, "HC_ROAD=", HC_ROAD
       print *, "TC_ROAD=", TC_ROAD
       print *, "D_ROAD=", D_ROAD
       print *, "HC_WALL=", HC_WALL
       print *, "TC_WALL=", TC_WALL
       print *, "D_WALL=", D_WALL
       print *, "NURBTYPE=", NURBTYPE
       print *, "ILEAFCOD=", ILEAFCOD
       print *, "Z0_TOWN=", Z0_TOWN
       print *, "BLD=", BLD
       print *, "BLD_HEIGHT=", BLD_HEIGHT
       print *, "BLD_HL_RATIO=", BLD_HL_RATIO
       print *, "AROOF=", AROOF
       print *, "EROOF=", EROOF
       print *, "AROAD=", AROAD
       print *, "EROAD=", EROAD
       print *, "AWALL=", AWALL
       print *, "EWALL=", EWALL
       print *, "HTRAF=", HTRAF
       print *, "HINDU=", HINDU
       print *, "PLETRAF=", PLETRAF
       print *, "PLEINDU=", PLEINDU
       call fatal_error(h//" reading namelist")
    else
       oneNamelistFile%teb_spm=teb_spm
       oneNamelistFile%fusfiles=fusfiles
       oneNamelistFile%ifusflg=ifusflg
       oneNamelistFile%ifusfn=ifusfn
       oneNamelistFile%ichemi=ichemi
       oneNamelistFile%ichemi_in=ichemi_in
       oneNamelistFile%chemdata_in=chemdata_in
       oneNamelistFile%isource=isource
       oneNamelistFile%weekdayin=weekdayin
       oneNamelistFile%rushh1=rushh1
       oneNamelistFile%rushh2=rushh2
       oneNamelistFile%daylight=daylight
       oneNamelistFile%efsat=efsat
       oneNamelistFile%efsun=efsun
       oneNamelistFile%eindno=eindno
       oneNamelistFile%eindno2=eindno2
       oneNamelistFile%eindpm=eindpm
       oneNamelistFile%eindco=eindco
       oneNamelistFile%eindso2=eindso2
       oneNamelistFile%eindvoc=eindvoc
       oneNamelistFile%eveino=eveino
       oneNamelistFile%eveino2=eveino2
       oneNamelistFile%eveipm=eveipm
       oneNamelistFile%eveico=eveico
       oneNamelistFile%eveiso2=eveiso2
       oneNamelistFile%eveivoc=eveivoc
       oneNamelistFile%iteb=iteb
       oneNamelistFile%tminbld=tminbld
       oneNamelistFile%nteb=nteb
       oneNamelistFile%hc_roof=hc_roof
       oneNamelistFile%tc_roof=tc_roof
       oneNamelistFile%d_roof=d_roof
       oneNamelistFile%hc_road=hc_road
       oneNamelistFile%tc_road=tc_road
       oneNamelistFile%d_road=d_road
       oneNamelistFile%hc_wall=hc_wall
       oneNamelistFile%tc_wall=tc_wall
       oneNamelistFile%d_wall=d_wall
       oneNamelistFile%nurbtype=nurbtype
       oneNamelistFile%ileafcod=ileafcod
       oneNamelistFile%z0_town=z0_town
       oneNamelistFile%bld=bld
       oneNamelistFile%bld_height=bld_height
       oneNamelistFile%bld_hl_ratio=bld_hl_ratio
       oneNamelistFile%aroof=aroof
       oneNamelistFile%eroof=eroof
       oneNamelistFile%aroad=aroad
       oneNamelistFile%eroad=eroad
       oneNamelistFile%awall=awall
       oneNamelistFile%ewall=ewall
       oneNamelistFile%htraf=htraf
       oneNamelistFile%hindu=hindu
       oneNamelistFile%pletraf=pletraf
       oneNamelistFile%pleindu=pleindu
    end if

    read (iunit, iostat=err, NML=MODEL_FILE_INFO)

    if (err /= 0 .and. err /= 106) then
       write(*,"(a)") h//"**(ERROR)** reading section MODEL_FILE_INFO "//&
            &"of namelist file "//trim(oneNamelistFile%fileName)
       write(*,"(a)") h//" compare values read with file contents:"
       write (*, "(a)") " namelist MODEL_FILE_INFO: "
       write (*,*) "initial=", initial
       write (*,*) "varfpfx=", trim(varfpfx)
       write (*,*) "nudlat=", nudlat
       write (*,*) "tnudlat=", tnudlat
       write (*,*) "tnudcent=", tnudcent
       write (*,*) "tnudtop=", tnudtop
       write (*,*) "znudtop=", znudtop
       write (*,*) "ioutput=", ioutput
       write (*,*) "hfilout=", trim(hfilout)
       write (*,*) "afilout=", trim(afilout)
       write (*,*) "frqhis=", frqhis
       write (*,*) "frqanl=", frqanl
       write (*,*) "ipos=", ipos
       write (*,*) "topfiles=", trim(topfiles)
       write (*,*) "sfcfiles=", trim(sfcfiles)
       write (*,*) "sstfpfx=", trim(sstfpfx)
       write (*,*) "ndvifpfx=", trim(ndvifpfx)
       write (*,*) "itoptfn=", (trim(itoptfn(i))//";", i =1,size(itoptfn))
       write (*,*) "isstfn=", (trim(isstfn(i))//";", i=1,size(isstfn))
       write (*,*) "ivegtfn=", (trim(ivegtfn(i))//";", i = 1, size(ivegtfn))
       write (*,*) "isoilfn=", (trim(isoilfn(i))//";", i = 1, size(isoilfn))
       write (*,*) "ndvifn=", (trim(ndvifn(i))//";", i=1,size(ndvifn))
       call fatal_error(h//" reading namelist")
    else
       oneNamelistFile%initial=initial
       oneNamelistFile%varfpfx=varfpfx
       oneNamelistFile%nudlat=nudlat
       oneNamelistFile%tnudlat=tnudlat
       oneNamelistFile%tnudcent=tnudcent
       oneNamelistFile%tnudtop=tnudtop
       oneNamelistFile%znudtop=znudtop
       oneNamelistFile%ioutput=ioutput
       oneNamelistFile%hfilout=hfilout
       oneNamelistFile%afilout=afilout
       oneNamelistFile%frqhis=frqhis
       oneNamelistFile%frqanl=frqanl
       oneNamelistFile%ipos=ipos
       oneNamelistFile%topfiles=topfiles
       oneNamelistFile%sfcfiles=sfcfiles
       oneNamelistFile%sstfpfx=sstfpfx
       oneNamelistFile%ndvifpfx=ndvifpfx
       oneNamelistFile%iupdndvi=iupdndvi
       oneNamelistFile%iupdsst=iupdsst
       oneNamelistFile%itoptfn=itoptfn
       oneNamelistFile%isstfn=isstfn
       oneNamelistFile%ivegtfn=ivegtfn
       oneNamelistFile%isoilfn=isoilfn
       oneNamelistFile%ndvifn=ndvifn
    end if

    read (iunit2, iostat=err, NML=MODEL_FILE_INFO2)

    if (err /= 0 .and. err /= 106) then
       write(*,"(a)") h//"**(ERROR)** reading section MODEL_FILE_INFO "//&
            &"of namelist file "//trim(advanced_ramsin)
       write(*,"(a)") h//" compare values read with file contents:"
       write (*, "(a)") " namelist MODEL_FILE_INFO: "
       write (*,*) "nud_type=", nud_type
       write (*,*) "vwait1=", vwait1
       write (*,*) "vwaittot=", vwaittot
       write (*,*) "nud_hfile=", trim(nud_hfile)
       write (*,*) "nudlat=", nudlat
       write (*,*) "timeWindowIAU=", timeWindowIAU
       write (*,*) "ramp=", ramp
       write (*,*) "wt_nudge_grid=", wt_nudge_grid
       write (*,*) "wt_nudge_uv=", wt_nudge_uv
       write (*,*) "wt_nudge_th=", wt_nudge_th
       write (*,*) "wt_nudge_pi=", wt_nudge_pi
       write (*,*) "wt_nudge_rt=", wt_nudge_rt
       write (*,*) "applyIAU=", applyIAU
       write (*,*) "fileNameIAU=", trim(fileNameIAU)
       write (*,*) "nud_cond=", nud_cond
       write (*,*) "cond_hfile=", trim(cond_hfile)
       write (*,*) "tcond_beg=", tcond_beg
       write (*,*) "tcond_end=", tcond_end
       write (*,*) "t_nudge_rc=", t_nudge_rc
       write (*,*) "wt_nudgec_grid=", wt_nudgec_grid
       write (*,*) "if_oda=", if_oda
       write (*,*) "oda_upaprefix=", trim(oda_upaprefix)
       write (*,*) "oda_sfcprefix=", trim(oda_sfcprefix)
       write (*,*) "frqoda=", frqoda
       write (*,*) "todabeg=", todabeg
       write (*,*) "todaend=", todaend
       write (*,*) "tnudoda=", tnudoda
       write (*,*) "wt_oda_grid=", wt_oda_grid
       write (*,*) "wt_oda_uv=", wt_oda_uv
       write (*,*) "wt_oda_th=", wt_oda_th
       write (*,*) "wt_oda_pi=", wt_oda_pi
       write (*,*) "wt_oda_rt=", wt_oda_rt
       write (*,*) "roda_sfce=", roda_sfce
       write (*,*) "roda_sfc0=", roda_sfc0
       write (*,*) "roda_upae=", roda_upae
       write (*,*) "roda_upa0=", roda_upa0
       write (*,*) "roda_hgt=", roda_hgt
       write (*,*) "roda_zfact=", roda_zfact
       write (*,*) "oda_sfc_til=", oda_sfc_til
       write (*,*) "oda_sfc_tel=", oda_sfc_tel
       write (*,*) "oda_upa_til=", oda_upa_til
       write (*,*) "oda_upa_tel=", oda_upa_tel
       write (*,*) "if_cuinv=", if_cuinv
       write (*,*) "cu_prefix=", trim(cu_prefix)
       write (*,*) "tnudcu=", tnudcu
       write (*,*) "wt_cu_grid=", wt_cu_grid
       write (*,*) "tcu_beg=", tcu_beg
       write (*,*) "tcu_end=", tcu_end
       write (*,*) "cu_tel=", cu_tel
       write (*,*) "cu_til=", cu_til
       write (*,*) "timstr=", timstr
       write (*,*) "hfilin=", trim(hfilin)
       write (*,*) "ipastin=", ipastin
       write (*,*) "pastfn=", trim(pastfn)
       write (*,*) "iclobber=", iclobber
       write (*,*) "ihistdel=", ihistdel
       write (*,*) "frqlite=", frqlite
       write (*,*) "xlite=", xlite
       write (*,*) "ylite=", ylite
       write (*,*) "zlite=", zlite
       write (*,*) "nlite_vars=", nlite_vars
       write (*,*) "lite_vars=", (trim(lite_vars(i))//";", i=1,size(lite_vars))
       write (*,*) "avgtim=", avgtim
       write (*,*) "frqmean=", frqmean
       write (*,*) "frqboth=", frqboth
       write (*,*) "kwrite=", kwrite
       write (*,*) "frqprt=", frqprt
       write (*,*) "initfld=", initfld
       write (*,*) "prtcputime", prtcputime
       write (*,*) "itoptflg=", itoptflg
       write (*,*) "isstflg=", isstflg
       write (*,*) "ivegtflg=", ivegtflg
       write (*,*) "isoilflg=", isoilflg
       write (*,*) "ndviflg=", ndviflg
       write (*,*) "nofilflg=", nofilflg
       write (*,*) "iupdndvi=", iupdndvi
       write (*,*) "iupdsst=", iupdsst
       write (*,*) "itopsflg=", itopsflg
       write (*,*) "toptenh=", toptenh
       write (*,*) "toptwvl=", toptwvl
       write (*,*) "iz0flg=", iz0flg
       write (*,*) "z0max=", z0max
       write (*,*) "z0fact=", z0fact
       write (*,*) "mkcoltab=", mkcoltab
       write (*,*) "coltabfn=", trim(coltabfn)
       write (*,*) "mapaotfile=", trim(mapaotfile)
       write (*,*) "JulesIn=",trim(julesin)
       call fatal_error(h//" reading namelist")
    else
       oneNamelistFile%nud_type=nud_type
       oneNamelistFile%vwait1=vwait1
       oneNamelistFile%vwaittot=vwaittot
       oneNamelistFile%nud_hfile=nud_hfile
       oneNamelistFile%nudlat=nudlat
       oneNamelistFile%timeWindowIAU=timeWindowIAU
       oneNamelistFile%ramp=ramp
       oneNamelistFile%wt_nudge_grid=wt_nudge_grid
       oneNamelistFile%wt_nudge_uv=wt_nudge_uv
       oneNamelistFile%wt_nudge_th=wt_nudge_th
       oneNamelistFile%wt_nudge_pi=wt_nudge_pi
       oneNamelistFile%wt_nudge_rt=wt_nudge_rt
       
       oneNamelistFile%applyIAU=applyIAU
       oneNamelistFile%fileNameIAU=fileNameIAU       
       
       oneNamelistFile%nud_cond=nud_cond
       oneNamelistFile%cond_hfile=cond_hfile
       oneNamelistFile%tcond_beg=tcond_beg
       oneNamelistFile%tcond_end=tcond_end
       oneNamelistFile%t_nudge_rc=t_nudge_rc
       oneNamelistFile%wt_nudgec_grid=wt_nudgec_grid
       oneNamelistFile%if_oda=if_oda
       oneNamelistFile%oda_upaprefix=oda_upaprefix
       oneNamelistFile%oda_sfcprefix=oda_sfcprefix
       oneNamelistFile%frqoda=frqoda
       oneNamelistFile%todabeg=todabeg
       oneNamelistFile%todaend=todaend
       oneNamelistFile%tnudoda=tnudoda
       oneNamelistFile%wt_oda_grid=wt_oda_grid
       oneNamelistFile%wt_oda_uv=wt_oda_uv
       oneNamelistFile%wt_oda_th=wt_oda_th
       oneNamelistFile%wt_oda_pi=wt_oda_pi
       oneNamelistFile%wt_oda_rt=wt_oda_rt
       oneNamelistFile%roda_sfce=roda_sfce
       oneNamelistFile%roda_sfc0=roda_sfc0
       oneNamelistFile%roda_upae=roda_upae
       oneNamelistFile%roda_upa0=roda_upa0
       oneNamelistFile%roda_hgt=roda_hgt
       oneNamelistFile%roda_zfact=roda_zfact
       oneNamelistFile%oda_sfc_til=oda_sfc_til
       oneNamelistFile%oda_sfc_tel=oda_sfc_tel
       oneNamelistFile%oda_upa_til=oda_upa_til
       oneNamelistFile%oda_upa_tel=oda_upa_tel
       oneNamelistFile%if_cuinv=if_cuinv
       oneNamelistFile%cu_prefix=cu_prefix
       oneNamelistFile%tnudcu=tnudcu
       oneNamelistFile%wt_cu_grid=wt_cu_grid
       oneNamelistFile%tcu_beg=tcu_beg
       oneNamelistFile%tcu_end=tcu_end
       oneNamelistFile%cu_tel=cu_tel
       oneNamelistFile%cu_til=cu_til
       oneNamelistFile%timstr=timstr
       oneNamelistFile%hfilin=hfilin
       oneNamelistFile%ipastin=ipastin
       oneNamelistFile%pastfn=pastfn
       oneNamelistFile%iclobber=iclobber
       oneNamelistFile%ihistdel=ihistdel
       oneNamelistFile%frqlite=frqlite
       oneNamelistFile%xlite=xlite
       oneNamelistFile%ylite=ylite
       oneNamelistFile%zlite=zlite
       oneNamelistFile%nlite_vars=nlite_vars
       oneNamelistFile%lite_vars=lite_vars
       oneNamelistFile%avgtim=avgtim
       oneNamelistFile%frqmean=frqmean
       oneNamelistFile%frqboth=frqboth
       oneNamelistFile%kwrite=kwrite
       oneNamelistFile%frqprt=frqprt
       oneNamelistFile%initfld=initfld
       oneNamelistFile%prtcputime=prtcputime
       oneNamelistFile%itoptflg=itoptflg
       oneNamelistFile%isstflg=isstflg
       oneNamelistFile%ivegtflg=ivegtflg
       oneNamelistFile%isoilflg=isoilflg
       oneNamelistFile%ndviflg=ndviflg
       oneNamelistFile%nofilflg=nofilflg
       oneNamelistFile%iupdndvi=iupdndvi
       oneNamelistFile%iupdsst=iupdsst
       oneNamelistFile%itopsflg=itopsflg
       oneNamelistFile%toptenh=toptenh
       oneNamelistFile%toptwvl=toptwvl
       oneNamelistFile%iz0flg=iz0flg
       oneNamelistFile%z0max=z0max
       oneNamelistFile%z0fact=z0fact
       oneNamelistFile%mkcoltab=mkcoltab
       oneNamelistFile%coltabfn=coltabfn
       oneNamelistFile%mapaotfile=mapaotfile
       oneNamelistFile%julesin=julesin
    end if

    read (iunit, iostat=err, NML=MODEL_OPTIONS)
    if (err /= 0 .and. err /= 106) then
       write(*,"(a)") h//"**(ERROR)** reading section MODEL_OPTIONS "//&
            &"of namelist file "//trim(oneNamelistFile%fileName)
       write(*,"(a)") h//" compare values read with file contents:"
       write (*, *) "sfire=",sfire
       write (*, *) "iswrtyp=",iswrtyp
       write (*, *) "ilwrtyp=",ilwrtyp
       write (*, *) "radfrq=",radfrq
       write (*, *) "nnqparm=",nnqparm
       write (*, *) "closure_type=",closure_type
       write (*, *) "nnshcu=",nnshcu
       write (*, *) "confrq=",confrq
       write (*, *) "shcufrq=",shcufrq
       write (*, *) "isfcl=",isfcl
       write (*, *) "isfcl_ocean=",isfcl_ocean
       write (*, *) "soil_moist_fail=",soil_moist_fail
       write (*, *) "usdata_in=",trim(usdata_in)
       write (*, *) "usmodel_in=",trim(usmodel_in)
       write (*, *) "mcphys_type=",mcphys_type
       write (*, *) "level=",level
       call fatal_error(h//" reading namelist")
    else
 !--(DMK-CCATT-INI)-----------------------------------------------------------
       oneNamelistFile%sfire = sfire
       oneNamelistFile%iswrtyp=iswrtyp
       oneNamelistFile%ilwrtyp=ilwrtyp
       oneNamelistFile%radfrq=radfrq
       oneNamelistFile%nnqparm=nnqparm
       oneNamelistFile%closure_type=closure_type
       oneNamelistFile%nnshcu=nnshcu
       oneNamelistFile%confrq=confrq
       oneNamelistFile%shcufrq=shcufrq
       oneNamelistFile%isfcl=isfcl
       oneNamelistFile%isfcl_ocean=isfcl_ocean
       oneNamelistFile%soil_moist_fail=soil_moist_fail
       oneNamelistFile%usdata_in=usdata_in
       oneNamelistFile%usmodel_in=usmodel_in
       oneNamelistFile%mcphys_type=mcphys_type
       oneNamelistFile%level=level
    end if


    read (iunit2, iostat=err, NML=MODEL_OPTIONS2)
    if (err /= 0 .and. err /= 106) then
       write(*,"(a)") h//"**(ERROR)** reading section MODEL_OPTIONS2 "//&
            &"of namelist file "//trim(advanced_ramsin)
       write(*,"(a)") h//" compare values read with file contents:"
       write (*,*) "naddsc=",naddsc
       write (*, *) "icorflg=",icorflg
       write (*, *) "dyncore_flag=",dyncore_flag
       write (*, *) "pd_or_mnt_constraint=", pd_or_mnt_constraint
       write (*, *) "order_h=", order_h
       write (*, *) "order_v=",order_v
       write (*, *) "iexev=",iexev
       write (*, *) "imassflx=",imassflx
       write (*, *) "vveldamp=",vveldamp
       write (*, *) "ibnd=",ibnd
       write (*, *) "jbnd=",jbnd
       write (*, *) "cphas=",cphas
       write (*, *) "lsflg=",lsflg
       write (*, *) "nfpt=",nfpt
       write (*, *) "distim=",distim
       write (*, *) "raddatfn=", RADDATFN
       write (*, *) "radtun=",radtun
       write (*, *) "lonrad=",lonrad
       write (*, *) "wcldbs=",wcldbs
       write (*, *) "g3d_spread=",g3d_spread
       write (*, *) "g3d_smoothh=",g3d_smoothh
       write (*, *) "g3d_smoothv=",g3d_smoothv
       write (*, *) "npatch=",npatch
       write (*, *) "nvegpat=",nvegpat
       write (*, *) "nvgcon=",nvgcon
       write (*, *) "pctlcon=",pctlcon
       write (*, *) "nslcon=",nslcon
       write (*, *) "drtcon=",drtcon
       write (*, *) "zrough=",zrough
       write (*, *) "albedo=",albedo
       write (*, *) "seatmp=",seatmp
       write (*, *) "dthcon=",dthcon
       write (*, *) "soil_moist=",soil_moist
       write (*, *) "slz=",slz
       write (*, *) "slmstr=",slmstr
       write (*, *) "stgoff=",stgoff
       write (*, *) "if_urban_canopy=",if_urban_canopy
       write (*, *) "idiffk=",idiffk
       write (*, *) "ihorgrad=",ihorgrad
       write (*, *) "csx=",csx
       write (*, *) "csz=",csz
       write (*, *) "xkhkm=",xkhkm
       write (*, *) "zkhkm=",zkhkm
       write (*, *) "akmin=",akmin
       write (*, *) "mcphys_type=",mcphys_type
       write (*, *) "level=",level
       write (*, *) "icloud=",icloud
       write (*, *) "idriz=",idriz
       write (*, *) "irime=",irime
       write (*, *) "iplaws=",iplaws
       write (*, *) "iccnlev=",iccnlev
       write (*, *) "irain=",irain
       write (*, *) "ipris=",ipris
       write (*, *) "isnow=",isnow
       write (*, *) "iaggr=",iaggr
       write (*, *) "igraup=",igraup
       write (*, *) "ihail=",ihail
       write (*, *) "cparm=",cparm
       write (*, *) "rparm=",rparm
       write (*, *) "pparm=",pparm
       write (*, *) "sparm=",sparm
       write (*, *) "aparm=",aparm
       write (*, *) "gparm=",gparm
       write (*, *) "hparm=",hparm
       write (*, *) "dparm=",dparm
       write (*, *) "cnparm=",cnparm
       write (*, *) "gnparm=",gnparm
       write (*, *) "epsil=",epsil
       write (*, *) "gnu=",gnu
       write (*, *) "windfarm=",windfarm
       write (*, *) "wfFile=",wfFile
       write (*, *) "damModule=",damModule
       write (*, *) "frqPrecip=",frqPrecip
       write (*, *) "damOutPrefix=",damOutPrefix
       write (*, *) "evaluate=",evaluate
       write (*, *) "evaluatePrefix=",evaluatePrefix
       write (*,*)  "advmnt=", advmnt
       write (*,*)  "GhostZoneLength", GhostZoneLength
       call fatal_error(h//" reading namelist")
    else

       oneNamelistFile%dyncore_flag        =dyncore_flag
       oneNamelistFile%pd_or_mnt_constraint=pd_or_mnt_constraint
       oneNamelistFile%order_h             =order_h
       oneNamelistFile%order_v             =order_v
       oneNamelistFile%advmnt              =advmnt
       oneNamelistFile%GhostZoneLength     =GhostZoneLength
       oneNamelistFile%naddsc              =naddsc
       oneNamelistFile%icorflg             =icorflg
       oneNamelistFile%iexev               =iexev
       oneNamelistFile%imassflx            =imassflx
       oneNamelistFile%vveldamp            =vveldamp
       oneNamelistFile%ibnd                =ibnd
       oneNamelistFile%jbnd                =jbnd
       oneNamelistFile%cphas               =cphas
       oneNamelistFile%lsflg               =lsflg
       oneNamelistFile%nfpt                =nfpt
       oneNamelistFile%distim              =distim
       oneNamelistFile%raddatfn            =raddatfn
       oneNamelistFile%radtun              =radtun
       oneNamelistFile%lonrad              =lonrad
       oneNamelistFile%wcldbs=wcldbs
       oneNamelistFile%g3d_spread=g3d_spread
       oneNamelistFile%g3d_smoothh=g3d_smoothh
       oneNamelistFile%g3d_smoothv=g3d_smoothv
       oneNamelistFile%npatch=npatch
       oneNamelistFile%nvegpat=nvegpat
       oneNamelistFile%nvgcon=nvgcon
       oneNamelistFile%pctlcon=pctlcon
       oneNamelistFile%nslcon=nslcon
       oneNamelistFile%drtcon=drtcon
       oneNamelistFile%zrough=zrough
       oneNamelistFile%albedo=albedo
       oneNamelistFile%seatmp=seatmp
       oneNamelistFile%dthcon=dthcon
       oneNamelistFile%soil_moist=soil_moist
       oneNamelistFile%slz=slz
       oneNamelistFile%slmstr=slmstr
       oneNamelistFile%stgoff=stgoff
       oneNamelistFile%if_urban_canopy=if_urban_canopy
       oneNamelistFile%idiffk=idiffk
       oneNamelistFile%ihorgrad=ihorgrad
       oneNamelistFile%csx=csx
       oneNamelistFile%csz=csz
       oneNamelistFile%xkhkm=xkhkm
       oneNamelistFile%zkhkm=zkhkm
       oneNamelistFile%akmin=akmin
       oneNamelistFile%icloud=icloud
       oneNamelistFile%idriz=idriz
       oneNamelistFile%iccnlev=iccnlev
       oneNamelistFile%iplaws=iplaws
       oneNamelistFile%irime=irime
       oneNamelistFile%irain=irain
       oneNamelistFile%ipris=ipris
       oneNamelistFile%isnow=isnow
       oneNamelistFile%iaggr=iaggr
       oneNamelistFile%igraup=igraup
       oneNamelistFile%ihail=ihail
       oneNamelistFile%cparm=cparm
       oneNamelistFile%rparm=rparm
       oneNamelistFile%pparm=pparm
       oneNamelistFile%sparm=sparm
       oneNamelistFile%aparm=aparm
       oneNamelistFile%gparm=gparm
       oneNamelistFile%hparm=hparm
       oneNamelistFile%dparm=dparm
       oneNamelistFile%cnparm=cnparm
       oneNamelistFile%gnparm=gnparm
       oneNamelistFile%epsil=epsil
       oneNamelistFile%gnu=gnu
       oneNamelistFile%windfarm=windfarm
       oneNamelistFile%wfFile=wfFile
       oneNamelistFile%damModule=damModule
       oneNamelistFile%frqPrecip=frqPrecip
       oneNamelistFile%damOutPrefix=damOutPrefix
       oneNamelistFile%evaluate=evaluate
       oneNamelistFile%evaluatePrefix=evaluatePrefix
    end if


    read (iunit2, iostat=err, NML=MODEL_SOUND)
    if (err /= 0 .and. err /= 106 .and. err/=111) then
       write(*,"(a)") h//"**(ERROR)** reading section MODEL_SOUND "//&
            &"of namelist file "//trim(advanced_ramsin)
       write(*,"(a)") h//" compare values read with file contents:"
       write (*, *) "ipsflg=",ipsflg
       write (*, *) "itsflg=",itsflg
       write (*, *) "irtsflg=",irtsflg
       write (*, *) "iusflg=",iusflg
       write (*, *) "hs=",hs
       write (*, *) "ps=",ps
       write (*, *) "ts=",ts
       write (*, *) "rts=",rts
       write (*, *) "us=",us
       write (*, *) "vs=",vs
       call fatal_error(h//" reading namelist")
    else
       oneNamelistFile%ipsflg=ipsflg
       oneNamelistFile%itsflg=itsflg
       oneNamelistFile%irtsflg=irtsflg
       oneNamelistFile%iusflg=iusflg
       oneNamelistFile%hs=hs
       oneNamelistFile%ps=ps
       oneNamelistFile%ts=ts
       oneNamelistFile%rts=rts
       oneNamelistFile%us=us
       oneNamelistFile%vs=vs
    end if

    read (iunit2, iostat=err, NML=MODEL_PRINT)
    if (err /= 0 .and. err /= 106) then
       write(*,"(a)") h//"**(ERROR)** reading section MODEL_PRINT "//&
            &"of namelist file "//trim(advanced_ramsin)
       write(*,"(a)") h//" compare values read with file contents:"
       write (*, *) "nplt=",nplt
       write (*, *) "iplfld=",(trim(iplfld(i))//";", i=1,size(iplfld))
       write (*, *) "ixsctn=",ixsctn
       write (*, *) "isbval=",isbval
       call fatal_error(h//" reading namelist")
    else
       oneNamelistFile%nplt=nplt
       oneNamelistFile%iplfld=iplfld
       oneNamelistFile%ixsctn=ixsctn
       oneNamelistFile%isbval=isbval
    end if

    read (iunit, iostat=err, NML=ISAN_CONTROL)
    if (err /= 0 .and. err /= 106) then
       write(*,"(a)") h//"**(ERROR)** reading section ISAN_CONTROL "//&
            &"of namelist file "//trim(oneNamelistFile%fileName)
       write(*,"(a)") h//" compare values read with file contents:"
       write (*, *) "isan_inc=",isan_inc
       write (*, *) "iapr=",trim(iapr)
       write (*, *) "varpfx=",trim(varpfx)
       call fatal_error(h//" reading namelist")
    else
       oneNamelistFile%isan_inc=isan_inc
       oneNamelistFile%iapr=iapr
       oneNamelistFile%varpfx=varpfx
    end if

    read (iunit2, iostat=err, NML=ISAN_CONTROL2)
    if (err /= 0 .and. err /= 106) then
       write(*,"(a)") h//"**(ERROR)** reading section ISAN_CONTROL2 "//&
            &"of namelist file "//trim(advanced_ramsin)
       write(*,"(a)") h//" compare values read with file contents:"
       write (*, *) "iszstage=",iszstage
       write (*, *) "ivrstage=",ivrstage
       write (*, *) "guess1st=",guess1st
       write (*, *) "i1st_flg=",i1st_flg
       write (*, *) "iupa_flg=",iupa_flg
       write (*, *) "isfc_flg=",isfc_flg
       write (*, *) "iarawi=",trim(iarawi)
       write (*, *) "iasrfce=",trim(iasrfce)
       write (*, *) "ioflgisz=",ioflgisz
       write (*, *) "ioflgvar=",ioflgvar
       call fatal_error(h//" reading namelist")
    else
       oneNamelistFile%iszstage=iszstage
       oneNamelistFile%ivrstage=ivrstage
       oneNamelistFile%guess1st=guess1st
       oneNamelistFile%i1st_flg=i1st_flg
       oneNamelistFile%iupa_flg=iupa_flg
       oneNamelistFile%isfc_flg=isfc_flg
       oneNamelistFile%iarawi=iarawi
       oneNamelistFile%iasrfce=iasrfce
       oneNamelistFile%ioflgisz=ioflgisz
       oneNamelistFile%ioflgvar=ioflgvar
    end if

    read (iunit, iostat=err, NML=ISAN_ISENTROPIC)
    if (err /= 0 .and. err /= 106) then
       write(*,"(a)") h//"**(ERROR)** reading section ISAN_ISENTROPIC "//&
            &"of namelist file "//trim(oneNamelistFile%fileName)
       write(*,"(a)") h//" compare values read with file contents:"
       write (*,*)  "icFileType=", icFileType             
       write (*,*)  "icPrefix=",trim(icPrefix)
       write (*,*)  "wind_u_varname=",trim(wind_u_varname)
       write (*,*)  "wind_v_varname=",trim(wind_v_varname)
       write (*,*)  "temperature_varname=",trim(temperature_varname)
       write (*,*)  "geo_varname=",trim(geo_varname)
       write (*,*)  "ur_varname=",trim(ur_varname)
       write (*,*)  "initial_latitude=",initial_latitude
       write (*,*)  "final_latitude=",final_latitude
       write (*,*)  "initial_longitude=",initial_longitude
       write (*,*)  "final_longitude=",final_longitude
       write (*,*)  "z_max_level=",z_max_level
       write (*,*)  "scale_factor=",scale_factor


       call fatal_error(h//" reading namelist")
    else
       oneNamelistFile%icFileType= icFileType             
       oneNamelistFile%icPrefix=trim(icPrefix)
       oneNamelistFile%wind_u_varname=trim(wind_u_varname)
       oneNamelistFile%wind_v_varname=trim(wind_v_varname)
       oneNamelistFile%temperature_varname=trim(temperature_varname)
       oneNamelistFile%geo_varname=trim(geo_varname)
       oneNamelistFile%ur_varname=trim(ur_varname)
       oneNamelistFile%initial_latitude=initial_latitude
       oneNamelistFile%final_latitude=final_latitude
       oneNamelistFile%initial_longitude=initial_longitude
       oneNamelistFile%final_longitude=final_longitude
       oneNamelistFile%z_max_level=z_max_level
       oneNamelistFile%scale_factor=scale_factor
    end if



    read (iunit2, iostat=err, NML=ISAN_ISENTROPIC2)
    if (err /= 0 .and. err /= 106) then
       write(*,"(a)") h//"**(ERROR)** reading section ISAN_ISENTROPIC2 "//&
            &"of namelist file "//trim(advanced_ramsin)
       write(*,"(a)") h//" compare values read with file contents:"
       write (*, *) "nisn=",nisn
       write (*, *) "levth=",levth
       write (*, *) "nigrids=",nigrids
       write (*, *) "topsigz=",topsigz
       write (*, *) "hybbot=",hybbot
       write (*, *) "hybtop=",hybtop
       write (*, *) "sfcinf=",sfcinf
       write (*, *) "sigzwt=",sigzwt
       write (*, *) "nfeedvar=",nfeedvar
       write (*, *) "maxsta=",maxsta
       write (*, *) "maxsfc=",maxsfc
       write (*, *) "notsta=",notsta
       write (*, *) "notid=",(trim(notid(i))//";", i=1,size(notid))
       write (*, *) "iobswin=",iobswin
       write (*, *) "stasep=",stasep
       write (*, *) "igridfl=",igridfl
       write (*, *) "gridwt=",gridwt
       write (*, *) "gobsep=",gobsep
       write (*, *) "gobrad=",gobrad
       write (*, *) "wvlnth=",wvlnth
       write (*, *) "swvlnth=",swvlnth
       write (*, *) "respon=",respon
       write (*,*)  "dlimit=",dlimit
       write (*,*)  "ulimit=",ulimit
       write (*,*)  "ccGradsWrite=",ccGradsWrite
       write (*,*)  "icGradsPrefix=",trim(icGradsPrefix)


       call fatal_error(h//" reading namelist")
    else
       oneNamelistFile%nisn=nisn
       oneNamelistFile%levth=levth
       oneNamelistFile%nigrids=nigrids
       oneNamelistFile%topsigz=topsigz
       oneNamelistFile%hybbot=hybbot
       oneNamelistFile%hybtop=hybtop
       oneNamelistFile%sfcinf=sfcinf
       oneNamelistFile%sigzwt=sigzwt
       oneNamelistFile%nfeedvar=nfeedvar
       oneNamelistFile%maxsta=maxsta
       oneNamelistFile%maxsfc=maxsfc
       oneNamelistFile%notsta=notsta
       oneNamelistFile%notid=notid
       oneNamelistFile%iobswin=iobswin
       oneNamelistFile%stasep=stasep
       oneNamelistFile%igridfl=igridfl
       oneNamelistFile%gridwt=gridwt
       oneNamelistFile%gobsep=gobsep
       oneNamelistFile%gobrad=gobrad
       oneNamelistFile%wvlnth=wvlnth
       oneNamelistFile%swvlnth=swvlnth
       oneNamelistFile%respon=respon
       oneNamelistFile%dlimit=dlimit
       oneNamelistFile%ulimit=ulimit
       oneNamelistFile%ccGradsWrite=ccGradsWrite
       oneNamelistFile%icGradsPrefix=icGradsPrefix
    end if

    ! namelist POST
    read (iunit, iostat=err, NML=POST)
    if (err /= 0 .and. err /= 106) then
       write(*,"(a)") h//"**(ERROR)** reading section POST "//&
            &"of namelist file "//trim(oneNamelistFile%fileName)
       write (*, *) "nvp="         ,nvp
       write (*, *) "vp="          ,vp
       write (*, *) "gprefix="     ,gprefix
       write (*, *) "csvFile="     ,csvFile
       write (*, *) "anl2gra="     ,anl2gra
       write (*, *) "proj="        ,proj
       write (*, *) "mean_type="   ,mean_type
       write (*, *) "lati="        ,lati
       write (*, *) "loni="        ,loni
       write (*, *) "latf="        ,latf
       write (*, *) "lonf="        ,lonf
       write (*, *) "zlevmax="     ,zlevmax
       write (*, *) "ipresslev="   ,ipresslev
       write (*, *) "inplevs="     ,inplevs
       write (*, *) "iplevs="      ,iplevs
       write (*, *) "mechanism="   ,mechanism
       write (*, *) "ascii_data="  ,ascii_data
       write (*, *) "site_lat="    ,site_lat
       write (*, *) "site_lon="    ,site_lon
       call fatal_error(h//" reading namelist")
    else
       oneNamelistFile%nvp       =nvp
       oneNamelistFile%vp        =vp
       oneNamelistFile%gprefix   =gprefix
       oneNamelistFile%csvFile   =csvFile
       oneNamelistFile%anl2gra   =anl2gra
       oneNamelistFile%proj      =proj
       oneNamelistFile%mean_type =mean_type
       oneNamelistFile%lati      =lati
       oneNamelistFile%loni      =loni
       oneNamelistFile%latf      =latf
       oneNamelistFile%lonf      =lonf
       oneNamelistFile%zlevmax   =zlevmax
       oneNamelistFile%ipresslev =ipresslev
       oneNamelistFile%inplevs   =inplevs
       oneNamelistFile%iplevs    =iplevs
       oneNamelistFile%mechanism =mechanism
       oneNamelistFile%ascii_data=ascii_data
       oneNamelistFile%site_lat  =site_lat
       oneNamelistFile%site_lon  =site_lon
    end if

    ! namelist DIGITAL FILTER
    read (iunit2, iostat=err, NML=DIGITALFILTER)
    if (err /= 0 .and. err /= 106) then
       write(*,"(a)") h//"**(ERROR)** reading section DIGITAL FILTER "//&
            &"of namelist file "//trim(advanced_ramsin)

    else
	oneNamelistFile%applyDigitalFilter = applyDigitalFilter
      	oneNamelistFile%digitalFilterTimeWindow=digitalFilterTimeWindow
    end if

       ! namelist METEOGRAM
    read (iunit2, iostat=err, NML=METEOGRAM)
    if (err /= 0 .and. err /= 106) then
       write(*,"(a)") h//"**(ERROR)** reading section METEOGRAM "//&
            &"of namelist file "//trim(advanced_ramsin)

    else
	oneNamelistFile%applyDigitalFilter = applyDigitalFilter
      	oneNamelistFile%digitalFilterTimeWindow=digitalFilterTimeWindow

	oneNamelistFile%applyMeteogram = applyMeteogram
     	oneNamelistFile%meteogramFreq  = meteogramFreq
     	oneNamelistFile%meteogramMap   = meteogramMap
     	oneNamelistFile%meteogramDir   = meteogramDir

    end if

    close(iunit, iostat=err)
    if (err /= 0) then
       write(c0,"(i10)") err
       call fatal_error(h//" closing file "//&
            trim(oneNamelistFile%fileName)//" returned iostat="//&
            trim(adjustl(c0)))
    end if
    close(iunit2, iostat=err)
    if (err /= 0) then
       write(c0,"(i10)") err
       call fatal_error(h//" closing file "//&
            trim(advanced_ramsin)//" returned iostat="//&
            trim(adjustl(c0)))
    end if

  end subroutine ReadNamelistFile







!!$  subroutine StoreNamelistFile(oneNamelistFile)
!!$
!!$    use io_params, only: frqboth, &
!!$         afilout,                 &
!!$         avgtim,                  &
!!$         frqanl,                  &
!!$         frqhis,                  &
!!$         frqlite,                 &
!!$         frqmean,                 &
!!$         frqprt,                  &
!!$         hfilin,                  &
!!$         hfilout,                 &
!!$         iclobber,                &
!!$         ihistdel,                &
!!$         initfld,                 &
!!$         prtcputime,              &
!!$         ioutput,                 &
!!$         ipastin,                 &
!!$         iplfld,                  &
!!$         isbval,                  &
!!$         isoilflg,                &
!!$         isoilfn,                 &
!!$         isstflg,                 &
!!$         isstfn,                  &
!!$         itopsflg,                &
!!$         itoptflg,                &
!!$         itoptfn,                 &
!!$         iupdndvi,                &
!!$         iupdsst,                 &
!!$         ivegtflg,                &
!!$         ivegtfn,                 &
!!$         ixsctn,                  &
!!$         iz0flg,                  &
!!$         kwrite,                  &
!!$         lite_vars,               &
!!$         ndviflg,                 &
!!$         ndvifn,                  &
!!$         ndvifpfx,                &
!!$         nlite_vars,              &
!!$         nofilflg,                &
!!$         nplt,                    &
!!$         pastfn,                  &
!!$         sfcfiles,                &
!!$         sstfpfx,                 &
!!$         timstr,                  &
!!$         topfiles,                &
!!$         toptenh,                 &
!!$         toptwvl,                 &
!!$         xlite,                   &
!!$         ylite,                   &
!!$         z0fact,                  &
!!$         z0max,                   &
!!$         zlite,                   &
!!$                                ! TEB
!!$         ifusflg,                 &
!!$         ifusfn,                  &
!!$         fusfiles
!!$    use isan_coms, only: gobrad, &
!!$         gobsep, &
!!$         gridwt, &
!!$         guess1st, &
!!$         hybbot, &
!!$         hybtop, &
!!$         i1st_flg, &
!!$         iapr, &
!!$         iarawi, &
!!$         iasrfce, &
!!$         igridfl, &
!!$         iobswin, &
!!$         ioflgisz, &
!!$         ioflgvar, &
!!$         isan_inc, &
!!$         isfc_flg, &
!!$         iszstage, &
!!$         iupa_flg, &
!!$         ivrstage, &
!!$         levth, &
!!$         maxsfc, &
!!$         maxsta, &
!!$         nfeedvar, &
!!$         nigrids, &
!!$         nisn, &
!!$         notid, &
!!$         notsta, &
!!$         respon, &
!!$         sfcinf, &
!!$         sigzwt, &
!!$         stasep, &
!!$         swvlnth, &
!!$         topsigz, &
!!$         varpfx, &
!!$         wvlnth
!!$    use mem_cuparm, only: confrq, &
!!$         cu_prefix, &
!!$         cu_tel, &
!!$         cu_til, &
!!$         if_cuinv, &
!!$         nnqparm, &
!!$         tcu_beg, &
!!$         tcu_end, &
!!$         tnudcu, &
!!$         wcldbs, &
!!$         wt_cu_grid
!!$    use mem_globrad, only: raddatfn
!!$    use mem_grell_param, only: closure_type
!!$    use mem_grid, only: centlat, &
!!$         centlon, &
!!$         cphas, &
!!$         deltax, &
!!$         deltay, &
!!$         deltaz, &
!!$         distim, &
!!$         dtlong, &
!!$         dzmax, &
!!$         dzrat, &
!!$         expnme, &
!!$         gridu, &
!!$         gridv, &
!!$         ibnd, &
!!$         icorflg, &
!!$         idate1, &
!!$         ideltat, &
!!$         if_adap, &
!!$         ihtran, &
!!$         imonth1, &
!!$         initial, &
!!$         itime1, &
!!$         iyear1, &
!!$         jbnd, &
!!$         lsflg, &
!!$         nacoust, &
!!$         naddsc, &
!!$         nestz1, &
!!$         nestz2, &
!!$         nfpt, &
!!$         ngrids, &
!!$         ninest, &
!!$         njnest, &
!!$         nknest, &
!!$         nndtrat, &
!!$         nnstbot, &
!!$         nnsttop, &
!!$         nnxp, &
!!$         nnyp, &
!!$         nnzp, &
!!$         npatch, &
!!$         nstratx, &
!!$         nstraty, &
!!$         nstratz1, &
!!$         nstratz2, &
!!$         nxtnest, &
!!$         nzg, &
!!$         nzs, &
!!$         polelat, &
!!$         polelon, &
!!$         runtype, &
!!$         timeunit, &
!!$         timmax, &
!!$         zz
!!$    use mem_leaf, only: albedo, &
!!$         drtcon, &
!!$         dthcon, &
!!$         isfcl, &
!!$         nslcon, &
!!$         nvegpat, &
!!$         nvgcon, &
!!$         pctlcon, &
!!$         seatmp, &
!!$         slmstr, &
!!$         slz, &
!!$         stgoff, &
!!$         zrough
!!$    use mem_oda, only: frqoda, &
!!$         if_oda, &
!!$         oda_sfc_tel, &
!!$         oda_sfc_til, &
!!$         oda_sfcprefix, &
!!$         oda_upa_tel, &
!!$         oda_upa_til, &
!!$         oda_upaprefix, &
!!$         roda_hgt, &
!!$         roda_sfc0, &
!!$         roda_sfce, &
!!$         roda_upa0, &
!!$         roda_upae, &
!!$         roda_zfact, &
!!$         tnudoda, &
!!$         todabeg, &
!!$         todaend, &
!!$         wt_oda_grid, &
!!$         wt_oda_pi, &
!!$         wt_oda_rt, &
!!$         wt_oda_th, &
!!$         wt_oda_uv
!!$    use mem_radiate, only: ilwrtyp, &
!!$         iswrtyp, &
!!$         lonrad, &
!!$         radfrq
!!$    use soilMoisture, only: soil_moist, &
!!$         soil_moist_fail, &
!!$         usdata_in, &
!!$         usmodel_in
!!$    use mem_turb, only: akmin, &
!!$         csx, &
!!$         csz, &
!!$         idiffk, &
!!$         if_urban_canopy, &
!!$         ihorgrad, &
!!$         xkhkm, &
!!$         zkhkm
!!$    use mem_varinit, only: cond_hfile, &
!!$         nud_cond, &
!!$         nud_hfile, &
!!$         nud_type, &
!!$         nudlat, &
!!$         t_nudge_rc, &
!!$         tcond_beg, &
!!$         tcond_end, &
!!$         tnudcent, &
!!$         tnudlat, &
!!$         tnudtop, &
!!$         varfpfx, &
!!$         vwait1, &
!!$         vwaittot, &
!!$         wt_nudge_grid, &
!!$         wt_nudge_pi, &
!!$         wt_nudge_rt, &
!!$         wt_nudge_th, &
!!$         wt_nudge_uv, &
!!$         wt_nudgec_grid, &
!!$         znudtop
!!$    use micphys, only: &
!!$         aparm, &
!!$         coltabfn, &
!!$         cparm, &
!!$         gnu, &
!!$         gparm, &
!!$         hparm, &
!!$         iaggr, &
!!$         icloud, &
!!$         igraup, &
!!$         ihail, &
!!$         ipris, &
!!$         irain, &
!!$         isnow, &
!!$         level, &
!!$         mkcoltab, &
!!$         pparm, &
!!$         rparm, &
!!$         sparm
!!$    use node_mod, only: &
!!$         load_bal
!!$    use ref_sounding, only: &
!!$         hs, &
!!$         ipsflg, &
!!$         irtsflg, &
!!$         itsflg, &
!!$         iusflg, &
!!$         ps, &
!!$         rts, &
!!$         ts, &
!!$         us, &
!!$         vs
!!$    use shcu_vars_const, only: &
!!$         nnshcu, &
!!$         shcufrq
!!$    use sib_vars, only: &
!!$         co2_init, &
!!$         n_co2
!!$
!!$    use catt_start, only: &
!!$         CATT
!!$
!!$    use emission_source_map, only: &
!!$         firemapfn, &
!!$         plumerise,                           &
!!$         define_proc
!!$
!!$    use plume_utils, only: &
!!$         prfrq
!!$
!!$    use mem_scalar, only: &
!!$         recycle_tracers
!!$
!!$    use teb_spm_start, only: &
!!$         teb_spm
!!$
!!$    use mem_emiss, only : &
!!$         ichemi,          &
!!$         ichemi_in,       &
!!$         chemdata_in,     &
!!$         isource,         &
!!$         weekdayin,       &
!!$         efsat,           &
!!$         efsun,           &
!!$         eindno,          &
!!$         eindno2,         &
!!$         eindpm,          &
!!$         eindco,          &
!!$         eindso2,         &
!!$         eindvoc,         &
!!$         eveino,          &
!!$         eveino2,         &
!!$         eveipm,          &
!!$         eveico,          &
!!$         eveiso2,         &
!!$         eveivoc
!!$
!!$    use teb_vars_const, only : &
!!$         rushh1,               &
!!$         rushh2,               &
!!$         daylight,             &
!!$         iteb,                 &
!!$         tminbld,              &
!!$         nteb,                 &
!!$         hc_roof,              &
!!$         tc_roof,              &
!!$         d_roof,               &
!!$         hc_road,              &
!!$         d_road,               &
!!$         tc_road,              &
!!$         d_wall,               &
!!$         tc_wall,              &
!!$         hc_wall,              &
!!$         nurbtype,             &
!!$         ileafcod,             &
!!$         z0_town,              &
!!$         bld,                  &
!!$         bld_height,           &
!!$         bld_hl_ratio,         &
!!$         aroof,                &
!!$         eroof,                &
!!$         aroad,                &
!!$         eroad,                &
!!$         awall,                &
!!$         ewall,                &
!!$         htraf,                &
!!$         hindu,                &
!!$         pletraf,              &
!!$         pleindu
!!$
!!$    ! Explicit domain decomposition
!!$    use domain_decomp, only: &
!!$         domain_fname
!!$
!!$    implicit none
!!$    type(namelistFile), pointer :: oneNamelistFile
!!$
!!$
!!$
!!$    frqboth = oneNamelistFile%frqboth
!!$    afilout = oneNamelistFile%afilout
!!$    avgtim = oneNamelistFile%avgtim
!!$    frqanl = oneNamelistFile%frqanl
!!$    frqhis = oneNamelistFile%frqhis
!!$    frqlite = oneNamelistFile%frqlite
!!$    frqmean = oneNamelistFile%frqmean
!!$    frqprt = oneNamelistFile%frqprt
!!$    hfilin = oneNamelistFile%hfilin
!!$    hfilout = oneNamelistFile%hfilout
!!$    iclobber = oneNamelistFile%iclobber
!!$    ihistdel = oneNamelistFile%ihistdel
!!$    initfld = oneNamelistFile%initfld
!!$    prtcputime = oneNamelistFile%prtcputime
!!$    ioutput = oneNamelistFile%ioutput
!!$    ipastin = oneNamelistFile%ipastin
!!$    iplfld = oneNamelistFile%iplfld
!!$    isbval = oneNamelistFile%isbval
!!$    isoilflg = oneNamelistFile%isoilflg
!!$    isoilfn = oneNamelistFile%isoilfn
!!$    isstflg = oneNamelistFile%isstflg
!!$    isstfn = oneNamelistFile%isstfn
!!$    itopsflg = oneNamelistFile%itopsflg
!!$    itoptflg = oneNamelistFile%itoptflg
!!$    itoptfn = oneNamelistFile%itoptfn
!!$    iupdndvi = oneNamelistFile%iupdndvi
!!$    iupdsst = oneNamelistFile%iupdsst
!!$    ivegtflg = oneNamelistFile%ivegtflg
!!$    ivegtfn = oneNamelistFile%ivegtfn
!!$    ixsctn = oneNamelistFile%ixsctn
!!$    iz0flg = oneNamelistFile%iz0flg
!!$    kwrite = oneNamelistFile%kwrite
!!$    lite_vars = oneNamelistFile%lite_vars
!!$    ndviflg = oneNamelistFile%ndviflg
!!$    ndvifn = oneNamelistFile%ndvifn
!!$    ndvifpfx = oneNamelistFile%ndvifpfx
!!$    nlite_vars = oneNamelistFile%nlite_vars
!!$    nofilflg = oneNamelistFile%nofilflg
!!$    nplt = oneNamelistFile%nplt
!!$    pastfn = oneNamelistFile%pastfn
!!$    sfcfiles = oneNamelistFile%sfcfiles
!!$    sstfpfx = oneNamelistFile%sstfpfx
!!$    timstr = oneNamelistFile%timstr
!!$    topfiles = oneNamelistFile%topfiles
!!$    toptenh = oneNamelistFile%toptenh
!!$    toptwvl = oneNamelistFile%toptwvl
!!$    xlite = oneNamelistFile%xlite
!!$    ylite = oneNamelistFile%ylite
!!$    z0fact = oneNamelistFile%z0fact
!!$    z0max = oneNamelistFile%z0max
!!$    zlite = oneNamelistFile%zlite
!!$    ifusflg = oneNamelistFile%ifusflg
!!$    ifusfn = oneNamelistFile%ifusfn
!!$    fusfiles = oneNamelistFile%fusfiles
!!$    gobrad = oneNamelistFile%gobrad
!!$    gobsep = oneNamelistFile%gobsep
!!$    gridwt = oneNamelistFile%gridwt
!!$    guess1st = oneNamelistFile%guess1st
!!$    hybbot = oneNamelistFile%hybbot
!!$    hybtop = oneNamelistFile%hybtop
!!$    i1st_flg = oneNamelistFile%i1st_flg
!!$    iapr = oneNamelistFile%iapr
!!$    iarawi = oneNamelistFile%iarawi
!!$    iasrfce = oneNamelistFile%iasrfce
!!$    igridfl = oneNamelistFile%igridfl
!!$    iobswin = oneNamelistFile%iobswin
!!$    ioflgisz = oneNamelistFile%ioflgisz
!!$    ioflgvar = oneNamelistFile%ioflgvar
!!$    isan_inc = oneNamelistFile%isan_inc
!!$    isfc_flg = oneNamelistFile%isfc_flg
!!$    iszstage = oneNamelistFile%iszstage
!!$    iupa_flg = oneNamelistFile%iupa_flg
!!$    ivrstage = oneNamelistFile%ivrstage
!!$    levth = oneNamelistFile%levth
!!$    maxsfc = oneNamelistFile%maxsfc
!!$    maxsta = oneNamelistFile%maxsta
!!$    nfeedvar = oneNamelistFile%nfeedvar
!!$    nigrids = oneNamelistFile%nigrids
!!$    nisn = oneNamelistFile%nisn
!!$    notid = oneNamelistFile%notid
!!$    notsta = oneNamelistFile%notsta
!!$    respon = oneNamelistFile%respon
!!$    sfcinf = oneNamelistFile%sfcinf
!!$    sigzwt = oneNamelistFile%sigzwt
!!$    stasep = oneNamelistFile%stasep
!!$    swvlnth = oneNamelistFile%swvlnth
!!$    topsigz = oneNamelistFile%topsigz
!!$    varpfx = oneNamelistFile%varpfx
!!$    wvlnth = oneNamelistFile%wvlnth
!!$    confrq = oneNamelistFile%confrq
!!$    cu_prefix = oneNamelistFile%cu_prefix
!!$    cu_tel = oneNamelistFile%cu_tel
!!$    cu_til = oneNamelistFile%cu_til
!!$    if_cuinv = oneNamelistFile%if_cuinv
!!$    nnqparm = oneNamelistFile%nnqparm
!!$    tcu_beg = oneNamelistFile%tcu_beg
!!$    tcu_end = oneNamelistFile%tcu_end
!!$    tnudcu = oneNamelistFile%tnudcu
!!$    wcldbs = oneNamelistFile%wcldbs
!!$    wt_cu_grid = oneNamelistFile%wt_cu_grid
!!$    raddatfn = oneNamelistFile%raddatfn
!!$    closure_type = oneNamelistFile%closure_type
!!$    centlat = oneNamelistFile%centlat
!!$    centlon = oneNamelistFile%centlon
!!$    cphas = oneNamelistFile%cphas
!!$    deltax = oneNamelistFile%deltax
!!$    deltay = oneNamelistFile%deltay
!!$    deltaz = oneNamelistFile%deltaz
!!$    distim = oneNamelistFile%distim
!!$    dtlong = oneNamelistFile%dtlong
!!$    dzmax = oneNamelistFile%dzmax
!!$    dzrat = oneNamelistFile%dzrat
!!$    expnme = oneNamelistFile%expnme
!!$    gridu = oneNamelistFile%gridu
!!$    gridv = oneNamelistFile%gridv
!!$    ibnd = oneNamelistFile%ibnd
!!$    icorflg = oneNamelistFile%icorflg
!!$    idate1 = oneNamelistFile%idate1
!!$    ideltat = oneNamelistFile%ideltat
!!$    if_adap = oneNamelistFile%if_adap
!!$    ihtran = oneNamelistFile%ihtran
!!$    imonth1 = oneNamelistFile%imonth1
!!$    initial = oneNamelistFile%initial
!!$    itime1 = oneNamelistFile%itime1
!!$    iyear1 = oneNamelistFile%iyear1
!!$    jbnd = oneNamelistFile%jbnd
!!$    lsflg = oneNamelistFile%lsflg
!!$    nacoust = oneNamelistFile%nacoust
!!$    naddsc = oneNamelistFile%naddsc
!!$    nestz1 = oneNamelistFile%nestz1
!!$    nestz2 = oneNamelistFile%nestz2
!!$    nfpt = oneNamelistFile%nfpt
!!$    ngrids = oneNamelistFile%ngrids
!!$    ninest = oneNamelistFile%ninest
!!$    njnest = oneNamelistFile%njnest
!!$    nknest = oneNamelistFile%nknest
!!$    nndtrat = oneNamelistFile%nndtrat
!!$    nnstbot = oneNamelistFile%nnstbot
!!$    nnsttop = oneNamelistFile%nnsttop
!!$    nnxp = oneNamelistFile%nnxp
!!$    nnyp = oneNamelistFile%nnyp
!!$    nnzp = oneNamelistFile%nnzp
!!$    npatch = oneNamelistFile%npatch
!!$    nstratx = oneNamelistFile%nstratx
!!$    nstraty = oneNamelistFile%nstraty
!!$    nstratz1 = oneNamelistFile%nstratz1
!!$    nstratz2 = oneNamelistFile%nstratz2
!!$    nxtnest = oneNamelistFile%nxtnest
!!$    nzg = oneNamelistFile%nzg
!!$    nzs = oneNamelistFile%nzs
!!$    polelat = oneNamelistFile%polelat
!!$    polelon = oneNamelistFile%polelon
!!$    runtype = oneNamelistFile%runtype
!!$    timeunit = oneNamelistFile%timeunit
!!$    timmax = oneNamelistFile%timmax
!!$    zz = oneNamelistFile%zz
!!$    albedo = oneNamelistFile%albedo
!!$    drtcon = oneNamelistFile%drtcon
!!$    dthcon = oneNamelistFile%dthcon
!!$    isfcl = oneNamelistFile%isfcl
!!$    nslcon = oneNamelistFile%nslcon
!!$    nvegpat = oneNamelistFile%nvegpat
!!$    nvgcon = oneNamelistFile%nvgcon
!!$    pctlcon = oneNamelistFile%pctlcon
!!$    seatmp = oneNamelistFile%seatmp
!!$    slmstr = oneNamelistFile%slmstr
!!$    slz = oneNamelistFile%slz
!!$    stgoff = oneNamelistFile%stgoff
!!$    zrough = oneNamelistFile%zrough
!!$    frqoda = oneNamelistFile%frqoda
!!$    if_oda = oneNamelistFile%if_oda
!!$    oda_sfc_tel = oneNamelistFile%oda_sfc_tel
!!$    oda_sfc_til = oneNamelistFile%oda_sfc_til
!!$    oda_sfcprefix = oneNamelistFile%oda_sfcprefix
!!$    oda_upa_tel = oneNamelistFile%oda_upa_tel
!!$    oda_upa_til = oneNamelistFile%oda_upa_til
!!$    oda_upaprefix = oneNamelistFile%oda_upaprefix
!!$    roda_hgt = oneNamelistFile%roda_hgt
!!$    roda_sfc0 = oneNamelistFile%roda_sfc0
!!$    roda_sfce = oneNamelistFile%roda_sfce
!!$    roda_upa0 = oneNamelistFile%roda_upa0
!!$    roda_upae = oneNamelistFile%roda_upae
!!$    roda_zfact = oneNamelistFile%roda_zfact
!!$    tnudoda = oneNamelistFile%tnudoda
!!$    todabeg = oneNamelistFile%todabeg
!!$    todaend = oneNamelistFile%todaend
!!$    wt_oda_grid = oneNamelistFile%wt_oda_grid
!!$    wt_oda_pi = oneNamelistFile%wt_oda_pi
!!$    wt_oda_rt = oneNamelistFile%wt_oda_rt
!!$    wt_oda_th = oneNamelistFile%wt_oda_th
!!$    wt_oda_uv = oneNamelistFile%wt_oda_uv
!!$    ilwrtyp = oneNamelistFile%ilwrtyp
!!$    iswrtyp = oneNamelistFile%iswrtyp
!!$    lonrad = oneNamelistFile%lonrad
!!$    radfrq = oneNamelistFile%radfrq
!!$    soil_moist = oneNamelistFile%soil_moist
!!$    soil_moist_fail = oneNamelistFile%soil_moist_fail
!!$    usdata_in = oneNamelistFile%usdata_in
!!$    usmodel_in = oneNamelistFile%usmodel_in
!!$    akmin = oneNamelistFile%akmin
!!$    csx = oneNamelistFile%csx
!!$    csz = oneNamelistFile%csz
!!$    idiffk = oneNamelistFile%idiffk
!!$    if_urban_canopy = oneNamelistFile%if_urban_canopy
!!$    ihorgrad = oneNamelistFile%ihorgrad
!!$    xkhkm = oneNamelistFile%xkhkm
!!$    zkhkm = oneNamelistFile%zkhkm
!!$    cond_hfile = oneNamelistFile%cond_hfile
!!$    nud_cond = oneNamelistFile%nud_cond
!!$    nud_hfile = oneNamelistFile%nud_hfile
!!$    nud_type = oneNamelistFile%nud_type
!!$    nudlat = oneNamelistFile%nudlat
!!$    t_nudge_rc = oneNamelistFile%t_nudge_rc
!!$    tcond_beg = oneNamelistFile%tcond_beg
!!$    tcond_end = oneNamelistFile%tcond_end
!!$    tnudcent = oneNamelistFile%tnudcent
!!$    tnudlat = oneNamelistFile%tnudlat
!!$    tnudtop = oneNamelistFile%tnudtop
!!$    varfpfx = oneNamelistFile%varfpfx
!!$    vwait1 = oneNamelistFile%vwait1
!!$    vwaittot = oneNamelistFile%vwaittot
!!$    wt_nudge_grid = oneNamelistFile%wt_nudge_grid
!!$    wt_nudge_pi = oneNamelistFile%wt_nudge_pi
!!$    wt_nudge_rt = oneNamelistFile%wt_nudge_rt
!!$    wt_nudge_th = oneNamelistFile%wt_nudge_th
!!$    wt_nudge_uv = oneNamelistFile%wt_nudge_uv
!!$    wt_nudgec_grid = oneNamelistFile%wt_nudgec_grid
!!$    znudtop = oneNamelistFile%znudtop
!!$    aparm = oneNamelistFile%aparm
!!$    coltabfn = oneNamelistFile%coltabfn
!!$    cparm = oneNamelistFile%cparm
!!$    gnu = oneNamelistFile%gnu
!!$    gparm = oneNamelistFile%gparm
!!$    hparm = oneNamelistFile%hparm
!!$    iaggr = oneNamelistFile%iaggr
!!$    icloud = oneNamelistFile%icloud
!!$    igraup = oneNamelistFile%igraup
!!$    ihail = oneNamelistFile%ihail
!!$    ipris = oneNamelistFile%ipris
!!$    irain = oneNamelistFile%irain
!!$    isnow = oneNamelistFile%isnow
!!$    level = oneNamelistFile%level
!!$    mkcoltab = oneNamelistFile%mkcoltab
!!$    pparm = oneNamelistFile%pparm
!!$    rparm = oneNamelistFile%rparm
!!$    sparm = oneNamelistFile%sparm
!!$    load_bal = oneNamelistFile%load_bal
!!$    hs = oneNamelistFile%hs
!!$    ipsflg = oneNamelistFile%ipsflg
!!$    irtsflg = oneNamelistFile%irtsflg
!!$    itsflg = oneNamelistFile%itsflg
!!$    iusflg = oneNamelistFile%iusflg
!!$    ps = oneNamelistFile%ps
!!$    rts = oneNamelistFile%rts
!!$    ts = oneNamelistFile%ts
!!$    us = oneNamelistFile%us
!!$    vs = oneNamelistFile%vs
!!$    nnshcu = oneNamelistFile%nnshcu
!!$    shcufrq = oneNamelistFile%shcufrq
!!$    co2_init = oneNamelistFile%co2_init
!!$    n_co2 = oneNamelistFile%n_co2
!!$    catt = oneNamelistFile%catt
!!$    firemapfn = oneNamelistFile%firemapfn
!!$    plumerise = oneNamelistFile%plumerise
!!$    define_proc = oneNamelistFile%define_proc
!!$    prfrq = oneNamelistFile%prfrq
!!$    recycle_tracers = oneNamelistFile%recycle_tracers
!!$    teb_spm = oneNamelistFile%teb_spm
!!$    ichemi = oneNamelistFile%ichemi
!!$    ichemi_in = oneNamelistFile%ichemi_in
!!$    chemdata_in = oneNamelistFile%chemdata_in
!!$    isource = oneNamelistFile%isource
!!$    weekdayin = oneNamelistFile%weekdayin
!!$    efsat = oneNamelistFile%efsat
!!$    efsun = oneNamelistFile%efsun
!!$    eindno = oneNamelistFile%eindno
!!$    eindno2 = oneNamelistFile%eindno2
!!$    eindpm = oneNamelistFile%eindpm
!!$    eindco = oneNamelistFile%eindco
!!$    eindso2 = oneNamelistFile%eindso2
!!$    eindvoc = oneNamelistFile%eindvoc
!!$    eveino = oneNamelistFile%eveino
!!$    eveino2 = oneNamelistFile%eveino2
!!$    eveipm = oneNamelistFile%eveipm
!!$    eveico = oneNamelistFile%eveico
!!$    eveiso2 = oneNamelistFile%eveiso2
!!$    eveivoc = oneNamelistFile%eveivoc
!!$    rushh1 = oneNamelistFile%rushh1
!!$    rushh2 = oneNamelistFile%rushh2
!!$    daylight = oneNamelistFile%daylight
!!$    iteb = oneNamelistFile%iteb
!!$    tminbld = oneNamelistFile%tminbld
!!$    nteb = oneNamelistFile%nteb
!!$    hc_roof = oneNamelistFile%hc_roof
!!$    tc_roof = oneNamelistFile%tc_roof
!!$    d_roof = oneNamelistFile%d_roof
!!$    hc_road = oneNamelistFile%hc_road
!!$    d_road = oneNamelistFile%d_road
!!$    tc_road = oneNamelistFile%tc_road
!!$    d_wall = oneNamelistFile%d_wall
!!$    tc_wall = oneNamelistFile%tc_wall
!!$    hc_wall = oneNamelistFile%hc_wall
!!$    nurbtype = oneNamelistFile%nurbtype
!!$    ileafcod = oneNamelistFile%ileafcod
!!$    z0_town = oneNamelistFile%z0_town
!!$    bld = oneNamelistFile%bld
!!$    bld_height = oneNamelistFile%bld_height
!!$    bld_hl_ratio = oneNamelistFile%bld_hl_ratio
!!$    aroof = oneNamelistFile%aroof
!!$    eroof = oneNamelistFile%eroof
!!$    aroad = oneNamelistFile%aroad
!!$    eroad = oneNamelistFile%eroad
!!$    awall = oneNamelistFile%awall
!!$    ewall = oneNamelistFile%ewall
!!$    htraf = oneNamelistFile%htraf
!!$    hindu = oneNamelistFile%hindu
!!$    pletraf = oneNamelistFile%pletraf
!!$    pleindu = oneNamelistFile%pleindu
!!$    domain_fname = oneNamelistFile%domain_fname
!!$  end subroutine StoreNamelistFile


  !**********************************************************************





  subroutine BroadcastNamelistFile(oneNamelistFile, oneParallelEnvironment)
    use ParLib, only: parf_bcast
    use ModParallelEnvironment, only: parallelEnvironment
    implicit none
    type(namelistFile), pointer :: oneNamelistFile
    type(parallelEnvironment), pointer :: oneParallelEnvironment

    include "i8.h"

    ! MODEL_GRIDS
    call parf_bcast(oneNamelistFile%expnme,&
         int(len(oneNamelistFile%expnme),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%runtype,&
         int(len(oneNamelistFile%runtype),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%timeunit,&
         int(len(oneNamelistFile%timeunit),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%timmax,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%load_bal,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%imonth1,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%idate1,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%iyear1,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%itime1,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%ngrids,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%nnxp,&
         int(size(oneNamelistFile%nnxp,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%nnyp,&
         int(size(oneNamelistFile%nnyp,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%nnzp,&
         int(size(oneNamelistFile%nnzp,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%nzg,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%nzs,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%nxtnest,&
         int(size(oneNamelistFile%nxtnest,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%domain_fname,&
         int(len(oneNamelistFile%domain_fname),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%if_adap,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%ihtran,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%deltax,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%deltay,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%deltaz,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%dzrat,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%dzmax,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%fixLevels,&
         oneParallelEnvironment%master_num)    
    call parf_bcast(oneNamelistFile%zz,&
         int(size(oneNamelistFile%zz,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%dtlong,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%nacoust,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%ideltat,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%nstratx,&
         int(size(oneNamelistFile%nstratx,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%nstraty,&
         int(size(oneNamelistFile%nstraty,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%nndtrat,&
         int(size(oneNamelistFile%nndtrat,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%nestz1,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%nstratz1,&
         int(size(oneNamelistFile%nstratz1,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%nestz2,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%nstratz2,&
         int(size(oneNamelistFile%nstratz2,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%polelat,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%polelon,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%centlat,&
         int(size(oneNamelistFile%centlat,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%centlon,&
         int(size(oneNamelistFile%centlat,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%ninest,&
         int(size(oneNamelistFile%ninest,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%njnest,&
         int(size(oneNamelistFile%njnest,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%nknest,&
         int(size(oneNamelistFile%nknest,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%nnsttop,&
         int(size(oneNamelistFile%nnsttop,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%nnstbot,&
         int(size(oneNamelistFile%nnstbot,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%gridu,&
         int(size(oneNamelistFile%gridu,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%gridv,&
         int(size(oneNamelistFile%gridv,1),i8),&
         oneParallelEnvironment%master_num)

    ! CCATT_INFO
    call parf_bcast(oneNamelistFile%ccatt,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%chemistry,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%split_method,&
         int(len(oneNamelistFile%split_method),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%chem_timestep,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%chemistry_aq,&
         oneParallelEnvironment%master_num)
!    call parf_bcast(oneNamelistFile%aerosol,&
!         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%chem_assim,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%srcmapfn,&
         int(len(oneNamelistFile%srcmapfn),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%recycle_tracers,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%def_proc_src,&
         int(len(oneNamelistFile%def_proc_src),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%diur_cycle,&
         int(size(oneNamelistFile%diur_cycle,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%na_extra2d,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%na_extra3d,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%plumerise,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%prfrq,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%volcanoes,&
         oneParallelEnvironment%master_num)
!Matrix
    call parf_bcast(oneNamelistFile%aerosol,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%aer_timestep,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%aer_assim,&
         oneParallelEnvironment%master_num)

    call parf_bcast(oneNamelistFile%mech,&
         oneParallelEnvironment%master_num)
!--(DMK-CCATT-OLD)-----------------------------------------------------------
!    ! CATT_INFO
!    call parf_bcast(oneNamelistFile%catt,&
!         oneParallelEnvironment%master_num)
!    call parf_bcast(oneNamelistFile%firemapfn,&
!         int(len(oneNamelistFile%firemapfn),i8),&
!         oneParallelEnvironment%master_num)
!    call parf_bcast(oneNamelistFile%recycle_tracers,&
!         oneParallelEnvironment%master_num)
!    call parf_bcast(oneNamelistFile%plumerise,&
!         oneParallelEnvironment%master_num)
!    call parf_bcast(oneNamelistFile%define_proc,&
!         oneParallelEnvironment%master_num)
!    call parf_bcast(oneNamelistFile%prfrq,&
!         oneParallelEnvironment%master_num)
!--(DMK-CCATT-FIM)-----------------------------------------------------------
    ! TEB_SPM_INFO
    call parf_bcast(oneNamelistFile%teb_spm,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%fusfiles,&
         int(len(oneNamelistFile%fusfiles),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%ifusflg,&
         int(size(oneNamelistFile%ifusflg,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%IFUSFN,&
         int(len(oneNamelistFile%IFUSFN(1)),i8), &
         int(size(oneNamelistFile%IFUSFN,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%ichemi,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%ichemi_in,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%chemdata_in,&
         int(len(oneNamelistFile%chemdata_in),i8), &
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%isource,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%weekdayin,&
         int(len(oneNamelistFile%weekdayin),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%rushh1,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%rushh2,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%daylight,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%efsat,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%efsun,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%eindno,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%eindno2,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%eindpm,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%eindco,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%eindso2,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%eindvoc,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%eveino,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%eveino2,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%eveipm,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%eveico,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%eveiso2,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%eveivoc,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%iteb,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%tminbld,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%nteb,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%hc_roof,&
         int(size(oneNamelistFile%hc_roof,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%tc_roof,&
         int(size(oneNamelistFile%tc_roof,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%d_roof,&
         int(size(oneNamelistFile%d_roof,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%hc_road,&
         int(size(oneNamelistFile%hc_road,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%tc_road,&
         int(size(oneNamelistFile%tc_road,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%d_road,&
         int(size(oneNamelistFile%d_road,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%hc_wall,&
         int(size(oneNamelistFile%hc_wall,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%tc_wall,&
         int(size(oneNamelistFile%tc_wall,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%d_wall,&
         int(size(oneNamelistFile%d_wall,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%nurbtype,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%ileafcod,&
         int(size(oneNamelistFile%ileafcod,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%z0_town,&
         int(size(oneNamelistFile%z0_town,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%bld,&
         int(size(oneNamelistFile%bld,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%bld_height,&
         int(size(oneNamelistFile%bld_height,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%bld_hl_ratio,&
         int(size(oneNamelistFile%bld_hl_ratio,1),i8), &
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%aroof,&
         int(size(oneNamelistFile%aroof,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%eroof,&
         int(size(oneNamelistFile%eroof,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%aroad,&
         int(size(oneNamelistFile%aroad,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%eroad,&
         int(size(oneNamelistFile%eroad,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%awall,&
         int(size(oneNamelistFile%awall,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%ewall,&
         int(size(oneNamelistFile%ewall,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%htraf,&
         int(size(oneNamelistFile%htraf,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%hindu,&
         int(size(oneNamelistFile%hindu,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%pletraf,&
         int(size(oneNamelistFile%pletraf,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%pleindu,&
         int(size(oneNamelistFile%pleindu,1),i8),&
         oneParallelEnvironment%master_num)
    ! MODEL_FILE_INFO
    call parf_bcast(oneNamelistFile%initial,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%nud_type,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%varfpfx,&
         int(len(oneNamelistFile%varfpfx),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%vwait1,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%vwaittot,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%nud_hfile,&
         int(len(oneNamelistFile%nud_hfile),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%nudlat,&
         oneParallelEnvironment%master_num)

    call parf_bcast(oneNamelistFile%timeWindowIAU,&
         oneParallelEnvironment%master_num)    
    call parf_bcast(oneNamelistFile%ramp,&
         oneParallelEnvironment%master_num)    
    call parf_bcast(oneNamelistFile%tnudlat,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%tnudcent,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%tnudtop,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%znudtop,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%wt_nudge_grid,&
         int(size(oneNamelistFile%wt_nudge_grid,1),i8), &
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%wt_nudge_uv,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%wt_nudge_th,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%wt_nudge_pi,&
         oneParallelEnvironment%master_num)

    call parf_bcast(oneNamelistFile%wt_nudge_rt,&
         oneParallelEnvironment%master_num)

    call parf_bcast(oneNamelistFile%applyIAU,&
         oneParallelEnvironment%master_num)
    
    call parf_bcast(oneNamelistFile%fileNameIAU,&
         int(len(oneNamelistFile%fileNameIAU),i8),&
         oneParallelEnvironment%master_num)


    call parf_bcast(oneNamelistFile%nud_cond,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%cond_hfile,&
         int(len(oneNamelistFile%cond_hfile),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%tcond_beg,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%tcond_end,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%t_nudge_rc,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%wt_nudgec_grid,&
         int(size(oneNamelistFile%wt_nudgec_grid,1),i8), &
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%if_oda,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%oda_upaprefix,&
         int(len(oneNamelistFile%oda_upaprefix),i8), &
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%oda_sfcprefix,&
         int(len(oneNamelistFile%oda_sfcprefix),i8), &
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%frqoda,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%todabeg,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%todaend,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%tnudoda,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%wt_oda_grid,&
         int(size(oneNamelistFile%wt_oda_grid,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%wt_oda_uv,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%wt_oda_th,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%wt_oda_pi,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%wt_oda_rt,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%roda_sfce,&
         int(size(oneNamelistFile%roda_sfce,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%roda_sfc0,&
         int(size(oneNamelistFile%roda_sfc0,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%roda_upae,&
         int(size(oneNamelistFile%roda_upae,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%roda_upa0,&
         int(size(oneNamelistFile%roda_upa0,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%roda_hgt,&
         int(size(oneNamelistFile%roda_hgt,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%roda_zfact,&
         int(size(oneNamelistFile%roda_zfact,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%oda_sfc_til,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%oda_sfc_tel,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%oda_upa_til,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%oda_upa_tel,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%if_cuinv,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%cu_prefix,&
         int(len(oneNamelistFile%cu_prefix),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%tnudcu,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%wt_cu_grid,&
         int(size(oneNamelistFile%wt_cu_grid,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%tcu_beg,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%tcu_end,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%cu_tel,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%cu_til,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%timstr,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%hfilin,&
         int(len(oneNamelistFile%hfilin),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%ipastin,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%pastfn,&
         int(len(oneNamelistFile%pastfn),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%ioutput,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%hfilout,&
         int(len(oneNamelistFile%hfilout),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%afilout,&
         int(len(oneNamelistFile%afilout),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%iclobber,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%ihistdel,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%frqhis,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%frqanl,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%frqlite,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%ipos,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%xlite,&
         int(len(oneNamelistFile%xlite),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%ylite,&
         int(len(oneNamelistFile%ylite),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%zlite,&
         int(len(oneNamelistFile%zlite),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%nlite_vars,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%lite_vars,&
         int(len(oneNamelistFile%lite_vars(1)),i8), &
         int(size(oneNamelistFile%lite_vars,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%avgtim,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%frqmean,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%frqboth,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%kwrite,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%frqprt,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%initfld,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%prtcputime,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%topfiles,&
         int(len(oneNamelistFile%topfiles),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%sfcfiles,&
         int(len(oneNamelistFile%sfcfiles),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%sstfpfx,&
         int(len(oneNamelistFile%sstfpfx),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%ndvifpfx,&
         int(len(oneNamelistFile%ndvifpfx),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%itoptflg,&
         int(size(oneNamelistFile%itoptflg,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%isstflg,&
         int(size(oneNamelistFile%isstflg,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%ivegtflg,&
         int(size(oneNamelistFile%ivegtflg,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%isoilflg,&
         int(size(oneNamelistFile%isoilflg,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%ndviflg,&
         int(size(oneNamelistFile%ndviflg,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%nofilflg,&
         int(size(oneNamelistFile%nofilflg,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%iupdndvi,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%iupdsst,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%itoptfn,&
         int(len(oneNamelistFile%itoptfn(1)),i8), &
         int(size(oneNamelistFile%itoptfn,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%isstfn,&
         int(len(oneNamelistFile%isstfn(1)),i8), &
         int(size(oneNamelistFile%isstfn,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%ivegtfn,&
         int(len(oneNamelistFile%ivegtfn(1)),i8), &
         int(size(oneNamelistFile%ivegtfn,1),i8), &
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%isoilfn,&
         int(len(oneNamelistFile%isoilfn(1)),i8), &
         int(size(oneNamelistFile%isoilfn,1),i8), &
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%ndvifn,&
         int(len(oneNamelistFile%ndvifn(1)),i8), &
         int(size(oneNamelistFile%ndvifn,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%itopsflg,&
         int(size(oneNamelistFile%itopsflg,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%toptenh,&
         int(size(oneNamelistFile%toptenh,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%toptwvl,&
         int(size(oneNamelistFile%toptwvl,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%iz0flg,&
         int(size(oneNamelistFile%iz0flg,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%z0max,&
         int(size(oneNamelistFile%z0max,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%z0fact,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%mkcoltab,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%coltabfn,&
         int(len(oneNamelistFile%coltabfn),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%mapaotfile,&
         int(len(oneNamelistFile%mapaotfile),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%julesin,&
         int(len(oneNamelistFile%julesin),i8),&
         oneParallelEnvironment%master_num)
    ! MODEL_OPTIONS
!--(DMK-CCATT-INI)-----------------------------------------------------------
    call parf_bcast(oneNamelistFile%dyncore_flag,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%pd_or_mnt_constraint, &
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%order_h, &
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%order_v,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%advmnt,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%ghostzonelength,&
         oneParallelEnvironment%master_num)

    call parf_bcast(oneNamelistFile%naddsc,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%icorflg,&
         oneParallelEnvironment%master_num)

!--(DMK-CCATT-INI)-----------------------------------------------------------
    call parf_bcast(oneNamelistFile%iexev,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%imassflx,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%vveldamp,&
         oneParallelEnvironment%master_num)
!--(DMK-CCATT-FIM)-----------------------------------------------------------
    call parf_bcast(oneNamelistFile%ibnd,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%jbnd,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%cphas,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%lsflg,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%nfpt,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%distim,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%sfire,&
         oneParallelEnvironment%master_num)    
    call parf_bcast(oneNamelistFile%iswrtyp,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%ilwrtyp,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%raddatfn,&
         int(len(oneNamelistFile%raddatfn),i8),&
         oneParallelEnvironment%master_num)

    call parf_bcast(oneNamelistFile%radfrq,&
         oneParallelEnvironment%master_num)

    call parf_bcast(oneNamelistFile%radtun,&
         oneParallelEnvironment%master_num)

    call parf_bcast(oneNamelistFile%lonrad,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%nnqparm,&
         int(size(oneNamelistFile%nnqparm,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%closure_type,&
         int(len(oneNamelistFile%closure_type),i8), &
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%nnshcu,&
         int(size(oneNamelistFile%nnshcu,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%confrq,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%shcufrq,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%wcldbs,&
         oneParallelEnvironment%master_num)

     call parf_bcast(oneNamelistFile%g3d_spread,&
         oneParallelEnvironment%master_num)
     call parf_bcast(oneNamelistFile%g3d_smoothh,&
         oneParallelEnvironment%master_num)
     call parf_bcast(oneNamelistFile%g3d_smoothv,&
         oneParallelEnvironment%master_num)

    call parf_bcast(oneNamelistFile%npatch,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%nvegpat,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%isfcl,&
         oneParallelEnvironment%master_num)

    call parf_bcast(oneNamelistFile%isfcl_ocean,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%nvgcon,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%pctlcon,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%nslcon,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%drtcon,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%zrough,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%albedo,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%seatmp,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%dthcon,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%soil_moist,&
         int(len(oneNamelistFile%soil_moist),i8), &
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%soil_moist_fail,&
         int(len(oneNamelistFile%soil_moist_fail),i8), &
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%usdata_in,&
         int(len(oneNamelistFile%usdata_in),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%usmodel_in,&
         int(len(oneNamelistFile%usmodel_in),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%slz,&
         int(size(oneNamelistFile%slz,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%slmstr,&
         int(size(oneNamelistFile%slmstr,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%stgoff,&
         int(size(oneNamelistFile%stgoff,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%if_urban_canopy,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%idiffk,&
         int(size(oneNamelistFile%idiffk,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%ihorgrad,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%csx,&
         int(size(oneNamelistFile%csx,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%csz,&
         int(size(oneNamelistFile%csz,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%xkhkm,&
         int(size(oneNamelistFile%xkhkm,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%zkhkm,&
         int(size(oneNamelistFile%zkhkm,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%akmin,&
         int(size(oneNamelistFile%akmin,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%mcphys_type,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%level,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%icloud,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%idriz,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%irime,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%iplaws,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%iccnlev,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%irain,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%ipris,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%isnow,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%iaggr,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%igraup,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%ihail,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%cparm,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%rparm,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%pparm,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%sparm,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%aparm,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%gparm,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%hparm,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%dparm,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%cnparm,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%gnparm,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%epsil,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%gnu,&
         int(size(oneNamelistFile%gnu,1),i8),&
         oneParallelEnvironment%master_num)
    !-------
    call parf_bcast(oneNamelistFile%windfarm,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%wffile,&
         int(len(oneNamelistFile%wffile),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%damModule,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%frqPrecip,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%damOutPrefix,&
         int(len(oneNamelistFile%damOutPrefix),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%evaluatePrefix,&
         int(len(oneNamelistFile%evaluatePrefix),i8),&
         oneParallelEnvironment%master_num)
     call parf_bcast(oneNamelistFile%evaluate,&
         oneParallelEnvironment%master_num)   
    !-------
    ! MODEL_SOUND
    call parf_bcast(oneNamelistFile%ipsflg,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%itsflg,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%irtsflg,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%iusflg,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%hs,&
         int(size(oneNamelistFile%hs,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%ps,&
         int(size(oneNamelistFile%hs,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%ts,&
         int(size(oneNamelistFile%ts,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%rts,&
         int(size(oneNamelistFile%rts,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%us,&
         int(size(oneNamelistFile%us,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%vs,&
         int(size(oneNamelistFile%vs,1),i8),&
         oneParallelEnvironment%master_num)
    ! MODEL_PRINT
    call parf_bcast(oneNamelistFile%nplt,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%iplfld,&
         int(len(oneNamelistFile%iplfld(1)),i8), &
         int(size(oneNamelistFile%iplfld,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%ixsctn,&
         int(size(oneNamelistFile%ixsctn,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%isbval,&
         int(size(oneNamelistFile%isbval,1),i8),&
         oneParallelEnvironment%master_num)
    ! ISAN_CONTROL
    call parf_bcast(oneNamelistFile%iszstage,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%ivrstage,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%isan_inc,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%guess1st,&
         int(len(oneNamelistFile%guess1st),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%i1st_flg,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%iupa_flg,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%isfc_flg,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%iapr,&
         int(len(oneNamelistFile%iapr),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%iarawi,&
         int(len(oneNamelistFile%iarawi),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%iasrfce,&
         int(len(oneNamelistFile%iasrfce),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%varpfx,&
         int(len(oneNamelistFile%varpfx),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%ioflgisz,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%ioflgvar,&
         oneParallelEnvironment%master_num)
    ! ISAN_ISENTROPIC
    call parf_bcast(oneNamelistFile%nisn,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%levth,&
         int(size(oneNamelistFile%levth,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%nigrids,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%topsigz,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%hybbot,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%hybtop,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%sfcinf,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%sigzwt,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%nfeedvar,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%maxsta,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%maxsfc,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%notsta,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%notid,&
         int(len(oneNamelistFile%notid(1)),i8), &
         int(size(oneNamelistFile%notid,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%iobswin,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%stasep,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%igridfl,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%gridwt,&
         int(size(oneNamelistFile%gridwt,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%gobsep,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%gobrad,&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%wvlnth,&
         int(size(oneNamelistFile%wvlnth,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%swvlnth,&
         int(size(oneNamelistFile%swvlnth,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%respon,&
         int(size(oneNamelistFile%respon,1),i8),&
         oneParallelEnvironment%master_num)
    !
    call parf_bcast(oneNamelistFile%icFileType, &
        oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%icPrefix, &
         int(len(oneNamelistFile%icPrefix),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%wind_u_varname , & 
         int(len(oneNamelistFile%wind_u_varname),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%wind_v_varname , &
         int(len(oneNamelistFile%wind_v_varname),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%temperature_varname, &
         int(len(oneNamelistFile%temperature_varname),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%geo_varname, &
         int(len(oneNamelistFile%geo_varname),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%ur_varname  , &  
         int(len(oneNamelistFile%ur_varname),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%initial_latitude, &
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%final_latitude, &
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%initial_longitude, &
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%final_longitude, &
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%z_max_level, &
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%dlimit, &
         int(size(oneNamelistFile%dlimit,1),i8),&
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%ulimit, &   
         int(size(oneNamelistFile%ulimit,1),i8),& 
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%scale_factor, &   
         int(size(oneNamelistFile%scale_factor,1),i8),& 
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%ccGradsWrite, &
         oneParallelEnvironment%master_num)    
    call parf_bcast(oneNamelistFile%icGradsPrefix, &
         int(len(oneNamelistFile%icGradsPrefix),i8),&
         oneParallelEnvironment%master_num)
    ! POST
    call parf_bcast(oneNamelistFile%nvp, &
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%vp, &
         int(len(oneNamelistFile%vp(1)),i8), &
         int(size(oneNamelistFile%vp,1),i8), &
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%gprefix, &
         int(len(oneNamelistFile%gprefix),i8), &
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%csvFile, &
         int(len(oneNamelistFile%csvFile),i8), &
         oneParallelEnvironment%master_num)    
    call parf_bcast(oneNamelistFile%anl2gra, &
         int(len(oneNamelistFile%anl2gra),i8), &
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%proj, &
         int(len(oneNamelistFile%proj),i8), &
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%mean_type, &
         int(len(oneNamelistFile%mean_type),i8), &
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%lati, &
         int(maxgrds,i8), &
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%loni, &
         int(maxgrds,i8), &
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%latf, &
         int(maxgrds,i8), &
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%lonf, &
         int(maxgrds,i8), &
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%zlevmax, &
         int(maxgrds,i8), &
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%ipresslev, &
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%inplevs, &
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%iplevs, &
         int(nzpmax,i8), &
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%mechanism, &
         int(len(oneNamelistFile%mechanism),i8), &
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%ascii_data, &
         int(len(oneNamelistFile%ascii_data),i8), &
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%site_lat, &
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%site_lon, &
         oneParallelEnvironment%master_num)

 !digital filter

    call parf_bcast(oneNamelistFile%applyDigitalFilter, &
         oneParallelEnvironment%master_num)
    call parf_bcast(oneNamelistFile%digitalFilterTimeWindow, &
         oneParallelEnvironment%master_num)

 !meteogram

   call parf_bcast(oneNamelistFile%applyMeteogram, &
                   oneParallelEnvironment%master_num)

 call parf_bcast(oneNamelistFile%meteogramFreq, &
         oneParallelEnvironment%master_num)

  call parf_bcast(oneNamelistFile%meteogramMap, &
         int(len(oneNamelistFile%meteogramMap),i8), &
         oneParallelEnvironment%master_num)

  call parf_bcast(oneNamelistFile%meteogramDir, &
         int(len(oneNamelistFile%meteogramDir),i8), &
         oneParallelEnvironment%master_num)

  end subroutine BroadcastNamelistFile



  subroutine TimeUnitsToSeconds(oneNamelistFile)
    type(namelistFile), pointer :: oneNamelistFile

    real :: tfact
    character(len=*), parameter :: h="**(TimeUnitsToSeconds)**"
    character(len=*), parameter :: header="**(TimeUnitsToSeconds)**"
    character(len=*), parameter :: version=""

    select case (oneNamelistFile%timeunit)
    case ("d","D")
       tfact = 86400.
    case ("h","H")
       tfact =  3600.
    case ("m","M")
       tfact =    60.
    case ("s","S")
       tfact =     1.
    case default
       call fatal_error(h//" Namelist timeunit ="//&
            &trim(oneNamelistFile%timeunit)//" has to be d, h, m or s")
    end select

    oneNamelistFile%timmax=oneNamelistFile%timmax*tfact
    oneNamelistFile%timstr=oneNamelistFile%timstr*tfact
  end subroutine TimeUnitsToSeconds




  subroutine DumpNamelistFile(oneNamelistFile,nmachs,mchnum,master_num)
    use modPrintInitial

    implicit none
    type(namelistFile), pointer :: oneNamelistFile

    character(len=*),parameter :: revision='6.0'
    character(len=*),parameter :: license='CC Attribution-ShareAlike 4.0 International'

    include "constants.f90"

    integer, intent(in) :: nmachs,mchnum,master_num

    integer :: ng,np,k,m,ngrids
    integer :: lenLine, lenThis,sofname
    character(len=1) :: ck

    ! This routine prints out a listing of the values of all variables
    ! in the NAMELISTS


    iErrNumber=bramsHeader(revision,license,nmachs,trim(oneNamelistFile%runtype),oneNamelistFile%fileName &
                          ,mchnum,master_num,oneNamelistFile%iyear1,oneNamelistFile%imonth1,oneNamelistFile%idate1 &
                          ,oneNamelistFile%itime1,oneNamelistFile%expnme,oneNamelistFile%ngrids &
                          ,oneNamelistFile%timmax,oneNamelistFile%timeunit)

    ngrids=oneNamelistFile%ngrids

    write(*,fmt='(A)') 
    write(*,"(a)") "+-------------------------------------------- Variables from RAMSIN ---------------------------------------------------"
    write(*,fmt='(A)') 

    iErrNumber=printGridHeader(ngrids)
    iErrNumber=printOneVarGrid(ngrids,'NNXP       ',oneNamelistFile%nnxp    ,"I14.4")
    iErrNumber=printOneVarGrid(ngrids,'NNYP       ',oneNamelistFile%nnyp    ,"I14.4")
    iErrNumber=printOneVarGrid(ngrids,'NNZP       ',oneNamelistFile%nnzp    ,"I14.4")
    iErrNumber=printOneVarGrid(ngrids,'NXTNEST    ',oneNamelistFile%nxtnest ,"I14.4")
    iErrNumber=printOneVarGrid(ngrids,'NSTRATX    ',oneNamelistFile%nstratx ,"I14.4")
    iErrNumber=printOneVarGrid(ngrids,'NSTRATY    ',oneNamelistFile%nstraty ,"I14.4")
    iErrNumber=printOneVarGrid(ngrids,'NSTRATZ1   ',oneNamelistFile%nstratz1 ,"I14.4")
    iErrNumber=printOneVarGrid(ngrids,'NSTRATZ2   ',oneNamelistFile%nstratz2 ,"I14.4")
    iErrNumber=printOneVarGrid(ngrids,'IDIFFK     ',oneNamelistFile%idiffk  ,"I14.4")
    iErrNumber=printOneVarGrid(ngrids,'NNDTRAT    ',oneNamelistFile%nndtrat ,"I14.4")
    iErrNumber=printOneVarGrid(ngrids,'NINEST     ',oneNamelistFile%ninest  ,"I14.4")
    iErrNumber=printOneVarGrid(ngrids,'NJNEST     ',oneNamelistFile%njnest  ,"I14.4")
    iErrNumber=printOneVarGrid(ngrids,'NKNEST     ',oneNamelistFile%nknest  ,"I14.4")
    iErrNumber=printOneVarGrid(ngrids,'NNSTTOP    ',oneNamelistFile%nnsttop ,"I14.4")
    iErrNumber=printOneVarGrid(ngrids,'NNSTBOT    ',oneNamelistFile%nnstbot ,"I14.4")      
    iErrNumber=printOneVarGrid(ngrids,'CENTLAT    ',oneNamelistFile%centlat ,"F14.6")
    iErrNumber=printOneVarGrid(ngrids,'CENTLON    ',oneNamelistFile%centlon ,"F14.6")
    iErrNumber=printOneVarGrid(ngrids,'NNQPARM    ',oneNamelistFile%nnqparm ,"I14.1")
    iErrNumber=printOneVarGrid(ngrids,'NNSHCU     ',oneNamelistFile%nnshcu  ,"I14.1") 
    iErrNumber=printOneVarGrid(ngrids,'ITOPSFLG   ',oneNamelistFile%itopsflg,"I14.1")
    iErrNumber=printOneVarGrid(ngrids,'TOPTENH    ',oneNamelistFile%toptenh,"F14.1") 
    iErrNumber=printOneVarGrid(ngrids,'TOPTWVL    ',oneNamelistFile%toptwvl,"F14.1") 
    iErrNumber=printOneVarGrid(ngrids,'IZ0FLG     ',oneNamelistFile%iz0flg ,"I14") 
    iErrNumber=printOneVarGrid(ngrids,'Z0MAX      ',oneNamelistFile%z0max  ,"F14.1")
    iErrNumber=printOneVarGrid(ngrids,'WT_NUDGE_GR',oneNamelistFile%wt_nudge_grid,"F14.1")
    iErrNumber=printOneVarGrid(ngrids,'CSX        ',oneNamelistFile%csx    ,"F14.2") 
    iErrNumber=printOneVarGrid(ngrids,'CSZ        ',oneNamelistFile%csz    ,"F14.2") 
    iErrNumber=printOneVarGrid(ngrids,'XKHKM      ',oneNamelistFile%xkhkm  ,"F14.2")
    iErrNumber=printOneVarGrid(ngrids,'ZKHKM      ',oneNamelistFile%zkhkm  ,"F14.2")
    iErrNumber=printOneVarGrid(ngrids,'AKMIN      ',oneNamelistFile%akmin  ,"F14.2")
    iErrNumber=printOneVarGrid(ngrids,'GRIDU      ',oneNamelistFile%gridu  ,"E14.5")
    iErrNumber=printOneVarGrid(ngrids,'GRIDV      ',oneNamelistFile%gridv  ,"E14.5")      
    iErrNumber=printOneVarGrid(ngrids,'ITOPTFLG',oneNamelistFile%itoptflg,"I14.1") 
    iErrNumber=printOneVarGrid(ngrids,'ISSTFLG ',oneNamelistFile%isstflg,"I14.1") 
    iErrNumber=printOneVarGrid(ngrids,'IVEGTFLG',oneNamelistFile%ivegtflg,"I14.1") 
    iErrNumber=printOneVarGrid(ngrids,'ISOILFLG',oneNamelistFile%isoilflg,"I14.1") 
    iErrNumber=printOneVarGrid(ngrids,'NDVIFLG ',oneNamelistFile%ndviflg,"I14.1") 
    iErrNumber=printOneVarGrid(ngrids,'NOFILFLG',oneNamelistFile%nofilflg,"I14.1") 
    if(oneNamelistFile%ipos>0) then
      iErrNumber=printOneVarGrid(ngrids,'ZLEVMAX',oneNamelistFile%zlevmax,"I14.3") 
      iErrNumber=printOneVarGrid(ngrids,'LATI',oneNamelistFile%lati,'F14.2')
      iErrNumber=printOneVarGrid(ngrids,'LATF',oneNamelistFile%latf,'F14.2')
      iErrNumber=printOneVarGrid(ngrids,'LONI',oneNamelistFile%loni,'F14.2')
      iErrNumber=printOneVarGrid(ngrids,'LONF',oneNamelistFile%lonf,'F14.2')
    endif
    iErrNumber=printGridTail(ngrids)

    iErrNumber=printVarHeader(4)
    !Block1
    iErrNumber=printOneLineVars(4,(/      &
       conv2String('POLELAT',"")          &
      ,conv2String('POLELON',"")          &
      ,conv2String('NZG    ',"")          &
      ,conv2String('NZS    ',"")/),(/     &
       conv2String(oneNamelistFile%polelat     ,"F14.3"    ) &
      ,conv2String(oneNamelistFile%polelon     ,"F14.3"    ) &
      ,conv2String(oneNamelistFile%nzg         ,"I14.2"      ) &
      ,conv2String(oneNamelistFile%nzs         ,"I14.2"      )/))
    !Block2
    iErrNumber=printOneLineVars(4,(/      &
       conv2String('DELTAX ',"")          &
      ,conv2String('DELTAY ',"")          &
      ,conv2String('DTLONG ',"")          &
      ,conv2String('IDELTAT',"")/),(/     &
       conv2String(oneNamelistFile%deltax ,"F14.1") &
      ,conv2String(oneNamelistFile%deltay ,"F14.1") &
      ,conv2String(oneNamelistFile%dtlong ,"F14.1") &
      ,conv2String(oneNamelistFile%ideltat,"I14.1"  )/))    
    !Block3
    iErrNumber=printOneLineVars(4,(/      &
       conv2String('DELTAZ   ',"")          &
      ,conv2String('DZRAT    ',"")          &
      ,conv2String('DZMAX    ',"")          &
      ,conv2String('FIXLEVELS',"")/),(/     &
       conv2String(oneNamelistFile%deltaz   ,"F14.1") &
      ,conv2String(oneNamelistFile%dzrat    ,"F14.2") &
      ,conv2String(oneNamelistFile%dzmax    ,"F14.1") &
      ,conv2String(oneNamelistFile%fixlevels,"I14"  )/))    
    !Block4
    iErrNumber=printOneLineVars(4,(/      &
       conv2String('NESTZ1      ',"")          &
      ,conv2String('NESTZ2      ',"")          &
      ,conv2String('IF_ADAP     ',"")          &
      ,conv2String('IHTRAN      ',"")/),(/     &
       conv2String(oneNamelistFile%nestz1       ,"I14.1") &
      ,conv2String(oneNamelistFile%nestz2       ,"I14.1") &
      ,conv2String(oneNamelistFile%if_adap      ,"I14.1") &
      ,conv2String(oneNamelistFile%ihtran       ,"I14.1"  )/))
    !Block5
    iErrNumber=printOneLineVars(4,(/      &
       conv2String('NACOUST'                    ,"")          &
      ,conv2String('DYNCORE_FLAG'               ,"")          &
      ,conv2String('ADVMNT     '                ,"")          &
      ,conv2String('PD_OR_MNT_CONSTRAINT'       ,"")/),(/     &
       conv2String(oneNamelistFile%nacoust                   ,"I14.1") &
      ,conv2String(oneNamelistFile%dyncore_flag              ,"I14.1") &
      ,conv2String(oneNamelistFile%advmnt                    ,"I14.1") &
      ,conv2String(oneNamelistFile%pd_or_mnt_constraint      ,"I14.1"  )/))
    !Block6
    iErrNumber=printOneLineVars(4,(/      &
       conv2String('ORDER_H    '       ,"")          &
      ,conv2String('ORDER_V    '       ,"")          &
      ,conv2String('GHOSTZONE  '       ,"")          &
      ,conv2String('IPOS       '       ,"")/),(/     &
       conv2String(oneNamelistFile%order_h        ,"I14.1") &
      ,conv2String(oneNamelistFile%order_v        ,"I14.1") &
      ,conv2String(oneNamelistFile%ghostzonelength,"I14.1") &
      ,conv2String(oneNamelistFile%ipos           ,"I14.1"  )/))
    !Block7
    iErrNumber=printOneLineVars(4,(/      &
       conv2String('IOUTPUT      '       ,"")          &
      ,conv2String('ICLOBBER     '       ,"")          &
      ,conv2String('IHISTDEL     '       ,"")          &
      ,conv2String('FRQHIS       '       ,"")/),(/     &
       conv2String(oneNamelistFile%ioutput ,"I14.1") &
      ,conv2String(oneNamelistFile%iclobber,"I14.1") &
      ,conv2String(oneNamelistFile%ihistdel,"I14.1") &
      ,conv2String(oneNamelistFile%frqhis  ,"F14.1"  )/))
    !Block8
    iErrNumber=printOneLineVars(4,(/      &
       conv2String('FRQANL    '       ,"")          &
      ,conv2String('APPLYIAU  '       ,"")          &
      ,conv2String('TIMEW  IAU'       ,"")          &
      ,conv2String('RAMP      '       ,"")/),(/     &
       conv2String(oneNamelistFile%frqanl       ,"F14.1") &
      ,conv2String(oneNamelistFile%applyIAU     ,"I14.1") &
      ,conv2String(oneNamelistFile%timeWindowIAU,"F14.1") &
      ,conv2String(oneNamelistFile%ramp         ,"F14.1"  )/))
    !Block9
    iErrNumber=printOneLineVars(4,(/      &
       conv2String('NUD_TYPE  '       ,"")          &
      ,conv2String('TNUDCENT  '       ,"")          &
      ,conv2String('TNUDLAT   '       ,"")          &
      ,conv2String('TNUDTOP   '       ,"")/),(/     &
       conv2String(oneNamelistFile%nud_type,"I14.1") &
      ,conv2String(oneNamelistFile%tnudcent,"F14.1") &
      ,conv2String(oneNamelistFile%tnudlat ,"F14.1") &
      ,conv2String(oneNamelistFile%tnudtop ,"F14.1"  )/))
    !Block10
    iErrNumber=printOneLineVars(4,(/      &
       conv2String('ZNUDTOP    '       ,"")          &
      ,conv2String('NUDLAT     '       ,"")          &
      ,conv2String('WT_NUDGE_UV'       ,"")          &
      ,conv2String('WT_NUDGE_TH'       ,"")/),(/     &
       conv2String(oneNamelistFile%znudtop    ,"F14.1") &
      ,conv2String(oneNamelistFile%nudlat     ,"I14.2") &
      ,conv2String(oneNamelistFile%wt_nudge_uv,"F14.1") &
      ,conv2String(oneNamelistFile%wt_nudge_th,"F14.1"  )/))
    !Block11
    iErrNumber=printOneLineVars(4,(/      &
       conv2String('WT_NUDGE_PI'       ,"")          &
      ,conv2String('WT_NUDGE_RT'       ,"")          &
      ,conv2String('INITIAL    '       ,"")          &
      ,conv2String('IPASTIN    '       ,"")/),(/     &
       conv2String(oneNamelistFile%wt_nudge_pi,"F14.1") &
      ,conv2String(oneNamelistFile%wt_nudge_rt,"F14.1") &
      ,conv2String(oneNamelistFile%initial    ,"I14.1") &
      ,conv2String(oneNamelistFile%ipastin    ,"I14.1"  )/))
    !Block12
    iErrNumber=printOneLineVars(4,(/      &
       conv2String('IOFLGISZ '       ,"")          &
      ,conv2String('IOFLGVAR '       ,"")          &
      ,conv2String('VWAIT1   '       ,"")          &
      ,conv2String('VWAITTOT '       ,"")/),(/     &
       conv2String(oneNamelistFile%ioflgisz,"I14.1") &
      ,conv2String(oneNamelistFile%ioflgvar,"I14.1") &
      ,conv2String(oneNamelistFile%vwait1  ,"F14.1") &
      ,conv2String(oneNamelistFile%vwaittot,"F14.1"  )/))
    !Block13
    iErrNumber=printOneLineVars(4,(/      &
       conv2String('TIMSTR  '       ,"")          &
      ,conv2String('KWRITE  '       ,"")          &
      ,conv2String('FRQPRT  '       ,"")          &
      ,conv2String('INITFLD '       ,"")/),(/     &
       conv2String(oneNamelistFile%timstr ,"F14.1") &
      ,conv2String(oneNamelistFile%kwrite ,"I14.1") &
      ,conv2String(oneNamelistFile%frqprt ,"F14.1") &
      ,conv2String(oneNamelistFile%initfld,"I14.1"  )/))
    !Block14
    iErrNumber=printOneLineVars(4,(/      &
       conv2String('TEB_SPM    ',"")          &
      ,conv2String('ICORFLG    ',"")          &
      ,conv2String('VVELDAMP   ',"")          &
      ,conv2String('IEXEV      ',"")/),(/     &
       conv2String(oneNamelistFile%teb_spm ,"I14.1") &
      ,conv2String(oneNamelistFile%icorflg ,"I14.1") &
      ,conv2String(oneNamelistFile%vveldamp,"I14.1") &
      ,conv2String(oneNamelistFile%iexev   ,"I14.1"  )/))
    !Block15
    iErrNumber=printOneLineVars(4,(/      &
       conv2String('IMASSFLX '       ,"")          &
      ,conv2String('IBND     '       ,"")          &
      ,conv2String('JBND     '       ,"")          &
      ,conv2String('CPHAS    '       ,"")/),(/     &
       conv2String(oneNamelistFile%imassflx,"I14.1") &
      ,conv2String(oneNamelistFile%ibnd    ,"I14.2") &
      ,conv2String(oneNamelistFile%jbnd    ,"I14.2") &
      ,conv2String(oneNamelistFile%cphas   ,"F14.1"  )/))
    !Block16
    iErrNumber=printOneLineVars(4,(/      &
       conv2String('LSFLG   '       ,"")          &
      ,conv2String('NFPT    '       ,"")          &
      ,conv2String('DISTIM  '       ,"")          &
      ,conv2String('ISWRTYP '       ,"")/),(/     &
       conv2String(oneNamelistFile%lsflg  ,"I14.1") &
      ,conv2String(oneNamelistFile%nfpt   ,"I14.1") &
      ,conv2String(oneNamelistFile%distim ,"F14.1") &
      ,conv2String(oneNamelistFile%iswrtyp,"I14.1"  )/))
    !Block17
    iErrNumber=printOneLineVars(4,(/      &
       conv2String('ILWRTYP     ',"")          &
      ,conv2String('LONRAD      ',"")          &
      ,conv2String('RADFRQ      ',"")          &
      ,conv2String('CLOSURE_TYPE',"")/),(/     &
       conv2String(oneNamelistFile%ilwrtyp     ,"I14.1") &
      ,conv2String(oneNamelistFile%lonrad      ,"I14.1") &
      ,conv2String(oneNamelistFile%radfrq      ,"F14.1") &
      ,conv2String(oneNamelistFile%closure_type,"A14"  )/))
    !Block18
    iErrNumber=printOneLineVars(4,(/      &
       conv2String('G3D_SPREAD '       ,"")          &
      ,conv2String('CONFRQ     '       ,"")          &
      ,conv2String('SHCUFRQ    '       ,"")          &
      ,conv2String('WCLDBS     '       ,"")/),(/     &
       conv2String(oneNamelistFile%g3d_spread,"I14.1") &
      ,conv2String(oneNamelistFile%confrq    ,"F14.1") &
      ,conv2String(oneNamelistFile%shcufrq   ,"F14.1") &
      ,conv2String(oneNamelistFile%wcldbs    ,"F14.1"  )/))
    !Block19
    iErrNumber=printOneLineVars(4,(/      &
       conv2String('Z0FACT     '       ,"")          &
      ,conv2String('MCPHYS_TYPE'       ,"")          &
      ,conv2String('M  LEVEL   '       ,"")          &
      ,conv2String('IRIME      '       ,"")/),(/     &
       conv2String(oneNamelistFile%z0fact     ,"F14.1") &
      ,conv2String(oneNamelistFile%mcphys_type,"I14.1") &
      ,conv2String(oneNamelistFile%level      ,"I14.1") &
      ,conv2String(oneNamelistFile%irime      ,"I14.1")/))
    !Block20
    iErrNumber=printOneLineVars(4,(/      &
       conv2String('IPLAWS'       ,"")          &
      ,conv2String('ICLOUD'       ,"")          &
      ,conv2String('IDRIZ '       ,"")          &
      ,conv2String('IRAIN '       ,"")/),(/     &
       conv2String(oneNamelistFile%iplaws,"I14.1") &
      ,conv2String(oneNamelistFile%icloud,"I14.1") &
      ,conv2String(oneNamelistFile%idriz ,"I14.1") &
      ,conv2String(oneNamelistFile%irain ,"I14.1"  )/))
    !Block21
    iErrNumber=printOneLineVars(4,(/      &
       conv2String('IPRIS '       ,"")          &
      ,conv2String('ISNOW '       ,"")          &
      ,conv2String('IAGGR '       ,"")          &
      ,conv2String('IGRAUP'       ,"")/),(/     &
       conv2String(oneNamelistFile%ipris ,"I14.1") &
      ,conv2String(oneNamelistFile%isnow ,"I14.1") &
      ,conv2String(oneNamelistFile%iaggr ,"I14.1") &
      ,conv2String(oneNamelistFile%igraup,"I14.1"  )/))
    !Block22
    iErrNumber=printOneLineVars(4,(/      &
       conv2String('IHAIL'       ,"")          &
      ,conv2String('CPARM'       ,"")          &
      ,conv2String('RPARM'       ,"")          &
      ,conv2String('PPARM'       ,"")/),(/     &
       conv2String(oneNamelistFile%ihail,"I14.1") &
      ,conv2String(oneNamelistFile%cparm,"E14.6") &
      ,conv2String(oneNamelistFile%rparm,"E14.6") &
      ,conv2String(oneNamelistFile%pparm,"E14.6")/))
    !Block23
    iErrNumber=printOneLineVars(4,(/      &
       conv2String('SPARM'       ,"")          &
      ,conv2String('APARM'       ,"")          &
      ,conv2String('GPARM'       ,"")          &
      ,conv2String('HPARM'       ,"")/),(/     &
       conv2String(oneNamelistFile%sparm,"E14.6") &
      ,conv2String(oneNamelistFile%aparm,"E14.6") &
      ,conv2String(oneNamelistFile%gparm,"E14.6") &
      ,conv2String(oneNamelistFile%hparm,"E14.6")/))
    !Block24
    iErrNumber=printOneLineVars(4,(/      &
       conv2String('DPARM   '       ,"")          &
      ,conv2String('EPSIL   '       ,"")          &
      ,conv2String('MKCOLTAB'       ,"")          &
      ,conv2String('IHORGRAD'       ,"")/),(/     &
       conv2String(oneNamelistFile%dparm   ,"E14.6") &
      ,conv2String(oneNamelistFile%epsil   ,"E14.6") &
      ,conv2String(oneNamelistFile%mkcoltab,"I14.1") &
      ,conv2String(oneNamelistFile%ihorgrad,"I14.1")/))
    !Block25
    iErrNumber=printOneLineVars(4,(/      &
       conv2String('DPARM   '       ,"")          &
      ,conv2String('EPSIL   '       ,"")          &
      ,conv2String('MKCOLTAB'       ,"")          &
      ,conv2String('IHORGRAD'       ,"")/),(/     &
       conv2String(oneNamelistFile%dparm   ,"F14.6") &
      ,conv2String(oneNamelistFile%epsil   ,"F14.6") &
      ,conv2String(oneNamelistFile%mkcoltab,"I14.1") &
      ,conv2String(oneNamelistFile%ihorgrad,"I14.1")/))
    !Block26
    iErrNumber=printOneLineVars(4,(/      &
       conv2String('NPATCH     '       ,"")          &
      ,conv2String('NVEGPAT    '       ,"")          &
      ,conv2String('ISFCL      '       ,"")          &
      ,conv2String('ISFCL_OCEAN'       ,"")/),(/     &
       conv2String(oneNamelistFile%npatch     ,"I14.1") &
      ,conv2String(oneNamelistFile%nvegpat    ,"I14.1") &
      ,conv2String(oneNamelistFile%isfcl      ,"I14.1") &
      ,conv2String(oneNamelistFile%isfcl_ocean,"I14.1")/))
    !Block27
    iErrNumber=printOneLineVars(4,(/      &
       conv2String('NVGCON '       ,"")          &
      ,conv2String('PCTLCON'       ,"")          &
      ,conv2String('NSLCON '       ,"")          &
      ,conv2String('ZROUGH '       ,"")/),(/     &
       conv2String(oneNamelistFile%iplaws ,"I14.1") &
      ,conv2String(oneNamelistFile%pctlcon,"F14.1") &
      ,conv2String(oneNamelistFile%nslcon ,"I14.1") &
      ,conv2String(oneNamelistFile%zrough ,"F14.2")/))
    !Block28
    iErrNumber=printOneLineVars(4,(/      &
       conv2String('ALBEDO'       ,"")          &
      ,conv2String('DTHCON'       ,"")          &
      ,conv2String('DRTCON'       ,"")          &
      ,conv2String('IPSFLG'       ,"")/),(/     &
       conv2String(oneNamelistFile%albedo,"F14.2") &
      ,conv2String(oneNamelistFile%dthcon,"F14.2") &
      ,conv2String(oneNamelistFile%drtcon,"F14.2") &
      ,conv2String(oneNamelistFile%ipsflg,"I14.1")/))
    !Block29
    iErrNumber=printOneLineVars(4,(/      &
       conv2String('ITSFLG  '       ,"")          &
      ,conv2String('IRTSFLG '       ,"")          &
      ,conv2String('IUSFLG  '       ,"")          &
      ,conv2String('WINDFARM'       ,"")/),(/     &
       conv2String(oneNamelistFile%itsflg  ,"I14.1") &
      ,conv2String(oneNamelistFile%irtsflg ,"I14.1") &
      ,conv2String(oneNamelistFile%iusflg  ,"I14.1") &
      ,conv2String(oneNamelistFile%windfarm,"I14.1")/))
    !Block30
    iErrNumber=printOneLineVars(4,(/      &
       conv2String('DAMMODULE'       ,"")          &
      ,conv2String('FRQPRECIP'       ,"")          &
      ,conv2String('EVALUATE '       ,"")          &
      ,conv2String('CHEMISTRY'       ,"")/),(/     &
       conv2String(oneNamelistFile%dammodule,"I14.1") &
      ,conv2String(oneNamelistFile%frqPrecip,"F14.4") &
      ,conv2String(oneNamelistFile%evaluate ,"I14.1") &
      ,conv2String(oneNamelistFile%chemistry,"I14.1")/))
    !Block31
    iErrNumber=printOneLineVars(4,(/      &
       conv2String('SPLIT_METHOD'       ,"")          &
      ,conv2String('CHEM_TIMESTEP'       ,"")          &
      ,conv2String('CHEMISTRY_AQ '       ,"")          &
      ,conv2String('CHEM_ASSIM'       ,"")/),(/     &
       conv2String(oneNamelistFile%split_method,"A14") &
      ,conv2String(oneNamelistFile%chem_timestep,"F14.4") &
      ,conv2String(oneNamelistFile%chemistry_aq ,"I14.1") &
      ,conv2String(oneNamelistFile%chem_assim,"I14.1")/))
    !Block32
    iErrNumber=printOneLineVars(4,(/      &
       conv2String('RECYCLE_TRACERS'       ,"")          &
      ,conv2String('DEF_PROC_SRC'       ,"")          &
      ,conv2String('NA_EXTRA2D '       ,"")          &
      ,conv2String('NA_EXTRA3D'       ,"")/),(/     &
       conv2String(oneNamelistFile%recycle_tracers,"I14.1") &
      ,conv2String(oneNamelistFile%def_proc_src,"A14") &
      ,conv2String(oneNamelistFile%na_extra2d ,"I14.4") &
      ,conv2String(oneNamelistFile%na_extra3d,"I14.4")/))
    !Block33
    iErrNumber=printOneLineVars(4,(/      &
       conv2String('PLUMERISE'       ,"")          &
      ,conv2String('PRFRQ'       ,"")          &
      ,conv2String('VOLCANOES '       ,"")          &
      ,conv2String('AEROSOL'       ,"")/),(/     &
       conv2String(oneNamelistFile%plumerise,"I14.1") &
      ,conv2String(oneNamelistFile%prfrq,"F14.4") &
      ,conv2String(oneNamelistFile%volcanoes ,"I14.1") &
      ,conv2String(oneNamelistFile%aerosol,"I14.1")/))
    !Block34
    iErrNumber=printOneLineVars(4,(/      &
       conv2String('AER_TIMESTEP'       ,"")          &
      ,conv2String('CCATT'       ,"")          &
      ,conv2String('ISZSTAGE '       ,"")          &
      ,conv2String('IVRSTAGE'       ,"")/),(/     &
       conv2String(oneNamelistFile%aer_timestep,"F14.4") &
      ,conv2String(oneNamelistFile%catt,"I14.1") &
      ,conv2String(oneNamelistFile%iszstage ,"I14.1") &
      ,conv2String(oneNamelistFile%ivrstage,"I14.1")/))
    !Block35
    iErrNumber=printOneLineVars(4,(/      &
       conv2String('ISAN_INC    '       ,"")          &
      ,conv2String('GUESS1ST    '       ,"")          &
      ,conv2String('I1ST_FLG    '       ,"")          &
      ,conv2String('IUPA_FLG    '       ,"")/),(/     &
       conv2String(oneNamelistFile%isan_inc       ,"I14.6") &
      ,conv2String(oneNamelistFile%guess1st,"A14") &
      ,conv2String(oneNamelistFile%i1st_flg,"I14.1") &
      ,conv2String(oneNamelistFile%iupa_flg,"I14.1"  )/))
    !Block36
    iErrNumber=printOneLineVars(4,(/      &
       conv2String('ISFC_FLG    '       ,"")          &
      ,conv2String('IUPDNDVI    '       ,"")          &
      ,conv2String('IUPDSST     '       ,"")          &
      ,conv2String('NOT USED    '       ,"")/),(/     &
       conv2String(oneNamelistFile%isfc_flg,"I14.1") &
      ,conv2String(oneNamelistFile%iupdndvi,"I14.1") &
      ,conv2String(oneNamelistFile%iupdsst,"I14.1") &
      ,conv2String('',"A14"  )/))
    !Block37
    iErrNumber=printOneLineVars(4,(/      &
       conv2String('NISN        '       ,"")          &
      ,conv2String('NIGRIDS     '       ,"")          &
      ,conv2String('TOPSIGZ     '       ,"")          &
      ,conv2String('HYBBOT      '       ,"")/),(/     &
       conv2String(oneNamelistFile%nisn   ,"I14.4") &
      ,conv2String(oneNamelistFile%nigrids,"I14.4") &
      ,conv2String(oneNamelistFile%topsigz,"F14.6") &
      ,conv2String(oneNamelistFile%hybbot ,"F14.6"  )/))
    !Block38
    iErrNumber=printOneLineVars(4,(/      &
       conv2String('HYBTOP      '       ,"")          &
      ,conv2String('SFCINF      '       ,"")          &
      ,conv2String('SIGZWT      '       ,"")          &
      ,conv2String('NFEEDVAR    '       ,"")/),(/     &
       conv2String(oneNamelistFile%hybtop   ,"F14.6") &
      ,conv2String(oneNamelistFile%sfcinf   ,"F14.6") &
      ,conv2String(oneNamelistFile%sigzwt   ,"F14.6") &
      ,conv2String(oneNamelistFile%nfeedvar ,"I14.1"  )/))
    !Block40
    iErrNumber=printOneLineVars(4,(/      &
       conv2String('MAXSTA      '       ,"")          &
      ,conv2String('MAXSFC      '       ,"")          &
      ,conv2String('NOTSTA      '       ,"")          &
      ,conv2String('NOT USED     '       ,"")/),(/     &
       conv2String(oneNamelistFile%maxsta   ,"I14.6") &
      ,conv2String(oneNamelistFile%maxsfc   ,"I14.6") &
      ,conv2String(oneNamelistFile%notsta   ,"I14.1") &
      ,conv2String('' ,"A14"  )/))
    !Block41
    iErrNumber=printOneLineVars(4,(/      &
       conv2String('IOBSWIN      '       ,"")          &
      ,conv2String('STASEP      '       ,"")          &
      ,conv2String('IGRIDFL      '       ,"")          &
      ,conv2String('NOT USED'       ,"")/),(/     &
       conv2String(oneNamelistFile%iobswin  ,"I14.6") &
      ,conv2String(oneNamelistFile%stasep   ,"F14.2") &
      ,conv2String(oneNamelistFile%igridfl  ,"I14.1") &
      ,conv2String(''   ,"A14"  )/))
    !Block42
    iErrNumber=printOneLineVars(4,(/      &
       conv2String('GRIDWT(1)   '       ,"")          &
      ,conv2String('GRIDWT(2)   '       ,"")          &
      ,conv2String('WVLNTH(1)      '       ,"")          &
      ,conv2String('WVLNTH(2)       '       ,"")/),(/     &
       conv2String(oneNamelistFile%gridwt(1)  ,"F14.6") &
      ,conv2String(oneNamelistFile%gridwt(2)  ,"F14.6") &
      ,conv2String(oneNamelistFile%wvlnth(1)  ,"F14.6") &
      ,conv2String(oneNamelistFile%wvlnth(2)  ,"F14.6"  )/))
    !Block43
    iErrNumber=printOneLineVars(4,(/      &
       conv2String('SWVLNTH(1)   '       ,"")          &
      ,conv2String('SWVLNTH(2)   '       ,"")          &
      ,conv2String('RESPON(1)      '       ,"")          &
      ,conv2String('RESPON(2)       '       ,"")/),(/     &
       conv2String(oneNamelistFile%swvlnth(1)  ,"F14.6") &
      ,conv2String(oneNamelistFile%swvlnth(2)  ,"F14.6") &
      ,conv2String(oneNamelistFile%respon(1)  ,"F14.6") &
      ,conv2String(oneNamelistFile%respon(2)  ,"F14.6"  )/))
    !Block44
    iErrNumber=printOneLineVars(4,(/      &
       conv2String('GOBSEP   '       ,"")          &
      ,conv2String('GOBRAD   '       ,"")          &
      ,conv2String('ICFILETYPE    '       ,"")          &
      ,conv2String('WIND_U_VARNAME'       ,"")/),(/     &
       conv2String(oneNamelistFile%gobsep  ,"F14.6") &
      ,conv2String(oneNamelistFile%gobrad  ,"F14.6") &
      ,conv2String(oneNamelistFile%icFileType  ,"I14.1") &
      ,conv2String(oneNamelistFile%wind_u_varname  ,"A14"  )/))
    !Block45
    iErrNumber=printOneLineVars(4,(/      &
       conv2String('WIND_V_VARNAME'       ,"")          &
      ,conv2String('TEMP_VARNAME  '       ,"")          &
      ,conv2String('GEO_VARNAME   '       ,"")          &
      ,conv2String('UR_VARNAME    '       ,"")/),(/     &
       conv2String(oneNamelistFile%wind_v_varname  ,"A14") &
      ,conv2String(oneNamelistFile%temperature_varname  ,"A14") &
      ,conv2String(oneNamelistFile%geo_varname  ,"A14") &
      ,conv2String(oneNamelistFile%ur_varname  ,"A14"  )/))
    !Block46
    iErrNumber=printOneLineVars(4,(/      &
       conv2String('DLIMIT U'       ,"")          &
      ,conv2String('DLIMIT V'       ,"")          &
      ,conv2String('ULIMIT U'       ,"")          &
      ,conv2String('ULIMIT V'       ,"")/),(/     &
       conv2String(oneNamelistFile%dlimit(1)  ,"F14.4") &
      ,conv2String(oneNamelistFile%dlimit(2)  ,"F14.4") &
      ,conv2String(oneNamelistFile%ulimit(1)  ,"F14.4") &
      ,conv2String(oneNamelistFile%ulimit(2)  ,"F14.4"  )/))
    !Block47
    iErrNumber=printOneLineVars(4,(/      &
       conv2String('DLIMIT T'       ,"")          &
      ,conv2String('DLIMIT Z'       ,"")          &
      ,conv2String('ULIMIT T'       ,"")          &
      ,conv2String('ULIMIT Z'       ,"")/),(/     &
       conv2String(oneNamelistFile%dlimit(3)  ,"F14.4") &
      ,conv2String(oneNamelistFile%dlimit(4)  ,"F14.4") &
      ,conv2String(oneNamelistFile%ulimit(3)  ,"F14.4") &
      ,conv2String(oneNamelistFile%ulimit(4)  ,"F14.4"  )/))
    !Block48
    iErrNumber=printOneLineVars(4,(/      &
       conv2String('DLIMIT R'       ,"")          &
      ,conv2String('ULIMIT R'       ,"")          &
      ,conv2String('SCALE_FACTOR U'       ,"")          &
      ,conv2String('SCALE_FACTOR V'       ,"")/),(/     &
       conv2String(oneNamelistFile%dlimit(5)  ,"F14.4") &
      ,conv2String(oneNamelistFile%ulimit(5)  ,"F14.4") &
      ,conv2String(oneNamelistFile%scale_factor(1)  ,"F14.4") &
      ,conv2String(oneNamelistFile%scale_factor(2)  ,"F14.4"  )/))
    !Block49
    iErrNumber=printOneLineVars(4,(/      &
       conv2String('SCALE_FACTOR T'       ,"")          &
      ,conv2String('SCALE_FACTOR Z'       ,"")          &
      ,conv2String('SCALE_FACTOR R'       ,"")          &
      ,conv2String('NPLT          '       ,"")/),(/     &
       conv2String(oneNamelistFile%scale_factor(3)  ,"F14.4") &
      ,conv2String(oneNamelistFile%scale_factor(4)  ,"F14.4") &
      ,conv2String(oneNamelistFile%scale_factor(5)  ,"F14.4") &
      ,conv2String(oneNamelistFile%nplt  ,"I14.2"  )/))
    !Block50
    iErrNumber=printOneLineVars(4,(/      &
       conv2String('INITIAL_LATITUDE'       ,"")          &
      ,conv2String('FINAL_LATITUDE  '       ,"")          &
      ,conv2String('INITIAL_LONGITUDE   '       ,"")          &
      ,conv2String('FINAL_LONGITUDE    '       ,"")/),(/     &
       conv2String(oneNamelistFile%initial_latitude  ,"F14.4") &
      ,conv2String(oneNamelistFile%final_latitude  ,"F14.4") &
      ,conv2String(oneNamelistFile%initial_longitude  ,"F14.4") &
      ,conv2String(oneNamelistFile%final_longitude  ,"F14.4"  )/))
    !Block51
    iErrNumber=printOneLineVars(4,(/      &
       conv2String('Z_MAX_LEVEL'       ,"")          &
      ,conv2String('CCGRADSWRITE  '       ,"")          &
      ,conv2String('NVP   '       ,"")          &
      ,conv2String('ANL2GRA    '       ,"")/),(/     &
       conv2String(oneNamelistFile%z_max_level  ,"I14.3") &
      ,conv2String(oneNamelistFile%ccGradsWrite  ,"I14.1") &
      ,conv2String(oneNamelistFile%nvp  ,"I14.3") &
      ,conv2String(oneNamelistFile%anl2gra  ,"A14"  )/))
    !Block52
    iErrNumber=printOneLineVars(4,(/      &
       conv2String('PROJ'       ,"")          &
      ,conv2String('MEAN_TYPE  '       ,"")          &
      ,conv2String('IPRESSLEV   '       ,"")          &
      ,conv2String('INPLEVS    '       ,"")/),(/     &
       conv2String(oneNamelistFile%proj  ,"A14") &
      ,conv2String(oneNamelistFile%mean_type  ,"A14") &
      ,conv2String(oneNamelistFile%ipresslev  ,"I14.1") &
      ,conv2String(oneNamelistFile%inplevs  ,"I14.3"  )/))
    !Block53
    iErrNumber=printOneLineVars(4,(/      &
       conv2String('ASCII_DATA'       ,"")          &
      ,conv2String('SITE_LAT  '       ,"")          &
      ,conv2String('SITE_LON   '       ,"")          &
      ,conv2String('INPLEVS    '       ,"")/),(/     &
       conv2String(oneNamelistFile%ascii_data  ,"A14") &
      ,conv2String(oneNamelistFile%site_lat  ,"F14.4") &
      ,conv2String(oneNamelistFile%site_lon  ,"F14.4") &
      ,conv2String(oneNamelistFile%inplevs  ,"I14.3"  )/))
    !Block54
    iErrNumber=printOneLineVars(4,(/      &
       conv2String('APPLYDIGITALFILTER'       ,"")          &
      ,conv2String('DFILTERIMEWINDOW'       ,"")          &
      ,conv2String('APPLYMETEOGRAM   '       ,"")          &
      ,conv2String('METEOGRAMFREQ    '       ,"")/),(/     &
       conv2String(oneNamelistFile%applyDigitalFilter  ,"A14") &
      ,conv2String(oneNamelistFile%digitalFilterTimeWindow  ,"F14.1") &
      ,conv2String(oneNamelistFile%applyMeteogram  ,"A14") &
      ,conv2String(oneNamelistFile%meteogramFreq  ,"F14.4"  )/))
    iErrNumber=printVarTail(4)
    if(oneNamelistFile%deltaz==0) then 
      iErrNumber=dumpMessage(c_tty,c_yes,'','',c_noError,'ZZ       :',oneNamelistFile%zz,"F7.1")
    elseif(oneNamelistFile%fixlevels>0) then 
      iErrNumber=dumpMessage(c_tty,c_yes,'','',c_noError,'ZZ FIX   :',oneNamelistFile%zz(1:oneNamelistFile%fixlevels),"F7.1")
    endif
    iErrNumber=dumpMessage(c_tty,c_yes,'','',c_noError,'SLZ       :',oneNamelistFile%slz,"F8.2")
    iErrNumber=dumpMessage(c_tty,c_yes,'','',c_noError,'SLMSTR    :',oneNamelistFile%slmstr,"F8.2")
    iErrNumber=dumpMessage(c_tty,c_yes,'','',c_noError,'STGOFF    :',oneNamelistFile%stgoff,"F8.2")
    iErrNumber=dumpMessage(c_tty,c_yes,'','',c_noError,'LEVTH     :',oneNamelistFile%stgoff,"F8.1")
    !iErrNumber=dumpMessage(c_tty,c_yes,'','',c_noError,'IPLFLD    :',oneNamelistFile%iplfld(1:10))
    iErrNumber=dumpMessage(c_tty,c_yes,'','',c_noError,'IXSCTN    :',oneNamelistFile%ixsctn(1:10),"I1")
    iErrNumber=dumpMessage(c_tty,c_yes,'','',c_noError,'ISBVAL    :',oneNamelistFile%isbval,"I1")
    iErrNumber=dumpMessage(c_tty,c_yes,'','',c_noError,'DIUR_CYCLE:',oneNamelistFile%diur_cycle(1:nsrc),"I1")
    iErrNumber=dumpMessage(c_tty,c_yes,'','',c_noError,'GNU       :',oneNamelistFile%gnu,"F4.1")

    write(*,fmt='(A)') 
    write(*,"(a)") "+-------------------------------------------- Files from RAMSIN ---------------------------------------------------"
    write(*,fmt='(A)') 

    iErrNumber=printFileHeader()
    iErrNumber=printOneFile("USDATA_IN ",trim(oneNamelistFile%usdata_in)   ," I ")
    iErrNumber=printOneFile("USMODEL_IN",trim(oneNamelistFile%usmodel_in)  ," I ")
    iErrNumber=printOneFile('HFILOUT    ',trim(oneNamelistFile%hfilout)    ," O " )
    iErrNumber=printOneFile('AFILOUT    ',trim(oneNamelistFile%afilout)    ," O " )
    iErrNumber=printOneFile('FILENAMEIAU',trim(oneNamelistFile%fileNameIAU),"I/O" )
    iErrNumber=printOneFile('VARFPFX    ',trim(oneNamelistFile%varfpfx)    ," O " ) 
    iErrNumber=printOneFile('VARPFX     ',trim(oneNamelistFile%varpfx)     ," I " )
    iErrNumber=printOneFile('NUD_HFILE  ',trim(oneNamelistFile%nud_hfile)  ," I " ) 
    iErrNumber=printOneFile('HFILIN     ',trim(oneNamelistFile%hfilin)     ," I " ) 
    iErrNumber=printOneFile('PASTFN     ',trim(oneNamelistFile%pastfn)     ," I " ) 
    iErrNumber=printOneFile('RADDATFN   ',trim(oneNamelistFile%raddatfn)   ," I " )
    iErrNumber=printOneFile('JULESIN    ',trim(oneNamelistFile%julesin)    ," I " )
    iErrNumber=printOneFile('COLTABFN   ',trim(oneNamelistFile%coltabfn)   ," I " )       
    iErrNumber=printOneFile('SOIL_MOIST ',trim(oneNamelistFile%soil_moist) ," I " )
    iErrNumber=printOneFile('SOIL_M_FAIL',trim(oneNamelistFile%soil_moist_fail) ," I "  )
    iErrNumber=printOneFile('SRCMAPFN   ',trim(oneNamelistFile%srcmapfn)   ," I "       )
    iErrNumber=printOneFile('MAPAOTFILE    ',trim(oneNamelistFile%MapAOTFile),"I")
    iErrNumber=printOneFile('SSTFPFX       ',trim(oneNamelistFile%sstfpfx),"O")
    iErrNumber=printOneFile('NDVIFPFX      ',trim(oneNamelistFile%ndvifpfx),"O")
    iErrNumber=printOneFile('WFFILE        ',trim(oneNamelistFile%wfFile),"O")
    iErrNumber=printOneFile('DAMOUTPREFIX  ',trim(oneNamelistFile%damOutPrefix),"O")
    iErrNumber=printOneFile('EVALUATEPREFIX',trim(oneNamelistFile%evaluatePrefix),"O")
    do k=1,ngrids
      write(ck,fmt='(I1.1)') k
      iErrNumber=printOneFile('ITOPTFN('//ck//')',trim(oneNamelistFile%itoptfn(k)),"I")
      iErrNumber=printOneFile('ISSTFN ('//ck//')',trim(oneNamelistFile%isstfn(k)),"I")
      iErrNumber=printOneFile('IVEGTFN('//ck//')',trim(oneNamelistFile%ivegtfn(k)),"I")
      iErrNumber=printOneFile('ISOILFN('//ck//')',trim(oneNamelistFile%isoilfn(k)),"I")
      iErrNumber=printOneFile('NDVIFN ('//ck//')',trim(oneNamelistFile%ndvifn(k)),"I")
    enddo
    if (len_trim(oneNamelistFile%domain_fname) > 0) then
      iErrNumber=printOneFile('DOMAIN_FNAME ('//ck//')',trim(oneNamelistFile%domain_fname),"I")
    else
      iErrNumber=printOneFile('DOMAIN_FNAME ('//ck//')','Using Domain decomposition',"*")
    endif
    iErrNumber=printOneFile('ICPREFIX',trim(oneNamelistFile%icPrefix),"I")
    iErrNumber=printOneFile('ICGRADSPREFIX',trim(oneNamelistFile%icGradsPrefix),"O")
    iErrNumber=printOneFile('GPREFIX',trim(oneNamelistFile%gprefix),"O")
    iErrNumber=printOneFile('IAPR   ',trim(oneNamelistFile%iapr),"I")
    iErrNumber=printOneFile('IARAWI ',trim(oneNamelistFile%iarawi),"I")
    iErrNumber=printOneFile('IASRFCE',trim(oneNamelistFile%iasrfce),"I")
    iErrNumber=printOneFile('CSVFILE     ',trim(oneNamelistFile%csvFile),"I")
    iErrNumber=printOneFile('METEOGRAMMAP',trim(oneNamelistFile%meteogramMap),"I")
    iErrNumber=printOneFile('METEOGRAMDIR',trim(oneNamelistFile%meteogramDir),"I")
    iErrNumber=printFileTail()


    if(oneNamelistFile%ipos>0) then
      write(*,fmt='(A)') 
      write(*,"(a)") "+---------------------------------------Variables Selected to Output-----------------------------------------------"
      write(*,fmt='(A)') 
      iErrNumber=csvHeader(oneNamelistFile%csvFile)
      do k=1,oneNamelistFile%nvp
          iErrNumber=printOnecsv(oneNamelistFile%vp(k))
      enddo 
      iErrNumber=csvTail()
      iErrNumber=dumpMessage(c_tty,c_yes,'','',c_noError,'IPLEVS:' &
                             ,oneNamelistFile%iplevs(1:oneNamelistFile%inplevs),"I4")
    
    endif
    
    !iErrNumber=dumpMessage(c_tty,c_yes,'RAMSIN','Namelist',c_notice,'HS         :',oneNamelistFile%hs,"F4.1")
    !iErrNumber=dumpMessage(c_tty,c_yes,'RAMSIN','Namelist',c_notice,'PS         :',oneNamelistFile%ps,"F8.2")
    !iErrNumber=dumpMessage(c_tty,c_yes,'RAMSIN','Namelist',c_notice,'TS         :',oneNamelistFile%ts,"F8.2") 
    !iErrNumber=dumpMessage(c_tty,c_yes,'RAMSIN','Namelist',c_notice,'RTS        :',oneNamelistFile%rts,"F8.2") 
    !iErrNumber=dumpMessage(c_tty,c_yes,'RAMSIN','Namelist',c_notice,'US         :',oneNamelistFile%us,"F8.2") 
    !iErrNumber=dumpMessage(c_tty,c_yes,'RAMSIN','Namelist',c_notice,'VS         :',oneNamelistFile%vs,"F8.2")
    write(*,*)


  end subroutine DumpNamelistFile


end module ModNamelistFile
