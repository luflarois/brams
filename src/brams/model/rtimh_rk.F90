
!############################# Change Log ##################################
! 5.2.0
!
!MB: new file for the Runge-Kutta dynamical core
!    (after Wicker, Skamarock, 2002, MWR)
!    Saulo Freitas (INPE), Michael Baldauf (DWD)
!
!############################# Change Log ##################################

module ModTimestep_RK

  !MB: for testing only
  logical :: flag_mb_adv_test=.false.
  character(len=*), parameter :: h="**(rtimh_rk)**"
logical,parameter :: stepDebug=.true.

contains
subroutine timestep_rk(OneGrid,oneNamelistFile)

  use ModMessageSet, only: &
       PostRecvSendMsgs, &
       WaitRecvMsgs

  use ModAcoust, only:         &
       acoustic_new,            &
       init_div_damping_coeff,  &
       deallocate_alpha_div,    &
       apply_div_damping
!
  use ModGrid, only: &
       Grid

  use mem_basic, only: &
       basic_g  ! INTENT(INOUT)

  use node_mod, only: &
       mzp, mxp, myp,  & ! INTENT(IN)
       ia, iz, ja, jz, & ! INTENT(IN)
       i0, j0,         & ! INTENT(IN)
       izu, jzv,       & ! INTENT(IN)
       mynum,          & ! INTENT(IN)
       ibcon,          & ! INTENT(IN)
       nmachs, &
       mchnum,        &
       master_num

  use mem_cuparm, only: &
       NNQPARM, & ! INTENT(IN)
       IF_CUINV   ! INTENT(IN)

  use mem_varinit, only: &
       NUD_TYPE   &! INTENT(IN)
      ,varinit_g  ! INTENT(IN)
      
  use mem_turb, only: &
       IF_URBAN_CANOPY, & ! INTENT(IN)
       ihorgrad           ! INTENT(IN)

  use mem_oda,   only: &
       if_oda ! INTENT(IN)

  use micphys,   only: &
       mcphys_type,  &! INTENT(IN)
       level          ! INTENT(IN)

  use mem_grid, only: &
       ngrids,     & ! INTENT(IN)
       ngrid,      & ! INTENT(IN)
       npatch,     & ! INTENT(IN)
       time,       & ! INTENT(IN)
       dts,        &
       dtlong,     & ! INTENT(IN)
       dtlongn,    & ! INTENT(IN)
       ideltat,    &
       nnacoust,   &
       iyear1,     & ! INTENT(IN)
       imonth1,    & ! INTENT(IN)
       idate1,     & ! INTENT(IN)
       grid_g,     & ! INTENT(INOUT)
       nxtnest,    & ! INTENT(IN)
       if_adap,    & ! INTENT(IN)
       dtlt,       & ! INTENT(IN)
       istp,       & ! INTENT(IN)
       jdim,       & ! INTENT(IN)
       nzp,        & ! INTENT(IN)
       f_thermo_e, & ! INTENT(IN)
       f_thermo_w, & ! INTENT(IN)
       f_thermo_s, & ! INTENT(IN)
       f_thermo_n, &   ! INTENT(IN)
       zt,         &
       zm,         &
       dzt,        &
       nzpmax,     &
       itime1,     &
       vveldamp    ! INTENT(IN)

  use shcu_vars_const, only: & ! For Shallow Cumulus Paramet.
       NNSHCU ! INTENT(IN)

  use mem_scalar, only: & ! For SiB
       scalar_g ! INTENT(IN)

  use mem_leaf, only: & 
       ISFCL          & ! INTENT(IN)
      ,ISFCL_OCEAN      ! INTENT(IN)

  ! TEB_SPM
  use teb_spm_start, only: &
       TEB_SPM ! INTENT(IN)
  use mem_emiss, only: &
       ichemi,         & ! INTENT(IN)
       isource           ! INTENT(IN)

  !ALF- Necessary in new advection scheme
  use advect_kit, only : &
       calc_advec          ! Subroutine

  ! For specific optimization depending the type of machine
  use machine_arq, only: &
       machine ! INTENT(IN)

  use rconstants, only:  &
       g,                & ! (IN)
       cp,               & ! (IN)
       cpor,             & ! (IN)
       p00,              & ! (IN)
       rgas,             & ! (IN)
       pi180               ! (IN)

  use ccatt_start, only: &
       ccatt               ! (IN)

  use mem_tend, only: &
       tend

  use mem_chem1, only: &
       nvert_src=>chem1_src_z_dim_g, & ! (IN)
       chem1_g,                      & ! (INOUT)
       nsrc,                         & ! (IN)
       chem1_src_g,                  & ! %sc_src(INOUT)
       chemistry,                    & ! (IN)
       split_method,                 & ! (IN)
       n_dyn_chem,                   &
       ntimes_src

  use mem_aer1, only:                  &
       aerosol,                        &! (IN)
       aer1_g,                         &! %sc_src(INOUT)
       aer2_g,                         &! %sc_src(INOUT)
       aer_nvert_src=>aer1_src_z_dim_g  ! (IN)

  use mem_plume_chem1, only:  plume_mean_g &   ! %flam_frac(IN), %fire_size(IN)
                             ,plume_fre_g
  use mem_stilt, only: &
       iexev,          &  ! (IN)
       imassflx,       &  ! (IN)
       stilt_g            ! %dnp (IN)

  use mem_radiate, only: radiate_g

  use chem_sources, only :     &
       alloc_emiss_cycle,      &  ! Subroutine
       init_actual_time_index, &  ! Subroutine
       emiss_cycle,            &  ! (INOUT)
       emiss_cycle_alloc,      &
       srcmapfn                   ! (IN)

  use ChemSourcesDriver, only:  sources_driver            ! Subroutine

  use ChemDryDepDriver , only:  drydep_driver             ! Subroutine

  use module_chemistry_driver, only: chemistry_driver ! Subroutine

  use radiation, only: radiate ! Subroutine

  use ModTimeStamp, only: SynchronizedTimeStamp, TimeStamp

  use cuparm_grell3, only: cuparm_grell3_catt !&  ! subroutine
                          !,g3d_g

  use digitalfilter, only: 	        &
                    applyDigitalFilter, & ! subroutine
                    fileNameDF,		& ! intent(inout) - file control
                    dfVars,             &
                    applyDF

  USE monotonic_adv, only:                 &
                           advmnt_driver,  &        ! subroutine
                           advmnt

  USE DriverMatrix, ONLY: MatrixDriver  !Matrix Aerosol Model


  USE rams_microphysics_2M, only: micro_2M_rams60,negadj1_2M_rams60

  use rrtm_driv, only:  rrtm_driver

  USE mem_radiate, ONLY: &
         ilwrtyp, iswrtyp

  use mem_turb, only:    &
       turb_g

  use ModComm, only: commHaloAcou

  use wind_farm, only: wind_farm_driver,windfarm

  use optical, only: &
            aodDriver

  use aerClimMod, only: &
          no_months, &
          no_src_types, &
          specieName, &
          aercam

  use ModNamelistFile, only : namelistFile


  use modIau, only:    &
     CreateIauTendency &
    ,readIauTendency   &
    ,GetIauTendency    &
    ,timeWindowIAU     &
    ,applyIAU
 
  use mod_leaf3_ocean_only, only: sfclyr_ocean_only

  use sfireMod, only: sfclyr_sfire

  use mem_sfire, only: sfire

  implicit none

  type(Grid), pointer :: OneGrid
  type(namelistFile), pointer :: oneNamelistFile

  ! execution time instrumentation
  include "constants.f90"
  include "tsNames.h"
  include "i8.h"

  LOGICAL, parameter :: flag_Coriolis_in_every_RK_step = .FALSE.

  integer :: l_rk
  real    :: rk_beta(3)
  integer :: rk_nmbr_small_timesteps(3)
  integer, parameter ::  &!rk_order = 2 ! for testing purposes (but works only with upwind1 and upwind3 advection)
                           rk_order = 3 ! for Wicker, Skamarock (2002) MWR-scheme

  integer :: n
  character(len=2) :: crk
  logical :: singleProcRun
  character(len=2) :: cnzp
  integer :: nv,itime,i,j
  character(len=256) :: julesFile

  !MB: only for testing
  !integer :: nmbr_gpts
  !real    :: pm,tm

  singleProcRun = nmachs == 1
  julesFile=oneNamelistFile%julesin

  !        +------------------------------------------------------------------+
  !        |   Timestep driver for the Runge-Kutta non-hydrostatic time-split |
  !        |      model.                                                      |
  !        +------------------------------------------------------------------+
  !
  !  Zero out all tendency arrays.
  !--------------------------------
  call TEND0()  

!------------------TMP 
!------------------TMP 
!------------------TMP 
! if(applyIAU == 1 ) then
!    if(mynum==1) print*,"timeIAU=",time,timeWindowIAU*0.5,abs ( time - dtlt - timeWindowIAU*0.5),applyIAU
!    call flush(6)
!    
!    call CreateIauTendency(ngrid, mzp*mxp*myp, mzp, mxp, myp,ia,iz,ja,jz&
!          ,varinit_g(ngrid)%varup(:,:,:),varinit_g(ngrid)%varvp(:,:,:)  &
!          ,varinit_g(ngrid)%varpp(:,:,:),varinit_g(ngrid)%vartp(:,:,:)  &
!          ,varinit_g(ngrid)%varrp(:,:,:)                                &
!	  
!          ,varinit_g(ngrid)%varuf(:,:,:),varinit_g(ngrid)%varvf(:,:,:)  &
!          ,varinit_g(ngrid)%varpf(:,:,:),varinit_g(ngrid)%vartf(:,:,:)  &
!          ,varinit_g(ngrid)%varrf(:,:,:)                                &
!
!          ,basic_g(ngrid)%up     (:,:,:)   ,basic_g(ngrid)%vp  (:,:,:)  &
!          ,basic_g(ngrid)%theta  (:,:,:)   ,basic_g(ngrid)%rtp (:,:,:)  &
!          ,basic_g(ngrid)%pp     (:,:,:)                                )
!  !RETURN
!  endif
!------------------TMP 
!------------------TMP 
!------------------TMP 
 
  ! Implements the Incremental Analysis Update procedure -
  ! phase 2: add the IAU tendencies
  !-------------------------------
  if(applyIAU == 2 ) then
    call readIauTendency(ngrid, mzp*mxp*myp, mzp, mxp, myp,ia,iz,ja,jz ,time,dtlt)
    call GetIauTendency (ngrid, mzp*mxp*myp, mzp, mxp, myp,ia,iz,ja,jz ,time,dtlt)      
  endif

  !  Thermodynamic diagnosis
  !--------------------------------
  if (mcphys_type <= 1 .and. level/=3) then
     call THERMO(mzp,mxp,myp,ia,iz,ja,jz,'SUPSAT'); !if(stepDebug) print *,mynum,THERMO'
  endif

  if (iexev == 2) &
     ! evolution of the Exner pressure: compression term
     call exevolve(mzp,mxp,myp,ngrid,ia,iz,ja,jz,izu,jzv,jdim,mynum,dtlt,'ADV')

  if (CCATT==1 .and. chemistry >= 0) call aodDriver(mzp,mxp,myp,ia,iz,ja,jz,ngrids)

  !  Radiation parameterization
  !--------------------------------
  call RADIATE(mzp,mxp,myp,ia,iz,ja,jz,mynum)

  !  Surface layer, soil and veggie model
  !----------------------------------------
  if (isfcl<=2) then
     call SFCLYR(mzp,mxp,myp,ia,iz,ja,jz,ibcon)
#ifdef JULES
  elseif (isfcl == 5) then
     if (time==0.) &
     call SFCLYR      (mzp,mxp,myp,ia,iz,ja,jz,ibcon)

     call SFCLYR_JULES(mzp,mxp,myp,ia,iz,ja,jz,jdim,julesFile,i0,j0)
     !--- this combines the JULES land + LEAF ocean models.
     if (isfcl_ocean == 1) & 
     call sfclyr_ocean_only  (mzp,mxp,myp,ia,iz,ja,jz,ibcon)
     
!if(stepDebug) print *,mynum,SFCLYR_JULES'
#endif
  endif

  !- Sea salt Aerossol inline source
     call SeaSaltDriver(ia,iz,ja,jz,ngrid,mxp,myp)

  !-- emission/deposition for CCATT chemistry models
  if (CCATT==1 .and. chemistry >= 0) then

     !- emission
     !- allocation for diurnal cycle of emission arrays   and
     !- actual_time_index on nodes
     if( (.not. emiss_cycle_alloc) .and. (chemistry >= 0) .and. &
          trim(srcmapfn) .ne.  'NONE' .and. trim(srcmapfn) .ne.  'none') then

	  call alloc_emiss_cycle(mxp,myp,ngrids,nsrc)

          call init_actual_time_index(nsrc,ntimes_src)

     endif
     !

     !plume_mean_g(:,:) instead of plume_mean_g(:,ngrid) to avoid memory errors.
     !emiss_cycle(:,:)  instead of emiss_cycle(:,ngrid)  to avoid memory errors.
     !the same for the others var
     call sources_driver(ngrid, mzp,mxp,myp,ia,iz,ja,jz,                          &
                         g,cp,cpor,p00,rgas,pi180,                                &
                         radiate_g(ngrid)%cosz,basic_g(ngrid)%theta,              &
                         basic_g(ngrid)%pp,basic_g(ngrid)%pi0,basic_g(ngrid)%rv,  &
                         basic_g(ngrid)%dn0,basic_g(ngrid)%up,basic_g(ngrid)%vp,  &
	   		 time,iyear1,imonth1,idate1,itime1,dtlt,                  &
                         grid_g(ngrid)%rtgt,grid_g(ngrid)%lpw,grid_g(ngrid)%glat, &
                         grid_g(ngrid)%glon,zt,zm,dzt,nzpmax,                     &
                         nvert_src    (:,:),                                  &
                         chem1_g      (:,:),                                  &
                         chem1_src_g  (:,:,:,:),                              &
                         aer1_g       (:,:,:),                                &
                         aer_nvert_src(:,:),                                  &
                         plume_mean_g (:,:),                                  &
                         stilt_g(ngrid)%dnp,                                  &
                         emiss_cycle  (:,:),                                  &
                         aer2_g       (:,:),                                  &
			 plume_fre_g  (:,:)                                   )


     !- call dry deposition and sedimentation routines
     call drydep_driver(mzp,mxp,myp,ia,iz,ja,jz)

     !- call Matrix Aerosol Model
     !----------------------------------------
     !if(AEROSOL==2) then
     !print*,"not doing matrix";call flush(6)
     !    CALL MatrixDriver(ia,iz,ja,jz,mzp,mxp,myp)
     !endif

  endif

  !!$  call TimeStamp(TS_PHYSICS)

  !  Send boundaries to adjoining nodes
  !-------------------------------------------
  if (nmachs > 1) then
     call PostRecvSendMsgs(OneGrid%SelectedGhostZoneSend, OneGrid%SelectedGhostZoneRecv)
  endif

  !  Coriolis terms
  !  ----------------------------------------
  if ( .NOT. flag_Coriolis_in_every_RK_step ) then
    call CORLOS(mzp,mxp,myp,i0,j0,ia,iz,ja,jz,izu,jzv, tend%ut, tend%vt)
  end if

  !  Cumulus parameterization version 1
  !----------------------------------------
  if (NNQPARM(ngrid)==1 .or. IF_CUINV==1) then
     call cuparm()
  end if

  !  Urban canopy parameterization
  !----------------------------------------
  if (IF_URBAN_CANOPY==1) call urban_canopy()
  !if(stepDebug) print *,mynum,urban_canopy'
  !  Analysis nudging and boundary condition
  !------------------------------------------
  !!$  call TimeStamp(TS_PHYSICS)
  
  if (NUD_TYPE>0) call DATASSIM()
  !if(stepDebug) print *,mynum,DATASSIM'
  !!$  call TimeStamp(TS_DINAMICS)

  !  Observation data assimilation
  !----------------------------------------
  if (IF_ODA==1) call oda_nudge()
  !if(stepDebug) print *,mynum,oda_nudge'
  !!$  call TimeStamp(TS_PHYSICS)

  !  Nested grid boundaries
  !----------------------------------------
  if (nxtnest(ngrid)>=1) call nstbdriv()
  !if(stepDebug) print *,mynum,nstbdriv'
  !  Rayleigh friction for theta
  !----------------------------------------
  call RAYFT()
  !if(stepDebug) print *,mynum,RAYFT'
  !  Get the overlap region between parallel nodes
  !---------------------------------------------------
  if (nmachs > 1) then
     call WaitRecvMsgs(OneGrid%SelectedGhostZoneSend, OneGrid%SelectedGhostZoneRecv)
  endif

  if (iexev == 2) &
     call exevolve(mzp,mxp,myp,ngrid,ia,iz,ja,jz,izu,jzv,jdim,mynum,dtlt,'THA')
  !if(stepDebug) print *,mynum,exevolve'
  !- task 2:  NO production by "eclair"
  if (ccatt == 1) &
     call chemistry_driver(mzp,mxp,myp,ia,iz,ja,jz,2,50)

  !- CATT & Chemistry == CCATT
  !----------------------------------------
  if (ccatt==1 .and. split_method== 'PARALLEL' .and. n_dyn_chem==1) then
     ! task 3 : production/loss by chemical processes and inclusion of the
     ! chemistry tendency at the total tendency
     CALL chemistry_driver(mzp,mxp,myp,ia,iz,ja,jz,3,50)
  endif
  if (ccatt==1 ) then
     ! task 4 : mass transfer between gas and liquid
     call chemistry_driver(mzp,mxp,myp,ia,iz,ja,jz,4,50)
  endif

  !---------------------------------------------------
  ! Shallow  cumulus parameterization by Souza
  if (NNSHCU(ngrid)==1) call SHCUPA()
  !---------------------------------------------------
  !if(stepDebug) print *,mynum,SHCUPA'
  
  if (TEB_SPM==1) then
     ! Update urban emissions
     !----------------------------------------
     if (isource==1) then
        call sources_teb(mzp, mxp, myp, ia, iz, ja, jz, ngrid, ngrids)
     endif
     !  Update chemistry
     !----------------------------------------
     if (ichemi==1) then
        call ozone(mzp, mxp, myp, ia, iz, ja, jz, ngrid, dtlt)
     endif
  endif

  if (iexev == 2) &
  call exevolve(mzp,mxp,myp,ngrid,ia,iz,ja,jz,izu,jzv,jdim,mynum,dtlt,'THS')


  !  Sub-grid diffusion terms
  !----------------------------------------
  if ((if_adap==0) .and. (ihorgrad==2)) then
     call diffuse_brams31() !call optimized subroutine
  else
     call diffuse()
  endif !;    call dumpVarAllLatLonk(tend%wt, 'WT'  ,412,0,0,1,mxp,1,myp,1,mzp,0.0,0.0)

  !- STILT-BRAMS coupling (ML)
  if (imassflx == 1) call prep_advflx_to_stilt(mzp,mxp,myp,ia,iz,ja,jz,ngrid)

  !- large and subgrid scale forcing for shallow and deep cumulus
  IF( NNQPARM(ngrid) >=2  ) CALL prepare_lsf(NNQPARM(ngrid), NNSHCU(ngrid),1)

  !- cumulus parameterizations options: G3d - GD-FIM and GF
  IF(NNQPARM(ngrid)>=3) CALL CUPARM_GRELL3_CATT(OneGrid,1,NNQPARM(ngrid),NNSHCU(ngrid))

  !------------------------------------------------------------------------------
  ! init preparations for Runge-Kutta  -loop
  !------------------------------------------------------------------------------

  if ( rk_order == 2 ) then
    ! Wicker, Skamarock (1998)-RK-scheme
    rk_beta(1) = 1.0 / 2.0    ! = beta(2,1) of Butcher tableau
    rk_beta(2) = 1.0          ! = beta(3,2) of Butcher tableau

    if ( MOD( nnacoust(ngrid), 2) /= 0 ) then
       call fatal_error("ERROR in timestep_rk: nnacoust(ngrid) must be an integer multiple of 2")
    end if

    rk_nmbr_small_timesteps(1) = nnacoust(ngrid) / 2
    rk_nmbr_small_timesteps(2) = nnacoust(ngrid)

  else if ( rk_order == 3 ) then
    ! Wicker, Skamarock (2002)-RK-scheme
    rk_beta(1) = 1.0 / 3.0    ! = beta(2,1) of Butcher tableau
    rk_beta(2) = 1.0 / 2.0    ! = beta(3,2) of Butcher tableau
    rk_beta(3) = 1.0          ! = beta(4,3) of Butcher tableau

    if ( MOD( nnacoust(ngrid), 6) /= 0 .and.  MOD( nnacoust(ngrid), 4) /= 0 ) then
       call fatal_error("ERROR in timestep_rk: nnacoust(ngrid) must be an integer multiple of 6 or 4")
    end if
    rk_nmbr_small_timesteps(1) = nnacoust(ngrid) / 3
    rk_nmbr_small_timesteps(2) = nnacoust(ngrid) / 2
    rk_nmbr_small_timesteps(3) = nnacoust(ngrid)
    !------
    !- optimization : rk_nmbr_small_timesteps(1)=1 e dts=2* dtlt / nnacoust(ngrid) - uncoment the line below
    rk_nmbr_small_timesteps(1) = 1.
    !------

  else
     call fatal_error("ERROR in timestep_rk: false value for rk_order")
  end if

  dts = dtlt / nnacoust(ngrid)

  if ( apply_div_damping ) then
     ! calculation of divergence damping coeff:
     !MB: ATTENTION: dts different for different nests?!!
     if ( ideltat == 0 ) then
        ! it is sufficient to calculate alpha_div only once:
        if ( istp == 1 ) then
           call init_div_damping_coeff( dts )
        end if
     else
        call init_div_damping_coeff( dts )
     end if
  end if


  ! begin of Runge-Kutta loop
  !---------------------------------------------------

  !MB>>
  !nmbr_gpts = mxp * myp * mzp    !MB: only for testing!!!
  !write(*,"(A,I2,4F15.8)") "uc  ", 0, minval( basic_g(ngrid)%uc  ), maxval( basic_g(ngrid)%uc  ), &
  !	  sum( basic_g(ngrid)%uc  )/nmbr_gpts, basic_g(ngrid)%uc(7,8,9)
  !write(*,"(A,I2,4F15.8)") "vc  ", 0, minval( basic_g(ngrid)%vc  ), maxval( basic_g(ngrid)%vc  ), &
  !	  sum( basic_g(ngrid)%vc  )/nmbr_gpts, basic_g(ngrid)%vc(7,8,9)
  !write(*,"(A,I2,4F15.8)") "wc  ", 0, minval( basic_g(ngrid)%wc  ), maxval( basic_g(ngrid)%wc  ), &
  !	  sum( basic_g(ngrid)%wc  )/nmbr_gpts, basic_g(ngrid)%wc(7,8,9)
  !write(*,"(A,I2,4F15.8)") "pc  ", 0, minval( basic_g(ngrid)%pc  ), maxval( basic_g(ngrid)%pc  ), &
  !	  sum( basic_g(ngrid)%pc  )/nmbr_gpts, basic_g(ngrid)%pc(7,8,9)
  !write(*,"(A,I2,4F15.8)") "thc ", 0, minval( basic_g(ngrid)%thc ), maxval( basic_g(ngrid)%thc ), &
  !	  sum( basic_g(ngrid)%thc )/nmbr_gpts, basic_g(ngrid)%thc(7,8,9)
  !
  !write(*,"(A,I2,3F15.8)") "ut  ", 0, minval( tend%ut  ), maxval( tend%ut  ), sum( tend%ut  )/nmbr_gpts
  !write(*,"(A,I2,3F15.8)") "vt  ", 0, minval( tend%vt  ), maxval( tend%vt  ), sum( tend%vt  )/nmbr_gpts
  !write(*,"(A,I2,3F15.8)") "wt  ", 0, minval( tend%wt  ), maxval( tend%wt  ), sum( tend%wt  )/nmbr_gpts
  !write(*,"(A,I2,3F15.8)") "pt  ", 0, minval( tend%pt  ), maxval( tend%pt  ), sum( tend%pt  )/nmbr_gpts
  !write(*,"(A,I2,3F15.8)") "tht ", 0, minval( tend%tht ), maxval( tend%tht ), sum( tend%tht )/nmbr_gpts
  !!MB<<

  !MB:
  !if ( flag_mb_adv_test) then
  !  !! ascii output for horizontal scalar advection test:
  !  !do i=1, mxp
  !  !  write(*,*) "thc  ", i, basic_g(ngrid)%thc(6,i,5)
  !  !end do
  !  !! ascii output for vertical scalar advection test:
  !  do k=1, mzp
  !    write(*,*) "thc  ", k, zt(k), basic_g(ngrid)%thc(k,10,5)
  !  end do
  !end if

  !  Lateral velocity boundaries - radiative
  !-------------------------------------------
  call LATBND()

  do l_rk = 1, rk_order
     !write(crk,fmt='(I2.2)') l_rk
     !print*,"rk:",rk_beta(l_rk),rk_nmbr_small_timesteps(l_rk),  rk_nmbr_small_timesteps(l_rk)*dts,dtlt

     !initialize the tendencies with the physics tendencies
     tend%ut_rk (:) =tend%ut (:)
     tend%vt_rk (:) =tend%vt (:)
     tend%wt_rk (:) =tend%wt (:)
     tend%pt_rk (:) =tend%pt (:) !; call dumpVarAllLatLonk(tend%tht,'THT',515,l_rk,0,1,mxp,1,myp,1,mzp,0.0,0.0)
     tend%tht_rk(:) =tend%tht(:) !; call dumpVarAllLatLonk(tend%tht_rk,'THT_RK',516,l_rk,0,1,mxp,1,myp,1,mzp,0.0,0.0)

     !.. call dumpVarAllLatLonk3P(tend%ut_rk , 'UT'  ,528,l_rk,0,1,mxp,1,myp,1,mzp,0.0,0.0,h)
     !.. call dumpVarAllLatLonk3P(tend%vt_rk , 'VT'  ,528,l_rk,0,1,mxp,1,myp,1,mzp,0.0,0.0,h)
     !.. call dumpVarAllLatLonk3P(tend%wt_rk , 'WT'  ,528,l_rk,0,1,mxp,1,myp,1,mzp,0.0,0.0,h)
     !.. call dumpVarAllLatLonk3P(tend%pt_rk , 'PT'  ,528,l_rk,0,1,mxp,1,myp,1,mzp,0.0,0.0,h)
     !.. call dumpVarAllLatLonk3P(tend%tht_rk,'THT'  ,528,l_rk,0,1,mxp,1,myp,1,mzp,0.0,0.0,h)

     ! advection should give back tendencies
     ! ut_rk, vt_rk, wt_rk, pt_rk, tht_rk = physics tend + advection tendency

     !  Velocity advection
     !----------------------------------------
     call advectc_rk('V',mzp,mxp,myp,ia,iz,ja,jz,izu,jzv,mynum,l_rk)!;  call dumpVarAllLatLonk(tend%wt_rk , 'WT'  ,530,l_rk,0,1,mxp,1,myp,1,mzp,0.0,0.0)

     !----------------------------------------
     !call dumpVarAllLatLonk(basic_g(ngrid)%thp,'aTHPe'//crk,1,mxp,1,myp,1,mzp,0.0,0.0) !ok na 1
     !call dumpVarAllLatLonk(basic_g(ngrid)%thc,'aTHCe'//crk,1,mxp,1,myp,1,mzp,0.0,0.0) !ok na 1

     !  advection of pi and theta_il
     call advectc_rk('THETAIL',mzp,mxp,myp,ia,iz,ja,jz,izu,jzv,mynum,l_rk)!;  call dumpVarAllLatLonk(tend%wt_rk , 'WT'  ,536,l_rk,0,1,mxp,1,myp,1,mzp,0.0,0.0)
     call advectc_rk('PI'     ,mzp,mxp,myp,ia,iz,ja,jz,izu,jzv,mynum,l_rk)!;  call dumpVarAllLatLonk(tend%wt_rk , 'WT'  ,537,l_rk,0,1,mxp,1,myp,1,mzp,0.0,0.0)

     if ( flag_Coriolis_in_every_RK_step ) then
       !if ( .not. flag_mb_adv_test)
          call CORLOS(mzp,mxp,myp,i0,j0,ia,iz,ja,jz,izu,jzv, tend%ut_rk, tend%vt_rk)
       !end if
     end if

     !  Buoyancy term for w equation
     !----------------------------------------
     !if ( .not. flag_mb_adv_test)
        call BUOYANCY( tend%wt_rk )
     !end if

     if ( l_rk > 1 ) then
       ! (not necessary in the first RK substep)
       basic_g(ngrid)%uc (:,:,:) = basic_g(ngrid)%up (:,:,:)
       basic_g(ngrid)%vc (:,:,:) = basic_g(ngrid)%vp (:,:,:)
       basic_g(ngrid)%wc (:,:,:) = basic_g(ngrid)%wp (:,:,:)
       basic_g(ngrid)%pc (:,:,:) = basic_g(ngrid)%pp (:,:,:)
       basic_g(ngrid)%thc(:,:,:) = basic_g(ngrid)%thp(:,:,:)
     end if

     !-  Acoustic small timesteps
     if ( .not. flag_mb_adv_test) then
       call acoustic_new(OneGrid, rk_nmbr_small_timesteps(l_rk),l_rk )
     else
      ! just for testing:
      ! without the acoustic subroutine, tendencies must be added "manually"
      ! update pp -> pc (similar for u, v, w, if needed.)
      call update_long_rk(int(mxp*myp*mzp,i8),dtlt,rk_beta(l_rk) &
                        ,basic_g(ngrid)%pc,basic_g(ngrid)%pp  &
                        ,tend%pt_rk)
     endif
    
     !- update thp -> thc (theta_il is not contained in acoustic_new)
     if (.not. singleProcRun) then
        call commHaloAcou(tend%tht_rk,mxp,myp,mzp,myNum,'tht_rk')
     endif
     call update_long_rk(int(mxp*myp*mzp,i8),dtlt,rk_beta(l_rk) &
                        ,basic_g(ngrid)%thc,basic_g(ngrid)%thp  &
                        ,tend%tht_rk)

     !- determine theta (dry potential temp.) for the buoyancy term:
     call theta_thp_rk(mzp,mxp,myp,ia,iz,ja,jz,"get_theta")

     !-damping on vertical velocity to keep stability
     !MB: does this act on wc???
     if(vveldamp == 1) call w_damping(mzp,mxp,myp,ia,iz,ja,jz,mynum)


  end do
  ! end of Runge-Kutta loop
  ! i.e. the fields
  !    basic_g%uc, ..%vc, ..%wc, ..%pc, ..%thc
  ! contain the fields at the new time level n+1
  !---------------------------------------------------
  !
  if ( apply_div_damping ) then
     if ( ideltat /= 0 ) then
        call deallocate_alpha_div
     end if
  end if

  !
  !
  !  water species, tke and tracers advection
  !----------------------------------------
  IF(advmnt == 1) THEN
      !- monotonic advection scheme
      CALL advmnt_driver('SCALAR',mzp,mxp,myp,ia,iz,ja,jz,izu,jzv,mynum)
  ELSEIF(advmnt == 0) THEN
      !- using the 2nd order forward upstream
      CALL ADVECTc   ('SCALAR',mzp,mxp,myp,ia,iz,ja,jz,izu,jzv,mynum)
  ELSEIF(advmnt == 3) THEN
      !- using the WS advection
      CALL advectc_rk('SCALAR',mzp,mxp,myp,ia,iz,ja,jz,izu,jzv,mynum,l_rk)

  ENDIF

  !  Update scalars (water species, tke and tracers)
  !----------------------------------------
   call PREDTR()

  !-  copy current time into past time (u,v,w,pi,thetail)
  !---> thp   must be changed to THC for microphysics/bc/theta update
  !---> pp    must be changed to PC  for microphysics
  !---> wp    must be changed to WC  for microphysics
  !---> up,vp must be changed to UC,VC for output
   basic_g(ngrid)%up (:,:,:) = basic_g(ngrid)%uc (:,:,:)
   basic_g(ngrid)%vp (:,:,:) = basic_g(ngrid)%vc (:,:,:)
   basic_g(ngrid)%wp (:,:,:) = basic_g(ngrid)%wc (:,:,:)
   basic_g(ngrid)%pp (:,:,:) = basic_g(ngrid)%pc (:,:,:)
   basic_g(ngrid)%thp(:,:,:) = basic_g(ngrid)%thc(:,:,:)
  !---->
  !---->
  !
  !  Moisture variables positive definite
  !----------------------------------------
  if     (mcphys_type == 0) then
           call negadj1(mzp,mxp,myp)

   elseif(mcphys_type == 1) then
           call negadj1_2M_rams60(mzp,mxp,myp)
  endif

  !  Microphysics (applied on THP, just updated)
  !----------------------------------------
  if (mcphys_type == 0 .and. level==3) then
     if (machine==1 .and. TEB_SPM==0) then
        !- optimized version only for SX-6
        call micro_opt()
     else
        !- original Version used in a Generic IA32 machine
        call micro()
     endif
  endif
  
  if (mcphys_type == 1 .and. level==3) then
        !- 2M rams microphysics
        call micro_2M_rams60()

  elseif (mcphys_type == 2 .or. mcphys_type == 3 ) then
        !- G. Thompson microphysics
        call micro_thompson()

  elseif(mcphys_type == 4 ) then
        call micro_gfdl()
  
  elseif(mcphys_type >= 5 ) then
        call micro_wsm()
  endif
  !=====================================================

  if (sfire > 0) call sfclyr_sfire(mzp,mxp,myp,ia,iz,ja,jz)

  !=====================================================

  !- Thermodynamic diagnosis
  if (mcphys_type <= 1 .and. level==3)  then
     call THERMO(mzp,mxp,myp,1,mxp,1,myp,'MICRO')
  endif

  !  Apply scalar b.c.'s (THP is changed here)
  !----------------------------------------
  call TRSETS()

  !---> THC must be changed to THP to include microphysics/trsets changes
  !---> for the next timestep
        basic_g(ngrid)%thc(:,:,:) = basic_g(ngrid)%thp(:,:,:)
  !--->

  !  Lateral velocity boundaries - radiative
  !-------------------------------------------
  !srf  call LATBND()

  !  Velocity/pressure boundary conditions
  !----------------------------------------
  call VPSETS()

  !- call THERMO on the boundaries
  call thermo_boundary_driver((time+dtlongn(ngrid)), dtlong, &
       f_thermo_e, f_thermo_w, f_thermo_s, f_thermo_n, &
       nzp, mxp, myp, jdim)

  if (iexev == 2) call get_true_air_density(mzp,mxp,myp,ia,iz,ja,jz)

  !----------------------------------------
  !- chemistry - microphysics tranfers - sedimentation and tranfer from clouds to rain
  if (ccatt==1) then
     ! task 5 : sedimentation and mass transfer between clouds and rain
     call chemistry_driver(mzp,mxp,myp,ia,iz,ja,jz,5,50)
  endif
  !
  !----------------------------------------
  !- chemistry/aerosol solvers
  if (ccatt==1) THEN
    if ( (split_method== 'PARALLEL' .AND. N_DYN_CHEM > 1) .OR.  &
         (split_method== 'SEQUENTIAL'                   ) .OR.  &
 	 (split_method== 'SYMMETRIC'                    )       )THEN

       ! task 3 : production/loss by chemical processes and final updated
       !  	of each specie
       call chemistry_driver(mzp,mxp,myp,ia,iz,ja,jz,3,50)
    endif

    !- call Matrix Aerosol Model
    !- using symmetric/sequential spliting operator
    if(AEROSOL==2) then

       CALL MatrixDriver(ia,iz,ja,jz,mzp,mxp,myp)
    endif
  endif
  if (ccatt==1 .and. aerosol == 1) THEN
    call aer_background(ngrid,mzp,mxp,myp,ia,iz,ja,jz)
  endif
  !----------------------------------------

  if (TEB_SPM==1) then
     !EDF  emission module
     if (isource==1) then
        ! Apply only for last finner grid
        if (ngrid==ngrids) then
           CALL le_fontes(ngrid, mzp, mxp, myp, &
                npatch, ia, iz, ja, jz, (time+dtlongn(1)))
        endif
     endif
     !EDF
  endif

 !- windfarm
  call wind_farm_driver(ngrid,mzp,mxp,myp,ia,iz,ja,jz)

 !- apply digital filter
 if(applyDF) call applyDigitalFilter(fileNameDF, dfVars)


 ! Implements the Incremental Analysis Update procedure -
 !  phase 1: create and store the IAU tendencies
 !----------------------------------------
 !
 if(applyIAU == 1 .and. abs ( time + dtlt - timeWindowIAU*0.5) < 0.01) then
    
    if(mynum==1) print*,"timeIAU=",time,timeWindowIAU*0.5,&
           abs ( time - dtlt - timeWindowIAU*0.5),applyIAU
    if(mynum==1)call flush(6)
    
    call CreateIauTendency(ngrid, mzp*mxp*myp, mzp, mxp, myp,ia,iz,ja,jz&
          ,varinit_g(ngrid)%varup(:,:,:),varinit_g(ngrid)%varvp(:,:,:)  &
          ,varinit_g(ngrid)%varpp(:,:,:),varinit_g(ngrid)%vartp(:,:,:)  &
          ,varinit_g(ngrid)%varrp(:,:,:)                                &
	  
          ,varinit_g(ngrid)%varuf(:,:,:),varinit_g(ngrid)%varvf(:,:,:)  &
          ,varinit_g(ngrid)%varpf(:,:,:),varinit_g(ngrid)%vartf(:,:,:)  &
          ,varinit_g(ngrid)%varrf(:,:,:)                                &

          ,basic_g(ngrid)%up     (:,:,:)   ,basic_g(ngrid)%vp  (:,:,:)  &
          ,basic_g(ngrid)%theta  (:,:,:)   ,basic_g(ngrid)%rtp (:,:,:)  &
          ,basic_g(ngrid)%pp     (:,:,:)                                )
  endif

!!$  call TimeStamp(TS_PHYSICS)
!if(stepDebug) print *,mynum,end of timestep'
end subroutine timestep_rk


end module ModTimestep_RK

!*************************************************************************
!-srf-  temp routine only for testing - not being used
!
subroutine adv_p_driver(m1,m2,m3,ifm,ia,iz,ja,jz,izu,jzv,jdim,mynum,edt,key)


   use mem_basic,   only: basic_g
   use mem_grid,    only: grid_g, itopo,dyncore_flag
   use mem_stilt,    only: stilt_g
   use mem_tend,    only: tend
   use mem_scratch, only: scratch
   use micphys,     only: level           !if(vapour_on)  use therm_lib,   only: vapour_on

   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   character(len=*) , intent(in) :: key
   integer          , intent(in) :: m1,m2,m3,ifm,ia,iz,ja,jz,izu,jzv,jdim,mynum
   real             , intent(in) :: edt
   !----- Local variables -----------------------------------------------------------------!
   integer :: i,j,k
   !---------------------------------------------------------------------------------------!
       call adv_p(m1,m2,m3,ia,iz,ja,jz,izu,jzv,jdim,itopo                                 &
                  ,grid_g(ifm)%rtgu                       ,grid_g(ifm)%fmapui              &
                  ,grid_g(ifm)%rtgv                       ,grid_g(ifm)%fmapvi              &
                  ,grid_g(ifm)%f13t                       ,grid_g(ifm)%f23t                &
                  ,grid_g(ifm)%rtgt                       ,grid_g(ifm)%fmapt               &
                  ,grid_g(ifm)%dxt                        ,grid_g(ifm)%dyt                 &
                  ,basic_g(ifm)%uc                        ,basic_g(ifm)%dn0u               &
                  ,basic_g(ifm)%vc                        ,basic_g(ifm)%dn0v               &
                  ,basic_g(ifm)%dn0                       ,basic_g(ifm)%wc                 &
                  ,basic_g(ifm)%pc                        ,tend%pt                         )

end subroutine adv_p_driver

subroutine adv_p(m1,m2,m3,ia,iz,ja,jz,izu,jzv,jdim,itopo,rtgu,fmapui,rtgv,fmapvi,f13t    &
                  ,f23t,rtgt,fmapt,dxt,dyt,uc,dn0u,vc,dn0v,dn0,wc,pc,pt)

   use mem_grid , only : hw4 & ! intent(in)
                       , dzt ! ! intent(in)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer , intent(in)                         :: m1,m2,m3,ia,iz,ja,jz,izu,jzv,itopo,jdim
   real    , intent(in)   , dimension(m2,m3)    :: rtgu,fmapui,rtgv,fmapvi,f13t,f23t
   real    , intent(in)   , dimension(m2,m3)    :: rtgt,fmapt,dxt,dyt
   real    , intent(in)   , dimension(m1,m2,m3) :: uc,dn0u,vc,dn0v,dn0,wc,pc
   real    , intent(inout), dimension(m1,m2,m3) :: pt
   !----- Local variables -----------------------------------------------------------------!
   integer                                      :: i,j,k,im,jm
   real                                         :: c1z,c1x,c1y
   real    ,                dimension(m1,m2,m3) :: flxu,flxv,flxw
   !---------------------------------------------------------------------------------------!



   !----- Compute momentum fluxes flxu, flxv, flxw ----------------------------------------!
   do j = 1,m3
      do i = 1,m2
         do k = 1,m1
            flxu(k,i,j) = uc(k,i,j) * dn0u(k,i,j) * rtgu(i,j) * fmapui(i,j)
            flxv(k,i,j) = vc(k,i,j) * dn0v(k,i,j) * rtgv(i,j) * fmapvi(i,j)
         enddo
      enddo
   enddo

   if(itopo == 0) then
      do j = 1,m3
         do i = 1,m2
            do k = 1,m1-1
               flxw(k,i,j) = wc(k,i,j) * .5 * (dn0(k,i,j) + dn0(k+1,i,j))
            end do
         end do
      end do
   else
      do j = 1,m3
         jm = max(j-1,1)
         do i = 1,m2
            im = max(i-1,1)
            do k = 1,m1-1
               flxw(k,i,j) = wc(k,i,j) * .5 * (dn0(k,i,j) + dn0(k+1,i,j))                  &
                           + hw4(k) * ( ( flxu(k,i,j) + flxu(k+1,i,j)                      &
                                        + flxu(k,im,j) + flxu(k+1,im,j) ) * f13t(i,j)      &
                                      + ( flxv(k,i,j) + flxv(k+1,i,j)                      &
                                        + flxv(k,i,jm) + flxv(k+1,i,jm) ) * f23t(i,j) )
            end do
         end do
      end do
   end if

   !---------------------------------------------------------------------------------------!
   !  Compute advection contribution of zonal gradient to Exner function tendency.         !
   !---------------------------------------------------------------------------------------!
   do j = ja,jz
      do i = ia,izu
         c1x = 0.5 / rtgt(i,j) * fmapt(i,j) * dxt(i,j)
         do k = 2,m1-1
            pt(k,i,j) = pt(k,i,j)                                                          &
                      - c1x / dn0(k,i,j)                                                   &
                      * ( flxu(k,i,j)   * (pc(k,i,j) + pc(k,i+1,j))                        &
                        - flxu(k,i-1,j) * (pc(k,i,j) + pc(k,i-1,j))                        &
                        - (flxu(k,i,j) - flxu(k,i-1,j)) * 2.* pc(k,i,j) )
          end do
       end do
    end do

   !---------------------------------------------------------------------------------------!
   !  Compute advection contribution of meridional gradient to Exner function tendency.    !
   !---------------------------------------------------------------------------------------!
   do j=ja,jzv
      do i=ia,iz
         c1y = 0.5 / rtgt(i,j) * fmapt(i,j) * dyt(i,j)
         do k=2,m1-1
            pt(k,i,j) = pt(k,i,j)                                                          &
                      - c1y /dn0(k,i,j)                                                    &
                      * ( flxv(k,i,j)      * (pc(k,i,j)+pc(k,i,j+jdim))                    &
                        - flxv(k,i,j-jdim) * (pc(k,i,j)+pc(k,i,j-jdim))                    &
                        - (flxv(k,i,j)-flxv(k,i,j-jdim)) * 2.* pc(k,i,j) )
         end do
      end do
   end do

   !---------------------------------------------------------------------------------------!
   !  Compute advection contribution of vertical gradient to Exner function tendency.      !
   !---------------------------------------------------------------------------------------!
   do j=ja,jz
      do i=ia,iz
         c1z = 0.5 / rtgt(i,j)
         do k=2,m1-1
            pt(k,i,j) = pt(k,i,j)                                                          &
                      - c1z * dzt(k) /dn0(k,i,j)                                           &
                      * ( flxw(k,i,j)   * (pc(k,i,j)+pc(k+1,i,j))                          &
                        - flxw(k-1,i,j) * (pc(k,i,j)+pc(k-1,i,j))                          &
                        -  (flxw(k,i,j)-flxw(k-1,i,j)) * 2. * pc(k,i,j) )
         end do
      end do
   end do
   return
end subroutine adv_p
!==========================================================================================!
