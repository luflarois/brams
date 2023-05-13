!###########################################################################
! CCATT- B - Regional Atmospheric Modeling System - RAMS
!###########################################################################
MODULE chem_sources

     
  USE mem_chem1, ONLY:       &
       maxsrcfiles,          & ! Parameter
       chem1_vars,           & ! Type
       chem1_src_vars          ! Type

  USE mem_aer1, ONLY:        &
       aer1_vars
       
  USE mem_plume_chem1, ONLY: &
       plume_vars,           & ! Type
       plume_mean_vars,      & ! Type
       plume_fre_vars,       &
       iflam_frac ,          &
       imean_frp  ,	     &
       istd_frp   ,	     &
       imean_size ,	     &
       istd_size 	 

  USE mem_volc_chem1, ONLY:  &
       volc_mean_vars          ! Type

  USE ModNamelistFile, ONLY: &
       namelistFile

  USE ModDateUtils, ONLY: &
       date_abs_secs2,    &    ! Subroutine
       date_add_to,       &    ! Subroutine
       date_make_big,     &    ! Subroutine
       julday                  ! Function

  USE ReadBcst, ONLY: &
       ReadStoreFullFieldAndOwnChunk, & ! Subroutine
       Broadcast                        ! Subroutine

  use ParLib, only: parf_bcast ! Subroutine

  USE io_params, ONLY:    &
       srctime1,          &
       srctime2,          &
       fnames_src,        &
       itotdate_src,      &
       src_times,         &
       actual_time_index, &
       nsrcfiles,         &
       next_srcfile


  IMPLICIT NONE

  INCLUDE "ranks.h"


  CHARACTER(len=256)  :: srcmapfn

  CHARACTER(len=32)   :: def_proc_src ! 'last_sources' 
                                      !    if does not exist => keep the current sources
                                      ! 'stop' 
                                      !    if does not exist => stop the execution

  LOGICAL             :: got_srcfiles_inv = .FALSE.

  INTEGER             :: src_swap=0

  !- arrays for the diurnal cycle of emission
  TYPE, PUBLIC ::  cycle_emission
     DOUBLE PRECISION, POINTER, DIMENSION(:,:) :: dcnorma_inv,   &
                                                  emission_rate, &
                                                  emission_rate_NOX
  END TYPE cycle_emission
  REAL, PARAMETER                                   :: emiss_cycle_time = 86400.
  LOGICAL                                           :: emiss_cycle_alloc = .FALSE.
  TYPE(cycle_emission), ALLOCATABLE, DIMENSION(:,:) :: emiss_cycle

  PRIVATE

  TYPE GlobalEmissData
     REAL, POINTER :: sc_src(:,:)        ! chem1
     REAL, POINTER :: src_dummy_2d(:,:)  ! chem1
     REAL, POINTER :: flam_frac(:,:)     ! plume_chem1
     REAL, POINTER :: fire_size(:,:)     ! plume_chem1
     REAL, POINTER :: plum_heigth(:,:)   ! volc_chem1
     REAL, POINTER :: vent_elev(:,:)     ! volc_chem1
     REAL, POINTER :: duration (:,:)     ! volc_chem1
     REAL, POINTER :: mean_frp (:,:)     ! plume_chem1
     REAL, POINTER :: std_frp  (:,:)     ! plume_chem1
     REAL, POINTER :: mean_size(:,:)     ! plume_chem1
     REAL, POINTER :: std_size (:,:)     ! plume_chem1
  END TYPE GlobalEmissData

  TYPE(GlobalEmissData), ALLOCATABLE, TARGET :: oneGlobalEmissData(:)

!!$  INTEGER, PARAMETER                              :: maxsrcfiles   = 1500 ! Movido para mem_chem1
!!$  REAL                                            :: srctime1=0.          ! Movido para io_params.f90
!!$  REAL                                            :: srctime2=0.          ! Movido para io_params.f90
!!$  CHARACTER(len=256), ALLOCATABLE, DIMENSION(:,:) :: fnames_src           ! Movido para io_params.f90
!!$  CHARACTER(len=14),  ALLOCATABLE, DIMENSION(:,:) :: itotdate_src         ! Movido para io_params.f90
!!$  REAL,               ALLOCATABLE, DIMENSION(:,:) :: src_times            ! Movido para io_params.f90
!!$  INTEGER,            ALLOCATABLE, DIMENSION(:,:) :: actual_time_index    ! Movido para io_params.f90
!!$  INTEGER,            ALLOCATABLE, DIMENSION(:)   :: nsrcfiles            ! Movido para io_params.f90
!!$  INTEGER,            ALLOCATABLE, DIMENSION(:)   :: next_srcfile         ! Movido para io_params.f90

  INTEGER,PARAMETER :: nucle = 1 ! nucleation mode
  INTEGER,PARAMETER :: accum = 2 ! accumulation mode
  INTEGER,PARAMETER :: coarse = 3 ! coarse mode
  INTEGER,PARAMETER :: aer_bburn=002
  INTEGER,PARAMETER :: urban=003
  INTEGER,PARAMETER :: bioge=004
  INTEGER,PARAMETER :: marin=005
  INTEGER,PARAMETER :: v_ash=006

  PUBLIC :: read_sourcemaps,                & ! Subroutine
            sources,                        & ! Subroutine
            init_actual_time_index,         & ! Subroutine
            alloc_emiss_cycle,              & ! Subroutine
            get_diurnal_cycle_normalized,   & ! Subroutine
            oneGlobalEmissData,             & ! Type(GlobalEmissData)
            StoreNamelistFileAtChemSources, & ! Subroutine
            alloc_GlobalEmissData,          & ! Subroutine
            nullify_GlobalEmissData,        & ! Subroutine
            dealloc_GlobalEmissData,        & ! Subroutine
            def_proc_src,                   &
            srcmapfn,                       &
            emiss_cycle_alloc,              &
            emiss_cycle_time,               &
            emiss_cycle

CONTAINS



  !---------------------------------------------------------------------------- 
  SUBROUTINE alloc_GlobalEmissData(global, n2, n3)
    IMPLICIT NONE
    TYPE (GlobalEmissData), INTENT (inout) :: global
    INTEGER,                INTENT(in)     :: n2
    INTEGER,                INTENT(in)     :: n3

    ALLOCATE (global%sc_src(n2,n3))
    ALLOCATE (global%src_dummy_2d(n2,n3))
    ALLOCATE (global%flam_frac(n2,n3))
    ALLOCATE (global%fire_size(n2,n3))
    ALLOCATE (global%plum_heigth(n2,n3))
    ALLOCATE (global%vent_elev(n2,n3))
    ALLOCATE (global%duration(n2,n3))
    ALLOCATE (global%mean_frp(n2,n3))
    ALLOCATE (global%std_frp(n2,n3))
    ALLOCATE (global%mean_size(n2,n3))
    ALLOCATE (global%std_size(n2,n3))

  END SUBROUTINE alloc_GlobalEmissData



  !---------------------------------------------------------------------------- 
  SUBROUTINE nullify_GlobalEmissData(global)
    IMPLICIT NONE
    TYPE (GlobalEmissData), INTENT (inout) :: global

    NULLIFY (global%sc_src)
    NULLIFY (global%src_dummy_2d)
    NULLIFY (global%flam_frac)
    NULLIFY (global%fire_size)
    NULLIFY (global%plum_heigth)
    NULLIFY (global%vent_elev)
    NULLIFY (global%duration)
    NULLIFY (global%mean_frp )
    NULLIFY (global%std_frp  )
    NULLIFY (global%mean_size)
    NULLIFY (global%std_size )

  END SUBROUTINE nullify_GlobalEmissData



  !---------------------------------------------------------------------------- 
  SUBROUTINE dealloc_GlobalEmissData(global)
    IMPLICIT NONE
    TYPE (GlobalEmissData), INTENT (inout) :: global

    IF (ASSOCIATED(global%sc_src))      DEALLOCATE (global%sc_src)
    IF (ASSOCIATED(global%src_dummy_2d)) DEALLOCATE (global%src_dummy_2d)
    IF (ASSOCIATED(global%flam_frac))   DEALLOCATE (global%flam_frac)
    IF (ASSOCIATED(global%fire_size))   DEALLOCATE (global%fire_size)
    IF (ASSOCIATED(global%plum_heigth)) DEALLOCATE (global%plum_heigth)
    IF (ASSOCIATED(global%vent_elev))   DEALLOCATE (global%vent_elev)
    IF (ASSOCIATED(global%duration))   DEALLOCATE (global%duration)
    IF (ASSOCIATED(global%mean_frp))   DEALLOCATE (global%mean_frp)
    IF (ASSOCIATED(global%std_frp))   DEALLOCATE (global%std_frp)
    IF (ASSOCIATED(global%mean_size))   DEALLOCATE (global%mean_size)
    IF (ASSOCIATED(global%std_size))   DEALLOCATE (global%std_size)

  END SUBROUTINE dealloc_GlobalEmissData



!  !---------------------------------------------------------------------------- 
!  SUBROUTINE alloc_emiss_cycle(mmxp,mmyp,ngrids,nsrc)
!
!    ! node_mod
!    INTEGER , INTENT(IN) :: mmxp(ngrids)
!    INTEGER , INTENT(IN) :: mmyp(ngrids)
!    
!    ! mem_grid
!    INTEGER , INTENT(IN) :: ngrids
!    
!    ! mem_chem1
!    INTEGER , INTENT(IN) :: nsrc
!    
!    INTEGER ng,isrc
!
!    ALLOCATE (emiss_cycle(nsrc,ngrids))
!    DO ng=1,ngrids
!       DO isrc=1,nsrc      
!          ALLOCATE(emiss_cycle(isrc,ng)%dcnorma_inv( mmxp(ng),mmyp(ng) ) )    
!          emiss_cycle(isrc,ng)%dcnorma_inv = 0.
!          ALLOCATE(emiss_cycle(isrc,ng)%emission_rate( mmxp(ng),mmyp(ng) ) )    
!          emiss_cycle(isrc,ng)%emission_rate = 0.
!          ALLOCATE(emiss_cycle(isrc,ng)%emission_rate_NOX( mmxp(ng),mmyp(ng) ) )    
!          emiss_cycle(isrc,ng)%emission_rate_NOX = 0.
!       ENDDO
!    ENDDO
!
!    emiss_cycle_alloc=.TRUE.
!
!  END SUBROUTINE alloc_emiss_cycle


  !---------------------------------------------------------------------------- 
  !--(DMK-BRAMS 5.0)-----------------------------------------------------------
  ! Somente para BRAMS 5.0, 1 grade
  SUBROUTINE alloc_emiss_cycle(mxp,myp,ngrids,nsrc)

    ! node_mod
    INTEGER , INTENT(IN) :: mxp
    INTEGER , INTENT(IN) :: myp
    
    ! mem_grid
    INTEGER , INTENT(IN) :: ngrids
    
    ! mem_chem1
    INTEGER , INTENT(IN) :: nsrc
    
    INTEGER ng,isrc

    ALLOCATE (emiss_cycle(nsrc,ngrids))
    DO ng=1,ngrids
       DO isrc=1,nsrc      
          ALLOCATE(emiss_cycle(isrc,ng)%dcnorma_inv(mxp,myp) )    
          emiss_cycle(isrc,ng)%dcnorma_inv = 0.
          ALLOCATE(emiss_cycle(isrc,ng)%emission_rate(mxp,myp) )    
          emiss_cycle(isrc,ng)%emission_rate = 0.
          ALLOCATE(emiss_cycle(isrc,ng)%emission_rate_NOX(mxp,myp) )    
          emiss_cycle(isrc,ng)%emission_rate_NOX = 0.
       ENDDO
    ENDDO

    emiss_cycle_alloc=.TRUE.

  END SUBROUTINE alloc_emiss_cycle


  !---------------------------------------------------------------------------- 
  SUBROUTINE read_sourcemaps(ng,m1,m2,m3,ia,iz,ja,jz, &
                             time,iyear1,imonth1,idate1,itime1,&
                             ngrids,timmax,chem_nspecies,spc_chem_alloc,   &
                             src,off,nsrc,nvert_src,chem1_src_g,bburn,     & 
                             antro, bioge,geoge,                           &
                             spc_chem_name,on,chemical_mechanism,          &
                             emiss_ajust,co,aer_nspecies,spc_aer_alloc,    &
                             spc_aer_name,src_name,      &
                             chemistry,ntimes_src,aer1_g,nmodes,aerosol,   &
                             plumerise,nveg_agreg,plume_mean_g,nzpmax,dzt, &
                             rtgt,topt,transport,plume_g,        &
                             tropical_forest,boreal_forest,savannah,       &
                             grassland,diur_cycle,volcanoes,volc_mean_g,   &
                             dn0,zt,zm,mchnum,master_num,mass_bin_dist,CO2,ISFCL,&
			     aerosol_mechanism,plume_fre_g,emiss_ajust_aer)

    ! original
    INTEGER , INTENT(IN) :: ng
    INTEGER , INTENT(IN) :: m1
    INTEGER , INTENT(IN) :: m2
    INTEGER , INTENT(IN) :: m3
    INTEGER , INTENT(IN) :: ia
    INTEGER , INTENT(IN) :: iz
    INTEGER , INTENT(IN) :: ja
    INTEGER , INTENT(IN) :: jz
    INTEGER , INTENT(IN) :: isfcl,CO2  !DSM
    
    ! grid_dims
    INTEGER, INTENT(IN) :: nzpmax

    ! mem_grid
    REAL    , INTENT(IN) :: time
    INTEGER , INTENT(IN) :: iyear1
    INTEGER , INTENT(IN) :: imonth1
    INTEGER , INTENT(IN) :: idate1
    INTEGER , INTENT(IN) :: itime1
    INTEGER , INTENT(IN) :: ngrids
    REAL    , INTENT(IN) :: timmax
    REAL    , INTENT(IN) :: dzt(nzpmax)
    REAL    , INTENT(IN) :: zt(nzpmax)
    REAL    , INTENT(IN) :: zm(nzpmax)
    REAL    , INTENT(IN) :: rtgt(m2,m3)
    REAL    , INTENT(IN) :: topt(m2,m3)
    REAL    , INTENT(IN) :: dn0(m1,m2,m3)

    ! chem1_list
    INTEGER          , INTENT(IN) :: chem_nspecies
    INTEGER          , INTENT(IN) :: spc_chem_alloc(6,chem_nspecies)
    INTEGER          , INTENT(IN) :: src
    CHARACTER(LEN=8) , INTENT(IN) :: spc_chem_name(chem_nspecies)
    INTEGER          , INTENT(IN) :: on
    INTEGER          , INTENT(IN) :: off
    CHARACTER(LEN=24), INTENT(IN) :: chemical_mechanism
    REAL             , INTENT(IN) :: emiss_ajust(chem_nspecies)
    INTEGER          , INTENT(IN) :: CO
    INTEGER          , INTENT(IN) :: transport

    ! mem_chem1
    INTEGER              , INTENT(IN)    :: nsrc
    INTEGER              , INTENT(IN)    :: nvert_src(nsrc,ngrids)
    TYPE(chem1_src_vars) , INTENT(INOUT) :: chem1_src_g(2,nsrc,chem_nspecies,ngrids)
    INTEGER              , INTENT(IN)    :: bburn,antro, bioge,  geoge
    CHARACTER(LEN=20)    , INTENT(IN)    :: src_name(nsrc)
    INTEGER              , INTENT(IN)    :: chemistry
    INTEGER              , INTENT(IN)    :: ntimes_src(nsrc)
    INTEGER              , INTENT(IN)    :: diur_cycle(nsrc)

    ! aer1_list
    INTEGER          , INTENT(IN) :: aer_nspecies
    INTEGER          , INTENT(IN) :: spc_aer_alloc(6,nmodes,aer_nspecies)
    CHARACTER(LEN=8) , INTENT(IN) :: spc_aer_name(nmodes,aer_nspecies)
    INTEGER          , INTENT(IN) :: nmodes
    REAL             , INTENT(INOUT) :: mass_bin_dist(nmodes)
    CHARACTER(LEN=24), INTENT(IN) :: aerosol_mechanism
    REAL             , INTENT(IN) :: emiss_ajust_aer(nmodes,aer_nspecies) 
    
    ! mem_aer1
    TYPE (aer1_vars) , INTENT(INOUT) :: aer1_g(nmodes,aer_nspecies,ngrids)
    INTEGER          , INTENT(IN)    :: aerosol

    ! mem_plume_chem1
    INTEGER               , INTENT(IN)    :: nveg_agreg
    INTEGER               , INTENT(IN)    :: plumerise
    TYPE (plume_mean_vars), INTENT(INOUT) :: plume_mean_g(nveg_agreg,ngrids)
    TYPE (plume_fre_vars) , INTENT(INOUT) :: plume_fre_g(5,ngrids)
    TYPE (plume_vars)     , INTENT(INOUT) :: plume_g(nveg_agreg,chem_nspecies,ngrids)
    INTEGER               , INTENT(IN)    :: tropical_forest
    INTEGER               , INTENT(IN)    :: boreal_forest
    INTEGER               , INTENT(IN)    :: savannah
    INTEGER               , INTENT(IN)    :: grassland

    ! mem_volc_chem1
    TYPE (volc_mean_vars), INTENT(INOUT) :: volc_mean_g(ngrids)
    INTEGER              , INTENT(IN)    :: volcanoes    

    ! node_mod
    INTEGER, INTENT(IN)   :: mchnum
    INTEGER, INTENT(IN)   :: master_num

    !- local var
    INTEGER :: iunit,isrc,k1,k2,isrctime,iwtdo
    CHARACTER(len=256)  :: fname
    CHARACTER(len=2)   :: cgrid


    IF(srcmapfn(1:LEN_TRIM(srcmapfn)) == 'NONE' .OR. srcmapfn(1:LEN_TRIM(srcmapfn)) == 'none') RETURN

    IF( .NOT. got_srcfiles_inv ) THEN 

       ! Out: nsrcfiles(maxgrds), 
       !      fnames_src(maxsrcfiles,maxgrds), 
       !      itotdate_src(maxsrcfiles,maxgrds),
       !      src_times(maxsrcfiles,maxgrds), 
       !      next_srcfile(maxgrds), 
       !      srctime1, 
       !      srctime2, 
       !      got_srcfiles_inv

       CALL src_file_inv(srcmapfn,iyear1,imonth1,idate1,itime1,ngrids,time,timmax,nsrc,diur_cycle, &
                         mchnum,master_num)

       ! Out: actual_time_index(max_ntimes_src,nsrc)

       CALL init_actual_time_index(nsrc,ntimes_src)

    ENDIF

    ! In:  next_srcfile(maxgrds), 
    !      nsrcfiles(maxgrds)
    !
    ! Out: iwtdo, 
    !      srctime1, 
    !      srctime2
       
    !- check if files required does exist and decide what to do in case not
    IF(mchnum==master_num) CALL check_src_files(next_srcfile(ng),nsrcfiles(ng),iwtdo)
    
    !PRINT *, 'Aqui 1: LFR',mchnum,iwtdo; call flush(6)
    CALL Broadcast(iwtdo, master_num, "iwtdo")
    !PRINT *, 'Aqui 2: LFR',mchnum,iwtdo;call flush(6)
    

    IF(iwtdo == 0) RETURN 

    CALL Broadcast(nsrcfiles, master_num, "nsrcfiles")

    !- loop at source files (isrctime will be > 1, only if will use linterp for any source)
    DO isrctime=1,MAXVAL(ntimes_src(1:nsrc))


       !- swap sources:  copy time level 2->1 and read the next data for the time level 2
       IF(src_swap == 1 .AND. MAXVAL(ntimes_src(1:nsrc)) > 1 .AND. isrctime==1) THEN 

          ! Out: chem1_src_g
          CALL swap_sources(m1,m2,m3,time,chem_nspecies,spc_chem_alloc, &
                            src,off,nsrc,nvert_src(:,ng),chem1_src_g(:,:,:,ng),bburn,geoge)

          next_srcfile(ng) = next_srcfile(ng) + 1 
          CYCLE
       ENDIF

     
       !- next file to open
       fname=fnames_src(next_srcfile(ng),ng)

       ! InOut: chem1_src_g, 
       !        aer1_g,
       !        plume_mean_g
       !        volc_mean_g

       !- read emission dataset using V-format
       CALL read_sources_vfm(ng,m1,m2,m3,iyear1,imonth1,idate1,isrctime,    &
                             fname(1:LEN_TRIM(fname)),chem_nspecies,        &
                             spc_chem_alloc,spc_chem_name,src,on,off,       &
                             chemical_mechanism,emiss_ajust,co,aer_nspecies,&
                             spc_aer_alloc,spc_aer_name,urban,nucle,accum,  &
                             nsrc,chem1_src_g(:,:,:,ng),src_name,chemistry,aer1_g(:,:,ng), &
                             nmodes,aerosol,plumerise,nveg_agreg,plume_mean_g(:,ng),   &
			     volc_mean_g(ng),volcanoes,mchnum,master_num,mass_bin_dist,v_ash, &
			     aerosol_mechanism,plume_fre_g(:,ng),emiss_ajust_aer)

       !- split bburn emissions into flaming/smoldering parts
       IF(plumerise /= 0) &
          CALL emis_flam_smold(m1,m2,m3,isrctime,                   &
                              nsrc,chem1_src_g(:,:,:,ng),bburn,chem_nspecies, &
                              spc_chem_alloc,src,on,off,transport,  &
                              aer1_g(:,:,ng),aerosol,aer_nspecies,          &
                             
			      spc_aer_alloc,nmodes,aer_bburn,       &
                             
			      plume_mean_g(:,ng),plume_g(:,:,ng),tropical_forest, &
                              boreal_forest,savannah,grassland,     &
                              nveg_agreg,spc_aer_name,plume_fre_g(:,ng),plumerise)

       !-----
       !- to use the mass conservation fix
       !- convert from [kg/m^2]to mixing ratio expressed in [ppbm = 1.e9 kg/kg]
       !call convert_to_mixing_ratio(ng,m1,m2,m3,isrctime)

       !- convert from [kg/m^2]to tracer density expressed in [ 1.e9 kg/m3]
       CALL convert_to_tracer_density(m1,m2,m3,ia,iz,ja,jz,isrctime,        &
                                      nzpmax,dzt,rtgt,nsrc,                 &
                                      chem1_src_g(:,:,:,ng),bburn,chem_nspecies,  &
                                      spc_chem_alloc,src,off,transport,     &
                                      aer1_g(:,:,ng),aerosol,aer_nspecies,  &
                                      spc_aer_alloc,nmodes,aer_bburn,geoge, &
                                      volcanoes,bioge,CO2,ISFCL,spc_aer_name) 

       !- for now volcanic emissions only works with SIMPLE aerosol model
       IF(AEROSOL == 1 .and. volcanoes == 1) THEN
	  CALL vert_dist_of_volcanic_emission(m1,m2,m3,ia,iz,ja,jz,isrctime,&
				      nzpmax,dzt,zt,zm,rtgt,topt,dn0,nsrc,&
				      chem1_src_g(:,:,:,ng),geoge,chem_nspecies,  &
				      spc_chem_alloc,src,off,transport, &
				      aer1_g(:,:,ng),aerosol,aer_nspecies,	&
				      spc_aer_alloc,nmodes,volc_mean_g(ng), v_ash)
       ENDIF
       !-----

       !- next time (don't change the position of this line)
       IF(isrctime==1) next_srcfile(ng) = next_srcfile(ng) + 1 

    ENDDO

    IF(mchnum==master_num) THEN
      !- update srctime for the next time of reading/update    
      srctime1	 =   src_times(next_srcfile(ng)-1,1)
      srctime2	 =   src_times(next_srcfile(ng)  ,1)
    END IF
    if (mchnum == master_num) print *, mchnum,'LFR-DBG - Comunicando srctime2 = ',srctime2; call flush(6)
    !CALL Broadcast(srctime1, master_num, "srctime1") !LFR
    CALL Broadcast(srctime2, master_num, "srctime2") 
    print *, mchnum, 'LFR-DBG - recebido srctime2 = ',srctime2, src_times(next_srcfile(ng)  ,1); call flush(6)
    !******* debug *****
    !write(*,fmt='(I3.3,1X,A,F12.6,1X,F12.6)') &
    !    mchnum,'LFR:srctime1 and 2',srctime1,srctime2
    !CALL flush(6)
    !*******************
    IF(ng==ngrids) src_swap = 1


    !PRINT*,'next_srcfile,srctime1-2=',next_srcfile(ng),srctime1,srctime2
    !CALL flush(6)


  END SUBROUTINE read_sourcemaps




  !----------------------------------------------------------------------
  SUBROUTINE src_file_inv(srcpref,iyear1,imonth1,idate1,itime1,ngrids, &
                          time,timmax,nsrc,diur_cycle,mchnum,master_num)
			  
    ! original
    CHARACTER(len=*) , INTENT(IN) :: srcpref
    INTEGER          , INTENT(IN) :: iyear1
    INTEGER          , INTENT(IN) :: imonth1
    INTEGER          , INTENT(IN) :: idate1
    INTEGER          , INTENT(IN) :: itime1
    INTEGER          , INTENT(IN) :: ngrids
    REAL             , INTENT(IN) :: time
    REAL             , INTENT(IN) :: timmax

    ! mem_chem1
    INTEGER , INTENT(IN) :: nsrc
    INTEGER , INTENT(IN) :: diur_cycle(nsrc)

    ! node_mod
    INTEGER, INTENT(IN) :: mchnum
    INTEGER, INTENT(IN) :: master_num

    INTEGER :: nc,nf,lnf,nvftot, ng,it,isrc,nfstart
    INTEGER :: inyear,inmonth,indate,inhour
    INTEGER :: index
    REAL    :: localInc
    LOGICAL :: there

    CHARACTER(len=256), DIMENSION(maxsrcfiles) :: fnames
    CHARACTER(len=256) :: vpref
    CHARACTER(len=14)  :: itotdate,itotdate_current
    CHARACTER(len=2)   :: cgrid

    INTEGER                :: lenFnames_src
    INTEGER                :: lenItot_src
    INTEGER                :: sizeCharVec
    INTEGER                :: sizeIntVec
    INTEGER                :: sizeRealVec
    INTEGER                :: lastChar
    INTEGER                :: ierr2
    CHARACTER(len=8)       :: c0, c1
    INTEGER,   ALLOCATABLE :: intVec(:)
    CHARACTER, ALLOCATABLE :: charVec(:)
    REAL,      ALLOCATABLE :: realVec(:)
    REAL(kind=r8)          :: secs_init,secs_src
    CHARACTER(len=*), PARAMETER :: h="**(src_file_inv)**"
!    REAL(kind=8) :: secs_init,secs_src ! CCATT-BRAMS 4.3

    CHARACTER(len=256) :: sVarName


    sizeIntVec     = 2*ngrids
    lenFnames_src  = LEN(fnames_src)
    lenItot_src    = LEN(itotdate_src)


    IF (mchnum==master_num) THEN

       DO ng=1,ngrids
  	
          fnames(1:maxsrcfiles)= 'XXXXXXXXXXXXXXXX'
    
          !Get abs seconds of run start
          CALL date_abs_secs2(iyear1,imonth1,idate1,itime1*100,secs_init)
	
          ! get the current time  
          CALL date_add_to(iyear1,imonth1,idate1,itime1*100  &
        	           ,time,'s',inyear,inmonth,indate,inhour)

          CALL date_make_big (inyear,inmonth,indate,inhour,itotdate_current)

          ! Go through src files and make inventory

          nc=LEN_TRIM(srcpref)
          nvftot=-1
          vpref=srcpref
    	             
          WRITE(cgrid,'(a1,i1)') 'g',ng

          !### TEMPORARY CALL SYSTEM ALTERNATIVE ##### RMF 
          !RMF: 
          !with this struct you can only select between daily or hourly emissions.
	
          index = 0
          there = .FALSE.
	
          SELECT CASE(MINVAL(diur_cycle))
          CASE(0)
             nvftot = CEILING(((timmax/3600))) + 1 ; localInc = 3600.   
          CASE(1)
             nvftot = CEILING(((timmax/3600))/24.) ; localInc = 86400.   
          END SELECT
	
          DO nf = 1, nvftot	
		
             CALL makefnam (sVarName,srcpref,0,inyear,inmonth,indate,inhour,'T',cgrid,'vfm')
	          print*,'Looking for sources files -->: ',TRIM(sVarName)
             INQUIRE(file=TRIM(sVarName),exist=there)
		
             IF (there) THEN  
                index = index + 1 
                fnames(index) = TRIM(sVarName)           		
             ENDIF
		
             IF(nvftot > maxsrcfiles .OR. index > maxsrcfiles) THEN
                CALL fatal_error('Too many sources files')
             ENDIF
			
             CALL date_add_to(inyear,inmonth,indate,inhour  &
                              ,localInc,'s',inyear,inmonth,indate,inhour)	
          END DO
  
          nvftot = index

          IF (nvftot .EQ. 0) THEN
             CALL fatal_error('Sources files not found!')
          ENDIF
  
          !### TEMPORARY CALL SYSTEM ALTERNATIVE #####  RMF

          nsrcfiles(ng)=0
          DO nf=1,nvftot
             lnf=LEN_TRIM(fnames(nf))
             !print*,lnf,fnames(nf)

             READ(fnames(nf)(lnf-23:lnf-6),20) inyear,inmonth,indate,inhour
20           FORMAT(i4,1x,i2,1x,i2,1x,i6)

             CALL date_make_big(inyear,inmonth,indate,inhour,itotdate)

             nsrcfiles(ng)=nsrcfiles(ng)+1
             fnames_src(nsrcfiles(ng),ng)=fnames(nf)
             itotdate_src(nsrcfiles(ng),ng)=itotdate

             CALL date_abs_secs2(inyear,inmonth,indate,inhour,secs_src)
             src_times(nsrcfiles(ng),ng)=secs_src - secs_init

          ENDDO
       
          CALL RAMS_dintsort(nsrcfiles(ng),itotdate_src(:,ng),fnames_src(:,ng))

          !  start printing section
          !--------------------------------------------------------------
       
          PRINT*,' '
          PRINT*,' '
          PRINT*,' '
          PRINT*,'-------------------------------------------------------------'
          PRINT*,'-----------  Sources File Input Inventory --for --- GRID=', ng
          PRINT*,'-------------------------------------------------------------'
          DO nf=1,nsrcfiles(ng)
             PRINT 8,  nf, itotdate_src(nf,ng),src_times(nf,ng) ,TRIM(fnames_src(nf,ng))
          ENDDO
8         FORMAT(i4,1x,a16,1x,f10.0,2x,a)
          PRINT*,'------------------------------------------------------'
       
       ENDDO ! ngrids

       IF(ngrids > 1 .AND. INT((SUM(nsrcfiles(1:ngrids)))/ nsrcfiles(1)) .NE. ngrids) THEN
          CALL fatal_error("The number of src files for each grid, MUST be the same")
       ENDIF

       !- Are there enough src files available(for now only 1 grid is considered)

       ng=1
       IF(src_times(nsrcfiles(ng),ng) < timmax) THEN 
          PRINT*,'=============================================================================='
          PRINT*,'Warning:'
          PRINT*,'Not enough source files for the entire time integration were found, model will continue.'
          PRINT*,'=============================================================================='
       ENDIF

       !- perform some initializations (for now only 1 grid is considered)
       ng = 1 
       next_srcfile(1:ngrids) = 0
       DO nf=1,nvftot
          IF(itotdate_src(nf,ng) == itotdate_current) THEN
             next_srcfile(1:ngrids) = nf
             EXIT
          ELSEIF(itotdate_src(nf,ng) > itotdate_current) THEN
             next_srcfile(1:ngrids) = nf-1
             EXIT
          ENDIF
       ENDDO

       IF (next_srcfile(ng) < 1) THEN
            CALL fatal_error(' next_srcfile < 1 ')
       ENDIF

       IF (next_srcfile(ng) > nvftot-1 .AND. SUM(diur_cycle(1:nsrc)) < 4 ) THEN
            CALL fatal_error('next_srcfile > nvftot-1')
       ENDIF

       srctime1     =   src_times(next_srcfile(ng)  ,ng)
       srctime2     =   src_times(next_srcfile(ng)+1,ng)

       IF (next_srcfile(ng) > nvftot-1 .AND. SUM(diur_cycle(1:nsrc)) == 4 ) &
            srctime2 = srctime1 + 86400. 

       !- fill src_times arrays above nvftot with valid numbers (to use in 
       !- case of forecast or not more available data)
       DO ng=1,ngrids
          DO nf=nsrcfiles(ng)+1,maxsrcfiles
             src_times(nf,ng) = src_times(nf-1,ng) + (srctime2-srctime1)
          ENDDO
       ENDDO

       got_srcfiles_inv = .TRUE.

       ng=1
       PRINT*,'next_srcfile=',next_srcfile(ng),srctime1,srctime2,itotdate_current
       PRINT*,'next_srcfile= ',TRIM(fnames_src(next_srcfile(ng),ng))
       CALL flush(6)

    END IF ! IF (mchnum == master_num)

    ! allocate 'int' broadcast area
    ALLOCATE(intVec(sizeIntVec), stat=ierr2)
    IF (ierr2/=0) THEN
       WRITE(c0,"(i8)") ierr2
       WRITE(c1,"(i8)") sizeIntVec
       CALL fatal_error(h//" allocate intVec("//TRIM(ADJUSTL(c1))// &
            ") fails with stat="//TRIM(ADJUSTL(c0)))
    END IF

    ! master process gathers data for broadcasting
    
    IF (mchnum==master_num) THEN
       intVec(1:ngrids) = next_srcfile(1:ngrids)
       intVec(ngrids+1:sizeIntVec) = nsrcfiles(1:ngrids)
    END IF

    ! broadcast integer data to remaining processes
    ! envia: nsrcfiles(1:ngrids), 
    !        next_srcfile(1:ngrids)
    CALL Broadcast(intVec, master_num, "intVec")

    ! scatter broadcasted data
    next_srcfile(1:ngrids) = intVec(1:ngrids)
    nsrcfiles(1:ngrids) = intVec(ngrids+1:sizeIntVec)

    ! deallocate broadcast area
    DEALLOCATE(intVec, stat=ierr2)
    IF (ierr2/=0) THEN
       WRITE(c0,"(i8)") ierr2
       CALL fatal_error(h//" deallocate intVec fails with stat="//&
            TRIM(ADJUSTL(c0)))
    END IF

    ! allocate broadcast area
    sizeCharVec = SUM(nsrcfiles(1:ngrids)*(lenFnames_src + lenItot_src))
    ALLOCATE(charVec(sizeCharVec), stat=ierr2)
    IF (ierr2/=0) THEN
        WRITE(c0,"(i8)") ierr2
        WRITE(c1,"(i8)") sizeCharVec
        CALL fatal_error(h//" allocate charVec("//TRIM(ADJUSTL(c1))//&
             ") fails with stat="//TRIM(ADJUSTL(c0)))
     END IF

     ! master process prepares broadcast data

     IF (mchnum==master_num) THEN
        lastChar = 0
        DO ng=1,ngrids
           DO nf=1,nsrcfiles(ng)
              DO nc=1,lenFnames_src
                 charVec(lastChar+nc) = fnames_src(nf,ng)(nc:nc)
              END DO
              lastChar = lastChar + lenFnames_src
              DO nc=1,lenItot_src
                 charVec(lastChar+nc) = itotdate_src(nf,ng)(nc:nc)
              END DO
              lastChar = lastChar + lenItot_src
           END DO
        END DO
     END IF

     ! broadcast character data to remaining processes
     ! envia fnames_src(1:nsrcfiles(ng),1:ngrids),  
     !       itotdate_src(1:nsrcfiles(ng),1:ngrids)
     CALL Broadcast(charVec, master_num, "charVec")

     ! scatter broadcasted data
     lastChar=0
     DO ng=1,ngrids
        DO nf=1,nsrcfiles(ng)
           DO nc=1,lenFnames_src
              fnames_src(nf,ng)(nc:nc) = charVec(lastChar+nc)
           END DO
           lastChar = lastChar + lenFnames_src
           DO nc=1,lenItot_src
              itotdate_src(nf,ng)(nc:nc) = charVec(lastChar+nc)
           END DO
           lastChar = lastChar + lenItot_src
        END DO
     END DO

     ! deallocate broadcast area
     DEALLOCATE(charVec, stat=ierr2)
     IF (ierr2/=0) THEN
        WRITE(c0,"(i8)") ierr2
        CALL fatal_error(h//" deallocate charVec fails with stat="//&
             TRIM(ADJUSTL(c0)))
     END IF
     got_srcfiles_inv = .TRUE.   !Modificacao do Massaru


  END SUBROUTINE src_file_inv
  




  !--------------------------------------------------------------

  SUBROUTINE init_actual_time_index(nsrc,ntimes_src)

    INTEGER , INTENT(IN) :: nsrc
    INTEGER , INTENT(IN) :: ntimes_src(nsrc)
    

    INTEGER :: it

    !- index to control memory access of src arrays, because some
    !- arrays have 2 time levels and others only 1 time level
    !- for time 1, the memory position is always allocated for all sources types
    it = 1
    actual_time_index(it,1:nsrc) = 1
    !- for the second, will depend on if linterp is wanted or not.
    !- if linterp is not desired for any sources, actual_time_index = 1
    !- and will have actually one only memory position 
    it = 2
    actual_time_index(it,1:nsrc) = ntimes_src(1:nsrc) 

  END SUBROUTINE init_actual_time_index




  !--------------------------------------------------------------

  SUBROUTINE check_src_files(next_srcfile,nsrcfiles,iwtdo)

    ! original
    INTEGER , INTENT(INOUT) :: next_srcfile
    INTEGER , INTENT(IN)    :: nsrcfiles
    INTEGER , INTENT(INOUT) :: iwtdo

    iwtdo = 1

    !-check if is not greater the max number defined
    IF(next_srcfile+1 > maxsrcfiles) THEN
       CALL fatal_error("next_srcfile(ng)+1 > maxsrcfiles")
    ENDIF
    !PRINT*,' next_srcfile(ng) , nsrcfiles(ng)',next_srcfile, nsrcfiles
    !CALL flush(6)

    IF(next_srcfile > nsrcfiles) THEN 
!    IF(next_srcfile(ng)+1 > nsrcfiles(ng)) THEN 

       !- situation 1 : stop model execution
       IF(def_proc_src(1:LEN_TRIM(def_proc_src)) == 'STOP' ) THEN
          CALL fatal_error("Not src files available!")

       !- situation 2 :keep the current sources
       ELSEIF(def_proc_src(1:LEN_TRIM(def_proc_src)) == 'LAST_SOURCES' )THEN

	       iwtdo=0

	       next_srcfile=next_srcfile+1

          !- update srctime for the next time of reading/update    
          srctime1	 =   src_times(next_srcfile-1,1)
          srctime2	 =   src_times(next_srcfile  ,1)

          !print*,'-----------------------------------------------------'
          PRINT*,'Not src files available:'
	        PRINT*,'using previous day sources, model will continue ...'
	        PRINT*,'srctime1 and 2=',srctime1,srctime2
	        CALL flush(6)

          RETURN
       ENDIF

    ENDIF


  END SUBROUTINE check_src_files




  !--------------------------------------------------------------

  SUBROUTINE swap_sources(m1,m2,m3,time,chem_nspecies,spc_chem_alloc, &
                          src,off,nsrc,nvert_src,chem1_src_g,bburn,geoge)

    ! original
    INTEGER , INTENT(IN) :: m1
    INTEGER , INTENT(IN) :: m2
    INTEGER , INTENT(IN) :: m3
    REAL    , INTENT(IN) :: time

    ! chem1_list
    INTEGER , INTENT(IN) :: chem_nspecies
    INTEGER , INTENT(IN) :: spc_chem_alloc(6,chem_nspecies)
    INTEGER , INTENT(IN) :: src
    INTEGER , INTENT(IN) :: off

    ! mem_chem1
    INTEGER             , INTENT(IN)    :: nsrc
    INTEGER             , INTENT(IN)    :: nvert_src(nsrc)
    TYPE(chem1_src_vars), INTENT(INOUT) :: chem1_src_g(2,nsrc,chem_nspecies)
    INTEGER             , INTENT(IN)    :: bburn,geoge

    INTEGER ispc,isrc
    INTEGER, PARAMETER :: it1=1, & ! time level 1
                          it2=2    ! time level 2

    DO isrc=1,nsrc  

       !- if the time level 2 uses the same memory allocation area
       !- of time level 1 => nothing to do
       IF(actual_time_index(it2,isrc) == 1) CYCLE
       IF(isrc == bburn) STOP 444
       IF(isrc == geoge) STOP 444

       !- else: make the swap     
       DO ispc=1,chem_nspecies
          IF(spc_chem_alloc(src,ispc) == off) CYCLE

          chem1_src_g(it1,isrc,ispc)%sc_src(1:nvert_src(isrc),1:m2,1:m3) = &
               chem1_src_g(it2,isrc,ispc)%sc_src(1:nvert_src(isrc),1:m2,1:m3) 
       ENDDO
    ENDDO

    !- aerosol section still need to be done

    PRINT*,'--> source swapped done at time (h)=',time/3600.


  END SUBROUTINE swap_sources





  !--------------------------------------------------------------

  SUBROUTINE read_sources_vfm(ng,m1,m2,m3,iyear,imon,iday,isrctime,fname, &
                              chem_nspecies,spc_chem_alloc,spc_chem_name, &
                              src,on,off,chemical_mechanism,emiss_ajust,  &
			      
                              co,aer_nspecies,spc_aer_alloc,spc_aer_name, &
                              
			      urban,nucle,accum,nsrc,chem1_src_g,src_name,&
                              
			      chemistry,aer1_g,nmodes,aerosol,plumerise,  &
                              nveg_agreg,plume_mean_g,volc_mean_g,volcanoes, &
                              mchnum,master_num,mass_bin_dist,v_ash,aerosol_mechanism,&
			      plume_fre_g,emiss_ajust_aer)


  use mem_grid, only: grid_g

    ! original
    INTEGER       , INTENT(IN) :: ng
    INTEGER       , INTENT(IN) :: m1
    INTEGER       , INTENT(IN) :: m2
    INTEGER       , INTENT(IN) :: m3
    INTEGER       , INTENT(IN) :: iyear
    INTEGER       , INTENT(IN) :: imon
    INTEGER       , INTENT(IN) :: iday
    INTEGER       , INTENT(IN) :: isrctime
    CHARACTER*(*) , INTENT(IN) :: fname

    ! chem1_list
    INTEGER          , INTENT(IN) :: chem_nspecies
    INTEGER          , INTENT(IN) :: spc_chem_alloc(6,chem_nspecies)
    CHARACTER(LEN=8) , INTENT(IN) :: spc_chem_name(chem_nspecies)
    INTEGER          , INTENT(IN) :: src
    INTEGER          , INTENT(IN) :: on
    INTEGER          , INTENT(IN) :: off
    CHARACTER(LEN=24), INTENT(IN) :: chemical_mechanism
    REAL             , INTENT(IN) :: emiss_ajust(chem_nspecies)
    INTEGER          , INTENT(IN) :: CO

    ! aer1_list
    INTEGER          , INTENT(IN) :: aer_nspecies
    INTEGER          , INTENT(IN) :: spc_aer_alloc(6,nmodes,aer_nspecies)
    CHARACTER(LEN=8) , INTENT(IN) :: spc_aer_name(nmodes,aer_nspecies)
    INTEGER          , INTENT(IN) :: urban
    INTEGER          , INTENT(IN) :: nucle
    INTEGER          , INTENT(IN) :: accum
    INTEGER          , INTENT(IN) :: nmodes
    INTEGER          , INTENT(IN) :: v_ash
    REAL             , INTENT(INOUT) :: mass_bin_dist(nmodes)
    CHARACTER(LEN=24), INTENT(IN) :: aerosol_mechanism
    REAL             , INTENT(IN) :: emiss_ajust_aer(nmodes,aer_nspecies)

    ! mem_chem1
    INTEGER              , INTENT(IN)    :: nsrc
    TYPE(chem1_src_vars) , INTENT(INOUT) :: chem1_src_g(2,nsrc,chem_nspecies)
    CHARACTER(LEN=20)    , INTENT(IN)    :: src_name(nsrc)
    INTEGER              , INTENT(IN)    :: chemistry

    ! mem_aer1
    TYPE (aer1_vars) , INTENT(INOUT) :: aer1_g(nmodes,aer_nspecies)
    INTEGER          , INTENT(IN)    :: aerosol
    INTEGER          , INTENT(IN)    :: plumerise

    ! mem_plume_chem1
    INTEGER               , INTENT(IN)    :: nveg_agreg
    TYPE (plume_mean_vars), INTENT(INOUT) :: plume_mean_g(nveg_agreg)
    TYPE (plume_fre_vars) , INTENT(INOUT) :: plume_fre_g(5)

    ! mem_volc_chem1
    TYPE (volc_mean_vars), INTENT(INOUT) :: volc_mean_g
    INTEGER              , INTENT(IN)    :: volcanoes    

    ! node_mod
    INTEGER, INTENT(in) :: mchnum
    INTEGER, INTENT(in) :: master_num

    !- local var
    INTEGER :: iunit,ispc,isrc,iveg_agreg,nvert,ihour,itim
    CHARACTER(len=20) read_spc_name,read_src_name,section,read_units,date
    CHARACTER(len=20) read_mean,read_veg_name
    INTEGER read_ident_chem_mec,read_ident_src,dummy,read_ident_veg
    INTEGER read_ident_aer,read_aer_mode, imode
    INTEGER nspecies,nxp,nyp,i,j,ii,jj,iax,ipx,ivx
    REAL dep_glon(2), dep_glat(2)
!    REAL, ALLOCATABLE, DIMENSION(:,:) :: src_dummy_2d
    REAL, POINTER, DIMENSION(:,:) :: src_dummy_2d
    CHARACTER(len=32) :: chemical_mechanism_test,aerosol_mechanism_test
    LOGICAL :: there

    INTEGER :: exDo


    INTEGER                :: nc
    INTEGER                :: ierr2
    CHARACTER(len=8)       :: c0, c1, c2
    INTEGER                :: sizeIntVec
    INTEGER,   ALLOCATABLE :: intVec(:)
    REAL, POINTER          :: dummy_sc_src(:,:)
    INTEGER                :: lastChar
    INTEGER                :: sizeCharVec
    CHARACTER, ALLOCATABLE :: charVec(:)
    CHARACTER(len=*), PARAMETER :: h="**(read_sources_vfm)**"

    CHARACTER(len=64) :: cdummy
    CHARACTER(len=12) :: type_volc_process

      integer :: recn,recordLen

    !- initialization of local vars
    type_volc_process="XXXXXXXX"
    read_spc_name="XXXXXXXX"
    read_src_name="XXXXXXXX"
    section      ="XXXXXXXX"
    read_units   ="XXXXXXXX"
    date         ="XXXXXXXX"
    read_mean    ="XXXXXXXX"
    read_veg_name="XXXXXXXX"
    iax=0;ipx=0;ivx=0
    !- initial attributions/allocations
    nspecies = chem_nspecies + aer_nspecies*nmodes

    ihour = 0
    nvert = 1
    iunit = 2
    !- This routine does not allow negative fluxes.
    ALLOCATE (src_dummy_2d(m2,m3));src_dummy_2d=0.

!################### LFR
!  recordLen=4*m2*m3
!  open(unit=33,file='aer.gra',&
!      action='WRITE',status='REPLACE',form='UNFORMATTED',access='DIRECT', &
!      recl=recordLen)
!  recn=1
!
!
!################### LFR

    !-need to zerout aerosol sources, before reading the emission of the new day
    IF(AEROSOL > 0) THEN
       mass_bin_dist = 0.0
       DO ispc=1,aer_nspecies
          DO imode=1,nmodes
             IF(spc_aer_alloc(src,imode,ispc) == on) aer1_g(imode,ispc)%sc_src = 0.
          ENDDO
       ENDDO
    ENDIF

    if (mchnum==master_num) then
      print*,'---------------------------------------------------------------------------'
      print*,'opening emission file= ',fname(1:len_trim(fname))
      inquire(file=fname,exist=there)
      if(.not.there) then
         call fatal_error("emission file not found!")
      endif
      !- open the source file
      open(unit=iunit,file=fname,form='formatted',status='old') 
      ispc = 0 ; isrc = 0; iveg_agreg = 0
      read(iunit,*) nxp,(dep_glon(i),i=1,2)
      read(iunit,*) nyp,(dep_glat(i),i=1,2)
      read(iunit,*) date 
      write (*,fmt='(a)') '=== sources header (nxpoints,nypoint,lon,lat,date) ==='
      write (*,fmt='(4x,2(i4.4,1x),4(f8.3,1x),a)') &
                     nxp,nyp,(dep_glon(i),i=1,2),(dep_glat(i),i=1,2),date
      !- test if the source data is for the chemical mechanism that will be used:
      read(iunit,*)  chemical_mechanism_test,aerosol_mechanism_test
      if(trim( chemical_mechanism_test ) /=  trim(chemical_mechanism)) then
         call fatal_error("wrong chem mechanism at chem_sources. expected="// &
                          trim(chemical_mechanism(1:len_trim(chemical_mechanism)))//" read="// &
                          trim(chemical_mechanism_test(1:len_trim(chemical_mechanism_test))))
      else
         print*,'   chem mechanism= ',trim(chemical_mechanism(1:len_trim(chemical_mechanism)))
      endif
      if(trim( aerosol_mechanism_test ) /=  trim(aerosol_mechanism)) then
         call fatal_error("wrong aer mechanism at chem_sources. expected="// &
                          trim(aerosol_mechanism(1:len_trim(aerosol_mechanism)))//" read="// &
                          trim(aerosol_mechanism_test(1:len_trim(aerosol_mechanism_test))))
      else
         print*,'   aer  mechanism= ',trim(aerosol_mechanism(1:len_trim(aerosol_mechanism)))
         if(trim(aerosol_mechanism) == "matrix" ) then
           print*,"   level of matrix aer model is not checked in the read sources routine"
         endif
      endif
    end if ! (mchnum==master_num)

    do i=1,nspecies*nsrc*5
      if (mchnum==master_num) then
        read(iunit,*,iostat=exdo) section
        if(trim(section) == 'aerosol' .and. aerosol == off .and. iax==0) then
          print *,'warning: aerosol is present in emission file but turned off in ramsin.'
          print *,'the related emission fields will be ignored.'
          iax=1
        end if
        if( (trim(section) == 'plume' .or. trim(section) == 'plumefre') &
             .and. plumerise == 0 .and. ipx==0) then
          print *,'warning: plume information is present in emission file but plumerise' 
          print *,'is turned off in ramsin. the data will be ignored.'
          ipx=1
        end if
        if(trim(section) == 'volcanoes' .and. volcanoes /= on .and. ivx==0) then
          print *,'warning: volcanic emission data is present but turned off in ramsin.'
          print *, 'the data will be ignored.'
          ivx=0
        end if
      end if
      !
      call broadcast(exdo, master_num,"exdo")
      call broadcast(section, master_num, "section")
      !
      if (exdo<0) exit ! eof

      !emission section  --------------------------------

      !%%%%%%%%%%     chemistry section  %%%%%%%%%%%%%
      if(trim(section) == 'chemistry') then
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if (mchnum==master_num) then
          read(iunit,*)   read_spc_name &
                         , read_ident_chem_mec &
                         , read_src_name       &
                         , read_ident_src      &
                         , read_units
          ispc = read_ident_chem_mec
          isrc = read_ident_src
          itim = actual_time_index(isrctime,isrc)
          print*,"   chem spc/src: ", read_spc_name , read_src_name 
        endif
        sizeintvec=3
        !allocate 'int' broadcast area
        allocate(intvec(sizeintvec), stat=ierr2)
        if (ierr2/=0) then
          write(c0,"(i8)") ierr2
          write(c1,"(i8)") sizeintvec
          call fatal_error(h//" allocate intvec("//trim(adjustl(c1))// &
                                ") fails with stat="//trim(adjustl(c0)))
        endif
        ! master process gathers data for broadcasting   
        if (mchnum==master_num) then
          intvec(1) = ispc
          intvec(2) = isrc
          intvec(3) = itim
        endif
        ! broadcast integer data to remaining processes
        ! envia: ispc, isrc, itim
        call broadcast(intvec, master_num, "intvec")
        ! scatter broadcasted data
        ispc = intvec(1)
        isrc = intvec(2)
        itim = intvec(3)
        ! deallocate broadcast area
        deallocate(intvec, stat=ierr2)
        if (ierr2/=0) then
          write(c0,"(i8)") ierr2
          call fatal_error(h//" deallocate intvec fails with stat=" &
            //trim(adjustl(c0)))
        endif
        if(.not. associated( chem1_src_g(itim,isrc,ispc)%sc_src)) then
          call fatal_error("chem source memory not allocated for specie "//trim(read_spc_name)//&
                           " and source "//trim(read_src_name))
        endif
        allocate(dummy_sc_src(m2,m3), stat=ierr2)
        if (ierr2/=0) then
          write(c0,"(i8)") ierr2
          call fatal_error(h//" allocate dummy_sc_src fails with stat="//trim(adjustl(c0)))
        endif

        call readstorefullfieldandownchunk(ng,iunit,oneglobalemissdata(ng)%sc_src, &
                                             dummy_sc_src,                         &
                                             trim(spc_chem_name(ispc))//"_"//trim(src_name(isrc)))
!       call vfirec(iunit,chem1_src_g(itim,isrc,ispc)%sc_src(1,:,:),m2*m3,'lin')    
        where( dummy_sc_src(:,:) < 0.) dummy_sc_src(:,:)  = 0.0
        chem1_src_g(itim,isrc,ispc)%sc_src(1,1:m2,1:m3)=dummy_sc_src(1:m2,1:m3)
        if (mchnum==master_num) then 
          if (mchnum==master_num) write(*,fmt='("    Maxval,MinVal: ",2(E13.6,1X))') & 
            maxval(chem1_src_g(itim,isrc,ispc)%sc_src(1,1:m2,1:m3)),&
            minval(chem1_src_g(itim,isrc,ispc)%sc_src(1,1:m2,1:m3))
        endif
        deallocate(dummy_sc_src, stat=ierr2)
        if (ierr2/=0) then
          write(c0,"(i8)") ierr2
          call fatal_error(h//" deallocate dummy_sc_src fails with stat=" &
               //trim(adjustl(c0)))
        endif
        !- biogenic co only for tracer runs         
        if(chemistry > 0 .and.  ispc == co .and. spc_chem_alloc(3  ,co) == on ) then           
           chem1_src_g(itim,3,co)%sc_src(1,:,:) = 0.
        endif
        !
        !- ajust emissions 
        chem1_src_g(itim,isrc,ispc)%sc_src(1,:,:) = emiss_ajust(ispc) &
                *chem1_src_g(itim,isrc,ispc)%sc_src(1,:,:)
        !srf- especial para o barca - voc x 0.6 somente para urbano:
        !- ajustando emissoes antropicas
        !if(emiss_ajust(ispc) < 1. .and. trim(src_name(isrc))=='antro') then
        !     chem1_src_g(itim,isrc,ispc,ng)%sc_src(1,:,:) = emiss_ajust(ispc)*&
        !     chem1_src_g(itim,isrc,ispc,ng)%sc_src(1,:,:)
        !     print*,'barca voc urbanos=',trim(src_name(isrc)),emiss_ajust(ispc),'x',trim(read_spc_name)
        !endif
        !srf- especial para o barca - end
      
      !%%%%%%%%%%     aerosol section  %%%%%%%%%%%%%
      elseif(trim(section) == 'aerosol' ) then     
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if (mchnum==master_num) then
          read(iunit,*)   read_spc_name       &
                        , read_ident_aer      &
                        , read_aer_mode       &
                        , read_src_name       &
                        , read_ident_src      &
                        , read_units
          ispc = read_ident_aer
          imode= read_aer_mode
          isrc = read_ident_src
          print*,"   aer  spc/src: ",read_spc_name,read_src_name,ispc,imode    
      
          if(trim(read_spc_name)=='v_ash1' .and. ispc==v_ash) then
            print*,'-> reading size distr for ', trim(spc_aer_name(imode,ispc)), ' source: ' &
                  , trim(src_name(isrc))
            read(iunit,*)   cdummy
            read(iunit,*)   mass_bin_dist(:)
            print*,'-> reading volcanic mass_bin_dist = ',mass_bin_dist(:)
          endif
        endif
        sizeintvec=3
        ! allocate 'int' broadcast area
        allocate(intvec(sizeintvec), stat=ierr2)
        if (ierr2/=0) then
          write(c0,"(i8)") ierr2
          write(c1,"(i8)") sizeintvec
          call fatal_error(h//" allocate intvec("//trim(adjustl(c1)) &
               //") fails with stat="//trim(adjustl(c0)))
        endif
        ! master process gathers data for broadcasting   
        if (mchnum==master_num) then
          intvec(1) = ispc
          intvec(2) = imode
          intvec(3) = isrc
        endif
        ! broadcast integer data to remaining processes
        ! envia: ispc, imode, isrc
        call broadcast(intvec, master_num, "intvec")
        ! scatter broadcasted data
        ispc  = intvec(1)
        imode = intvec(2)
        isrc  = intvec(3)
        ! deallocate broadcast area
        deallocate(intvec, stat=ierr2)
        if (ierr2/=0) then
          write(c0,"(i8)") ierr2
          call fatal_error(h//" deallocate intvec fails with stat=" &
               //trim(adjustl(c0)))
        endif
        !- broadcast of the mass distribution of volcanic ash
        call broadcast(mass_bin_dist, master_num, "mass_bin_dist")
        write(c2,"(i8)") mchnum      
        call readstorefullfieldandownchunk(ng,iunit,oneglobalemissdata(ng)%src_dummy_2d &
                                          ,src_dummy_2d &
                                          ,trim(spc_aer_name(imode,ispc))//"_" &
                                          //trim(src_name(isrc)))

        !call vfirec(iunit,src_dummy_2d,m2*m3,'lin')
        if( aerosol > 0) then
          if (spc_aer_alloc(src,imode,ispc) == on) then
            where( src_dummy_2d(:,:) < 0.) src_dummy_2d(:,:)  = 0.
	    
!-19OCT2020 srf - check this later for MATRIX aer mechanism 
!           aer1_g(imode,ispc)%sc_src(1,1:m2,1:m3) = src_dummy_2d(1:m2,1:m3)
            aer1_g(imode,ispc)%sc_src(1,1:m2,1:m3) = emiss_ajust_aer(imode,ispc)*src_dummy_2d(1:m2,1:m3)
!
	    if (mchnum==master_num) write(*,fmt='("    Maxval,MinVal: ",2(E13.6,1X))') &
                maxval(aer1_g(imode,ispc)%sc_src(1,1:m2,1:m3)) &
               ,minval(aer1_g(imode,ispc)%sc_src(1,1:m2,1:m3))
!################### LFR
!            write(33,rec=recn) aer1_g(imode,ispc)%sc_src(1,1:m2,1:m3)
!            recn=recn+1
!################### LFR


          endif
       
        else
            call fatal_error('aer memory not allocated')
        endif  
      !%%%%%%%%%%     plume section for frp methodology %%%%%%%%%%%%%
      elseif(trim(section) == 'plumefre') then
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if(plumerise == 1.and. mchnum==master_num) then
          print*,"emission file prepared for plumerise 2 but plumerise flag in ramsin is not 2"
          call fatal_error('stop at routine read_sources 2')
        endif
        if (mchnum==master_num) then 
          print*,"   reading fre related products for plumerise model version 2"
          read(iunit,*)   read_mean      
        endif
        call broadcast(read_mean, master_num, "read_mean")
        if(trim(read_mean) == 'mean_fct' ) then     
          call readstorefullfieldandownchunk(ng,iunit,oneGlobalEmissData(ng)%flam_frac &
               ,src_dummy_2d,TRIM(read_mean))
          if(plumerise == 2) then
            where(src_dummy_2d < 0.0) src_dummy_2d  = 0.0
            plume_fre_g(iflam_frac)%pvar(1:m2,1:m3)=src_dummy_2d(1:m2,1:m3)
            if (mchnum==master_num) then 
              print*,"   frp max/min: ",read_mean,maxval(plume_fre_g(iflam_frac)%pvar) &
                                                 ,minval(plume_fre_g(iflam_frac)%pvar)
            endif
          endif
        elseif(trim(read_mean) == 'mean_frp' ) then     
          call readstorefullfieldandownchunk(ng,iunit,oneglobalemissdata(ng)%mean_frp &
                ,src_dummy_2d,trim(read_mean))
          if(plumerise == 2) then
            where(src_dummy_2d < 0.0) src_dummy_2d  = 0.0
            plume_fre_g(imean_frp)%pvar(1:m2,1:m3)=src_dummy_2d(1:m2,1:m3)*1.e+6 ! convert from mw to w
            if (mchnum==master_num) then 
              print*,"   frp max/min: ",read_mean,maxval(plume_fre_g(imean_frp)%pvar) &
                                                ,minval(plume_fre_g(imean_frp)%pvar)
            endif
          endif
        elseif(trim(read_mean) == 'std_frp' ) then
          call readstorefullfieldandownchunk(ng,iunit,oneglobalemissdata(ng)%std_frp &
                ,src_dummy_2d,trim(read_mean))
          if(plumerise == 2) then
            where(src_dummy_2d < 0.0) src_dummy_2d  = 0.0
            plume_fre_g(istd_frp)%pvar(1:m2,1:m3)=src_dummy_2d(1:m2,1:m3)*1.e+6! convert from mw to w 
            if (mchnum==master_num) then 
              print*,"   frp max/min: ",read_mean,maxval(plume_fre_g(istd_frp)%pvar) &
                                               ,minval(plume_fre_g(istd_frp)%pvar)
            endif
          endif  
        elseif(trim(read_mean) == 'mean_size' ) then
          call readstorefullfieldandownchunk(ng,iunit,oneglobalemissdata(ng)%mean_size &
               ,src_dummy_2d,trim(read_mean))
          if(plumerise == 2) then
            where(src_dummy_2d < 0.0) src_dummy_2d  = 0.0
            plume_fre_g(imean_size)%pvar(1:m2,1:m3)=src_dummy_2d(1:m2,1:m3)*1.e+6! convert from km2 to m2
            if (mchnum==master_num) then 
              print*,"   frp max/min: ",read_mean,maxval(plume_fre_g(imean_size)%pvar)& 
                                                 ,minval(plume_fre_g(imean_size)%pvar)
            endif
          endif
        elseif(trim(read_mean) == 'std_size' ) then
          call readstorefullfieldandownchunk(ng,iunit,oneglobalemissdata(ng)%std_size &
               ,src_dummy_2d,TRIM(read_mean))
          if(plumerise == 2) then
            where(src_dummy_2d < 0.0) src_dummy_2d  = 0.0
            plume_fre_g(istd_size)%pvar(1:m2,1:m3)=src_dummy_2d(1:m2,1:m3)*1.e+6! convert from km2 to m2
            if (mchnum==master_num) then 
              print*,"   frp max/min: ",read_mean,maxval(plume_fre_g(istd_size)%pvar) &
                                                 ,minval(plume_fre_g(istd_size)%pvar)
            endif
          endif
        else
          call fatal_error('unknow error in frp methodology')
        endif
      !%%%%%%%%%%     plume section   %%%%%%%%%%%%%
      elseif(trim(section) == 'plume') then
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if(plumerise == 2 .and. mchnum==master_num) then
          print*,"emission file prepared for plumerise 1 but plumerise flag in ramsin is not 1"
          call fatal_error('stop at routine read_sources 1')
        endif    
	      if (mchnum==master_num) then
          read(iunit,*)   read_mean &
                           , dummy &
                           , read_veg_name       &
                           , read_ident_veg      &
                           , read_units
          iveg_agreg = read_ident_veg       
        endif
        call broadcast(iveg_agreg, master_num, "iveg_agreg")
        ! allocate broadcast area
        sizecharvec = len(read_mean)+len(read_veg_name)
        allocate(charvec(sizecharvec), stat=ierr2)
        if (ierr2/=0) then
          write(c0,"(i8)") ierr2
          write(c1,"(i8)") sizecharvec
          call fatal_error(h//" allocate charvec("//trim(adjustl(c1))//&
                              ") fails with stat="//trim(adjustl(c0)))
        endif
        ! master process prepares broadcast data
        if (mchnum==master_num) then
          lastchar = 0
          do nc=1,len(read_mean)
            charvec(lastchar+nc) = read_mean(nc:nc)
          enddo
          lastchar = lastchar+len(read_mean)
          do nc=1,len(read_veg_name)
            charvec(lastchar+nc) = read_veg_name(nc:nc)
          enddo
        endif
        ! broadcast character data to remaining processes
        ! envia read_mean
        call broadcast(charvec, master_num, "charvec")
        ! scatter broadcasted data
        lastchar=0
        do nc=1,len(read_mean)
          read_mean(nc:nc) = charvec(lastchar+nc)
        end do
        lastchar = lastchar+len(read_mean)
        do nc=1,len(read_veg_name)
          read_veg_name(nc:nc) = charvec(lastchar+nc)
        end do
        ! deallocate broadcast area
        deallocate(charvec, stat=ierr2)
        if (ierr2/=0) then
          write(c0,"(i8)") ierr2
          call fatal_error(h//" deallocate charvec fails with stat="//&
                            trim(adjustl(c0)))
        end if
        if(trim(read_mean) == 'mean_fct' ) then
           call readstorefullfieldandownchunk(ng,iunit,oneglobalemissdata(ng)%flam_frac, &
                                                src_dummy_2d,                              &
!                                                plume_mean_g(iveg_agreg)%flam_frac,      &
                                                trim(read_mean)//"_"//trim(read_veg_name))
!             call vfirec(iunit,plume_mean_g(iveg_agreg)%flam_frac(:,:),m2*m3,'lin')

             if(plumerise == 1) then
           where(src_dummy_2d < 0.0) src_dummy_2d  = 0.0
                 plume_mean_g(iveg_agreg)%flam_frac(1:m2,1:m3)=src_dummy_2d(1:m2,1:m3)
             endif
    
    elseif(trim(read_mean) == 'firesize' ) then

             call readstorefullfieldandownchunk(ng,iunit,oneglobalemissdata(ng)%fire_size, &
                                                src_dummy_2d,                              &
!                                               plume_mean_g(iveg_agreg)%fire_size,      &
                                                trim(read_mean)//"_"//trim(read_veg_name))
!             call vfirec(iunit,plume_mean_g(iveg_agreg)%fire_size(:,:),m2*m3,'lin')
             if(plumerise == 1) then
           where(src_dummy_2d < 0.0) src_dummy_2d  = 0.0
                 plume_mean_g(iveg_agreg)%fire_size(1:m2,1:m3)=src_dummy_2d(1:m2,1:m3)
             endif
          
    else
                call fatal_error('model not yet prepared for flam-frac of each one of species')
          endif


       !- volcanoes section  --------------------------------
       elseif(   trim(section) == 'volcanic-eruption' .or. &
                 trim(section) == 'volcanic-degassing'     ) then

          if (mchnum==master_num) then
             read(iunit,*)  read_mean , &
                            read_units
             if(trim(section) == 'volcanic-eruption' ) type_volc_process="eruption"
             if(trim(section) == 'volcanic-degassing') type_volc_process="degassing"
          endif

          call broadcast(read_mean, master_num, "read_mean")
          call broadcast(read_units, master_num, "read_units")
          call broadcast(type_volc_process, master_num, "type_volc_process")

          !write(c2,"(i8)") mchnum
          !print *,'proc ',trim(adjustl(c2)),' volc: reading ',trim(read_mean),' units= ',trim(read_units)
          
          if(trim(read_mean) == 'inject_height' ) then

             call readstorefullfieldandownchunk(ng,iunit,oneglobalemissdata(ng)%plum_heigth, &
                                                src_dummy_2d,                                &
!                                               volc_mean_g%plum_heigth,                    &
                                                trim(read_mean)//"_"//trim(read_units))
!             call vfirec(iunit,volc_mean_g%plum_heigth(:,:),m2*m3,'lin')
             if(volcanoes == on) then
                
    where( src_dummy_2d(:,:) < 0.) src_dummy_2d(:,:)  = 0.
          volc_mean_g%plum_heigth(1:m2,1:m3)=src_dummy_2d(1:m2,1:m3)
                
    if (mchnum==master_num)print*,'max inject_height (m) =',maxval(oneglobalemissdata(ng)%plum_heigth)
             endif
          elseif(trim(read_mean) == 'vent_elevation' ) then   

             call readstorefullfieldandownchunk(ng,iunit,oneglobalemissdata(ng)%vent_elev, &
                                                src_dummy_2d,                              &
!                                               volc_mean_g%vent_elev,                    &
                                                trim(read_mean)//"_"//trim(read_units))
!             call vfirec(iunit,volc_mean_g%vent_elev(:,:),m2*m3,'lin')

             if(volcanoes == on) then                
    where( src_dummy_2d(:,:) < 0.) src_dummy_2d(:,:)  = 0.    
    volc_mean_g%vent_elev(1:m2,1:m3)=src_dummy_2d(1:m2,1:m3)
                
    if (mchnum==master_num)print*,'max vent_elevation (m) =',maxval(oneglobalemissdata(ng)%vent_elev(:,:))

             endif
          elseif(trim(read_mean) == 'time_duration' ) then   

             call readstorefullfieldandownchunk(ng,iunit,oneglobalemissdata(ng)%duration, &
                                                src_dummy_2d,                             &
!                                               volc_mean_g%duration,                    &
                                                trim(read_mean)//"_"//trim(read_units))
!             call vfirec(iunit,volc_mean_g%vent_elev(:,:),m2*m3,'lin')

           
             if(volcanoes == on) then                
    where( src_dummy_2d(:,:) < 0.) src_dummy_2d(:,:)  = 0.    
       
          volc_mean_g%duration(1:m2,1:m3)=src_dummy_2d(1:m2,1:m3)
                if (mchnum==master_num)print*,'max time duration (m) =',maxval(oneglobalemissdata(ng)%duration(:,:))
             endif
          else
             call fatal_error('volc sources not prepared yet for: '//trim(read_mean))
          endif

       else 
          write (*,fmt='(a)') 'read_mean="'//trim(read_mean)//'"'
          call fatal_error('unknown error!')
       endif

    enddo

100 continue   
    
    if (mchnum==master_num) then 
       print*,'---------------------------------------------------------------------------'
       print*,'---------------------------------------------------------------------------'
      close (iunit) 
    endif

    !- volcanic ash mass distribution:
    IF(AEROSOL == 1 .and. trim(type_volc_process)=="eruption") THEN
       !- inverting the array mass_bin_dist:
       !- bin 1 will be the smallest in size, bin 10 will be the bigger
       !
       !!do i=1,int(nbins/2)
       !!  rdummy                          =mass_bin_dist(nbins-(i-1))
       !!  mass_bin_dist(nbins-(i-1))=mass_bin_dist(i          )
       !!  mass_bin_dist(i          )=rdummy
       !!enddo 
       !
       ! 2nd way
       mass_bin_dist(:)=mass_bin_dist(nmodes:1:-1)
       print*, ' volc ash total mass = ', maxval(aer1_g(1,v_ash)%sc_src(1,:,:)),&
                mass_bin_dist(:); call flush(6)
       
       !
       !
       ! -- assuming the src array of ash bin = 1 (V_ASH1) contains the _total_ ash mass
       do imode=2,nmodes
         aer1_g(imode,v_ash)%sc_src(1,:,:)= mass_bin_dist(imode)*aer1_g(1,v_ash)%sc_src(1,:,:)
                
       enddo
       ! - corrected ash mass for bin 1
       aer1_g(1,v_ash)%sc_src(1,:,:)= mass_bin_dist(1)*aer1_g(1,v_ash)%sc_src(1,:,:)
       
       !print*," volc ash=",v_ash&
       !       ,1,mass_bin_dist(1),maxval(aer1_g(1,v_ash)%sc_src(1,:,:))

! -alterando a distribuicao de concentracao de cinzas 29/12/2013
! -aplicado ao caso do vulcano chileno PYHUE
!       do imode=1,nmodes
!         if(imode .le. 5) then
!	   aer1_g(imode,v_ash)%sc_src(1,:,:)=aer1_g(imode,v_ash)%sc_src(1,:,:)/60.
!	  Else
!	   aer1_g(imode,v_ash)%sc_src(1,:,:)=aer1_g(imode,v_ash)%sc_src(1,:,:)*1.65
!         ENDIF         
!       enddo
!
    ENDIF


    IF(trim(type_volc_process)=="degassing" .and. volcanoes == on) THEN
     ! converting the plume_heigth of degass volcanic from meters above sea level
     ! to meters above the vent, to allow an unified treatment with the eruption
     ! volcanoes at the vertical mass distribution routine:
     volc_mean_g%plum_heigth(:,:) = volc_mean_g%plum_heigth(:,:) &
                                  - volc_mean_g%vent_elev  (:,:)
    ENDIF 

!--(DMK-BRAMS-5.0-INI)-----------------------------------------------------------------------------------
!    IF(AEROSOL > 0) THEN
!       !- special section for aerosol sulfate 
!       !- mode Aitken = 50% mode Accumulation
!       !* ref: Stier et al., The aerosol-climate model ECHAM5-HAM. Atmos.
!       !  Chem. Phys., 5,1125-1156,2005.
!       IF(spc_aer_alloc(src,nucle,urban) == on) &
!            aer1_g(nucle,urban)%sc_src(1,:,:)=0.5*aer1_g(accum,urban)%sc_src(1,:,:)
!
!       IF(spc_aer_alloc(src,accum,urban) == on) &
!            aer1_g(accum,urban)%sc_src(1,:,:)=0.5*aer1_g(accum,urban)%sc_src(1,:,:)
!
!    ENDIF
!--(DMK-BRAMS-5.0-INI)-----------------------------------------------------------------------------------

    DEALLOCATE (src_dummy_2d)


!################### LFR
!
!  close(33)
!
!  open(unit=33,file='aer.ctl' &
!       ,action='WRITE',status='replace',form='FORMATTED')
!
!  !writing the name of grads file
!  write(33,*) 'dset ^aer.gra'
!  !writing others infos to ctl
!  write(33,*) 'undef -0.9990000E+34'
!  write(33,*) 'title AerTest'
!  write(33,*) 'xdef ',m2,' linear ',grid_g(1)%glon(1,1),grid_g(1)%glon(2,1)-grid_g(1)%glon(1,1)
!  write(33,*) 'ydef ',m3,' linear ',grid_g(1)%glat(1,1),grid_g(1)%glat(1,2)-grid_g(1)%glat(1,1)
!  write(33,*) 'zdef ',1,'levels',1000
!  write(33,*) 'tdef 1 linear 00:00z01jan2018     1mo'
!  write(33,*) 'vars ',4
!  write(33,*) 'urban2 1 99 urban2'
!  write(33,*) 'urban3 1 99 urban3' 
!  write(33,*) 'bburn2 1 99 bburn2'
!  write(33,*) 'bburn3 1 99 bburn3'
!  write(33,*) 'endvars'!
!
!  close(33)
!
!################### LFR   


END SUBROUTINE read_sources_vfm

  !-------------------------------------------------------------
  ! (DMK) NOT USED
  !
  ! subroutine convert_to_mixing_ratio(ng,m1,m2,m3,isrctime) 
  !  use mem_basic, only: basic_g
  !  use mem_grid, only: dzt,grid_g
  !  use mem_chem1 
  !  use chem1_list, only : chem_nspecies=>nspecies &
  !                        ,spc_chem_alloc=>spc_alloc&
  !		          ,spc_chem_name =>spc_name,src,on,off,offline&
  !			  ,transport
  !
  !  use mem_aer1
  !  use aer1_list, only : aer_nspecies=>nspecies &
  !                       ,spc_aer_alloc=>spc_alloc, nmodes&
  !                       ,spc_aer_name =>aer_name &
  !			 ,aer_bburn => bburn 
  !  implicit none
  !
  !  integer,intent(IN) :: ng,m1,m2,m3,isrctime
  !  
  !  
  !  ! local var
  !  !Fator de conversao de unidades    
  !  real,parameter :: fcu =1.e+9 !=> ppbm 
  !  !LPCE
  !  real :: rhodz_inv,dz
  !  integer :: i,j,ksrc,isrc,ispc,imode,itim
  !  
  !  ksrc=2 !surface level of emission in the model
  !  
  !  do j=1,m3
  !     do i=1,m2
  !        
  !        ! Todas as unidades estao em kg/m2/dia => use 'dz' em vez de 'vol'
  !        !    vol = 1./(dxt(i,j)*dyt(i,j)*dzt(k))*rtgt(i,j)
  !
  !        dz        = grid_g(ng)%rtgt(i,j)/dzt(ksrc) ! dzt=1/(z(k)-z(k-1))	
  !	  rhodz_inv = 1./(dz*basic_g(ng)%dn0(ksrc,i,j))
  !	  do ispc=1,chem_nspecies
  !           if(spc_chem_alloc(src,ispc) == off) cycle
  !
  !                !- convert from kg/m^2  to  density (kg[gas]/m^3*1.e9)
  ! 	     	  do isrc=1,nsrc
  !
  !		      itim = actual_time_index(isrctime,isrc)
  !
  !	              chem1_src_g(itim,isrc,ispc,ng)%sc_src(1,i,j) = &
  !		      chem1_src_g(itim,isrc,ispc,ng)%sc_src(1,i,j) * fcu * rhodz_inv
  !                enddo
  !                ! copy smoldering emission from  bburn to ksrc
  !		  itim = actual_time_index(isrctime,bburn)
  !                chem1_src_g(itim,bburn,ispc,ng)%sc_src(ksrc,i,j) = &
  !		  chem1_src_g(itim,bburn,ispc,ng)%sc_src(1   ,i,j)
  !
  !	  enddo
  !	enddo 
  !  enddo
  !  
  ! !- aerosol section 
  ! IF(AEROSOL > 0 ) then
  !  do j=1,m3
  !     do i=1,m2
  !        
  !        ! Todas as unidades estao em kg/m2/dia => use 'dz' em vez de 'vol'
  !        !    vol = 1./(dxt(i,j)*dyt(i,j)*dzt(k))*rtgt(i,j)
  !
  !        dz        = grid_g(ng)%rtgt(i,j)/dzt(ksrc) ! dzt=1/(z(k)-z(k-1))	
  !	  rhodz_inv = 1./(dz*basic_g(ng)%dn0(ksrc,i,j))
  ! 	  do ispc=1,aer_nspecies
  !	    do imode=1,nmodes
  !	     
  !	      if(spc_aer_alloc(src,imode,ispc)       == off .or. &
  !               spc_aer_alloc(transport,imode,ispc) == off) cycle
  !            
  !	      !- convert from kg/m^2  to  density (kg[aer]/m^3*1.e9)
  !           
  !	     
  !	      aer1_g(imode,ispc,ng)%sc_src(1,i,j) = &
  !	      aer1_g(imode,ispc,ng)%sc_src(1,i,j) * fcu * rhodz_inv
  !          
  !          enddo
  !        enddo
  !
  !	 
  !	  ! copy smoldering emission from aerosol bburn to ksrc
  !	  do imode=1,nmodes
  !	     
  !	     if(spc_aer_alloc(src,imode,aer_bburn)       == off .or. &
  !               spc_aer_alloc(transport,imode,aer_bburn) == off) cycle
  !         
  !	     aer1_g(imode,aer_bburn,ng)%sc_src(ksrc,i,j) = &
  !	     aer1_g(imode,aer_bburn,ng)%sc_src(1   ,i,j)
  !        enddo
  !
  !	enddo 
  !  enddo
  ! ENDIF
  !
  ! end subroutine convert_to_mixing_ratio
  !-------------------------------------------------------------

  SUBROUTINE convert_to_tracer_density(m1,m2,m3,ia,iz,ja,jz,isrctime,      &
                                       nzpmax,dzt,rtgt,nsrc,chem1_src_g,   &
                                       bburn,chem_nspecies,spc_chem_alloc, &
                                       src,off,transport,aer1_g,aerosol,   &
                                       aer_nspecies,spc_aer_alloc,nmodes,  &
                                       aer_bburn,geoge,volcanoes,bioge,CO2,&
				       ISFCL,spc_aer_name) 

    ! original
    INTEGER              , INTENT(IN)    :: bioge
    INTEGER , INTENT(IN) :: m1
    INTEGER , INTENT(IN) :: m2
    INTEGER , INTENT(IN) :: m3
    INTEGER , INTENT(IN) :: ia
    INTEGER , INTENT(IN) :: iz
    INTEGER , INTENT(IN) :: ja
    INTEGER , INTENT(IN) :: jz
    INTEGER , INTENT(IN) :: isrctime
    INTEGER , INTENT(IN) :: isfcl,CO2  !DSM

    ! grid_dims
    INTEGER , INTENT(IN) :: nzpmax

    ! mem_grid
    REAL    , INTENT(IN) :: dzt(nzpmax)
    REAL    , INTENT(IN) :: rtgt(m2,m3)

    ! chem1_list
    INTEGER , INTENT(IN) :: chem_nspecies
    INTEGER , INTENT(IN) :: spc_chem_alloc(6,chem_nspecies)
    INTEGER , INTENT(IN) :: src
    INTEGER , INTENT(IN) :: off
    INTEGER , INTENT(IN) :: transport

    ! aer1_list
    INTEGER , INTENT(IN) :: aer_nspecies
    INTEGER , INTENT(IN) :: spc_aer_alloc(6,nmodes,aer_nspecies)
    INTEGER , INTENT(IN) :: nmodes
    INTEGER , INTENT(IN) :: aer_bburn
    CHARACTER(LEN=8) , INTENT(IN) :: spc_aer_name(nmodes,aer_nspecies)

    ! mem_chem1
    INTEGER              , INTENT(IN)    :: nsrc
    TYPE(chem1_src_vars) , INTENT(INOUT) :: chem1_src_g(2,nsrc,chem_nspecies)
    INTEGER              , INTENT(IN)    :: bburn,geoge

    ! mem_aer1
    TYPE (aer1_vars) , INTENT(INOUT) :: aer1_g(nmodes,aer_nspecies)
    INTEGER          , INTENT(IN)    :: aerosol

    ! mem_volc_chem1
    INTEGER          , INTENT(IN)    :: volcanoes    


    ! local var
    !Fator de conversao de unidades    
    REAL,PARAMETER :: fcu =1.e+9 !=> ppbm 
    !LPCE
    REAL :: dz_inv,dz
    INTEGER :: i,j,ksrc,isrc,ispc,imode,itim

    ksrc=2 !surface level of emission in the model

    !- chemistry section

    DO j=ja,jz
       DO i=ia,iz

          ! Todas as unidades estao em kg/m2/dia => use 'dz' em vez de 'vol'
          !    vol = 1./(dxt(i,j)*dyt(i,j)*dzt(k))*rtgt(i,j)

!         dz        = grid_g(ng)%rtgt(i,j)/dzt(ksrc) ! dzt=1/(z(k)-z(k-1))	
          dz        = rtgt(i,j)/dzt(ksrc) ! dzt=1/(z(k)-z(k-1))	
	  dz_inv    = 1./dz
	  DO ispc=1,chem_nspecies
	  
         IF(spc_chem_alloc(src,ispc) == off) CYCLE

         !- convert from kg/m^2  to  density (kg[gas]/m^3*1.e9)
         DO isrc=1,nsrc
	       
	   !- only call geoge emissions if volcanoes is ON.
	    IF(isrc==geoge .AND. volcanoes == off) CYCLE

            if(ISFCL == 5 .and. ispc == CO2 .and. isrc == bioge) cycle  !DSM e Saulo

            itim = actual_time_index(isrctime,isrc)
            chem1_src_g(itim,isrc,ispc)%sc_src(1,i,j) = &
            chem1_src_g(itim,isrc,ispc)%sc_src(1,i,j) * fcu * dz_inv
         ENDDO
         ! copy smoldering emission from  bburn to ksrc
         itim = actual_time_index(isrctime,bburn)

         chem1_src_g(itim,bburn,ispc)%sc_src(ksrc,i,j) = &
         chem1_src_g(itim,bburn,ispc)%sc_src(1   ,i,j)

     ENDDO
    ENDDO
    ENDDO

    !- aerosol section 
    IF(AEROSOL > 0 ) THEN
       DO j=ja,jz
          DO i=ia,iz

             ! Todas as unidades estao em kg/m2/dia => use 'dz' em vez de 'vol'
             !    vol = 1./(dxt(i,j)*dyt(i,j)*dzt(k))*rtgt(i,j)

!            dz        = grid_g(ng)%rtgt(i,j)/dzt(ksrc) ! dzt=1/(z(k)-z(k-1))	
             dz        = rtgt(i,j)/dzt(ksrc) ! dzt=1/(z(k)-z(k-1))	
             dz_inv    = 1./dz
             DO ispc=1,aer_nspecies
                DO imode=1,nmodes

                   IF(spc_aer_alloc(src      ,imode,ispc) == off .OR. &
                      spc_aer_alloc(transport,imode,ispc) == off) CYCLE

                   !- convert from kg/m^2  to  density (kg[aer]/m^3*1.e9)

                   aer1_g(imode,ispc)%sc_src(1,i,j) = &
                   aer1_g(imode,ispc)%sc_src(1,i,j) * fcu * dz_inv



                ENDDO
             ENDDO


             ! copy smoldering emission from aerosol bburn to ksrc
             IF(AEROSOL==1) THEN
	     
	       DO imode=1,nmodes

                  IF(spc_aer_alloc(src      ,imode,aer_bburn) == off .OR. &
                     spc_aer_alloc(transport,imode,aer_bburn) == off) CYCLE

                     aer1_g(imode,aer_bburn)%sc_src(ksrc,i,j) = aer1_g(imode,aer_bburn)%sc_src(1,i,j)
               ENDDO
	     ELSEIF(AEROSOL==2) THEN ! For MATRIX
              DO ispc=1,aer_nspecies
	        DO imode=1,nmodes

                  IF(spc_aer_alloc(src      ,imode,ispc) == off .OR. &
                     spc_aer_alloc(transport,imode,ispc) == off) CYCLE

	          IF(spc_aer_name(imode,ispc)=="boc_bcar" .or. &
	             spc_aer_name(imode,ispc)=="boc_ocar" )    THEN
                      aer1_g(imode,ispc)%sc_src(ksrc,i,j) = aer1_g(imode,ispc)%sc_src(1,i,j)
                  ENDIF
	        ENDDO
              ENDDO
	     ENDIF

          ENDDO
       ENDDO
    ENDIF

  END SUBROUTINE convert_to_tracer_density

  !----------------------------------------------------------------------

  SUBROUTINE  emis_flam_smold(n1,n2,n3,isrctime,                 &
                              nsrc,chem1_src_g,bburn,chem_nspecies, &
                              spc_chem_alloc,src,on,off,transport,  &
                              aer1_g,aerosol,aer_nspecies,          &
                              spc_aer_alloc,nmodes,aer_bburn,       &
                              plume_mean_g,plume_g,tropical_forest, &
                              boreal_forest,savannah,grassland,     &
                              nveg_agreg,spc_aer_name,plume_fre_g,plumerise)

    ! original
    INTEGER , INTENT(IN) :: n1
    INTEGER , INTENT(IN) :: n2
    INTEGER , INTENT(IN) :: n3
    INTEGER , INTENT(IN) :: isrctime,plumerise

    ! chem1_list
    INTEGER          , INTENT(IN) :: chem_nspecies
    INTEGER          , INTENT(IN) :: spc_chem_alloc(6,chem_nspecies)
    INTEGER          , INTENT(IN) :: src
    INTEGER          , INTENT(IN) :: on
    INTEGER          , INTENT(IN) :: off
    INTEGER          , INTENT(IN) :: transport

    ! aer1_list
    INTEGER          , INTENT(IN) :: nmodes
    INTEGER          , INTENT(IN) :: aer_nspecies
    INTEGER          , INTENT(IN) :: spc_aer_alloc(6,nmodes,aer_nspecies)
    INTEGER          , INTENT(IN) :: aer_bburn

    ! mem_chem1
    INTEGER              , INTENT(IN)    :: nsrc
    TYPE(chem1_src_vars) , INTENT(INOUT) :: chem1_src_g(2,nsrc,chem_nspecies)
    INTEGER              , INTENT(IN)    :: bburn


    ! mem_aer1
    TYPE (aer1_vars) , INTENT(INOUT) :: aer1_g(nmodes,aer_nspecies)
    INTEGER          , INTENT(IN)    :: aerosol

    ! mem_plume_chem1
    INTEGER               , INTENT(IN)    :: nveg_agreg
    TYPE (plume_mean_vars), INTENT(INOUT) :: plume_mean_g(nveg_agreg)
    TYPE (plume_fre_vars) , INTENT(INOUT) :: plume_fre_g(5)
    TYPE (plume_vars)     , INTENT(INOUT) :: plume_g(nveg_agreg,chem_nspecies)
    INTEGER               , INTENT(IN)    :: tropical_forest
    INTEGER               , INTENT(IN)    :: boreal_forest
    INTEGER               , INTENT(IN)    :: savannah
    INTEGER               , INTENT(IN)    :: grassland
    CHARACTER(LEN=8) , INTENT(IN) :: spc_aer_name(nmodes,aer_nspecies)

    REAL,DIMENSION(n2,n3) :: smold_frac 
    INTEGER iv,ispc,i,j,imode,itim
    INTEGER:: imean_plume
    imean_plume = 1 !change this at alloc_plume_chem1 routine also


    !- time index of memory allocation position        
    itim = actual_time_index(isrctime,bburn)
    IF(itim > 1) THEN
        CALL fatal_error('Time level 2 not allowed when plumerise is used!')
    ENDIF

    IF(imean_plume == on) THEN
       !-----  
       !- calcula a emissao smoldering e fatores para obtencao da fracao
       !- flaming em funcao da emissao smoldering
       if(plumerise == 1) then
         smold_frac(1:n2,1:n3) = 1.- ( plume_mean_g(tropical_forest)%flam_frac(1:n2,1:n3) + &
                                       plume_mean_g(boreal_forest)%flam_frac(1:n2,1:n3) + &
                                       plume_mean_g(savannah     )%flam_frac(1:n2,1:n3) + &
                                       plume_mean_g(grassland    )%flam_frac(1:n2,1:n3)	)    
       elseif(plumerise == 2) then
       
         smold_frac(1:n2,1:n3) = 1.- plume_fre_g(iflam_frac)%pvar(1:n2,1:n3)
       
       endif

       !- chemistry section (only for bburn source)
       DO ispc = 1,chem_nspecies

  	  IF(spc_chem_alloc(src,ispc) /= on) CYCLE 

          !- convert from 'total' emisson to 'smoldering' part
  	  chem1_src_g(itim,bburn,ispc)%sc_src(1,:,:) = smold_frac(:,:) * &  
                                                       chem1_src_g(itim,bburn,ispc)%sc_src(1,:,:)
       ENDDO
       IF(AEROSOL == 1 ) THEN
          !- aerosol section (only for bburn aerosols)
          DO imode=1,nmodes

             IF(spc_aer_alloc(src,imode,aer_bburn)       == off .OR. &
                spc_aer_alloc(transport,imode,aer_bburn) == off) CYCLE

             !- convert from 'total' emisson to 'smoldering' part
             aer1_g(imode,aer_bburn)%sc_src(1,:,:) = smold_frac(:,:) * & 
                                                     aer1_g(imode,aer_bburn)%sc_src(1,:,:)

          ENDDO
       !- this for MATRIX
       ELSEIF(AEROSOL ==  2 ) THEN
          !- aerosol section (only for bburn aerosols)
          DO ispc = 1,aer_nspecies
	  
	    DO imode=1,nmodes

             IF(spc_aer_alloc(src,imode,ispc)       == off .OR. &
                spc_aer_alloc(transport,imode,ispc) == off) CYCLE
             
	     !-only for bburn aerosols) 
	     IF(spc_aer_name(imode,ispc)=="boc_bcar" .or. &
	        spc_aer_name(imode,ispc)=="boc_ocar" )    THEN
	     !print*,"bburn=",spc_aer_name(imode,ispc),maxval(aer1_g(imode,ispc)%sc_src(1,:,:))
             
	     !- convert from 'total' emisson to 'smoldering'
                aer1_g(imode,ispc)%sc_src(1,:,:) = smold_frac(:,:) * & 
                                                   aer1_g(imode,ispc)%sc_src(1,:,:)
             ENDIF
            ENDDO
          ENDDO
       ENDIF

      !- convert from flaming fraction to relationship with phase smoldering emission
      if(plumerise == 1) then
         do iv = 1, nveg_agreg
           plume_mean_g(iv)%flam_frac(:,:) = plume_mean_g(iv)%flam_frac(:,:)/ &
                                            (1.e-8+smold_frac(:,:))
         enddo
       elseif(plumerise == 2) then
       
	 plume_fre_g(iflam_frac)%pvar(:,:) = plume_fre_g(iflam_frac)%pvar(:,:)/ &
                                            (1.e-8+smold_frac(:,:))       
       endif
       !-----

    ELSE

       !-----      case where each specie has his own flaming fraction ----------------
       CALL fatal_error('Aerosol emission not ready for this option!')

       DO ispc = 1,chem_nspecies
          IF(spc_chem_alloc(src,ispc) /= on) CYCLE 
          smold_frac(:,:) = 1.- ( plume_g(tropical_forest,ispc)%fct(:,:) + &
                                  plume_g(boreal_forest,ispc)%fct(:,:) + &
                                  plume_g(savannah	    ,ispc)%fct(:,:) + &
                                  plume_g(grassland    ,ispc)%fct(:,:)   )

          !- convert from 'total' emisson to 'smoldering'
          chem1_src_g(itim,bburn,ispc)%sc_src(1,:,:) = smold_frac(:,:) * &  
                                                       chem1_src_g(itim,bburn,ispc)%sc_src(1,:,:)

          !- convert from flaming fraction to relationship with phase smoldering emission
          DO iv = 1, nveg_agreg
             plume_g(iv,ispc)%fct(:,:) = plume_g(iv,ispc)%fct(:,:)/ &
                                    (1.e-8+smold_frac(:,:))
             !- flamming emission =  plume_g(iv,iscp,ng)%fct(:,:) * &
             !		             chem1_src_g(itim,bburn,ispc,ng)%sc_src(1,:,:)
          ENDDO

       ENDDO
       !-----
    ENDIF
  END SUBROUTINE  emis_flam_smold
  !----------------------------------------------------------------------


  !----------------------------------------------------------------------

  SUBROUTINE sources(m1,m2,m3,ia,iz,ja,jz,itime1,time,imonth1,idate1,iyear1,glon,    &
                     chem_nspecies,spc_chem_alloc,src,on,off,transport,nsrc,nvert_src,  &
                     chem1_src_g,chem1_g,bburn,bioge,antro,geoge,diur_cycle,aer_nvert_src,    &
                     aer1_g,aerosol,aer_nspecies,spc_aer_alloc,nmodes,aer_bburn,        &
                     aer_sdust,aer_urban,aer_bioge,aer_marin,dnp,iexev,dn0,cosz,spc_chem_name, &
                     emiss_cycle,volcanoes,aer_v_ash,CO2,isfcl,spc_aer_name,matrix_level,&
		     aer2_g,SO2)

    use aer1_list, ONLY :  akk,sulf,acc,bc1,bcar,occ,ocar,dd1,dust,dd2,boc, numb_alloc  
    
    implicit none
    
    ! original
    INTEGER , INTENT(IN) :: m1
    INTEGER , INTENT(IN) :: m2
    INTEGER , INTENT(IN) :: m3
    INTEGER , INTENT(IN) :: ia
    INTEGER , INTENT(IN) :: iz
    INTEGER , INTENT(IN) :: ja
    INTEGER , INTENT(IN) :: jz
    REAL    , INTENT(IN) :: time
    INTEGER , INTENT(IN) :: imonth1
    INTEGER , INTENT(IN) :: idate1
    INTEGER , INTENT(IN) :: iyear1
    INTEGER , INTENT(IN) :: itime1
    INTEGER , INTENT(IN) :: isfcl,CO2  !DSM
    INTEGER , INTENT(IN) :: SO2  ! for matrix

    ! chem1_list
    INTEGER , INTENT(IN) :: chem_nspecies
    INTEGER , INTENT(IN) :: spc_chem_alloc(6,chem_nspecies)
    INTEGER , INTENT(IN) :: src
    INTEGER , INTENT(IN) :: on
    INTEGER , INTENT(IN) :: off
    INTEGER , INTENT(IN) :: transport
    CHARACTER(LEN=8), INTENT(IN),DIMENSION(chem_nspecies) :: spc_chem_name

    ! aer1_list
    INTEGER , INTENT(IN) :: aer_nspecies
    INTEGER , INTENT(IN) :: spc_aer_alloc(6,nmodes,aer_nspecies)
    INTEGER , INTENT(IN) :: nmodes
    INTEGER , INTENT(IN) :: aer_bburn
    INTEGER , INTENT(IN) :: aer_sdust
    INTEGER , INTENT(IN) :: aer_urban
    INTEGER , INTENT(IN) :: aer_bioge
    INTEGER , INTENT(IN) :: aer_marin
    INTEGER , INTENT(IN) :: aer_v_ash
    CHARACTER(LEN=8) , INTENT(IN) :: spc_aer_name(nmodes,aer_nspecies)
    CHARACTER(LEN=1 ), INTENT(IN) :: matrix_level

    ! mem_chem1
    INTEGER              , INTENT(IN)    :: nsrc
    INTEGER              , INTENT(IN)    :: nvert_src(nsrc)
    TYPE(chem1_src_vars) , INTENT(INOUT) :: chem1_src_g(2,nsrc,chem_nspecies)
    TYPE(chem1_vars)     , INTENT(INOUT) :: chem1_g(chem_nspecies)
    INTEGER              , INTENT(IN)    :: bburn
    INTEGER              , INTENT(IN)    :: bioge
    INTEGER              , INTENT(IN)    :: geoge
    INTEGER              , INTENT(IN)    :: antro
    INTEGER              , INTENT(IN)    :: diur_cycle(nsrc)

    ! mem_aer1
    INTEGER         , INTENT(IN)    :: aer_nvert_src(aer_nspecies)
    TYPE(aer1_vars) , INTENT(INOUT) :: aer1_g(nmodes,aer_nspecies)
    TYPE(aer1_vars) , INTENT(INOUT) :: aer2_g(nmodes)
    INTEGER         , INTENT(IN)    :: aerosol

    ! mem_stilt
    REAL    , POINTER    :: dnp(:,:,:) ! in
    INTEGER , INTENT(IN) :: iexev

    ! mem_basic
    REAL , POINTER :: dn0(:,:,:) ! in

    ! mem_radiate
    REAL , INTENT(IN) :: cosz(m2,m3)

    TYPE(cycle_emission), INTENT(INOUT) :: emiss_cycle(nsrc)

    ! mem_volc_chem1
    INTEGER, INTENT(IN) :: volcanoes

    INTEGER :: k_src,k_tend,it1,it2

    DOUBLE PRECISION :: tlinterp
    !-MFA
    REAL :: timeq2,timeq3,gglon,fuso,alfa(nsrc)
    !-for Lagrange
    REAL :: g,no,src1,src2,srcn1,srcn2,time1,time2,ztime1,ztime2,local_hour,htime1,htime2
    DOUBLE PRECISION :: tlinterp2
    !- ANTRO FROM CETESB-
    INTEGER, PARAMETER :: diur_cetesb_flag=0
    !-MFA
    
    INTEGER :: iweek,idays,j,i,k,ispc,isrc,k2,imode
    REAL :: tign,strtim,timeq,r_q,r_antro,real_time,jd
    REAL, DIMENSION(7) :: week_CYCLE
    !                     dia da semana:    SEG   TER	QUA   QUI   SEX   SAB  DOM  
    !                            iweek=     1      2	  3	4     5     6	 7
    !- dados cetesb/campinas/2005
    DATA (week_CYCLE(iweek),iweek=1,7) /1.1, 1.1, 1.1, 1.1, 1.1, 0.83,0.67/ !total = 7

    REAL rt(nsrc),rt_aer(aer_nspecies)
    REAL, PARAMETER :: bx_bburn  = 18.041288 * 3600., & !- pico em 18 UTC
         cx        = 2.5 * 3600.,       & ! 2.184936 * 3600., &
         rinti     = 0.8941851* 2.1813936e-8    , & ! 1/integral
         ax        = 2000.6038        , &
         bx_antro  = 9. *3600.  ,&  ! local time of peak 1 
         cx_antro  = 16.*3600.  ,&  ! local time of peak 2
         rsum      = 1.3848466E-05  !2.1311416E-05   ! 1/integral

    REAL, POINTER, SAVE,DIMENSION(:,:,:) :: rho_air
    REAL, DIMENSION(m2,m3) :: glon

    REAL local_cosz(m2,m3), local_emiss_bioge_diur_cycle(m2,m3)

    !- nocturnal/background/constant emission for biogenic/urban-industrial-transp processes
    REAL, PARAMETER :: f_nct=0.15                  &! 15% per day
         , f_nct_dvd86400=f_nct/86400. &
         , um86400=1./86400.

    !- parameters for converting emission in mass to emission in number 
    !- only for matrix 
    integer, parameter :: nemis_spcs        = 10 
    real(4), parameter :: pi                = 3.141592654
    real(4), parameter :: pi6               = pi/6.0
    real(4), parameter :: emis_dens_sulf    = 1.770e+00 !g/cm^3]
    real(4), parameter :: emis_dens_bcar_bb = 1.390e+00 ![g/cm^3]
    real(4), parameter :: emis_dens_bcar_ur = 1.700e+00 ![g/cm^3]
    real(4), parameter :: emis_dens_ocar    = 1.000e+00 ![g/cm^3]
    real(4), parameter :: emis_dens_seas    = 2.165e+00 ![g/cm^3]
    real(4), parameter :: emis_dens_dust    = 2.600e+00 ![g/cm^3]
    real, dimension(nemis_spcs) :: emis_dens = (/ emis_dens_sulf, emis_dens_sulf, &
         emis_dens_bcar_ur, emis_dens_bcar_bb,emis_dens_bcar_bb, emis_dens_ocar,  &
         emis_dens_seas   , emis_dens_seas   , emis_dens_dust  , emis_dens_dust /)

                                                 

    real(4), parameter :: dg_akk_sulf = 0.013   ![um]
    real(4), parameter :: dg_acc_sulf = 0.068   ![um]
    real(4), parameter :: dg_bc1_bcar = 0.030   ![um]
    real(4), parameter :: dg_boc_bcar = 0.021   ![um]
    real(4), parameter :: dg_boc_ocar = 0.021   ![um]
    real(4), parameter :: dg_occ_ocar = 0.030   ![um]
    real(4), parameter :: dg_ssa_seas = 0.370   ![um]
    real(4), parameter :: dg_ssc_seas = 3.930   ![um]
    real(4), parameter :: dg_dd1_dust = 0.580   ![um]
    real(4), parameter :: dg_dd2_dust = 5.400   ![um]


    real, dimension(nemis_spcs) :: dgn0_emis = (/ dg_akk_sulf, dg_acc_sulf, &
         dg_bc1_bcar, dg_boc_bcar, dg_boc_ocar, dg_occ_ocar, dg_ssa_seas, &
	 dg_ssc_seas, dg_dd1_dust, dg_dd2_dust /)

    real(4), parameter :: sg_akk_sulf = 1.600d+00
    real(4), parameter :: sg_acc_sulf = 1.800d+00
    real(4), parameter :: sg_bc1_bcar = 1.800d+00
    real(4), parameter :: sg_boc_bcar = 1.550d+00
    real(4), parameter :: sg_boc_ocar = 1.555d+00
    real(4), parameter :: sg_occ_ocar = 1.800d+00
    real(4), parameter :: sg_ssa_seas = 1.800d+00
    real(4), parameter :: sg_ssc_seas = 2.000d+00
    real(4), parameter :: sg_dd1_dust = 1.800d+00
    real(4), parameter :: sg_dd2_dust = 1.800d+00 
 
    real, dimension(nemis_spcs) :: sig0_emis = (/ sg_akk_sulf, sg_acc_sulf, &
         sg_bc1_bcar, sg_boc_bcar, sg_boc_ocar, sg_occ_ocar, sg_ssa_seas, &
	 sg_ssc_seas, sg_dd1_dust, sg_dd2_dust /)

    real, dimension(nemis_spcs) ::recip_part_mass,dp0_emis
    REAL, DIMENSION(m1,m2,m3) :: dummy_src
    REAL :: factor
    !-end of parameters for matrix.

    dummy_src=0.0
   
    !- if using mass conservation fix : air dens changes with the time evolution
    !   IF( iexev == 2 )  rho_air => stilt_g(ng)%dnp(:,:,:) 
    IF( iexev == 2 )  rho_air => dnp(:,:,:) 

    !- if not, air dens = air dens of basic state and  need to define for each when 
    !- have nested grids
    !   IF( iexev == 1  ) rho_air => basic_g(ng)%dn0(:,:,:) 
    IF( iexev == 1  ) rho_air => dn0(:,:,:) 
    !
    !number of days of simulation
    idays = INT(( float(itime1)/100. + time/3600.)/24.+.00001)  

    
    
    IF(diur_cycle(bburn) == ON) THEN
       !-------------biomass burning diurnal cycle --------------------
       tign  = REAL(idays)*24.*3600.

       ! Modulacao da queimada media durante o ciclo diurno(unidade: 1/s)
       ! com a int( r_q dt) (0 - 24h)= 1.
       timeq= ( time + float(itime1)*0.01*3600. - tign )

       r_q  = rinti*( ax * EXP( -(timeq-bx_bburn)**2/(2.*cx**2) ) + 100. -  &
            5.6712963e-4*( timeq ))

       emiss_cycle(bburn)%emission_rate(:,:)=r_q
       rt(bburn)= r_q
       alfa(bburn) = 0.
    ELSE
       emiss_cycle(bburn)%emission_rate(:,:)= um86400 ! = 1./86400.
       rt(bburn)= um86400 ! = 1./86400.
       alfa(bburn) = 1.
    ENDIF


    IF(diur_cycle(antro) == ON) THEN

       !------------- anthropogenic diurnal cycle (industrial,urban, ...)
       real_time = time + float(itime1)*0.01*3600. !UTC

       ! weekly cycle
       ! week day
       !v1
       iweek= INT(((float(julday(imonth1,idate1,iyear1))/7. - &
            INT(julday(imonth1,idate1,iyear1)/7))*7.)) + 3
       IF(iweek.GT.7) iweek = iweek-7
       !v2
       !call Greg2Jul(0, idate1, imonth1, iyear1, jd)
       !jd=jd+float(nint(time/86400.))
       !iweek = int(AMOD(jd, 7.)) 
       !if(iweek < 1 .or. iweek > 7) stop 315

       !- diurnal cycle
       IF (diur_cetesb_flag == 0) THEN
         DO j=ja,jz
     	  DO i=ia,iz

             gglon = glon(i,j)
             fuso = INT(gglon/15)
             !-to better keep continuity
             idays = INT(( real_time +fuso*3600. )/86400.+.00001)  
             tign  = REAL(idays)*86400.

             timeq2= (real_time  -tign + fuso*3600.) - bx_antro
             timeq3= (real_time  -tign + fuso*3600.) - cx_antro
             !r_antro  =(exp(-((timeq2)**2)/((3.*3600) **2))+(exp(-((timeq3)**2)/((3.*3600) **2))+0.1))*rsum
             r_antro=(EXP(-(timeq2**2)/(18400.**2))+EXP(-(timeq3**2)/(18500.**2))+0.1)*rsum

             !- weekly + diurnal cycle
             r_antro = r_antro * week_CYCLE(iweek)

             emiss_cycle(antro)%emission_rate(i,j)=r_antro 
             alfa(antro) = 0.
    	  ENDDO
         ENDDO
       
       ENDIF
       IF (diur_cetesb_flag == 1) THEN
          
           DO j=ja,jz
              DO i=ia,iz

        	 gglon = glon(i,j)
        	 fuso  = INT(gglon/15)
        	 
		 !-to better keep continuity
        	 
		 idays = INT(( real_time )/86400.+.00001)  
        	 tign  = REAL(idays)*86400.
        	 
		 
		 local_hour = ((real_time - tign)/3600)+fuso
        	 
        	 
		 IF (local_hour.LT.0) THEN
        	 local_hour = (24) + local_hour
        	 ENDIF
        	 
        	
		 CALL interpolation_antro(local_hour,src1,src2,srcn1,srcn2,time1,time2)
        	 
        	 
		 htime1 = time1 - fuso
        	 IF (htime1.GE.24) THEN
        	 htime1 = htime1 - 24
        	 ENDIF
        	 
        	 htime2 = time2 - fuso
        	 IF (htime2.GE.24) THEN
        	 htime2 = htime2 - 24
        	 ENDIF
        	 
		 
		 ztime1 = ((htime1*3600) + tign)
        	 ztime2 = ((htime2*3600) + tign)
        	 
		 
		 tlinterp=DBLE(src1 + (((time - ztime1)*(src2 - src1))/(ztime2 - ztime1)))
        	 tlinterp2=DBLE(srcn1 + (((time - ztime1)*(srcn2 - srcn1))/(ztime2 - ztime1)))
        	 
        	 !PRINT*,'TESTE',time,local_hour,time1,time2,ztime1,ztime2,tign,src2,src1,tlinterp,fuso
        	 
        	 emiss_cycle(antro)%emission_rate    (:,:)=tlinterp/3600
        	 emiss_cycle(antro)%emission_rate_NOX(:,:)=tlinterp2/3600
        	 
        	 alfa(antro)= 0.    
        	 r_antro=tlinterp/3600.
              ENDDO
           ENDDO
       ENDIF
        IF (diur_cetesb_flag.NE.1.AND.diur_cetesb_flag.NE.0) THEN
           CALL fatal_error('No definition for diur_cetesb_flag!')
       ENDIF
	
    ELSE
       !----------- sources linearly time interpolated
       IF(srctime2<=srctime1) THEN
           CALL fatal_error('srctime2<=srctime1! linterp')
       ENDIF
       tlinterp=DBLE(time-srctime1)/DBLE(srctime2-srctime1)
       emiss_cycle(antro)%emission_rate(:,:)=tlinterp
       alfa(antro)= 1.
    ENDIF

    IF(diur_cycle(bioge) == ON) THEN

       !---------- sources with diurnal cycle and spatial dependence
       ! - using zenital angle from radiate routine/including constant background emission (f_cnt%)
!       emiss_cycle(bioge,ng)%emission_rate(:,:)= f_nct_dvd86400 + emiss_cycle(bioge,ng)%dcnorma_inv(:,:)&    
!            *MAX(0.,radiate_g(ng)%cosz(:,:)) * (1.-f_nct) 
       emiss_cycle(bioge)%emission_rate(:,:)= f_nct_dvd86400 + emiss_cycle(bioge)%dcnorma_inv(:,:)&    
            *MAX(0.,cosz(:,:)) * (1.-f_nct) 
       alfa(bioge) = 0.
    ELSE
       !----------- sources linearly time interpolated    
       IF(srctime2<=srctime1) THEN
            CALL fatal_error('srctime2==srctime1! linterp')
       ENDIF
       tlinterp=DBLE(time-srctime1)/DBLE(srctime2-srctime1)
       emiss_cycle(bioge)%emission_rate(:,:) = tlinterp
       alfa(bioge) = 1.
    ENDIF

    IF(diur_cycle(geoge) == ON  ) THEN
       emiss_cycle(geoge)%emission_rate(:,:)=um86400 ! 1./86400 sec
       alfa(geoge) = 0.
    ELSE
       !----------- sources linearly time interpolated    
       IF(srctime2<=srctime1) THEN
            CALL fatal_error('srctime2==srctime1! linterp')
       ENDIF
       tlinterp=DBLE(time-srctime1)/DBLE(srctime2-srctime1)
       emiss_cycle(geoge)%emission_rate(:,:) = tlinterp
       alfa(geoge) = 1.
    ENDIF

    !-srf temporary array to store emission cycle for biogenic species.
    !-    (CO2 from JULES is already instantaneous while others biogenic emissions 
    !-    from MEGAN is per day)
    IF(ISFCL == 5) local_emiss_bioge_diur_cycle(:,:)=emiss_cycle(bioge)%emission_rate(:,:)
    !

    ! print*,'emiss_cycle(antro,ng)%emission_rate=',emiss_cycle(antro,ng)%emission_rate(ia,ja)&
    ! ,srctime1,srctime2,tlinterp

    !-------------------------- perform emissions
    !- chemistry section 
    DO ispc=1,chem_nspecies

       IF(spc_chem_alloc(src,ispc) == off .OR. spc_chem_alloc(transport,ispc) == off)  CYCLE

       DO isrc=1,nsrc

          !- DSM/SRF - special treatment for biogenic/jules (CO2) emissions
          IF(ISFCL == 5 .and. isrc == bioge) then
             IF( ispc == CO2 ) THEN
	        emiss_cycle(bioge)%emission_rate(:,:)=1.
             ELSE
	        emiss_cycle(bioge)%emission_rate(:,:)=local_emiss_bioge_diur_cycle(:,:)
             END IF
	  END IF
          
          !- only call geoge emissions if volcanoes is ON.
	  IF(isrc==geoge .AND. volcanoes == off) CYCLE

          !- memory position of source array for each time level
          it1 = 1 ! always 1
          it2 = actual_time_index(2,isrc)!might be 1 or 2, will be 1 if linterp = OFF
                                         !=>chem1_src_g(it1,:,:,:)% = chem1_src_g(it2,:,:,:)%

          k_src = 1 
          !- for anthropogenic and biogenic emissions
	  k2 = 2                                          ! control for 2-dim src emission field
        
	  !- for biomass burning and volcanic emissions          
	  IF(isrc == bburn .OR. isrc == geoge)  k2 = m1-1  ! control for 3-dim src emission field

          DO k=2,k2

             IF(isrc == bburn .OR. isrc == geoge) k_src=k  ! control for 3-dim source emission field


	     IF (diur_cetesb_flag == 1 .AND. isrc == antro) THEN

	    	CALL source_to_tend_cycle_cetesb ( m1,m2,m3,ia,iz,ja,jz   &
            	     ,chem1_src_g(it1,isrc,ispc)%sc_src &! source data at time level 1
            	     ,chem1_src_g(it2,isrc,ispc)%sc_src &! source data at time level 2
            	     ,chem1_g	 (	   ispc)%sc_t	&! tendency array
            	     ,nvert_src(isrc)	        	&! vertical size of source array
            	     ,emiss_cycle(isrc)%emission_rate   &! diurnal cycle of emission
            	     ,alfa(isrc)			&! alfa cte 
            	     ,k_src				&! vertical level of source array (where data is stored)
            	     ,k 				&! vertical level of tendency
            	     ,rho_air				&! air density (to convert to mixing ratio tendency)
	   	     ,emiss_cycle(isrc)%emission_rate_NOX &! peso hora atual NOX     
	   	     ,ispc,spc_chem_name,chem_nspecies  )				      

             ELSE

            	CALL source_to_tend_cycle ( m1,m2,m3,ia,iz,ja,jz   &
            	     ,chem1_src_g(it1,isrc,ispc)%sc_src	   &! source data at time level 1
            	     ,chem1_src_g(it2,isrc,ispc)%sc_src	   &! source data at time level 2
            	     ,chem1_g	 (	   ispc)%sc_t      &! tendency array
            	     ,nvert_src(isrc)			   &! vertical size of source array
            	     ,emiss_cycle(isrc)%emission_rate	   &! diurnal cycle of emission
            	     ,alfa(isrc)			   &! alfa cte 
            	     ,k_src				   &! vertical level of source array (where data is stored)
            	     ,k 				   &! vertical level of tendency
            	     ,rho_air)  			    ! air density (to convert to mixing ratio tendency)
             ENDIF
 
          ENDDO
       ENDDO
    ENDDO

    !- aerosol section ----------------------------------------
    IF(AEROSOL == 1)  THEN
       !- still need implementation of the emission cycle with space dependence
       rt_aer(aer_bburn) = rt(bburn)
       rt_aer(aer_sdust) = 1.0 ! on-line emission
       rt_aer(aer_urban) = 1.157407e-5 !rt(antro)  ! < must be fixed later with actual diurnal cycle
       rt_aer(aer_bioge) = 1.157407e-5 
       rt_aer(aer_marin) = 1.0 ! on-line emission 
       rt_aer(aer_v_ash) = 1.157407e-5 ! should be 1./duration

       DO ispc=1,aer_nspecies

          DO imode=1,nmodes

             IF(spc_aer_alloc(src      ,imode,ispc) == off .OR. &
                spc_aer_alloc(transport,imode,ispc) == off) CYCLE


             k_src = 1;		  k2 = 2     ! control for 2-dim aerosol (bioge, antro, marin, sdust)
             IF(ispc == aer_bburn .or. ispc == aer_v_ash)  k2 = m1-1  ! control for 3-dim bburn aerosol

             DO k=2,k2

                IF(ispc == aer_bburn .or. ispc == aer_v_ash) k_src=k     ! control  for 3-dim bburn aerosol

                CALL source_to_tend (m1,m2,m3,ia,iz,ja,jz &
                     ,aer1_g(imode,ispc)%sc_src           &! source data
                     ,aer1_g(imode,ispc)%sc_t		  &! tendency array
                     ,aer_nvert_src(ispc)		  &! vertical size of source array
                     ,rt_aer(ispc)			  &! diurnal cycle of emission
                     ,k_src 				  &! vertical level of source array (where data is stored)
                     ,k					  &! vertical level of tendency
                     ,rho_air)				   ! air density (to convert to mixing ratio tendency)

             ENDDO
          ENDDO
       ENDDO
    ELSEIF(AEROSOL == 2)  THEN
       
       if(matrix_level .ne. "1") stop "Emissions are only configured for matrix_level=1"
       
       !- section for emission of mass --------------------
       DO ispc=1,aer_nspecies
          DO imode=1,nmodes
             IF(spc_aer_alloc(src      ,imode,ispc) == off .OR. &
                spc_aer_alloc(transport,imode,ispc) == off) CYCLE

                !- case anthropogenic emssions
	         IF(   spc_aer_name(imode,ispc)=="occ_ocar" .or. &
	               spc_aer_name(imode,ispc)=="bc1_bcar" )	 THEN
		    k2  	 = 2	 ! control for 2-dim aerosol (bioge, antro, marin, sdust)
                    rt_aer(ispc) = r_antro
	            k_src = 1
                    !print*,"anthro",maxval(aer1_g(imode,ispc)%sc_src),spc_aer_name(imode,ispc), rt_aer(ispc)
		
                    !DO k=2,k2
                        CALL source_to_tend (m1,m2,m3,ia,iz,ja,jz &
                        ,aer1_g(imode,ispc)%sc_src         	  &! source data
                        ,aer1_g(imode,ispc)%sc_t		  &! tendency array
                        ,aer_nvert_src(ispc)		          &! vertical size of source array
                        ,rt_aer(ispc)				  &! diurnal cycle of emission
                        ,k_src 				          &! vertical level of source array (where data is stored)
                        ,k2					  &! vertical level of tendency
                        ,rho_air)				   ! air density (to convert to mixing ratio tendency)
		   ! ENDDO
               
	        !- case sea salt emssions
	        ELSEIF( spc_aer_name(imode,ispc)=="ssa_seas" .or. &
	                spc_aer_name(imode,ispc)=="ssc_seas" )    THEN
		    k2  	 = 2	 ! control for 2-dim aerosol (bioge, antro, marin, sdust)
                    rt_aer(ispc) = 1.0
	            k_src = 1
                    !print*,"seas",maxval(aer1_g(imode,ispc)%sc_src),spc_aer_name(imode,ispc)
		
                    !DO k=2,k2
                        CALL source_to_tend (m1,m2,m3,ia,iz,ja,jz &
                        ,aer1_g(imode,ispc)%sc_src         	  &! source data
                        ,aer1_g(imode,ispc)%sc_t		  &! tendency array
                        ,aer_nvert_src(ispc)		          &! vertical size of source array
                        ,rt_aer(ispc)				  &! diurnal cycle of emission
                        ,k_src 				          &! vertical level of source array (where data is stored)
                        ,k2					  &! vertical level of tendency
                        ,rho_air)				   ! air density (to convert to mixing ratio tendency)
		    !ENDDO
		
		!- case biomass burning emissions !the same must be for volcanic ASH
		 ELSEIF(spc_aer_name(imode,ispc)=="boc_bcar" .or. &
	                spc_aer_name(imode,ispc)=="boc_ocar" )    THEN
                    !print*,"bb",maxval(aer1_g(imode,ispc)%sc_src),spc_aer_name(imode,ispc), rt(bburn)
		    rt_aer(ispc) = rt(bburn)
                    k2 = m1-1  !
                    DO k=2,k2
                        k_src=k
                        
			CALL source_to_tend (m1,m2,m3,ia,iz,ja,jz &
                        ,aer1_g(imode,ispc)%sc_src         	  &! source data
                        ,aer1_g(imode,ispc)%sc_t		  &! tendency array
                        ,aer_nvert_src(ispc)		          &! vertical size of source array
                        ,rt_aer(ispc)				  &! diurnal cycle of emission
                        ,k_src 				          &! vertical level of source array (where data is stored)
                        ,k					  &! vertical level of tendency
                        ,rho_air)				   ! air density (to convert to mixing ratio tendency)
		    ENDDO
                 ENDIF
          ENDDO
       ENDDO
       !- section for emission of mass 
       !- special section for on-line emissions (not from prep-chem-src)
       !- this section will provide emissions for aerosol species akk_sulf and acc_sulf
       !- for antro and bio-burn processes. These emissions are based on the SO2 gas emission
       !===> be sure emissions for the SO2 gas is provided (arrays chem1_src_g(*,isrc,SO2)%sc_src)
       IF(spc_chem_alloc(src,SO2) == ON)  then
       
         DO ispc=sulf,sulf
            DO imode=1,nmodes
                 IF(spc_aer_alloc(transport,imode,ispc) == off) CYCLE
                
		
	         IF(   spc_aer_name(imode,ispc)=="akk_sulf" .or. &
		       spc_aer_name(imode,ispc)=="acc_sulf"      ) THEN
		     !print*,"src emissions for aer=",imode, ispc,spc_aer_name(imode,ispc)
                 
		     IF(   spc_aer_name(imode,ispc)=="akk_sulf" ) factor = 0.01 * 0.025* 96.07/62.66  
		 
		     IF(   spc_aer_name(imode,ispc)=="acc_sulf" ) factor = 0.99 * 0.025* 96.07/62.66  
 		 
		 
		      ! 1) case anthropogenic emissions
		      isrc  = antro
                      dummy_src(1:m1,1:m2,1:m3) = 0.5*( chem1_src_g(it1,isrc,SO2)%sc_src(1:m1,1:m2,1:m3) + &
            	                                        chem1_src_g(it2,isrc,SO2)%sc_src(1:m1,1:m2,1:m3)  )
		
		      k2  	 = 2	 ! control for 2-dim aerosol (bioge, antro, marin, sdust)
                      rt_aer(ispc) = r_antro * factor
	              k_src = 1
                      CALL source_to_tend (m1,m2,m3,ia,iz,ja,jz &
                      ,dummy_src		                &! source data
                      ,aer1_g(imode,ispc)%sc_t  		&! tendency array
                      ,aer_nvert_src(ispc)			&! vertical size of source array
                      ,rt_aer(ispc)				&! diurnal cycle of emission
                      ,k_src					&! vertical level of source array (where data is stored)
                      ,k2					&! vertical level of tendency
                      ,rho_air) 				 ! air density (to convert to mixing ratio tendency)

		      ! 2) case biomass burning emissions
		      isrc  = bburn
                      dummy_src(1:m1,1:m2,1:m3) = 0.5*( chem1_src_g(it1,isrc,SO2)%sc_src(1:m1,1:m2,1:m3) + &
            	                                        chem1_src_g(it2,isrc,SO2)%sc_src(1:m1,1:m2,1:m3)  )
		      rt_aer(ispc) = rt(bburn)* factor
                      k2 = m1-1  !
                      DO k=2,k2
                        k_src=k
                        
			CALL source_to_tend (m1,m2,m3,ia,iz,ja,jz &
                        ,dummy_src       	                  &! source data
                        ,aer1_g(imode,ispc)%sc_t		  &! tendency array
                        ,aer_nvert_src(ispc)		          &! vertical size of source array
                        ,rt_aer(ispc)				  &! diurnal cycle of emission
                        ,k_src 				          &! vertical level of source array (where data is stored)
                        ,k					  &! vertical level of tendency
                        ,rho_air)				   ! air density (to convert to mixing ratio tendency)
		      ENDDO

                 ENDIF

            ENDDO
         ENDDO
       ELSE
        stop "Emissions for mass of aer akk/acc are required, but SO2 emissions are not set"
       
       ENDIF
       !
       !print*,"done emissions for mass of akk acc"; call flush(6)
       !
       !--- section for emission of number concentration
       !- A) getting parameter recip_part_mass
       DO i =1,nemis_spcs
       
         dp0_emis(i) = 1.0e-06 * dgn0_emis(i) * exp( 1.5e+00 * ( log(sig0_emis(i)) )**2 )  ! convert from [um] to [m]
         
	 !recip_part_mass(i) = 1.0e+00 / ( 1.0e+12 * emis_dens(i) * pi6 * dp0_emis(i)**3 )
	 !- unit 1/kg
	 recip_part_mass(i) = 1.0e+00 / ( 1.0e+3 * emis_dens(i) * pi6 * dp0_emis(i)**3 )
       
       ENDDO
       !
       !- B) getting number emissions in terms of mass emissions and performing the emissions
       !
       DO i =1,nemis_spcs
        if(i==1 ) then ;imode=akk ;   ispc =sulf; endif  !   
        if(i==2 ) then ;imode=acc ;   ispc =sulf; endif    
        if(i==3 ) then ;imode=bc1 ;   ispc =bcar; endif    
        if(i==4 ) then ;imode=boc ;   ispc =bcar; endif  ! boc bcar + boc ocar
        if(i==5 ) then ; cycle                  ; endif  ! DUMMY  
        if(i==6 ) then ;imode=occ ;   ispc =ocar; endif   
        if(i==7 ) then ; cycle                  ; endif  ! sea salt has not number emission  
        if(i==8 ) then ; cycle                  ; endif  ! sea salt has not number emission   
        if(i==9 ) then ;imode=dd1 ;   ispc =dust; endif    
        if(i==10) then ;imode=dd2 ;   ispc =dust; endif    
        !print*,"nemis=",i,imode,ispc;call flush(6) 
	IF(numb_alloc(transport,imode) == off) cycle

          !-special treatment for akk/acc - sulf	 
	  IF( (imode == akk .or. imode==acc) .and. ispc==sulf) THEN
	           
	     if(  imode == akk ) factor = 0.01 * 0.025* 96.07/62.66 * 1.e-9*recip_part_mass(i)   !-unit  #/m3/s	
	     if(  imode == acc ) factor = 0.99 * 0.025* 96.07/62.66 * 1.e-9*recip_part_mass(i)   !-unit  #/m3/s
		 
	     ! 1) case anthropogenic emissions
	     isrc  = antro
             dummy_src(1:m1,1:m2,1:m3) = 0.5*( chem1_src_g(it1,isrc,SO2)%sc_src(1:m1,1:m2,1:m3) + &
             				       chem1_src_g(it2,isrc,SO2)%sc_src(1:m1,1:m2,1:m3)  )
	
	     k2 	= 2	! control for 2-dim aerosol (bioge, antro, marin, sdust)
             rt_aer(ispc) = r_antro * factor 
	     k_src = 1
             CALL source_to_tend (m1,m2,m3,ia,iz,ja,jz &
             ,dummy_src 			       &! source data
             ,aer2_g(imode)%sc_t		       &! tendency array
             ,aer_nvert_src(ispc)		       &! vertical size of source array
             ,rt_aer(ispc)			       &! diurnal cycle of emission
             ,k_src				       &! vertical level of source array (where data is stored)
             ,k2				       &! vertical level of tendency
             ,rho_air)  				! air density (to convert to mixing ratio tendency)
 	      

	     ! 2) case biomass burning emissions
	     isrc  = bburn
             dummy_src(1:m1,1:m2,1:m3) = 0.5*( chem1_src_g(it1,isrc,SO2)%sc_src(1:m1,1:m2,1:m3) + &
                                               chem1_src_g(it2,isrc,SO2)%sc_src(1:m1,1:m2,1:m3)  )
	     rt_aer(ispc) = rt(bburn)* factor
             k2 = m1-1  !
             DO k=2,k2
                k_src=k
                
		CALL source_to_tend (m1,m2,m3,ia,iz,ja,jz &
                ,dummy_src       	                  &! source data
                ,aer2_g(imode)%sc_t		          &! tendency array
                ,aer_nvert_src(ispc)		          &! vertical size of source array
                ,rt_aer(ispc)				  &! diurnal cycle of emission
                ,k_src 				          &! vertical level of source array (where data is stored)
                ,k					  &! vertical level of tendency
                ,rho_air)				   ! air density (to convert to mixing ratio tendency)
	      ENDDO
	  !-special treatment for boc - bcar 
	  ELSEIF( imode == boc .and. ispc==bcar) THEN
	  
	     !-unit  #/m3/s
	     dummy_src(:,:,:) =1.e-9*(recip_part_mass(4) * aer1_g(boc,bcar)%sc_src(:,:,:) +&
	                              recip_part_mass(5) * aer1_g(boc,ocar)%sc_src(:,:,:) )
	     ! only biomass burning emissions
	     isrc  = bburn
	     rt_aer(ispc) = rt(bburn)
             k2 = m1-1  !
             DO k=2,k2
                k_src=k
		CALL source_to_tend (m1,m2,m3,ia,iz,ja,jz &
                ,dummy_src       	                  &! source data
                ,aer2_g(imode)%sc_t		          &! tendency array
                ,aer_nvert_src(ispc)		          &! vertical size of source array
                ,rt_aer(ispc)				  &! diurnal cycle of emission
                ,k_src 				          &! vertical level of source array (where data is stored)
                ,k					  &! vertical level of tendency
                ,rho_air)				   ! air density (to convert to mixing ratio tendency)
	     ENDDO
	  
	  
	  ELSE
	      IF(spc_aer_alloc(src      ,imode,ispc) == off) then
	        print*," not aloc for",imode,ispc
	        stop 7777
	      endif
	  
	     !-unit #/m3/s
	     dummy_src(:,:,:) =1.e-9*recip_part_mass(i) * aer1_g(imode,ispc)%sc_src(:,:,:)
             !       
	     ! only anthropogenic emissions
	     isrc  = antro	
	     k2	     = 2     ! control for 2-dim aerosol (antro, marin, sdust)
             rt_aer(ispc) = r_antro 
	     k_src = 1
             CALL source_to_tend (m1,m2,m3,ia,iz,ja,jz &
              ,dummy_src			    &! source data
              ,aer2_g(imode)%sc_t		    &! tendency array
              ,aer_nvert_src(ispc)  		    &! vertical size of source array
              ,rt_aer(ispc) 			    &! diurnal cycle of emission
              ,k_src				    &! vertical level of source array (where data is stored)
              ,k2				    &! vertical level of tendency
              ,rho_air)				     ! air density (to convert to mixing ratio tendency)
	  
          ENDIF	
       ENDDO
    ENDIF ! end of aerosol section (aerosol model) 
  !
  END SUBROUTINE sources
  !------------------------------------------------------------------

  SUBROUTINE source_to_tend(m1,m2,m3,ia,iz,ja,jz,sc_src,sc_t,nvert,rt,k_src,k_tend,rho_air)
   
    ! original
    INTEGER , INTENT(IN)    :: m1
    INTEGER,  INTENT(IN)    :: m2
    INTEGER,  INTENT(IN)    :: m3
    INTEGER,  INTENT(IN)    :: ia
    INTEGER,  INTENT(IN)    :: iz
    INTEGER,  INTENT(IN)    :: ja
    INTEGER,  INTENT(IN)    :: jz
    REAL,     INTENT(IN)    :: sc_src(nvert,m2,m3)
    REAL,     INTENT(INOUT) :: sc_t(m1,m2,m3)
    INTEGER,  INTENT(IN)    :: nvert
    REAL,     INTENT(IN)    :: rt
    INTEGER,  INTENT(IN)    :: k_src
    INTEGER,  INTENT(IN)    :: k_tend
    REAL,     INTENT(IN)    :: rho_air(m1,m2,m3)

    !- also this routine assumes that the source term is expressed in terms of  density (kg/m3)
    !- and, then, the conversion to mixing ratio (division by air density) is needed.
    !- k_src express the level where the src data is stored in the sc_src array
    sc_t(k_tend,ia:iz,ja:jz) = sc_t(k_tend,ia:iz,ja:jz) +  sc_src(k_src,ia:iz,ja:jz)*rt /  &
                                !- the level of air density
                                !- must be the same of tendency
                               rho_air(k_tend,ia:iz,ja:jz) 

  END SUBROUTINE source_to_tend


  !------------------------------------------------------------------

  SUBROUTINE source_to_tend_cycle(m1,m2,m3,ia,iz,ja,jz,sc_src1,sc_src2,sc_t,nvert, &
                                  rt,alfa,k_src,k_tend,rho_air)

    ! original
    INTEGER          , INTENT(IN)    :: m1
    INTEGER          , INTENT(IN)    :: m2
    INTEGER          , INTENT(IN)    :: m3
    INTEGER          , INTENT(IN)    :: ia
    INTEGER          , INTENT(IN)    :: iz
    INTEGER          , INTENT(IN)    :: ja
    INTEGER          , INTENT(IN)    :: jz
    REAL             , INTENT(IN)    :: sc_src1(nvert,m2,m3)
    REAL             , INTENT(IN)    :: sc_src2(nvert,m2,m3)
    REAL             , INTENT(INOUT) :: sc_t(m1,m2,m3)
    INTEGER          , INTENT(IN)    :: nvert
    DOUBLE PRECISION , INTENT(IN)    :: rt(m2,m3)
    REAL             , INTENT(IN)    :: alfa
    INTEGER          , INTENT(IN)    :: k_src
    INTEGER          , INTENT(IN)    :: k_tend
    REAL             , INTENT(IN)    :: rho_air(m1,m2,m3)

    INTEGER i,j

    !- also this routine assumes that the source term is expressed in terms of  density (kg/m3)
    !- and, then, the conversion to mixing ratio (division by air density) is needed.
    !- k_src express the level where the src data is stored in the sc_src array
    ! if(nint(alfa) == 0) then

    !   do j=ja,jz ; do i=ia,iz
    !  
    !
    !   sc_t(k_tend,i,j) = sc_t(k_tend,i,j) +  (                        &
    !                      sc_src1(k_src,i,j)* ( 1.0D0- rt(i,j) )*alfa  +  &
    !                      sc_src2(k_src,i,j)* rt(i,j)               )/ &
    !                      rho_air(k_tend,i,j) !- the level of air density
    !	 	                         !- must be the same of tendency
    !                     
    !   enddo;enddo

    !original
    IF(alfa > 0.) THEN 


       sc_t(k_tend,ia:iz,ja:jz) = sc_t   (k_tend ,ia:iz,ja:jz) + (  &
            sc_src1(k_src  ,ia:iz,ja:jz) * ( 1.0D0- rt(ia:iz,ja:jz) )*alfa +     &
            sc_src2(k_src  ,ia:iz,ja:jz) *          rt(ia:iz,ja:jz)        )  /  &
            rho_air(k_tend ,ia:iz,ja:jz) !- the level of air density
                                         !- must be the same of tendency
    ELSE

       sc_t(k_tend,ia:iz,ja:jz) = sc_t   (k_tend ,ia:iz,ja:jz) + (  &
!!!    sc_src1(k_src  ,ia:iz,ja:jz) * ( 1.0D0- rt(ia:iz,ja:jz) )*alfa  +  &
            sc_src2(k_src  ,ia:iz,ja:jz) *          REAL(rt(ia:iz,ja:jz))         )  /  &
            rho_air(k_tend ,ia:iz,ja:jz) !- the level of air density
                                         !- must be the same of tendency

    ENDIF

  END SUBROUTINE source_to_tend_cycle
  !------------------------------------------------------------------
  SUBROUTINE source_to_tend_cycle_cetesb(m1,m2,m3,ia,iz,ja,jz,sc_src1,sc_src2,sc_t,nvert, &
                                  rt,alfa,k_src,k_tend,rho_air,rt2,ispc,spc_chem_name,chem_nspecies)

    !USE chem1_list, ONLY: spc_name
    
    ! original
    INTEGER          , INTENT(IN)    :: chem_nspecies
    INTEGER          , INTENT(IN)    :: m1
    INTEGER          , INTENT(IN)    :: m2
    INTEGER          , INTENT(IN)    :: m3
    INTEGER          , INTENT(IN)    :: ia
    INTEGER          , INTENT(IN)    :: iz
    INTEGER          , INTENT(IN)    :: ja
    INTEGER          , INTENT(IN)    :: jz
    REAL             , INTENT(IN)    :: sc_src1(nvert,m2,m3)
    REAL             , INTENT(IN)    :: sc_src2(nvert,m2,m3)
    REAL             , INTENT(INOUT) :: sc_t(m1,m2,m3)
    INTEGER          , INTENT(IN)    :: nvert
    DOUBLE PRECISION , INTENT(IN)    :: rt(m2,m3)
    REAL             , INTENT(IN)    :: alfa
    INTEGER          , INTENT(IN)    :: k_src
    INTEGER          , INTENT(IN)    :: k_tend
    REAL             , INTENT(IN)    :: rho_air(m1,m2,m3)
    DOUBLE PRECISION , INTENT(IN)    :: rt2(m2,m3) 
    INTEGER          , INTENT(IN)    :: ispc
    CHARACTER(LEN=8), INTENT(IN),DIMENSION(chem_nspecies) :: spc_chem_name
    INTEGER i,j

         IF (spc_chem_name(ispc)=='NO'.OR.spc_chem_name(ispc)=='NO2') THEN
      sc_t(k_tend,ia:iz,ja:jz) = sc_t   (k_tend ,ia:iz,ja:jz) + (  &
!!!                           sc_src1(k_src  ,ia:iz,ja:jz) * ( 1.0D0- rt(ia:iz,ja:jz) )*alfa  +  &
            sc_src2(k_src  ,ia:iz,ja:jz) *          REAL(rt2(ia:iz,ja:jz))         )  /  &
            rho_air(k_tend ,ia:iz,ja:jz) !- the level of air density
       !- must be the same of tendency
     
	 ELSE

      sc_t(k_tend,ia:iz,ja:jz) = sc_t   (k_tend ,ia:iz,ja:jz) + (  &
!!!                           sc_src1(k_src  ,ia:iz,ja:jz) * ( 1.0D0- rt(ia:iz,ja:jz) )*alfa  +  &
            sc_src2(k_src  ,ia:iz,ja:jz) *          REAL(rt(ia:iz,ja:jz))         )  /  &
            rho_air(k_tend ,ia:iz,ja:jz) !- the level of air density
       !- must be the same of tendency
     
         ENDIF

  END SUBROUTINE source_to_tend_cycle_cetesb
  !------------------------------------------------------------------
  ! NOT USED
  !
  !  subroutine invert_air_dens(m1,rho_air)
  !   implicit none
  !   integer m1
  !   real rho_air(m1)
  !   rho_air(:)=1./rho_air(:)
  ! 
  !  end subroutine invert_air_dens
  !------------------------------------------------------------------


  !------------------------------------------------------------------
  SUBROUTINE get_diurnal_cycle_normalized(m2,m3,ia,iz,ja,jz,dtlt,glat, &
                                          glon,imonth1,idate1,iyear1,itime1, &
                                          antro,bioge,pi180,nsrc,emiss_cycle)
      
    ! original
    INTEGER , INTENT(IN) :: m2
    INTEGER , INTENT(IN) :: m3
    INTEGER , INTENT(IN) :: ia
    INTEGER , INTENT(IN) :: iz
    INTEGER , INTENT(IN) :: ja
    INTEGER , INTENT(IN) :: jz
    REAL    , INTENT(IN) :: dtlt
    REAL    , INTENT(IN) :: glat(m2,m3)
    REAL    , INTENT(IN) :: glon(m2,m3)
    
    ! mem_grid
    INTEGER , INTENT(IN) :: imonth1
    INTEGER , INTENT(IN) :: idate1
    INTEGER , INTENT(IN) :: iyear1
    INTEGER , INTENT(IN) :: itime1
      
    ! mem_chem1
    INTEGER , INTENT(IN) :: antro
    INTEGER , INTENT(IN) :: bioge

    ! rconstants
    REAL    , INTENT(IN) :: pi180


    INTEGER , INTENT(IN) :: nsrc
    TYPE(cycle_emission), INTENT(INOUT) :: emiss_cycle(nsrc)

    INTEGER :: jday,i,j

!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
!    INTEGER :: julday
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------

    REAL :: solfac,tdec,sdec,cdec,declin,d0,d02,dayhr,radlat,cslcsd,snlsnd, &
            gglon,dayhrr,hrangl

    !- local var
    INTEGER it
    REAL time_x,dt_x,cosz
    REAL dcnorma(m2,m3)
    REAL xxx

    !print*,'-----------------------------------------------------' 
    !print*,'getting diurnal/space variable emission rates - time=',time/3600. ;call flush(6)
    !print*,'-----------------------------------------------------'


    dt_x=0.
    time_x=0.
    dcnorma(:,:)=0.


    DO it=1,NINT(86400./dtlt) 

       jday = julday(imonth1,idate1,iyear1) ! don't change the position of this line

       jday = jday + NINT(time_x/86400.)
       !      sdec - sine of declination, cdec - cosine of declination
       declin = -23.5 * COS(6.283 / 365. * (jday + 9)) * pi180
       sdec = SIN(declin)
       cdec = COS(declin)

       ! Find the factor, solfac, to multiply the solar constant to correct
       ! for Earth's varying distance to the sun.

       !d0 = 6.2831853 * float(jday-1) / 365.
       !d02 = d0 * 2.
       !solfac = 1.000110 + 0.034221 * cos (d0) + 0.001280 * sin(d0)  &
       !     + 0.000719 * cos(d02) + 0.000077 * sin(d02)

       ! Find the hour angle, THEN get cosine of zenith angle.

       dayhr = time_x / 3600. + float(itime1/100) + float(MOD(itime1,100)) / 60.

       DO j = ja,jz
          DO i = ia,iz
             radlat = glat(i,j) * pi180
             !IF (lonrad .eq. 0) radlat = centlat(1) * pi180
             !IF (radlat .eq. declin) radlat = radlat + 1.e-5
             cslcsd = COS(radlat) * cdec
             snlsnd = SIN(radlat) * sdec
             !gglon = glon(i,j)
             !IF (lonrad .eq. 0) gglon = centlon(1)
             dayhrr = MOD(dayhr+glon(i,j)/15.+24.,24.)
             hrangl = 15. * (dayhrr - 12.) * pi180
             cosz = snlsnd + cslcsd * COS(hrangl)
             !cosz = min(cosz+1.0E-10, 1.0) 

             dcnorma(i,j)=dcnorma(i,j)+MAX(0.,cosz)

          END DO
       END DO
       time_x=time_x+dtlt
    END DO
    DO j = ja,jz
       DO i = ia,iz
    !- invert dcnorma to save computation time
    !      dcnorma(:,:)=1./(dcnorma(:,:)*dtlt)
           dcnorma(i,j)=1./(dcnorma(i,j)*dtlt)
    ENDDO; ENDDO
    !- transfer the emission cycle to bioge and antro arrays only, for now  
    emiss_cycle(bioge)%dcnorma_inv( :,:) = dcnorma(:,:)
    emiss_cycle(antro)%dcnorma_inv( :,:) = dcnorma(:,:)

    !print*,'dcnorma_INV=',dcnorma(3,4); call flush(6)

  END SUBROUTINE get_diurnal_cycle_normalized
 !------------------------------------------------------------------
  SUBROUTINE interpolation_antro(h_hour,src1,src2,srcn1,srcn2,time1,time2)

	INTEGER :: i, j
	
	REAL               :: newX, g, gTemp, no
	REAL               :: src1,src2,srcn1,srcn2,time1,time2
	REAL, DIMENSION(25) :: f, x, c
	INTEGER  :: INDICE
        REAL :: h_hour,diff

!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
!- Marcelo - curva para SaoPaulo
!	DATA f /0.027,0.016,0.019,0.023,0.027,0.031,0.035,0.040,0.044,0.048,0.051,0.054,0.056,0.057,0.057,&
!                0.057,0.056,0.054,0.051,0.048,0.044,0.040,0.035,0.031,0.027/
!	DATA c /0.018,0.015,0.020,0.025,0.031,0.036,0.041,0.045,0.048,0.050,0.050,0.049,0.047,0.047,0.048,&
!                0.050,0.053,0.056,0.058,0.056,0.052,0.045,0.035,0.026,0.018/ 
!	DATA x /0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24/
!	
!- Madeleine - curva para Lima
	DATA f /0.007, 0.007, 0.003, 0.003, 0.003, 0.010, 0.041, 0.051, 0.082, 0.071, 0.055, 0.049, 0.058, &
	        0.068, 0.063, 0.070, 0.081, 0.082, 0.063, 0.063, 0.026, 0.022, 0.014, 0.010, 0.007/
	DATA c /0.007, 0.007, 0.000, 0.000, 0.000, 0.009, 0.021, 0.034, 0.093, 0.071, 0.059, 0.057, 0.062, &
	        0.076, 0.071, 0.077, 0.077, 0.074, 0.065, 0.062, 0.028, 0.026, 0.017, 0.009, 0.007/ 
	DATA x /0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24/
!--(DMK-CCATT-BRAMS-5.0-OLD)------------------------------------------------------------------
!	DATA f /0.027,0.016,0.019,0.023,0.027,0.031,0.035,0.040,0.044,0.048,0.051,0.054,0.056,0.057,0.057,0.057,0.056,0.054,0.051,0.048,0.044,0.040,0.035,0.031,0.027/
!	DATA c /0.018,0.015,0.020,0.025,0.031,0.036,0.041,0.045,0.048,0.050,0.050,0.049,0.047,0.047,0.048,0.050,0.053,0.056,0.058,0.056,0.052,0.045,0.035,0.026,0.018/ 
!	DATA x /0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24/
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------     
			
                CALL sortZ(24,x,h_hour,1,INDICE)
		
		
		diff = x(INDICE) - h_hour
		
				
		
		IF (diff.GT.0) THEN
		time2 = x(INDICE)
		time1 = x(INDICE-1)
		src2 = c(INDICE)
		src1 = c(INDICE-1)
		srcn2 = f(INDICE)
		srcn1 = f(INDICE-1)
		
		ELSE
		
		time1 = x(INDICE)
		time2 = x(INDICE+1)
		src1 = c(INDICE)
		src2 = c(INDICE+1)
		srcn1 = f(INDICE)
		srcn2 = f(INDICE+1)
		
		ENDIF
	        
		
		
  END SUBROUTINE interpolation_antro
 !------------------------------------------------------------------
  SUBROUTINE sortZ(qtLevels, levels, pLevel, order, outIndex)
	
		INTEGER, INTENT(IN)                    :: qtLevels
		REAL, DIMENSION(qtLevels), INTENT(IN)  :: levels
		REAL, INTENT(IN)                       :: pLevel
		INTEGER, INTENT(IN)                    :: order
		INTEGER, INTENT(OUT) :: outIndex
		
		INTEGER                   :: i,     &
			                     j,     &
			                     minIdx
		REAL                      :: currMin
		REAL, DIMENSION(qtLevels) :: ztemp
	
			
			ztemp = levels
		
			!A SIMPLE SORT THAT RETURNS THE 'order' INDEXES NEAREST OF 'pLevel'
		
			DO j = 1, order
		
				currMin = ABS(ztemp(1) - pLevel)
				minIdx  = 1
			
				DO i = 2, qtLevels
			
					IF( ABS(ztemp(i) - pLevel) .LE. currMin)THEN
				
						currMin = ABS(ztemp(i) - pLevel)
						minIdx  = i
									
					END IF
				END DO				

				outIndex = minIdx
				ztemp(minIdx) = -9E5
					
			END DO	
			 
  END SUBROUTINE sortZ
 !------------------------------------------------------------------
  SUBROUTINE vert_dist_of_volcanic_emission(m1,m2,m3,ia,iz,ja,jz,isrctime,   &
                         nzpmax,dzt,zt,zm,rtgt,topt,dn0,nsrc,     &
			 chem1_src_g,                             &
                         geoge,chem_nspecies,spc_chem_alloc, &
                         src,off,transport,aer1_g,aerosol,   &
                         aer_nspecies,spc_aer_alloc,nmodes,  &
                         volc_mean_g,v_ash) 

    INTEGER , INTENT(IN) :: m1
    INTEGER , INTENT(IN) :: m2
    INTEGER , INTENT(IN) :: m3
    INTEGER , INTENT(IN) :: ia
    INTEGER , INTENT(IN) :: iz
    INTEGER , INTENT(IN) :: ja
    INTEGER , INTENT(IN) :: jz
    INTEGER , INTENT(IN) :: isrctime

    ! grid_dims
    INTEGER , INTENT(IN) :: nzpmax

    ! mem_grid
    REAL    , INTENT(IN) :: dzt(nzpmax),zt(nzpmax),zm(nzpmax)
    REAL    , INTENT(IN) :: rtgt(m2,m3)
    REAL    , INTENT(IN) :: topt(m2,m3)
    REAL    , INTENT(IN) :: dn0(m1,m2,m3)

    ! chem1_list
    INTEGER , INTENT(IN) :: chem_nspecies
    INTEGER , INTENT(IN) :: spc_chem_alloc(6,chem_nspecies)
    INTEGER , INTENT(IN) :: src
    INTEGER , INTENT(IN) :: off
    INTEGER , INTENT(IN) :: transport

    ! aer1_list
    INTEGER , INTENT(IN) :: aer_nspecies
    INTEGER , INTENT(IN) :: spc_aer_alloc(6,nmodes,aer_nspecies)
    INTEGER , INTENT(IN) :: nmodes
    INTEGER , INTENT(IN) :: v_ash
    

    ! mem_chem1
    INTEGER              , INTENT(IN)    :: nsrc
    TYPE(chem1_src_vars) , INTENT(INOUT) :: chem1_src_g(2,nsrc,chem_nspecies)
    INTEGER              , INTENT(IN)    :: geoge

    ! mem_aer1
    TYPE (aer1_vars) , INTENT(INOUT) :: aer1_g(nmodes,aer_nspecies)
    INTEGER          , INTENT(IN)    :: aerosol

    ! mem_volc
    TYPE (volc_mean_vars), INTENT(INOUT) :: volc_mean_g

    ! local var
    REAL :: dz_inv,dz,x1,actual_vent_elev
    INTEGER :: i,j,ksrc,isrc,ispc,imode,itim,k_initial,k_final,kk4,ko,kl,k
    INTEGER :: k_init_umbr
!
    REAL, PARAMETER :: percen_mass_umbrel=.75 
    REAL, PARAMETER :: base_umbrel=.25 
    REAL, ALLOCATABLE :: vert_mass_dist(:)


    ksrc=2          !surface level of emission in the model
    isrc = geoge    !this routine is only for volcanoes
    itim = actual_time_index(isrctime,isrc)


    ALLOCATE(vert_mass_dist(m1))

    !- chemistry section
    DO j=ja,jz
       DO i=ia,iz

       !- if there is not mass to distribute = > cycle
       x1=0.
       DO ispc=1,chem_nspecies
             IF(spc_chem_alloc(src,ispc) == off) CYCLE
              x1=MAX(x1,chem1_src_g(itim,isrc,ispc)%sc_src(1,i,j))
       ENDDO
       DO imode=1,nmodes         
             IF(spc_aer_alloc(src,imode,v_ash) == off ) CYCLE
             x1=MAX(x1,aer1_g(imode,v_ash)%sc_src(1,i,j))
       ENDDO
       !if(maxval( chem1_src_g(itim,isrc,1:chem_nspecies,ng)%sc_src(1,i,j)) <1.e-10) cycle
       IF(x1<1.e-16) CYCLE
       PRINT*,' ==============================================================================' 
       PRINT*,' -> active volcano found at grid point=',i,j
       PRINT*,' -> vent topo plum-heigth=',volc_mean_g%vent_elev(i,j),topt(i,j),volc_mean_g%plum_heigth(i,j)
      
       !- convert units back to 1.e+9 kg/m^2
       !- (this routine is called after 'convert_to_tracer_density' routine
       !-  where the emission field was converted to kg/m^3*1.e9)
       !- so, we need to convert back to create the vertical emission field
       dz	 = rtgt(i,j)/dzt(ksrc) ! dzt=1/(z(k)-z(k-1)) 
       DO ispc=1,chem_nspecies
             IF(spc_chem_alloc(src,ispc) == off) CYCLE
              
	      chem1_src_g(itim,isrc,ispc)%sc_src(1,i,j) = &
              chem1_src_g(itim,isrc,ispc)%sc_src(1,i,j) * dz
       ENDDO
       !-do the same for aerosols
       !-aerosol section (only for volcanoes)
       DO imode=1,nmodes
           IF(spc_aer_alloc(src,imode,v_ash) == off ) CYCLE
           aer1_g(imode,v_ash)%sc_src(1,i,j) =  &
           aer1_g(imode,v_ash)%sc_src(1,i,j) * dz
	ENDDO

       !-
       !- performs the vertical distribution
       
       ! - initial settings 
       k_final    =2     !- volc cloud top    
       k_init_umbr=2     !- volc cloud base
       k_initial  =2     !- volc vent 
       
       !- check vertical position with the actual model topography
       !
       volc_mean_g%vent_elev  (i,j)=MAX( volc_mean_g%vent_elev  (i,j),topt(i,j) )
       
       !- here plum_height is converted to height above sea level:
       !
       !volc_mean_g%plum_heigth(i,j)=MAX( volc_mean_g%plum_heigth(i,j)                             ,topt(i,j) )
        volc_mean_g%plum_heigth(i,j)=MAX( volc_mean_g%plum_heigth(i,j)+volc_mean_g%vent_elev  (i,j),topt(i,j) )
       
       !- 
       IF( ABS(volc_mean_g%plum_heigth(i,j)-volc_mean_g%vent_elev(i,j)) .LE. dz) THEN
             !- case of pseudo 'surface emission', but possibly elevated due mismatch between
             !- model grid-scale topography and the actual elevation of the volcano vent
             DO k=m1-1,2,-1
             !print*,'x',zm(k)*rtgt(i,j)+topt(i,j) , (1.-base_umbrel),volc_mean_g(ng)%plum_heigth(i,j)
        	 IF(zm(k)*rtgt(i,j)+topt(i,j) < volc_mean_g%plum_heigth(i,j))THEN
        	   k_initial=k
        	   EXIT
        	 ENDIF
             ENDDO   
             k_initial=MAX(2,k_initial)
             k_final  =k_initial
	     vert_mass_dist(1:m1) = 0.; vert_mass_dist(k_final) = 1.
	     !print*,'pseudo sfc emis=',k_initial,k_final
       
       ELSE
             !- possibly more than 1 vert emission layer, we will use the 'umbrella' shape for this case:
             !- top volc cloud above sea level   (meters)     
             DO k=m1,2,-1
          
        	  IF(zm(k)*rtgt(i,j)+topt(i,j) < volc_mean_g%plum_heigth(i,j))THEN
        	    k_final=k+1
        	    EXIT
        	  ENDIF
             ENDDO
             k_final=MAX(2,k_final)
             
             IF( k_final == 2 ) THEN                     ! surface emission
        	 k_init_umbr = 2
        	 k_initial = 2
		 vert_mass_dist(1:m1) = 0.;vert_mass_dist(k_final)=1.
		
             ELSE					 ! aerial emission
	     
	    	 !- bottom of volc cloud above sea level
            	 DO k=m1-1,2,-1
            	 !print*,'x',zm(k)*rtgt(i,j)+topt(i,j) , (1.-base_umbrel),volc_mean_g(ng)%plum_heigth(i,j)
            	       IF(zm(k)*rtgt(i,j)+topt(i,j) < (1.-base_umbrel)*volc_mean_g%plum_heigth(i,j))THEN
            	 	k_init_umbr=k
            	 	EXIT
            	      ENDIF
            	 ENDDO   
            	 k_init_umbr=MAX(2,k_init_umbr)
            	 !print*,'kin-kfin=' ,k_init_umbr,k_final
        	 
        	 !- vertical mass distribution (initialization)
        	 vert_mass_dist(1:m1) = 0.
        	 
        	 !- part 1
        	 !- parabolic vertical distribution between k_init_umbr and k_final
        	 kk4 = k_final-k_init_umbr+2
        	 DO ko=2,kk4-1
        	     kl=ko+k_init_umbr-1
        	     vert_mass_dist(kl) = 6.*percen_mass_umbrel* float(ko)/float(kk4)**2 * (1. - float(ko)/float(kk4))
        	 ENDDO
        	 
        	 IF(SUM(vert_mass_dist(1:m1)) .NE. percen_mass_umbrel) THEN
        	     x1= ( percen_mass_umbrel- SUM(vert_mass_dist(1:m1)) )/float(k_final-k_init_umbr+1)
        	     DO ko=k_init_umbr,k_final
        	       vert_mass_dist(ko) = vert_mass_dist(ko)+ x1 !- values between 0 and 1.
        	     ENDDO
        	     !print*,'new mass=',sum(vmd)*100.,x1
        	     !pause
        	 ENDIF
        	 !
        	 !- part 2
        	 !- determine the actual level of the volcanoe vent
		 DO k=m1-1,2,-1
          	     IF(zm(k)*rtgt(i,j)+topt(i,j) < volc_mean_g%vent_elev(i,j))THEN
        	        k_initial=k
        	        EXIT
        	     ENDIF
                 ENDDO   
                 k_initial=MAX(2,k_initial)

		 !- linear detrainment from vent to base of umbrella
        	  DO ko=k_initial,k_init_umbr
        	     !vert_mass_dist(ko)=float(ko)/float(k_init_umbr)
        	     vert_mass_dist(ko)=(log((float(ko)/float((k_initial-1)))))**2./float(k_init_umbr)
                  ENDDO
        	  x1=SUM(vert_mass_dist(2:k_init_umbr))
        	  DO ko=2,k_init_umbr
        	      vert_mass_dist(ko)=(1.-percen_mass_umbrel)*vert_mass_dist(ko)/x1
        	  ENDDO
        	  
        	 ! if(k_init_umbr > 2 ) vert_mass_dist(k_init_umbr)=max(vert_mass_dist(k_init_umbr),vert_mass_dist(k_init_umbr-1))
        	 !- check final mass conservation
        	  IF(SUM(vert_mass_dist(1:m1)) .NE. 1) THEN
        	     x1= ( 1- SUM(vert_mass_dist(1:m1)) )/float(k_final-k_init_umbr+1)
        	     DO ko=k_init_umbr,k_final
        	       vert_mass_dist(ko) = vert_mass_dist(ko)+ x1 !- values between 0 and 1.
        	     ENDDO
        	     !print*,'new mass=',sum(vmd)*100.,x1
        	     !pause
        	  ENDIF
        	 
              ENDIF
        ENDIF
	
	PRINT*,'vert mass distr=',SUM(vert_mass_dist)*100.,i,j,volc_mean_g%plum_heigth(i,j),volc_mean_g%vent_elev(i,j)
      
        DO ispc=1,chem_nspecies
	  
           IF(spc_chem_alloc(src,ispc) == off) CYCLE

           DO k=k_initial,k_final
            
	     !dz= rtgt(i,j)/dzt(ksrc) ! dzt=1/(z(k)-z(k-1)) 
	      dz_inv = dzt(k)/rtgt(i,j)
	     
	     !- convert from kg/m^2  to  density (kg[gas]/m^3*1.e9)
             
	      chem1_src_g(itim,isrc,ispc)%sc_src(k,i,j) = &
              chem1_src_g(itim,isrc,ispc)%sc_src(1,i,j) * dz_inv * vert_mass_dist(k)
	     
	      !IF(chem1_src_g(itim,isrc,ispc)%sc_src(k,i,j)>0.) &
	      !  PRINT*,'ijk sc chem=',i,j,k,vert_mass_dist(k),chem1_src_g(itim,isrc,ispc)%sc_src(k,i,j)
	      
	      
           ENDDO
	ENDDO
	
        !-aerosol section (only for volcanoes)
	DO imode=1,nmodes
         
           IF(spc_aer_alloc(src,imode,v_ash) == off ) CYCLE

       !!    if(imode==1) print*, 'total ash mass kg/m2= ',sum( aer1_g(1:nmodes,v_ash)%sc_src(1,i,j) )*1.e-9
           print*,'----------- imode------ k ----- vert dist ------- mass of ash kg/m^3*1.e9 --------' 
           DO k=k_initial,k_final
            
	     !dz= rtgt(i,j)/dzt(ksrc) ! dzt=1/(z(k)-z(k-1)) 
	      dz_inv = dzt(k)/rtgt(i,j)
	     
	     !- convert from kg/m^2  to  density (kg[gas]/m^3*1.e9)
              aer1_g(imode,v_ash)%sc_src(k,i,j) =  &
              aer1_g(imode,v_ash)%sc_src(1,i,j) * dz_inv * vert_mass_dist(k)

	      !IF(aer1_g(imode,v_ash)%sc_src(k,i,j)>0. .and. imode==1) &
	      PRINT*,imode,k,vert_mass_dist(k),aer1_g(imode,v_ash)%sc_src(k,i,j)
	      
	      
           ENDDO
	ENDDO
	
!tmp-----------
!	PRINT*,' k,vert_mass_dist(k),aer1_g(1,v_ash)%sc_src(k,i,j),height asl' 
!	DO imode=1,nmodes
!           DO k=m1,2,-1
!	        aer1_g(imode,v_ash)%sc_p(k,i,j)=aer1_g(imode,v_ash)%sc_src(k,i,j)/dn0(k,i,j)
!	      
!	        !PRINT*,k,vert_mass_dist(k),aer1_g(1,v_ash)%sc_src(k,i,j),zm(k)*rtgt(i,j)+topt(i,j)
!	        if(imode==nmodes)PRINT*,k,vert_mass_dist(k),aer1_g(imode,v_ash)%sc_p(k,i,j),zm(k)*rtgt(i,j)+topt(i,j)
!	        aer1_g(imode,v_ash)%sc_src(k,i,j) =0.
!           ENDDO
!	ENDDO
!	print*,'sum=' ,sum(vert_mass_dist)
!tmp-----------
           PRINT*,' ==============================================================================' 
       
    !- big-loop X-Y
       ENDDO
    ENDDO
    DEALLOCATE(vert_mass_dist)

 END SUBROUTINE vert_dist_of_volcanic_emission



 SUBROUTINE StoreNamelistFileAtChemSources(oneNamelistFile)
   TYPE(namelistFile), POINTER :: oneNamelistFile
   srcmapfn = oneNamelistFile%srcmapfn
   def_proc_src = oneNamelistFile%def_proc_src
 END SUBROUTINE StoreNamelistFileAtChemSources



  !------------------------------------------------------------------
  ! (DMK) NOT USED
  !
  ! SUBROUTINE get_cosz(m2,m3,ia,iz,ja,jz,glat,glon,cosz)
  !
  !  use mem_grid   , only: imonth1,idate1,iyear1,time,itime1,centlat, &
  !       centlon
  !  use mem_radiate, only: lonrad
  !  use rconstants , only: pi180
  !
  !  implicit none
  !
  !  integer :: m2,m3,ia,iz,ja,jz,jday,i,j,julday
  !
  !  real :: solfac,tdec,sdec,cdec,declin,d0,d02,dayhr,radlat,cslcsd,snlsnd  &
  !       ,gglon,dayhrr,hrangl,dtlt
  !  real, dimension(m2,m3) :: glat,glon,cosz
  !
  !  jday = julday(imonth1,idate1,iyear1)
  !  
  !  jday = jday + nint(time/86400.)
  !  !      sdec - sine of declination, cdec - cosine of declination
  !  declin = -23.5 * cos(6.283 / 365. * (jday + 9)) * pi180
  !  sdec = sin(declin)
  !  cdec = cos(declin)
  !  
  !  ! Find the factor, solfac, to multiply the solar constant to correct
  !  ! for Earth's varying distance to the sun.
  !  
  !  !d0 = 6.2831853 * float(jday-1) / 365.
  !  !d02 = d0 * 2.
  !  !solfac = 1.000110 + 0.034221 * cos (d0) + 0.001280 * sin(d0)  &
  !  !     + 0.000719 * cos(d02) + 0.000077 * sin(d02)
  !  
  !  ! Find the hour angle, THEN get cosine of zenith angle.
  !  
  !  dayhr = time / 3600. + float(itime1/100) + float(mod(itime1,100)) / 60.
  !  
  !  DO j = ja,jz
  !     DO i = ia,iz
  !        radlat = glat(i,j) * pi180
  !        !IF (lonrad .eq. 0) radlat = centlat(1) * pi180
  !        !IF (radlat .eq. declin) radlat = radlat + 1.e-5
  !        cslcsd = cos(radlat) * cdec
  !        snlsnd = sin(radlat) * sdec
  !        !gglon = glon(i,j)
  !        !IF (lonrad .eq. 0) gglon = centlon(1)
  !        dayhrr = mod(dayhr+glon(i,j)/15.+24.,24.)
  !        hrangl = 15. * (dayhrr - 12.) * pi180
  !        cosz(i,j) = snlsnd + cslcsd * cos(hrangl)
  !        
  !     END DO
  !  END DO
  !  
  ! END SUBROUTINE get_cosz

  !------------------------------------------------------------------
  ! (DMK) NOT USED
  ! SUBROUTINE Greg2Jul(h, d, m, y, jd)
  !  
  !   INTEGER , INTENT(IN)  :: d
  !   INTEGER , INTENT(IN)  :: m
  !   INTEGER , INTENT(IN)  :: y
  !   INTEGER , INTENT(IN)  :: h
  !   REAL    , INTENT(OUT) :: jd
  !
  !   jd = ( 1461 * ( y + 4800 + ( m - 14 ) / 12 ) ) / 4 +        &
  !        ( 367 * ( m - 2 - 12 * ( ( m - 14 ) / 12 ) ) ) / 12 -  &
  !        ( 3 * ( ( y + 4900 + ( m - 14 ) / 12 ) / 100 ) ) / 4 + &
  !        d - 32075
  !  
  !   jd = jd + (h/24.)
  !  
  ! END SUBROUTINE Greg2Jul
  !------------------------------------------------------------------

END MODULE chem_sources

!----------------------------------------------------------------------------
! (DMK) NOT USED
!
! !- to be used for interpolation from some data (lat/lon gridded) to model 
! SUBROUTINE interp_init_concentration(ng,n1,n2,n3,glat,glon,dn0,rtgt,dxt,dyt)
!  IMPLICIT NONE
!
!  INTEGER,INTENT(IN) :: ng,n1,n2,n3
!
!  REAL, DIMENSION(n2,n3)    :: dxt,dyt,rtgt,glat,glon
!  REAL, DIMENSION(n1,n2,n3) :: s3p,dn0
!  !variaveis locais
!  REAL ::  latnf,lonnf,dlatr,dlonr  !,ilats,lonn,lats,ilons,latn,lons
!  REAL ::  undef  
!  INTEGER :: i,j,k,qi1,qi2,qj1,qj2,ncount,ii,jj,jc,ic,i1,j1,i2,j2,nzvar !,kk
!  REAL, ALLOCATABLE :: DATAf(:,:)	 ! model concentration
!  REAL, ALLOCATABLE :: tmp(:,:,:)	 !dado dummy temporario
!  !- informacoes da grade dos dados que serao interpolados (input data)
!  REAL, ALLOCATABLE :: api_us(:,:,:),prlat(:,:),prlon(:,:),usdum(:) !api_us= dado inicial a ser interpolado
!  INTEGER,PARAMETER ::nlon=546,nlat=596,n4us=1  !n4us eh usado no caso do dado nao ser somente de um nivel vertical
!  REAL,PARAMETER :: ilatn= 0.090497! delta lat
!  REAL,PARAMETER :: ilonn= 0.091743! delta lon
!  REAL,PARAMETER :: latni=-38.! lat mais ao sul
!  REAL,PARAMETER :: lonni=-85.! lon mais a oeste
!  !real,parameter :: fcu =1.e+6  !=> mg [gas/part] /kg [ar]
!
!  latnf=latni + (nlat-1)*ilatn
!  lonnf=lonni + (nlon-1)*ilonn
!
!  !valor indefinido para o modelo
!  undef = 0.
!
!
!  !grade do dado a ser interpolado
!  ALLOCATE(prlat(nlon,nlat),prlon(nlon,nlat))
!  CALL api_prlatlon(nlon,nlat,prlat,prlon,ilatn,ilonn,latni,lonni)
!  ALLOCATE(api_us(nlon,nlat,n4us),usdum(n4us))
!
!  !- need to adapt the lines below in order to read the data:
!  PRINT*,'--------------------------------------'
!  PRINT*,'Opening initial data='
!  nzvar=22 !number of vert levels of data to be ingested
!  ALLOCATE(tmp(nlon,nlat,nzvar))
!
!  OPEN(2,status='OLD',form='unformatted',access='direct', &
!       recl=4*nlat*nlon*nzvar,file='9_19_2002_262.bin')
!  READ(2) tmp
!  CLOSE(2)
!  !	   call swap32(tmp,nlat*nlon*nzvar)
!  DO i=1,nlon
!     DO j=1,nlat
!	api_us(i,j,1) = tmp (i,j,18)
!	!	   if(tmp(i,j,18).gt. 0.)print*,i,j,tmp (i,j,18)
!     ENDDO
!  ENDDO
!  PRINT*,'--------------------------------------'
!
!  ! dado final        na grade do modelo
!  ALLOCATE(DATAf(n2,n3))
!
!
!  ! loop no dominio do modelo
!  DO j = 1,n3
!     DO i = 1,n2
!
!	! evite pontos fora do dominio da grade de umidade de solo
!	IF(glat(i,j) .LT. latni .OR. glat(i,j) .GT. latnf .OR. &
!	     glon(i,j) .LT. lonni .OR. glon(i,j) .GT. lonnf) THEN
!
!	   DATAf(i,j) = undef
!	   GOTO 1111
!	ENDIF
!
!	CALL interpolacao (glon(i,j),glat(i,j),nlon,nlat,prlat,prlon,i1,i2,ic,j1,j2,jc)
!	IF(ic.GE.0 .AND. jc .GE. 0) THEN
!	   !print*,ic,jc,i,j,ifm
!	   dlonr=0.5*(glon(n2,j)-glon(1,j))/float(n2-1)
!	   dlatr=0.5*(glat(i,n3)-glat(i,1))/float(n3-1)
!	   qi1=INT(dlonr/ilonn+0.5)
!	   qi2=INT(dlonr/ilonn+0.5)
!	   qj1=INT(dlatr/ilatn+0.5)
!	   qj2=INT(dlatr/ilatn+0.5)
!	   !print*,jc-qj1,jc+qj2,ic-qi1,ic+qi2
!	   DO k=1,n4us
!	      ncount = 0
!	      usdum(k)=0.
!
!	      ! print*,'======================================== i j=',i,j
!	      ! print*,jc-qj1,jc+qj2,ic-qi1,ic+qi2
!	      DO jj =MAX(1,jc-qj1),MIN(nlat,jc+qj2)
!		 DO ii = MAX(1,ic-qi1),MIN(nlon,ic+qi2)
!
!		    IF (api_us(ii,jj,k).GE.1.e-5) THEN
!		       ncount = ncount + 1
!		       usdum(k) = usdum(k) + api_us(ii,jj,k)
!
!		       ! print*,'ii jj=',ii,jj
!		       ! print*,'ncount datai soma=',ncount,api_us(ii,jj,k),usdum(k)
!
!		    ENDIF
!		 ENDDO
!	      ENDDO
!	      IF(ncount .GT. 0 ) THEN
!		 DATAf(i,j) = usdum(k) / (float(ncount))
!	      ELSE
!		 DATAf(i,j) = undef
!	      ENDIF
!
!	      ! print*,'FINAL=',dataf(i,j)
!	      ! print*,'========================================'
!
!	      !if(usdum(k) .lt. 1.e-5) then
!	      ! print*,glon(i,j),glat(i,j),usdum(k)
!	      ! stop
!	      !endif
!	   ENDDO
!	   !
!	ENDIF
! 1111	CONTINUE
!
!     ENDDO
!  ENDDO
!
!  DEALLOCATE(api_us,usdum,prlat,prlon)
!
!  OPEN(2,status='UNKNOWN',form='unformatted',access='direct', &
!       recl=4*n2*n3,file='new.bin')
!  WRITE(2) DATAf
!  CLOSE(2)
!
!  
!  CALL azero(n1*n2*n3,s3p)   !zera s3p
!  DO j=1,n3
!     DO i=1,n2
!	DO k=1,7
!	   s3p(k,i,j)  = DATAf(i,j)/dn0(k,i,j)   
!	ENDDO
!     ENDDO
!  ENDDO
!  !
!  DEALLOCATE(tmp,DATAf)
! END SUBROUTINE interp_init_concentration
!-----------------------------------------------------------------------
