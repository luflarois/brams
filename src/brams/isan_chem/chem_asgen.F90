!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

!--(DMK-CCATT-INI)----------------------------------------------------------------
SUBROUTINE chem_isan_driver (name_name)
  !--(DMK-CCATT-OLD)----------------------------------------------------------------
  !SUBROUTINE isan_driver (name_name)
  !--(DMK-CCATT-FIM)----------------------------------------------------------------
  use ModDateUtils
  use isan_coms,only: &
       dnref,        &
       guess1st,     &
       hybbot,       &
       hybtop,       &
       idate,        &
       ihour,        &
       iproc_dates,  &
       iproc_flag,   &
       iproc_names,  &
       imonth,       &
       innpr,        &
       ioflgisz,     &
       ioflgvar,     &
       iszstage,     &
       is_grids,     &
       ivrstage,     &
       iyear,        &
       maxisfiles,   &
       maxiy,        &
       maxix,        & 
       maxiz,        &
       maxsigz,      &
       natime,       &
       nfeedvar,     &
       nigrids,      &
       nisn,         &
       npry,         &
       nprx,         &
       nprz,         &
       npdates,      &
       nsigz,        &
       piref,        &
       pi_p,         &
       pi_r,         &
       pi_s,         &
       pi_scrb,      &
       pi_scra,      &
       pi_u,         &
       pi_v,         &
       ps_p,         &
       ps_r,         &
       ps_scrb,      &
       ps_scra,      &
       ps_t,         &
       ps_u,         &
       ps_v,         &
       p_u,          &
       rr_scr1,      &
       rr_scr2,      &
       rr_vt2da,     &
       rs_p,         &
       rs_r,         &
       rs_s,         &
       rs_t,         &
       rs_top,       &
       rs_u,         &
       rs_v,         &
       rs_qual,      &
       rs_sfp,       &
       rs_sft,       &
       rs_slp,       &
       rs_snow,      &
       rs_sst,       &
       rtref,        &
       sigz,         &
       thref,        &
       topsigz,      &
       varpfx,       &
       icFileType

  !--(DMK-CCATT-INI)----------------------------------------------------------------
  !srf-chem
  use chem1_list, only : chemical_mechanism, & ! intent(in)
       spc_name
  use chem_isan_coms                        
  use mem_chem1, only:  CHEM_ASSIM, CHEMISTRY     ! intent(in)
  use aer1_list, only: aer_name
  !srf-chem-end
  !--(DMK-CCATT-FIM)----------------------------------------------------------------

  use mem_grid
  use io_params
  use node_mod, only:  &
       mynum, &
       mchnum, &
       master_num, &
       nmachs
  use mem_aer1, only: aer_assim
  use dump

  implicit none

  include "files.h"
  include "constants.f90"
  include "UseVfm.h"
  
  character(len=*),parameter :: sourceName='chem_asgen.F90'
  character(len=*),parameter :: procedureName='chem_isan_driver'

  character(len=*), intent(IN) :: name_name

  character(len=3) :: csuff
  character(len=f_name_length) :: locfn, locfna
  character(len=16) :: assSpecieName(50),assAerSpecieName(50)

!!$  ! fnames dimensioned by a max number of files possible
!!$  character(len=128) :: fnames(maxisfiles)

  integer, dimension(maxgrds) :: itoptn,iglatn,iglonn
  integer :: ifm,icm,ng,i,k,ifileok

  integer :: recordLen,nvar,irec,j
  character(len=256) :: aerFName

  !--(DMK-CCATT-INI)----------------------------------------------------------------
  !srf-chem
  integer ::  nspc,nm
  real, allocatable :: scratch2D(:,:)
  integer :: ios, ierr
  character(len=8) :: str(10)
  character(len=*), parameter :: h="**(chem_isan_driver)**"
  !srf-chem-end
  !--(DMK-CCATT-FIM)----------------------------------------------------------------

  ! Not necessary because all information is read by Namelist

!!$  ! Read input ISAN namelists
!!$
!!$  call rams_f_open(1,name_name(1:len_trim(name_name)),'FORMATTED','OLD','READ',0)
!!$  call namein_isan(1,'$ISAN_CONTROL',1)
!!$  call namein_isan(1,'$ISAN_ISENTROPIC',1)
!!$  close(1)
  call opspec4

  !  Allocate grid data type since almost all model memory is not used for
  !     ISAN.

  allocate(grid_g(ngrids))
  do ng=1,ngrids
     write (*,fmt='("Proc #",I4.4,", start ISAN grid alloc, grid=",I1,", nnzp=",I3.3,", nnxp=",I3.3,", nnyp=",I3.3)') mynum,ng,nnzp(ng),nnxp(ng),nnyp(ng)
     call nullify_grid(grid_g(ng))
     call alloc_grid(grid_g(ng),nnzp(ng),nnxp(ng),nnyp(ng),ng,if_adap) 
  enddo

  ! Allocate global grid variables data type.
  allocate(oneGlobalGridData(ngrids))
  do ng=1,ngrids
     call nullify_GlobalGridData(oneGlobalGridData(ng))
     call alloc_GlobalGridData(oneGlobalGridData(ng), nnxp(ng), nnyp(ng))
  end do

  ! ISAN runs all the time in sigma-z vertical coordinate. Reset ADAP flag...
  if_adap=0 

  ! Topo is read from surface files.

  do ng=1,ngrids
     call newgrid(ng)
     call TopReadStoreOwnChunk(ng)
  enddo

  !  Setup RAMS hroizontal and vertical grid structure.
  call gridSetup(3) !(2)

  !  Allocate RAMS grid arrays where data analysis will be put 
  !     for output and feedback.  Old way was to use "A"

  !--(DMK-CCATT-INI)----------------------------------------------------------------
  !srf-chem  
  ! calculate the number of species in the 4dda scheme
  nspecies=0
  nspecies_aer_in=0
  if(CHEMISTRY>=0) then
     do nspc=1,total_nspecies
        if(spc_alloc(fdda,nspc) == 1 .and. CHEM_ASSIM == 1) then
           nspecies = nspecies + 1
           assSpecieName(nSpecies)=spc_name(nspc)
        endif
     enddo
     if(CHEM_ASSIM == 1) write(*,fmt='("Proc #",I4.4,", number of species in the 4DDA scheme=",I4.4)') mynum,nspecies 

     do nspc=1,total_speciesAer
        do nm=1,nmodes
           print *,'LFR####### ',nspc,nm,spc_alloc_aer(fdda,nm,nspc)
           if(spc_alloc_aer(fdda,nm,nspc) == 1 .and. aer_ASSIM == 1) then
              nspecies_aer_in=nspecies_aer_in+1
              assAerSpecieName(nspecies_aer_in)=aer_name(nm,nspc)
           endif
        enddo
     enddo
     if(AER_ASSIM == 1) write(*,fmt='("Proc #",I4.4,", number of Aer species in the 4DDA scheme=",I4.4)') mynum,nspecies_aer_in 
  endif
  !--(DMK-CCATT-FIM)----------------------------------------------------------------

  maxix=maxval(nnxp(1:nigrids))   
  maxiy=maxval(nnyp(1:nigrids))   
  maxiz=maxval(nnzp(1:nigrids))   

  do ngrid=1,nigrids  
     allocate(is_grids(ngrid)%rr_u (nnzp(ngrid),nnxp(ngrid),nnyp(ngrid)))
     allocate(is_grids(ngrid)%rr_v (nnzp(ngrid),nnxp(ngrid),nnyp(ngrid)))
     allocate(is_grids(ngrid)%rr_t (nnzp(ngrid),nnxp(ngrid),nnyp(ngrid)))
     allocate(is_grids(ngrid)%rr_r (nnzp(ngrid),nnxp(ngrid),nnyp(ngrid)))
     allocate(is_grids(ngrid)%rr_p (nnzp(ngrid),nnxp(ngrid),nnyp(ngrid)))
     allocate(is_grids(ngrid)%rr_ug (nnzp(ngrid),nnxp(ngrid),nnyp(ngrid)))
     allocate(is_grids(ngrid)%rr_vg (nnzp(ngrid),nnxp(ngrid),nnyp(ngrid)))
     allocate(is_grids(ngrid)%rr_tg (nnzp(ngrid),nnxp(ngrid),nnyp(ngrid)))
     allocate(is_grids(ngrid)%rr_rg (nnzp(ngrid),nnxp(ngrid),nnyp(ngrid)))
     allocate(is_grids(ngrid)%rr_pg (nnzp(ngrid),nnxp(ngrid),nnyp(ngrid)))
     allocate(is_grids(ngrid)%rr_pi0 (nnzp(ngrid),nnxp(ngrid),nnyp(ngrid)))
     allocate(is_grids(ngrid)%rr_th0 (nnzp(ngrid),nnxp(ngrid),nnyp(ngrid)))
     allocate(is_grids(ngrid)%rr_dn0 (nnzp(ngrid),nnxp(ngrid),nnyp(ngrid)))
     allocate(is_grids(ngrid)%rr_dn0u(nnzp(ngrid),nnxp(ngrid),nnyp(ngrid)))
     allocate(is_grids(ngrid)%rr_dn0v(nnzp(ngrid),nnxp(ngrid),nnyp(ngrid)))

     allocate(is_grids(ngrid)%rr_slp (nnxp(ngrid),nnyp(ngrid)))
     allocate(is_grids(ngrid)%rr_sfp (nnxp(ngrid),nnyp(ngrid)))
     allocate(is_grids(ngrid)%rr_sft (nnxp(ngrid),nnyp(ngrid)))
     allocate(is_grids(ngrid)%rr_snow(nnxp(ngrid),nnyp(ngrid)))
     allocate(is_grids(ngrid)%rr_sst (nnxp(ngrid),nnyp(ngrid)))

     !--(DMK-CCATT-INI)----------------------------------------------------------------
     !srf-chem 
     if(CHEM_ASSIM == 1 .and. nspecies>0) then   
        allocate(chem_is_grids(ngrid)%rr_sc (nnzp(ngrid),nnxp(ngrid),nnyp(ngrid),nspecies))
        allocate(chem_is_grids(ngrid)%rr_scg(nnzp(ngrid),nnxp(ngrid),nnyp(ngrid),nspecies))
        allocate(chem_is_grids(ngrid)%rr_sc0(nnzp(ngrid),nnxp(ngrid),nnyp(ngrid),nspecies))
     endif
     if(AER_ASSIM == 1 .and. nspecies_aer_in>0) then   
        allocate(aer_is_grids(ngrid)%rr_sc (nnzp(ngrid),nnxp(ngrid),nnyp(ngrid),nspecies_aer_in))
        allocate(aer_is_grids(ngrid)%rr_scg(nnzp(ngrid),nnxp(ngrid),nnyp(ngrid),nspecies_aer_in))
        allocate(aer_is_grids(ngrid)%rr_sc0(nnzp(ngrid),nnxp(ngrid),nnyp(ngrid),nspecies_aer_in))
     endif


     !srf-chem-end
     !--(DMK-CCATT-FIM)----------------------------------------------------------------
  enddo
  allocate(rr_scr1(maxix*maxiy*maxiz))
  allocate(rr_scr2(maxix*maxiy*maxiz))
  allocate(rr_vt2da(maxix*maxiy))

  ! Do inventory of input file names, determine which times to process
  call chem_ISAN_file_inv(iyear1,imonth1,idate1,itime1,timmax,CHEM_ASSIM, CHEMISTRY)

  write(*,fmt='("Proc #",I4.4,", Total of dates to process=",I4.4," using ",I4.4," procs ",I2.2)') mynum,npdates,nmachs,nhemgrd2

  do natime=1,npdates
     !if(mynum/=natime) cycle 
     if (iproc_flag(natime,1) == 0) cycle
     call date_unmake_big(iyear,imonth,idate,ihour,iproc_dates(natime))
     print*
     print*,'================================================================================================'
     print*,'ISAN processing time: ',natime,' ',iyear,imonth,idate,ihour
     print*,'================================================================================================'
     print*
     ihour=ihour/100
     if(icFileType==1 .or. icFileType==2 .or. icFileType==3 .or. icFileType==4) then
        innpr=iproc_names(natime,5)(1:len_trim(iproc_names(natime,5)))
     else
        innpr=iproc_names(natime,1)
     endif
     if(iproc_flag(natime,2).eq.1.and.guess1st.eq.'PRESS') then

        if(icFileType==1) then
           if(nhemgrd2 == 0) then
              call chem_pressure_stage_grib2(nnxp(1),nnyp(1),nhemgrd2  &
                   ,grid_g(1)%glat(1,1),grid_g(1)%glon(1,1)  &
                   ,grid_g(1)%glat(1,1),grid_g(1)%glon(1,1))
           else
              call chem_pressure_stage_grib2(nnxp(1),nnyp(1),nhemgrd2  &
                   ,grid_g(1)%glat(1,1),grid_g(1)%glon(1,1)  &
                   ,grid_g(nhemgrd2)%glat(1,1)  &
                   ,grid_g(nhemgrd2)%glon(1,1))
           endif
           print*,'after pressure stagep_u(nprx/2,npry/2,1:nprz):',nprx,npry,p_u(nprx/2,npry/2,1:nprz)
        elseif(icFileType==2 .or. icFileType==3) then
#ifdef cdf
           if(nhemgrd2 == 0) then
              call chem_pressure_stage_netCDF(nnxp(1),nnyp(1),nhemgrd2  &
                   ,grid_g(1)%glat(1,1),grid_g(1)%glon(1,1)  &
                   ,grid_g(1)%glat(1,1),grid_g(1)%glon(1,1))
           else
              call chem_pressure_stage_netCDF(nnxp(1),nnyp(1),nhemgrd2  &
                   ,grid_g(1)%glat(1,1),grid_g(1)%glon(1,1)  &
                   ,grid_g(nhemgrd2)%glat(1,1)  &
                   ,grid_g(nhemgrd2)%glon(1,1))
           endif
#endif
           print*,'after pressure stagep_u(nprx/2,npry/2,1:nprz):',nprx,npry,p_u(nprx/2,npry/2,1:nprz)
        elseif(icFileType==4) then
           if(nhemgrd2 == 0) then
              call chem_pressure_stage_grads(nnxp(1),nnyp(1),nnzp(1),nhemgrd2  &
                   ,oneGlobalGridData(1)%global_glat(1,1),oneGlobalGridData(1)%global_glon(1,1)  &
                   ,oneGlobalGridData(1)%global_glat(1,1),oneGlobalGridData(1)%global_glon(1,1))
           else
              call chem_pressure_stage_grads(nnxp(1),nnyp(1),nnzp(1),nhemgrd2  &
                   ,oneGlobalGridData(1)%global_glat(1,1),oneGlobalGridData(1)%global_glon(1,1)  &
                   ,oneGlobalGridData(nhemgrd2)%global_glat(1,1),oneGlobalGridData(nhemgrd2)%global_glon(1,1))
           endif
           print*,'after pressure stagep_u(nprx/2,npry/2,1:nprz):',nprx,npry,p_u(nprx/2,npry/2,1:nprz)
        else
           if(nhemgrd2 == 0) then
              call chem_pressure_stage(nnxp(1),nnyp(1),nhemgrd2  &
                   ,grid_g(1)%glat(1,1),grid_g(1)%glon(1,1)  &
                   ,grid_g(1)%glat(1,1),grid_g(1)%glon(1,1))
           else
              call chem_pressure_stage(nnxp(1),nnyp(1),nhemgrd2  &
                   ,grid_g(1)%glat(1,1),grid_g(1)%glon(1,1)  &
                   ,grid_g(nhemgrd2)%glat(1,1)  &
                   ,grid_g(nhemgrd2)%glon(1,1))
           endif
           print*,'after pressure stagep_u(nprx/2,npry/2,1:nprz):',nprx,npry,p_u(nprx/2,npry/2,1:nprz)
        endif
     endif

     ! Isentropic/sigma-z analysis to all requested RAMS grids

     do ngrid=1,nigrids

        ! Find number of sigma-z levels

        if(guess1st == 'RAMS') then
           topsigz=50.e3 ;   hybtop =50.e3;   hybbot =50.e3
           print*,'************FIRST GUESS RAMS*************'
           print*,'*resetting topsigz,hybbot,hybtop to 50km*'
           print*,'****************************************'
        endif
        do k=1,nnzp(ngrid)
           if(ztn(k,ngrid) > topsigz) goto 100
           sigz(k)=ztn(k,ngrid)
        enddo
        k=nnzp(ngrid)+1
100     continue
        nsigz=k-1


        !         Allocate memory for isentropic analysis
        !         --------------------------------------------------------
        print*,'Allocating RAMS polar/isentropic grid-'  &
             ,ngrid,nnxp(ngrid),nnyp(ngrid),nisn

        allocate(pi_u(nnxp(ngrid),nnyp(ngrid),nisn))
        allocate(pi_v(nnxp(ngrid),nnyp(ngrid),nisn))
        allocate(pi_p(nnxp(ngrid),nnyp(ngrid),nisn))
        allocate(pi_s(nnxp(ngrid),nnyp(ngrid),nisn))
        allocate(pi_r(nnxp(ngrid),nnyp(ngrid),nisn))
        allocate(pi_scra(nnxp(ngrid),nnyp(ngrid),nisn))
        allocate(pi_scrb(nnxp(ngrid),nnyp(ngrid),nisn))

        if(CHEM_ASSIM == 1 .and. nspecies>0) &   
             allocate(pi_sc(nnxp(ngrid),nnyp(ngrid),nisn,nspecies))

        if(aer_ASSIM == 1 .and. nspecies_aer_in>0) &   
             allocate(pi_aer_sc(nnxp(ngrid),nnyp(ngrid),nisn,nspecies_aer_in))

        !         Allocate memory for sigma-z analysis 
        !         --------------------------------------------------------
        print*,'Allocating RAMS polar/sigmaz grid-'  &
             ,ngrid,nnxp(ngrid),nnyp(ngrid),nsigz

        allocate(ps_u(nnxp(ngrid),nnyp(ngrid),nsigz))
        allocate(ps_v(nnxp(ngrid),nnyp(ngrid),nsigz))
        allocate(ps_p(nnxp(ngrid),nnyp(ngrid),nsigz))
        allocate(ps_t(nnxp(ngrid),nnyp(ngrid),nsigz))
        allocate(ps_r(nnxp(ngrid),nnyp(ngrid),nsigz))
        allocate(ps_scra(nnxp(ngrid),nnyp(ngrid),nsigz))
        allocate(ps_scrb(nnxp(ngrid),nnyp(ngrid),nsigz))

        if(CHEM_ASSIM == 1 .and. nspecies>0)&
	     allocate(ps_sc(nnxp(ngrid),nnyp(ngrid),nsigz,nspecies))

        if(AER_ASSIM == 1 .and. nspecies_aer_in>0)&
             allocate(ps_aer_sc(nnxp(ngrid),nnyp(ngrid),nsigz,nspecies_aer_in))   


        !         Allocate memory for surface analysis 
        !         --------------------------------------------------------
        print*,'Allocating RAMS polar/surface grid-'  &
             ,ngrid,nnxp(ngrid),nnyp(ngrid)

        allocate(rs_u(nnxp(ngrid),nnyp(ngrid)))
        allocate(rs_v(nnxp(ngrid),nnyp(ngrid)))
        allocate(rs_p(nnxp(ngrid),nnyp(ngrid)))
        allocate(rs_t(nnxp(ngrid),nnyp(ngrid)))
        allocate(rs_r(nnxp(ngrid),nnyp(ngrid)))
        allocate(rs_s(nnxp(ngrid),nnyp(ngrid)))
        allocate(rs_top(nnxp(ngrid),nnyp(ngrid)))
        allocate(rs_qual(nnxp(ngrid),nnyp(ngrid)))
        allocate(rs_slp(nnxp(ngrid),nnyp(ngrid)))
        allocate(rs_sfp(nnxp(ngrid),nnyp(ngrid)))
        allocate(rs_sft(nnxp(ngrid),nnyp(ngrid)))
        allocate(rs_snow(nnxp(ngrid),nnyp(ngrid)))
        allocate(rs_sst(nnxp(ngrid),nnyp(ngrid)))

        ! !LFR Abrindo grads para escrita
        !            PRINT *,nprx,npry,NPRZ

        !            recordLen=4*nprx*npry
        !            irec=1
        !            write(aerFName,fmt='("aer1_",I4.4,I2.2,I2.2,I6.6)') iyear,imonth,idate  &
        !                 ,ihour*100
        !            open(unit=33,file=trim(aerFName)//'_bci.gra',&
        !                   action='WRITE',status='REPLACE',form='UNFORMATTED',access='DIRECT', &
        !                   recl=recordLen)
        !           do nspc=1,nspecies_aer_in

        ! !LFR Escrevendo Grads
        !              do k=1,nprz
        !               write(33,rec=irec) p_aer_sc(:,:,K,nspc)
        !               irec=irec+1
        !               enddo
        !           enddo
        !           close(33)
        !            open(unit=33,file=trim(aerFName)//'_bci.ctl' &
        !             ,action='WRITE',status='replace',form='FORMATTED')

        !            !writing the name of grads file
        !            write(33,*) 'dset ^'//trim(aerFName)//'_bci.gra'
        !            !writing others infos to ctl
        !            write(33,*) 'undef -0.9990000E+34'
        !            write(33,*) 'title Aerossols In'
        !            write(33,*) 'xdef ',nprx,' linear ',oneGlobalGridData(1)%global_glon(1,1) &
        !              ,oneGlobalGridData(1)%global_glon(2,1)-oneGlobalGridData(1)%global_glon(1,1)
        !            write(33,*) 'ydef ',npry,' linear ',oneGlobalGridData(1)%global_glat(1,1) &
        !              ,oneGlobalGridData(1)%global_glat(1,2)-oneGlobalGridData(1)%global_glat(1,1)
        !            write(33,*) 'zdef ',nprz,'levels',(k,k=1,nprz)
        !            write(33,*) 'tdef 1 linear 00:00z01jan2018     1mo'
        !            write(33,*) 'vars ',nspecies_aer_in
        !            do nvar=1,nspecies_aer_in
        !              write(33,*) assAerSpecieName(nvar),nprz,'99 ',assAerSpecieName(nvar)
        !            enddo
        !            write(33,*) 'endvars'

        !            close(33)

        ! !LFR ---


        ! Do isentropic and sigma-z analysis

        if(iszstage == 1) then
           call chem_isnstage ()	   

           ! !LFR Abrindo grads para escrita
           !            PRINT *,nnxp(ngrid),nnyp(ngrid),nprz,nspecies_aer_in
           !            print *,size(pp_aer_sc,1),size(pp_aer_sc,2),size(pp_aer_sc,3),size(pp_aer_sc,4)

           !            recordLen=4*nnxp(ngrid)*nnyp(ngrid)
           !            irec=1
           !            write(aerFName,fmt='("aer1_",I4.4,I2.2,I2.2,I6.6)') iyear,imonth,idate  &
           !                 ,ihour*100
           !            open(unit=33,file=trim(aerFName)//'_aci.gra',&
           !                   action='WRITE',status='REPLACE',form='UNFORMATTED',access='DIRECT', &
           !                   recl=recordLen)
           !           do nspc=1,nspecies_aer_in

           ! !LFR Escrevendo Grads
           !              do k=1,nprz
           !                 print *,k,nspc,size(pp_aer_sc,1),size(pp_aer_sc,2),size(pp_aer_sc,3),size(pp_aer_sc,4)
           !                 ! if(k<=2) then
           !                 !   do i=1,nnxp(ngrid)
           !                 !     do j=1,nnyp(ngrid)
           !                 !       write(44,fmt='(3(I2.2,1X),E18.6)') i,j,k,pp_aer_sc(i,j,K,nspc); call flush(44)
           !                 !     enddo
           !                 !   enddo
           !                 ! endif
           !                 ! write(33,rec=irec) pp_aer_sc(:,:,K,nspc)
           !                 irec=irec+1
           !              enddo
           !           enddo
           !           close(33)
           !            open(unit=33,file=trim(aerFName)//'_aci.ctl' &
           !             ,action='WRITE',status='replace',form='FORMATTED')

           !            !writing the name of grads file
           !            write(33,*) 'dset ^'//trim(aerFName)//'_aci.gra'
           !            !writing others infos to ctl
           !            write(33,*) 'undef -0.9990000E+34'
           !            write(33,*) 'title Aerossols In'
           !            write(33,*) 'xdef ',nnxp(ng),' linear ',oneGlobalGridData(1)%global_glon(1,1) &
           !              ,oneGlobalGridData(1)%global_glon(2,1)-oneGlobalGridData(1)%global_glon(1,1)
           !            write(33,*) 'ydef ',nnyp(ng),' linear ',oneGlobalGridData(1)%global_glat(1,1) &
           !              ,oneGlobalGridData(1)%global_glat(1,2)-oneGlobalGridData(1)%global_glat(1,1)
           !            write(33,*) 'zdef ',nprz,'levels',(k,k=1,nprz)
           !            write(33,*) 'tdef 1 linear 00:00z01jan2018     1mo'
           !            write(33,*) 'vars ',nspecies_aer_in
           !            do nvar=1,nspecies_aer_in
           !              write(33,*) assAerSpecieName(nvar),nprz,'99 ',assAerSpecieName(nvar)
           !            enddo
           !            write(33,*) 'endvars'

           !            close(33)

           ! !LFR ---




           ! Output isentropic file if desired
           if(CHEMISTRY >= 0 .and. CHEM_ASSIM == 1 .and. ioflgisz == 1) then
              print*,'CHEM Assimilation is not for ready isentropic output.'
              print*,'Please, use ioflgisz=0'
              stop 3548
           endif

           if(ioflgisz == 1)then
              if (useVfm) then
                 write(csuff,'(a1,i1)') 'g',ngrid
                 call makefnam (locfn,varpfx,0,iyear,imonth,idate,  &
                      ihour*100,'I',csuff,'vfm')
                 call rams_f_open(13,locfn(1:len_trim(locfn)),'FORMATTED','REPLACE','READ',iclobber)
                else
                   call makefnam (locfn,varpfx,0,iyear,imonth,idate,  &
                        ihour*100,'I',csuff,'bin')
                   open(13, action="read", file=trim(locfn), form="unformatted", iostat=ios)
                   if (ios /= 0) then
                      write(str(1),"(i8)") ios
                      call fatal_error(h//" opening file "//trim(locfn)//" returns "//trim(adjustl(str(1))))
                   end if
                end if
                call isenio ('OUT',13,nnxp(ngrid),nnyp(ngrid))
                call sigzio ('OUT',13,nnxp(ngrid),nnyp(ngrid))
              close(13)
           endif

        elseif(ivrstage == 1) then
           write(csuff,'(a1,i1)') 'g',ngrid
           if (useVfm) then
              call makefnam(locfn,varpfx,0,iyear,imonth,idate,  &
                   ihour*100,'I',csuff,'vfm')
              call rams_f_open(13,locfn(1:len_trim(locfn)),'FORMATTED','OLD','READ',iclobber)
           else
              call makefnam(locfn,varpfx,0,iyear,imonth,idate,  &
                   ihour*100,'I',csuff,'bin')
              open(13, action="read", file=trim(locfn), form="unformatted", iostat=ios)
              if (ios /= 0) then
                 write(str(1),"(i8)") ios
                 call fatal_error(h//" opening file "//trim(locfn)//" returns "//trim(adjustl(str(1))))
              end if
           end if
           call isenio('IN',13,nnxp(ngrid),nnyp(ngrid))
           call sigzio('IN',13,nnxp(ngrid),nnyp(ngrid))
           close(13)
        endif

        if(ivrstage == 1) then
           call chem_makevarf(ngrid)
        endif

        deallocate(pi_u,pi_v,pi_p,pi_s,pi_r,pi_scra,pi_scrb)
        deallocate(ps_u,ps_v,ps_p,ps_t,ps_r,ps_scra,ps_scrb)
        deallocate(rs_u,rs_v,rs_p,rs_t,rs_r,rs_s,rs_top,rs_qual)
        deallocate(rs_slp,rs_sfp,rs_sft,rs_snow,rs_sst)

        if(CHEM_ASSIM == 1 .and. nspecies>0) deallocate(pi_sc,ps_sc)
        if(AER_ASSIM == 1 .and. nspecies_aer_in>0) deallocate(pi_aer_sc,ps_aer_sc)

     enddo

     ! Do the nesting feedback and write out the "varfiles"
     if(ivrstage == 1) then

        if(nfeedvar == 1 .and. nigrids > 1) then

           ! fill reference states for all grids

           call varfile_refstate(nnzp(1),nnxp(1),nnyp(1)  &
                ,is_grids(1)%rr_t(1,1,1),is_grids(1)%rr_p(1,1,1)  &
                ,is_grids(1)%rr_pi0(1,1,1),is_grids(1)%rr_th0(1,1,1)  &
                ,is_grids(1)%rr_r(1,1,1),is_grids(1)%rr_dn0(1,1,1)  &
                ,is_grids(1)%rr_dn0u(1,1,1),is_grids(1)%rr_dn0v(1,1,1)  &
                ,grid_g(1)%topt(1,1),grid_g(1)%rtgt(1,1)  &
                ,ztn(1,1),ztop,piref(1,1),thref(1,1),dnref(1,1),rtref(1,1))
           is_grids(1)%rr_p(1:nnzp(1),1:nnxp(1),1:nnyp(1)) =  &
                is_grids(1)%rr_p  (1:nnzp(1),1:nnxp(1),1:nnyp(1)) &
                -is_grids(1)%rr_pi0(1:nnzp(1),1:nnxp(1),1:nnyp(1))

           do ifm=2,nigrids
              icm=nxtnest(ifm)
              call fmrefs1d_isan(ifm,icm,maxsigz,nnzp(ifm)  &
                   ,piref,thref,dnref,rtref)
              call fmrefs3d_isan(ifm,icm,nnzp(ifm),nnxp(ifm),nnyp(ifm) &
                   ,nnzp(icm),nnxp(icm),nnyp(icm),maxiz,maxix,maxiy  &
                   ,nnstbot(ifm),nnsttop(ifm),jdim  &
                   ,rr_scr1(1),rr_scr2(1),rr_vt2da(1)  &
                   ,grid_g(ifm)%topt(1,1),grid_g(icm)%topt(1,1) &
                   ,is_grids(icm)%rr_dn0(1,1,1),is_grids(ifm)%rr_dn0(1,1,1) &
                   ,is_grids(icm)%rr_th0(1,1,1),is_grids(ifm)%rr_th0(1,1,1) &
                   ,is_grids(ifm)%rr_pi0(1,1,1),is_grids(ifm)%rr_dn0u(1,1,1) &
                   ,is_grids(ifm)%rr_dn0v(1,1,1),ztn(1,ifm),ztop )

              call fmdn0_isan(ifm,icm,nnzp(ifm),nnxp(ifm),nnyp(ifm) &
                   ,nnzp(icm),nnxp(icm),nnyp(icm),maxiz,maxix,maxiy &
                   ,rr_scr1(1),rr_scr2(1)  &
                   ,grid_g(ifm)%topt(1,1),grid_g(icm)%topt(1,1) &
                   ,is_grids(ifm)%rr_dn0(1,1,1),is_grids(ifm)%rr_dn0u(1,1,1) &
                   ,is_grids(ifm)%rr_dn0v(1,1,1),ztn(1,ifm),ztop )

              is_grids(ifm)%rr_p(1:nnzp(ifm),1:nnxp(ifm),1:nnyp(ifm)) =  &
                   is_grids(ifm)%rr_p  (1:nnzp(ifm),1:nnxp(ifm),1:nnyp(ifm)) &
                   -is_grids(ifm)%rr_pi0(1:nnzp(ifm),1:nnxp(ifm),1:nnyp(ifm))

           enddo

           ! Feed back u,v,pi,t,and rt fields

           do ifm=nigrids,2,-1
              icm=nxtnest(ifm)
              if (icm == 0) cycle
              call chem_varfile_nstfeed(ifm,icm,nnzp(ifm),nnxp(ifm),nnyp(ifm) &
                   ,nnzp(icm),nnxp(icm),nnyp(icm)  &
                   ,nnstbot(icm),nnsttop(icm))
           enddo

           do ifm=1,nigrids
              is_grids(ifm)%rr_p(1:nnzp(ifm),1:nnxp(ifm),1:nnyp(ifm)) =  &
                   is_grids(ifm)%rr_p  (1:nnzp(ifm),1:nnxp(ifm),1:nnyp(ifm)) &
                   +is_grids(ifm)%rr_pi0(1:nnzp(ifm),1:nnxp(ifm),1:nnyp(ifm))
           enddo

        endif


        if(ioflgvar == 1) then
           do ng=1,nigrids
              nxyzp=nnxp(ng)*nnyp(ng)*nnzp(ng)
              nxyp =nnxp(ng)*nnyp(ng)
              write(csuff,'(a1,i1)') 'g',ng

              if (useVfm) then

                 ! vfm file
                 
                 call makefnam (locfn,varpfx,0,iyear,imonth,idate  &
                      ,ihour*100,'V',csuff,'vfm')

                 print *,'*** Open  '//locfn(1:len_trim(locfn))//' for write *** '
                 call rams_f_open (2,locfn(1:len_trim(locfn)),'FORMATTED','REPLACE','WRITE',iclobber)
                 write(2,11) 999999,2
11               format(i7,i3)
                 write(2,10) iyear,imonth,idate,ihour  &
                      ,nnxp(ng),nnyp(ng),nnzp(ng)  &
                      ,platn(ng),plonn(ng),deltaxn(ng)  &
                      ,deltayn(ng),deltazn(ng),dzrat,dzmax
10               format(7i5,7f14.5)
              else

                 ! vfm file
                 
                 call makefnam (locfn,varpfx,0,iyear,imonth,idate  &
                      ,ihour*100,'V',csuff,'bin')

                 print *,'*** Open  '//locfn(1:len_trim(locfn))//' for write *** '
                 open(2, action="write", file=trim(locfn), form="unformatted", iostat=ios)
                 if (ios /= 0) then
                    write(str(1),"(i8)") ios
                    call fatal_error(h//" opening file "//trim(locfn)//" returns "//trim(adjustl(str(1))))
                 end if
                 write(2) 999999,2
                 write(2) iyear,imonth,idate,ihour  &
                      ,nnxp(ng),nnyp(ng),nnzp(ng)  &
                      ,platn(ng),plonn(ng),deltaxn(ng)  &
                      ,deltayn(ng),deltazn(ng),dzrat,dzmax
              end if
                 
              !- chemicam mechanism that this file is prepared for
              if (useVfm) then

                 ! vfm file
                 
                 if(CHEM_ASSIM == 1 .and. nspecies>0) then
                    write(2,*)  trim(chemical_mechanism(1:len_trim(chemical_mechanism)))
                 endif
                 write(*,fmt='("Writing Meteo assimilation for ",I3.3," vars.")') 5 
                 
                 call vforec(2,is_grids(ng)%rr_u(1,1,1),nxyzp,18  &
                      ,rr_scr1(1),'LIN')
                 call vforec(2,is_grids(ng)%rr_v(1,1,1),nxyzp,18  &
                      ,rr_scr1(1),'LIN')
                 call vforec(2,is_grids(ng)%rr_p(1,1,1),nxyzp,18  &
                      ,rr_scr1(1),'LIN')
                 call vforec(2,is_grids(ng)%rr_t(1,1,1),nxyzp,18  &
                      ,rr_scr1(1),'LIN')
                 call vforec(2,is_grids(ng)%rr_r(1,1,1),nxyzp,18  &
                      ,rr_scr1(1),'LIN')

              else
                 
                 if(CHEM_ASSIM == 1 .and. nspecies>0) then
                    write(2)  trim(chemical_mechanism(1:len_trim(chemical_mechanism)))
                 endif
                 write(*,fmt='("Writing Meteo assimilation for ",I3.3," vars.")') 5 

                 write (2) is_grids(ng)%rr_u
                 write (2) is_grids(ng)%rr_v
                 write (2) is_grids(ng)%rr_p
                 write (2) is_grids(ng)%rr_t
                 write (2) is_grids(ng)%rr_r
              end if

              write(*,fmt='("Var: ",I3.3,1X,A16," MaxVal: ",E13.4," MinVal: ",E13.4)') 1, &
                   'U',maxval(is_grids(ng)%rr_u(:,:,:)), &
                   minval(is_grids(ng)%rr_u(:,:,:))
              write(*,fmt='("Var: ",I3.3,1X,A16," MaxVal: ",E13.4," MinVal: ",E13.4)') 2, &
                   'V',maxval(is_grids(ng)%rr_v(:,:,:)), &
                   minval(is_grids(ng)%rr_v(:,:,:))
              write(*,fmt='("Var: ",I3.3,1X,A16," MaxVal: ",E13.4," MinVal: ",E13.4)') 3, &
                   'P',maxval(is_grids(ng)%rr_p(:,:,:)), &
                   minval(is_grids(ng)%rr_p(:,:,:))
              write(*,fmt='("Var: ",I3.3,1X,A16," MaxVal: ",E13.4," MinVal: ",E13.4)') 4, &
                   'T',maxval(is_grids(ng)%rr_t(:,:,:)), &
                   minval(is_grids(ng)%rr_t(:,:,:))
              write(*,fmt='("Var: ",I3.3,1X,A16," MaxVal: ",E13.4," MinVal: ",E13.4)') 5, &
                   'R',maxval(is_grids(ng)%rr_r(:,:,:)), &
                   minval(is_grids(ng)%rr_r(:,:,:))

              if (CHEM_ASSIM == 1 .and. nspecies>0) then
                 print*,'-------------------------------------------------'                   
                 write(*,fmt='("Writing chemistry assimilation for ",I3.3," Species.")') nspecies                   
	         do nspc=1,nspecies
                    if (useVfm) then
                       call vforec(2,chem_is_grids(ng)%rr_sc(1,1,1,nspc),nxyzp,18  &
                            ,rr_scr1(1),'LIN')
                    else
                       write (2) chem_is_grids(ng)%rr_sc(:,:,:,nspc)
                    end if
                    write(*,fmt='("Spc: ",I3.3,1X,A16," MaxVal: ",E13.4," MinVal: ",E13.4)') nspc, &
                         trim(assSpecieName(nspc)),maxval(chem_is_grids(ng)%rr_sc(:,:,:,nspc)), &
                         minval(chem_is_grids(ng)%rr_sc(:,:,:,nspc))

                    if(maxval(chem_is_grids(ng)%rr_sc(:,:,:,nspc)) < 1.e-10) &
                         iErrNumber=dumpMessage(c_tty,c_yes,sourceName,procedureName &
                         ,c_fatal,'wrong dpchem file. Maxval < 1.0e-10')	      

                 end do
	      end if
              if(AER_ASSIM == 1 .and. nspecies_aer_in>0) then

                 print*,'-------------------------------------------------'                   
                 write(*,fmt='("Writing aerosol assimilation for ",I3.3," Species.")') nspecies_aer_in                   
                 do nspc=1,nspecies_aer_in

                    if (useVfm) then
                       call vforec(2,aer_is_grids(ng)%rr_sc(1,1,1,nspc),nxyzp,18  &
                            ,rr_scr1(1),'LIN')
                    else
                       write (2) aer_is_grids(ng)%rr_sc(:,:,:,nspc)
                    end if
                    
                    write(*,fmt='("Spc: ",I3.3,1X,A16," MaxVal: ",E13.4," MinVal: ",E13.4)') nspc, &
                         trim(assAerSpecieName(nspc)),maxval(aer_is_grids(ng)%rr_sc(:,:,:,nspc)), &
                         minval(aer_is_grids(ng)%rr_sc(:,:,:,nspc))

                    if(maxval(aer_is_grids(ng)%rr_sc(:,:,:,nspc)) < 1.e-10) &
                         iErrNumber=dumpMessage(c_tty,c_yes,sourceName,procedureName &
                         ,c_fatal,'wrong dpchem file. Maxval < 1.0e-10')        

                 end do
              end if

              if (useVfm) then
                 call vmissw(is_grids(ng)%rr_slp(1,1),nxyp,rr_vt2da(1),1E30,-1.)
                 call vforec(2,rr_vt2da(1),nxyp,18,rr_scr1(1),'LIN')
                 call vmissw(is_grids(ng)%rr_sfp(1,1),nxyp,rr_vt2da(1),1E30,-1.)
                 call vforec(2,rr_vt2da(1),nxyp,18,rr_scr1(1),'LIN')
                 call vmissw(is_grids(ng)%rr_sft(1,1),nxyp,rr_vt2da(1),1E30,-1.)
                 call vforec(2,rr_vt2da(1),nxyp,18,rr_scr1(1),'LIN')
                 call vmissw(is_grids(ng)%rr_snow(1,1),nxyp,rr_vt2da(1),1E30,-1.)
                 call vforec(2,rr_vt2da(1),nxyp,18,rr_scr1(1),'LIN')
                 call vmissw(is_grids(ng)%rr_sst(1,1),nxyp,rr_vt2da(1),1E30,-1.)
                 call vforec(2,rr_vt2da(1),nxyp,18,rr_scr1(1),'LIN')
              else
                 allocate(scratch2D(nnxp(ng),nnyp(ng)), stat=ierr)
                 if (ierr /= 0) then
                    write(str(1),"(i8)") ierr
                    write(str(2),"(i8)") nnxp(ng)
                    write(str(3),"(i8)") nnyp(ng)
                    call fatal_error(h//" allocate(scratch2d("//&
                         trim(adjustl(str(2)))//","//trim(adjustl(str(3)))//&
                         ") fails with ierr="//trim(adjustl(str(1))))
                 end if
                 call vmissw(is_grids(ng)%rr_slp(1,1),nxyp,scratch2D,1E30,-1.)
                 write (2) scratch2D
                 call vmissw(is_grids(ng)%rr_sfp(1,1),nxyp,scratch2D,1E30,-1.)
                 write (2) scratch2D
                 call vmissw(is_grids(ng)%rr_sft(1,1),nxyp,scratch2D,1E30,-1.)
                 write (2) scratch2D
                 call vmissw(is_grids(ng)%rr_snow(1,1),nxyp,scratch2D,1E30,-1.)
                 write (2) scratch2D
                 call vmissw(is_grids(ng)%rr_sst(1,1),nxyp,scratch2D,1E30,-1.)
                 write (2) scratch2D
                 deallocate(scratch2D, stat=ierr)
                 if (ierr /= 0) then
                    write(str(1),"(i8)") ierr
                    call fatal_error(h//" deallocate scratch2d"//&
                         " fails with ierr="//trim(adjustl(str(1))))
                 end if
              end if
              close(2)
           end do
           call makefnam (locfn,varpfx,0,iyear,imonth,idate  &
                ,ihour*100,'V','$','tag')
           call rams_f_open (2,locfn(1:len_trim(locfn)),'FORMATTED','REPLACE','WRITE',iclobber)
           write(2,*) nigrids
           close(2)

        end if

     end if

  end do

end SUBROUTINE chem_isan_driver
