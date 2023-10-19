!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine MakeSfcFiles()

  use grid_dims, only: &
       nxpmax,         &
       nypmax

  use mem_mksfc, only: &
       sfcfile_p,      &
       mksfc_scr1,     &
       mksfc_scr2,     &
       mksfc_vt2da,    &
       mksfc_vt2db,    &
       scrx,           &
       npq,            &
       glatp,          &
       glonp,          &
       datq_patch,     &
       datp,           &
       nvndvif,        &
       nvsstf,         &
       glatp,          &
       alloc_sfcfile,  &
       dealloc_sfcfile

  use mem_grid, only:  &
       ngrids,         &
       nnxp,           &
       nnyp,           &
       nnxyp,          &
       npatch,         &
       nzg,            &
       ngridsh,        &
       nxtnest,        &
       deltaxn,        &
       deltayn,        &
       platn,          &
       plonn,          &
       xtn,            &
       ytn,            &
       runtype

  use io_params, only: &
       iupdsst,        &
       isstflg,        &
       itoptflg,       &
       iupdndvi,       &
       ndvifn,         &
       ndviflg,        &
       isoilfn,        &
       isoilflg,       &
       ivegtfn,        &
       ivegtflg

  use node_mod, only: &
       mynum,         &
       mchnum,        &
       master_num,    &
       nodemxp,       &
       nodemyp

  use ReadBcst, only: &
       Broadcast

  ! TEB_SPM
  use teb_spm_start, only: TEB_SPM !INTENT(IN)
  !

  implicit none
  
  integer :: ifm,icm,nvtime,ivtime,ng1,ng2,ng1t,ng2t,ng1s,ng2s
  integer :: isfcerr,itoperr,issterr,indvierr

  ! TEB_SPM
  integer :: ng1f, ng2f
  integer :: ifusoerr
  !

  ! This subroutine makes sure that all surface, topo, sst, and ndvi
  ! files required for the present run exist and are correct. 

  ! The logical choices made are as follows:

  ! For runtype = 'MAKESFC':
  !    Make all surface, topo, sst, and ndvi files for grids 1:NGRIDS.
  ! For runtype = 'INITIAL' or 'MAKEVFILE':
  !    Check for existence and correctness of all surface, topo, sst,
  !       and ndvi files for grids 1:NGRIDS.  Remake all grids if the files
  !       are incorrect for the set of files: topo, sfc/ndvi, sst
  ! For runtype = 'HISTORY':
  !    If NGRIDS > NGRIDSH, check for correctness of all sfcfiles
  !       for grids NGRIDSH+1:NGRIDS.  For any grid (ifm) with an 
  !       incorrect sfcfile, remake the file, requiring that
  !       itoptflg(ifm) = 0.
  !    Check for existence and correctness of all sst and ndvi 
  !       files for grids 1:NGRIDS.  For any grid (ifm) with
  !       incorrect files, make sst and ndvi files.

  integer :: ierr
  integer :: intVec(3)
  character(len=8) :: c0, c1, c2
  character(len=*), parameter :: h="**(MakeSfcFiles)**"

  ! initialize error flags

  isfcerr = 0
  itoperr = 0
  issterr = 0
  indvierr = 0
  ifusoerr = 0

  ! Allocate memory needed for initializing sfcfiles

  if (allocated(sfcfile_p)) then
     deallocate(sfcfile_p, stat=ierr)
     if (ierr /= 0) then
        write(c0,"(i8)") ierr
        call fatal_error(h//" deallocate sfcfile_p fails with stat="//trim(adjustl(c0)))
     end if
  end if
  allocate(sfcfile_p(ngrids), stat=ierr)
  if (ierr /= 0) then
     write(c0,"(i8)") ierr
     write(c1,"(i8)") ngrids
     call fatal_error(h//" allocate sfcfile_p("//trim(adjustl(c1))//&
          ") fails with stat="//trim(adjustl(c0)))
  end if
  do ifm = 1,ngrids
!!$     call alloc_sfcfile(sfcfile_p(ifm), &
!!$          nodemxp(mynum,ifm), nodemyp(mynum,ifm), nzg, npatch)
     call alloc_sfcfile(sfcfile_p(ifm), nnxp(ifm), nnyp(ifm), nzg, npatch)
  end do

  allocate (mksfc_scr1(nxpmax,nypmax), stat=ierr)
  if (ierr /= 0) then
     write(c0,"(i8)") ierr
     write(c1,"(i8)") nxpmax
     write(c2,"(i8)") nypmax
     call fatal_error(h//" allocate scr1("//trim(adjustl(c1))//&
          ","//trim(adjustl(c2))//") fails with stat="//trim(adjustl(c0)))
  end if
  allocate (mksfc_scr2(nxpmax,nypmax), stat=ierr)
  if (ierr /= 0) then
     write(c0,"(i8)") ierr
     write(c1,"(i8)") nxpmax
     write(c2,"(i8)") nypmax
     call fatal_error(h//" allocate scr2("//trim(adjustl(c1))//&
          ","//trim(adjustl(c2))//") fails with stat="//trim(adjustl(c0)))
  end if
  allocate (mksfc_vt2da(nxpmax,nypmax), stat=ierr)
  if (ierr /= 0) then
     write(c0,"(i8)") ierr
     write(c1,"(i8)") nxpmax
     write(c2,"(i8)") nypmax
     call fatal_error(h//" allocate vt2da("//trim(adjustl(c1))//&
          ","//trim(adjustl(c2))//") fails with stat="//trim(adjustl(c0)))
  end if
  allocate (mksfc_vt2db(nxpmax,nypmax), stat=ierr)
  if (ierr /= 0) then
     write(c0,"(i8)") ierr
     write(c1,"(i8)") nxpmax
     write(c2,"(i8)") nypmax
     call fatal_error(h//" allocate vt2db("//trim(adjustl(c1))//&
          ","//trim(adjustl(c2))//") fails with stat="//trim(adjustl(c0)))
  end if
  allocate (scrx(maxval(nnxyp(1:ngrids))*nzg) , stat=ierr)
  if (ierr /= 0) then
     write(c0,"(i8)") ierr
     write(c1,"(i8)") maxval(nnxyp(1:ngrids))*nzg
     call fatal_error(h//" allocate scrx("//trim(adjustl(c1))//&
          ") fails with stat="//trim(adjustl(c0)))
  end if

  if (runtype(1:7) == 'MAKESFC') then

     itoperr = 1
     isfcerr = 1
     issterr = 1
     indvierr = 1

     ng1=1 ; ng2=ngrids      ! sst,ndvi grid bounds
     ng1t=1 ; ng2t=ngrids    ! topo grid bounds
     ng1s=1 ; ng2s=ngrids    ! sfc grid bounds

     if (TEB_SPM==1) then
        ifusoerr = 1
        ng1f=1 ; ng2f=ngrids    ! fuso grid bounds
     end if

  elseif (runtype(1:7) == 'INITIAL' .or. runtype(1:9) == 'MAKEVFILE') then

     ! One process checks sfc files (no side effects)

     if (mchnum == master_num) then
        do ifm = 1,ngrids
           call sfc_check(ifm,isfcerr)
           if(isfcerr == 1) exit
        end do
     end if

     ! One process checks topo files (no side effects)

     if (mchnum == master_num) then
        do ifm = 1,ngrids
           call top_check(ifm,itoperr)
           if(itoperr == 1) exit
        end do
     end if

     ! Check sst files: builds sst file data structure at module io_params;
     ! all processes should build the data structure and return issterr

!!$     if (runtype(1:9)=='MAKEVFILE') then
     do ifm = 1,ngrids
        call SstReadStoreOwnChunk(2,ifm,issterr)
        if(issterr == 1) exit
     end do

     ! Check ndvi files:builds ndvi file data structure at module io_params;
     ! all processes should build the data structure and return indvierr

     do ifm = 1,ngrids
        if(ndviflg(ifm) >= 0) then
           call NdviReadStoreOwnChunk(2,ifm,indvierr)
        end if
        if(indvierr == 1) exit
     end do
!!$     endif

     ! One process checks fuso files 

     if (TEB_SPM==1 .and. mchnum == master_num) then
        do ifm = 1,ngrids
           call fuso_check(ifm,ifusoerr)
           if(ifusoerr == 1) exit
        end do
     end if

     ! broadcast checking results done by a single process

     intVec=(/isfcerr, itoperr, ifusoerr/)

     call Broadcast(intVec(1:3), master_num, "isfcerr,itoperr,ifusoerr")

     isfcerr = intVec(1)
     itoperr = intVec(2)
     ifusoerr = intVec(3)

     ! If we are making ndvi files, we must also make the sfc files (and vice versa)

     if(indvierr == 1) isfcerr = 1
     if(isfcerr == 1) indvierr = 1

     ! return if no check error

     if (isfcerr==0 .and. issterr==0 .and.  &
          itoperr==0 .and.indvierr==0 .and. &
          ifusoerr==0) then
        return
     else
        if (mchnum == master_num) then
           print*, h//' Nonexistent or incorrect surface files for'
           print*, h//'   INITIAL or MAKEVFILE runtype (re)making'
           if(isfcerr == 1 ) print*, h//'   sfc files'
           if(itoperr == 1 ) print*, h//'   top files'
           if(issterr == 1 ) print*, h//'   sst files'
           if(indvierr == 1) print*, h//'   ndvi files'
           if (TEB_SPM==1) then
              if(ifusoerr == 1) print*, h//'   fuso files'
           end if
        end if
     end if

     ng1=1 ; ng2=ngrids
     ng1t=1 ; ng2t=ngrids
     ng1s=1 ; ng2s=ngrids
     if (TEB_SPM==1) then
        ng1f=1 ; ng2f=ngrids
     end if

  elseif (runtype(1:7) == 'HISTORY') then

     !   We can't do this set of checks until we have done a full history start
     !     since we may have to interpolate added grids from coarse grid
     !     fields. This is call from INITLZ on a history start.

     ! Check topo files for added grids. If these are not correct, only
     !   allow remakes if topo is to be interpolated from existing grids.

     do ifm = ngridsh+1,ngrids
        call top_check(ifm,itoperr)
        if (itoperr == 1) then
           if (itoptflg(ifm) /= 0) then
              print*, 'Nonexistent or incorrect TOPFILE for grid ',ifm
              print*, '   which is being added on a history start:'
              print*, '   ITOPTFLG must be set to zero in this case.'
              stop 'added grid - itoptflg'
           end if
        end if
     end do
     ng1t=ngridsh+1 ; ng2t=ngrids

     ! Check sfc files for added grids. Existing grid info is read from
     !   history file.
     do ifm = ngridsh+1,ngrids
        call sfc_check(ifm,isfcerr)
        if(isfcerr == 1) exit
     end do
     ng1s=ngridsh+1 ; ng2s=ngrids

     if (TEB_SPM==1) then
        ! Check fuso files for added grids. Existing grid info is read from
        !   history file.
        do ifm = ngridsh+1,ngrids
           call fuso_check(ifm,ifusoerr)
           if(ifusoerr == 1) exit
        end do
        ng1f=ngridsh+1 ; ng2f=ngrids
     end if

     ! Check sst files for all grids. This is a potentially time-dependent field,
     !   so we need to check if all files are there. If there is a set of files
     !   for a grid that is incomplete, remake all files.

     !   This is somewhat dangerous, since the namelist could have changed 
     !   since the previous run to
     !   specify different data files and we really have no way of 
     !   knowing. We will make the assumption that this didn't occur.
     do ifm = 1,ngrids
        call SstReadStoreOwnChunk(2, ifm, issterr)
        if(issterr == 1) exit
     end do

     ! Do same for NDVI files
     do ifm = 1,ngrids
        if(ndviflg(ifm) >= 0) call NdviReadStoreOwnChunk(2, ifm, indvierr)
        if(indvierr == 1) exit
     end do

     ! If we are making ndvi files, we must also make the sfc files. We will
     !    remake all grids since we may be interpolating from coarser grid.
     if(indvierr == 1 .or. isfcerr == 1) ng1s=1

     ng1=1 ; ng2=ngrids

  end if

  ! If we got here, at least one set of files are bad. Re-make the bad ones.

  if (mchnum==master_num) then

     !------------------------------------------
     !  TOP (topo and roughness) file creation
     if(itoperr == 1) then
        ! do topography, topo roughness on all grids
        call toptnest(ng1t,ng2t)

        do ifm = 1,ngrids
           call top_write(ifm)
        end do
     endif

     ! TEB_SPM
     if (TEB_SPM==1) then
        !------------------------------------------
        !  FUSO  file creation
        if(ifusoerr == 1) then
           ! do FUSO (Local Time) on all grids
           call fusonest(ng1f,ng2f)
           do ifm = 1,ngrids
              call fuso_write(ifm)
           end do
        end if
     end if

     !------------------------------------------
     !  SFC (veg class, patch area, soil type) and NDVI file creation
     !  If we are making ndvi files, we must also make the sfc files
     !  (and vice versa)

     if (isfcerr==1 .or. indvierr==1) then

        ! If iupdndvi = 1, require that:
        !    (1) ndviflg = 1 for a grid that has nxtnest = 0,
        !    (2) ndviflg /= 2 for all other grids.

        if (iupdndvi==1) then
           do ifm=1,ngrids
              if (ndviflg(ifm)/=1 .and. nxtnest(ifm)==0) then
                 print*, 'iupdndvi = 1 and ndviflg /= 1 for grid ', ifm
                 call fatal_error('iupdndvi error')
              end if

              if (ndviflg(ifm)==2) then
                 print*, 'iupdndvi = 1 and ndviflg = 2 for grid ', ifm
                 call fatal_error('iupdndvi error')
              end if
           end do
        end if

        do ifm=ng1s,ng2s

           if (ivegtflg(ifm)==1 .or. isoilflg(ifm)==1 .or. &
                ndviflg(ifm)==1) then

              ! Find size of patch arrays
              call patch_array_size(npq, (xtn(2,ifm)-xtn(1,ifm)),            &
                   ivegtflg(ifm), ivegtfn(ifm), isoilflg(ifm), isoilfn(ifm), &
                   ndviflg(ifm), ndvifn(ifm))

              ! Allocate arrays that need npq dimension
              allocate (glatp(npq,npq,nnxp(ifm),nnyp(ifm)),  &
                   glonp(npq,npq,nnxp(ifm),nnyp(ifm)),       &
                   datq_patch(npq,npq,nnxp(ifm),nnyp(ifm)),  &
                   datp(npq,npq,nnxp(ifm),nnyp(ifm)))

              ! Fill lat-lon
              call patch_latlon(nnxp(ifm), nnyp(ifm),                  &
                   xtn(1,ifm), ytn(1,ifm), deltaxn(ifm), deltayn(ifm), &
                   platn(ifm), plonn(ifm))
           end if

           ! do sfcfile
           call geonest_file(ifm)
           call sfc_write(ifm)

           ! do ndvifile
           if (ndviflg(ifm)==1) then
              call ndvi_read_dataheader(ifm)
              nvtime = nvndvif(ifm)
           elseif (isstflg(ifm)==0) then
              nvtime = nvndvif(nxtnest(ifm))
           else
              nvtime = 1
           end if

           do ivtime=1,nvtime
              call ndvinest(ifm, ivtime)
              call ndvi_write(ifm, ivtime)
           end do

           ! Deallocate arrays that needed npq dimension
           if(allocated(glatp)) deallocate (glatp,glonp,datq_patch,datp)

        end do

     end if

     !------------------------------------------
     !  SST file creation
     if(issterr == 1) then

        ! If iupdsst = 1, require that:
        !    (1) isstflg = 1 for a grid that has nxtnest = 0,
        !    (2) isstflg /= 2 for all other grids.

        if (iupdsst == 1) then
           do ifm = 1,ngrids
              if (isstflg(ifm) /= 1 .and. nxtnest(ifm) == 0) then
                 print*, 'iupdsst = 1 and isstflg /= 1 for grid ', ifm
                 stop 'iupdsst'
              end if

              if (isstflg(ifm) == 2) then
                 print*, 'iupdsst = 1 and isstflg = 2 for grid ', ifm
                 stop 'iupdsst'
              end if
           end do
        end if

        do ifm = ng1,ng2
           ! do sstfile
           if (isstflg(ifm) == 1) then
              call sst_read_dataheader(ifm)
              nvtime = nvsstf(ifm)
           elseif (isstflg(ifm) == 0) then
              nvtime = nvsstf(nxtnest(ifm))
           else
              nvtime = 1
           end if

           do ivtime = 1,nvtime
              call sstnest (ifm,ivtime)
              call sst_write (ifm,ivtime)
           end do
        end do

     end if

  endif

  ! Deallocate memory needed for initializing sfcfiles

  do ifm = 1,ngrids
     call dealloc_sfcfile(sfcfile_p(ifm))
  end do

  deallocate (sfcfile_p, stat=ierr)
  if (ierr /= 0) then
     write(c0,"(i8)") ierr
     write(c1,"(i8)") ngrids
     call fatal_error(h//" deallocate sfcfile_p("//trim(adjustl(c1))//&
          ") fails with stat="//trim(adjustl(c0)))
  end if

  deallocate (mksfc_scr1, stat=ierr)
  if (ierr /= 0) then
     write(c0,"(i8)") ierr
     write(c1,"(i8)") nxpmax
     write(c2,"(i8)") nypmax
     call fatal_error(h//" allocate scr1("//trim(adjustl(c1))//&
          ","//trim(adjustl(c2))//") fails with stat="//trim(adjustl(c0)))
  end if
  deallocate (mksfc_scr2, stat=ierr)
  if (ierr /= 0) then
     write(c0,"(i8)") ierr
     write(c1,"(i8)") nxpmax
     write(c2,"(i8)") nypmax
     call fatal_error(h//" allocate scr2("//trim(adjustl(c1))//&
          ","//trim(adjustl(c2))//") fails with stat="//trim(adjustl(c0)))
  end if
  deallocate (mksfc_vt2da, stat=ierr)
  if (ierr /= 0) then
     write(c0,"(i8)") ierr
     write(c1,"(i8)") nxpmax
     write(c2,"(i8)") nypmax
     call fatal_error(h//" allocate vtda("//trim(adjustl(c1))//&
          ","//trim(adjustl(c2))//") fails with stat="//trim(adjustl(c0)))
  end if
  deallocate (mksfc_vt2db, stat=ierr)
  if (ierr /= 0) then
     write(c0,"(i8)") ierr
     write(c1,"(i8)") nxpmax
     write(c2,"(i8)") nypmax
     call fatal_error(h//" allocate vtdb("//trim(adjustl(c1))//&
          ","//trim(adjustl(c2))//") fails with stat="//trim(adjustl(c0)))
  end if
  deallocate (scrx, stat=ierr)
  if (ierr /= 0) then
     write(c0,"(i8)") ierr
     write(c1,"(i8)") maxval(nnxyp(1:ngrids))*nzg
     call fatal_error(h//" deallocate scrx("//trim(adjustl(c1))//&
          ") fails with stat="//trim(adjustl(c0)))
  end if

end subroutine MakeSfcFiles
