!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine sfc_check(ifm,ierr)

  use mem_grid, only: &
       platn, plonn, xtn, ytn, deltaxn, deltayn, &
       nnxp, nnyp, nzg, npatch

  use io_params, only: &
       sfcfiles, ivegtflg, isoilflg, nofilflg

  ! This subroutine checks for the existence of a surface file for
  ! grid number ifm, and if it exists, also checks for agreement of
  ! grid configuration between the file and the current model run.
  ! If the file does not exist or does not match grid configuration,
  ! the flag ifileok is returned with a value of 0.  If the file
  ! exists and is ok, ifileok is returned with a value of 1.

  implicit none

  integer :: ifm,ierr

  integer :: lc,isfc_marker,isfc_ver,nsfx,nsfy,nsfzg  &
       ,nsivegtflg,nsisoilflg,nsnofilflg,nspatch
  real ::  sfdx,sfdy,sfplat,sfplon,sflat,sflon,glatr,glonr

  include "files.h"
  include "UseVfm.h"

  character(len=f_name_length) :: flnm
  character(len=2) :: cgrid
  logical there
  integer :: ios
  
  lc=len_trim(sfcfiles)
  write(cgrid,'(a1,i1)') 'g',ifm
  if (useVfm) then
     flnm=trim(sfcfiles)//'-S-'//cgrid//'.vfm'
  else
     flnm=trim(sfcfiles)//'-S-'//cgrid//'.bin'
  end if

  inquire(file=flnm(1:len_trim(flnm)),exist=there)

  if(.not.there) then
     ierr = 1
     print*,'SFCfile:', trim(flnm), ' for grid ',ifm,' not there.'
     print*,'------------------------------------------------'
     return
  endif

  call xy_ll(glatr,glonr,platn(ifm),plonn(ifm),xtn(1,ifm),ytn(1,ifm))

  if (useVfm) then

     call rams_f_open(25,flnm(1:len_trim(flnm)),'FORMATTED','OLD','READ',0)
     read (25,*) isfc_marker,isfc_ver
     read (25,100) nsfx,nsfy,nsfzg,nspatch  &
          ,sfdx,sfdy,sfplat,sfplon,sflat,sflon
     read (25,101) nsivegtflg,nsisoilflg,nsnofilflg
100  format(4i5,2f15.5,4f11.5)
101  format(5i5,2f11.5,i5,2f11.5)
     close (25)
  else

     open(25, action="read", file=trim(flnm), form="unformatted", iostat=ios)
     read (25) isfc_marker,isfc_ver
     read (25) nsfx,nsfy,nsfzg,nspatch  &
          ,sfdx,sfdy,sfplat,sfplon,sflat,sflon
     read (25) nsivegtflg,nsisoilflg,nsnofilflg
     close (25)
  end if

  if (nsfx                       .ne. nnxp(ifm)     .or.  &
       nsfy                       .ne. nnyp(ifm)     .or.  &
       nsfzg                      .ne. nzg           .or.  &
       nspatch                    .ne. npatch        .or.  &
       abs(sfdx-deltaxn(ifm))     .gt. .001          .or.  &
       abs(sfdy-deltayn(ifm))     .gt. .001          .or.  &
       abs(sfplat-platn(ifm))     .gt. .001          .or.  &
       abs(sfplon-plonn(ifm))     .gt. .001          .or.  &
       abs(sflat-glatr)           .gt. .001          .or.  &
       abs(sflon-glonr)           .gt. .001          .or.  &
       nsivegtflg                 .ne. ivegtflg(ifm) .or.  &
       nsisoilflg                 .ne. isoilflg(ifm) .or.  &
       nsnofilflg                 .ne. nofilflg(ifm) ) then

     ierr = 1

     print*,'SFCfile mismatch on grid:',ifm
     print*,'Values: model, file'
     print*,'-------------------'
     print*,'nnxp:',nnxp(ifm),nsfx
     print*,'nnyp:',nnyp(ifm),nsfy
     print*,'deltax:',deltaxn(ifm),sfdx
     print*,'deltay:',deltayn(ifm),sfdy
     print*,'platn:',platn(ifm),sfplat
     print*,'plonn:',plonn(ifm),sfplon
     print*,'SW lat:',glatr,sflat
     print*,'SW lon:',glonr,sflon
     print*,'ivegtflg:',ivegtflg(ifm),nsivegtflg
     print*,'isoilflg:',isoilflg(ifm),nsisoilflg
     print*,'nofilflg:',nofilflg(ifm),nsnofilflg
     print*,'-------------------'

  else

     ierr = 0

  endif

  return
end subroutine sfc_check

!*****************************************************************************

subroutine sfc_write(ifm)

  use mem_mksfc, only: &
       sfcfile_p, scrx

  use mem_grid, only: &
       platn, plonn, xtn, ytn, deltaxn, deltayn, &
       nnxp, nnyp, nzg, npatch, nnxyp

  use io_params, only: &
       sfcfiles, ivegtflg, isoilflg, nofilflg

  implicit none

  include "files.h"
  include "UseVfm.h"
  
  integer :: ifm,ip,k,i,j
  integer :: ios
  real :: glatr,glonr
  character(len=f_name_length) :: flnm
  character(len=2) :: cgrid
  character(len=*), parameter :: h="**(sfc_write)**"
  
  !     write surface characteristics, one file for each grid


  write(cgrid,'(a1,i1)') 'g',ifm

  if (useVfm) then
     flnm=trim(sfcfiles)//'-S-'//cgrid//'.vfm'
  else
     flnm=trim(sfcfiles)//'-S-'//cgrid//'.bin'
  end if


  call xy_ll(glatr,glonr,platn(ifm),plonn(ifm),xtn(1,ifm),ytn(1,ifm))


    if (useVfm) then
       call rams_f_open (25,flnm(1:len_trim(flnm)),'FORMATTED','REPLACE','WRITE',1)
       rewind 25
       
       write(25,99) 999999,3
99     format(2i8)

       write(25,100) nnxp(ifm),nnyp(ifm),nzg,npatch  &
            ,deltaxn(ifm),deltayn(ifm),platn(ifm),plonn(ifm)  &
            ,glatr,glonr
100    format(4i5,2f15.5,4f11.5)

       write(25,101) ivegtflg(ifm),isoilflg(ifm),nofilflg(ifm)
101    format(3i5)


       do ip = 1,npatch
          call vforec(25,sfcfile_p(ifm)%patch_area(1,1,ip),nnxyp(ifm),24,scrx,'LIN')
       enddo

       do ip = 1,npatch
          call vforec(25,sfcfile_p(ifm)%leaf_class(1,1,ip),nnxyp(ifm),24,scrx,'LIN')
       enddo
          
       do ip = 1,npatch
          call vforec(25,sfcfile_p(ifm)%soil_text(1,1,1,ip),nzg*nnxyp(ifm),24,scrx,'LIN')
       enddo

       close(25)

    else

       open(25, action="write", file=trim(flnm), form="unformatted", iostat=ios)
       if (ios /= 0) then
          call fatal_error(h//" opening file "//trim(flnm))
       end if
       
       write(25) 999999,3

       write(25) nnxp(ifm),nnyp(ifm),nzg,npatch  &
            ,deltaxn(ifm),deltayn(ifm),platn(ifm),plonn(ifm)  &
            ,glatr,glonr

       write(25) ivegtflg(ifm),isoilflg(ifm),nofilflg(ifm)

       do ip = 1,npatch
          write(25) sfcfile_p(ifm)%patch_area(:,:,ip)
       enddo

       do ip = 1,npatch
          write(25) sfcfile_p(ifm)%leaf_class(:,:,ip)
       enddo
       
       do ip = 1,npatch
          write(25) sfcfile_p(ifm)%soil_text(:,:,:,ip)
       enddo

       close(25)
    end if

  return
end subroutine sfc_write



subroutine SfcReadStoreOwnChunk(ifm)
  use dump, only: &
    dumpMessage

  use mem_grid, only: &
       npatch, &
       nnxp, &
       nnyp, &
       nzg

  use mem_leaf, only: &
       leaf_g

  use io_params, only: &
       sfcfiles

  use node_mod, only:      &
       mchnum,             &
       mynum,              &
       master_num, nodemxp, nodemyp

  use ReadBcst, only: &
       ReadStoreOwnChunk

  implicit none

  include "files.h"
  include "constants.f90"
  include "UseVfm.h"
  
  integer :: ifm
  integer :: ipat
  integer :: ios
  logical :: there
  character(len=f_name_length) :: flnm
  character(len=2) :: cgrid
  character(len=1) :: dummy
  integer :: ierr
  character(len=8) :: c0
  character(len=*), parameter :: h="**(SfcReadStoreOwnChunk)**"
  character(len=2) :: cipat

  real, pointer :: p2D(:,:), p3D(:,:,:)

  ! master process opens file and skips headers

  if (mchnum == master_num) then
     write(cgrid,'(a1,i1)') 'g',ifm
    if (useVfm) then
       flnm=trim(sfcfiles)//'-S-'//cgrid//'.vfm'
    else
       flnm=trim(sfcfiles)//'-S-'//cgrid//'.bin'
    end if

     inquire(file=flnm(1:len_trim(flnm)),exist=there)

     if(.not.there) then
        !call fatal_error(h//" file "//trim(flnm)//" not there")
        iErrNumber=dumpMessage(c_tty,c_yes,h,modelVersion,c_fatal, &
          " file "//trim(flnm)//" not there")
     end if

    if (useVfm) then

       ! use VFM format

       call rams_f_open(25,flnm(1:len_trim(flnm)),'FORMATTED','OLD','READ',0)
       rewind 25

       ! Skip header records (already been checked)
       read(25,*) dummy
       read(25,*) dummy
       read(25,*) dummy

    else

       ! use binary format
    
       open(25, action="read", file=trim(flnm), form="unformatted", iostat=ios)
       rewind 25

       ! Skip header records (already been checked)
       read(25) dummy
       read(25) dummy
       read(25) dummy

    end if


  end if

  ! deals with patch_area

  do ipat = 1,npatch
     write(cipat,fmt='(I2.2)') ipat
     p2D => leaf_g(ifm)%patch_area(:,:,ipat)
     call ReadStoreOwnChunk(ifm, 25, p2D, "patch_area"//cipat)
  end do
  ! deals with leaf_class

  do ipat = 1,npatch
     p2D => leaf_g(ifm)%leaf_class(:,:,ipat)
     call ReadStoreOwnChunk(ifm, 25, p2D, "leaf_class")
  end do

  ! deals with soil_text

  do ipat = 1,npatch
     p3D => leaf_g(ifm)%soil_text(:,:,:,ipat)
     call ReadStoreOwnChunk(ifm, 25, p3D, nzg, "soil_text")
  end do


  if (mchnum == master_num) then
     close (25)
  end if

end subroutine SfcReadStoreOwnChunk
