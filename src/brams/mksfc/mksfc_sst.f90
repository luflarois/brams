!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


subroutine sst_read_dataheader(ifm)

  use io_params,only: &
      isstfn
  use mem_mksfc,only: &
      idatevs,        &
      ihourvs,        &
      imonthvs,       &
      iyearvs,        &
      nvsstf,         &
      vsstfil 

  implicit none

  include "files.h"

  integer :: ifm

  integer :: itime,nc
  character(len=f_name_length) :: flnm,line,line2
  character(len=1) :: dummy
  logical :: there
!!$  integer, external :: lastslash
!!$  integer :: ind

  ! Read header file for all sstdata files (all times and locations).  The header
  ! file contains:
  ! first line:       geographic block size (degrees), file size, geographic
  !                   starting point for the dataset (south pole), and offsets
  ! second line:      number of data times (NDTIM)
  ! next NDTIM lines: file prefix, year, month, day, and hour for one data time
  ! last 3 lines:     comments describing the first, second, and following lines
  !                   above

  ! Construct header file name

  flnm=trim(isstfn(ifm))//'HEADER'

  inquire(file=flnm(1:len_trim(flnm)),exist=there)
  if (.not.there) then
     print*,'SSTDATA header file:', trim(flnm), ' for grid ',ifm,' not there.'
     stop 'sst_read_fileheader-1'
  endif

  ! Read this header file

  call rams_f_open(25,flnm(1:len_trim(flnm)),'FORMATTED','OLD','READ',0)
  rewind 25

  ! read number of data times in dataset

  read(25,*) dummy
  read(25,*) nvsstf(ifm)

  if (nvsstf(ifm) <= 0) then
     print*, 'No SST input files found with specified prefix or incorrect header'
     close(25)
     stop 'sst_read_fileheader-2'
  endif

  ! read prefix list and times 

  do itime = 1,nvsstf(ifm)
     read(25,'(A80)') line
     call char_strip_var(line, flnm, line2)
!!$     read (line, *) flnm
!!$     flnm = trim(adjustl(adjustr(flnm)))
!!$     ind = index(line, flnm)
!!$     line2 = line(ind+len(flnm):)
     read (line2,*) iyearvs(itime,ifm), imonthvs(itime,ifm),  &
          idatevs(itime,ifm), ihourvs(itime,ifm)

     vsstfil(itime,ifm)=trim(isstfn(ifm))//trim(flnm)

  enddo

  close(25)

  return
end subroutine sst_read_dataheader

!******************************************************************************

subroutine sstnest(ifm, ivtime)

  use mem_mksfc, only: &
       SFCFILE_P, &
       MKSFC_SCR1, &
       MKSFC_SCR2, &
       NVSSTF, &
       IYEARVS, &
       IMONTHVS, &
       IDATEVS, &
       IHOURVS, &
       VSSTFIL, &
       MKSFC_VT2DA, &
       MKSFC_VT2DB

  use mem_grid, only: &
       NXTNEST, &
       NNXP, &
       NNYP, &
       NXPMAX, &
       NYPMAX, &
       IYEAR1, &
       IMONTH1, &
       IDATE1, &
       IHOUR1

  use io_params, only: &
       ISSTFLG,        &
       ISSTFN

  implicit none
  ! Arguments:
  integer, intent(IN) :: ifm, ivtime
  ! Local Variables:
  integer :: icm

  icm = nxtnest(ifm)

  ! Initialize SEATP and SEATF in subroutine sstinit

  call sstinit(nnxp(ifm),nnyp(ifm),ifm, sfcfile_p(ifm)%seatf(1,1))

  if (icm>=1 .and. isstflg(ifm)==0) then

     ! Interpolate SEATF from coarser grid

     call fillscr(1, nxpmax, nypmax, 1, nnxp(icm), nnyp(icm), 1, 1,  &
          mksfc_scr1, sfcfile_p(icm)%seatf(1,1))
     call eintp(mksfc_scr1, mksfc_scr2, 1, nxpmax, nypmax, 1,        &
          nnxp(ifm), nnyp(ifm), ifm, 2, 't', 0, 0)
     call fillvar(1, nxpmax, nypmax, 1, nnxp(ifm), nnyp(ifm), 1, 1,  &
          mksfc_scr2, sfcfile_p(ifm)%seatf(1,1))

     nvsstf(ifm)                 = nvsstf(icm)
     iyearvs (1:nvsstf(ifm),ifm) = iyearvs (1:nvsstf(ifm),icm)
     imonthvs(1:nvsstf(ifm),ifm) = imonthvs(1:nvsstf(ifm),icm)
     idatevs (1:nvsstf(ifm),ifm) = idatevs (1:nvsstf(ifm),icm)
     ihourvs (1:nvsstf(ifm),ifm) = ihourvs (1:nvsstf(ifm),icm)

  elseif (isstflg(ifm)==1) then

     ! Interpolate SEATF from standard dataset

     call geodat(nnxp(ifm), nnyp(ifm), sfcfile_p(ifm)%seatf(1,1),     &
          isstfn(ifm)(1:len_trim(isstfn(ifm))), vsstfil(ivtime,ifm)(1:len_trim(vsstfil(ivtime,ifm))), mksfc_vt2da, mksfc_vt2db, &
          ifm, 'SST')

  else

     iyearvs(1,ifm) = iyear1; imonthvs(1,ifm) = imonth1
     idatevs(1,ifm) = idate1; ihourvs(1,ifm)  = ihour1        

  endif

  ! If desired, override current values of SEATF with user-defined
  ! changes to subroutine sstinit_user.

  call sstinit_user(nnxp(ifm), nnyp(ifm), ifm, sfcfile_p(ifm)%seatf(1,1))

end subroutine sstnest

! ****************************************************************************

subroutine sstinit(n2,n3,ifm,seatf)

  use mem_leaf

  implicit none
  integer :: n2,n3,ifm,i,j
  real, dimension(n2,n3) :: seatf

  ! Fill the SEATF array with a default value of seatmp.  This 
  ! default is used only when a standard RAMS sst dataset is not used and when 
  ! no overrides to sea temperature are defined in subroutine sstinit_user 
  ! in the file ruser.f90.

  do j = 1,n3
     do i = 1,n2
        seatf(i,j) = seatmp
     enddo
  enddo
  return
end subroutine sstinit


!****************************************************************************

subroutine sst_write(ifm,ivt)

  use mem_grid
  use mem_mksfc,only: &
      idatevs,        &
      ihourvs,        &
      imonthvs,       &
      iyearvs,        &
      nvsstf,         &
      scrx,           &
      sfcfile_p,      &
      vsstfil 
  use io_params,only: &
      sstfpfx

  implicit none

  include "files.h"
  include "UseVfm.h"
  integer :: ifm,ivt,i,j
  integer :: ios
  
  real :: glatr,glonr
  character(len=f_name_length) :: flnm
  character(len=2) :: cgrid
  character(len=*), parameter :: h="**(sst_write)**"
  
  ! Write sst data to sst file for one grid and one time
  write(cgrid,'(a1,i1)') 'g',ifm

    if (useVfm) then
       call makefnam(flnm,sstfpfx,0.,iyearvs(ivt,ifm),imonthvs(ivt,ifm) &
            ,idatevs(ivt,ifm),ihourvs (ivt,ifm)*10000,'W',cgrid,'vfm')
    else
       call makefnam(flnm,sstfpfx,0.,iyearvs(ivt,ifm),imonthvs(ivt,ifm) &
            ,idatevs(ivt,ifm),ihourvs (ivt,ifm)*10000,'W',cgrid,'bin')
    end if

  call xy_ll(glatr,glonr,platn(ifm),plonn(ifm),xtn(1,ifm),ytn(1,ifm))

    if (useVfm) then

       call rams_f_open(25,flnm(1:len_trim(flnm)),'FORMATTED','REPLACE','WRITE',1)
       rewind 25
       
       write(25,99) 999999,2
       write(25,100) iyearvs(ivt,ifm),imonthvs(ivt,ifm) &
            ,idatevs(ivt,ifm),ihourvs (ivt,ifm)
       write(25,101) nnxp(ifm),nnyp(ifm)
       write(25,102) deltaxn(ifm),deltayn(ifm),platn(ifm),plonn(ifm)  &
            ,glatr,glonr
       
99     format(2i8)
100    format(1x,i4.4,2(1x,i2.2),1x,i4.4)
101    format(4i5)
102    format(6f16.5)
       
       call vforec(25,sfcfile_p(ifm)%seatf,nnxp(ifm)*nnyp(ifm),24,scrx,'LIN')

       close(25)


    else

       open(25, action="write", file=trim(flnm), form="unformatted", iostat=ios)
       if (ios /= 0) then
          call fatal_error(h//" opening file "//trim(flnm))
       end if
       rewind 25
       
       write(25) 999999,2
       write(25) iyearvs(ivt,ifm),imonthvs(ivt,ifm) &
            ,idatevs(ivt,ifm),ihourvs (ivt,ifm)
       write(25) nnxp(ifm),nnyp(ifm)
       write(25) deltaxn(ifm),deltayn(ifm),platn(ifm),plonn(ifm)  &
            ,glatr,glonr
       
       write(25) sfcfile_p(ifm)%seatf
       
       close(25)
    end if
  return
end subroutine sst_write
