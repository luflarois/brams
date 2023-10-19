!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


subroutine ndvi_read_dataheader(ifm)

  use mem_mksfc, only: &
       iyearvn, imonthvn, idatevn, ihourvn, vndvifil, nvndvif

  use io_params, only: &
       ndvifn

  use dump, only: &
     dumpMessage

  implicit none

  include "files.h"
  include "constants.f90"

  integer, intent(IN) :: ifm
  integer             :: itime, nc
  character(len=f_name_length)  :: flnm, line, line2
  character(len=1)    :: dummy
  logical             :: there
  character(len=*),parameter :: header='***(ndvi_read_dataheader)***'
  character(len=*),parameter :: version=''
  integer :: ierr

!!$  integer, external :: lastslash
!!$  integer :: ind

  ! Read header file for all ndvidata files (all times and locations).  The header
  ! file contains:
  ! first line:       geographic block size (degrees), file size, geographic
  !                   starting point for the dataset (south pole), and offsets
  ! second line:      number of data times (NDTIM)
  ! next NDTIM lines: file prefix, year, month, day, and hour for one data time
  ! last 3 lines:     comments describing the first, second, and following lines
  !                   above

  ! Construct header file name

  flnm = trim(ndvifn(ifm))//'HEADER'

  inquire(file=flnm(1:len_trim(flnm)), exist=there)
  if (.not.there) then
     print *, 'ndvifn data header file:', trim(flnm), &
          ' for grid ', ifm ,' not there.'
     !call fatal_error('ndvi_read_fileheader-1')
     iErrNumber=dumpMessage(c_tty,c_yes,header,modelVersion,c_fatal,'ndvi_read_fileheader-1')
  endif

  ! Read this header file

  call rams_f_open(25, flnm(1:len_trim(flnm)), 'FORMATTED', 'OLD', 'READ', 0)
  rewind 25

  ! read number of data times in dataset

  read (25,*) dummy
  read (25,*) nvndvif(ifm)

  if (nvndvif(ifm)<=0) then
     print *, &
          'No ndvi input files found with specified prefix or incorrect header'
     close(25)
     !call fatal_error('ndvi_read_fileheader-2')
     iErrNumber=dumpMessage(c_tty,c_yes,header,modelVersion,c_fatal,'ndvi_read_fileheader-2')
  endif

  ! read prefix list and times

  do itime = 1,nvndvif(ifm)
     read (25,'(A80)') line
     call char_strip_var(line, flnm, line2)
!!$     read (line, *) flnm
!!$     flnm = trim(adjustl(adjustr(flnm)))
!!$     ind = index(line, flnm)
!!$     line2 = line(ind+len(flnm):)
!!$
!!$     print *, "DEBUG-ALF:ndvi_read_dataheader:line,line2=", line,line
!!$     call flush(6)

     read (line2,*) iyearvn(itime,ifm), imonthvn(itime,ifm),  &
          idatevn(itime,ifm), ihourvn(itime,ifm)

!!$     nc=lastslash(ndvifn(ifm))
     nc = index(ndvifn(ifm), '/', .true.)
     vndvifil(itime,ifm)=ndvifn(ifm)(1:nc)//trim(flnm)

  enddo

  close(25)

end subroutine ndvi_read_dataheader

!******************************************************************************

subroutine ndvinest(ifm, ivtime)

  use mem_mksfc, only: &
       sfcfile_p, nvndvif, iyearvn, imonthvn, idatevn, ihourvn, vndvifil

  use mem_grid, only: &
       nxtnest, nnxp, nnyp, npatch, nzp, nzg, ipm, jpm, platn, plonn, &
       iyear1, imonth1, idate1, ihour1

  use mem_leaf, only: &
       nvegpat

  use io_params, only: &
       ndviflg, ivegtflg, ivegtfn, isoilflg, isoilfn, ndvifn

  implicit none
  ! Arguments:
  integer, intent(IN) :: ifm
  integer, intent(IN) :: ivtime
  ! Local Variables:
  integer :: icm, i, j, k, ic, jc, ipat

  icm = nxtnest(ifm)

  ! Initialize SEATP and SEATF in subroutine ndviinit

  call ndviinit(nnxp(ifm), nnyp(ifm), npatch, ifm, &
       sfcfile_p(ifm)%veg_ndvif(1,1,1))

  if (icm>=1 .and. ndviflg(ifm)==0) then

     ! Assign NDVIF from coarser grid

     do ipat=2,npatch
        do k=1,nzg
           do j=1,nnyp(ifm)
              do i=1,nnxp(ifm)
                 ic = ipm(i,ifm)
                 jc = jpm(j,ifm)

                 sfcfile_p(ifm)%veg_ndvif(i,j,ipat) = &
                      sfcfile_p(icm)%veg_ndvif(ic,jc,ipat)

              enddo
           enddo
        enddo
     enddo

     nvndvif(ifm)                 = nvndvif(icm)
     iyearvn (1:nvndvif(ifm),ifm) = iyearvn (1:nvndvif(ifm),icm)
     imonthvn(1:nvndvif(ifm),ifm) = imonthvn(1:nvndvif(ifm),icm)
     idatevn (1:nvndvif(ifm),ifm) = idatevn (1:nvndvif(ifm),icm)
     ihourvn (1:nvndvif(ifm),ifm) = ihourvn (1:nvndvif(ifm),icm)

  elseif (ndviflg(ifm)==1) then

     ! Assign NDVIF from standard dataset:

     call landuse_opqr(nnxp(ifm), nnyp(ifm), nzg, npatch, nvegpat,  &
          ivegtflg(ifm), ivegtfn(ifm), isoilflg(ifm), isoilfn(ifm), &
          ndviflg(ifm), ndvifn(ifm), vndvifil(ivtime,ifm),  &
          'ndvi', platn(ifm), plonn(ifm),        &
          sfcfile_p(ifm)%soil_text(1,1,1,1),  &
          sfcfile_p(ifm)%patch_area(1,1,1),   &
          sfcfile_p(ifm)%leaf_class(1,1,1),   &
          sfcfile_p(ifm)%veg_ndvif(1,1,1))
  else

     iyearvn (1,ifm) = iyear1; imonthvn(1,ifm) = imonth1
     idatevn (1,ifm) = idate1; ihourvn (1,ifm) = ihour1

  endif

  ! If desired, override current values of NDVIP and NDVIF in user-specified
  ! changes to subroutine ndviinit_user in the file ruser.f.

  call ndviinit_user(nnxp(ifm), nnyp(ifm), npatch, ifm,  &
       sfcfile_p(ifm)%veg_ndvif(1,1,1))

end subroutine ndvinest

! ****************************************************************************

subroutine ndviinit(n2, n3, npat, ifm, veg_ndvif)
  implicit none
  ! Arguments:
  integer, intent(IN) :: n2, n3, npat, ifm
  real, intent(OUT)   :: veg_ndvif(n2,n3,npat)
  ! Local Variables:
  integer ::i, j, ipat

  ! Fill the veg_ndvif array with a default value of .5.  This default is
  ! used only when a standard RAMS ndvi dataset is not used and when no
  ! overrides to ndvi values are defined in subroutine ndviinit_user in the
  ! file ruser.f.

  do j=1,n3
     do i=1,n2

        veg_ndvif(i,j,1) = 0.0
        veg_ndvif(i,j,2) = 0.5

        do ipat=3,npat
           veg_ndvif(i,j,ipat) = veg_ndvif(i,j,2)
        enddo

     enddo
  enddo

end subroutine ndviinit

!****************************************************************************

subroutine ndvi_write(ifm,ivt)

  use mem_mksfc, only: &
       iyearvn, imonthvn, idatevn, ihourvn, sfcfile_p, scrx

  use mem_grid, only: &
       platn, plonn, xtn, ytn, deltaxn, deltayn, &
       nnxp, nnyp, npatch, nnxyp

  use io_params, only: &
       ndvifpfx

  implicit none

  include "files.h"
  include "UseVfm.h"
  
  integer :: ifm,ivt,ip
  integer :: ios
  real :: glatr,glonr
  character(len=f_name_length) :: flnm
  character(len=2) :: cgrid
  character(len=*), parameter :: h="**(ndvi_write)**"
  
  ! Write ndvi data to ndvi file for one grid and one time

  write(cgrid,'(a1,i1)') 'g',ifm

  if (useVfm) then
     call makefnam(flnm,ndvifpfx,0.,iyearvn(ivt,ifm),imonthvn(ivt,ifm) &
          ,idatevn(ivt,ifm),ihourvn (ivt,ifm)*10000,'N',cgrid,'vfm')
  else
     call makefnam(flnm,ndvifpfx,0.,iyearvn(ivt,ifm),imonthvn(ivt,ifm) &
          ,idatevn(ivt,ifm),ihourvn (ivt,ifm)*10000,'N',cgrid,'bin')
  end if

  call xy_ll(glatr,glonr,platn(ifm),plonn(ifm),xtn(1,ifm),ytn(1,ifm))

  if (useVfm) then
     call rams_f_open(25,flnm(1:len_trim(flnm)),'FORMATTED','REPLACE','WRITE',1)
     rewind 25
     
     write(25,99) 999999,2
99   format(2i8)
     
     write(25,100) iyearvn(ivt,ifm),imonthvn(ivt,ifm) &
          ,idatevn(ivt,ifm),ihourvn (ivt,ifm)
100  format(1x,i4.4,2(1x,i2.2),1x,i4.4)
     
     write(25,101) nnxp(ifm),nnyp(ifm),npatch
101  format(4i5)
     
     write(25,102) deltaxn(ifm),deltayn(ifm),platn(ifm),plonn(ifm)  &
          ,glatr,glonr
102  format(6f16.5)
     
     do ip = 1,npatch
        call vforec(25,sfcfile_p(ifm)%veg_ndvif(1,1,ip),nnxyp(ifm),24,scrx,'LIN')
     enddo
     
     close(25)
     
  else
     open(25, action="write", file=trim(flnm), form="unformatted", iostat=ios)
     if (ios /= 0) then
        call fatal_error(h//" opening file "//trim(flnm))
     end if
     rewind 25
     
     write(25) 999999,2
     write(25) iyearvn(ivt,ifm),imonthvn(ivt,ifm) &
          ,idatevn(ivt,ifm),ihourvn (ivt,ifm)
     write(25) nnxp(ifm),nnyp(ifm),npatch
     write(25) deltaxn(ifm),deltayn(ifm),platn(ifm),plonn(ifm)  &
          ,glatr,glonr
     do ip = 1,npatch
        write(25) sfcfile_p(ifm)%veg_ndvif(:,:,ip)
     enddo
     
     close(25)
  end if
  return
end subroutine ndvi_write
