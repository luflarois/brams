!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine ndvi_check_header(ifm, flnm, ierr)

  use mem_grid, only: &
       platn, plonn, deltaxn, deltayn, xtn, ytn, &
       nnxp, nnyp, npatch

  ! This subroutine checks for the existence of an ndvi file for
  ! grid number ifm, and if it exists, also checks for agreement of
  ! grid configuration between the file and the current model run.
  ! If the file does not exist or does not match grid configuration,
  ! the flag ierr is returned with a value of 1.  If the file 
  ! exists and is ok, ierr is returned with a value of 0.

  implicit none

  include "files.h"
  include "UseVfm.h"
  
  ! Arguments:
  integer, intent(IN)       :: ifm
  integer, intent(OUT)      :: ierr
  character(len=*), intent(IN) :: flnm
  ! Local Variables:
  integer :: iiyear, iimonth, iidate, iihour,  &
       isfc_marker, isfc_ver, nsfx, nsfy, nsfpat
  integer :: ios
  real :: sfdx, sfdy, sfplat, sfplon, sflat, sflon, glatr, glonr
  character(len=*), parameter :: h="**(ndvi_check_header)**"
  
  ierr = 0

  call xy_ll(glatr, glonr, platn(ifm), plonn(ifm), xtn(1,ifm), ytn(1,ifm))

      if (useVfm) then

       call rams_f_open(25, flnm(1:len_trim(flnm)), 'FORMATTED', 'OLD', 'READ', 0)
       read (25,*) isfc_marker, isfc_ver
       read (25,*) iiyear, iimonth, iidate, iihour
       read (25,*) nsfx, nsfy, nsfpat
       read (25,*) sfdx, sfdy, sfplat, sfplon, sflat, sflon
       close (25)

    else

       open(25, action="read", file=trim(flnm), form="unformatted", iostat=ios)
       if (ios /= 0) then
          call fatal_error(h//" error opening file "//trim(flnm))
       end if
       read (25) isfc_marker, isfc_ver
       read (25) iiyear, iimonth, iidate, iihour
       read (25) nsfx, nsfy, nsfpat
       read (25) sfdx, sfdy, sfplat, sfplon, sflat, sflon
       close (25)

    end if

    if (nsfx/=nnxp(ifm) .or. nsfy/=nnyp(ifm) .or. nsfpat/=npatch .or. &
       abs(sfdx-deltaxn(ifm))>0.001 .or. &
       abs(sfdy-deltayn(ifm))>0.001 .or. &
       abs(sfplat-platn(ifm))>0.001 .or. &
       abs(sfplon-plonn(ifm))>0.001 .or. &
       abs(sflat-glatr)      >0.001 .or. &
       abs(sflon-glonr)      >0.001)     then

     ierr = 1

     print *, 'ndvifile mismatch on grid:',ifm
     print *, 'Values: model, file'
     print *, '-------------------'
     print *, 'nnxp:',   nnxp(ifm), nsfx
     print *, 'nnyp:',   nnyp(ifm), nsfy
     print *, 'npatch:', npatch, nsfpat
     print *, 'deltax:', deltaxn(ifm), sfdx
     print *, 'deltay:', deltayn(ifm), sfdy
     print *, 'platn:',  platn(ifm), sfplat
     print *, 'plonn:',  plonn(ifm), sfplon
     print *, 'SW lat:', glatr, sflat
     print *, 'SW lon:', glonr, sflon
     print *, '-------------------'

   else

     ierr = 0

  endif

end subroutine ndvi_check_header


!=========================================================

subroutine NdviReadStoreOwnChunk(runflag, ifm, ierr)
  use ModDateUtils
  use mem_grid, only: &
       iyear1, imonth1, idate1, itime1, time, timmax, &
       runtype, ngrids, npatch, nnxp, nnyp

  use mem_leaf, only: &
       leaf_g

  use io_params, only: &
       itotdate_ndvi, &
       indviflp, &
       indviflf, &
       indvicyclic, &
       indvicycdata, &
       iupdndvi, &
       ndvifpfx, &
       ndvitime1, &
       ndvitime2, &
       nndvifiles

  use node_mod, only:      &
       mchnum,             &
       master_num

  implicit none
  include "i8.h"
  ! Arguments:
  integer, intent(IN)  :: runflag, ifm
  integer, intent(OUT) :: ierr
  ! Local Variables:
  character(len=14), save  :: totdate_start, totdate_init
  character(len=14)  :: totdate, totdatem
  integer :: iyears, imonths, idates, ihours, nf, ng, i, j, ip
  real    :: timefac_ndvi
  real(kind=8) :: secs_init, secs1, secs2
  character(len=*), parameter :: h="**(NdviReadStoreOwnChunk)**"

  ierr = 0

  if (runflag==1 .or. runflag==2) then   ! Initialization(1) or file check(2)

     ! Inventory all ndvi surface files. 

     call NdviFileInv(ndvifpfx(1:len_trim(ndvifpfx)), ierr)

     if (ierr==1) then
        if (runflag==2) return
        if (runflag==1) then
           call fatal_error(h//' ndvi_read: error on init')
        end if
     end if

     ! Find init and start date

     call date_make_big(iyear1, imonth1, idate1, itime1*100, totdate_init)
     if (runtype=='HISTORY') then
        call date_add_to_big(totdate_init, time, 's', totdate_start)
     elseif (runtype=='INITIAL') then
        totdate_start = totdate_init
     elseif (runtype=='MAKEVFILE') then
        totdate_start = totdate_init
     end if

     ! Do some checks on times

     do ng=1,ngrids
        if (itotdate_ndvi(1,ng)>totdate_start .and. indvicycdata==0) then
           if (mchnum==master_num) then
              print*, h//' initial ndvi file time for grid ', ng, &
                   ' later than beginning of run'
           end if
           ierr = 1
           return
        end if

        if (iupdndvi==1 .and. nndvifiles(ng)==1) then
           if (mchnum==master_num) then
              print*, h//' updating ndvi values but only one ndvi file', &
                   ' for grid', ng
           end if
           ierr = 1
           return
        end if

        call date_add_to_big(totdate_init, timmax, 's', totdatem)
        if (iupdndvi==1 .and. indvicycdata==0 .and. &
             itotdate_ndvi(nndvifiles(ng),ng)<totdatem) then
           if (mchnum==master_num) then
              print*, h//' final ndvi file time for grid ', ng, &
                   'earlier than end of run - making new ndvi files'
           end if
           ierr = 1
           return
        end if
     end do

     if (ierr==1) then
        if (runflag==1) then
           call fatal_error(h//' ndvi_read: error on init time checks')
        end if
     end if

     ! If we are only checking, we're done.

     if (runflag==2) return 

     do ng=1,ngrids

        ! Change the ndvi file dates to the current year. We will increment
        !   when we need to read a future file.
        if (indvicyclic==1) then
           do nf=1,nndvifiles(ng)
              itotdate_ndvi(nf,ng)(1:4) = totdate_start(1:4)
           end do
        end if

        ! Find past time file. The files are ordered in time, so we only 
        !    need to find when start time is greater than a file.
        indviflp(ng) = 0
        do nf=nndvifiles(ng),1,-1
           if (totdate_start>=itotdate_ndvi(nf,ng)) then
              indviflp(ng) = nf
              exit
           end if
        end do

        indviflf(ng) = indviflp(ng) + 1

        ! If we are cyclic, we possibly didn't find a start time later than
        !   a file time. I think it is safe to assume that we will use the 
        !   last file for the start file. Also see if future time will be on 
        !   next cycle.
        if (indvicycdata==1 ) then
           if (indviflp(ng)==0)             indviflp(ng) = nndvifiles(ng)
           if (indviflf(ng)>nndvifiles(ng)) indviflf(ng) = 1
        end if

        ! Read past/future time ndvi field

        call NdviUpdate(0, indviflp(ng))
        if (iupdndvi==1) then
           ! Read future time ndvi field if updating
           call NdviUpdate(1, indviflf(ng))

           ! Compute times as number of seconds past 1 Jan 1900
           call date_abs_secs(totdate_init, secs_init)

           ! Get model time of past file
           totdatem = itotdate_ndvi(indviflp(ng),ng)
           if (indvicyclic==1) then
              ! If month of past file > current month, subtract a year
              if (totdatem(5:6)>totdate_start(5:6))   &
                   call date_add_to_big(totdatem, -365., 'd', totdatem)
           endif
           call date_abs_secs(totdatem, secs1)

           totdatem = itotdate_ndvi(indviflf(ng), ng)
           if (indvicyclic==1) then
              if (totdatem<totdate_start) then
                 ! Future file is in next year. Update all file names for a new year
                 call date_add_to_big(totdatem, 365., 'd', totdatem)
                 do nf=1,nndvifiles(ng)
                    itotdate_ndvi(nf,ng)(1:4) = totdatem(1:4)
                 enddo
              endif
           endif
           call date_abs_secs(totdatem, secs2)

           ndvitime1(ng) = secs1 - secs_init
           ndvitime2(ng) = secs2 - secs_init
        else
           do ip=1,npatch
              leaf_g(ng)%veg_ndvip(:,:,ip) = leaf_g(ng)%veg_ndvif(:,:,ip)
           enddo
           ndvitime1(ng) = 0.
           ndvitime2(ng) = 0.
        endif

        ! Fill "current" time ndvi values
        if (iupdndvi==0) then
           timefac_ndvi = 0.
        else
           timefac_ndvi = (time - ndvitime1(ng))/(ndvitime2(ng) - ndvitime1(ng))
        endif

        do ip=1,npatch
           leaf_g(ng)%veg_ndvic(:,:,ip) = leaf_g(ng)%veg_ndvip(:,:,ip) + &
                (leaf_g(ng)%veg_ndvif(:,:,ip) - &
                leaf_g(ng)%veg_ndvip(:,:,ip))*timefac_ndvi
        end do

     enddo

     return

  elseif (runflag==3) then   ! Runtime file increment

     call date_add_to_big(totdate_init, time, 's', totdate)

     if (time>=ndvitime2(ifm) ) then

        ! Update ndvi fields
        indviflp(ifm) = indviflf(ifm)
        indviflf(ifm) = indviflp(ifm) + 1

        ! Compute times as number of seconds past 1 Jan 1900
        !   If cyclic, modify file date 
        call date_abs_secs(totdate_init, secs_init)
        call date_add_to_big(totdate_init, time, 's', totdate)

        totdatem = itotdate_ndvi(indviflp(ifm), ifm)
        call date_abs_secs(totdatem, secs1)

        ! Need to deal with cyclic. 
        if (indvicyclic==1 .and. indviflf(ifm)>nndvifiles(ifm)) then
           indviflf(ifm) = 1
           ! Update all file names for a new year
           totdatem = itotdate_ndvi(indviflf(ifm), ifm)
           call date_add_to_big(totdatem, 365., 'd', totdatem)
           do nf=1,nndvifiles(ifm)
              itotdate_ndvi(nf,ifm)(1:4) = totdatem(1:4)
           enddo
        endif

        totdatem = itotdate_ndvi(indviflf(ifm), ifm)
        call date_abs_secs(totdatem, secs2)

        ndvitime1(ifm) = secs1 - secs_init
        ndvitime2(ifm) = secs2 - secs_init

        ! Finally read the actual field           
        call NdviUpdate(1, indviflf(ifm))

        if (mchnum==master_num) then
           print *, h//' Switched ndvi files:',ndvitime1(ifm),ndvitime2(ifm), &
                secs1, secs2, secs_init
        end if

     else
        return
     endif
  endif
end subroutine NdviReadStoreOwnChunk







subroutine NdviFileInv (sfilin, ierr)
  use ModDateUtils
  use mem_grid, only: &
       iyear1, imonth1, idate1, itime1, ihour1, ngrids, &
       maxndvifiles, nxtnest

  use mem_mksfc, only: &
      nvsstf,          &
      idatevn,         &
      ihourvn,         &
      imonthvn,        &
      iyearvn,         &
      nvndvif     

  use io_params, only: &
       indvicyclic,    & ! intent(out)
       indvicycdata,   & ! intent(out)
       nndvifiles,     & ! intent(out)
       fnames_ndvi,    & ! intent(out)
       itotdate_ndvi,  & ! intent(out)
       iupdndvi,       & ! intent(in)
       ndviflg,        & 
       ISSTFLG

  use node_mod, only:      &
       mchnum,             &
       master_num

  use ReadBcst, only: &
       Broadcast

  implicit none

  include "files.h"
  include "UseVfm.h"
  
  ! Arguments:
  character(len=*), intent(IN) :: sfilin
  integer, intent(OUT)         :: ierr
  ! Local Variables:
  integer :: nc, nf, lnf, nftot, ng, icm, ifm
  integer :: inyear, inmonth, indate, inhour
  character(len=f_name_length), dimension(maxndvifiles) :: fnames
  character(len=f_name_length) :: rams_filelist_arg
  character(len=14)  :: itotdate
  character(len=1)  :: cgrid
  character(len=2)  :: cgrid2
  real(kind=8) :: secs_init, secs_file
  integer :: lastChar
  integer :: lenFnames
  integer :: lenItot
  integer :: sizeIntVec
  integer :: sizeCharVec
  integer :: ierr2
  integer :: nvtime,ivtime,indice
  integer, allocatable :: intVec(:)
  character, allocatable :: charVec(:)
  character(len=8) :: c0, c1
  character(len=*), parameter :: h="**(NdviFileInv)**"
  character(len=f_name_length) :: flnm
  logical :: there

  ! size and length of data to be broadcasted

  sizeIntVec = ngrids+3
  lenFnames  = len(fnames_ndvi)
  lenItot    = len(itotdate_ndvi)

  ! only one process does file inventory

  if (mchnum==master_num) then

     ! io_params initialization

     nndvifiles   = 0
!     fnames_ndvi= ""
!     itotdate_ndvi= ""
     indvicyclic  = 0
     indvicycdata = 0
     ierr         = 0

     ! Get abs seconds of run start

     call date_abs_secs2(iyear1, imonth1, idate1, itime1*100, secs_init)

     ! Go through ndvi files and make inventory. We unfortunately have to do this
     !   for all grids.

     indvicyclic  = 0
     indvicycdata = 0

     indice = 1

     do ng=1,ngrids

        write(cgrid2,'(a1,i1)') 'g',ng
        nftot = -1
        write (cgrid,'(i1)') ng

        !!$ rams_filelist_arg = trim(sfilin)//'-N-*-g'//cgrid//'.vfm'

        if (ndviflg(ng)==1) then
           call ndvi_read_dataheader(ng)
           nvtime = nvndvif(ng)
        elseif (isstflg(ng)==0) then
           nvtime = nvndvif(nxtnest(ng))
        else
           nvtime = 1
        end if
        
        ifm = ng
        icm = nxtnest(ifm)

        do ivtime=1,nvtime
           !!call ndvinest(ng, ivtime)
           if (ng>=1 .and. ndviflg(ifm)==0) then
              nvndvif(ifm)                 = nvndvif(icm)
              iyearvn (1:nvndvif(ifm),ifm) = iyearvn (1:nvndvif(ifm),icm)
              imonthvn(1:nvndvif(ifm),ifm) = imonthvn(1:nvndvif(ifm),icm)
              idatevn (1:nvndvif(ifm),ifm) = idatevn (1:nvndvif(ifm),icm)
              ihourvn (1:nvndvif(ifm),ifm) = ihourvn (1:nvndvif(ifm),icm)
           endif
           
             if (useVfm) then
                call makefnam(flnm, sfilin, 0., &
                     iyearvn(ivtime,ng), imonthvn(ivtime,ng), &
                     idatevn(ivtime,ng), ihourvn (ivtime,ng)*10000, &
                     'N', cgrid2, 'vfm')
             else
                call makefnam(flnm, sfilin, 0., &
                     iyearvn(ivtime,ng), imonthvn(ivtime,ng), &
                     idatevn(ivtime,ng), ihourvn (ivtime,ng)*10000, &
                     'N', cgrid2, 'bin')
             end if

           inquire(file=flnm(1:len_trim(flnm)), exist=there)

           if (there) then
              fnames(indice) = trim(flnm)
              indice = indice + 1 
           endif
        end do

        nftot = indice - 1
        
        !!$ call RAMS_filelist(fnames, rams_filelist_arg, nftot)
        
        if (nftot<=0) then
           print *, 'No ndvi files for grid '//cgrid
           ierr = 1
           exit
        end if

        if (nftot>maxndvifiles) then
           call fatal_error(h//' too many ndvi files')
        end if


        ! We need to see if the data is to be considered "cyclic". Count how many files
        !  have the year 0000. If all of them do not, we probably have an error.
        !  indvicycdata = 1 when data is cyclic
        !  indvicyclic =1 when we are have cyclic data and will be updating in time
        
        nndvifiles(ng) = 0
        do nf=1,nftot
           lnf = len_trim(fnames(nf))
           read (fnames(nf)(lnf-23:lnf-7), "(i4,1x,i2,1x,i2,1x,i6)") &
                inyear, inmonth, indate, inhour

           ! Check the file headers
           call ndvi_check_header(ng,fnames(nf)(1:len_trim(fnames(nf))),ierr)
           if (ierr/=0) exit

           if (inyear==0) indvicycdata = indvicycdata + 1

           call date_make_big(inyear,inmonth,indate,inhour,itotdate)

           nndvifiles(ng)                   = nndvifiles(ng) + 1
           fnames_ndvi(nndvifiles(ng),ng)   = fnames(nf)
           itotdate_ndvi(nndvifiles(ng),ng) = itotdate

        end do

        if (ierr==0) then
           call RAMS_dintsort(nndvifiles(ng), itotdate_ndvi(1,ng), &
                fnames_ndvi(1,ng))

           !  start printing section
           open(unit=22,file='brams.log',position='append',action='write')
           write(unit=22,fmt='(A)')  ' '
           write(unit=22,fmt='(A)')  ' '
           write(unit=22,fmt='(A)')  ' '
           write(unit=22,fmt='(A)')  '-------------------------------------------------------------'
           write(unit=22,fmt='(A)')  '-----------  NDVI Input File Inventory: Grid '//cgrid
           write(unit=22,fmt='(A)')  '-------------------------------------------------------------'
           do nf=1,nndvifiles(ng)
              write(unit=22,fmt='(A,A,A)')   itotdate_ndvi(nf,ng), '   ', trim(fnames_ndvi(nf,ng))
           end do
           write(unit=22,fmt='(A)') '------------------------------------------------------'
           close(unit=22)
        else
           exit
        end if
     end do

     ! Check the cyclic condition. Only relevant if we are updating in time.
     !   WE ARE ONLY ALLOWING CYCLIC ON ALL GRIDS

     if (ierr==0) then
        if (indvicycdata>0) then
           if (indvicycdata/=sum(nndvifiles(1:ngrids))) then
              call fatal_error(h// ' All ndvi surface files do not have year 0000;'//&
                   ' This confuses the gods and can not occur')
           end if
           indvicycdata = 1
        else
           indvicycdata = 0
        end if
     end if

     ! Set the main cyclic flag. Only relevant if we are updating in time.

     if (iupdndvi==1 .and. indvicycdata==1) indvicyclic = 1
  end if

  ! allocate broadcast area

  allocate(intVec(sizeIntVec), stat=ierr2)
  if (ierr2 /= 0) then
     write(c0,"(i8)") ierr2
     write(c1,"(i8)") sizeIntVec
     call fatal_error(h//" allocate intVec("//trim(adjustl(c1))//&
          ") fails with stat="//trim(adjustl(c0)))
  end if

  ! master process gathers data for broadcasting

  if (mchnum==master_num) then
     intVec(1:3)          = (/ierr, indvicyclic, indvicycdata /)
     intVec(4:sizeIntVec) = nndvifiles(1:ngrids)
  end if

  ! broadcast integer data to remaining processes

  call Broadcast(intVec, master_num, "intVec")

  ! scatter broadcasted data

  ierr         = intVec(1)
  indvicyclic  = intVec(2)
  indvicycdata = intVec(3)

  nndvifiles(:)        = 0
  nndvifiles(1:ngrids) = intVec(4:sizeIntVec)
!!$  nndvifiles(ngrids+1:)=0

  ! deallocate broadcast area

  deallocate(intVec, stat=ierr2)
  if (ierr2 /= 0) then
     write(c0,"(i8)") ierr2
     call fatal_error(h//" deallocate intVec fails with stat="//trim(adjustl(c0)))
  end if

  ! broadcast remaining data only if inventory succeeds

  if (ierr==0) then

     ! allocate broadcast area

     sizeCharVec = sum(nndvifiles(1:ngrids))*(lenFnames + lenItot)
     allocate(charVec(sizeCharVec), stat=ierr2)
     if (ierr2/=0) then
        write(c0,"(i8)") ierr2
        write(c1,"(i8)") sizeCharVec
        call fatal_error(h//" allocate charVec("//trim(adjustl(c1))//&
             ") fails with stat="//trim(adjustl(c0)))
     end if

     ! master process prepares broadcast data

     if (mchnum==master_num) then
        lastChar = 0
        do ng=1,ngrids
           do nf=1,nndvifiles(ng)
              do nc=1,lenFnames
                 charVec(lastChar+nc) = fnames_ndvi(nf,ng)(nc:nc)
              end do
              lastChar = lastChar+lenFnames
              do nc=1,lenItot
                 charVec(lastChar+nc) = itotdate_ndvi(nf,ng)(nc:nc)
              end do
              lastChar = lastChar + lenItot
           end do
        end do
     end if

     ! broadcast character data to remaining processes

     call Broadcast(charVec, master_num, "charVec")

     ! scatter broadcasted data

     lastChar = 0
     do ng=1,ngrids
        do nf=1,nndvifiles(ng)
           do nc=1,lenFnames
              fnames_ndvi(nf,ng)(nc:nc) = charVec(lastChar+nc)
           end do
           lastChar = lastChar + lenFnames
           do nc = 1, lenItot
              itotdate_ndvi(nf,ng)(nc:nc)=charVec(lastChar+nc)
           end do
           lastChar = lastChar + lenItot
        end do
     end do

     ! deallocate broadcast area

     deallocate(charVec, stat=ierr2)
     if (ierr2 /= 0) then
        write(c0,"(i8)") ierr2
        call fatal_error(h//" deallocate charVec fails with stat="//trim(adjustl(c0)))
     end if
  end if
end subroutine NdviFileInv








subroutine NdviUpdate(iswap,nfile)

  use mem_leaf, only: &
       leaf_g

  use mem_grid, only: &
       ngrids, npatch, nnxp, nnyp

  use io_params, only: &
       fnames_ndvi

  use node_mod, only:      &
       mchnum,             &
       mynum,              &
       master_num

  use ReadBcst, only: &
       ReadStoreOwnChunk

  implicit none

  include "files.h"
  include "UseVfm.h"

  integer :: iswap,nfile
  integer :: ios
  integer :: ng,nc,ip
  character(len=1) :: cgrid
  character(len=f_name_length) :: flnm
  character(len=1) :: dummy
  real, pointer :: p2D(:,:)
  character(len=*), parameter :: h="**(NdviUpdate)**"

  ! Put new fields into future arrays. If iswap == 1, 
  !     swap future into past first

  if (iswap == 1) then
     do ng=1,ngrids
        do ip = 1,npatch
           leaf_g(ng)%veg_ndvip(:,:,ip) = leaf_g(ng)%veg_ndvif(:,:,ip)
        enddo
     enddo
  endif


  ! Open the input file for each grid and read field.
  do ng=1,ngrids

     ! master process builds filename, opens file and skips headers

     if (mchnum == master_num) then
        write(cgrid, '(i1)') ng
        flnm=fnames_ndvi(nfile,ng)
        nc=len_trim(flnm)-4
        flnm(nc:nc)=cgrid
        if (useVfm) then
           call rams_f_open(25,flnm(1:len_trim(flnm)),'FORMATTED','OLD','READ',0)
           read(25,*) dummy
           read(25,*) dummy
           read(25,*) dummy
           read(25,*) dummy
        else
           open(25, action="read", file=trim(flnm), form="unformatted", iostat=ios)
           if (ios /= 0) then
              call fatal_error(h//" error opening file "//trim(flnm))
           end if
           read(25) dummy
           read(25) dummy
           read(25) dummy
           read(25) dummy
        end if
        
     end if

     ! deals with ndivf

     do ip = 1,npatch
        p2D => leaf_g(ng)%veg_ndvif(:,:,ip)
        call ReadStoreOwnChunk(ng, 25, p2D, "ndvif")
     end do

     ! master process close the input file

     if (mchnum == master_num) then
        close(25)
     end if
  end do
end subroutine NdviUpdate

