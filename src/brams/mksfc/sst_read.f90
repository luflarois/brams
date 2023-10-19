!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine sst_check_header(ifm,flnm,ierr)

  use mem_grid, only: &
       platn, &
       plonn, &
       xtn, &
       ytn, &
       nnxp, &
       nnyp, &
       deltaxn, &
       deltayn

  ! This subroutine checks for the existence of an sst file for
  ! grid number ifm, and if it exists, also checks for agreement of
  ! grid configuration between the file and the current model run.
  ! If the file does not exist or does not match grid configuration,
  ! the flag ierr is returned with a value of 1.  If the file
  ! exists and is ok, ierr is returned with a value of 0.

  implicit none

  include "files.h"
  include "UseVfm.h"

  integer :: ifm,ierr
  character(len=*),intent(in) :: flnm

  integer :: iiyear,iimonth,iidate,iihour  &
       ,isfc_marker,isfc_ver,nsfx,nsfy
  integer :: ios
  real :: sfdx,sfdy,sfplat,sfplon,sflat,sflon,glatr,glonr
  character(len=*), parameter :: h="**(sst_check_header)**"
  
  ierr = 0

  call xy_ll(glatr,glonr,platn(ifm),plonn(ifm),xtn(1,ifm),ytn(1,ifm))

  if (useVfm) then

     call rams_f_open(25,flnm(1:len_trim(flnm)),'FORMATTED','OLD','READ',0)
     read (25,*) isfc_marker,isfc_ver
     read (25,*) iiyear,iimonth,iidate,iihour
     read (25,*) nsfx,nsfy
     read (25,*) sfdx,sfdy,sfplat,sfplon,sflat,sflon
     close (25)

  else

     open(25, action="read", file=trim(flnm), form="unformatted", iostat=ios)
     if (ios /= 0) then
        call fatal_error(h//" error opening file "//trim(flnm))
     end if
     read (25) isfc_marker,isfc_ver
     read (25) iiyear,iimonth,iidate,iihour
     read (25) nsfx,nsfy
     read (25) sfdx,sfdy,sfplat,sfplon,sflat,sflon
     close (25)
  end if

  if (nsfx                   .ne. nnxp(ifm) .or.  &
       nsfy                   .ne. nnyp(ifm) .or.  &
       abs(sfdx-deltaxn(ifm)) .gt. .001      .or.  &
       abs(sfdy-deltayn(ifm)) .gt. .001      .or.  &
       abs(sfplat-platn(ifm)) .gt. .001      .or.  &
       abs(sfplon-plonn(ifm)) .gt. .001      .or.  &
       abs(sflat-glatr)       .gt. .001      .or.  &
       abs(sflon-glonr)       .gt. .001) then

     ierr = 1

     print*,'SSTfile mismatch on grid:',ifm
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
     print*,'-------------------'

  else

     ierr = 0
  endif
end subroutine sst_check_header

!==========================================================

subroutine SstReadStoreOwnChunk(runflag, ifm, ierr)
  use ModDateUtils
  use mem_grid, only: &
       iyear1, imonth1, idate1, itime1, time, &
       runtype, ngrids, timmax, nnxp, nnyp

  use mem_leaf, only: &
       leaf_g

  use io_params, only: &
       itotdate_sst, &
       isstflp, &
       isstflf, &
       ssttime1, &
       ssttime2, &
       sstfpfx, &
       isstcycdata, &
       isstcyclic, &
       iupdsst, &
       nsstfiles

  use node_mod, only: &
       mchnum,        &
       master_num

  implicit none

  integer, intent(in)  :: runflag
  integer, intent(in)  :: ifm
  integer, intent(out) :: ierr

  character(len=14), save :: totdate_start, totdate_init
  character(len=14)  :: totdate, totdatem
  integer :: iyears, imonths, idates, ihours, nf, ng
  real(kind=8) :: secs_init, secs1, secs2
  character(len=*), parameter :: h="**(SstReadStoreOwnChunk)**"

  ierr = 0

  if (runflag==1 .or. runflag==2) then   ! Initialization(1) or file check(2)

     ! Inventory all sst surface files.

     call SstFileInv(sstfpfx(1:len_trim(sstfpfx)), ierr)
     if (ierr==1) then
        if (runflag==2) then
           return
        else if (runflag==1) then
           call fatal_error(h//' sst_read: error on init')
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
     endif

     ! Do some checks on times

     do ng=1,ngrids
        if (itotdate_sst(1,ng)>totdate_start .and. isstcycdata==0) then
           if (mchnum==master_num) then
              print *, h//' initial sst file time for grid ', ng, &
                   ' later than beginning of run'
           end if
           ierr = 1
           return
        endif

        if (iupdsst==1 .and. nsstfiles(ng)==1) then
           if (mchnum==master_num) then
              print*, h//' updating SST values but only one SST file',  &
                   ' for grid', ng
           end if
           ierr = 1
           return
        endif

        call date_add_to_big(totdate_init, timmax, 's', totdatem)
        if (iupdsst==1 .and. isstcycdata==0 .and. &
             itotdate_sst(nsstfiles(ng),ng)<totdatem) then
           if (mchnum==master_num) then
              print*, h//' final sst file time for grid ', ng,  &
                   'earlier than end of run - making new sst files'
           end if
           ierr = 1
           return
        endif
     enddo

     if (ierr==1) then
        if (runflag==1) then
           call fatal_error(h//' sst_read: error on init time checks')
        end if
     end if

     ! If we are only checking, we're done.
     if (runflag==2) return

     do ng=1,ngrids

        ! Change the sst file dates to the current year. We will increment
        !   when we need to read a future file.

        if (isstcycdata==1) then
           do nf=1,nsstfiles(ng)
              itotdate_sst(nf,ng)(1:4) = totdate_start(1:4)
           end do
        end if

        ! Find past time file. The files are ordered in time, so we only
        !    need to find when start time is greater than a file.

        isstflp(ng) = 0
        do nf=nsstfiles(ng),1,-1
           if (totdate_start>=itotdate_sst(nf,ng)) then
              isstflp(ng) = nf
              exit
           end if
        end do

        isstflf(ng) = isstflp(ng) + 1

        ! If we are cyclic, we possibly didn't find a start time later than
        !   a file time. I think it is safe to assume that we will use the
        !   last file for the start file. Also see if future time will be on
        !   next cycle.

        if (isstcycdata==1) then
           if (isstflp(ng)==0)            isstflp(ng) = nsstfiles(ng)
           if (isstflf(ng)>nsstfiles(ng)) isstflf(ng) = 1
        endif

        ! Read past time sst field

        call SstUpdate(0, isstflp(ng))

        if (iupdsst==1) then

           ! Read future time sst field if updating

           call SstUpdate(1, isstflf(ng))

           ! Compute times as number of seconds past 1 Jan 1900

           call date_abs_secs(totdate_init, secs_init)

           ! Get model time of past file

           totdatem = itotdate_sst(isstflp(ng),ng)
           if (isstcyclic==1) then

              ! If month of past file > current month, subtract a year

              if (totdatem(5:6)>totdate_start(5:6)) &
                   call date_add_to_big(totdatem, -365., 'd', totdatem)
           endif
           call date_abs_secs(totdatem, secs1)

           totdatem = itotdate_sst(isstflf(ng),ng)
           if (isstcyclic==1) then
              if (totdatem<totdate_start) then

                 ! Future file is in next year. Update all file names for a new year

                 call date_add_to_big(totdatem, 365., 'd', totdatem)
                 do nf=1,nsstfiles(ng)
                    itotdate_sst(nf,ng)(1:4) = totdatem(1:4)
                 enddo
              endif
           endif
           call date_abs_secs(totdatem, secs2)

           ssttime1(ng) = secs1 - secs_init
           ssttime2(ng) = secs2 - secs_init
        else
           leaf_g(ng)%seatp(:,:) = leaf_g(ng)%seatf(:,:)
           ssttime1(ng)          = 0.
           ssttime2(ng)          = 0.
        endif

     enddo

     return

  elseif (runflag==3) then   ! Runtime file increment

     if (time>=ssttime2(ifm)) then

        ! Update sst fields
        isstflp(ifm) = isstflf(ifm)
        isstflf(ifm) = isstflp(ifm) + 1

        ! Compute times as number of seconds past 1 Jan 1900
        !   If cyclic, modify file date
        call date_abs_secs(totdate_init, secs_init)
        call date_add_to_big(totdate_init,time, 's', totdate)

        totdatem = itotdate_sst(isstflp(ifm),ifm)
        call date_abs_secs(totdatem, secs1)

        ! Need to deal with cyclic.
        if (isstcyclic==1 .and. isstflf(ifm)>nsstfiles(ifm)) then
           isstflf(ifm) = 1
           ! Update all file names for a new year
           totdatem     = itotdate_sst(isstflf(ifm),ifm)
           call date_add_to_big(totdatem, 365., 'd', totdatem)
           do nf=1,nsstfiles(ifm)
              itotdate_sst(nf,ifm)(1:4) = totdatem(1:4)
           enddo
        endif

        totdatem = itotdate_sst(isstflf(ifm),ifm)
        call date_abs_secs(totdatem, secs2)

        ssttime1(ifm) = secs1 - secs_init
        ssttime2(ifm) = secs2 - secs_init

        ! Finally read the actual field

        call SstUpdate(1, isstflf(ifm))

        if (mchnum==master_num) then
           print*, h//' Switched sst files:',ssttime1(ifm),ssttime2(ifm), &
                secs1, secs2, secs_init
        end if

     else

        return

     endif

  endif
end subroutine SstReadStoreOwnChunk








subroutine SstFileInv(sfilin, ierr)
  use ModDateUtils
  use mem_grid, only: &
       iyear1, &
       imonth1, &
       idate1, &
       itime1, &
       ngrids, &
       maxsstfiles, & 
       nxtnest

  use mem_mksfc, only: &
       nvsstf,          &
       idatevs,         &
       ihourvs,         &
       imonthvs,        &
       iyearvs       

  use io_params, only: &
       isstcyclic, &     ! intent(out)
       isstcycdata, &    ! intent(out) 
       nsstfiles, &      ! intent(out)
       fnames_sst, &     ! intent(out)
       itotdate_sst, &   ! intent(out)
       iupdsst,&         ! intent(in)
       isstflg


  use node_mod, only:      &
       mchnum,             &
       master_num

  use ReadBcst, only: &
       Broadcast

  use grid_dims, only: maxfiles


  implicit none

  include "files.h"
  include "UseVfm.h"

  character(len=*), intent(IN) :: sfilin
  integer, intent(OUT)           :: ierr

  integer :: nc, nf, lnf, nftot, ng, icm, ifm
  integer :: inyear, inmonth, indate, inhour

!!$  integer, parameter :: maxfiles=1000
  character(len=f_name_length) :: fnames(maxfiles)
  character(len=f_name_length) :: rams_filelist_arg
  character(len=014) :: itotdate
  character(len=001) :: cgrid

  logical :: there

  integer :: lastChar
  integer :: lenFnames
  integer :: lenItot
  integer :: sizeIntVec
  integer :: sizeCharVec
  integer :: ierr2
  integer, allocatable :: intVec(:)
  character, allocatable :: charVec(:)

  real(kind=8) :: secs_init, secs_file
  character(len=8) :: c0, c1
  integer :: nvtime,ivtime,indice
  character(len=2) :: cgrid2
  character(len=f_name_length) :: flnm
  character(len=*), parameter :: h="**(SstFileInv)**"

  ! size and length of data to be broadcasted

  sizeIntVec = ngrids+3
  lenFnames  = len(fnames_sst)
  lenItot    = len(itotdate_sst)

  ! only one process does file inventory

  if (mchnum==master_num) then

     ! io_params initialization

     nsstfiles   = 0
     !     fnames_sst= ""
     !     itotdate_sst= ""
     isstcyclic  = 0
     isstcycdata = 0
     ierr        = 0

     ! Get abs seconds of run start

     call date_abs_secs2(iyear1, imonth1, idate1, itime1*100, secs_init)

     ! Go through sst files and make inventory. We unfortunately have to do this
     !   for all grids.


     do ng=1,ngrids

        write(cgrid2,'(a1,i1)') 'g',ng

        nftot = -1
        write(cgrid,'(i1)') ng
        rams_filelist_arg = trim(sfilin)//'-W-*-g'//cgrid//'.vfm'

!!$        print *, "DEBUG-ALF:SstFileInv:trim(sfilin)=", trim(sfilin)
!!$        print *, "DEBUG-ALF:SstFileInv:cgrid=", cgrid
!!$        print *, "DEBUG-ALF:SstFileInv:argumento=", rams_filelist_arg
!!$        call flush(6)


        if (isstflg(ng) == 1) then
           call sst_read_dataheader(ng)
           nvtime = nvsstf(ng)
        elseif (isstflg(ng) == 0) then
           nvtime = nvsstf(nxtnest(ng))
        else
           nvtime = 1
        end if

        indice = 1

        ifm = ng
        icm = nxtnest(ifm)

        do ivtime = 1,nvtime
           !!call sstnest (ng,ivtime)
           if (icm>=1 .and. isstflg(ifm)==0) then
              nvsstf(ifm)                 = nvsstf(icm)
              iyearvs (1:nvsstf(ifm),ifm) = iyearvs (1:nvsstf(ifm),icm)
              imonthvs(1:nvsstf(ifm),ifm) = imonthvs(1:nvsstf(ifm),icm)
              idatevs (1:nvsstf(ifm),ifm) = idatevs (1:nvsstf(ifm),icm)
              ihourvs (1:nvsstf(ifm),ifm) = ihourvs (1:nvsstf(ifm),icm)
           endif

           if (useVfm) then
              call makefnam(flnm, sfilin, 0., &
                   iyearvs(ivtime,ng), imonthvs(ivtime,ng), &
                   idatevs(ivtime,ng), ihourvs(ivtime,ng)*10000, &
                   'W',cgrid2,'vfm')
           else
              call makefnam(flnm, sfilin, 0., &
                   iyearvs(ivtime,ng), imonthvs(ivtime,ng), &
                   idatevs(ivtime,ng), ihourvs(ivtime,ng)*10000, &
                   'W',cgrid2,'bin')
           end if


           inquire(file=flnm(1:len_trim(flnm)), exist=there)

           if (there) then
              fnames(indice) = trim(flnm)
              indice = indice + 1 
           endif

        end do

        nftot = indice - 1


!!$     call RAMS_filelist(fnames, rams_filelist_arg, nftot)

        if (nftot<=0) then
           print *, 'No sst files for grid '//cgrid
           ierr = 1
           exit
        end if

        if (nftot>maxsstfiles) then
           call fatal_error(h//' too many sst files')
        end if

        ! We need to see if the data is to be considered "cyclic". Count how many files
        !  have the year 0000. If all of them do not, we probably have an error.
        !  isstcycdata = 1 when data is cyclic
        !  isstcyclic =1 when we are have cyclic data and will be updating in time

        nsstfiles(ng) = 0

        do nf=1,nftot
           lnf = len_trim(fnames(nf))
           read(fnames(nf) (lnf-23:lnf-7), "(i4,1x,i2,1x,i2,1x,i6)") &
                inyear, inmonth, indate, inhour

           ! Check the file headers

           call sst_check_header(ng, fnames(nf)(1:len_trim(fnames(nf))), ierr)
           if (ierr/=0) exit

           if (inyear==0) isstcycdata = isstcycdata + 1

           call date_make_big(inyear, inmonth, indate, inhour, itotdate)

           nsstfiles(ng)                  = nsstfiles(ng) + 1
           fnames_sst(nsstfiles(ng),ng)   = fnames(nf)
           itotdate_sst(nsstfiles(ng),ng) = itotdate
        end do

        if (ierr==0) then
           call RAMS_dintsort(nsstfiles(ng), itotdate_sst(1,ng), &
                fnames_sst(1,ng))

           !  start printing section

           open(unit=22,file='brams.log',position='append',action='write')
           write(unit=22,fmt='(A)') ' '
           write(unit=22,fmt='(A)') ' '
           write(unit=22,fmt='(A)') ' '
           write(unit=22,fmt='(A)') '-------------------------------------------------------------'
           write(unit=22,fmt='(A)') '-----------  SST Input File Inventory: Grid '//cgrid
           write(unit=22,fmt='(A)') '-------------------------------------------------------------'
           do nf=1,nsstfiles(ng)
              write(unit=22,fmt='(A,A,A)') itotdate_sst(nf,ng), '   ', trim(fnames_sst(nf,ng))
           enddo
           write(unit=22,fmt='(A)') '------------------------------------------------------'
           close(unit=22)
        end if
        if (ierr/=0) exit
     end do

     ! Check the cyclic data condition.
     ! WE ARE ONLY ALLOWING CYCLIC ON ALL GRIDS

     if (ierr==0) then
        if (isstcycdata>0) then
           if (isstcycdata/=sum(nsstfiles(1:ngrids))) then
              call fatal_error(h//&
                   ' All sst surface files do not have year 0000;'//&
                   ' This confuses the gods and can not occur.')
           endif
           isstcycdata = 1
        else
           isstcycdata = 0
        endif

        ! Set the main cyclic flag. Only relevant if we are updating in time.

        if (iupdsst==1 .and. isstcycdata==1) then
           isstcyclic = 1
        end if
     end if
  end if

  ! allocate broadcast area

  allocate(intVec(sizeIntVec), stat=ierr2)
  if (ierr2/=0) then
     write(c0,"(i8)") ierr2
     write(c1,"(i8)") sizeIntVec
     call fatal_error(h//" allocate intVec("//trim(adjustl(c1))//&
          ") fails with stat="//trim(adjustl(c0)))
  end if

  ! master process gathers data for broadcasting

  if (mchnum==master_num) then
     intVec(1:3)          = (/ierr, isstcyclic, isstcycdata /)
     intVec(4:sizeIntVec) = nsstfiles(1:ngrids)
  end if

  ! broadcast integer data to remaining processes

  call Broadcast(intVec, master_num, "intVec")

  ! scatter broadcasted data

  ierr        = intVec(1)
  isstcyclic  = intVec(2)
  isstcycdata = intVec(3)

  nsstfiles(:) = 0
  nsstfiles(1:ngrids) = intVec(4:sizeIntVec)
!!$  nsstfiles(ngrids+1:)=0

  ! deallocate broadcast area

  deallocate(intVec, stat=ierr2)
  if (ierr2/=0) then
     write(c0,"(i8)") ierr2
     call fatal_error(h//" deallocate intVec fails with stat="//&
          trim(adjustl(c0)))
  end if

  ! broadcast remaining data only if inventory succeeds

  if (ierr==0) then

     ! allocate broadcast area

     sizeCharVec = sum(nsstfiles(1:ngrids))*(lenFnames + lenItot)
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
           do nf=1,nsstfiles(ng)
              do nc=1,lenFnames
                 charVec(lastChar+nc) = fnames_sst(nf,ng)(nc:nc)
              end do
              lastChar = lastChar + lenFnames
              do nc=1,lenItot
                 charVec(lastChar+nc) = itotdate_sst(nf,ng)(nc:nc)
              end do
              lastChar = lastChar + lenItot
           end do
        end do
     end if

     ! broadcast character data to remaining processes

     call Broadcast(charVec, master_num, "charVec")

     ! scatter broadcasted data

     lastChar=0
     do ng=1,ngrids
        do nf=1,nsstfiles(ng)
           do nc=1,lenFnames
              fnames_sst(nf,ng)(nc:nc) = charVec(lastChar+nc)
           end do
           lastChar = lastChar + lenFnames
           do nc=1,lenItot
              itotdate_sst(nf,ng)(nc:nc) = charVec(lastChar+nc)
           end do
           lastChar = lastChar + lenItot
        end do
     end do

     ! deallocate broadcast area

     deallocate(charVec, stat=ierr2)
     if (ierr2/=0) then
        write(c0,"(i8)") ierr2
        call fatal_error(h//" deallocate charVec fails with stat="//&
             trim(adjustl(c0)))
     end if
  end if
end subroutine SstFileInv







subroutine SstUpdate(iswap, nfile)

  use mem_grid, only: &
       ngrids

  use mem_leaf, only: &
       leaf_g

  use io_params, only: &
       fnames_sst

  use node_mod, only:      &
       mchnum,             &
       mynum,              &
       master_num

  use ReadBcst, only: &
       ReadStoreOwnChunk

  implicit none

  include "files.h"
  include "UseVfm.h"

  integer, intent(IN) :: iswap, nfile

  integer            :: ng, nc
  integer :: ios
  character(len=001) :: cgrid
  character(len=f_name_length) :: flnm
  character(len=001) :: dummy
  character(len=*), parameter :: h="**(SstUpdate)**"

  ! Put new fields into future arrays. If iswap == 1,
  !     swap future into past first

  if (iswap==1) then
     do ng=1,ngrids
        leaf_g(ng)%seatp(:,:) = leaf_g(ng)%seatf(:,:)
     enddo
  endif


  ! Open the input file for each grid and read field.

  do ng=1,ngrids

     ! master process builds filename, opens file and skips headers

     if (mchnum==master_num) then
        write(cgrid, '(i1)') ng
        flnm        = trim(fnames_sst(nfile,ng))
        nc          = len_trim(flnm) - 4
        flnm(nc:nc) = cgrid

        if (useVfm) then
           call rams_f_open(25, flnm(1:len_trim(flnm)), 'FORMATTED', 'OLD', 'READ', 0)

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

     ! deals with seatf

     call ReadStoreOwnChunk(ng, 25, leaf_g(ng)%seatf, "seatf")

     ! master process closes file

     if (mchnum==master_num) then
        close(25)
     end if
  enddo

end subroutine SstUpdate
