module ModVarfFile

   use ModDateUtils

   use mem_scratch, only: &
      vctr2

   use ParLib, only: &
      parf_bcast, &
      parf_barrier, &
      parf_allreduce_sum

   use rconstants, only: &
      g, &
      cp, &
      cpor, &
      p00, &
      rgas

   use ref_sounding, only: &
      topref, &
      iref, &
      jref, &
      u01dn, &
      v01dn, &
      rt01dn, &
      th01dn, &
      pi01dn, &
      dn01dn

   use mem_varinit, only: &
      vwait1, &
      vwaittot, &
      varfpfx, &
      lastdate_iv, &
      vtime1, &
      vtime2, &
      fnames_varf, &
      maxnudfiles, &
      nvarffiles, &
      itotdate_varf, &
      varf_times, &
      tnudcent, &
      tnudtop, &
      nvarffl, &
      varinit_g, &
      wt_nudge_uv, &
      wt_nudge_th, &
      wt_nudge_pi, &
      wt_nudge_rt
   use node_mod, only: &
      i0, &
      j0, &
      mchnum, &
      master_num, &
      nmachs, &
      mxp, &
      myp, &
      mzp, &
      nodei0, &
      nodej0, &
      nodemxp, &
      nodemyp, &
      mynum

   use ReadBcst, only: &
      gatherData, &
      ProcWithMin, &
      ReadStoreOwnChunk, &
      Broadcast, &
      storeOwnChunk_3D, &
      DumpFullField

   use isan_coms, only: &
      isan_inc

   use mem_grid, only: &
      iyear1, &
      imonth1, &
      itime1, &
      idate1, &
      time, &
      ngrids, &
      nxtnest, &
      nnzp, &
      nnxp, &
      nnyp, &
      grid_g, &
      timmax, &
      ngrid, &
      nxp, &
      nyp, &
      nzp, &
      nnxyzp, &
      nxyp, &
      nxyzp, &
      deltax, &
      deltay, &
      deltaz, &
      dzrat, &
      dzmax, &
      if_adap, &
      nzg, &
      nzs, &
      npatch, &
      platn, &
      plonn, &
      ztop, &
      ztn, &
      zt, &
      dzm

   use ModGridTree, only: &
      GridTree, &
      FetchGrid

   use ModGrid, only: &
      Grid

   use mem_chem1, only: &
      chem1_g, &
      chem_assim, &
      chemistry

   use mem_leaf, only: &
      leaf_g

   use mem_basic, only: &
      basic_g

   use micphys, only: &
      level

   use chem1_list, only: &
      nspecies, &
      spc_alloc, &
      fdda, &
      chemical_mechanism, &
      init_ajust, &
      spc_name

   use ModMessageSet, only: &
      PostRecvSendMsgs, &
      WaitRecvMsgs

   implicit none

   include "files.h"
   include "i8.h"

   private
   public :: VarfReadStoreOwnChunk

contains

   subroutine VarfReadStoreOwnChunk(AllGrids, ivflag, initialFlag)

      type(GridTree), pointer :: AllGrids
      integer, intent(IN) :: ivflag
      integer, intent(in) :: initialFlag

      character(len=14)  :: itotdate_current
      integer :: iyears, imonths, idates, ihours, nf, ifm, &
                 ivar_wait, nwaits, nw, ivwait1, irsleep, irslp, ifileok, icm
      type(Grid), pointer :: OneGrid => null()
      character(len=8) :: c0
      character(len=*), parameter :: h = "**(VarfReadStoreOwnChunk)**"
      !srf - special weights for pressure (only for operations)
      real, dimension(mzp, mxp, myp) :: varwts_for_operations_only

      ! See if we want to possibly wait for files to be available.
      !    This will control some logic...

      ivar_wait = 0
      if (vwait1 > 0.0 .and. vwaittot > 0.0) ivar_wait = 1
      if (ivar_wait == 0) then
         nwaits = 1
      elseif (ivar_wait == 1) then
         if (mchnum == master_num) then
            print *, ' Will wait for varfiles if not present'
            print *, '  Check interval:', vwait1, ' Fail time:', vwaittot
         end if
         nwaits = int(vwaittot/vwait1) + 1
      end if

      if (ivflag == 0) then ! Initialization of initial fields

         call date_make_big(iyear1, imonth1, idate1, itime1*100, itotdate_current)

         wait: do nw = 1, nwaits

            ! Inventory all varf files.

            call VarfFileInv(varfpfx, iyear1, imonth1, idate1, itime1, initialFlag)

            ! The initial time must have an exact time match.
            nvarffl = 0
            do nf = 1, nvarffiles

               if (itotdate_current == itotdate_varf(nf)) then
                  nvarffl = nf
                  exit wait
               end if
            end do
            if (ivar_wait == 0) then
               call fatal_error(h//' No initial varfiles found with prefix '// &
                                trim(varfpfx))
            elseif (ivar_wait == 1 .and. nw < nwaits) then
               if (mchnum == master_num) then
                  print *, 'No initial varfiles found: ', trim(varfpfx)
                  print *, '    Waiting:', vwait1, ' seconds.   Total wait:', &
                     nw*vwait1
               end if
               ivwait1 = nint(vwait1)
               irslp = irsleep(ivwait1)
            elseif (ivar_wait == 1 .and. nw == nwaits) then
               call fatal_error(h// &
                                " Waited too long; no initial varfiles found with prefix "// &
                                trim(varfpfx))
            end if

         end do wait

         ! Now do actual initialization for the coarse grid
         !     and find 1D reference state
         call newgrid(1)

         OneGrid => FetchGrid(AllGrids, 1)

         call VarfUpdate(0, OneGrid, ifileok, 1, initialFlag)

         !  On all fine grids, initialize the 1-D reference state arrays,
         !  the 3-D reference state arrays,
         !  and the prognostic atmospheric fields by interpolation.

         call FineMeshRefSound1D(2, ngrids)

         do ifm = 2, ngrids

            !**(JP)** not converted

            call fatal_error(h//"**(JP)** not converted for multiple grids")

            icm = nxtnest(ifm)

            ! Get 3D reference state for this grid
            call fmrefs3d(ifm)

            ! Interpolate prognostic fields. These will be overwritten if the varfile
            !    exists
            call prgintrp(nnzp(icm), nnxp(icm), nnyp(icm), &
                          nnzp(icm), nnxp(icm), nnyp(icm), 0, 0, ifm, 1, 0)

            call newgrid(ifm)

            ! See if this grid's varfile is created.

            OneGrid => FetchGrid(AllGrids, ifm)

            call VarfUpdate(0, OneGrid, ifileok, 1, initialFlag)

            open (unit=22, file='brams.log', position='append', action='write')
            if (ifileok == 1) then
               ! Everything's cool...
               if (mchnum == master_num) then
                  write (c0, "(i8)") ifm
                  write (unit=22, fmt="(a)") &
                     " Initial varfile read of grid-"//trim(adjustl(c0))
               end if
            else
               ! Using interpolated nudging arrays from parent grid.
               if (mchnum == master_num) then
                  write (c0, "(i8)") ifm
                  write (unit=22, fmt="(a)") &
                     " Initial interpolation of grid-"//trim(adjustl(c0))
               end if
            end if
            close (unit=22)

            call fmdn0(ifm)

         end do

         ! ALF
         lastdate_iv = itotdate_varf(nvarffl)

         return

      elseif (ivflag == 1) then ! Fill nudging arrays and compute weights

         ! If a history start, we will assume a past time file is there.
         !  Take closest past time.

         call date_add_to(iyear1, imonth1, idate1, itime1*100, &
                          time, 's', iyears, imonths, idates, ihours)
         call date_make_big(iyears, imonths, idates, ihours, itotdate_current)

         call VarfFileInv(varfpfx, iyear1, imonth1, idate1, itime1, initialFlag)

         nvarffl = 0
         do nf = nvarffiles, 1, -1
            if (itotdate_varf(nf) <= itotdate_current) then
               nvarffl = nf
               exit
            end if
         end do
         if (nvarffl == 0) then
            call fatal_error(h// &
                             "No past varfiles found on nudge fill:"//trim(varfpfx))
         end if

         ! Compute weighting factors for grid 1
         call VariableWeight(nnzp(1), nodemxp(mynum, 1), nodemyp(mynum, 1), nnxp(1), &
                             nnyp(1), nodei0(mynum, 1), nodej0(mynum, 1), &
                             grid_g(1)%topt(1, 1), grid_g(1)%rtgt(1, 1), varinit_g(1)%varwts(1, 1, 1), &
                             !srf
                             varwts_for_operations_only)

         if (chem_assim == 1 .and. chemistry >= 0) &
            call VariableWeightChem(nnzp(1), nodemxp(mynum, 1), nodemyp(mynum, 1), nnxp(1), &
                                    nnyp(1), nodei0(mynum, 1), nodej0(mynum, 1), &
                                    grid_g(1)%topt(1, 1), grid_g(1)%rtgt(1, 1), &
                                    varinit_g(1)%varwts_chem(1, 1, 1))

         ! Read files

         do ifm = 1, ngrids
            icm = nxtnest(ifm)

            call newgrid(ifm)

            ! Interpolate weights to all other grids
            if (ifm > 1) then
               call VarfIntrp(ifm, 1)
            end if

            ! See if this grid's varfile is created.

            OneGrid => FetchGrid(AllGrids, ifm)

            call VarfUpdate(0, OneGrid, ifileok, 0, initialFlag)
            open (unit=22, file='brams.log', position='append', action='write')
            if (ifileok == 1) then
               ! Everything's cool...
               if (mchnum == master_num) then
                  write (unit=22, fmt='(A,1X,I2.2)') 'Varfile read of grid-', ifm
               end if
            else
               ! Using interpolated nudging arrays from parent grid.
               call VarfIntrp(ifm, 2)
               if (mchnum == master_num) then
                  write (unit=22, fmt='(A,1X,I2.2)') 'Interpolation of grid-', ifm
               end if
            end if

         end do
         vtime2 = varf_times(nvarffl)
         if (mchnum == master_num) then
            write (unit=22, fmt='(A,1X,I2.2,1X,F12.2)') 'New varfile times:', nvarffl, vtime2
         end if
         close (unit=22)

         ! ALF
         lastdate_iv = itotdate_varf(nvarffl)

      elseif (ivflag == 2) then ! Runtime file increment

         ! Find current date/time

         call date_add_to(iyear1, imonth1, idate1, itime1*100, &
                          time, 's', iyears, imonths, idates, ihours)
         call date_make_big(iyears, imonths, idates, ihours, itotdate_current)

      end if

      ! Find the next varfile in the list, waiting for it if necessary

      wait2: do nw = 1, nwaits

         ! Redo the inventory in case new files showed up

         call VarfFileInv(varfpfx, iyear1, imonth1, idate1, itime1, initialFlag)

         nvarffl = 0
         do nf = 1, nvarffiles
            if (itotdate_varf(nf) > itotdate_current) then
               nvarffl = nf
               exit wait2
            end if
         end do
         if (ivar_wait == 0) then
            call fatal_error(h//' No future varfiles found with prefix '// &
                             trim(varfpfx))
         elseif (ivar_wait == 1 .and. nw < nwaits) then
            if (mchnum == master_num) then
               print *, 'No initial varfiles found: ', trim(varfpfx)
               print *, '    Waiting:', vwaittot, ' seconds'
            end if
            ivwait1 = nint(vwait1)
            irslp = irsleep(ivwait1)
         elseif (ivar_wait == 1 .and. nw == nwaits) then
            call fatal_error(h// &
                             " Waited too long; no future varfiles found with prefix "// &
                             trim(varfpfx))
         end if

      end do wait2

      ! Read future files

      do ifm = 1, ngrids
         icm = nxtnest(ifm)

         call newgrid(ifm)

         ! See if this grid's varfile is created.

         OneGrid => FetchGrid(AllGrids, ifm)

         call VarfUpdate(1, OneGrid, ifileok, 0, initialFlag)
         open (unit=22, file='brams.log', position='append', action='write')
         if (ifileok == 1) then
            ! Everything's cool...
            if (mchnum == master_num) then
               write (unit=22, fmt='(A,1X,I2.2)') 'Future varfile read of grid-', ifm
            end if
         else
            call VarfIntrp(ifm, 2)
            if (mchnum == master_num) then
               write (unit=22, fmt='(A,1X,I2.2)') 'Future interpolation of grid-', ifm
            end if
         end if
         close (unit=22)

      end do

      vtime1 = vtime2
      vtime2 = varf_times(nvarffl)

      if (mchnum == master_num) then
         open (unit=22, file='brams.log', position='append', action='write')
         write (unit=22, fmt='(A,1X,I2.2,2(F12.2,1X))') 'New varfile times:', nvarffl, vtime1, vtime2
         close (unit=22)
      end if

      ! ALF
      lastdate_iv = itotdate_varf(nvarffl)
   end subroutine VarfReadStoreOwnChunk

   subroutine VarfFileInv(varpref, iyear1, imonth1, idate1, itime1, flag)
      character(len=f_name_length), intent(in) :: varpref
      integer, intent(in) :: iyear1
      integer, intent(in) :: imonth1
      integer, intent(in) :: idate1
      integer, intent(in) :: itime1
      integer, intent(in) :: flag

      character(len=14)  :: itotdate_current
      integer :: iyears, imonths, idates, ihours
      integer :: nc, nf, lnf
      integer :: inyear, inmonth, indate, inhour
      character(len=f_name_length) :: sVarName

      integer :: localTime

      character(len=f_name_length), dimension(maxnudfiles) :: fnames
!!$  character(len=f_name_length) :: vpref
      character(len=f_name_length)::rams_filelist_arg
      character(len=14)  :: itotdate
      real(kind=8) :: secs_init, secs_varf

      integer :: lastChar
      integer :: lenFnames
      integer :: sizeCharVec
      integer :: ierr2
      character, allocatable :: charVec(:)

      integer :: indice

      character(len=8) :: c0, c1
      character(len=*), parameter :: h = "**(VarfFileInv)**"
      logical there
      real :: time_inc_sec

      ! Get abs seconds of run start

      nvarffiles = -1
      call date_abs_secs2(iyear1, imonth1, idate1, itime1*100, secs_init)

      ! one process goes through history files and make inventory

      if (mchnum == master_num) then
!!$     vpref=varpref

         nvarffiles = ceiling(((timmax/3600)/(isan_inc/100)) + 1)

!print*,"A===|" ,iyear1, imonth1, idate1, itime1*100,time

         call date_add_to(iyear1, imonth1, idate1, itime1*100, &
                          time, 's', iyears, imonths, idates, ihours)

         localTime = itime1

!print*,"B===|" , nvarffiles,iyears, imonths, idates, ihours,time

         call makefnam(sVarName, varpref, 0, iyears, imonths, idates, ihours, 'V', '$', 'tag')
         inquire (file=sVarName(1:len_trim(sVarName)), exist=there)
         indice = 1
         if (there) then
            fnames(indice) = trim(sVarName)
            indice = indice + 1
         else
            !     print*,"IVAR file not found:",sVarName(1:len_trim(sVarName))
         end if

         !print*,"================================================================"
         !print*,"Making Varfile Input Inventory for the remaining time simulation"

         do nf = 1, nvarffiles

!print*,"===|A" , nf,iyears, imonths, idates, ihours, localTime*100

            time_inc_sec = isan_inc/100*3600.

!          call date_add_to(iyears, imonths, idates, localTime*100,  &

!print*,"===|B" , nf,iyears, imonths, idates, ihours

            !--(DMK-CCATT-INI)-----------------------------------------------------
            if (flag .eq. 2) then

               call makefnam(sVarName, varpref, 0, iyears, imonths, idates, ihours, 'V', '$', 'tag')
               !print*,"2 sVarName",trim(sVarName), iyears, imonths, idates, ihours
               call date_add_to(iyears, imonths, idates, ihours, &
                                time_inc_sec, 's', iyears, imonths, idates, ihours)
            else if (flag .eq. 4) then

               call makefnam(sVarName, varpref, 0, iyears, imonths, idates, ihours, 'A', 'head', 'txt')
               call date_add_to(iyears, imonths, idates, ihours, &
                                time_inc_sec, 's', iyears, imonths, idates, ihours)
            end if

            !
            inquire (file=sVarName(1:len_trim(sVarName)), exist=there)

            if (there) then
               fnames(indice) = trim(sVarName)
               indice = indice + 1
            else
               !print*,"IVAR file not found:",sVarName(1:len_trim(sVarName))
            end if

            !if (localTime<=2100) then
!         ! if (localTime<=1800) then
            !   localTime = localTime + isan_inc
            !else
            !   localTime = 000000
            !   localTime = localTime + isan_inc
            !endif
         end do

         nvarffiles = indice - 1

!!$     rams_filelist_arg = trim(varpref)//'*.tag'

!!$     call RAMS_filelist(fnames,rams_filelist_arg , nvarffiles)

         if (nvarffiles > maxnudfiles) then
            call fatal_error(h//' too many varf files')
         end if
      end if

      ! broadcast number of files

      call Broadcast(nvarffiles, master_num, "nvftot")

      ! allocate broadcast area

      lenFnames = len(fnames)
      sizeCharVec = nvarffiles*lenFnames
      allocate (charVec(sizeCharVec), stat=ierr2)
      if (ierr2 /= 0) then
         write (c0, "(i8)") ierr2
         write (c1, "(i8)") sizeCharVec
         call fatal_error(h//" allocate charVec("//trim(adjustl(c1))// &
                          ") failed with stat="//trim(adjustl(c0)))
      end if

      ! master process prepares broadcast data

      if (mchnum == master_num) then
         lastChar = 0
         do nf = 1, nvarffiles
            do nc = 1, lenFnames
               charVec(lastChar + nc) = fnames(nf) (nc:nc)
            end do
            lastChar = lastChar + lenFnames
         end do
      end if

      ! broadcast character data to remaining processes

      call Broadcast(charVec, master_num, "charVec")

      ! scatter broadcasted data

      lastChar = 0

      do nf = 1, nvarffiles
         do nc = 1, lenFnames
            fnames(nf) (nc:nc) = charVec(lastChar + nc)
         end do
         lastChar = lastChar + lenFnames
      end do

      ! deallocate broadcast area

      deallocate (charVec, stat=ierr2)
      if (ierr2 /= 0) then
         write (c0, "(i8)") ierr2
         call fatal_error(h//" deallocate charVec fails with stat="// &
                          trim(adjustl(c0)))
      end if

      !

      do nf = 1, nvarffiles
         lnf = len_trim(fnames(nf))

         !--(DMK-CCATT-INI)-----------------------------------------------------
         if (flag .eq. 2) then
            read (fnames(nf) (lnf - 20:lnf - 4), "(i4,1x,i2,1x,i2,1x,i6)") inyear, inmonth, indate, inhour
         else if (flag .eq. 4) then
            read (fnames(nf) (lnf - 25:lnf - 9), "(i4,1x,i2,1x,i2,1x,i6)") inyear, inmonth, indate, inhour
         end if
         !--(DMK-CCATT-OLD)-----------------------------------------------------
         !     read(fnames(nf)(lnf-20:lnf-4), "(i4,1x,i2,1x,i2,1x,i6)") &
         !          inyear, inmonth, indate, inhour
         !--(DMK-CCATT-FIM)-----------------------------------------------------

         call date_make_big(inyear, inmonth, indate, inhour, itotdate)

         fnames_varf(nf) = fnames(nf)

         itotdate_varf(nf) = itotdate

         call date_abs_secs2(inyear, inmonth, indate, inhour, secs_varf)
         varf_times(nf) = secs_varf - secs_init

      end do

      call RAMS_dintsort(nvarffiles, itotdate_varf, fnames_varf)

      !  start printing section

      if (mchnum == master_num) then
         open (unit=22, file='brams.log', position='append', action='write')
         write (unit=22, fmt='(A)') ' '
         write (unit=22, fmt='(A)') ' '
         write (unit=22, fmt='(A)') ' '
         write (unit=22, fmt='(A)') '-------------------------------------------------------------'
         write (unit=22, fmt='(A)') '-----------       Varfile Input Inventory       -------------'
         write (unit=22, fmt='(A)') '-------------------------------------------------------------'
         do nf = 1, nvarffiles
            write (unit=22, fmt='(i4,1x,a16,1x,f10.0,2x,a)') nf, itotdate_varf(nf) &
               , varf_times(nf), trim(fnames_varf(nf))
         end do
         write (unit=22, fmt='(A)') '------------------------------------------------------'
         close (unit=22)
      end if

   end subroutine VarfFileInv

   subroutine VarfUpdate(iswap, OneGrid, ifileok, initflag, initialFlag)

      use chem1_list, only: spc_name, spc_alloc, on, fdda, nspecies
      use mem_chem1, only: CHEM_ASSIM, CHEMISTRY
      use mem_aer1, only: AEROSOL, aer1_g, AER_ASSIM
      use aer1_list, only: aer_name => spc_name, &
                           aer_mod_spc => aer_name, &
                           aer_nspecies => nspecies, &
                           aer_alloc => spc_alloc, &
                           aer_fdda => fdda, &
                           aer_on => on, &
                           aer_nmodes => nmodes

      integer, intent(in) :: iswap
      type(Grid), pointer :: OneGrid
      integer, intent(out) :: ifileok
      integer, intent(in) :: initflag
      integer, intent(in) :: initialFlag

      real :: scratch(nxyp), scratch2(nxyzp)

      !---------------------------------------------------------------+
      !    "Variable initialization"  initialization routines
      !---------------------------------------------------------------+
      logical there
      integer :: iver_var, nc, iyearx, imonthx, idatex, ihourx &
                 , nxpx, nypx, nzpx, imarker, iun
      real :: rlatx, wlon1x, deltaxx, deltayx, deltazx, dzratx, dzmaxx
      character(len=7)   :: cgrid

      !--(DMK-CCATT-INI)-----------------------------------------------------
      integer :: nspc, k
      character(len=32) :: chemical_mechanism_test
      !--(DMK-CCATT-END)-----------------------------------------------------

      character(len=f_name_length) :: flnm
      character(len=*), parameter :: h = "**(VarfUpdate)**"

      !--(DMK-CCATT-INI)-----------------------------------------------------
      !initial 4 related
      character(len=2)                        :: cng
      character(len=f_name_length)                :: fileName

      logical                                :: sameGrid
      integer, external                        :: cio_i
      integer, external                        :: cio_f
      integer, external                        :: RAMS_getvar
      integer                             :: imode
      integer                                :: ngr
      integer                                :: ie
      integer                                :: iunhd
      integer                                :: inhunt
      integer, dimension(1)                :: ngrids1
      integer                                 :: nzg1
      integer                                 :: nzs1
      integer                                 :: npatch1
      integer                                 :: ngrid1
      integer                                :: maxarr
      integer                                :: maxarr2
      integer                                :: maxx1
      integer                                :: maxy1
      integer                                :: maxz1
      real                                 :: time1
      real                                :: ztop1
      character(len=2)                        :: cmode
      integer :: iTime1
      integer :: iZtop1

      integer, allocatable, dimension(:)        :: nnxp1
      integer, allocatable, dimension(:)        :: nnyp1
      integer, allocatable, dimension(:)        :: nnzp1
      integer :: nm

      real, allocatable, dimension(:)         :: platn1
      real, allocatable, dimension(:)         :: plonn1
      real, allocatable, dimension(:)         :: scr
      real, allocatable, dimension(:)         :: scr2
      real, allocatable, dimension(:)         :: scr3
      real, allocatable, dimension(:, :)         :: topt1
      real, allocatable, dimension(:, :)        :: xmn1
      real, allocatable, dimension(:, :)        :: xtn1
      real, allocatable, dimension(:, :)        :: ymn1
      real, allocatable, dimension(:, :)        :: ytn1
      real, allocatable, dimension(:, :)        :: zmn1
      real, allocatable, dimension(:, :)        :: ztn1
      integer :: recordLen, irec
      integer, external :: outRealSize
      character(len=2) :: cmynum
      character(len=6) :: ctime
      !--(DMK-CCATT-FIM)-----------------------------------------------------

      !      Check and see what we are doing. If it is initial time, read
      !        fields into regular arrays. If not, see if nudging will be done
      !        on this grid if it is a nested grid.
      if (ngrid > 1 .and. tnudcent + tnudtop < .001 .and. initflag == 0) return

      ! Put new fields into varinit future arrays. If iswap == 1,
      !     swap future into past first

      ! print*,"====>In",trim(h),time,initialFlag

      if (iswap == 1) then
         varinit_g(ngrid)%varup(:, :, :) = varinit_g(ngrid)%varuf(:, :, :)
         varinit_g(ngrid)%varvp(:, :, :) = varinit_g(ngrid)%varvf(:, :, :)
         varinit_g(ngrid)%varpp(:, :, :) = varinit_g(ngrid)%varpf(:, :, :)
         varinit_g(ngrid)%vartp(:, :, :) = varinit_g(ngrid)%vartf(:, :, :)
         varinit_g(ngrid)%varrp(:, :, :) = varinit_g(ngrid)%varrf(:, :, :)

         !--(DMK-CCATT-INI)-----------------------------------------------------
         if (chem_assim == 1 .and. chemistry >= 0) then
            do nspc = 1, nspecies
               if (spc_alloc(fdda, nspc) == 1) &
                  chem1_g(nspc, ngrid)%sc_pp(:, :, :) = chem1_g(nspc, ngrid)%sc_pf(:, :, :)
            end do
         end if

         if (aer_assim == 1 .and. aerosol >= 1 .and. chemistry >= 0) then
            do nspc = 1, aer_nspecies
               do imode = 1, aer_nmodes
                  if (aer_alloc(aer_fdda, imode, nspc) == aer_on) then
                     aer1_g(imode, nspc, ngrid)%sc_pp(:, :, :) = aer1_g(imode, nspc, ngrid)%sc_pf(:, :, :)
                  end if
               end do
            end do
         end if

         !--(DMK-CCATT-FIM)-----------------------------------------------------

      end if

      ! master makes data file name from tag file name and checks existence

      if (mchnum == master_num) then

         !--(DMK-CCATT-INI)-----------------------------------------------------
         if (initialFlag .eq. 2) then
            !--(DMK-CCATT-FIM)-----------------------------------------------------

            write (cgrid, '(a2,i1,a4)') '-g', ngrid, '.vfm'
            nc = len_trim(fnames_varf(nvarffl))
            flnm = fnames_varf(nvarffl) (1:nc - 4)//trim(cgrid)

            inquire (file=flnm(1:len_trim(flnm)), exist=there)

            !print*,"====>In",trim(h),time,initialFlag,trim(flnm)

            ! Gotta have grid 1...
            if (.not. there .and. ngrid == 1) then
               print *
               call fatal_error(h//" No grid 1 varfile ("//trim(flnm)//") found")
            end if

            !--(DMK-CCATT-INI)-----------------------------------------------------
         else if (initialFlag .eq. 4) then
            iunhd = 11
            inhunt = 10

            fileName = fnames_varf(nvarffl) !fnames_varf(1)

            inquire (file=trim(fileName), exist=there)

            if (.not. there) then
               print *
               call fatal_error(h//" No grid 1 varfile ("//trim(fileName)//") found")
               stop ! # RMF: isnt the right way to stop the model for sure!
            end if

            call rams_read_header(fileName(1:len_trim(fileName) - 9))

            call rams_f_open(iunhd, fileName, 'FORMATTED', 'OLD', 'READ', 0)

            ie = cio_i(iunhd, 1, 'ngrids', ngrids1(1), 1)
            !print*,"====>In",trim(h),time,initialFlag,trim(flnm)

         end if
         !--(DMK-CCATT-FIM)-----------------------------------------------------

         if (there) then
            ifileok = 1
         else
            ifileok = 0
         end if
      end if

      ! all processes decide on continuing computing

      call Broadcast(ifileok, master_num, "ifileok")
      if (ifileok == 0) then
         return
      end if

      !--(DMK-CCATT-INI)-----------------------------------------------------
      if (initialFlag .eq. 2) then
         !--(DMK-CCATT-FIM)-----------------------------------------------------

         ! master opens var file to find out var file version
         ! and to check adequacy to specified problem

         if (mchnum == master_num) then
            iun = 22
            open (unit=42, file='brams.log', position='append', action='write')
            write (unit=42, fmt='(A)') "======================================================================="
            write (unit=42, fmt='(A,1X,A)') "in:", trim(h)
            write (unit=42, fmt='(a,2x,f10.0,2x,a)') ' time (s)', time, trim(flnm(1:len_trim(flnm)))
            write (unit=42, fmt='(A)') "======================================================================="
            close (unit=42)

            call rams_f_open(iun, flnm(1:len_trim(flnm)), 'FORMATTED', 'OLD', 'READ', 0)

            ! Find varfile "version"

            read (iun, *) imarker
            rewind (iun)

            if (imarker == 999999) then
               read (iun, *) imarker, iver_var
            else
               iver_var = 1
            end if

            read (iun, *) iyearx, imonthx, idatex, ihourx &
               , nxpx, nypx, nzpx, rlatx, wlon1x, deltaxx, deltayx, deltazx &
               , dzratx, dzmaxx

            if (nxp .ne. nxpx .or. &
                nyp .ne. nypx .or. &
                nzp .ne. nzpx .or. &
                abs(deltax - deltaxx) .gt. .001 .or. &
                abs(deltay - deltayx) .gt. .001 .or. &
                abs(deltaz - deltazx) .gt. .001 .or. &
                abs(dzrat - dzratx) .gt. .001 .or. &
                abs(dzmax - dzmaxx) .gt. .001 .or. &
                abs(platn(ngrid) - rlatx) .gt. .001 .or. &
                abs(plonn(ngrid) - wlon1x) .gt. .001) then

               call fatal_error(h//" grid mismatch between varfile ("//trim(flnm)//") and namelist")
            end if

            !--(DMK-CCATT-INI)-----------------------------------------------------
            !- test if the source data is for the chemical mechanism that will be used:
            if (chem_assim == 1 .and. chemistry >= 0) then
               read (iun, *) chemical_mechanism_test
               if (trim(chemical_mechanism_test) /= trim(chemical_mechanism)) then
                  call fatal_error(h//" Wrong Chem mechanism in Varfiles. Expected="// &
                                   trim(chemical_mechanism(1:len_trim(chemical_mechanism)))//" Read="// &
                                   trim(chemical_mechanism_test(1:len_trim(chemical_mechanism_test))))
               end if
            end if
            !--(DMK-CCATT-END)-----------------------------------------------------

            !  !-- SRF added to reuse ivar files with chemical fields but
            !  !-- only reading the meteo fields (meteo runs only).
            !  if (chem_assim == 0 .and. chemistry == -1 .and. maxval(spc_alloc(fdda,:))>0) then
            !     read(iun,*)  chemical_mechanism_test
            !     print*,"Reading IVAR file for chemical mechanism ", trim( chemical_mechanism_test )
            !     print*,"but model will only consider the meteo fields for assimilation"
            !  endif

         end if

         ! All processes need to know var file version

         call Broadcast(iver_var, master_num, "iver_var")

         ! deals with varfile fields into the "future" varinit arrays,
         ! to be swapped to the past arrays when needed.
         call ReadStoreOwnChunk(ngrid, iun, varinit_g(ngrid)%varuf, nzp, "varuf")
         call ReadStoreOwnChunk(ngrid, iun, varinit_g(ngrid)%varvf, nzp, "varvf")
         call ReadStoreOwnChunk(ngrid, iun, varinit_g(ngrid)%varpf, nzp, "varpf")
         call ReadStoreOwnChunk(ngrid, iun, varinit_g(ngrid)%vartf, nzp, "vartf")
         call ReadStoreOwnChunk(ngrid, iun, varinit_g(ngrid)%varrf, nzp, "varrf")

         varinit_g(ngrid)%varrf(:, :, :) = max(1.e-8, varinit_g(ngrid)%varrf(:, :, :))

         if (chem_assim == 1 .and. chemistry >= 0) then
            if (mchnum == master_num) then
               print *, '--------------------------------------------------------------------------'
               print *, ' 4DDA using chem mechanism= ', trim(chemical_mechanism(1:len_trim(chemical_mechanism)))
               print *, ' specie name - max - min values in assimilated dataset'
            end if

            do nspc = 1, nspecies
               if (spc_alloc(fdda, nspc) == 1) then

                  call ReadStoreOwnChunk(ngrid, iun, chem1_g(nspc, ngrid)%sc_pf, nzp, "sc_pf")

                  !no futuro cheque o limite inferior (e-23)
                  chem1_g(nspc, ngrid)%sc_pf(:, :, :) = init_ajust(nspc)*max(1.e-23, chem1_g(nspc, ngrid)%sc_pf(:, :, :))
                if (mchnum == master_num) print *, ' chem spc=', nspc, spc_name(nspc), maxval(chem1_g(nspc, ngrid)%sc_pf(:, :, :)) &
                     , minval(chem1_g(nspc, ngrid)%sc_pf(:, :, :))
               end if
            end do
            do nspc = 1, aer_nspecies
               do nm = 1, aer_nmodes
                  if (aer_alloc(aer_fdda, nm, nspc) == 1) then

                     call ReadStoreOwnChunk(ngrid, iun, aer1_g(nm, nspc, ngrid)%sc_pf, nzp, "sc_pf")

                     !no futuro cheque o limite inferior (e-23)
                     aer1_g(nm, nspc, ngrid)%sc_pf(:, :, :) = init_ajust(nspc)*max(1.e-23, aer1_g(nm, nspc, ngrid)%sc_pf(:, :, :))
                     if (mchnum == master_num) print *, ' aer spc=', nm, nspc, aer_mod_spc(nm, nspc) &
                                                , maxval(aer1_g(nm, nspc, ngrid)%sc_pf(:, :, :)) &
                        , minval(aer1_g(nm, nspc, ngrid)%sc_pf(:, :, :))
                  end if
               end do
            end do

            if (mchnum == master_num) &
               print *, '--------------------------------------------------------------------------'
         end if

         ! master skips a few fields and deals with snow_mass
         if (initflag == 1 .and. iver_var == 2) then
            ! Extract snow depth from the varfile. Ignore other 2D fields for now.
            if (mchnum == master_num) then
               call vfirec(iun, scratch, nxyp, 'LIN')
               call vfirec(iun, scratch, nxyp, 'LIN')
               call vfirec(iun, scratch, nxyp, 'LIN')
            end if
            call ReadStoreOwnChunk(ngrid, iun, leaf_g(ngrid)%snow_mass, "snow_mass")

            if (mchnum == master_num) then
               call vfirec(iun, scratch, nxyp, 'LIN')
               close (iun)
            end if
         end if

         !--(DMK-CCATT-INI)-----------------------------------------------------
      elseif (initialFlag .eq. 4) then

         call Broadcast(ngrids1(1), master_num, 'ngrids1(1)')

         allocate (nnxp1(ngrids1(1)))
         allocate (nnyp1(ngrids1(1)))
         allocate (nnzp1(ngrids1(1)))
         allocate (platn1(ngrids1(1)))
         allocate (plonn1(ngrids1(1)))

         if (mchnum == master_num) then

            ie = cio_i(iunhd, 1, 'nnxp', nnxp1, ngrids1(1))
            ie = cio_i(iunhd, 1, 'nnyp', nnyp1, ngrids1(1))
            ie = cio_i(iunhd, 1, 'nnzp', nnzp1, ngrids1(1))
            ie = cio_i(iunhd, 1, 'npatch', npatch1, 1)
            ie = cio_i(iunhd, 1, 'nzg', nzg1, 1)
            ie = cio_i(iunhd, 1, 'nzs', nzs1, 1)
            ie = cio_f(iunhd, 1, 'time', time1, 1)
            ie = cio_f(iunhd, 1, 'ztop', ztop1, 1)
            ie = cio_f(iunhd, 1, 'platn', platn1, ngrids1(1))
            ie = cio_f(iunhd, 1, 'plonn', plonn1, ngrids1(1))

            if (nzg .ne. nzg1 .or. nzs .ne. nzs1 .or. npatch .ne. npatch1) then
               print *, 'LEAF parameters must be same for initial-history start'
               print *, 'npatch: ', npatch, npatch1
               print *, 'nzg: ', nzg, nzg1
               print *, 'nzs: ', nzs, nzs1
               stop 'LEAF hist-init'
            end if

         end if

         iTime1 = int(time1)
         iZtop1 = int(ztop1)
         call Broadcast(ngrids1(1), master_num, 'ngrids1(1)')
         call Broadcast(nnxp1, master_num, 'nnxp1')
         call Broadcast(nnyp1, master_num, 'nnyp1')
         call Broadcast(nnzp1, master_num, 'nnzp1')
         call Broadcast(npatch1, master_num, 'npatch1')
         call Broadcast(nzg1, master_num, 'nzg1')
         call Broadcast(nzs1, master_num, 'nzs1')
         call Broadcast(platn1, master_num, 'platn1')
         call Broadcast(plonn1, master_num, 'plonn1')
         call Broadcast(itime1, master_num, 'time1')
         call Broadcast(iztop1, master_num, 'ztop1')
         time1 = real(iTime1)
         ztop1 = real(iztop1)

         maxarr = 0
         maxarr2 = 0
         maxx1 = 0
         maxy1 = 0
         maxz1 = 0

         ! Find maximum size of any array on history file. Allocate scratch array of
         ! this size.
         ! # RMF: is it really necessary? only interpolates grid1
         do ngr = 1, ngrids1(1)
            maxarr = max(maxarr, nnxp1(ngr)*nnyp1(ngr)*nnzp1(ngr), &
                         nnxp1(ngr)*nnyp1(ngr)*nzg1*npatch1, &
                         nnxp1(ngr)*nnyp1(ngr)*nzs1*npatch1) ! # RMF: this 'nwave' var is used by CCATT only -  nnxp1(ngr)*nnyp1(ngr)*nwave)!srf

            maxarr2 = max(maxarr2, nnxp1(ngr)*nnyp1(ngr))
            maxx1 = max(maxx1, nnxp1(ngr))
            maxy1 = max(maxy1, nnyp1(ngr))
            maxz1 = max(maxz1, nnzp1(ngr))
         end do

         allocate (xmn1(maxx1, ngrids1(1)))
         allocate (xtn1(maxx1, ngrids1(1)))
         allocate (ymn1(maxy1, ngrids1(1)))
         allocate (ytn1(maxy1, ngrids1(1)))
         allocate (zmn1(maxz1, ngrids1(1)))
         allocate (ztn1(maxz1, ngrids1(1)))

         ngrid1 = 1

         if (mchnum == master_num) then

            do ngr = 1, ngrids1(1)
               write (cng, '(i2.2)') ngr
               ie = cio_f(iunhd, 1, 'xmn'//cng, xmn1(1, ngr), nnxp1(ngr))
               ie = cio_f(iunhd, 1, 'xtn'//cng, xtn1(1, ngr), nnxp1(ngr))
               ie = cio_f(iunhd, 1, 'ymn'//cng, ymn1(1, ngr), nnyp1(ngr))
               ie = cio_f(iunhd, 1, 'ytn'//cng, ytn1(1, ngr), nnyp1(ngr))
               ie = cio_f(iunhd, 1, 'zmn'//cng, zmn1(1, ngr), nnzp1(ngr))
               ie = cio_f(iunhd, 1, 'ztn'//cng, ztn1(1, ngr), nnzp1(ngr))
            end do
            close (iunhd)
            close (inhunt)

         end if

         sameGrid = .false.

         if ((platn1(ngrid1) .eq. platn(ngrid1)) .and. (plonn1(ngrid1) .eq. &
                                                        plonn(ngrid1)) .and. (nnxp1(ngrid1) .eq. nnxp(ngrid1)) .and. &
             (nnyp1(ngrid1) .eq. nnyp(ngrid1)) .and. (nnzp1(ngrid1) .eq. &
                                                      nnzp(ngrid1)) .and. (npatch1 .eq. npatch) .and. &
             (ztop1 .eq. ztop1) .and. (nzg1 .eq. nzg) .and. (nzs1 .eq. nzs)) sameGrid = .true.

         call Broadcast(xmn1, master_num, 'xmn1')
         call Broadcast(xtn1, master_num, 'xtn1')
         call Broadcast(ymn1, master_num, 'ymn1')
         call Broadcast(ytn1, master_num, 'ytn1')
         call Broadcast(zmn1, master_num, 'zmn1')
         call Broadcast(ztn1, master_num, 'ztn1')
         call Broadcast(filename, master_num, 'filename')

         allocate (topt1(maxarr2, ngrids1(1)))
         allocate (scr(maxarr))
         allocate (scr2(maxarr))
         allocate (scr3(maxarr))

         if (mchnum == master_num) then

            print *, '% inithist-anl: ', trim(fileName(1:len_trim(fileName) - 9))

            ! # RMF: in initial 4 case must always be grid 1

            ie = RAMS_getvar('TOPT', 1, scr, scr2, trim(fileName(1:len_trim(fileName) - 9)))

            call unarrange(1, nnxp1(ngrid1), nnyp1(ngrid1), scr, topt1(:, 1))

         end if

         call Broadcast(topt1, master_num, 'topt1')

         ! ##### UP

         if (mchnum == master_num) then

            ie = RAMS_getvar('UP', 1, scr, scr2, trim(fileName(1:len_trim(fileName) - 9)))

         end if

         call Broadcast(scr, master_num, 'scr')

         if (.not. sameGrid) then

            !write(*,*) 'UP antes: ',mchnum,time,maxval(scr),minval(scr)

            call hi_interpInitial4(nnzp1(ngrid1), nnxp1(ngrid1), nnyp1(ngrid1), scr, &
                                   xmn1, xtn1, ymn1, ytn1, zmn1, ztn1, platn1(ngrid1), plonn1(ngrid1), &
                                   topt1, ztop1, mzp, mxp, myp, & !mxp,myp,  &
                                   varinit_g(ngrid1)%varuf, ngrid1, ngrid1, 'UP', 3)
! if(mchnum == 23) call writevar(varinit_g(ngrid1)%varuf,mxp,myp,mzp,'escr23-'//ctime)
!                write(*,*) 'varuf   :',mchnum,time,maxval(varinit_g(ngrid1)%varuf),minval(varinit_g(ngrid1)%varuf)
         else
            call unarrange(nnzp1(ngrid1), nnxp1(ngrid1), nnyp1(ngrid1), scr, scr3)
            call storeOwnChunk_3D(ngrid1, scr3, varinit_g(ngrid1)%varuf, nnzp(1), nnxp(1), nnyp(1), 'UP')
         end if

         if (mchnum == master_num) then

            ie = RAMS_getvar('VP', 1, scr, scr2, trim(fileName(1:len_trim(fileName) - 9)))

         end if

         call Broadcast(scr, master_num, 'scr')

         if (.not. sameGrid) then

            call hi_interpInitial4(nnzp1(ngrid1), nnxp1(ngrid1), nnyp1(ngrid1), scr, &
                                   xmn1, xtn1, ymn1, ytn1, zmn1, ztn1, platn1(ngrid1), plonn1(ngrid1), &
                                   topt1, ztop1, mzp, mxp, myp, &
                                   varinit_g(ngrid1)%varvf, ngrid1, ngrid1, 'VP', 3)
         else

            call unarrange(nnzp1(ngrid1), nnxp1(ngrid1), nnyp1(ngrid1), scr, scr3)
            call storeOwnChunk_3D(ngrid1, scr3, varinit_g(ngrid1)%varvf, nnzp(1), nnxp(1), nnyp(1), 'VP')
         end if

         ! ##### THETA

         if (mchnum == master_num) then

            ie = RAMS_getvar('THETA', 1, scr, scr2, trim(fileName(1:len_trim(fileName) - 9)))

         end if

         call Broadcast(scr, master_num, 'scr')

         if (.not. sameGrid) then

            call hi_interpInitial4(nnzp1(ngrid1), nnxp1(ngrid1), nnyp1(ngrid1), scr, &
                                   xmn1, xtn1, ymn1, ytn1, zmn1, ztn1, platn1(ngrid1), plonn1(ngrid1), &
                                   topt1, ztop1, mzp, mxp, myp, &
                                   varinit_g(ngrid1)%vartf, ngrid1, ngrid1, 'THETA', 3)
         else

            call unarrange(nnzp1(ngrid1), nnxp1(ngrid1), nnyp1(ngrid1), scr, scr3)
            call storeOwnChunk_3D(ngrid1, scr3, varinit_g(ngrid1)%vartf, nnzp(1), nnxp(1), nnyp(1), 'THETA')
         end if

         ! ##### PI

         if (mchnum == master_num) then

            ie = RAMS_getvar('PI', 1, scr, scr2, trim(fileName(1:len_trim(fileName) - 9)))

         end if

         call Broadcast(scr, master_num, 'scr')

         if (.not. sameGrid) then

            call hi_interpInitial4(nnzp1(ngrid1), nnxp1(ngrid1), nnyp1(ngrid1), scr, &
                                   xmn1, xtn1, ymn1, ytn1, zmn1, ztn1, platn1(ngrid1), plonn1(ngrid1), &
                                   topt1, ztop1, mzp, mxp, myp, &
                                   varinit_g(ngrid1)%varpf, ngrid1, ngrid1, 'PI', 3)
         else
      !!  write(*,*) 'LFR - DEBUG: ','VarfUpdate - 1154'
            call unarrange(nnzp1(ngrid1), nnxp1(ngrid1), nnyp1(ngrid1), scr, scr3)
            call storeOwnChunk_3D(ngrid1, scr3, varinit_g(ngrid1)%varpf, nnzp(1), nnxp(1), nnyp(1), 'PI')
         end if

         ! ##### RV

         if (mchnum == master_num) then

            ie = RAMS_getvar('RV', 1, scr, scr2, trim(fileName(1:len_trim(fileName) - 9)))

         end if

         call Broadcast(scr, master_num, 'scr')
         if (.not. sameGrid) then

            call hi_interpInitial4(nnzp1(ngrid1), nnxp1(ngrid1), nnyp1(ngrid1), scr, &
                                   xmn1, xtn1, ymn1, ytn1, zmn1, ztn1, platn1(ngrid1), plonn1(ngrid1), &
                                   topt1, ztop1, mzp, mxp, myp, &
                                   varinit_g(ngrid1)%varrf, ngrid1, ngrid1, 'RV', 3)
         else

            call unarrange(nnzp1(ngrid1), nnxp1(ngrid1), nnyp1(ngrid1), scr, scr3)
            call storeOwnChunk_3D(ngrid1, scr3, varinit_g(ngrid1)%varrf, nnzp(1), nnxp(1), nnyp(1), 'RV')
         end if

         varinit_g(ngrid)%varrf(:, :, :) = max(1.e-8, varinit_g(ngrid)%varrf(:, :, :))

         ! #### CHEM
         if (CHEM_ASSIM == on .and. CHEMISTRY >= 0) then
            do nspc = 1, nspecies
               !print*, trim(vtab_r(ni,ng)%name), '>>', trim(spc_name(nspc))//'P'
               if (spc_alloc(fdda, nspc) == on) then

                  if (mchnum == master_num) then
                     ie = RAMS_getvar(trim(spc_name(nspc))//'P', 1, scr, scr2, trim(fileName(1:len_trim(fileName) - 9)))
                  end if

                  call Broadcast(scr, master_num, 'scr')
                  if (.not. sameGrid) then

                     call hi_interpInitial4(nnzp1(ngrid1), nnxp1(ngrid1), nnyp1(ngrid1), scr, &
                                            xmn1, xtn1, ymn1, ytn1, zmn1, ztn1, platn1(ngrid1), plonn1(ngrid1), &
                                            topt1, ztop1, mzp, mxp, myp, &
                                            chem1_g(nspc, ngrid1)%sc_pf, ngrid1, ngrid1, trim(spc_name(nspc))//'P', 3)
                  else

                     call unarrange(nnzp1(ngrid1), nnxp1(ngrid1), nnyp1(ngrid1), scr, scr3)
                     call storeOwnChunk_3D(ngrid1, scr3, chem1_g(nspc, ngrid1)%sc_pf, nnzp(1), nnxp(1), nnyp(1) &
                                          , trim(spc_name(nspc))//'P')
                  end if

                  chem1_g(nspc, ngrid1)%sc_pf(:, :, :) = max(1.e-23, chem1_g(nspc, ngrid1)%sc_pf(:, :, :))

               end if
            end do
         end if

         ! ### AERO
         if (AER_ASSIM == on .and. AEROSOL >= 1 .and. CHEMISTRY >= 0) then
            do nspc = 1, aer_nspecies
               do imode = 1, aer_nmodes
                  write (cmode, '(BN, I2)') imode
                  cmode = adjustl(cmode)
                  if (aer_alloc(aer_fdda, imode, nspc) == aer_on) then

                     if (mchnum == master_num) then
                        ie = RAMS_getvar(trim(aer_name(nspc))//trim(cmode)//'P', 1, scr, scr2, &
                                         trim(fileName(1:len_trim(fileName) - 9)))
                     end if

                     call Broadcast(scr, master_num, 'scr')
                     if (.not. sameGrid) then

                        call hi_interpInitial4(nnzp1(ngrid1), nnxp1(ngrid1), nnyp1(ngrid1), scr, &
                                               xmn1, xtn1, ymn1, ytn1, zmn1, ztn1, platn1(ngrid1), plonn1(ngrid1), &
                                               topt1, ztop1, mzp, mxp, myp, &
                                       aer1_g(imode, nspc, ngrid1)%sc_pf, ngrid1, ngrid1, trim(aer_name(nspc))//trim(cmode)//'P', 3)
                     else

                        call unarrange(nnzp1(ngrid1), nnxp1(ngrid1), nnyp1(ngrid1), scr, scr3)
                        call storeOwnChunk_3D(ngrid1, scr3, aer1_g(imode, nspc, ngrid1)%sc_pf, nnzp(1), &
                                              nnxp(1), nnyp(1), trim(aer_name(nspc))//trim(cmode)//'P')
                     end if

                     aer1_g(imode, nspc, ngrid1)%sc_pf(:, :, :) = max(1.e-23, aer1_g(imode, nspc, ngrid1)%sc_pf(:, :, :))

                  end if
               end do
            end do
         end if

      end if !if(initialFlag .eq. 2)
      !--(DMK-CCATT-FIM)-----------------------------------------------------

      ! If running ADAP coord, do interpolation to Cartesian levels

      if (if_adap == 1) then

         !**(JP)** not converted

         call fatal_error(h//"**(JP)** varf_adap not converted")
         call varf_adap(nnzp(ngrid), nnxp(ngrid), nnyp(ngrid) &
                        , varinit_g(ngrid)%varuf, varinit_g(ngrid)%varvf &
                        , varinit_g(ngrid)%varpf, varinit_g(ngrid)%vartf &
                        , varinit_g(ngrid)%varrf, grid_g(ngrid)%topta)
      end if

      ! Find the reference state

      if (initflag == 1 .and. ngrid == 1) then
         call RefVar(OneGrid, mzp, mxp, myp &
                     , varinit_g(ngrid)%vartf, varinit_g(ngrid)%varpf &
                     , basic_g(ngrid)%pi0, basic_g(ngrid)%th0 &
                     , varinit_g(ngrid)%varrf, basic_g(ngrid)%dn0 &
                     , basic_g(ngrid)%dn0u, basic_g(ngrid)%dn0v &
                     , varinit_g(ngrid)%varuf, varinit_g(ngrid)%varvf &
                     , grid_g(ngrid)%topt, grid_g(ngrid)%rtgt &
                     , grid_g(ngrid)%topta, level)
      end if

      varinit_g(ngrid)%varpf(:, :, :) = &
         varinit_g(ngrid)%varpf(:, :, :) - basic_g(ngrid)%pi0(:, :, :)

      ! If this is an initialization, put data into regular arrays

      if (initflag == 1) then
         basic_g(ngrid)%uc(:, :, :) = varinit_g(ngrid)%varuf(:, :, :)
         basic_g(ngrid)%vc(:, :, :) = varinit_g(ngrid)%varvf(:, :, :)
         basic_g(ngrid)%pc(:, :, :) = varinit_g(ngrid)%varpf(:, :, :)
         basic_g(ngrid)%thp(:, :, :) = varinit_g(ngrid)%vartf(:, :, :)
         basic_g(ngrid)%rtp(:, :, :) = varinit_g(ngrid)%varrf(:, :, :)

         !--(DMK-CCATT-INI)------------------------------------------------------
         if (chem_assim == 1 .and. chemistry >= 0) then
            do nspc = 1, nspecies
               if (spc_alloc(fdda, nspc) == 1) &
                  chem1_g(nspc, ngrid)%sc_p(:, :, :) = chem1_g(nspc, ngrid)%sc_pf(:, :, :)
            end do
         end if
         !--(DMK-CCATT-FIM)------------------------------------------------------
         if (aer_assim == 1 .and. aerosol >= 1 .and. chemistry >= 0) then
            do nspc = 1, aer_nspecies
               do imode = 1, aer_nmodes
                  if (aer_alloc(aer_fdda, imode, nspc) == aer_on) then
                     aer1_g(imode, nspc, ngrid)%sc_p(:, :, :) = aer1_g(imode, nspc, ngrid)%sc_pf(:, :, :)
                  end if
               end do
            end do
         end if
      end if

   end subroutine VarfUpdate

   subroutine RefVar(OneGrid, &
                     n1, n2, n3, thp, pc, pi0, th0, rtp, dn0, dn0u, dn0v, uc, &
                     vc, topt, rtgt, topta, level)

      ! Arguments:
      type(Grid), pointer :: OneGrid
      integer, intent(in)  :: n1
      integer, intent(in)  :: n2
      integer, intent(in)  :: n3
      integer, intent(in)  :: level
      real, intent(in)  :: thp(n1, n2, n3)
      real, intent(in)  :: pc(n1, n2, n3)
      real, intent(out) :: pi0(n1, n2, n3)
      real, intent(in)  :: rtp(n1, n2, n3)
      real, intent(out) :: dn0(n1, n2, n3)
      real, intent(out) :: dn0u(n1, n2, n3)
      real, intent(out) :: dn0v(n1, n2, n3)
      real, intent(in)  :: uc(n1, n2, n3)
      real, intent(in)  :: vc(n1, n2, n3)
      real, intent(in)  :: topt(n2, n3)
      real, intent(in)  :: rtgt(n2, n3)
      real, intent(out) :: th0(n1, n2, n3)
      real, intent(in)  :: topta(n2, n3)
      ! local Variables:
      integer :: i, j, k
      integer :: rankMin
      real :: vctr1(nzp)
      real :: vctr2(nzp)
      real :: vctr4(nzp)
      integer :: indRefVar(2)
      logical, parameter :: dumpLocal = .false.
      character(len=16) :: c0, c1, c2, c3
      character(len=*), parameter :: h = "**(RefVar)**"

      real :: global_data(nnzp(1), nnxp(1), nnyp(1))
      real :: global_data2d(nnxp(1), nnyp(1))

      real, allocatable :: thp_global(:, :, :), uc_global(:, :, :), vc_global(:, :, :), &
                           rtp_global(:, :, :), pc_global(:, :, :), top_global(:, :)

      real :: thp_ref(n1), uc_ref(n1), vc_ref(n1), rtp_ref(n1), pc_ref(n1), top_ref
      integer :: nxyp, ierr
      character(LEN=16)  :: varn

      ! Reference sounding is point with lowest topography
      ! All processes have lowest topography point (topref);
      ! Only process that owns the point (rankMin) has point coordinates

      ! Calculating in serial mode to keep binary reprodutibility in real sum and
      ! average procedures
      ! Allocating global data

      nxyp = nnxp(1)*nnyp(1)

      ! master process gathers global data and computes average per level

      ! for THP

      varn = "THP"
      call gatherData(3, varn, 1, nnzp(1), nnxp(1), nnyp(1), &
                      nmachs, mchnum, mynum, master_num, &
                      thp, global_data)
      if (mchnum == master_num) then
         do k = 1, n1
            thp_ref(k) = sum(global_data(k, :, :))/real(nxyp)
         end do
      end if
      call parf_barrier(0)

      ! for UC

      varn = "UC"
      call gatherData(3, varn, 1, nnzp(1), nnxp(1), nnyp(1), &
                      nmachs, mchnum, mynum, master_num, &
                      uc, global_data)
      if (mchnum == master_num) then
         do k = 1, n1
            uc_ref(k) = sum(global_data(k, :, :))/real(nxyp)
         end do
      end if
      call parf_barrier(0)

      ! for VC

      varn = "VC"
      call gatherData(3, varn, 1, nnzp(1), nnxp(1), nnyp(1), &
                      nmachs, mchnum, mynum, master_num, &
                      vc, global_data)
      if (mchnum == master_num) then
         do k = 1, n1
            vc_ref(k) = sum(global_data(k, :, :))/real(nxyp)
         end do
      end if
      call parf_barrier(0)

      ! for RTP

      varn = "RTP"
      call gatherData(3, varn, 1, nnzp(1), nnxp(1), nnyp(1), &
                      nmachs, mchnum, mynum, master_num, &
                      rtp, global_data)
      if (mchnum == master_num) then
         do k = 1, n1
            rtp_ref(k) = sum(global_data(k, :, :))/real(nxyp)
         end do
      end if
      call parf_barrier(0)

      ! for PC

      varn = "PC"
      call gatherData(3, varn, 1, nnzp(1), nnxp(1), nnyp(1), &
                      nmachs, mchnum, mynum, master_num, &
                      pc, global_data)
      if (mchnum == master_num) then
         do k = 1, n1
            pc_ref(k) = sum(global_data(k, :, :))/real(nxyp)
         end do
      end if
      call parf_barrier(0)

      ! for TOPTA

      varn = "TOPTA"
      call gatherData(2, varn, 1, nnxp(1), nnyp(1), &
                      nmachs, mchnum, mynum, master_num, &
                      topta, global_data2d)
      if (mchnum == master_num) then
         top_ref = sum(global_data2d(:, :))/real(nxyp)
      end if
      call parf_barrier(0)

      ! Broadcast reference data:

      call parf_bcast(thp_ref, int(n1, i8), master_num)
      call parf_bcast(uc_ref, int(n1, i8), master_num)
      call parf_bcast(vc_ref, int(n1, i8), master_num)
      call parf_bcast(rtp_ref, int(n1, i8), master_num)
      call parf_bcast(pc_ref, int(n1, i8), master_num)
      call parf_bcast(top_ref, master_num)

      topref = top_ref ! Global module information

      if (dumpLocal) then
         print *, h//"mynum,topref=", mynum, topref
         print *, h//"mynum,thp_ref=", mynum, thp_ref
         print *, h//"mynum,uc_ref=", mynum, uc_ref
         print *, h//"mynum,vc_ref=", mynum, vc_ref
         print *, h//"mynum,rtp_ref=", mynum, rtp_ref
         print *, h//"mynum,pc_ref=", mynum, pc_ref
      end if

      if (if_adap == 0) then
         do k = 1, nzp
            vctr2(k) = ztn(k, ngrid)*(1.-topref/ztop) + topref
         end do
         if (dumpLocal) then
            print *, h//"mynum, vctr2=", mynum, vctr2
         end if
         call htint2(nzp, thp_ref(1), vctr2, nzp, vctr1, zt)
         if (dumpLocal) then
            print *, h//"new mynum, vctr1=", mynum, vctr1
         end if
         call htint2(nzp, uc_ref(1), vctr2, nzp, u01dn(1, ngrid), zt)
         call htint2(nzp, vc_ref(1), vctr2, nzp, v01dn(1, ngrid), zt)
         if (level >= 1) then
            call htint2(nzp, rtp_ref(1), vctr2, nzp, rt01dn(1, ngrid), zt)
         else
            rt01dn(1:nzp, ngrid) = 0.
         end if
      else
         ! *********Check latter for Shaved-ETA:
         vctr2(1:nzp) = ztn(1:nzp, ngrid)
         vctr1(1:nzp) = thp_ref(1:nzp) !thp(1:nzp, iref, jref)
         u01dn(1:nzp, ngrid) = uc_ref(1:nzp) !uc(1:nzp,iref,jref)
         v01dn(1:nzp, ngrid) = vc_ref(1:nzp) !vc(1:nzp,iref,jref)
         rt01dn(1:nzp, ngrid) = 0.
         if (level >= 1) rt01dn(1:nzp, ngrid) = rtp_ref(1:nzp) !rtp(1:nzp,iref,jref)
      end if

      do k = 1, nzp
         th01dn(k, ngrid) = vctr1(k)*(1.+0.61*rt01dn(k, ngrid))
      end do
      u01dn(1, ngrid) = u01dn(2, ngrid)
      v01dn(1, ngrid) = v01dn(2, ngrid)
      rt01dn(1, ngrid) = rt01dn(2, ngrid)
      th01dn(1, ngrid) = th01dn(2, ngrid)

      pi01dn(1, ngrid) = pc_ref(1) + g*(vctr2(1) - zt(1))/ &
                         (0.5*(th01dn(1, ngrid) + thp_ref(1)*(1.+0.61*rtp_ref(1))))
      do k = 2, nzp
         if (dumpLocal) then
            print *, h//"mynum,nzp,k,dzm(k-1),th01dn(k,ngrid),th01dn(k-1,ngrid)=", &
               mynum, nzp, k, dzm(k - 1), th01dn(k, ngrid), th01dn(k - 1, ngrid)
         end if
         pi01dn(k, ngrid) = pi01dn(k - 1, ngrid) - g/(dzm(k - 1)*0.5* &
                                                      (th01dn(k, ngrid) + th01dn(k - 1, ngrid)))
      end do

      do k = 1, nzp
         vctr4(k) = (pi01dn(k, ngrid)/cp)**cpor*p00
         dn01dn(k, ngrid) = cp*vctr4(k)/(rgas*th01dn(k, ngrid)*pi01dn(k, ngrid))
      end do

      ! Compute 3-D reference state from 1-D reference state

      call PropRefSound3D(OneGrid, nzp, n2, n3, pi0, dn0, dn0u, dn0v, th0, topt, rtgt)

   end subroutine RefVar

   subroutine PropRefSound3D(OneGrid, n1, n2, n3, pi0, dn0, dn0u, dn0v, th0, topt, rtgt)
      ! Arguments:
      type(Grid), pointer :: OneGrid
      integer, intent(in)  :: n1
      integer, intent(in)  :: n2
      integer, intent(in)  :: n3
      real, intent(out) :: pi0(n1, n2, n3)
      real, intent(out) :: dn0(n1, n2, n3)
      real, intent(out) :: dn0u(n1, n2, n3)
      real, intent(out) :: dn0v(n1, n2, n3)
      real, intent(out) :: th0(n1, n2, n3)
      real, intent(in)  :: topt(n2, n3)
      real, intent(in)  :: rtgt(n2, n3)
      ! Local Variables:
      integer :: i, j, k
      real :: c1, c2, c3
      character(len=*), parameter :: h = "**(PropRefSound3D)**"

      ! +---------------------------------------------------------------------
      ! _    This routine initializes the 3-D reference state arrays from the
      !        1-D reference state.
      ! +---------------------------------------------------------------------

      do j = 1, n3
         do i = 1, n2

            if (if_adap == 1) then

               do k = 1, n1
                  pi0(k, i, j) = pi01dn(k, ngrid)
                  th0(k, i, j) = th01dn(k, ngrid)
               end do
               c1 = g*2.

            else

               do k = 1, n1
                  vctr2(k) = zt(k)*rtgt(i, j) + topt(i, j)
               end do
               call htint(nzp, pi01dn(1, ngrid), zt, nzp, pi0(1, i, j), vctr2)
               call htint(nzp, th01dn(1, ngrid), zt, nzp, th0(1, i, j), vctr2)
               c1 = g*2.*(1.-topt(i, j)/ztop)

            end if

            c2 = 1.-cpor
            c3 = cp**c2
            do k = n1 - 1, 1, -1
               pi0(k, i, j) = pi0(k + 1, i, j) + c1/((th0(k, i, j) + th0(k + 1, i, j))*dzm(k))
            end do

            do k = 1, n1
               dn0(k, i, j) = (c3*p00)/(rgas*th0(k, i, j)*pi0(k, i, j)**c2)
            end do

         end do
      end do

      call FillDn0uv(OneGrid, n1, n2, n3, dn0, dn0u, dn0v)

   end subroutine PropRefSound3D

   subroutine FillDn0uv(OneGrid, n1, n2, n3, dn0, dn0u, dn0v)
      ! Arguments:
      type(Grid), pointer :: OneGrid
      integer, intent(in)  :: n1
      integer, intent(in)  :: n2
      integer, intent(in)  :: n3
      real, intent(in)  :: dn0(n1, n2, n3)
      real, intent(out) :: dn0u(n1, n2, n3)
      real, intent(out) :: dn0v(n1, n2, n3)
      ! Local variables:
      integer :: i, j, k, i1, j1
      integer :: ng
      character(len=*), parameter :: h = "**(FillDn0uv)**"

      ! on parallel runs, at processes where n2, n3 are not real borders,
      ! i upper border of dn0u ([n2,:]) is wrong
      ! j upper border of dn0v ([:,n3]) is wrong
      ! since dn0u and dn0v expressions are build for the full domain,
      ! but are computed on local chunks

      do j = 1, n3
         j1 = min(j + 1, n3)
         do i = 1, n2
            i1 = min(i + 1, n2)
            do k = 1, n1
               dn0u(k, i, j) = 0.5*(dn0(k, i, j) + dn0(k, i1, j))
               dn0v(k, i, j) = 0.5*(dn0(k, i, j) + dn0(k, i, j1))
            end do
         end do
      end do

      ! exchange borders to correct wrong values

      call PostRecvSendMsgs(OneGrid%SendDn0u, OneGrid%RecvDn0u)
      call PostRecvSendMsgs(OneGrid%SendDn0v, OneGrid%RecvDn0v)
      call WaitRecvMsgs(OneGrid%SendDn0u, OneGrid%RecvDn0u)
      call WaitRecvMsgs(OneGrid%SendDn0v, OneGrid%RecvDn0v)
   end subroutine FillDn0uv
end module ModVarfFile
