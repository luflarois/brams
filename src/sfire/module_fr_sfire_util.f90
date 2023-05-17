
module module_fr_sfire_util

   implicit none

   integer, save:: &
      fire_print_msg = 1, &
      fire_print_file = 1, &
      fuel_left_method = 1, &
      fuel_left_irl = 2, &
      fuel_left_jrl = 2, &
      boundary_guard = -1, &
      fire_grows_only = 1, &
      sfire_upwinding = 3, &
      fire_test_steps = 0, &
      fire_topo_from_atm = 1, &
      fire_advection = 0, &
      fire_wind_log_interp = 4, &
      fire_use_windrf = 0, &
      !fire_fmc_read = 1, & !LFR-Isilda
      fire_fmc_read = 0, &
      fire_ignition_clamp = 0, &
      fire_hfx_given = 0, &
      fire_hfx_num_lines = 1, &
      fire_update_fuel_frac = 1, &
      fndwi_from_ndwi = 1, &
      kfmc_ndwi = 0, &
      fire_can_top_read = 1

   real, save:: &
      fire_perimeter_time = 0., &
      fire_tign_in_time = 0., &
      fire_atm_feedback = 1., &
      fire_back_weight = 0.5, &
      fire_viscosity = 0.4, &
      fire_lfn_ext_up = 1, &
      fire_hfx_value = 0., &
      fire_hfx_latent_part = 0.084
!!!INTRODUZIDO POR ISILDA CUNHA MENEZES
   logical, save:: &
      crown = .true.
!!!!!

   integer, parameter:: REAL_SUM = 10, REAL_MAX = 20, REAL_MIN = 21, REAL_AMAX = 22, RNRM_SUM = 30, RNRM_MAX = 40

   type line_type
      real ros, &
         stop_time, &
         wind_red, &
         wrdist, &
         wrupwind, &
         start_x, &
         start_y, &
         end_x, &
         end_y, &
         start_time, &
         end_time, &
         trans_time, &
         radius, &
         hfx_value
   end type line_type

   integer, parameter:: fire_max_lines = 5

   integer:: stat_lev = 1

   type lines_type
      type(line_type):: line(fire_max_lines)
      integer:: num_lines, &
                max_lines, &
                longlat
      real::  unit_fxlong, unit_fxlat
   end type lines_type

contains

!logical function isnan(a)
!real, intent(in):: a
!isnan= (a.ne.a)
!return
!end function isnan

   logical function isnotfinite(aa)
      real, intent(in)::aa
      isnotfinite = (aa .ne. aa .or. .not. aa .le. huge(aa) .or. .not. aa .ge. -huge(aa))
   end function isnotfinite

   subroutine interpolate_z2fire(id, &
                                 istrip, &
                                 ids, ide, jds, jde, &
                                 ims, ime, jms, jme, &
                                 ips, ipe, jps, jpe, &
                                 its, ite, jts, jte, &
                                 ifds, ifde, jfds, jfde, &
                                 ifms, ifme, jfms, jfme, &
                                 ifts, ifte, jfts, jfte, &
                                 ir, jr, &
                                 zs, &
                                 zsf)

      implicit none

      integer, intent(in)::id, &
                            istrip, &
                            ids, ide, jds, jde, &
                            ims, ime, jms, jme, &
                            ips, ipe, jps, jpe, &
                            its, ite, jts, jte, &
                            ifds, ifde, jfds, jfde, &
                            ifms, ifme, jfms, jfme, &
                            ifts, ifte, jfts, jfte, &
                            ir, jr
      real, intent(in), dimension(ims:ime, jms:jme):: zs
      real, intent(out), dimension(ifms:ifme, jfms:jfme):: &
         zsf

      real, dimension(its - 2:ite + 2, jts - 2:jte + 2):: za
      integer:: i, j, jts1, jte1, its1, ite1, jfts1, jfte1, ifts1, ifte1, itso, jtso, iteo, jteo

      !print *, 'LFR-DBG: calling crash'; call flush (6)
      if (istrip .gt. 1) call crash('interpolate_z2fire: istrip should be 0 or 1 or less')

      jts1 = max(jts - 1, jds)
      its1 = max(its - 1, ids)
      jte1 = min(jte + 1, jde)
      ite1 = min(ite + 1, ide)
      !print *, 'LFR-DBG: ',its1, ite1, jts1, jte1, size(za,1),size(za,2),size(zs,1),size(zs,2)
      do j = jts1, jte1
         do i = its1, ite1
         !print *, 'LFR-DBG: ',i,j,size(zs,1)<i,size(zs,2)<j
         !print *, 'LFR-DBG: ',isnan(zs(i, j))
            za(i, j) = zs(i, j)
         end do
      end do

      !print *, 'LFR-DBG: calling continue_at_boundary'; call flush (6)
      call continue_at_boundary(1, 1, 0., &
                                its - 2, ite + 2, jts - 2, jte + 2, &
                                ids, ide, jds, jde, &
                                ips, ipe, jps, jpe, &
                                its1, ite1, jts1, jte1, &
                                itso, jtso, iteo, jteo, &
                                za)

      jfts1 = snode(jfts, jfds, -istrip)
      ifts1 = snode(ifts, ifds, -istrip)
      jfte1 = snode(jfte, jfde, +istrip)
      ifte1 = snode(ifte, ifde, +istrip)

      !print *, 'LFR-DBG: calling interpolate_2d'; call flush (6)
      call interpolate_2d( &
         its - 2, ite + 2, jts - 2, jte + 2, &
         its1 - 1, ite1 + 1, jts1 - 1, jte1 + 1, &
         ifms, ifme, jfms, jfme, &
         ifts1, ifte1, jfts1, jfte1, &
         ir, jr, &
         real(ids), real(jds), ifds + (ir - 1)*0.5, jfds + (jr - 1)*0.5, &
         za, &
         zsf)

   end subroutine interpolate_z2fire

   subroutine crash(s)
!use module_wrf_error
      implicit none
      character(len=*), intent(in)::s
!character(len=128)msg
      character(len=128):: msg !INTRODUZIDO POR ISILDA
      msg = 'crash: '//s
      call message(msg, level=0)
      print *,'Fatal-Error: ',msg
      stop 'Crash wrf'
!$OMP CRITICAL(SFIRE_MESSAGE_CRIT)
!call wrf_error_fatal3("<stdin>",204,&
!msg)
!$OMP END CRITICAL(SFIRE_MESSAGE_CRIT)
   end subroutine crash

   subroutine warning(s, level)
      implicit none

      character(len=*), intent(in)::s
      character(len=128)::msg
      integer, intent(in), optional::level
      msg = 'WARNING:'//s
      if (present(level)) then
         call message(msg, level=level)
      else
         call message(msg, level=0)
      end if
   end subroutine warning

   subroutine message(s, level)
!use module_wrf_error
      implicit none

      character(len=*), intent(in)::s
      integer, intent(in), optional::level

      character(len=128)::msg
      character(len=118)::t
      integer m, mlevel
      logical op

      if (present(level)) then
         mlevel = level
      else
         mlevel = 2
      end if
      if (fire_print_msg .ge. mlevel) then
         m = 0
!$OMP CRITICAL(SFIRE_MESSAGE_CRIT)
         msg = 'SFIRE:'//s
         !call wrf_message(msg)

!$OMP END CRITICAL(SFIRE_MESSAGE_CRIT)
      end if
   end subroutine message

   subroutine time_start
!use module_timing, only:start_timing
      implicit none
!call start_timing
   end subroutine time_start

   subroutine time_end(string)
!use module_timing, only:end_timing
      implicit none
!character(len=*)string
      character(len=*)::string
!call end_timing(string)
   end subroutine time_end

   integer function open_text_file(filename, rw)
      implicit none
      character(len=*), intent(in):: filename, rw

      character(len=128):: msg
      character(len=1)::act
      integer::iounit, ierr
      logical::op
      character(LEN=1024) :: entire_line

      do iounit = 19, 99
         inquire (iounit, opened=op)
         if (.not. op) goto 1
      end do
      call crash('open_text_file: Cannot find any available I/O unit')
1     continue
      act = rw(1:1)
      select case (act)
      case ('r', 'R')
         open (iounit, FILE=trim(filename), FORM='FORMATTED', STATUS='OLD', ACTION='READ', IOSTAT=ierr)
         if (ierr .ne. 0) goto 998

         read (UNIT=iounit, FMT='(A)', IOSTAT=ierr) entire_line
         if (ierr .ne. 0) then
            call message('open_text_file: Cannot read from file '//trim(filename), level=0)
            goto 999
         end if
         backspace (UNIT=iounit)
      case ('w', 'W')
         open (iounit, FILE=trim(filename), FORM='FORMATTED', STATUS='UNKNOWN', ACTION='WRITE', IOSTAT=ierr)
         if (ierr .ne. 0) goto 998
      case default
         write (msg, *) 'open_text_file: bad mode ', trim(rw), ' for file ', trim(filename)
         call crash(msg)
      end select
      open_text_file = iounit
      return

998   call message('open_text_file: Cannot open file '//trim(filename), level=0)
999   write (msg, *) 'unit ', iounit, ' error ', ierr
      call crash(msg)

   end function open_text_file

   subroutine error_namelist(nml_read_unit)

      integer, intent(in)::nml_read_unit
      character(LEN=1024) :: entire_line
      character(LEN=128) :: msg
      backspace (UNIT=nml_read_unit)
      backspace (UNIT=nml_read_unit)
      read (UNIT=nml_read_unit, FMT='(A)') entire_line
      call message("Maybe here?: "//trim(entire_line), level=-1)
      read (UNIT=nml_read_unit, FMT='(A)') entire_line
      call message("Maybe here?: "//trim(entire_line), level=-1)
      write (msg, *) "Error in namelist, unit=", nml_read_unit
      call crash(msg)
   end

   subroutine set_ideal_coord(dxf, dyf, &
                              ifds, ifde, jfds, jfde, &
                              ifms, ifme, jfms, jfme, &
                              ifts, ifte, jfts, jfte, &
                              fxlong, fxlat &
                              )
      implicit none

      real, intent(in)::dxf, dyf
      integer, intent(in):: &
         ifds, ifde, jfds, jfde, &
         ifms, ifme, jfms, jfme, &
         ifts, ifte, jfts, jfte
      real, intent(out), dimension(ifms:ifme, jfms:jfme)::fxlong, fxlat

      integer::i, j

      do j = jfts, jfte
         do i = ifts, ifte

            fxlong(i, j) = (i - ifds + 0.5)*dxf
            fxlat(i, j) = (j - jfds + 0.5)*dyf
         end do
      end do
   end subroutine set_ideal_coord

   subroutine continue_at_boundary(ix, iy, bias, &
                                   ims, ime, jms, jme, &
                                   ids, ide, jds, jde, &
                                   ips, ipe, jps, jpe, &
                                   its, ite, jts, jte, &
                                   itso, iteo, jtso, jteo, &
                                   lfn)
      implicit none

      integer, intent(in)::ix, iy
      real, intent(in)::bias
      integer, intent(in)::ims, ime, jms, jme, &
                            ids, ide, jds, jde, &
                            ips, ipe, jps, jpe, &
                            its, ite, jts, jte
      integer, intent(out)::itso, jtso, iteo, jteo
      real, intent(inout), dimension(ims:ime, jms:jme)::lfn

      integer i, j
      character(len=128)::msg
      integer::its1, ite1, jts1, jte1
      integer, parameter::halo = 1

      call check_mesh_2dim(its - 1, ite + 1, jts - 1, jte + 1, ims, ime, jms, jme)

      itso = its
      jtso = jts
      iteo = ite
      jteo = jte

      its1 = its
      jts1 = jts
      ite1 = ite
      jte1 = jte
      if (its .eq. ips .and. .not. its .eq. ids) its1 = its - halo
      if (jts .eq. jps .and. .not. jts .eq. jds) jts1 = jts - halo
      if (ite .eq. ipe .and. .not. ite .eq. ide) ite1 = ite + halo
      if (jte .eq. jpe .and. .not. jte .eq. jde) jte1 = jte + halo
!$OMP CRITICAL(SFIRE_UTIL_CRIT)
      write (msg, '(a,2i5,a,f5.2)') 'continue_at_boundary: directions', ix, iy, ' bias ', bias
      call message(msg, level=3)
!$OMP END CRITICAL(SFIRE_UTIL_CRIT)
      if (ix .ne. 0) then
         if (its .eq. ids) then
            do j = jts1, jte1
               lfn(ids - 1, j) = EX(lfn(ids, j), lfn(ids + 1, j))
            end do
            itso = ids - 1
         end if
         if (ite .eq. ide) then
            do j = jts1, jte1
               lfn(ide + 1, j) = EX(lfn(ide, j), lfn(ide - 1, j))
            end do
            iteo = ide + 1
         end if
!$OMP CRITICAL(SFIRE_UTIL_CRIT)
         write (msg, '(8(a,i5))') 'continue_at_boundary: x:', its, ':', ite, ',', jts, ':', jte, ' ->', itso, ':', iteo, ',', jts1 &
               , ':', jte1
         call message(msg, level=3)
!$OMP END CRITICAL(SFIRE_UTIL_CRIT)
      end if
      if (iy .ne. 0) then
         if (jts .eq. jds) then
            do i = its1, ite1
               lfn(i, jds - 1) = EX(lfn(i, jds), lfn(i, jds + 1))
            end do
            jtso = jds - 1
         end if
         if (jte .eq. jde) then
            do i = its1, ite1
               lfn(i, jde + 1) = EX(lfn(i, jde), lfn(i, jde - 1))
            end do
            jteo = jde + 1
         end if
!$OMP CRITICAL(SFIRE_UTIL_CRIT) 
        write (msg, '(8(a,i5))') 'continue_at_boundary: y:', its, ':', ite, ',', jts, ':', jte, ' ->', its1, ':', ite1, ',', jtso &
                                 , ':', jteo
!$OMP END CRITICAL(SFIRE_UTIL_CRIT)
         call message(msg, level=3)
      end if

      if (ix .ne. 0 .and. iy .ne. 0) then
         if (its .eq. ids .and. jts .eq. jds) lfn(ids - 1, jds - 1) = EX(lfn(ids, jds), lfn(ids + 1, jds + 1))
         if (its .eq. ids .and. jte .eq. jde) lfn(ids - 1, jde + 1) = EX(lfn(ids, jde), lfn(ids + 1, jde - 1))
         if (ite .eq. ide .and. jts .eq. jds) lfn(ide + 1, jds - 1) = EX(lfn(ide, jds), lfn(ide - 1, jds + 1))
         if (ite .eq. ide .and. jte .eq. jde) lfn(ide + 1, jde + 1) = EX(lfn(ide, jde), lfn(ide - 1, jde - 1))
      end if
      return
   contains
      real function EX(a, b)

         real a, b
         EX = (1.-bias)*(2.*a - b) + bias*max(2.*a - b, a, b)
      end function EX
   end subroutine continue_at_boundary

   subroutine check_mesh_2dim(ids, ide, jds, jde, ims, ime, jms, jme)
      implicit none
      integer, intent(in)::ids, ide, jds, jde, ims, ime, jms, jme
!character(len=128)msg
      character(len=128)::msg
      if (ids < ims .or. ide > ime .or. jds < jms .or. jde > jme) then
!$OMP CRITICAL(SFIRE_UTIL_CRIT)
         write (msg, *) 'mesh dimensions:  ', ids, ide, jds, jde
         call message(msg, level=0)
         write (msg, *) 'memory dimensions:', ims, ime, jms, jme
!$OMP END CRITICAL(SFIRE_UTIL_CRIT)
         call message(msg, level=0)
         call crash('check_mesh_2dim: memory dimensions too small')
      end if
   end subroutine check_mesh_2dim

   subroutine check_mesh_3dim(ids, ide, kds, kde, jds, jde, ims, ime, kms, kme, jms, jme)
      integer, intent(in)::ids, ide, jds, jde, ims, ime, jms, jme, kds, kde, kms, kme
      if (ids < ims .or. ide > ime .or. jds < jms .or. jde > jme .or. kds < kms .or. kde > kme) then
         call crash('memory dimensions too small')
      end if
   end subroutine check_mesh_3dim

   subroutine sum_2d_cells( &
      ifms, ifme, jfms, jfme, &
      ifts, ifte, jtfs, jfte, &
      v2, &
      ims, ime, jms, jme, &
      its, ite, jts, jte, &
      v1)
      implicit none

      integer, intent(in)::its, ite, jts, jte, ims, ime, jms, jme
!real, intent(out)::v1(ims:ime,jms:jme)
      integer, intent(in)::ifts, ifte, jtfs, jfte, ifms, ifme, jfms, jfme
!real, intent(in)::v2(ifms:ifme,jfms:jfme)
      real, dimension(ims:ime, jms:jme), intent(out):: v1 ! INTRODUZIDO POR ISILDA CM
      real, dimension(ifms:ifme, jfms:jfme), intent(in):: v2 ! INTRODUZIDO POR ISILDA CM

      integer:: i, i_f, j, j_f, ir, jr, isz1, isz2, jsz1, jsz2, ioff, joff, ibase, jbase
      real :: t
!character(len=128)msg
      character(len=128):: msg !INTRODUZIDO POR ISILDA

      call check_mesh_2dim(its, ite, jts, jte, ims, ime, jms, jme)
      call check_mesh_2dim(ifts, ifte, jtfs, jfte, ifms, ifme, jfms, jfme)

      isz1 = ite - its + 1
      jsz1 = jte - jts + 1
      isz2 = ifte - ifts + 1
      jsz2 = jfte - jtfs + 1

      if (isz1 .le. 0 .or. jsz1 .le. 0 .or. isz2 .le. 0 .or. jsz2 .le. 0) then
         call message('all mesh sizes must be positive', level=0)
         goto 9
      end if

      ir = isz2/isz1
      jr = jsz2/jsz1

      if (isz2 .ne. isz1*ir .or. jsz2 .ne. jsz1*jr) then
         call message('input mesh size must be multiple of output mesh size', level=0)
         goto 9
      end if

      do j = jts, jte
         jbase = jtfs + jr*(j - jts)
         do i = its, ite
            ibase = ifts + ir*(i - its)
            t = 0.
            do joff = 0, jr - 1
               j_f = joff + jbase
               do ioff = 0, ir - 1
                  i_f = ioff + ibase
                  t = t + v2(i_f, j_f)
               end do
            end do
            v1(i, j) = t
         end do
      end do

      return

9     continue
!$OMP CRITICAL(SFIRE_UTIL_CRIT)
      write (msg, 91) ifts, ifte, jtfs, jfte, ifms, ifme, jfms, jfme
      call message(msg, level=0)
      write (msg, 91) its, ite, jts, jte, ims, ime, jms, jme
      call message(msg, level=0)
      write (msg, 92) 'input  mesh size:', isz2, jsz2
      call message(msg, level=0)
91    format('dimensions: ', 8i8)
      write (msg, 92) 'output mesh size:', isz1, jsz1
      call message(msg, level=0)
92    format(a, 2i8)
!$OMP END CRITICAL(SFIRE_UTIL_CRIT)
      call crash('sum_2d_cells: bad mesh sizes')

   end subroutine sum_2d_cells

   subroutine interpolate_2d( &
      ims2, ime2, jms2, jme2, &
      its2, ite2, jts2, jte2, &
      ims1, ime1, jms1, jme1, &
      its1, ite1, jts1, jte1, &
      ir, jr, &
      rip2, rjp2, rip1, rjp1, &
      v2, &
      v1)
      implicit none

      integer, intent(in)::its1, ite1, jts1, jte1, ims1, ime1, jms1, jme1
      integer, intent(in)::its2, ite2, jts2, jte2, ims2, ime2, jms2, jme2
      integer, intent(in)::ir, jr
      real, intent(in):: rjp1, rip1, rjp2, rip2
!real, intent(out)::v1(ims1:ime1,jms1:jme1)
!real, intent(in)::v2(ims2:ime2,jms2:jme2)
      real, dimension(ims1:ime1, jms1:jme1), intent(out):: v1 ! INTRODUZIDO POR ISILDA CM
      real, dimension(ims2:ime2, jms2:jme2), intent(in):: v2 ! INTRODUZIDO POR ISILDA CM

      integer:: i1, i2, j1, j2, is, ie, js, je
      real:: tx, ty, rx, ry
      real:: rio, rjo
      intrinsic::ceiling, floor

      call check_mesh_2dim(its1, ite1, jts1, jte1, ims1, ime1, jms1, jme1)
      call check_mesh_2dim(its2, ite2, jts2, jte2, ims2, ime2, jms2, jme2)

      rx = 1./ir
      ry = 1./jr

      do j2 = jts2, jte2 - 1
         rjo = rjp1 + jr*(j2 - rjp2)
         js = max(jts1, ceiling(rjo))
         je = min(jte1, floor(rjo) + jr)
         do i2 = its2, ite2 - 1
            rio = rip1 + ir*(i2 - rip2)
            is = max(its1, ceiling(rio))
            ie = min(ite1, floor(rio) + ir)
            do j1 = js, je
               ty = (j1 - rjo)*ry
               do i1 = is, ie

                  tx = (i1 - rio)*rx

                  v1(i1, j1) = &
                     (1 - tx)*(1 - ty)*v2(i2, j2) &
                     + (1 - tx)*ty*v2(i2, j2 + 1) &
                     + tx*(1 - ty)*v2(i2 + 1, j2) &
                     + tx*ty*v2(i2 + 1, j2 + 1)

               end do
            end do
         end do
      end do

   end subroutine interpolate_2d

   subroutine interpolate_2d_cells2cells( &
      ids2, ide2, jds2, jde2, ims2, ime2, jms2, jme2, v2, &
      ids1, ide1, jds1, jde1, ims1, ime1, jms1, jme1, v1)
      implicit none

      integer, intent(in)::ids1, ide1, jds1, jde1, ims1, ime1, jms1, jme1
!real, intent(out)::v1(ims1:ime1,jms1:jme1)
      integer, intent(in)::ids2, ide2, jds2, jde2, ims2, ime2, jms2, jme2
!real, intent(in)::v2(ims2:ime2,jms2:jme2)
      real, dimension(ims1:ime1, jms1:jme1), intent(out):: v1 ! INTRODUZIDO POR ISILDA CM
      real, dimension(ims2:ime2, jms2:jme2), intent(in):: v2 ! INTRODUZIDO POR ISILDA CM

      integer:: ir, jr, isz1, isz2, jsz1, jsz2, ip, jp, ih, jh
!character(len=128)msg
      character(len=128):: msg !INTRODUZIDO POR ISILDA

      call check_mesh_2dim(ids1, ide1, jds1, jde1, ims1, ime1, jms1, jme1)
      call check_mesh_2dim(ids2, ide2, jds2, jde2, ims2, ime2, jms2, jme2)

      isz1 = ide1 - ids1 + 1
      jsz1 = jde1 - jds1 + 1
      isz2 = ide2 - ids2 + 1
      jsz2 = jde2 - jds2 + 1

      if (isz1 .le. 0 .or. jsz1 .le. 0 .or. isz2 .le. 0 .or. jsz2 .le. 0) goto 9
      if (mod(isz1, isz2) .ne. 0 .or. mod(jsz1, jsz2) .ne. 0) goto 9

      ir = isz1/isz2
      jr = jsz1/jsz2

      ih = ir/2
      jh = jr/2

      ip = mod(ir + 1, 2)
      jp = mod(jr + 1, 2)

      call interpolate_2d_w(ip, jp, ih, jh, ir, jr, &
                            ids2, ide2, jds2, jde2, ims2, ime2, jms2, jme2, v2, &
                            ids1, ide1, jds1, jde1, ims1, ime1, jms1, jme1, v1)

      return

9     continue
!$OMP CRITICAL(SFIRE_UTIL_CRIT)
      write (msg, 91) ids2, ide2, jds2, jde2, ims2, ime2, jms2, jme2
      call message(msg, level=0)
      write (msg, 91) ids1, ide1, jds1, jde1, ims1, ime1, jms1, jme1
      call message(msg, level=0)
      write (msg, 92) 'input  mesh size:', isz2, jsz2
      call message(msg, level=0)
91    format('dimensions: ', 8i8)
      write (msg, 92) 'output mesh size:', isz1, jsz1
      call message(msg, level=0)
92    format(a, 2i8)
      call crash("module_fr_sfire_util:interpolate_2dmesh_cells: bad mesh sizes")
!$OMP END CRITICAL(SFIRE_UTIL_CRIT)
   end subroutine interpolate_2d_cells2cells

   subroutine interpolate_2d_cells2nodes( &
      ids2, ide2, jds2, jde2, ims2, ime2, jms2, jme2, v2, &
      ids1, ide1, jds1, jde1, ims1, ime1, jms1, jme1, v1)
      implicit none

      integer, intent(in)::ids1, ide1, jds1, jde1, ims1, ime1, jms1, jme1
!real, intent(out)::v1(ims1:ime1,jms1:jme1)
      integer, intent(in)::ids2, ide2, jds2, jde2, ims2, ime2, jms2, jme2
!real, intent(in)::v2(ims2:ime2,jms2:jme2)
      real, dimension(ims1:ime1, jms1:jme1), intent(out):: v1 ! INTRODUZIDO POR ISILDA CM
      real, dimension(ims2:ime2, jms2:jme2), intent(in):: v2 ! INTRODUZIDO POR ISILDA CM

      integer:: ir, jr, isz1, isz2, jsz1, jsz2, ip, jp, ih, jh
!character(len=128)msg

      character(len=128):: msg !INTRODUZIDO POR ISILDA

      call check_mesh_2dim(ids1, ide1 + 1, jds1, jde1 + 1, ims1, ime1, jms1, jme1)
      call check_mesh_2dim(ids2, ide2, jds2, jde2, ims2, ime2, jms2, jme2)

      isz1 = ide1 - ids1 + 1
      jsz1 = jde1 - jds1 + 1
      isz2 = ide2 - ids2 + 1
      jsz2 = jde2 - jds2 + 1

      if (isz1 .le. 0 .or. jsz1 .le. 0 .or. isz2 .le. 0 .or. jsz2 .le. 0) goto 9
      if (mod(isz1, isz2) .ne. 0 .or. mod(jsz1, jsz2) .ne. 0) goto 9

      ir = isz1/isz2
      jr = jsz1/jsz2

      ih = (ir + 1)/2
      jh = (jr + 1)/2

      ip = mod(ir, 2)
      jp = mod(jr, 2)

      call interpolate_2d_w(ip, jp, ih, jh, ir, jr, &
                            ids2, ide2, jds2, jde2, ims2, ime2, jms2, jme2, v2, &
                            ids1, ide1 + 1, jds1, jde1 + 1, ims1, ime1, jms1, jme1, v1)

      return
9     continue
!$OMP CRITICAL(SFIRE_UTIL_CRIT)
      write (msg, 91) ids2, ide2, jds2, jde2, ims2, ime2, jms2, jme2
      call message(msg, level=0)
      write (msg, 91) ids1, ide1, jds1, jde1, ims1, ime1, jms1, jme1
      call message(msg, level=0)
      write (msg, 92) 'input  mesh size:', isz2, jsz2
      call message(msg, level=0)
91    format('dimensions: ', 8i8)
      write (msg, 92) 'output mesh size:', isz1, jsz1
      call message(msg, level=0)
92    format(a, 2i8)
      call crash("module_fr_sfire_util:interpolate_2d_cells2nodes: bad mesh sizes")
!$OMP END CRITICAL(SFIRE_UTIL_CRIT)
   end subroutine interpolate_2d_cells2nodes

   subroutine interpolate_2d_w(ip, jp, ih, jh, ir, jr, &
                               ids2, ide2, jds2, jde2, ims2, ime2, jms2, jme2, v2, &
                               ids1, ide1, jds1, jde1, ims1, ime1, jms1, jme1, v1)
      implicit none

      integer, intent(in)::ip, jp, ih, jh, ir, jr
      integer, intent(in)::ids1, ide1, jds1, jde1, ims1, ime1, jms1, jme1
!real, intent(out)::v1(ims1:ime1,jms1:jme1)
      integer, intent(in)::ids2, ide2, jds2, jde2, ims2, ime2, jms2, jme2
!real, intent(in)::v2(ims2:ime2,jms2:jme2)
      real, dimension(ims1:ime1, jms1:jme1), intent(out):: v1 ! INTRODUZIDO POR ISILDA CM
      real, dimension(ims2:ime2, jms2:jme2), intent(in):: v2 ! INTRODUZIDO POR ISILDA CM
      real:: tx, ty, rx, ry, half, xoff, yoff
      integer:: i1, i2, j1, j2, ioff, joff
      parameter(half=0.5)

      rx = ir
      ry = jr

      xoff = ip*half
      yoff = jp*half

      do j2 = jds2, jde2 - 1
         do i2 = ids2, ide2 - 1
            do ioff = 0, ir - ip
               do joff = 0, jr - jp

                  i1 = ioff + (ih + ids1) + ir*(i2 - ids2)
                  j1 = joff + (jh + jds1) + jr*(j2 - jds2)

                  tx = (ioff + xoff)/rx
                  ty = (joff + yoff)/ry

                  v1(i1, j1) = &
                     (1 - tx)*(1 - ty)*v2(i2, j2) &
                     + (1 - tx)*ty*v2(i2, j2 + 1) &
                     + tx*(1 - ty)*v2(i2 + 1, j2) &
                     + tx*ty*v2(i2 + 1, j2 + 1)

               end do
            end do
         end do
      end do

      do ioff = 0, ih - 1
         do j2 = jds2, jde2 - 1
            do joff = 0, jr - jp
               j1 = joff + (jh + jds1) + jr*(j2 - jds2)

               ty = (joff + yoff)/ry

               v1(ids1 + ioff, j1) = (1 - ty)*v2(ids2, j2) + ty*v2(ids2, j2 + 1)
               v1(ide1 - ioff, j1) = (1 - ty)*v2(ide2, j2) + ty*v2(ide2, j2 + 1)
            end do
         end do
      end do
      do joff = 0, jh - 1
         do i2 = ids2, ide2 - 1
            do ioff = 0, ir - ip
               i1 = ioff + (ih + ids1) + ir*(i2 - ids2)

               tx = (ioff + xoff)/rx

               v1(i1, jds1 + joff) = (1 - tx)*v2(i2, jds2) + tx*v2(i2 + 1, jds2)
               v1(i1, jde1 - joff) = (1 - tx)*v2(i2, jde2) + tx*v2(i2 + 1, jde2)
            end do
         end do
      end do

      do ioff = 0, ih - 1
         do joff = 0, jh - 1
            v1(ids1 + ioff, jds1 + joff) = v2(ids2, jds2)
            v1(ide1 - ioff, jds1 + joff) = v2(ide2, jds2)
            v1(ids1 + ioff, jde1 - joff) = v2(ids2, jde2)
            v1(ide1 - ioff, jde1 - joff) = v2(ide2, jde2)
         end do
      end do
   end subroutine interpolate_2d_w

   real function interp(ids, ide, jds, jde, ims, ime, jms, jme, x, y, v)
      implicit none

      integer, intent(in)::ids, ide, jds, jde, ims, ime, jms, jme
      real, intent(in)::x, y
      real, dimension(ims:ime, jms:jme), intent(in):: v

      intrinsic floor, min, max

      integer:: i, j
      real:: tx, ty

      i = floor(x)
      i = max(min(i, ide), ids)
      j = floor(y)
      j = max(min(j, jde), jds)

      tx = x - real(i)
      ty = y - real(j)

      interp = &
         (1 - tx)*(1 - ty)*v(i, j) &
         + tx*(1 - ty)*v(i + 1, j) &
         + (1 - tx)*ty*v(i, j + 1) &
         + tx*ty*v(i + 1, j + 1)

   end function interp

   subroutine meshdiffc_2d(ids, ide, jds, jde, &
                           ims1, ime1, jms1, jme1, &
                           dx, dy, &
                           lfn, &
                           diffCx, diffCy)
      implicit none

      integer, intent(in)::ids, ide, jds, jde, ims1, ime1, jms1, jme1
      real, intent(in):: dx, dy
      real, intent(in), dimension(ims1:ime1, jms1:jme1):: lfn
      real, intent(out), dimension(ims1:ime1, jms1:jme1):: diffCx, diffCy

      integer:: i, j
      real, dimension(ims1:ime1, jms1:jme1):: diffLx, diffRx, diffLy, diffRy

      call meshdiff_2d(ids, ide, jds, jde, &
                       ims1, ime1, jms1, jme1, &
                       dx, dy, &
                       lfn, &
                       diffLx, diffRx, diffLy, diffRy)

      do j = jds, jde + 1
         do i = ids, ide + 1
            diffCx(i, j) = 0.5*(diffLx(i, j) + diffRx(i, j))
            diffCy(i, j) = 0.5*(diffLy(i, j) + diffRy(i, j))
         end do
      end do
   end subroutine meshdiffc_2d

   subroutine meshdiff_2d(ids, ide, jds, jde, &
                          ims1, ime1, jms1, jme1, &
                          dx, dy, &
                          lfn, &
                          diffLx, diffRx, diffLy, diffRy)
      implicit none

      integer, intent(in)::ids, ide, jds, jde, ims1, ime1, jms1, jme1
      real, intent(in):: dx, dy
      real, intent(in), dimension(ims1:ime1, jms1:jme1):: lfn
      real, intent(out), dimension(ims1:ime1, jms1:jme1):: diffLx, diffRx, diffLy, diffRy

      integer:: i, j
      real:: tmpx, tmpy

      call check_mesh_2dim(ids, ide + 1, jds, jde + 1, ims1, ime1, jms1, jme1)

      do j = jds, jde
         do i = ids, ide
            tmpx = (lfn(i + 1, j) - lfn(i, j))/dx
            diffLx(i + 1, j) = tmpx
            diffRx(i, j) = tmpx
            tmpy = (lfn(i, j + 1) - lfn(i, j))/dy
            diffLy(i, j + 1) = tmpy
            diffRy(i, j) = tmpy
         end do

         diffLx(ids, j) = diffLx(ids + 1, j)
         diffRx(ide + 1, j) = diffRx(ide, j)
      end do

      do i = ids, ide
         tmpx = (lfn(i + 1, j) - lfn(i, j))/dx
         diffLx(i + 1, j) = tmpx
         diffRx(i, j) = tmpx
      end do

      do j = jds, jde
         tmpy = (lfn(i, j + 1) - lfn(i, j))/dy
         diffLy(i, j + 1) = tmpy
         diffRy(i, j) = tmpy
      end do

      diffLx(ids, j) = diffLx(ids + 1, j)
      diffRx(ide + 1, j) = diffRx(ide, j)
      do i = ids, ide + 1
         diffLy(i, jds) = diffLy(i, jds + 1)
         diffRy(i, jde + 1) = diffRy(i, jde)
      end do

   end subroutine meshdiff_2d

   real pure function sum_2darray(its, ite, jts, jte, &
                                  ims, ime, jms, jme, &
                                  a) result(saida)
      integer, intent(in)::its, ite, jts, jte, ims, ime, jms, jme
      real, intent(in)::a(ims:ime, jms:jme)

      integer:: i, j
      real:: t

      t = 0.
      do j = jts, jte
         do i = its, ite
            t = t + a(i, j)
         end do
      end do
!sum_2darray = t
      saida = t
   end function sum_2darray

   real pure function max_2darray(its, ite, jts, jte, &
                                  ims, ime, jms, jme, &
                                  a) result(saida)
      integer, intent(in)::its, ite, jts, jte, ims, ime, jms, jme
      real, intent(in)::a(ims:ime, jms:jme)

      integer:: i, j
      real:: T

      t = 0.
      do j = jts, jte
         do i = its, ite
            t = max(t, a(i, j))
         end do
      end do
!max_2darray = t
      saida = t
   end function max_2darray

   subroutine print_2d_stats_vec(ips, ipe, jps, jpe, &
                                 ims, ime, jms, jme, &
                                 ax, ay, name)
      implicit none
      integer, intent(in)::ips, ipe, jps, jpe, ims, ime, jms, jme
      real, intent(in), dimension(ims:ime, jms:jme)::ax, ay
      character(len=*), intent(in)::name
      integer:: i, j
      real:: t
      real:: avg_a, max_a, min_a
      character(len=25)::id
      id = name
      call print_2d_stats(ips, ipe, jps, jpe, &
                          ims, ime, jms, jme, &
                          ax, id//'/x ')
      call print_2d_stats(ips, ipe, jps, jpe, &
                          ims, ime, jms, jme, &
                          ay, id//'/y ')
      avg_a = 0
      max_a = -huge(max_a)
      min_a = huge(min_a)
      do j = jps, jpe
         do i = ips, ipe
            ! print*,'THKM',j,i,jpe,ipe,ax(i,j),ay(i,j)
            ! call flush(6)
            t = sqrt(ax(i, j)**2 + ay(i, j)**2)
            !print*,'TTTTTTTA',t
            !call flush(6)
            max_a = max(max_a, t)
            min_a = min(min_a, t)
            avg_a = avg_a + t
         end do
      end do
      avg_a = avg_a/((ipe - ips + 1)*(jpe - jps + 1))
!print*,'YYYavg_a',avg_a
!call flush(6)
      call print_stat_line(id//'/sz', ips, ipe, jps, jpe, min_a, max_a, avg_a)
   end subroutine print_2d_stats_vec

   subroutine print_stat_line(name, ips, ipe, jps, jpe, min_a, max_a, avg_a)

      implicit none

      integer, intent(in)::ips, ipe, jps, jpe
      character(len=*), intent(in)::name
      real, intent(in)::min_a, max_a, avg_a

      character(len=128)::msg
      character(len=35)::id

      if (.not. avg_a .eq. avg_a) then
         !print*,'print_stat_line','NaN'
         msg = 'NaN detected in '//trim(name)
         call crash(msg)
      end if
      if (fire_print_msg .eq. 0) return
      id = name
!$OMP CRITICAL(SFIRE_UTIL_CRIT)
      write (msg, '(a,4i5,3g13.5)') id, ips, ipe, jps, jpe, min_a, max_a, avg_a
!$OMP END CRITICAL(SFIRE_UTIL_CRIT)
      call message(msg, level=2)
   end subroutine print_stat_line

   subroutine print_3d_stats_by_slice(ips, ipe, kps, kpe, jps, jpe, &
                                      ims, ime, kms, kme, jms, jme, &
                                      a, name)
      implicit none
      integer, intent(in)::ips, ipe, jps, jpe, ims, ime, jms, jme, kms, kme, kps, kpe
      real, intent(in)::a(ims:ime, kms:kme, jms:jme)
      character(len=*), intent(in)::name
      integer::k
      character(len=128)::msg
      do k = kps, kpe

!$OMP CRITICAL(SFIRE_UTIL_CRIT)
         write (msg, '(i2,1x,a)') k, name
!$OMP END CRITICAL(SFIRE_UTIL_CRIT)
         call print_3d_stats(ips, ipe, k, k, jps, jpe, &
                             ims, ime, kms, kme, jms, jme, &
                             a, msg)
      end do
!print*,"estou a sair da rotina print_3d_stats_by_slice"
   end subroutine print_3d_stats_by_slice

   subroutine print_3d_stats(ips, ipe, kps, kpe, jps, jpe, &
                             ims, ime, kms, kme, jms, jme, &
                             a, name)
      implicit none
      integer, intent(in)::ips, ipe, jps, jpe, ims, ime, jms, jme, kms, kme, kps, kpe
      real, intent(in)::a(ims:ime, kms:kme, jms:jme)
      character(len=*), intent(in)::name
      integer:: i, j, k
      real:: avg_a, max_a, min_a, t, aa, bb
      character(len=128)::msg

!print*,"estou na rotina print_3d_stats_by_slice dentro da print_3d_stats"
      !print *, 'LFR-DBG para o nome enviado "', trim(name), '"'
      !print *, 'LFR-DBG min e max:', minval(a), maxval(a); call flush (6)
      bb = 0.
      do j = jps, jpe
         do k = kps, kpe
            do i = ips, ipe
               bb = bb + a(i, k, j)
            end do
         end do
      end do
!print*,'BBBB=',bb
!call flush(6)
      if (bb .eq. bb .and. fire_print_msg .eq. 0) return
      avg_a = 0.
      max_a = -huge(max_a)
      min_a = huge(min_a)
      t = huge(t)
      do j = jps, jpe
         do k = kps, kpe
            do i = ips, ipe
               aa = a(i, k, j)
               if (aa .ne. aa .or. .not. aa .le. t .or. .not. aa .ge. -t) goto 9
               max_a = max(max_a, aa)
               min_a = min(min_a, aa)
               avg_a = avg_a + aa
            end do
         end do
      end do
      if (bb .ne. bb) goto 10
      if (fire_print_msg .le. 0) return
      avg_a = avg_a/((ipe - ips + 1)*(jpe - jps + 1)*(kpe - kps + 1))
!print*,'GGGGGGG',avg_a
!call flush(6)
!print*,"vou levar o valor para print_stat_line"
!call flush(6)
      call print_stat_line(name, ips, ipe, jps, jpe, min_a, max_a, avg_a)
!print*,"estou fora da print_stat_line"
!call flush(6)
      return
9     continue
!$OMP CRITICAL(SFIRE_UTIL_CRIT)
      write (msg, 1) name, i, k, j, aa
      call message(msg, level=0)
1     format(a30, '(', i6, ',', i6, ',', i6, ') = ', g13.5)
      write (msg, 2) 'patch dimensions ', ips, ipe, kps, kpe, jps, jpe
      call message(msg, level=0)
      write (msg, 2) 'memory dimensions', ims, ime, kms, kme, jms, jme
      call message(msg, level=0)
2     format(a, 6i8)
!$OMP END CRITICAL(SFIRE_UTIL_CRIT)
      call print_stat_line(name, ips, ipe, jps, jpe, aa, aa, aa)
      if (aa .ne. aa) goto 10
      msg = 'Invalid floating point number detected in '//name
      call crash(msg)
10    msg = 'NaN detected in '//name
      call crash(msg)
!print*,"estou na rotina print_3d_stats_by_slice vou sair da rotina print_3d_stats"
   end subroutine print_3d_stats

   subroutine print_2d_stats(ips, ipe, jps, jpe, &
                             ims, ime, jms, jme, &
                             a, name)
      implicit none
      integer, intent(in)::ips, ipe, jps, jpe, ims, ime, jms, jme
      real, intent(in)::a(ims:ime, jms:jme)
      character(len=*), intent(in)::name

      call print_3d_stats(ips, ipe, 1, 1, jps, jpe, &
                          ims, ime, 1, 1, jms, jme, &
                          a, name)

   end subroutine print_2d_stats

   real pure function avg_2darray(its, ite, jts, jte, &
                                  ims, ime, jms, jme, &
                                  a) result(saida)
      integer, intent(in)::its, ite, jts, jte, ims, ime, jms, jme
      real, intent(in)::a(ims:ime, jms:jme)

!avg_2darray = sum_2darray( its,ite,jts,jte,               &
!                           ims,ime,jms,jme,               &
!                           a)/((ite-its+1)*(jte-jts+1))

      saida = sum_2darray(its, ite, jts, jte, &
                          ims, ime, jms, jme, &
                          a)/((ite - its + 1)*(jte - jts + 1))

   end function avg_2darray

   real pure function avg_2darray_vec(its, ite, jts, jte, &
                                      ims, ime, jms, jme, &
                                      ax, ay) result(saida)
      integer, intent(in)::its, ite, jts, jte, ims, ime, jms, jme
      real, intent(in), dimension(ims:ime, jms:jme):: ax, ay

      integer:: i, j
      real:: t
      t = 0.
      do j = jts, jte
         do i = its, ite
            t = t + sqrt(ax(i, j)**2 + ay(i, j)**2)
         end do
      end do
      t = t/((ite - its + 1)*(jte - jts + 1))
!avg_2darray_vec = t
      saida = t
   end function avg_2darray_vec

   subroutine print_array(its, ite, jts, jte, &
                          ims, ime, jms, jme, &
                          a, name, id)

      integer, intent(in)::its, ite, jts, jte, ims, ime, jms, jme, id
      real, intent(in), dimension(ims:ime, jms:jme):: a
      character(len=*), intent(in)::name

      integer:: i, j
      character(len=128)::msg

!$OMP CRITICAL(SFIRE_UTIL_CRIT)
      write (msg, *) name, ' start ', id, ' dim ', its, ite, jts, jte
      call message(msg)
      do j = jts, jte
         do i = its, ite
            write (msg, *) i, j, a(i, j)
            call message(msg)
         end do
      end do
      write (msg, *) name, ' end ', id
      call message(msg)
!$OMP END CRITICAL(SFIRE_UTIL_CRIT)
   end subroutine print_array

   subroutine write_array_m(its, ite, jts, jte, &
                            ims, ime, jms, jme, &
                            a, name, id)

      integer, intent(in)::its, ite, jts, jte, ims, ime, jms, jme, id
      real, intent(in), dimension(ims:ime, jms:jme):: a
      character(len=*), intent(in)::name
      !print *, "estou na subroutine write_array_m vou chamar a write_array_m3"
      call write_array_m3(its, ite, 1, 1, jts, jte, &
                          ims, ime, 1, 1, jms, jme, &
                          a, name, id)
   end subroutine write_array_m

   subroutine write_array_m3(its, ite, kts, kte, jts, jte, &
                             ims, ime, kms, kme, jms, jme, &
                             a, name, id)
!use module_dm

      implicit none

      integer, intent(in)::its, ite, jts, jte, ims, ime, jms, jme, kts, kte, kms, kme, id
      real, intent(in), dimension(ims:ime, kms:kme, jms:jme):: a
      character(len=*), intent(in)::name

      integer:: i, j, k, iu, ilen, myproc, nprocs
      logical:: op
      character(len=128)::fname, msg
      !print *, "estou na subroutine write_array_m3"
      if (fire_print_file .eq. 0 .or. id .le. 0) return
      call check_mesh_2dim(its, ite, jts, jte, ims, ime, jms, jme)
!call wrf_get_nproc (nprocs)
!call wrf_get_myproc(myproc)
      !print *, "estou na subroutine write_array_m3 vou escrever o nome do ficheiro", name
!#####$OMP CRITICAL(SFIRE_UTIL_CRIT)
!COMENTADO POR ISILDA CUNHA MENEZES
!if(nprocs.eq.1)then
      write (fname, 3) name, '_', id, '.txt'
!else
!    write(fname,4)name,'_',id,'.',myproc,'.txt'
!endif
!######
      iu = 0
      do i = 6, 99
         inquire (unit=i, opened=op)
         if (.not. op .and. iu .le. 0) iu = i
      end do
      if (iu .gt. 0) open (iu, file=trim(fname), form='formatted', status='unknown')

      if (iu .le. 0) call crash('write_array_m: cannot find available fortran unit')
!print*,"estou na subroutine write_array_m3 vou escrever a variavel", name
      write (iu, 1) real(its)
      write (iu, 1) real(ite)
      write (iu, 1) real(jts)
      write (iu, 1) real(jte)
      write (iu, 1) real(kts)
      write (iu, 1) real(kte)
      write (iu, 1) (((a(i, k, j), i=its, ite), j=jts, jte), k=kts, kte)
      close (iu)
      write (msg, 2) name, '(', its, ':', ite, ',', jts, ':', jte, ',', &
         kts, ':', kte, ') -> ', trim(fname)
!print*,"vou sair da rotina write_array_m3"
!$OMP END CRITICAL(SFIRE_UTIL_CRIT)
      call message(msg)
      return

1     format(e20.12)
2     format(2a, 3(i5, a, i5, a), 2a)
3     format(a, a, i8.8, a)
4     format(a, a, i8.8, a, i4.4, a)

   end subroutine write_array_m3

   subroutine read_array_2d_real(filename, a, its, ite, jts, jte, ims, ime, jms, jme)
!use module_dm
      implicit none

      integer, intent(in)::its, ite, jts, jte, ims, ime, jms, jme
      real, intent(out), dimension(ims:ime, jms:jme):: a
      character(len=*), intent(in)::filename

      integer:: i, j, ni, nj, mi, mj, nprocs, myproc, mythread, iu
      logical:: op
      character(len=128)::fname, msg

!call wrf_get_nproc (nprocs)
!call wrf_get_myproc( myproc )
      mythread = 0
      if (nprocs .ne. 1 .or. myproc .ne. 0 .or. mythread .ne. 0) &
         call crash('read_array_2d: parallel execution not supported')

      mi = ite - its + 1
      mj = jte - jts + 1
      write (msg, 2) 'reading array size ', mi, mj, ' from file ', trim(filename)
2     format(a, 2i6, 2a)
      call message(msg, level=1)

      call check_mesh_2dim(its, ite, jts, jte, ims, ime, jms, jme)

      iu = 0
      do i = 11, 99
         inquire (unit=i, opened=op)
         if (.not. op .and. iu .le. 0) iu = i
      end do
      if (iu .le. 0) call crash('read_array_2d: cannot find available fortran unit')

      if (iu .gt. 0) open (iu, file=filename, form='formatted', status='old', err=9)
      rewind (iu, err=9)

      read (iu, *, err=10) ni, nj
      if (ni .ne. mi .or. nj .ne. mj) then
         write (msg, '(a,2i6,a,2i6)') 'Array dimensions', ni, nj, ' in the input file should be ', mi, mj
         call message(msg, level=0)
         goto 10
      end if
      do i = its, ite
         read (iu, *, err=10) (a(i, j), j=jts, jte)
      end do
      close (iu, err=11)
      call print_2d_stats(its, ite, jts, jte, &
                          ims, ime, jms, jme, &
                          a, filename)
      write (6, *) its, jts, a(its, jts), loc(a(its, jts))
      return

9     msg = 'Error opening file '//trim(filename)
      call crash(msg)
10    msg = 'Error reading file '//trim(filename)
      call crash(msg)
11    msg = 'Error closing file '//trim(filename)
      call crash(msg)
   end subroutine read_array_2d_real

   pure integer function ifval(l, i, j) result(saida)
      implicit none
      logical, intent(in)::l
      integer, intent(in)::i, j
      if (l) then
         !ifval=i
         saida = i
      else
         !ifval=j
         saida = j
      end if

   end function ifval

   pure integer function snode(t, d, i) result(saida)
      implicit none
      integer, intent(in)::t, d, i
      if (t .ne. d) then
         ! snode=t
         saida = t
      else
         ! snode=t+i
         saida = t
      end if
   end function snode

   subroutine print_chsum(id, &
                          ims, ime, kms, kme, jms, jme, &
                          ids, ide, kds, kde, jds, jde, &
                          ips, ipe, kps, kpe, jps, jpe, &
                          istag, kstag, jstag, &
                          a, name)

      !USE module_dm , only : wrf_dm_bxor_integer

      integer, intent(in):: id, &
                            ims, ime, kms, kme, jms, jme, &
                            ids, ide, kds, kde, jds, jde, &
                            ips, ipe, kps, kpe, jps, jpe, &
                            istag, kstag, jstag
      real, intent(in), dimension(ims:ime, kms:kme, jms:jme)::a
      character(len=*)::name

      integer::lsum
      integer::i, j, k, n, ipe1, jpe1, kpe1, iel, thread, is, js, ks
      integer, save::psum, gsum
      real::rel
      equivalence(rel, iel)
!character(len=256)msg
      character(len=256):: msg

      if (fire_print_msg .le. 0) return
      
      ipe1 = ifval(ipe .eq. ide .and. istag .ne. 0, ipe + 1, ipe)
      kpe1 = ifval(kpe .eq. kde .and. kstag .ne. 0, kpe + 1, kpe)
      jpe1 = ifval(jpe .eq. jde .and. jstag .ne. 0, jpe + 1, jpe)
      is = ifval(istag .ne. 0, 1, 0)
      ks = ifval(kstag .ne. 0, 1, 0)
      js = ifval(jstag .ne. 0, 1, 0)

if(size(a,3)<jpe1) write (*,fmt='(A,6(I3.3,1X))') 'LFR-DBG: ,jpe,jde,jstag,jpe1,jme,jms:',jpe,jde,jstag,jpe1,jme,jms

      lsum = 0
      do j = jps, jpe1
         do k = kps, kpe1 - 1
            do i = ips, ipe1
               rel = a(i, k, j)
               lsum = ieor(lsum, iel)
            end do
         end do
      end do

      thread = 0

      if (thread .eq. 0) psum = 0
!$OMP BARRIER
!$OMP CRITICAL(CHSUM)
      psum = ieor(psum, lsum)
!$OMP END CRITICAL(CHSUM)
!$OMP BARRIER

      if (thread .eq. 0) then
         !gsum = wrf_dm_bxor_integer ( psum )
         gsum = psum
         write (msg, 1) id, name, ids, ide + is, kds, kde + ks, jds, jde + js, gsum
1        format(i6, 1x, a10, ' dims', 6i5, ' chsum ', z8.8)
         call message(msg)
      end if

   end subroutine print_chsum

   real function fun_real(fun, &
                          ims, ime, kms, kme, jms, jme, &
                          ids, ide, kds, kde, jds, jde, &
                          ips, ipe, kps, kpe, jps, jpe, &
                          istag, kstag, jstag, &
                          a, b)

      !USE module_dm , only : wrf_dm_sum_real , wrf_dm_max_real

      integer, intent(in)::  fun, &
                            ims, ime, kms, kme, jms, jme, &
                            ids, ide, kds, kde, jds, jde, &
                            ips, ipe, kps, kpe, jps, jpe, &
                            istag, kstag, jstag
      real, intent(in), dimension(ims:ime, kms:kme, jms:jme)::a, b

      real::lsum, void
      integer::i, j, k, n, ipe1, jpe1, kpe1, iel, thread, is, js, ks
      real, save::psum, gsum
      real::rel
      logical:: dosum, domax, domin
      character(len=256):: msg

      ipe1 = ifval(ipe .eq. ide .and. istag .ne. 0, ipe + 1, ipe)
      kpe1 = ifval(kpe .eq. kde .and. kstag .ne. 0, kpe + 1, kpe)
      jpe1 = ifval(jpe .eq. jde .and. jstag .ne. 0, jpe + 1, jpe)
      is = ifval(istag .ne. 0, 1, 0)
      ks = ifval(kstag .ne. 0, 1, 0)
      js = ifval(jstag .ne. 0, 1, 0)

      if (fun .eq. REAL_SUM) then
         void = 0.
         lsum = void
         do j = jps, jpe1
            do k = kps, kpe1
               do i = ips, ipe1
                  lsum = lsum + a(i, k, j)
               end do
            end do
         end do
      elseif (fun .eq. RNRM_SUM) then
         void = 0.
         lsum = void
         do j = jps, jpe1
            do k = kps, kpe1
               do i = ips, ipe1
                  lsum = lsum + sqrt(a(i, k, j)*a(i, k, j) + b(i, k, j)*b(i, k, j))
               end do
            end do
         end do
      elseif (fun .eq. REAL_MAX) then
         void = -huge(lsum)
         lsum = void
         do j = jps, jpe1
            do k = kps, kpe1
               do i = ips, ipe1
                  lsum = max(lsum, a(i, k, j))
               end do
            end do
         end do
      elseif (fun .eq. REAL_AMAX) then
         void = -huge(lsum)
         lsum = void
         do j = jps, jpe1
            do k = kps, kpe1
               do i = ips, ipe1
                  lsum = max(lsum, abs(a(i, k, j)))
               end do
            end do
         end do
      elseif (fun .eq. REAL_MIN) then
         void = huge(lsum)
         lsum = void
         do j = jps, jpe1
            do k = kps, kpe1
               do i = ips, ipe1
                  lsum = min(lsum, a(i, k, j))
               end do
            end do
         end do
      elseif (fun .eq. RNRM_MAX) then
         void = 0.
         lsum = void
         do j = jps, jpe1
            do k = kps, kpe1
               do i = ips, ipe1
                  lsum = max(lsum, sqrt(a(i, k, j)*a(i, k, j) + b(i, k, j)*b(i, k, j)))
               end do
            end do
         end do
      else
         call crash('fun_real: bad fun')
      end if

      if (lsum .ne. lsum) call message('fun_real: WARNING: NaN detected')

      dosum = fun .eq. REAL_SUM .or. fun .eq. RNRM_SUM
      domax = fun .eq. REAL_MAX .or. fun .eq. REAL_AMAX .or. fun .eq. RNRM_MAX
      domin = fun .eq. REAL_MIN

!$OMP SINGLE

      psum = void
!$OMP END SINGLE
!$OMP BARRIER

!$OMP CRITICAL(RDSUM)

      if (dosum) psum = psum + lsum
      if (domax) psum = max(psum, lsum)
      if (domin) psum = min(psum, lsum)
!$OMP END CRITICAL(RDSUM)

!$OMP BARRIER

!$OMP SINGLE

      !if(dosum) gsum = wrf_dm_sum_real ( psum )
      !if(domax) gsum = wrf_dm_max_real ( psum )
      if (gsum .ne. gsum) call message('fun_real: WARNING: NaN detected')
!$OMP END SINGLE

!$OMP BARRIER

      fun_real = gsum

   end function fun_real

   subroutine sfire_debug_hook(fire_debug_hook_sec)
      integer, intent(in)::fire_debug_hook_sec
      integer, save:: go = -1
      external:: wrf_dm_bcast_integer
      if (go < 0) then
         go = fire_debug_hook_sec
      end if
!do while (go .ne. 0)
!    call sleep(go)

      !call wrf_dm_bcast_integer(abs(go),1)
!enddo
   end subroutine sfire_debug_hook

end module module_fr_sfire_util

