
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
                                 ifds, ifde, jfds, jfde, &
                                 ifms, ifme, jfms, jfme, &
                                 ir, jr, &
                                 zs, &
                                 zsf)

      implicit none

      integer, intent(in)::id, &
                            istrip, &
                            ids, ide, jds, jde, &
                            ims, ime, jms, jme, &
                            ifds, ifde, jfds, jfde, &
                            ifms, ifme, jfms, jfme, &
                            ir, jr
      real, intent(in), dimension(ims:ime, jms:jme):: zs
      real, intent(out), dimension(ifms:ifme, jfms:jfme):: &
         zsf
      !##########PEREIRA
    !  real, dimension(ids - 2:ide + 2, jds - 2:jde + 2):: za
      real, dimension(ims:ime, jms:jme):: za
      !##################
      integer:: i, j, jts1, jte1, its1, ite1, jfts1, jfte1, ifts1, ifte1, itso, jtso, iteo, jteo
      print *, 'LFR-DBG: "ESTOU DENTRO ROTINA interpolate_z2fire'
      !print *, 'LFR-DBG: calling crash'; call flush (6)
      if (istrip .gt. 1) call crash('interpolate_z2fire: istrip should be 0 or 1 or less')

      jts1 = jds
      its1 = ids
      jte1 = jde
      ite1 = ide
      print *, 'LFR-DBG: ',its1, ite1, jts1, jte1, size(za,1),size(za,2),size(zs,1),size(zs,2)
      do j = jts1, jte1
         do i = its1, ite1
         !print *, 'LFR-DBG: ',i,j,size(zs,1)<i,size(zs,2)<j
         !print *, 'LFR-DBG: ',isnan(zs(i, j))
            za(i, j) = zs(i, j)
         end do
      end do

      print *, 'LFR-DBG: calling continue_at_boundary'; call flush (6)
      !####PEREIRA
     ! call continue_at_boundary(1, 1, 0., &
     !                           ids-2, ide+2, jds-2, jde+2, &
     !                           ids, ide, jds, jde, &
     !                           itso, jtso, iteo, jteo, &
     !                           za)
      call continue_at_boundary(1, 1, 0., &
                                ims, ime, jms, jme, &
                                ids, ide, jds, jde, &
                                itso, jtso, iteo, jteo, &
                                za)
     !####################################
      jfts1 = jfds
      ifts1 = ifds
      jfte1 = jfde
      ifte1 = ifde

      print *, 'LFR-DBG: calling interpolate_2d'; call flush (6)
     !###### PEREIRA 
    !  call interpolate_2d( &
    !     ids - 2, ide + 2, jds - 2, jde + 2, &
    !     its1 - 1, ite1 + 1, jts1 - 1, jte1 + 1, &
    !     ifms, ifme, jfms, jfme, &
    !     ifts1, ifte1, jfts1, jfte1, &
    !     ir, jr, &
    !     real(ids), real(jds), ifds + (ir - 1)*0.5, jfds + (jr - 1)*0.5, &
    !     za, &
    !     zsf)
 
       call interpolate_2d( &
         ims, ime, jms, jme, &
         its1, ite1, jts1, jte1, &
         ifms, ifme, jfms, jfme, &
         ifts1, ifte1, jfts1, jfte1, &
         ir, jr, &
         real(ids), real(jds), ifds + (ir - 1)*0.5, jfds + (jr - 1)*0.5, &
         za, &
         zsf)
     !##############################

   print *, 'LFR-DBG: "ESTOU A SAIR DA ROTINA interpolate_z2fire Ja fiz os calculos'
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

 !  subroutine set_ideal_coord(dxf, dyf, &
 !                             ifds, ifde, jfds, jfde, &
 !                             ifms, ifme, jfms, jfme, &
 !                             fxlong, fxlat &
 !                             )
 !     implicit none
 !
 !     real, intent(in)::dxf, dyf
 !     integer, intent(in):: &
 !        ifds, ifde, jfds, jfde, &
 !        ifms, ifme, jfms, jfme, &
 !     real, intent(out), dimension(ifms:ifme, jfms:jfme)::fxlong, fxlat
 !
 !     integer::i, j
 !
 !     do j = jfds, jfde
 !        do i = ifds, ifde
 !
 !           fxlong(i, j) = (i - ifds + 0.5)*dxf
 !           fxlat(i, j) = (j - jfds + 0.5)*dyf
 !        end do
 !     end do
 !  end subroutine set_ideal_coord

   subroutine continue_at_boundary(ix, iy, bias, &
                                   ims, ime, jms, jme, &
                                   ids, ide, jds, jde, &
                                   itso, iteo, jtso, jteo, &
                                   lfn)
      implicit none

      integer, intent(in)::ix, iy
      real, intent(in)::bias
      integer, intent(in)::ims, ime, jms, jme, &
                           ids, ide, jds, jde
                           
      integer, intent(out)::itso, jtso, iteo, jteo
      real, intent(inout), dimension(ims:ime, jms:jme)::lfn

      integer i, j
      character(len=128)::msg
      integer:: its1, ite1, jte1, jts1 
      
      call check_mesh_2dim(ids - 1, ide + 1, jds - 1, jde + 1, ims, ime, jms, jme)

      itso = ids
      jtso = jds
      iteo = ide
      jteo = jde

      its1 = ids
      jts1 = jds
      ite1 = ide
      jte1 = jde


      write (msg, '(a,2i5,a,f5.2)') 'continue_at_boundary: directions', ix, iy, ' bias ', bias
      call message(msg, level=3)

      if (ix .ne. 0) then

            do j = jts1, jte1
               lfn(ids - 1, j) = EX(lfn(ids, j), lfn(ids + 1, j))
            end do
            itso = ids - 1

     
            do j = jts1, jte1
               lfn(ide + 1, j) = EX(lfn(ide, j), lfn(ide - 1, j))
            end do
            iteo = ide + 1
     

         write (msg, '(8(a,i5))') 'continue_at_boundary: x:', ids, ':', ide, ',', jds, ':', jde, ' ->', itso, ':', iteo, ',', jts1 &
               , ':', jte1
         call message(msg, level=3)

      end if
      
      if (iy .ne. 0) then
         if (jds .eq. jds) then
            do i = its1, ite1
               lfn(i, jds - 1) = EX(lfn(i, jds), lfn(i, jds + 1))
            end do
            jtso = jds - 1
         end if
         if (jde .eq. jde) then
            do i = its1, ite1
               lfn(i, jde + 1) = EX(lfn(i, jde), lfn(i, jde - 1))
            end do
            jteo = jde + 1
         end if
 
        write (msg, '(8(a,i5))') 'continue_at_boundary: y:', ids, ':', ide, ',', jds, ':', jde, ' ->', its1, ':', ite1, ',', jtso &
                                 , ':', jteo

         call message(msg, level=3)
      end if

      if (ix .ne. 0 .and. iy .ne. 0) then
         if (ids .eq. ids .and. jds .eq. jds) lfn(ids - 1, jds - 1) = EX(lfn(ids, jds), lfn(ids + 1, jds + 1))
         if (ids .eq. ids .and. jde .eq. jde) lfn(ids - 1, jde + 1) = EX(lfn(ids, jde), lfn(ids + 1, jde - 1))
         if (ide .eq. ide .and. jds .eq. jds) lfn(ide + 1, jds - 1) = EX(lfn(ide, jds), lfn(ide - 1, jds + 1))
         if (ide .eq. ide .and. jde .eq. jde) lfn(ide + 1, jde + 1) = EX(lfn(ide, jde), lfn(ide - 1, jde - 1))
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
      ifds, ifde, jfds, jfde, &
      v2, &
      ims, ime, jms, jme, &
      ids, ide, jds, jde, &
      v1)
      implicit none

      integer, intent(in):: ids, ide, jds, jde, ims, ime, jms, jme
      !real, intent(out):: v1(ims:ime,jms:jme)
      integer, intent(in):: ifds, ifde, jfds, jfde, ifms, ifme, jfms, jfme
      !real, intent(in):: v2(ifms:ifme,jfms:jfme)
      real, dimension(ims:ime, jms:jme), intent(out):: v1 ! INTRODUZIDO POR ISILDA CM
      real, dimension(ifms:ifme, jfms:jfme), intent(in):: v2 ! INTRODUZIDO POR ISILDA CM

      integer:: i, i_f, j, j_f, ir, jr, isz1, isz2, jsz1, jsz2, ioff, joff, ibase, jbase
      real :: t
!character(len=128)msg
      character(len=128):: msg !INTRODUZIDO POR ISILDA

      isz1 = ide - ids + 1
      jsz1 = jde - jds + 1
      isz2 = ifde - ifds + 1
      jsz2 = jfde - jfds + 1

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

      do j = jds, jde
         jbase = jfds + jr*(j - jds)
         do i = ids, ide
            ibase = ifds + ir*(i - ids)
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
      write (msg, 91) ifds, ifde, jfds, jfde, ifms, ifme, jfms, jfme
      call message(msg, level=0)
      write (msg, 91) ids, ide, jds, jde, ims, ime, jms, jme
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

   real pure function sum_2darray(ids, ide, jds, jde, &
                                  ims, ime, jms, jme, &
                                  a) result(saida)
      integer, intent(in)::ids, ide, jds, jde, ims, ime, jms, jme
      real, intent(in)::a(ims:ime, jms:jme)

      integer:: i, j
      real:: t

      t = 0.
      do j = jds, jde
         do i = ids, ide
            t = t + a(i, j)
         end do
      end do
!sum_2darray = t
      saida = t
   end function sum_2darray

   real pure function max_2darray(ids, ide, jds, jde, &
                                  ims, ime, jms, jme, &
                                  a) result(saida)
      integer, intent(in)::ids, ide, jds, jde, ims, ime, jms, jme
      real, intent(in)::a(ims:ime, jms:jme)

      integer:: i, j
      real:: T

      t = 0.
      do j = jds, jde
         do i = ids, ide
            t = max(t, a(i, j))
         end do
      end do
!max_2darray = t
      saida = t
   end function max_2darray

   subroutine print_2d_stats_vec(ids, ide, jds, jde, &
                                 ims, ime, jms, jme, &
                                 ax, ay, name)
      implicit none
      integer, intent(in)::ids, ide, jds, jde, ims, ime, jms, jme
      real, intent(in), dimension(ims:ime, jms:jme)::ax, ay
      character(len=*), intent(in)::name
      integer:: i, j
      real:: t
      real:: avg_a, max_a, min_a
      character(len=25)::id
      id = name
      call print_2d_stats(ids, ide, jds, jde, &
                          ims, ime, jms, jme, &
                          ax, id//'/x ')
      call print_2d_stats(ids, ide, jds, jde, &
                          ims, ime, jms, jme, &
                          ay, id//'/y ')
      avg_a = 0
      max_a = -huge(max_a)
      min_a = huge(min_a)
      do j = jds, jde
         do i = ids, ide
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
      avg_a = avg_a/((ide - ids + 1)*(jde - jds + 1))
!print*,'YYYavg_a',avg_a
!call flush(6)
      call print_stat_line(id//'/sz', ids, ide, jds, jde, min_a, max_a, avg_a)
   end subroutine print_2d_stats_vec

   subroutine print_stat_line(name, ids, ide, jds, jde, min_a, max_a, avg_a)

      implicit none

      integer, intent(in)::ids, ide, jds, jde
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
      write (msg, '(a,4i5,3g13.5)') id, ids, ide, jds, jde, min_a, max_a, avg_a
!$OMP END CRITICAL(SFIRE_UTIL_CRIT)
      call message(msg, level=2)
   end subroutine print_stat_line

   subroutine print_3d_stats_by_slice(ids, ide, kds, kde, jds, jde, &
                                      ims, ime, kms, kme, jms, jme, &
                                      a, name)
      implicit none
      integer, intent(in)::ids, ide, jds, jde, ims, ime, jms, jme, kms, kme, kds, kde
      real, intent(in)::a(ims:ime, kms:kme, jms:jme)
      character(len=*), intent(in)::name
      integer::k
      character(len=128)::msg
      
      do k = kds, kde

         write (msg, '(i2,1x,a)') k, name

         call print_3d_stats(ids, ide, k, k, jds, jde, &
                             ims, ime, kms, kme, jms, jme, &
                             a, msg)
      end do
!print*,"estou a sair da rotina print_3d_stats_by_slice"
   end subroutine print_3d_stats_by_slice

subroutine print_3d_stats(ids, ide, kds, kde, jds, jde, &
                             ims, ime, kms, kme, jms, jme, &
                             a, name)
      implicit none
      integer, intent(in)::ids, ide, jds, jde, ims, ime, jms, jme, kms, kme, kds, kde
      real, intent(in)::a(ims:ime, kms:kme, jms:jme)
      character(len=*), intent(in)::name
      integer:: i, j, k
      real:: avg_a, max_a, min_a, t, aa, bb
      character(len=128)::msg

!print*,"estou na rotina print_3d_stats_by_slice dentro da print_3d_stats"
      !print *, 'LFR-DBG para o nome enviado "', trim(name), '"'
      !print *, 'LFR-DBG min e max:', minval(a), maxval(a); call flush (6)
      bb = 0.
      do j = jds, jde
         do k = kds, kde
            do i = ids, ide
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
      do j = jds, jde
         do k = kds, kde
            do i = ids, ide
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
      avg_a = avg_a/((ide - ids + 1)*(jde - jds + 1)*(kde - kds + 1))
!print*,'GGGGGGG',avg_a
!call flush(6)
!print*,"vou levar o valor para print_stat_line"
!call flush(6)
      call print_stat_line(name, ids, ide, jds, jde, min_a, max_a, avg_a)
!print*,"estou fora da print_stat_line"
!call flush(6)
      return
9     continue

      write (msg, 1) name, i, k, j, aa
      call message(msg, level=0)
1     format(a30, '(', i6, ',', i6, ',', i6, ') = ', g13.5)
      write (msg, 2) 'patch dimensions ', ids, ide, kds, kde, jds, jde
      call message(msg, level=0)
      write (msg, 2) 'memory dimensions', ims, ime, kms, kme, jms, jme
      call message(msg, level=0)
2     format(a, 6i8)

      call print_stat_line(name, ids, ide, jds, jde, aa, aa, aa)
      if (aa .ne. aa) goto 10
      msg = 'Invalid floating point number detected in '//name
      call crash(msg)
10    msg = 'NaN detected in '//name
      call crash(msg)
!print*,"estou na rotina print_3d_stats_by_slice vou sair da rotina print_3d_stats"

   end subroutine print_3d_stats


   subroutine print_2d_stats(ids, ide, jds, jde, &
                             ims, ime, jms, jme, &
                             a, name)
      implicit none
      integer, intent(in)::ids, ide, jds, jde, ims, ime, jms, jme
      real, intent(in)::a(ims:ime, jms:jme)
      character(len=*), intent(in)::name

      call print_3d_stats(ids, ide, 1, 1, jds, jde, &
                          ims, ime, 1, 1, jms, jme, &
                          a, name)

   end subroutine print_2d_stats

   real pure function avg_2darray(ids, ide, jds, jde, &
                                  ims, ime, jms, jme, &
                                  a) result(saida)
      integer, intent(in)::ids, ide, jds, jde, ims, ime, jms, jme
      real, intent(in)::a(ims:ime, jms:jme)

!avg_2darray = sum_2darray( ids,ide,jds,jde,               &
!                           ims,ime,jms,jme,               &
!                           a)/((ide-ids+1)*(jde-jds+1))

      saida = sum_2darray(ids, ide, jds, jde, &
                          ims, ime, jms, jme, &
                          a)/((ide - ids + 1)*(jde - jds + 1))

   end function avg_2darray

   real pure function avg_2darray_vec(ids, ide, jds, jde, &
                                      ims, ime, jms, jme, &
                                      ax, ay) result(saida)
      integer, intent(in)::ids, ide, jds, jde, ims, ime, jms, jme
      real, intent(in), dimension(ims:ime, jms:jme):: ax, ay

      integer:: i, j
      real:: t
      t = 0.
      do j = jds, jde
         do i = ids, ide
            t = t + sqrt(ax(i, j)**2 + ay(i, j)**2)
         end do
      end do
      t = t/((ide - ids + 1)*(jde - jds + 1))
!avg_2darray_vec = t
      saida = t
   end function avg_2darray_vec

   subroutine print_array(ids, ide, jds, jde, &
                          ims, ime, jms, jme, &
                          a, name, id)

      integer, intent(in)::ids, ide, jds, jde, ims, ime, jms, jme, id
      real, intent(in), dimension(ims:ime, jms:jme):: a
      character(len=*), intent(in)::name

      integer:: i, j
      character(len=128)::msg

!$OMP CRITICAL(SFIRE_UTIL_CRIT)
      write (msg, *) name, ' start ', id, ' dim ', ids, ide, jds, jde
      call message(msg)
      do j = jds, jde
         do i = ids, ide
            write (msg, *) i, j, a(i, j)
            call message(msg)
         end do
      end do
      write (msg, *) name, ' end ', id
      call message(msg)
!$OMP END CRITICAL(SFIRE_UTIL_CRIT)
   end subroutine print_array

   subroutine write_array_m(ids, ide, jds, jde, &
                            ims, ime, jms, jme, &
                            a, name, id)

      integer, intent(in)::ids, ide, jds, jde, ims, ime, jms, jme, id
      real, intent(in), dimension(ims:ime, jms:jme):: a
      character(len=*), intent(in)::name
      !print *, "estou na subroutine write_array_m vou chamar a write_array_m3"
      call write_array_m3(ids, ide, 1, 1, jds, jde, &
                          ims, ime, 1, 1, jms, jme, &
                          a, name, id)
   end subroutine write_array_m
   

   subroutine write_array_m3(ids, ide, kds, kde, jds, jde, &
                             ims, ime, kms, kme, jms, jme, &
                             a, name, id)
!use module_dm

      implicit none

      integer, intent(in)::ids, ide, jds, jde, ims, ime, jms, jme, kds, kde, kms, kme, id
      real, intent(in), dimension(ims:ime, kms:kme, jms:jme):: a
      character(len=*), intent(in)::name

      integer:: i, j, k, iu, ilen, myproc, nprocs
      logical:: op
      character(len=128)::fname, msg
      !print *, "estou na subroutine write_array_m3"
      if (fire_print_file .eq. 0 .or. id .le. 0) return
      call check_mesh_2dim(ids, ide, jds, jde, ims, ime, jms, jme)
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
      write (iu, 1) real(ids)
      write (iu, 1) real(ide)
      write (iu, 1) real(jds)
      write (iu, 1) real(jde)
      write (iu, 1) real(kds)
      write (iu, 1) real(kde)
      write (iu, 1) (((a(i, k, j), i=ids, ide), j=jds, jde), k=kds, kde)
      close (iu)
      write (msg, 2) name, '(', ids, ':', ide, ',', jds, ':', jde, ',', &
         kds, ':', kde, ') -> ', trim(fname)
!print*,"vou sair da rotina write_array_m3"
!$OMP END CRITICAL(SFIRE_UTIL_CRIT)
      call message(msg)
      return

1     format(e20.12)
2     format(2a, 3(i5, a, i5, a), 2a)
3     format(a, a, i8.8, a)
4     format(a, a, i8.8, a, i4.4, a)

   end subroutine write_array_m3

   subroutine read_array_2d_real(filename, a, ids, ide, jds, jde, ims, ime, jms, jme)
!use module_dm
      implicit none

      integer, intent(in)::ids, ide, jds, jde, ims, ime, jms, jme
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

      mi = ide - ids + 1
      mj = jde - jds + 1
      write (msg, 2) 'reading array size ', mi, mj, ' from file ', trim(filename)
2     format(a, 2i6, 2a)
      call message(msg, level=1)

      call check_mesh_2dim(ids, ide, jds, jde, ims, ime, jms, jme)

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
      do i = ids, ide
         read (iu, *, err=10) (a(i, j), j=jds, jde)
      end do
      close (iu, err=11)
      call print_2d_stats(ids, ide, jds, jde, &
                          ims, ime, jms, jme, &
                          a, filename)
      write (6, *) ids, jds, a(ids, jds), loc(a(ids, jds))
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
   
   

    subroutine print_chsum(id,&
                 ims, ime, kms, kme, jms, jme, &
                ids, ide, kds, kde, jds, jde, a, name)

    ! Defini��o dos argumentos de entrada
    integer, intent(in) :: id
    integer, intent(in) ::ims, ime, kms, kme, jms, jme, & 
                         ids, ide, kds, kde, jds, jde  ! �ndices principais
    real, intent(in), dimension(ims:ime, kms:kme, jms:jme) :: a
    character(len=*) :: name

    ! Vari�veis locais
    integer :: i, j, k, iel, lsum
    real :: rel
    character(len=256) :: msg
    equivalence(rel, iel)
    
    if (fire_print_msg .le. 0) return
    
    ! Inicializa��o
    lsum = 0

    ! Loop sobre a matriz principal para calcular o checksum
    do j = jds, jde
        do k = kds, kde
            do i = ids, ide
                rel = a(i, k, j)
                lsum = ieor(lsum, iel)  ! Opera��o bit a bit XOR para c�lculo do checksum
            end do
        end do
    end do

    ! Impress�o dos resultados do checksum
    write(msg, '(i6, 1x, a10, " dims ", 6i5, " chsum ", z8.8)') id, name, ids, ide, kds, kde, jds, jde, lsum
    call message(msg)

    end subroutine print_chsum


    real function fun_real(fun,  ims, ime, kms, kme, jms, jme, &
                                 ids, ide, kds, kde, jds, jde, a, b)
    implicit none

    ! Argumentos de entrada
    integer, intent(in) :: fun,  ims, ime, kms, kme, jms, jme, &
                                 ids, ide, kds, kde, jds, jde
    real, intent(in), dimension(ims:ime, kms:kme, jms:jme) :: a, b

    ! Vari�veis locais
    real :: lsum, gsum, void, psum
    integer :: i, j, k
    logical :: dosum, domax, domin
    character(len=256) :: msg

    ! Inicializa��o de vari�veis e l�gica para diferentes fun��es
    if (fun .eq. REAL_SUM) then
        void = 0.
        lsum = void
        do j = jds, jde
            do k = kds, kde
                do i = ids, ide
                    lsum = lsum + a(i, k, j)
                end do
            end do
        end do
    elseif (fun .eq. RNRM_SUM) then
        void = 0.
        lsum = void
        do j = jds, jde
            do k = kds, kde
                do i = ids, ide
                    lsum = lsum + sqrt(a(i, k, j)**2 + b(i, k, j)**2)
                end do
            end do
        end do
    elseif (fun .eq. REAL_MAX) then
        void = -huge(lsum)
        lsum = void
        do j = jds, jde
            do k = kds, kde
                do i = ids, ide
                    lsum = max(lsum, a(i, k, j))
                end do
            end do
        end do
    elseif (fun .eq. REAL_AMAX) then
        void = -huge(lsum)
        lsum = void
        do j = jds, jde
            do k = kds, kde
                do i = ids, ide
                    lsum = max(lsum, abs(a(i, k, j)))
                end do
            end do
        end do
    elseif (fun .eq. REAL_MIN) then
        void = huge(lsum)
        lsum = void
        do j = jds, jde
            do k = kds, kde
                do i = ids, ide
                    lsum = min(lsum, a(i, k, j))
                end do
            end do
        end do
    elseif (fun .eq. RNRM_MAX) then
        void = 0.
        lsum = void
        do j = jds, jde
            do k = kds, kde
                do i = ids, ide
                    lsum = max(lsum, sqrt(a(i, k, j)**2 + b(i, k, j)**2))
                end do
            end do
        end do
    else
        call crash('fun_real: bad fun')
    end if

    ! Verifica��o de NaN
    if (lsum .ne. lsum) call message('fun_real: WARNING: NaN detected')

    ! Ajustes para soma, m�ximo e m�nimo
    dosum = fun .eq. REAL_SUM .or. fun .eq. RNRM_SUM
    domax = fun .eq. REAL_MAX .or. fun .eq. REAL_AMAX .or. fun .eq. RNRM_MAX
    domin = fun .eq. REAL_MIN

    ! Calcula o resultado final
    psum = void
    if (dosum) psum = psum + lsum
    if (domax) psum = max(psum, lsum)
    if (domin) psum = min(psum, lsum)

    gsum = psum  ! Simplifica��o para remover paralelismo

    ! Verifica��o final de NaN
    if (gsum .ne. gsum) call message('fun_real: WARNING: NaN detected')

    ! Retorno do valor final
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


   subroutine build_NFFL(ifms, ifme, jfms ,jfme,&
                              ifds, ifde, jfds ,jfde,&
                              fxlat, fxlong,tcomb)

  implicit none
  real :: con_grau
  integer, intent(IN) :: ifms, ifme, jfms, jfme
  integer, intent(IN) :: ifds, ifde, jfds, jfde
  real, dimension(ifms:ifme, jfms:jfme), intent(OUT) :: tcomb
  real, dimension(ifms:ifme, jfms:jfme),intent(IN) :: fxlat, fxlong
  real :: point_latitude, point_longitude, grid_latitude,&
            grid_longitude
  real ::  distance
  integer :: i, j, ia, ja, bb,px,py
  integer, parameter :: escrever = 1 ! FLAG para escrever no ficheiro para ser lido no grads (1 => sim)
  integer, parameter :: rows = 24286
  integer, parameter :: ncols = 16439
  real, parameter :: cellsize = 0.00022511079999887
  real, parameter :: first_lat = 36.855458166893
  real, parameter :: first_lon = -9.7545209749261
  real, parameter :: res_ant_init = 1000.0
 ! real, dimension(rows) :: vlat
 ! real, dimension(ncols) :: vlon
 ! integer, dimension(rows) :: jj
 ! integer, dimension(ncols) :: ii
  real, allocatable :: vlat(:)
  real, allocatable :: vlon(:)
  integer, allocatable :: jj(:)
  integer, allocatable :: ii(:)
 ! real, dimension(ncols,rows) :: idclass
  real, allocatable :: idclass(:,:)
  character(len=255) :: filename
  character(len=255) :: lixo
  logical :: file_exists
  logical :: value_unit
  character(len=80) :: fname
  integer :: irec,recl_size,pp
  real :: res_ant 
  real :: p1,p2,p3,p4,p5
  real :: lat1,lon1,lat2,lon2,delta_lat,delta_lon,a,c,R

  allocate(vlat(rows))
  allocate(vlon(ncols))
  allocate(jj(rows))
  allocate(ii(ncols))
  allocate(idclass(ncols,rows))

  
  filename = "./data_comb/comb_port1.txt"

  inquire (file=trim(filename), exist=file_exists)
  if (.not. file_exists) then
    print *, "nao existe o ficheiro de modelos de combustivel"
    stop
  end if
  

  open (32, file=trim(adjustl(filename)), status='old', &
        access='sequential', form='formatted', action='read')


  read (32, *) lixo, lixo
  read (32, *) lixo, lixo
  read (32, *) lixo, lixo
  read (32, *) lixo, lixo
  read (32, *) lixo, lixo
  read (32, *) lixo, lixo
  read (32, *) ((idclass(i, j), i=1, ncols),j=rows,1,-1 )
   

   py=0
   do j = rows,1,-1
     vlat(j) = first_lat + j*cellsize
     if((vlat(j) >= fxlat(ifds, jfds)) .and.&
          (vlat(j) <= fxlat(ifde,jfde)))then
    py=1+py
    jj(py)=j
    endif
   enddo

   px=0
    do i = 1, ncols
      vlon(i) = first_lon + i*cellsize
   if((vlon(i) >= fxlong(ifds, jfds)) .and.&
          (vlon(i) <= fxlong(ifde,jfde)))then
     px=px+1
     ii(px)=i
    endif
    enddo

  do j = 1, rows
    do i = 1, ncols
      if (idclass(i, j) == -9999) then
        idclass(i, j) = 14.
      end if
    end do
  end do


!******Vai achar a posicao do no na malha cartesiana*********

  con_grau = 3.14159265358979323846/180.0
  R=6372.795477598*1000.0
  tcomb(:,:) = 14.0
  do j = jfds-1,jfde+1
       do i = ifds-1, ifde+1
         res_ant=res_ant_init
         grid_latitude = fxlat(i, j)
         grid_longitude = fxlong(i, j)
          do ja = 1,py
            do ia=1,px
            point_latitude = vlat(jj(ja))
            point_longitude = vlon(ii(ia))
            if((abs(grid_latitude - point_latitude) <= 0.02) .and. &
                (abs(grid_longitude - point_longitude) <= 0.02)) then
                  lat1 = point_latitude*con_grau
                  lon1 = point_longitude*con_grau
                  lat2 = grid_latitude*con_grau
                  lon2 = grid_longitude*con_grau
                  delta_lat = (lat2 - lat1)
                  delta_lon = (lon2 - lon1)
                  a = sin(delta_lat/2.0)**2.0 + cos(lat1)*cos(lat2)*sin(delta_lon/2.0)**2.0
                  c = 2.0*atan2(sqrt(a), sqrt(1.0 - a))
                  distance = R*c
                  if(distance < res_ant) then
                     res_ant = distance
                     tcomb(i, j) = idclass(ii(ia),jj(ja))
                   end if
              endif
            enddo
           end do
         end do
       end do

  
  if (escrever == 1) then
    write (fname, '(a)') 'NFFL_High_resol.bin'
    recl_size=4*((ifde - ifds + 1)*(jfde - jfds + 1))
    open(90,file=trim(adjustl(fname)),form='unformatted',access='direct',&
         status='replace',recl=recl_size)
    irec = 1
     write(90, rec=irec) ((tcomb(i, j), i=ifds,ifde), j=jfds, jfde)
    close (90)
  end if


  deallocate(vlat)
  deallocate(vlon)
  deallocate(jj)
  deallocate(ii)
  deallocate(idclass)


  close(32)


   end subroutine build_NFFL

   subroutine build_ZSF(ifms, ifme, jfms ,jfme,&
                              ifds, ifde, jfds ,jfde,&
                              fxlat, fxlong,ZSF)

  implicit none
  real :: con_grau
  integer, intent(IN) :: ifms, ifme, jfms, jfme
  integer, intent(IN) :: ifds, ifde, jfds, jfde
  real, dimension(ifms:ifme, jfms:jfme), intent(OUT) :: ZSF
  real, dimension(ifms:ifme, jfms:jfme), intent(IN) :: fxlat, fxlong
  real :: point_latitude, point_longitude, grid_latitude,&
            grid_longitude
  real ::  distance
  integer :: i, j, ia, ja, bb,px,py
  integer, parameter :: escrever = 1 ! FLAG para escrever no ficheiro para ser lido no grads (1 => sim)
  integer, parameter :: rows = 18001
  integer, parameter :: ncols = 10801
  real, parameter :: cellsize = 0.00027777777777778
  real, parameter :: first_lat = 36.999861111111
  real, parameter :: first_lon = -10.000138888889
  real, parameter :: res_ant_init = 1000.0
 ! real, dimension(rows) :: vlat
 ! real, dimension(ncols) :: vlon
 ! integer, dimension(rows) :: jj
 ! integer, dimension(ncols) :: ii
  real, allocatable :: vlat(:)
  real, allocatable :: vlon(:)
  integer, allocatable :: jj(:)
  integer, allocatable :: ii(:)
 ! real, dimension(ncols,rows) :: topo_txt
  real, allocatable ::  topo_txt(:,:)
  character(len=255) :: filename
  character(len=255) :: lixo
  logical :: file_exists
  logical :: value_unit
  character(len=80) :: fname
  integer :: irec,recl_size,pp
  real :: res_ant 
  real :: p1,p2,p3,p4,p5
  real :: lat1,lon1,lat2,lon2,delta_lat,delta_lon,a,c,R
  real :: sta,fin
  
  allocate(vlat(rows))
  allocate(vlon(ncols))
  allocate(jj(rows))
  allocate(ii(ncols))
  allocate(topo_txt(ncols,rows))
  
  filename = "./data_TOPO/topo_25m_portugal.txt"

  inquire (file=trim(filename), exist=file_exists)
  if (.not. file_exists) then
    print *, "nao existe o ficheiro de modelos de Topografia"
    stop
  end if
  
  open (32, file=trim(adjustl(filename)), status='old', &
        access='sequential', form='formatted', action='read')

  read (32, *) lixo, lixo
  read (32, *) lixo, lixo
  read (32, *) lixo, lixo
  read (32, *) lixo, lixo
  read (32, *) lixo, lixo
  read (32, *) lixo, lixo
  read (32, *) ((topo_txt(i, j), i=1, ncols),j=rows,1,-1 )
   
  print *, 'LFR->Topo lido'
   py=0
   do j = rows,1,-1
     vlat(j) = first_lat + j*cellsize
     if((vlat(j) >= fxlat(ifds, jfds)) .and.&
          (vlat(j) <= fxlat(ifde,jfde)))then
    py=1+py
    jj(py)=j
    endif
   enddo

   px=0
    do i = 1, ncols
      vlon(i) = first_lon + i*cellsize
   if((vlon(i) >= fxlong(ifds, jfds)) .and.&
          (vlon(i) <= fxlong(ifde,jfde)))then
     px=px+1
     ii(px)=i
    endif
    enddo
  print *, 'LFR->PY-PX'

  do j = 1, rows
    do i = 1, ncols
      if (topo_txt(i, j) == -9999) then
        topo_txt(i, j) = 14.
      end if
    end do
  end do

  print *, 'LFR->topo 14'
!******Vai achar a posicao do no na malha cartesiana*********

  con_grau = 3.14159265358979323846/180.0
  R=6372.795477598*1000.0
  ZSF(:,:) = 0.0
  call cpu_time(sta)
  do j = jfds-1,jfde+1
       do i = ifds-1, ifde+1
         call cpu_time(fin)
         print *,'LFR J=',j," ate ",jfde+1," ","i=",i," ate ",ifde+1,fin-sta
         sta = fin
         res_ant=res_ant_init
         grid_latitude = fxlat(i, j)
         grid_longitude = fxlong(i, j)
          do ja = 1,py
            do ia=1,px
              point_latitude = vlat(jj(ja))
              point_longitude = vlon(ii(ia))
              if((abs(grid_latitude - point_latitude) <= 0.02) .and. &
                 (abs(grid_longitude - point_longitude) <= 0.02)) then
                  lat1 = point_latitude*con_grau
                  lon1 = point_longitude*con_grau
                  lat2 = grid_latitude*con_grau
                  lon2 = grid_longitude*con_grau
                  delta_lat = (lat2 - lat1)
                  delta_lon = (lon2 - lon1)
                  a = sin(delta_lat/2.0)**2.0 + cos(lat1)*cos(lat2)*sin(delta_lon/2.0)**2.0
                  c = 2.0*atan2(sqrt(a), sqrt(1.0 - a))
                  distance = R*c
                  if(distance < res_ant) then
                     res_ant = distance
                     ZSF(i, j) = topo_txt(ii(ia),jj(ja))
                   end if
                endif
               enddo
              end do
           end do
       end do
  

  if (escrever == 1) then
    write (fname, '(a)') 'ZSF_High_resol.bin'
    recl_size=4*((ifde - ifds + 1)*(jfde - jfds + 1))
    open(90,file=trim(adjustl(fname)),form='unformatted',access='direct',&
         status='replace',recl=recl_size)
    irec = 1
     write(90, rec=irec) ((ZSF(i, j), i=ifds,ifde), j=jfds, jfde)
    close (90)
  end if

  deallocate(vlat)
  deallocate(vlon)
  deallocate(jj)
  deallocate(ii)
  deallocate(topo_txt)


  close(32)

  
  end subroutine build_ZSF



end module module_fr_sfire_util

