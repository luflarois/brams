
module module_fr_sfire_core

   use module_fr_sfire_phys, only: fire_params, fire_ros, set_fire_crown_params !set_fire_crown_params introduzido por ISILDA CUNHA MENEZES
   use module_fr_sfire_util

!!!!! ACRESCENTADO POR ISILDA CUNHA MENEZES
!use module_fr_sfire_crown
   use ModNamelistsfireFile, only: grid_config_rec_type
!!!!!!

   implicit none

contains

   subroutine init_no_fire( &
      ifds, ifde, jfds, jfde, &
      ifms, ifme, jfms, jfme, &
      ifts, ifte, jfts, jfte, &
      fdx, fdy, time_start, dt, &
      fuel_frac, fire_area, lfn, tign_in, &
      tign)
      implicit none

      integer, intent(in):: ifds, ifde, jfds, jfde
      integer, intent(in):: ifts, ifte, jfts, jfte
      integer, intent(in):: ifms, ifme, jfms, jfme
      real, intent(in) :: fdx, fdy, time_start, dt
      real, intent(out), dimension(ifms:ifme, jfms:jfme) :: &
         fuel_frac, fire_area, lfn, &
         tign_in
      real, intent(inout), dimension(ifms:ifme, jfms:jfme) :: &
         tign

      intrinsic epsilon

      integer:: i, j
      real:: lfn_init, time_init, time_now
      character(len=128):: msg

      time_now = time_start + dt
      time_init = time_start + 2*dt
      lfn_init = 2*max((ifde - ifds + 1)*fdx, (jfde - jfds + 1)*fdy)

      call print_2d_stats(ifts, ifte, jfts, jfte, ifms, ifme, jfms, jfme, tign, 'init_no_fire start: tign')

      if (fire_perimeter_time > 0. .and. fire_tign_in_time > 0.) then
         call crash('fire_perimeter_time > 0 and fire_tign_in_time > 0')
      end if
      !print *, 'LFR-DBG: cond', fire_perimeter_time, fire_tign_in_time; call flush (6)
      if (fire_perimeter_time > 0.) then

         write (msg, '(a,g14.6)') 'init_no_fire: using given ignition time to replay history to fire_perimeter_time=' &
                                  , fire_perimeter_time
!$OMP CRITICAL(SFIRE_CORE_CRIT)
         call message(msg, level=1)
!$OMP END CRITICAL(SFIRE_CORE_CRIT)
         do j = jfts, jfte
            do i = ifts, ifte
               fuel_frac(i, j) = 1.
               fire_area(i, j) = 0.
               lfn(i, j) = tign(i, j) - time_now
            end do
         end do
      elseif (fire_tign_in_time > 0.) then
         call message('init_no_fire: using ignition from given max fire arrival time', level=1)
         do j = jfts, jfte
            do i = ifts, ifte
               tign_in(i, j) = tign(i, j)
               fuel_frac(i, j) = 1.
               fire_area(i, j) = 0.
               lfn(i, j) = tign_in(i, j) - time_now
            end do
         end do
         call print_2d_stats(ifts, ifte, jfts, jfte, ifms, ifme, jfms, jfme, tign_in, 'init_no_fire: tign_in')
      else
         call message('init_no_fire: setting state to no fire', level=1)
         do j = jfts, jfte
            do i = ifts, ifte
               fuel_frac(i, j) = 1.
               fire_area(i, j) = 0.
               tign(i, j) = time_init
               lfn(i, j) = lfn_init
            end do
         end do
      end if
      !print *, 'LFR-DBG: core', 'tign(3,3):', tign(3, 3)
      call check_lfn_tign('init_no_fire', time_now, ifts, ifte, jfts, jfte, ifms, ifme, jfms, jfme, lfn, tign)
      call print_2d_stats(ifts, ifte, jfts, jfte, ifms, ifme, jfms, jfme, tign, 'init_no_fire: tign')
      call print_2d_stats(ifts, ifte, jfts, jfte, ifms, ifme, jfms, jfme, tign, 'init_no_fire: lfn')

      do j = jfts, jfte
         do i = ifts, ifte
            if (.not. lfn(i, j) > 0.) then
               write (msg, *) 'i,j=', i, j, ' tign=', tign(i, j), ' <= time_now =', time_now
               call message(msg, level=0)
               call crash('init_no_fire: ignition time must be after the end of the first time step')
            end if
         end do
      end do

   end subroutine init_no_fire

   subroutine ignite_from_tign_in( &
      ifds, ifde, jfds, jfde, &
      ifms, ifme, jfms, jfme, &
      ifts, ifte, jfts, jfte, &
      start_ts, end_ts, &
      tign_in, &
      lfn, tign, ignited)
      implicit none

      integer, intent(in):: ifds, ifde, jfds, jfde
      integer, intent(in):: ifts, ifte, jfts, jfte
      integer, intent(in):: ifms, ifme, jfms, jfme
      real, intent(in):: start_ts, end_ts
      real, intent(in), dimension(ifms:ifme, jfms:jfme) :: tign_in
      real, intent(inout), dimension(ifms:ifme, jfms:jfme) :: &
         lfn, tign
      integer, intent(out)::ignited

      integer:: i, j
      real:: lfn_old

      ignited = 0
      if (.not. start_ts < fire_tign_in_time) return

      call check_lfn_tign('ignite_from_tign_in start', end_ts, ifts, ifte, jfts, jfte, ifms, ifme, jfms, jfme, lfn, tign)
      !print *, 'LFR-DBG: tign_in', jfts, jfte, ifts, ifte, size(tign_in, 1), size(tign_in, 2), end_ts; call flush (6)
      do j = jfts, jfte
         do i = ifts, ifte
            !print *, 'LFR-DBG: Loop tign_in:', i, j
            !print *,'LFR-DBG: Loop tign_in:',isnan(tign_in(i,j)); call flush(6)
            if (.not. tign_in(i, j) > end_ts) then
               lfn_old = lfn(i, j)
               tign(i, j) = min(tign(i, j), tign_in(i, j))
               lfn(i, j) = min(lfn(i, j), tign(i, j) - end_ts)
               call check_lfn_tign_ij(i, j, 'ignite_from_tign_in', end_ts, lfn(i, j), tign(i, j))
               if (lfn_old > 0 .and. .not. lfn(i, j) > 0.) ignited = ignited + 1
            end if
         end do
      end do

      call check_lfn_tign('ignite_from_tign_in end', end_ts, ifts, ifte, jfts, jfte, ifms, ifme, jfms, jfme, lfn, tign)
      call print_2d_stats(ifts, ifte, jfts, jfte, ifms, ifme, jfms, jfme, tign, 'init_no_fire: tign')
      call print_2d_stats(ifts, ifte, jfts, jfte, ifms, ifme, jfms, jfme, tign, 'init_no_fire: lfn')

   end subroutine ignite_from_tign_in

   subroutine ignite_fire(ifds, ifde, jfds, jfde, &
                          ifms, ifme, jfms, jfme, &
                          ifts, ifte, jfts, jfte, &
                          ignition_line, &
                          start_ts, end_ts, &
                          coord_xf, coord_yf, &
                          unit_xf, unit_yf, &
                          lfn, tign, ignited)
      implicit none

      integer, intent(in):: ifds, ifde, jfds, jfde
      integer, intent(in):: ifts, ifte, jfts, jfte
      integer, intent(in):: ifms, ifme, jfms, jfme
      type(line_type), intent(in):: ignition_line
      real, intent(in):: start_ts, end_ts
      real, dimension(ifms:ifme, jfms:jfme), intent(in):: &
         coord_xf, coord_yf
      real, intent(in):: unit_xf, unit_yf
      real, intent(inout), dimension(ifms:ifme, jfms:jfme) :: &
         lfn, tign
      integer, intent(out):: ignited

      integer:: i, j
      real::lfn_new, tign_new, time_ign, ax, ay, rels, rele, d

      real:: sx, sy
      real:: ex, ey
      real:: st, et
      character(len=128):: msg
      real::cx2, cy2, dmax, axmin, axmax, aymin, aymax, dmin
      real:: start_x, start_y
      real:: end_x, end_y
      real:: radius
      real:: start_time, end_time
      real:: ros, tos
      integer:: msglevel = 4, smsg = 2
      real:: lfn_min

      call check_lfn_tign('ignite_fire start', end_ts, ifts, ifte, jfts, jfte, ifms, ifme, jfms, jfme, lfn, tign)

      start_x = ignition_line%start_x
      start_y = ignition_line%start_y
      end_x = ignition_line%end_x
      end_y = ignition_line%end_y
      start_time = ignition_line%start_time
      end_time = ignition_line%end_time
      radius = ignition_line%radius
      ros = ignition_line%ros

      tos = radius/ros
      st = start_time
      et = min(end_ts, end_time)

      ignited = 0
      print *, 'LFR-DBG: start_time,end_time:',start_time,end_time,start_time < end_time
      print *, 'LFR-DBG: start_ts,et,tos,end_ts,st:',start_ts,et,tos,end_ts,st
      if (start_ts > et + tos .or. end_ts < st) return

      if (start_time < end_time) then
         print *, 'LFR-DBG: I am inside!'
         rels = 0.
         sx = start_x
         sy = start_y
         rele = (et - start_time)/(end_time - start_time)
         ex = start_x + rele*(end_x - start_x)
         ey = start_y + rele*(end_y - start_y)
      else

         rels = 0.
         rele = 1.
         sx = start_x
         sy = start_y
         ex = end_x
         ey = end_y
      end if

      cx2 = unit_xf*unit_xf
      cy2 = unit_yf*unit_yf

      axmin = coord_xf(ifts, jfts)
      aymin = coord_yf(ifts, jfts)
      axmax = coord_xf(ifte, jfte)
      aymax = coord_yf(ifte, jfte)
      if (fire_print_msg .ge. smsg) then
!$OMP CRITICAL(SFIRE_CORE_CRIT)
         write (msg, '(a,2f11.6,a,2f11.6)') 'IGN from ', sx, sy, ' to ', ex, ey
         call message(msg, level=smsg)
         write (msg, '(a,2f10.2,a,2f10.2,a)') 'IGN timestep [', start_ts, end_ts, '] in [', start_time, end_time, ']'
         call message(msg, level=smsg)
         write (msg, '(a,2g13.6,a,2g13.6)') 'IGN tile coord from  ', axmin, aymin, ' to ', axmax, aymax
         call message(msg, level=smsg)
!$OMP END CRITICAL(SFIRE_CORE_CRIT)
      end if
      ignited = 0
      dmax = 0
      dmin = huge(dmax)
11    format('IGN ', 6(a, g17.7, 1x))
12    format('IGN ', 4(a, 2g13.7, 1x))

      do j = jfts, jfte
         do i = ifts, ifte
            call check_lfn_tign_ij(i, j, 'ignite_fire start', end_ts, lfn(i, j), tign(i, j))
            ax = coord_xf(i, j)
            ay = coord_yf(i, j)

            call nearest(d, time_ign, ax, ay, sx, sy, st, ex, ey, et, cx2, cy2)

            dmax = max(d, dmax)
            dmin = min(d, dmin)

            if (radius < ros*(end_ts - time_ign)) then
               lfn_new = d - radius
               tign_new = time_ign + d/ros
            else
               lfn_new = d - ros*(end_ts - time_ign)
               tign_new = lfn_new/ros + end_ts
            end if

            if (fire_print_msg .ge. msglevel) then
!$OMP CRITICAL(SFIRE_CORE_CRIT)
               write (msg, *) 'IGN1 i,j=', i, j, ' lfn(i,j)=', lfn(i, j), ' tign(i,j)=', tign(i, j)
               call message(msg, level=0)
               write (msg, *) 'IGN2 i,j=', i, j, ' lfn_new= ', lfn_new, ' time_ign= ', time_ign, ' d=', d
               call message(msg, level=msglevel)
!$OMP END CRITICAL(SFIRE_CORE_CRIT)
            end if
            if (.not. lfn_new > 0.) then
               ignited = ignited + 1
            end if
            if (.not. lfn(i, j) < 0. .and. lfn_new < 0.) then
               tign(i, j) = tign_new
               if (fire_print_msg .ge. msglevel) then
!$OMP CRITICAL(SFIRE_CORE_CRIT)
                  write (msg, '(a,2i6,a,2g13.6,a,f10.2,a,2f10.2,a)') 'IGN ignited cell ', i, j, ' at', ax, ay, &
                     ' time', tign(i, j), ' in [', start_ts, end_ts, ']'
                  call message(msg, level=0)
             write (msg, '(a,g10.3,a,f10.2,a,2f10.2,a)') 'IGN distance', d, ' from ignition line at', time_ign, ' in [', st, et, ']'
                  call message(msg, level=msglevel)
!$OMP END CRITICAL(SFIRE_CORE_CRIT)
               end if
               if (tign(i, j) < start_ts .or. tign(i, j) > end_ts) then
!$OMP CRITICAL(SFIRE_CORE_CRIT)
                  write (msg, '(a,2i6,a,f11.6,a,2f11.6,a)') 'WARNING ', i, j, &
                     ' fixing ignition time ', tign(i, j), ' outside of the time step [', start_ts, end_ts, ']'
                  call message(msg)
!$OMP END CRITICAL(SFIRE_CORE_CRIT)
                  tign(i, j) = min(max(tign(i, j), start_ts), end_ts)
               end if
            end if
            lfn(i, j) = min(lfn(i, j), lfn_new)
            if (fire_print_msg .ge. msglevel) then
!$OMP CRITICAL(SFIRE_CORE_CRIT)
               write (msg, *) 'IGN3 i,j=', i, j, ' lfn(i,j)=', lfn(i, j), ' tign(i,j)=', tign(i, j)
               call message(msg, level=msglevel)
!$OMP END CRITICAL(SFIRE_CORE_CRIT)
            end if
            call check_lfn_tign_ij(i, j, 'ignite_fire end', end_ts, lfn(i, j), tign(i, j))
         end do
      end do

      call check_lfn_tign("ignite_fire end:", end_ts, ifts, ifte, jfts, jfte, ifms, ifme, jfms, jfme, lfn, tign)
      if (fire_print_msg .ge. smsg) then
         lfn_min = huge(lfn_min)
         do j = jfts, jfte
            do i = ifts, ifte
               lfn_min = min(lfn_min, lfn(i, j))
            end do
         end do
!$OMP CRITICAL(SFIRE_CORE_CRIT)
         write (msg, '(a,2g13.2,a,g10.2,a,g10.2)') 'IGN units ', unit_xf, unit_yf, ' m max dist ', dmax, ' min', dmin
         call message(msg, level=smsg)
         write (msg, '(a,f6.1,a,f8.1,a,i10,a,g10.2)') 'IGN radius ', radius, ' time of spread', tos, &
            ' ignited nodes', ignited, ' lfn min', lfn_min
         call message(msg, level=smsg)
!$OMP END CRITICAL(SFIRE_CORE_CRIT)
      end if
   end subroutine ignite_fire

   subroutine check_lfn_tign(s, time_now, ifts, ifte, jfts, jfte, ifms, ifme, jfms, jfme, lfn, tign)

      implicit none

      character(len=*), intent(in)::s
      real, intent(in)::time_now
      integer, intent(in):: ifts, ifte, jfts, jfte
      integer, intent(in):: ifms, ifme, jfms, jfme
      real, intent(in), dimension(ifms:ifme, jfms:jfme) :: &
         lfn, tign

      integer:: i, j
      character(len=128)::msg

      do j = jfts, jfte
         do i = ifts, ifte
            call check_lfn_tign_ij(i, j, s, time_now, lfn(i, j), tign(i, j))
         end do
      end do

! !$OMP CRITICAL(SFIRE_CORE_CRIT)

! !$OMP END CRITICAL(SFIRE_CORE_CRIT)

   end subroutine check_lfn_tign

   subroutine check_lfn_tign_ij(i, j, s, time_now, lfnij, tignij)

      implicit none

      integer, intent(in)::i, j
      character(len=*), intent(in)::s
      real, intent(in)::time_now
      real, intent(in):: lfnij, tignij

      character(len=128):: msg, msg1

      if (.not. lfnij < 0. .and. tignij < time_now) then
!$OMP CRITICAL(SFIRE_CORE_CRIT)
1        format(a, i5, ',', i5, a, g13.6, a, g15.8, a, g15.8)
         msg1 = trim(s)//':check_lfn_tign_ij: inconsistent state'
         if (lfnij > 0.) then
            write (msg, 1) 'i,j=', i, j, ' lfn=', lfnij, '>0 tign=', tignij, '< time_now=', time_now
            call message(msg, level=0)
            call crash(msg1)
         else
            write (msg, 1) 'i,j=', i, j, ' lfn=', lfnij, '>=0 tign=', tignij, '< time_now=', time_now
            call message(msg, level=0)
            call warning(msg1, level=0)
         end if
!$OMP END CRITICAL(SFIRE_CORE_CRIT)
      end if

   end subroutine check_lfn_tign_ij

!!!!!!DEC$ ATTRIBUTES FORCEINLINE
   subroutine nearest(d, t, ax, ay, sx, sy, st, ex, ey, et, cx2, cy2)
      implicit none

      real, intent(out):: d, t
      real, intent(in):: ax, ay, sx, sy, st, ex, ey, et, cx2, cy2

      real:: mx, my, dam2, dames, am_es, cos2, dmc2, mcrel, mid_t, dif_t, des2, cx, cy
      character(len=128):: msg
      integer::msglevel = 4

11    format('IGN ', 6(a, g17.7, 1x))
12    format('IGN ', 4(a, 2g13.7, 1x))

      mx = (sx + ex)*0.5
      my = (sy + ey)*0.5
      dam2 = (ax - mx)*(ax - mx)*cx2 + (ay - my)*(ay - my)*cy2
      des2 = (ex - sx)*(ex - sx)*cx2 + (ey - sy)*(ey - sy)*cy2
      dames = dam2*des2
      am_es = (ax - mx)*(ex - sx)*cx2 + (ay - my)*(ey - sy)*cy2
      if (dames > 0) then
         cos2 = (am_es*am_es)/dames
      else
         cos2 = 0.
      end if
      dmc2 = dam2*cos2
      if (4.*dmc2 < des2) then

         mcrel = sign(sqrt(4.*dmc2/des2), am_es)
      elseif (am_es > 0) then
         mcrel = 1.0
      else
         mcrel = -1.0
      end if
      cx = (ex + sx)*0.5 + mcrel*(ex - sx)*0.5
      cy = (ey + sy)*0.5 + mcrel*(ey - sy)*0.5
      d = sqrt((ax - cx)*(ax - cx)*cx2 + (ay - cy)*(ay - cy)*cy2)
      t = (et + st)*0.5 + mcrel*(et - st)*0.5
      if (fire_print_msg .ge. msglevel) then
!$OMP CRITICAL(SFIRE_CORE_CRIT)
         write (msg, 12) 'find nearest to [', ax, ay, '] from [', sx, sy, '] [', ex, ey, ']'
         call message(msg, level=msglevel)
         write (msg, 12) 'end times', st, et, ' scale squared', cx2, cy2
         call message(msg, level=msglevel)
         write (msg, 11) 'nearest at mcrel=', mcrel, 'from the midpoint, t=', t
         call message(msg, level=msglevel)
         write (msg, 12) 'nearest is [', cx, cy, '] d=', d
         call message(msg, level=msglevel)
         write (msg, 11) 'dam2=', dam2, 'des2=', des2, 'dames=', dames
         call message(msg, level=msglevel)
         write (msg, 11) 'am_es=', am_es, 'cos2=', cos2, 'dmc2=', dmc2
         call message(msg)
!$OMP END CRITICAL(SFIRE_CORE_CRIT)
      end if
   end subroutine nearest

   subroutine fuel_left( &
      ifds, ifde, jfds, jfde, &
      ims, ime, jms, jme, &
      its, ite, jts, jte, &
      ifs, ife, jfs, jfe, &
      lfn, tign, fuel_time, time_now, fuel_frac, fire_area)
      implicit none

      integer, intent(in) ::ifds, ifde, jfds, jfde, its, ite, jts, jte, ims, ime &
                             , jms, jme, ifs, ife, jfs, jfe
      real, intent(in), dimension(ims:ime, jms:jme)::lfn, tign, fuel_time
      real, intent(in):: time_now
      real, intent(out), dimension(ifs:ife, jfs:jfe)::fuel_frac
      real, intent(out), dimension(ims:ime, jms:jme):: fire_area

      integer::i, j, ir, jr, icl, jcl, isubcl, jsubcl, i2, j2, ii, jj, its1, jts1, ite1, jte1
      real::fmax, frat, helpsum1, helpsum2, fuel_left_ff, fire_area_ff, rx, ry, tignf(2, 2)
      real, dimension(3, 3)::tff, lff

      real::lffij, lffi1j, lffij1, lffi1j1, tifij, tifi1j, tifij1, tifi1j1, tx, ty, txx, tyy

      character(len=128)::msg
      integer::m, omp_get_thread_num

      call check_mesh_2dim(its - 1, ite + 1, jts - 1, jte + 1, ims, ime, jms, jme)
      call check_mesh_2dim(its, ite, jts, jte, ifs, ife, jfs, jfe)
      call check_lfn_tign('fuel_left start', time_now, its, ite, jts, jte, ims, ime, jms, jme, lfn, tign)

      ir = fuel_left_irl
      jr = fuel_left_jrl

      if ((ir .ne. 2) .or. (jr .ne. 2)) then
         call crash('fuel_left: ir.ne.2 or jr.ne.2 ')
      end if

      rx = 1./ir
      ry = 1./jr

      its1 = max(its, ifds + 1)
      ite1 = min(ite, ifde - 1)
      jts1 = max(jts, jfds + 1)
      jte1 = min(jte, jfde - 1)

      do j = jts, jte
         do i = its, ite
            fuel_frac(i, j) = 1.
            fire_area(i, j) = 0.
         end do
      end do

      call check_lfn_tign('fuel_left', time_now, its, ite, jts, jte, ims, ime, jms, jme, lfn, tign)

      do icl = its1, ite1
         do jcl = jts1, jte1
            helpsum1 = 0
            helpsum2 = 0

            call tign_lfn_interpolation(time_now, icl, jcl, ims, ime, jms, jme, &
                                        tign, lfn, tff, lff)

            do isubcl = 1, ir
               do jsubcl = 1, jr
                  if (fuel_left_method .eq. 1) then
                     call fuel_left_cell_1(fuel_left_ff, fire_area_ff, &
                               lff(isubcl, jsubcl), lff(isubcl, jsubcl + 1), lff(isubcl + 1, jsubcl), lff(isubcl + 1, jsubcl + 1), &
                               tff(isubcl, jsubcl), tff(isubcl, jsubcl + 1), tff(isubcl + 1, jsubcl), tff(isubcl + 1, jsubcl + 1), &
                                           time_now, fuel_time(icl, jcl))
                  elseif (fuel_left_method .eq. 2) then
                     call fuel_left_cell_2(fuel_left_ff, fire_area_ff, &
                               lff(isubcl, jsubcl), lff(isubcl, jsubcl + 1), lff(isubcl + 1, jsubcl), lff(isubcl + 1, jsubcl + 1), &
                               tff(isubcl, jsubcl), tff(isubcl, jsubcl + 1), tff(isubcl + 1, jsubcl), tff(isubcl + 1, jsubcl + 1), &
                                           time_now, fuel_time(icl, jcl))

                  else
                     call crash('fuel_left: unknown fuel_left_method')
                  end if

                  if (fire_area_ff .lt. -1e-6 .or. &
                      (fire_area_ff .eq. 0. .and. fuel_left_ff .lt. 1.-1e-6)) then
!$OMP CRITICAL(SFIRE_CORE_CRIT)
                     write (msg, '(a,2i6,2(a,f11.8))') 'fuel_left: at node', i, j, &
                        ' of refined mesh fuel burnt', 1 - fuel_left_ff, ' fire area', fire_area_ff
!$OMP END CRITICAL(SFIRE_CORE_CRIT)
                     call crash(msg)
                  end if

                  helpsum1 = helpsum1 + fuel_left_ff
                  helpsum2 = helpsum2 + fire_area_ff

               end do
            end do

            fuel_frac(icl, jcl) = helpsum1/(ir*jr)
            fire_area(icl, jcl) = helpsum2/(ir*jr)
         end do
      end do

!!$OMP CRITICAL(SFIRE_CORE_CRIT)

!!$OMP END CRITICAL(SFIRE_CORE_CRIT)

!$OMP CRITICAL(SFIRE_CORE_CRIT)

!$OMP END CRITICAL(SFIRE_CORE_CRIT)

      return

   end subroutine fuel_left

   subroutine tign_lfn_interpolation(time_now, icl, jcl, ims, ime, jms, jme, &
                                     tign, lfn, tff, lff)
      real, intent(in):: time_now
      integer, intent(in) :: icl, jcl
      integer, intent(in) :: ims, ime, jms, jme
      real, intent(in), dimension(ims:ime, jms:jme)::lfn, tign
      real, intent(out), dimension(3, 3)::tff, lff

      call check_lfn_tign('tign_lfn_interpolation', time_now, icl - 1, icl, jcl - 1, jcl, ims, ime, jms, jme, lfn, tign)

      call tign_lfn_four_pnts_interp(tign(icl - 1, jcl - 1), tign(icl - 1, jcl), tign(icl, jcl - 1), &
                                     tign(icl, jcl), lfn(icl - 1, jcl - 1), lfn(icl - 1, jcl), &
                                     lfn(icl, jcl - 1), lfn(icl, jcl), lff(1, 1), tff(1, 1), time_now)

      call tign_lfn_line_interp(tign(icl - 1, jcl), tign(icl, jcl), lfn(icl - 1, jcl), lfn(icl, jcl), &
                                lff(1, 2), tff(1, 2), time_now, icl - 1, jcl, icl, jcl)

      call tign_lfn_four_pnts_interp(tign(icl - 1, jcl), tign(icl - 1, jcl + 1), tign(icl, jcl), &
                                     tign(icl, jcl + 1), lfn(icl - 1, jcl), lfn(icl - 1, jcl + 1), &
                                     lfn(icl, jcl), lfn(icl, jcl + 1), lff(1, 3), tff(1, 3), time_now)

      call tign_lfn_line_interp(tign(icl, jcl - 1), tign(icl, jcl), lfn(icl, jcl - 1), lfn(icl, jcl), &
                                lff(2, 1), tff(2, 1), time_now, icl, jcl - 1, icl, jcl)

      lff(2, 2) = lfn(icl, jcl)
      tff(2, 2) = tign(icl, jcl)

      call tign_lfn_line_interp(tign(icl, jcl + 1), tign(icl, jcl), lfn(icl, jcl + 1), lfn(icl, jcl), &
                                lff(2, 3), tff(2, 3), time_now, icl, jcl + 1, icl, jcl)

      call tign_lfn_four_pnts_interp(tign(icl, jcl - 1), tign(icl, jcl), tign(icl + 1, jcl - 1), &
                                     tign(icl + 1, jcl), lfn(icl, jcl - 1), lfn(icl, jcl), &
                                     lfn(icl + 1, jcl - 1), lfn(icl + 1, jcl), lff(3, 1), tff(3, 1), time_now)

      call tign_lfn_line_interp(tign(icl + 1, jcl), tign(icl, jcl), lfn(icl + 1, jcl), lfn(icl, jcl), &
                                lff(3, 2), tff(3, 2), time_now, icl + 1, jcl, icl, jcl)

      call tign_lfn_four_pnts_interp(tign(icl, jcl), tign(icl, jcl + 1), tign(icl + 1, jcl), &
                                     tign(icl + 1, jcl + 1), lfn(icl, jcl), lfn(icl, jcl + 1), &
                                     lfn(icl + 1, jcl), lfn(icl + 1, jcl + 1), lff(3, 3), tff(3, 3), time_now)

   end subroutine tign_lfn_interpolation

   subroutine tign_lfn_line_interp(tign1, tign2, lfn1, lfn2, lfn_subcl, tign_subcl, time_now, i1, j1, i2, j2)

      real, intent(in) :: tign1, tign2
      real, intent(in) :: lfn1, lfn2
      real, intent(in) :: time_now
      real, intent(out) :: lfn_subcl, tign_subcl
      integer, intent(in):: i1, j1, i2, j2

      real :: c
      character(len=128)::msg

      call check_lfn_tign_ij(i1, j1, 'tign_lfn_line_interp', time_now, lfn1, tign1)
      call check_lfn_tign_ij(i2, j2, 'tign_lfn_line_interp', time_now, lfn2, tign2)

      lfn_subcl = 0.5*(lfn1 + lfn2)

      if (.not. lfn_subcl < 0.) then
         tign_subcl = time_now
      elseif ((lfn1 < 0.) .and. (lfn2 < 0.)) then
         tign_subcl = 0.5*(tign1 + tign2)
      else
         if (lfn1 < 0.) then
            c = (tign1 - time_now)/lfn1
         elseif (lfn2 < 0.) then
            c = (tign2 - time_now)/lfn2
         else
            call crash('tign_lfn_line_interp: one of lfn1 or lfn2 should be < 0')
         end if
         if (c < 0.) then
!$OMP CRITICAL(SFIRE_CORE_CRIT)
            write (msg, *) 'tign1,tign2,time_now ', tign1, tign2, time_now
            call message(msg, level=0)
            write (msg, *) 'tign1-time_now,lfn1 ', tign1 - time_now, lfn1
            call message(msg, level=0)
            write (msg, *) 'tign2-time_now,lfn2 ', tign2 - time_now, lfn2
            call message(msg, level=0)
            write (msg, *) 'c ', c
            call message(msg, level=0)
!$OMP END CRITICAL(SFIRE_CORE_CRIT)
            call crash('tign_lfn_line_interp: bad ignition times, c<0')
         end if
         tign_subcl = c*lfn_subcl + time_now; 
      end if
   end subroutine tign_lfn_line_interp

   subroutine tign_lfn_four_pnts_interp(tign1, tign2, tign3, tign4, &
                                        lfn1, lfn2, lfn3, lfn4, lfn_subcl, tign_subcl, time_now)

      real, intent(in) :: tign1, tign2, tign3, tign4
      real, intent(in) :: lfn1, lfn2, lfn3, lfn4
      real, intent(in) :: time_now
      real, intent(out) :: lfn_subcl, tign_subcl

      real :: a, b, c, scale
      character(len=128)::msg

      call check_lfn_tign_ij(0, 1, 'tign_lfn_four_pnts_interp', time_now, lfn1, tign1)
      call check_lfn_tign_ij(0, 2, 'tign_lfn_four_pnts_interp', time_now, lfn2, tign2)
      call check_lfn_tign_ij(0, 3, 'tign_lfn_four_pnts_interp', time_now, lfn3, tign3)
      call check_lfn_tign_ij(0, 4, 'tign_lfn_four_pnts_interp', time_now, lfn4, tign4)

      lfn_subcl = 0.25*(lfn1 + lfn2 + lfn3 + lfn4)

      if (.not. lfn_subcl < 0.) then

         tign_subcl = time_now
      elseif ((lfn1 < 0.) .and. (lfn2 < 0.) .and. (lfn3 < 0.) .and. (lfn4 < 0.)) then

         tign_subcl = 0.25*(tign1 + tign2 + tign3 + tign4)
      else

         scale = -minval((/lfn1, lfn2, lfn3, lfn4/))
         a = 0.
         b = 0.
         if (lfn1 < 0.) then
            a = a + (lfn1/scale)**2
            b = b + (tign1 - time_now)*(lfn1/scale)
         end if
         if (lfn2 < 0.) then
            a = a + (lfn2/scale)**2
            b = b + (tign2 - time_now)*(lfn2/scale)
         end if
         if (lfn3 < 0.) then
            a = a + (lfn3/scale)**2
            b = b + (tign3 - time_now)*(lfn3/scale)
         end if
         if (lfn4 < 0.) then
            a = a + (lfn4/scale)**2
            b = b + (tign4 - time_now)*(lfn4/scale)
         end if
         if (.not. a > 0. .or. .not. scale > 0.) then
!$OMP CRITICAL(SFIRE_CORE_CRIT)
            write (msg, *) 'lfn1,2,3,4,_subcl,ssq = ', lfn1, lfn2, lfn3, lfn3, lfn_subcl, a, ' scale=', scale
!$OMP END CRITICAL(SFIRE_CORE_CRIT)
            call message(msg, level=0)
            call crash('tign_lfn_four_pnts_interp: internal: not a>0, none of lfn < 0 ?')
         end if
         if (b < 0.) call crash('tign_lfn_four_pnts_interp: bad ignition times, b<0')
         c = b/a
         tign_subcl = c*(lfn_subcl/scale) + time_now; 
      end if

   end subroutine tign_lfn_four_pnts_interp

   subroutine fuel_left_cell_1(fuel_frac_left, fire_frac_area, &
                               lfn00, lfn01, lfn10, lfn11, &
                               tign00, tign01, tign10, tign11, &
                               time_now, fuel_time_cell)

      implicit none

      real, intent(out):: fuel_frac_left, fire_frac_area
      real, intent(in)::lfn00, lfn01, lfn10, lfn11
      real, intent(in)::tign00, tign01, tign10, tign11
      real, intent(in)::time_now
      real, intent(in)::fuel_time_cell

      intrinsic tiny

      real::ps, aps, area, ta, out
      real::t00, t01, t10, t11
      real, parameter::safe = tiny(aps)
      character(len=128)::msg

      t00 = tign00 - time_now
      if (lfn00 > 0. .or. t00 > 0.) t00 = 0.
      t01 = tign01 - time_now
      if (lfn01 > 0. .or. t01 > 0.) t01 = 0.
      t10 = tign10 - time_now
      if (lfn10 > 0. .or. t10 > 0.) t10 = 0.
      t11 = tign11 - time_now
      if (lfn11 > 0. .or. t11 > 0.) t11 = 0.

      ps = lfn00 + lfn01 + lfn10 + lfn11
      aps = abs(lfn00) + abs(lfn01) + abs(lfn10) + abs(lfn11)
      aps = max(aps, safe)
      area = (-ps/aps + 1.)/2.
      area = max(area, 0.)
      area = min(area, 1.)

      ta = 0.25*(t00 + t01 + t10 + t11)

      out = 1.

      if (area > 0) out = area*exp(ta/fuel_time_cell) + (1.-area)

      if (out > 1.) then
!$OMP CRITICAL(SFIRE_CORE_CRIT)
         write (msg, *) 'out=', out, '>1 area=', area, ' ta=', ta
         call message(msg)
         write (msg, *) 'tign=', tign00, tign01, tign10, tign11, ' time_now=', time_now
!$OMP END CRITICAL(SFIRE_CORE_CRIT)
         call message(msg)

         call crash('fuel_left_cell_1: fuel fraction > 1')
      end if

      fuel_frac_left = out
      fire_frac_area = area

   end subroutine fuel_left_cell_1

   subroutine fuel_left_cell_2(fuel_frac_left, fire_frac_area, &
                               lfn00, lfn01, lfn10, lfn11, &
                               tign00, tign01, tign10, tign11, &
                               time_now, fuel_time_cell)

      implicit none

      real, intent(out):: fuel_frac_left, fire_frac_area
      real, intent(in)::lfn00, lfn01, lfn10, lfn11
      real, intent(in)::tign00, tign01, tign10, tign11
      real, intent(in)::time_now
      real, intent(in)::fuel_time_cell

      intrinsic tiny

      real(kind=8)::ps, aps, area, ta, out
      real(kind=8)::t00, t01, t10, t11
      real(kind=8), parameter::safe = tiny(aps)
      real(kind=8)::dx, dy
      integer::i, j, k

      real(kind=8), dimension(3)::u

      real(kind=8)::tweight, tdist
      integer::kk, ll, ss
      real(kind=8)::rnorm
      real(kind=8), dimension(8, 2)::xylist, xytlist
      real(kind=8), dimension(8)::tlist, llist, xt
      real(kind=8), dimension(5)::xx, yy
      real(kind=8), dimension(5)::lfn, tign

      integer:: npoint
      real(kind=8)::tt, x0, y0, xts, xte, yts, yte, xt1, xt2
      real(kind=8)::lfn0, lfn1, dist, nr, s, errQ, ae, ce, ceae, a0, a1, a2, d, cet
      real(kind=8)::s1, s2, s3
      real(kind=8)::upper, lower, ah, ch, aa, cc, aupp, cupp, alow, clow
      real(kind=8), dimension(2, 2)::mQ
      real(kind=8), dimension(2)::ut
      character(len=128)::msg

      intrinsic epsilon

      real(kind=8), parameter:: zero = 0., one = 1., eps = epsilon(zero)

      real(kind=8)::tign_middle, dt_dx, dt_dy, lfn_middle, a, b, c
      real(kind=8):: alg_err

      call crash('fuel_left_method=2 not working, please use fuel_left_method=1')

      call check_lfn_tign_ij(0, 0, 'fuel_left_cell_2', time_now, lfn00, tign00)
      call check_lfn_tign_ij(0, 1, 'fuel_left_cell_2', time_now, lfn01, tign01)
      call check_lfn_tign_ij(1, 0, 'fuel_left_cell_2', time_now, lfn10, tign10)
      call check_lfn_tign_ij(1, 1, 'fuel_left_cell_2', time_now, lfn11, tign11)

      alg_err = 0
      dx = 1
      dy = 1
      t00 = time_now - tign00
      if (lfn00 >= 0. .or. t00 < 0.) t00 = 0.
      t01 = time_now - tign01
      if (lfn01 >= 0. .or. t01 < 0.) t01 = 0.
      t10 = time_now - tign10
      if (lfn10 >= 0. .or. t10 < 0.) t10 = 0.
      t11 = time_now - tign11
      if (lfn11 >= 0. .or. t11 < 0.) t11 = 0.

      ps = lfn00 + lfn01 + lfn10 + lfn11
      aps = abs(lfn00) + abs(lfn01) + abs(lfn10) + abs(lfn11)
      aps = max(aps, safe)
      area = (-ps/aps + 1.)/2.
      area = max(area, zero)
      area = min(area, one)

      if (lfn00 >= 0 .and. lfn10 >= 0 .and. lfn01 >= 0 .and. lfn11 >= 0) then
         out = 1.0
         area = 0.

      else if (lfn00 <= 0 .and. lfn10 <= 0 .and. lfn01 <= 0 .and. lfn11 <= 0) then

         tign_middle = (t00 + t01 + t10 + t11)/4

         dt_dx = (t10 - t00 + t11 - t01)/2
         dt_dy = (t01 - t00 + t11 - t10)/2

         u(1) = dt_dx
         u(2) = dt_dy
         u(3) = tign_middle - (dt_dx + dt_dy)/2

         u(1) = -u(1)/fuel_time_cell
         u(2) = -u(2)/fuel_time_cell
         u(3) = -u(3)/fuel_time_cell
         s1 = u(1)
         s2 = u(2)
         out = exp(u(3))*intexp(s1)*intexp(s2)
         area = 1
         if (out < 0 .or. out > 1.0) then
            call message('WARNING: fuel_left_cell: case all burning: out should be between 0 and 1')
         end if

      else

         xx(1) = -0.5
         xx(2) = 0.5
         xx(3) = 0.5
         xx(4) = -0.5
         xx(5) = -0.5
         yy(1) = -0.5
         yy(2) = -0.5
         yy(3) = 0.5
         yy(4) = 0.5
         yy(5) = -0.5
         lfn(1) = lfn00
         lfn(2) = lfn10
         lfn(3) = lfn11
         lfn(4) = lfn01
         lfn(5) = lfn00
         tign(1) = t00
         tign(2) = t10
         tign(3) = t11
         tign(4) = t01
         tign(5) = t00
         npoint = 0

         do k = 1, 4
            lfn0 = lfn(k)
            lfn1 = lfn(k + 1)
            if (lfn0 <= 0.0) then
               npoint = npoint + 1
               xylist(npoint, 1) = xx(k)
               xylist(npoint, 2) = yy(k)
               tlist(npoint) = -tign(k)
               llist(npoint) = lfn0
            end if
            if (lfn0*lfn1 < 0) then
               npoint = npoint + 1

               tt = lfn0/(lfn0 - lfn1)
               x0 = xx(k) + (xx(k + 1) - xx(k))*tt
               y0 = yy(k) + (yy(k + 1) - yy(k))*tt
               xylist(npoint, 1) = x0
               xylist(npoint, 2) = y0
               tlist(npoint) = 0
               llist(npoint) = 0
            end if
         end do

         tlist(npoint + 1) = tlist(1)
         llist(npoint + 1) = llist(1)
         xylist(npoint + 1, 1) = xylist(1, 1)
         xylist(npoint + 1, 2) = xylist(1, 2)

         lfn_middle = (lfn00 + lfn01 + lfn10 + lfn11)/4
         dt_dx = (lfn10 - lfn00 + lfn11 - lfn01)/2
         dt_dy = (lfn01 - lfn00 + lfn11 - lfn10)/2
         u(1) = dt_dx
         u(2) = dt_dy
         u(3) = lfn_middle - (dt_dx + dt_dy)/2

         a = 0
         b = 0

         if (lfn00 <= 0) then
            a = a + lfn00*lfn00
            if (t00 < 0) then
               call crash('fuel_burnt_fd: tign(i1) should be less then time_now')
            else
               b = b + t00*lfn00
            end if
         end if

         if (lfn01 <= 0) then
            a = a + lfn01*lfn01
            if (t01 < 0) then
               call crash('fuel_burnt_fd: tign(i1) should be less then time_now')
            else
               b = b + t01*lfn01
            end if
         end if

         if (lfn10 <= 0) then
            a = a + lfn10*lfn10
            if (t10 < 0) then
               call crash('fuel_burnt_fd: tign(i1) should be less then time_now')
            else
               b = b + t10*lfn10
            end if
         end if

         if (lfn11 <= 0) then
            a = a + lfn11*lfn11
            if (t11 < 0) then
               call crash('fuel_burnt_fd: tign(i1) should be less then time_now')
            else
               b = b + t11*lfn11
            end if
         end if

         if (a == 0) then
            call crash('fuel_burnt_fd: if c is on fire then one of cells should be on fire')
         end if
         c = b/a
         u(1) = u(1)*c
         u(2) = u(2)*c
         u(3) = u(3)*c

         nr = sqrt(u(1)**2 + u(2)**2)
         if (.not. nr .gt. eps) then
            out = 1.
            goto 900
         end if
         c = u(1)/nr
         s = u(2)/nr

         mQ(1, 1) = c
         mQ(1, 2) = s
         mQ(2, 1) = -s
         mQ(2, 2) = c

         call matvec(mQ, 2, 2, u, 3, ut, 2, 2, 2)
         errQ = ut(2)
         ae = -ut(1)/fuel_time_cell
         ce = -u(3)/fuel_time_cell
         cet = ce
         call matmatp(xylist, 8, 2, mQ, 2, 2, xytlist, 8, 2, npoint + 1, 2, 2)
         call sortxt(xytlist, 8, 2, xt, 8, npoint)
         out = 0.0
         aupp = 0.0
         cupp = 0.0
         alow = 0.0
         clow = 0.0
         do k = 1, npoint - 1
            xt1 = xt(k)
            xt2 = xt(k + 1)
            upper = 0
            lower = 0
            ah = 0
            ch = 0
            if (xt2 - xt1 > eps*100) then

               do ss = 1, npoint
                  xts = xytlist(ss, 1)
                  yts = xytlist(ss, 2)
                  xte = xytlist(ss + 1, 1)
                  yte = xytlist(ss + 1, 2)

                  if ((xts > xt1 .and. xte > xt1) .or. &
                      (xts < xt2 .and. xte < xt2)) then
                     aa = 0
                     cc = 0
                  else
                     aa = (yts - yte)/(xts - xte)
                     cc = (xts*yte - xte*yts)/(xts - xte)
                     if (xte < xts) then
                        aupp = aa
                        cupp = cc
                        ah = ah + aa
                        ch = ch + cc
                        upper = upper + 1
                     else
                        alow = aa
                        clow = cc
                        lower = lower + 1
                     end if
                  end if
               end do
               ce = cet

               if (ae*xt1 + ce > 0) then
                  ce = ce - (ae*xt1 + ce)
               end if
               if (ae*xt2 + ce > 0) then
                  ce = ce - (ae*xt2 + ce)
               end if

               ah = aupp - alow
               ch = cupp - clow

               s1 = (xt2 - xt1)*((1./2.)*ah*(xt2 + xt1) + ch)

               ceae = ce/ae; 
               s2 = -ch*exp(ae*(xt1 + ceae))*(xt2 - xt1)*intexp(ae*(xt2 - xt1))

               a2 = (xt1 - xt2)*((1./4.)*(xt1 + xt2)*ceae**2 + (1./3.)* &
                                 (xt1**2 + xt1*xt2 + xt2**2)*ceae + (1./8.)* &
                                 (xt1**3 + xt1*(xt2**2) + xt1**2*xt2 + xt2**3))
               d = (ae**4)*a2

               if (abs(d) > eps) then

                  s3 = (exp(ae*(xt1 + ceae))*(ae*xt1 - 1) - &
                        exp(ae*(xt2 + ceae))*(ae*xt2 - 1))/(ae**2)

               else

                  a1 = (xt1 - xt2)*((1./2.)*ceae*(xt1 + xt2) + (1./3.)* &
                                    (xt1**2 + xt1*xt2 + xt2**2))

                  a0 = (1./2.)*(xt1 - xt2)*(xt1 + xt2)
                  s3 = a0 + a1*ae + a2*ae**2; 
               end if

               s3 = ah*s3
               out = out + s1 + s2 + s3
               if (out < 0. .or. out > 1.) then
                  write (msg, '(a,g14.4,a)') 'WARNING::fuel_fraction ', out, ' should be between 0 and 1'
               end if
if (out .ne. out .or. .not. out .le. huge(out) .or. .not. out .ge. -huge(out)) call crash('fuel_fraction out is not a valid number')

            end if
         end do

         out = 1 - out
      end if

900   continue
      fuel_frac_left = out
      fire_frac_area = area
      if (isnotfinite(fuel_frac_left)) call crash('fuel_frac_left is not a valid number')
      if (isnotfinite(fire_frac_area)) call crash('fire_frac_area is not a valid number')
   end subroutine fuel_left_cell_2

   real function intexp(ab)
      implicit none
      real(kind=8)::ab

      intrinsic epsilon

      real, parameter:: zero = 0., one = 1., eps = epsilon(zero)

      if (eps < abs(ab)**3/6.) then
         intexp = (exp(ab) - 1)/ab
      else
         intexp = 1 + ab/2.
      end if
   end function intexp

   subroutine sortxt(xytlist, nrow, ncolumn, xt, nxt, nvec)
      implicit none
      integer::nrow, ncolumn, nxt, nvec
      real(kind=8), dimension(nrow, ncolumn)::xytlist
      real(kind=8), dimension(nxt)::xt

      integer::i, j
      real(kind=8)::temp

      do i = 1, nvec
         xt(i) = xytlist(i, 1)
      end do

      do i = 1, nvec - 1
         do j = i + 1, nvec
            if (xt(i) > xt(j)) then
               temp = xt(i)
               xt(i) = xt(j)
               xt(j) = temp
            end if
         end do
      end do

   end subroutine sortxt

   subroutine matvec(A, m, n, V, nv, out, nout, nrow, ncolumn)
      implicit none
      integer::m, n, nv, nout, nrow, ncolumn
      real(kind=8), dimension(m, n)::A
      real(kind=8), dimension(nv)::V
      real(kind=8), dimension(nout)::out

      integer::i, j

      do i = 1, nrow
         out(i) = 0.0
         do j = 1, ncolumn
            out(i) = out(i) + A(i, j)*V(j)
         end do
      end do
   end subroutine

   subroutine matmatp(A, mA, nA, B, mB, nB, C, mC, nC, nrow, ncolumn, nP)
      implicit none
      integer::mA, nA, mB, nB, mC, nC, nrow, ncolumn, nP
      real(kind=8), dimension(mA, nA)::A
      real(kind=8), dimension(mB, nB)::B
      real(kind=8), dimension(mC, nC)::C
      integer::i, j, k
      do i = 1, nrow
         do j = 1, ncolumn
            C(i, j) = 0.0
            do k = 1, nP
               C(i, j) = C(i, j) + A(i, k)*B(j, k)
            end do
         end do
      end do
   end subroutine

   subroutine prop_ls(ifun, id, &
                      ipart, &
                      ids, ide, jds, jde, &
                      ims, ime, jms, jme, &
                      ips, ipe, jps, jpe, &
                      its, ite, jts, jte, &
                      ts, dt, dx, dy, &
                      tbound, &
                      lfn_in, lfn_out, tign, ros, &
                      fp, config_flags &
                      ) !config_flags e ifun foi introd por ISILDA CM
      implicit none

      integer, intent(in)::id, ipart, ims, ime, jms, jme, ids, ide, jds, jde, its, ite, jts, jte, ips, ipe, jps, jpe
      real, dimension(ims:ime, jms:jme), intent(inout)::lfn_in, tign
      real, dimension(ims:ime, jms:jme), intent(out)::lfn_out, ros
      real, intent(in)::dx, dy, ts, dt
      real, intent(out)::tbound
!### COMENTADO POR ISILDA CUNHA MENEZES
!type(fire_params),intent(in)::fp
!### INTRODUZIDO POR ISILDA CUNHA MENEZES
      type(fire_params), intent(inout)::fp
      type(grid_config_rec_type), intent(IN)  :: config_flags
!#####

      real, dimension(its - 1:ite + 1, jts - 1:jte + 1):: tend, lfn1

      real::grad2, rr, tbound2, a, a1

      real::gradx, grady, aspeed, err, aerr, time_future, time_now
      real::tmp, t0, t1, t2
      integer::ihs, ihe, jhs, jhe
      integer::ihs2, ihe2, jhs2, jhe2
      integer::itso, iteo, jtso, jteo
      integer::i, j, its1, ite1, jts1, jte1, k, kk, id1
      character(len=128)::msg
      integer::nfirenodes, nfireline, ierrx
      real::sum_err, min_err, max_err, sum_aerr, min_aerr, max_aerr

      integer, parameter :: mstep = 1000, printl = 1
      real, parameter:: zero = 0., one = 1., eps = epsilon(zero), tol = 100*eps, &
                        safe = 2., rmin = safe*tiny(zero), rmax = huge(zero)/safe

      intrinsic max, min, sqrt, nint, epsilon, tiny, huge
      integer, intent(in)::ifun !introduzido por ISILDA CM
      !print *, "estou no module_fr_sfire_core dentro da prop_ls"
      !call flush (6)
!$OMP CRITICAL(SFIRE_CORE_CRIT)
      write (msg, '(a,i6,a,i2,4(a,i5))') 'prop_ls:', id, ' part', ipart, ' tile', its, ':', ite, ',', jts, ':', jte
!$OMP END CRITICAL(SFIRE_CORE_CRIT)
      call message(msg)

      if (ipart == 1) then

         call check_lfn_tign('prop_ls start', ts, its, ite, jts, jte, ims, ime, jms, jme, lfn_in, tign)

         a = fire_back_weight
         a1 = 1.-a

         ihs2 = max(its - 2, ids)
         ihe2 = min(ite + 2, ide)
         jhs2 = max(jts - 2, jds)
         jhe2 = min(jte + 2, jde)

         ihs = max(its - 1, ids)
         ihe = min(ite + 1, ide)
         jhs = max(jts - 1, jds)
         jhe = min(jte + 1, jde)

         call write_array_m(ihs, ihe, jhs, jhe, ims, ime, jms, jme, lfn_in, 'lfn_in', id)

         call check_mesh_2dim(ihs2, ihe2, jhs2, jhe2, ims, ime, jms, jme)
         call print_2d_stats(ihs2, ihe2, jhs2, jhe2, ims, ime, jms, jme, lfn_in, 'prop_ls: lfn in')

         id1 = id
         if (id1 .ne. 0) id1 = id1 + 1000
         call tend_ls(ifun, id1, &
                      ims, ime, jms, jme, &
                      its - 1, ite + 1, jts - 1, jte + 1, &
                      ids, ide, jds, jde, &
                      ips, ipe, jps, jpe, &
                      ihs, ihe, jhs, jhe, &
                      ims, ime, jms, jme, &
                      its, ite, jts, jte, &
                      ts, dt, dx, dy, &
                      lfn_in, &
                      tbound, &
                      tend, ros, &
                      fp, config_flags &
                      ) !config_flags e ifun INTROD ISILDA CUNHA MENEZES
         !   print*,"estou no core na prop_ls, fora da tend_ls_1, vou para write_array_m"
         !   call flush(6)
         !print*,"MHHHH tend",tend
         !call flush(6)
         !print*,"MEEEEE lfn_in",lfn_in
         !call flush(6)
         !call flush(6)
         call write_array_m(ihs, ihe, jhs, jhe, its - 1, ite + 1, jts - 1, jte + 1, tend, 'tend1', id)

         do j = jhs, jhe
            do i = ihs, ihe
               lfn1(i, j) = lfn_in(i, j) + dt*tend(i, j)
               !if (lfn1(i,j) .gt. 10000000000000000000000.0) then
               ! print*,"FFFFFFPPPP lfn_in",lfn_in(i,j),"dt",dt,"tend",tend(i,j),"i",i,"j",j
               ! call flush(6)
               ! endif
            end do
         end do

         !call print_2d_stats(ihs, ihe, jhs, jhe, its - 1, ite + 1, jts - 1, jte + 1, lfn1, 'prop_ls: lfn1')

         !  print*,"estou no core na prop_ls, vou entrar outr vez na tend_ls agora a 2"
         !  call flush(6)
         if (id1 .ne. 0) id1 = id1 + 1000
         call tend_ls(ifun, id1, &
                      its - 1, ite + 1, jts - 1, jte + 1, &
                      its - 1, ite + 1, jts - 1, jte + 1, &
                      ids, ide, jds, jde, &
                      ips, ipe, jps, jpe, &
                      its, ite, jts, jte, &
                      ims, ime, jms, jme, &
                      its, ite, jts, jte, &
                      ts + dt, dt, dx, dy, &
                      lfn1, &
                      tbound2, &
                      tend, ros, &
                      fp, config_flags &
                      ) !config_flags e ifun INTROD ISILDA CUNHA MENEZES
         ! print*,"estou no core na prop_ls, fora da tend_ls_2, vou para write_array_m"
         call write_array_m(its, ite, jts, jte, its - 1, ite + 1, jts - 1, jte + 1, tend, 'tend2', id)
         ! print*,"estou no core, vou para print_2d_stats"
         !LFR call print_2d_stats(its, ite, jts, jte, ims, ime, jms, jme, ros, 'prop_ls: ros')
         !LFR call print_2d_stats(its, ite, jts, jte, its - 1, ite + 1, jts - 1, jte + 1, tend, 'prop_ls: tend2')
         ! print*,"estou no core na prop_ls, vou gerar o tbound"
         tbound = min(tbound, tbound2)

!$OMP CRITICAL(SFIRE_CORE_CRIT)
         write (msg, '(a,f10.2,4(a,f7.2))') 'prop_ls: time', ts, ' dt=', dt, ' bound', min(tbound, 999.99), &
            ' dx=', dx, ' dy=', dy
!$OMP END CRITICAL(SFIRE_CORE_CRIT)
         call message(msg)
         if (dt > tbound) then
!$OMP CRITICAL(SFIRE_CORE_CRIT)
            write (msg, '(2(a,f10.2))') 'prop_ls: WARNING: time step ', dt, &
               ' > bound =', tbound
!$OMP END CRITICAL(SFIRE_CORE_CRIT)
            call message(msg)
         end if

         !  print*,"estou no core na prop_ls, vou calcular lfn_out"
         !  call flush(6)
         do j = jts, jte
            do i = its, ite
               ! print*,"membros do lfn_out1",j,i,a1,a, lfn1(i,j), lfn_in(i,j),dt,tend(i,j),a1*lfn1(i,j), a*lfn_in(i,j), dt*tend(i,j)
               ! call flush(6)
               !  print*,"m1",  lfn1(i,j); call flush(6)
               !  print*,"m2",  lfn_in(i,j) ; call flush(6)
               !  print*,"m3", tend(i,j); call flush(6)

               lfn_out(i, j) = a1*lfn1(i, j) + a*(lfn_in(i, j) + dt*tend(i, j))
               ! print*,"minimo do lfn_out1", lfn_out(i,j)
               ! call flush(6)
               lfn_out(i, j) = min(lfn_out(i, j), lfn_in(i, j))
               ! print*,"minimo do lfn_out1A", lfn_out(i,j)
               ! call flush(6)
            end do
         end do
         ! print*,"estou no core na prop_ls, vou calcular depois do lfn_out, vou sair do if ipart==1"
         ! call flush(6)
      elseif (ipart == 2) then
         ! print*,"estou no core na prop_ls, vou para o continue_at_boundary"
         ! call flush(6)
         call continue_at_boundary(1, 1, zero, &
                                   ims, ime, jms, jme, &
                                   ids, ide, jds, jde, &
                                   ips, ipe, jps, jpe, &
                                   its, ite, jts, jte, &
                                   itso, iteo, jtso, jteo, &
                                   lfn_out)

         ! print*,"estou no core na prop_ls, vou calcular lfn1"
         ! call flush(6)
         do j = jts, jte
            do i = its, ite
               lfn1(i, j) = lfn_out(i, j)
            end do
         end do

         ! print*,"estou no core na prop_ls, vou calcular t0 t1"
         ! call flush(6)
         do j = jts, jte
            do i = its, ite
               t0 = min(lfn_in(i + 1, j), lfn_in(i - 1, j), lfn_in(i, j + 1), lfn_in(i, j - 1))
               t1 = min(lfn1(i + 1, j), lfn1(i - 1, j), lfn1(i, j + 1), lfn1(i, j - 1))
               if (.not. t0 > lfn_in(i, j) .and. t1 > lfn_out(i, j)) then
                  if (fire_print_msg > 2) then
!$OMP CRITICAL(SFIRE_CORE_CRIT)
                     write (msg, '(a,2i6)') 'prop_ls: new local min', i, j
                     call message(msg)
                     write (msg, '((a,g13.5,1x))') &
                        'prop_ls: new local min', lfn_out(i, j), 'fixed to', t1, 'incr', t1 - lfn_in(i, j)
                     call message(msg)
!$OMP END CRITICAL(SFIRE_CORE_CRIT)
                  end if
                  lfn_out(i, j) = t1
               end if
            end do
         end do

         ! print*,"estou no core na prop_ls, vou calcular time"
         ! call flush(6)
         time_now = ts + dt
         time_future = ts + 2*dt
         do j = jts, jte
            do i = its, ite
               call check_lfn_tign_ij(i, j, 'prop_ls before', ts, lfn_in(i, j), tign(i, j))
               if (lfn_out(i, j) < 0.) then
                  if (.not. lfn_in(i, j) < 0) then
                     tign(i, j) = time_now + dt*lfn_out(i, j)/(lfn_in(i, j) - lfn_out(i, j))
                     !   print*,"GGGLLLL tign",tign(i,j)
                     !   call flush(6)
                  end if
               else
                  tign(i, j) = time_future
               end if
               call check_lfn_tign_ij(i, j, 'prop_ls after', time_now, lfn_out(i, j), tign(i, j))

            end do
         end do

         nfirenodes = 0
         nfireline = 0
         sum_err = 0.
         min_err = rmax
         max_err = rmin
         sum_aerr = 0.
         min_aerr = rmax
         max_aerr = rmin
         its1 = its + 1
         jts1 = jts + 1
         ite1 = ite - 1
         jte1 = jte - 1

         ! print*,"estou no core na prop_ls, vou calcular grads asleed rr sum etc"
         ! call flush(6)
         do j = jts1, jte1
            do i = its1, ite1
               if (lfn_out(i, j) > 0.0) then
                  if (lfn_out(i + 1, j) <= 0 .or. lfn_out(i, j + 1) <= 0 .or. &
                      lfn_out(i - 1, j) <= 0 .or. lfn_out(i, j - 1) <= 0) then
                     gradx = (lfn_out(i + 1, j) - lfn_out(i - 1, j))/(2.0*dx)
                     grady = (lfn_out(i, j + 1) - lfn_out(i, j - 1))/(2.0*dy)
                     grad2 = sqrt(gradx*gradx + grady*grady)
                     aspeed = (lfn_in(i, j) - lfn_out(i, j))/(dt*max(grad2, rmin))
                     !   print*,"RRRPPP aspeed",aspeed
                     !   call flush(6)
                     rr = speed_func(ifun, gradx, grady, dx, dy, i, j, fp, ierrx, msg, &
                                     config_flags) !config_flags e ifun INTRODUZIDO POR ISILDA CM
                     err = aspeed - rr
                     sum_err = sum_err + err
                     min_err = min(min_err, err)
                     max_err = max(max_err, err)
                     aerr = abs(err)
                     sum_aerr = sum_aerr + aerr
                     min_aerr = min(min_aerr, aerr)
                     max_aerr = max(max_aerr, aerr)
                     nfireline = nfireline + 1
                  end if
               else
                  nfirenodes = nfirenodes + 1
               end if
            end do
         end do
!$OMP CRITICAL(SFIRE_CORE_CRIT)
         write (msg, '(2(a,i6,f8.4))') 'prop_ls: nodes burning', nfirenodes, &
            (100.*nfirenodes)/((ite1 - its1 + 1)*(jte1 - jts1 + 1)), '% next to fireline', nfireline
!$OMP END CRITICAL(SFIRE_CORE_CRIT)
         call message(msg)

         !print*,"estou no core na pro_ls, vou calcular boundary_guard"
         ! call flush(6)
         do k = -1, 1, 2

            do kk = 1, boundary_guard
               i = ids + k*kk
               if (i .ge. its .and. i .le. ite) then
                  do j = jts, jte
                     if (lfn_out(i, j) <= 0.) goto 9
                  end do
               end if
            end do

            do kk = 1, boundary_guard
               j = jds + k*kk
               if (j .ge. jts .and. j .le. jte) then
                  do i = its, ite
                     if (lfn_out(i, j) <= 0.) goto 9
                  end do
               end if
            end do
         end do
         goto 10
9        continue
!$OMP CRITICAL(SFIRE_CORE_CRIT)
         write (msg, '(a,i2,a,2i8)') 'prop_ls: fire', boundary_guard, &
            ' cells from domain boundary at node ', i, j
!$OMP END CRITICAL(SFIRE_CORE_CRIT)
         call message(msg)
         call crash('prop_ls: increase the fire region')
10       continue

         call print_2d_stats(its, ite, jts, jte, ims, ime, jms, jme, lfn_out, 'prop_ls: lfn out')
         call print_2d_stats(its, ite, jts, jte, ims, ime, jms, jme, tign, 'prop_ls: tign out')

         call check_lfn_tign('prop_ls end', ts + dt, its, ite, jts, jte, ims, ime, jms, jme, lfn_out, tign)

      else
         call crash('prop_ls: ipart must be 1 or 2')
      end if

      !print *, "estou no core vou sair prop_ls"; call flush (6)
   end subroutine prop_ls

!!!!! COMENTADO POR ISILDA CUNHA MENEZES
!subroutine tend_ls( id, &
!    lims,lime,ljms,ljme, &
!    tims,time,tjms,tjme, &
!    ids,ide,jds,jde, &
!    ips,ipe,jps,jpe, &
!    ints,inte,jnts,jnte, &
!    ims,ime,jms,jme, &
!    its,ite,jts,jte, &
!    t,dt,dx,dy,      &
!    lfn, &
!    tbound, &
!    tend, ros,  &
!    fp &
!)

!implicit none

!integer,intent(in)::id,lims,lime,ljms,ljme,tims,time,tjms,tjme
!integer,intent(in)::ims,ime,jms,jme,its,ite,jts,jte
!integer, intent(in)::ids,ide,jds,jde,ints,inte,jnts,jnte,ips,ipe,jps,jpe
!real,intent(in)::t
!real,intent(in)::dt,dx,dy
!real,dimension(lims:lime,ljms:ljme),intent(inout)::lfn
!real,intent(out)::tbound
!real,dimension(tims:time,tjms:tjme),intent(out)::tend
!real,dimension(ims:ime,jms:jme),intent(out)::ros
!type(fire_params),intent(in)::fp

!real:: te,diffLx,diffLy,diffRx,diffRy, &
!   diffCx,diffCy,diff2x,diff2y,grad,rr, &
!   ros_back,ros_wind,ros_slope,advx,advy,scale,nvx,nvy,speed
!integer::i,j,itso,iteo,jtso,jteo,ierrx,nerr
!integer::upwind_err=0
!character(len=128)msg,msg2

!real, parameter:: eps=epsilon(0.0)

!real, parameter:: zero=0.,one=1.,tol=100*eps, &
!    safe=2.,rmin=safe*tiny(zero),rmax=huge(zero)/safe!

!intrinsic max,min,sqrt,nint,tiny,huge

!    call check_mesh_2dim(ints-1,inte+1,jnts-1,jnte+1,lims,lime,ljms,ljme)
!    call check_mesh_2dim(ints,inte,jnts,jnte,tims,time,tjms,tjme)

!    call continue_at_boundary(1,1,fire_lfn_ext_up, &
!    lims,lime,ljms,ljme, &
!    ids,ide,jds,jde, &
!    ips,ipe,jps,jpe, &
!    ints,inte,jnts,jnte, &
!    itso,iteo,jtso,jteo, &
!    lfn)

!    call print_2d_stats(itso,iteo,jtso,jteo,lims,lime,ljms,ljme, &
!                   lfn,'tend_ls: lfn cont')

!    call write_array_m(ints-1,inte+1,jnts-1,jnte+1,lims,lime,ljms,ljme,lfn,'tend_lfn_in',id)

!    nerr=0
!    tbound=0
!    do j=jnts,jnte
!        do i=ints,inte

!            diffRx = (lfn(i+1,j)-lfn(i,j))/dx
!            diffLx = (lfn(i,j)-lfn(i-1,j))/dx
!            diffRy = (lfn(i,j+1)-lfn(i,j))/dy
!            diffLy = (lfn(i,j)-lfn(i,j-1))/dy
!            diffCx = diffLx+diffRx
!            diffCy = diffLy+diffRy

!            select case(sfire_upwinding)
!            case(0)
!                grad=sqrt(diffCx**2 + diffCy**2)
!            case(1)
!                diff2x=select_upwind(diffLx,diffRx)
!                diff2y=select_upwind(diffLy,diffRy)
!                grad=sqrt(diff2x*diff2x + diff2y*diff2y)
!            case(2)
!                diff2x=select_godunov(diffLx,diffRx)
!                diff2y=select_godunov(diffLy,diffRy)
!                grad=sqrt(diff2x*diff2x + diff2y*diff2y)
!            case(3)
!                diff2x=select_eno(diffLx,diffRx)
!                diff2y=select_eno(diffLy,diffRy)
!                grad=sqrt(diff2x*diff2x + diff2y*diff2y)
!            case(4)
!                grad=sqrt(max(diffLx,0.)**2+min(diffRx,0.)**2   &
!                        + max(diffLy,0.)**2+min(diffRy,0.)**2)
!            case default
!                grad=0.
!                upwind_err=1
!            end select

!            scale=sqrt(diffCx*diffCx+diffCy*diffCy+eps)
!            nvx=diffCx/scale
!            nvy=diffCy/scale

!            call fire_ros(ros_back,ros_wind,ros_slope, &
!            nvx,nvy,i,j,fp,ierrx,msg2)
!            nerr = nerr + ierrx

!            rr=ros_back + ros_wind + ros_slope
!            if(fire_grows_only.gt.0)rr=max(rr,0.)

!            if(i.ge.its.and.i.le.ite.and.j.ge.jts.and.j.le.jte)ros(i,j)=rr

!            te = -rr*grad

!            if (grad > 0.) then
!                 tbound = max(tbound,rr*(abs(diff2x)/dx+abs(diff2y)/dy)/grad)
!            endif

!            te=te + fire_viscosity*abs(rr)*((diffRx-diffLx)+(diffRy-diffLy))

!            tend(i,j)=te
!        enddo
!    enddo

!    if(upwind_err>0)call crash('sfire_upwinding value not supported, 3 recommended')

!    if(nerr>0)then
!!$OMP CRITICAL(SFIRE_CORE_CRIT)
!        write(msg,'(a,i6,1x,a)')'tend_ls:',nerr,'messages in rate of spread computations. Last message:'
!!$OMP END CRITICAL(SFIRE_CORE_CRIT)
!        call warning(msg)
!        call warning(msg2)
!    endif

!    call print_2d_stats(its,ite,jts,jte,ims,ime,jms,jme, ros,'tend_ls: ros')
!    call print_2d_stats(ints,inte,jnts,jnte,tims,time,tjms,tjme, &
!                   tend,'tend_ls: tend out')

!    tbound = 1/(tbound+tol)

!end subroutine tend_ls
!!!!!!!

!!!!! #### ALTERADA E INTRODUZIDO POR ISILDA CUNHA MENEZES
!!!CALCULA A PROPAGACAO ASSOCIANDO A VELOCIDADE TOTAL SE O MODELO COPA ESTIVER ACTIVO OU PROPAGACAO DE SUPERFICIE COM O METODO DE PROPAGACAO NUMERICO SET LEVEL

   subroutine tend_ls(ifun, id, &
                      lims, lime, ljms, ljme, &
                      tims, time, tjms, tjme, &
                      ids, ide, jds, jde, &
                      ips, ipe, jps, jpe, &
                      ints, inte, jnts, jnte, &
                      ims, ime, jms, jme, &
                      its, ite, jts, jte, &
                      t, dt, dx, dy, &
                      lfn, &
                      tbound, &
                      tend, ros, &
                      fp, config_flags &
                      ) !config_flags e ifun INTROD ISILDA CUNHA MENEZES

      implicit none

      integer, intent(in)::id, lims, lime, ljms, ljme, tims, time, tjms, tjme
      integer, intent(in)::ims, ime, jms, jme, its, ite, jts, jte
      integer, intent(in)::ids, ide, jds, jde, ints, inte, jnts, jnte, ips, ipe, jps, jpe
      real, intent(in)::t
      real, intent(in)::dt, dx, dy
      real, dimension(lims:lime, ljms:ljme), intent(inout)::lfn
      real, intent(out)::tbound
      real, dimension(tims:time, tjms:tjme), intent(out)::tend
      real, dimension(ims:ime, jms:jme), intent(out)::ros

!type(fire_params),intent(in)::fp
!!!!INTRODUZIDO POR ISILDA CUNHA MENEZES
      type(fire_params), intent(inout)::fp
      real :: ros_total
      type(grid_config_rec_type), intent(IN)  :: config_flags
      logical:: crown
      integer, intent(in)::ifun
!!!!!!!!!

      real:: te, diffLx, diffLy, diffRx, diffRy, &
             diffCx, diffCy, diff2x, diff2y, grad, rr, &
             ros_back, ros_wind, ros_slope, advx, advy, scale, nvx, nvy, speed
      integer::i, j, itso, iteo, jtso, jteo, ierrx, nerr
      integer::upwind_err = 0
      character(len=128) msg, msg2

      real, parameter:: eps = epsilon(0.0)

      real, parameter:: zero = 0., one = 1., tol = 100*eps, &
                        safe = 2., rmin = safe*tiny(zero), rmax = huge(zero)/safe

      intrinsic max, min, sqrt, nint, tiny, huge

      !   print*,"estou no module_fr_sfire_core dentro da tend_ls"
      !   call flush(6)
      crown = config_flags%crown ! Introd por ISILDA DA CUNHA MENEZES

      call check_mesh_2dim(ints - 1, inte + 1, jnts - 1, jnte + 1, lims, lime, ljms, ljme)
      call check_mesh_2dim(ints, inte, jnts, jnte, tims, time, tjms, tjme)

      call continue_at_boundary(1, 1, fire_lfn_ext_up, &
                                lims, lime, ljms, ljme, &
                                ids, ide, jds, jde, &
                                ips, ipe, jps, jpe, &
                                ints, inte, jnts, jnte, &
                                itso, iteo, jtso, jteo, &
                                lfn)

      call print_2d_stats(itso, iteo, jtso, jteo, lims, lime, ljms, ljme, &
                          lfn, 'tend_ls: lfn cont')

      call write_array_m(ints - 1, inte + 1, jnts - 1, jnte + 1, lims, lime, ljms, ljme, lfn, 'tend_lfn_in', id)

      nerr = 0
      tbound = 0
      do j = jnts, jnte
         do i = ints, inte

            diffRx = (lfn(i + 1, j) - lfn(i, j))/dx
            diffLx = (lfn(i, j) - lfn(i - 1, j))/dx
            diffRy = (lfn(i, j + 1) - lfn(i, j))/dy
            diffLy = (lfn(i, j) - lfn(i, j - 1))/dy
            diffCx = diffLx + diffRx
            diffCy = diffLy + diffRy

            select case (sfire_upwinding)
            case (0)
               grad = sqrt(diffCx**2 + diffCy**2)
            case (1)
               diff2x = select_upwind(diffLx, diffRx)
               diff2y = select_upwind(diffLy, diffRy)
               grad = sqrt(diff2x*diff2x + diff2y*diff2y)
            case (2)
               diff2x = select_godunov(diffLx, diffRx)
               diff2y = select_godunov(diffLy, diffRy)
               grad = sqrt(diff2x*diff2x + diff2y*diff2y)
            case (3)
               diff2x = select_eno(diffLx, diffRx)
               diff2y = select_eno(diffLy, diffRy)
               grad = sqrt(diff2x*diff2x + diff2y*diff2y)
            case (4)
               grad = sqrt(max(diffLx, 0.)**2 + min(diffRx, 0.)**2 &
                           + max(diffLy, 0.)**2 + min(diffRy, 0.)**2)
            case default
               grad = 0.
               upwind_err = 1
            end select

            scale = sqrt(diffCx*diffCx + diffCy*diffCy + eps)
            nvx = diffCx/scale
            nvy = diffCy/scale

!!!!######### chama R_total - velocidade de Rothermel de superficie e de copa

            if (crown) then
               !    print*,"estou no module_fr_sfire_core na tend_ls vou chamar fire_ros com crown"
               !    call flush(6)
               call fire_ros(ifun, ros_back, ros_wind, ros_slope, &
                             nvx, nvy, i, j, fp, ierrx, msg2, config_flags, ros_total) !config_flags intorduzido por ISILDA CM

               rr = ros_total
               !   print*, "COPA_BBB estou no module_fr_sfire_core na tend_ls entreguei ros_total",rr
               !   call flush(6)

            else !(crown .eq. .FALSE.) then
            !!! INTRODUZIDO POR ISILDA CUNHA MENEZES
               !  print*,"estou no module_fr_sfire_core na tend_ls vou chamar fire_ros sem crown"
               !  call flush(6)
               call fire_ros(ifun, ros_back, ros_wind, ros_slope, &
                             nvx, nvy, i, j, fp, ierrx, msg2, config_flags) !config_flags introduzido por ISILDA CM
            !!!!
               rr = ros_back + ros_wind + ros_slope

            end if

            nerr = nerr + ierrx

            if (fire_grows_only .gt. 0) rr = max(rr, 0.)

            if (i .ge. its .and. i .le. ite .and. j .ge. jts .and. j .le. jte) ros(i, j) = rr

            !  print*,"KKKKNNN SURF rr",rr,"grad",grad
            !  call flush(6)
            te = -rr*grad

            if (grad > 0.) then
               tbound = max(tbound, rr*(abs(diff2x)/dx + abs(diff2y)/dy)/grad)
               !    print*,'tbound',tbound
               !     call flush(6)
            end if

            !   print*,"TTTHHH1 SURF fire_viscosity",fire_viscosity,"abs(rr)",abs(rr)
            !    call flush(6)
            !   print*,"TTTHHH2 diffRx",diffRx,"diffLx",diffLx,"diffRy",diffRy,"diffLy",diffLy
            !     call flush(6)
            te = te + fire_viscosity*abs(rr)*((diffRx - diffLx) + (diffRy - diffLy))
            !  print*,"MMMVVV SURF tend",te,"i=",i,"j=",j
            !   call flush(6)
            tend(i, j) = te

         end do
      end do

      if (upwind_err > 0) call crash('sfire_upwinding value not supported, 3 recommended')

      if (nerr > 0) then
!$OMP CRITICAL(SFIRE_CORE_CRIT)
         write (msg, '(a,i6,1x,a)') 'tend_ls:', nerr, 'messages in rate of spread computations. Last message:'
!$OMP END CRITICAL(SFIRE_CORE_CRIT)
         call warning(msg)
         call warning(msg2)
      end if

      !LFR call print_2d_stats(its, ite, jts, jte, ims, ime, jms, jme, ros, 'tend_ls: ros')
      !LFR call print_2d_stats(ints, inte, jnts, jnte, tims, time, tjms, tjme, &
      !LFR                    tend, 'tend_ls: tend out')

      tbound = 1/(tbound + tol)

      !  print*,"vou sair da tend_ls"
      !  call flush(6)

   end subroutine tend_ls

!!!!

   real function select_upwind(diffLx, diffRx)
      implicit none
      real, intent(in):: diffLx, diffRx
      real diff2x

      diff2x = 0
      if (diffLx > 0 .and. diffRx > 0.) diff2x = diffLx
      if (diffLx < 0 .and. diffRx < 0.) diff2x = diffRx

      select_upwind = diff2x
   end function select_upwind

   real function select_godunov(diffLx, diffRx)
      implicit none
      real, intent(in):: diffLx, diffRx
      real diff2x, diffCx

      diff2x = 0
      diffCx = diffRx + diffLx
      if (diffLx > 0 .and. .not. diffCx < 0) diff2x = diffLx
      if (diffRx < 0 .and. diffCx < 0) diff2x = diffRx

      select_godunov = diff2x
   end function select_godunov

   real function select_eno(diffLx, diffRx)
      implicit none
      real, intent(in):: diffLx, diffRx
      real diff2x

      if (.not. diffLx > 0 .and. .not. diffRx > 0) then
         diff2x = diffRx
      elseif (.not. diffLx < 0 .and. .not. diffRx < 0) then
         diff2x = diffLx
      elseif (.not. diffLx < 0 .and. .not. diffRx > 0) then
         if (.not. abs(diffRx) < abs(diffLx)) then
            diff2x = diffRx
         else
            diff2x = diffLx
         end if
      else
         diff2x = 0.
      end if

      select_eno = diff2x
   end function select_eno

   real function speed_func(ifun, diffCx, diffCy, dx, dy, i, j, fp, ierrx, msg, config_flags) !ifun INTRODUZIDO POR ISILDA

      implicit none

      real, intent(in)::diffCx, diffCy
      real, intent(in)::dx, dy
      integer, intent(in)::i, j
!!! COMENTADO POR ISILDA CUNHA MENEZES
!type(fire_params),intent(in)::fp
!!! INTRODUZIDO POR ISILDA CUNHA MENEZES
      type(fire_params), intent(inout)::fp
!!!!!
      integer, intent(out)::ierrx
      character(len=*), intent(out)::msg

      real::scale, nvx, nvy, r
      real::ros_back, ros_wind, ros_slope
      real:: ros_total !!! INTRODUZIDO POR ISILDA CUNHA MENEZES
      type(grid_config_rec_type) :: config_flags !!! INTRODUZIDO POR ISILDA CUNHA MENEZES
      real, parameter:: eps = epsilon(0.0)
      integer, intent(in)::ifun ! INTRODUZIDO POR ISILDA CUNHA MENEZES

      scale = sqrt(diffCx*diffCx + diffCy*diffCy + eps)
      nvx = diffCx/scale
      nvy = diffCy/scale

            !!!! INTRODUZIDO POR ISILDA CUNHA MENEZES
      call fire_ros(ifun, ros_back, ros_wind, ros_slope, &
                    nvx, nvy, i, j, fp, ierrx, msg, config_flags, ros_total) !config_flagsintroduzido por ISILDA CM
      r = ros_total
            !!!!!!!!

      ! call fire_ros(ros_back,ros_wind,ros_slope, &
      ! nvx,nvy,i,j,fp,ierrx,msg)

      ! r=ros_back + ros_wind + ros_slope
      if (fire_grows_only .gt. 0) r = max(r, 0.)
      speed_func = r

   end function speed_func

end module module_fr_sfire_core
