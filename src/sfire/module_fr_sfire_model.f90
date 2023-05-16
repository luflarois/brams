
module module_fr_sfire_model

   use module_fr_sfire_core
   use module_fr_sfire_util
   use module_fr_sfire_phys

!!!!! INCORPORADO POR ISILDA CUNHA MENEZES
!use module_fr_sfire_crown
   use ModNamelistsfireFile, only: grid_config_rec_type
!!!!!

   implicit none

contains

   subroutine sfire_model( &
      id, &
      ifun, &
      restart, replay, &
      run_fuel_moisture, &
      ifuelread, nfuel_cat0, &
      ifds, ifde, jfds, jfde, &
      ifms, ifme, jfms, jfme, &
      ifps, ifpe, jfps, jfpe, &
      ifts, ifte, jfts, jfte, &
      time_start, dt, &
      fdx, fdy, &
      ignition, hfx, &
      coord_xf, coord_yf, &
      fire_hfx, &
      tign_in, &
      lfn, lfn_out, tign, fuel_frac, fire_area, &
      fuel_frac_burnt, &
      fgrnhfx, fgrnqfx, &
      ros, flineint, flineint2, &
      f_ros0, f_rosx, f_rosy, f_ros, &
      f_int, f_lineint, f_lineint2, &
      nfuel_cat, &
      fuel_time, fwh, fz0, &
      fp, &
      flineint_total, config_flags) !!!! flineint_total,config_flags introduzido por ISILDA CUNHA MENEZES

      implicit none

      integer, intent(in) :: id
      integer, intent(in) :: ifun

      logical, intent(in):: restart
      logical, intent(in):: replay
      logical, intent(in)::run_fuel_moisture

      integer, intent(in) :: ifuelread, nfuel_cat0
      integer, intent(in) :: ifds, ifde, jfds, jfde, &
                             ifps, ifpe, jfps, jfpe
      integer, intent(in) :: ifts, ifte, jfts, jfte
      integer, intent(in) :: ifms, ifme, jfms, jfme
      real, intent(in) :: time_start, dt
      real, intent(in) :: fdx, fdy

      type(lines_type), intent(in):: ignition, hfx
      real, dimension(ifms:ifme, jfms:jfme), intent(in):: &
         coord_xf, coord_yf
      real, dimension(ifms:ifme, jfms:jfme), intent(inout):: &
         fire_hfx

      real, dimension(ifms:ifme, jfms:jfme), intent(inout):: tign_in

      real, intent(inout), dimension(ifms:ifme, jfms:jfme):: &
         lfn, &
         tign, &
         fuel_frac

      real, intent(out), dimension(ifms:ifme, jfms:jfme):: &
         fire_area, &
         fuel_frac_burnt

      real, intent(out), dimension(ifms:ifme, jfms:jfme):: &
         lfn_out, &
         fgrnhfx, fgrnqfx, &
         ros, flineint, flineint2, &
         f_ros0, f_rosx, f_rosy, f_ros, &
         f_int, f_lineint, f_lineint2, flineint_total !!!! flineint_total introduzido por ISILDA CUNHA MENEZES

      real, intent(inout), dimension(ifms:ifme, jfms:jfme)::nfuel_cat
      real, intent(inout), dimension(ifms:ifme, jfms:jfme):: fuel_time, fwh, fz0
      type(fire_params), intent(inout)::fp

      integer :: xifms, xifme, xjfms, xjfme
      real, dimension(ifts:ifte, jfts:jfte)::fuel_frac_end
      integer::ignited, ig, i, j, itso, iteo, jtso, jteo
      real::tbound, err, erri, errj, maxgrad, grad, tfa, thf, mhf, tqf, mqf, aw, mw, t
      character(len=128)::msg
      logical:: freeze_fire
      real::fireline_mask = 0.

!#### INTRODUZIDO POR ISILDA CUNHA MENEZES
      type(grid_config_rec_type), intent(IN)  :: config_flags
      logical :: crown
!####

      call check_mesh_2dim(ifts - 1, ifte + 1, jfts - 1, jfte + 1, ifms, ifme, jfms, jfme)

      xifms = ifms
      xifme = ifme
      xjfms = jfms
      xjfme = jfme

      freeze_fire = fire_hfx_given .ne. 0
      !print *, 'LFR-DBG: IFUN:', ifun, " => ", restart, fire_tign_in_time, tiny(dt)

      if (ifun .eq. 1) then
! !$OMP SINGLE

! !$OMP END SINGLE
         if (replay) then

            call message('replay, setting the level set function', level=1)
            do j = jfts, jfte
               do i = ifts, ifte
                  lfn(i, j) = tign(i, j) - time_start
               end do
            end do
            call check_lfn_tign('replay, lfn from tign', time_start, ifts, ifte, jfts, jfte, ifms, ifme, jfms, jfme, lfn, tign)
         end if
      elseif (ifun .eq. 2) then

         !print *, "estou dentro ifun=2 na rotina sfire_model e vou entrar na continue_at_boundary"
         !call flush (6)
         call continue_at_boundary(1, 1, 0., &
                                   ifms, ifme, jfms, jfme, &
                                   ifds, ifde, jfds, jfde, &
                                   ifps, ifpe, jfps, jfpe, &
                                   ifts, ifte, jfts, jfte, &
                                   itso, iteo, jtso, jteo, &
                                   fp%zsf)

         err = 0.
         maxgrad = 0.
         do j = jfts, jfte
            do i = ifts, ifte
               erri = fp%dzdxf(i, j) - (fp%zsf(i + 1, j) - fp%zsf(i - 1, j))/(2.*fdx)
               errj = fp%dzdyf(i, j) - (fp%zsf(i, j + 1) - fp%zsf(i, j - 1))/(2.*fdy)
               err = max(err, abs(erri), abs(errj))
               grad = sqrt(fp%dzdxf(i, j)**2 + fp%dzdyf(i, j)**2)
               maxgrad = max(maxgrad, grad)
            end do
         end do
!$OMP CRITICAL(SFIRE_MODEL_CRIT)
         write (msg, *) 'max gradient ', maxgrad, ' max error against zsf', err
!$OMP END CRITICAL(SFIRE_MODEL_CRIT)
         call message(msg)
         !print *, "estou dentro ifun=2 na rotina sfire_model e vou entrar na set_nfuel_cat"
         !call flush (6)
         call set_nfuel_cat( &
            ifms, ifme, jfms, jfme, &
            ifts, ifte, jfts, jfte, &
            ifuelread, nfuel_cat0, &
            fp%zsf, nfuel_cat)

         !print *, "estou dentro ifun=2 na rotina sfire_model e vou entrar na set_fire_params"
         !call flush (6)
         call set_fire_params( &
            ifds, ifde, jfds, jfde, &
            ifms, ifme, jfms, jfme, &
            ifts, ifte, jfts, jfte, &
            fdx, fdy, nfuel_cat0, &
            nfuel_cat, fuel_time, &
            fp)

!!!!! AUTOR ISILDA DA CUNHA MENEZES
         !print *, "estou dentro ifun=2 na rotina sfire_model e vou entrar na set_fire_crown_params, crow=",crown
         !call flush (6)
         if (crown) then
            !print *, "estou no module_fire_model ifun= 2 e vou entrar na set_fire_crown_params"
            !call flush (6)
            call set_fire_crown_params( &
               ifds, ifde, jfds, jfde, &
               ifms, ifme, jfms, jfme, &
               ifts, ifte, jfts, jfte, &
               fdx, fdy, nfuel_cat0, &
               nfuel_cat, fuel_time, &
               fp)
         end if
!!!!!!!
         !print *, "estou dentro ifun=2 na rotina sfire_model e vou entrar na init_no_fire"
         !call flush (6)
         if (.not. restart) then
            call init_no_fire( &
               ifds, ifde, jfds, jfde, &
               ifms, ifme, jfms, jfme, &
               ifts, ifte, jfts, jfte, &
               fdx, fdy, time_start, dt, &
               fuel_frac, fire_area, lfn, tign_in, tign)

         end if

         call check_lfn_tign('after initialization:', time_start, ifts, ifte, jfts, jfte, ifms, ifme, jfms, jfme, lfn, tign)
         !print *, "estou dentro ifun=2 na rotina sfire_model e vou entrar na fuel_left"
         !call flush (6)
         if (replay) then
            call message('replay, recomputing fuel fraction')
            call fuel_left( &
               ifds, ifde, jfds, jfde, &
               ifms, ifme, jfms, jfme, &
               ifts, ifte, jfts, jfte, &
               ifms, ifme, jfms, jfme, &
               lfn, tign, fuel_time, time_start, fuel_frac, fire_area)
            call check_lfn_tign('recomputed fuel fraction', time_start, ifts, ifte, jfts, jfte, ifms, ifme, jfms, jfme, lfn, tign)
         end if

      elseif (ifun .eq. 3) then

      elseif (ifun .eq. 4) then

         if (run_fuel_moisture) then

            !print *, "estou dentro ifun=4 na rotina sfire_model e vou entrar na set_fire_params com o fuel_moisture"
            !call flush (6)
            call set_fire_params( &
               ifds, ifde, jfds, jfde, &
               ifms, ifme, jfms, jfme, &
               ifts, ifte, jfts, jfte, &
               fdx, fdy, nfuel_cat0, &
               nfuel_cat, fuel_time, &
               fp)
         end if

!!!!! AUTOR ISILDA DA CUNHA MENEZES
         if (crown) then
            !print *, "estou dentro ifun=4 no module_fire_model ifun= 4 e vou entrar na set_fire_crown_params"
            !call flush (6)
            call set_fire_crown_params( &
               ifds, ifde, jfds, jfde, &
               ifms, ifme, jfms, jfme, &
               ifts, ifte, jfts, jfte, &
               fdx, fdy, nfuel_cat0, &
               nfuel_cat, fuel_time, &
               fp)
         end if
!!!!!!!

         if (fire_print_msg .ge. stat_lev) then
            aw = fun_real(RNRM_SUM, &
                          ifms, ifme, 1, 1, jfms, jfme, &
                          ifds, ifde, 1, 1, jfds, jfde, &
                          ifts, ifte, 1, 1, jfts, jfte, &
                          0, 0, 0, &
                          fp%vx, fp%vy)/((ifde - ifds + 1)*(jfde - jfds + 1))
            mw = fun_real(RNRM_MAX, &
                          ifms, ifme, 1, 1, jfms, jfme, &
                          ifds, ifde, 1, 1, jfds, jfde, &
                          ifts, ifte, 1, 1, jfts, jfte, &
                          0, 0, 0, &
                          fp%vx, fp%vy)
!$OMP MASTER
            write (msg, 91) time_start, 'Average surface wind', aw, 'm/s'
            call message(msg, stat_lev)
            write (msg, 91) time_start, 'Maximum surface wind', mw, 'm/s'
            call message(msg, stat_lev)
!$OMP END MASTER
         end if

         call print_2d_stats(ifts, ifte, jfts, jfte, &
                             ifms, ifme, jfms, jfme, &
                             fuel_frac, 'model: fuel_frac start')

         if (.not. freeze_fire) then
            if (.not. replay) then
               !print *, "estou dentro ifun=4 na rotina sfire_model e vou entrar na prop_ls 1"
               !Scall flush (6)
               call prop_ls(ifun, id, 1, &
                            ifds, ifde, jfds, jfde, &
                            ifms, ifme, jfms, jfme, &
                            ifps, ifpe, jfps, jfpe, &
                            ifts, ifte, jfts, jfte, &
                            time_start, dt, fdx, fdy, tbound, &
                            lfn, lfn_out, tign, ros, fp, config_flags &
                            ) !config_flags e ifun int pro ISILDA CM
               !print *, "estou no model,sai do prop_ls"
               !call flush (6)
            else
               do j = jfts, jfte
                  do i = ifts, ifte
                     lfn_out(i, j) = tign(i, j) - (time_start + dt)
                  end do
               end do
            end if
         else
            call message('sfire_model: EXPERIMENTAL: skipping fireline propagation')

         end if

         call print_2d_stats(ifts, ifte, jfts, jfte, ifms, ifme, jfms, jfme, ros, 'model: ros')

      elseif (ifun .eq. 5) then

         if (.not. freeze_fire) then
            if (.not. replay) then
               !print *, "estou dentro ifun=5 na rotina sfire_model e vou entrar na prop_ls 2"
               !call flush (6)
               call prop_ls(ifun, id, 2, &
                            ifds, ifde, jfds, jfde, &
                            ifms, ifme, jfms, jfme, &
                            ifps, ifpe, jfps, jfpe, &
                            ifts, ifte, jfts, jfte, &
                            time_start, dt, fdx, fdy, tbound, &
                            lfn, lfn_out, tign, ros, fp, config_flags &
                            ) !config_flags e ifun int pro ISILDA CM
               !print *, "estou no model sai do prop_ls 2"
               !call flush (6)
            end if
         end if

         call check_lfn_tign('time step end', time_start + dt, ifts, ifte, jfts, jfte, ifms, ifme, jfms, jfme, lfn_out, tign)

         if (.not. freeze_fire) then
         do j = jfts, jfte
            do i = ifts, ifte
               lfn(i, j) = lfn_out(i, j)

            end do
         end do

         call check_lfn_tign('before ignition', time_start + dt, ifts, ifte, jfts, jfte, ifms, ifme, jfms, jfme, lfn, tign)

         if (fire_tign_in_time > tiny(dt)) then
            if (ignition%num_lines > 0) then
               call crash('ignition from lines and from tign_in are not compatible')
            end if
            !print *, "estou dentro ifun=5 na rotina sfire_model e vou entrar na ignite_from_tign_is"
            !call flush (6)
            call ignite_from_tign_in( &
               ifds, ifde, jfds, jfde, &
               ifms, ifme, jfms, jfme, &
               ifts, ifte, jfts, jfte, &
               time_start, time_start + dt, &
               tign_in, &
               lfn, tign, ignited)
            call check_lfn_tign('after ignite_from_tign_in', time_start + dt, ifts, ifte, jfts, jfte, ifms &
                                , ifme, jfms, jfme, lfn, tign)
         end if

         do ig = 1, ignition%num_lines

            !print *, "estou dentro ifun=5 na rotina sfire_model e vou entrar na ignite_fire"
            !call flush (6)
            call ignite_fire( &
               ifds, ifde, jfds, jfde, &
               ifms, ifme, jfms, jfme, &
               ifts, ifte, jfts, jfte, &
               ignition%line(ig), &
               time_start, time_start + dt, &
               coord_xf, coord_yf, ignition%unit_fxlong, ignition%unit_fxlat, &
               lfn, tign, ignited)

            call write_array_m(ifts, ifte, jfts, jfte, ifms, ifme, jfms, jfme, lfn, 'lfn_ig', id)
            call write_array_m(ifts, ifte, jfts, jfte, ifms, ifme, jfms, jfme, coord_xf, 'coord_xf_ig', id)
            call write_array_m(ifts, ifte, jfts, jfte, ifms, ifme, jfms, jfme, coord_yf, 'coord_yf_ig', id)

         end do
         else
         call message('sfire_model: EXPERIMENTAL: skipping ignition')
         end if

         call print_2d_stats(ifts, ifte, jfts, jfte, ifms, ifme, jfms, jfme, &
                             lfn, 'sfire_model: lfn out')

         call check_lfn_tign('after ignition', time_start + dt, ifts, ifte, jfts, jfte, ifms, ifme, jfms, jfme, lfn, tign)

      elseif (ifun .eq. 6) then

         !print *, "estou dentro ifun=6 na rotina sfire_model e vou entrar na fire_intensity"
         !call flush (6)
         call fire_intensity(fp, &
                             ifms, ifme, jfms, jfme, &
                             ifts, ifte, jfts, jfte, &
                             ifms, ifme, jfms, jfme, &
                             ros, nfuel_cat, &
                             flineint, flineint2)

!!!!!! INTRODUZIDO POR ISILDA CUNHA MENEZES
         if (crown) then
            !print *, "estou dentro ifun=6 no module_fire_model e vou entrar na fire_total_intensity"
            !call flush (6)
            call fire_total_intensity(fp, &
                                      ifms, ifme, jfms, jfme, &
                                      ifts, ifte, jfts, jfte, &
                                      ifms, ifme, jfms, jfme, &
                                      fp%fgip, ros, &
                                      flineint_total)
         end if
!!!!!!!!

         if (fireline_mask < 0.) then
            do j = jfts, jfte
               do i = ifts, ifte

                  if ((lfn(i - 1, j - 1) > 0. .and. lfn(i - 1, j) > 0. .and. lfn(i, j - 1) > 0. .and. lfn(i, j) > 0. .and. &
                       lfn(i + 1, j + 1) > 0. .and. lfn(i + 1, j) > 0. .and. lfn(i, j + 1) > 0.) .or. &
                      (lfn(i - 1, j - 1) < 0. .and. lfn(i - 1, j) < 0. .and. lfn(i, j - 1) < 0. .and. lfn(i, j) < 0. .and. &
                       lfn(i + 1, j + 1) < 0. .and. lfn(i + 1, j) < 0. .and. lfn(i, j + 1) < 0.) .or. &
                      i .eq. ifds .or. i .eq. ifde .or. j .eq. jfds .or. j .eq. jfde) then
                     ros(i, j) = fireline_mask
                     flineint(i, j) = fireline_mask
                     flineint2(i, j) = fireline_mask
                  end if
               end do
            end do
         end if

         !print *, "estou dentro ifun=6 na rotina sfire_model e vou entrar na fire_risk"
         !call flush (6)
         call fire_risk(fp, &
                        ifms, ifme, jfms, jfme, &
                        ifts, ifte, jfts, jfte, &
                        nfuel_cat, &
                        f_ros0, f_rosx, f_rosy, f_ros, &
                        f_int, f_lineint, f_lineint2)
                        
         select case (fire_hfx_given)

         case (0)

!$OMP CRITICAL(SFIRE_MODEL_CRIT)
            write (msg, *) 'time_start=', time_start, ' dt=', dt, ' before fuel_left'
!$OMP END CRITICAL(SFIRE_MODEL_CRIT)
            call message(msg)
            call print_2d_stats(ifts, ifte, jfts, jfte, ifms, ifme, jfms, jfme, lfn, 'model: lfn')
            call print_2d_stats(ifts, ifte, jfts, jfte, ifms, ifme, jfms, jfme, tign, 'model: tign')
            call check_lfn_tign('before fuel_left', time_start + dt, ifts, ifte, jfts, jfte, ifms, ifme, jfms, jfme, lfn, tign)

            call print_2d_stats(ifts, ifte, jfts, jfte, ifms, ifme, jfms, jfme, fuel_time, 'model: fuel_time')
            print *,'LFR-DBG: fire_update_fuel_frac=',fire_update_fuel_frac
            if (fire_update_fuel_frac .eq. 1) then
               !print *, "estou dentro ifun=6 na rotina sfire_model e vou entrar na fuel_left"
               !call flush (6)
               call fuel_left( &
                  ifds, ifde, jfds, jfde, &
                  ifms, ifme, jfms, jfme, &
                  ifts, ifte, jfts, jfte, &
                  ifts, ifte, jfts, jfte, &
                  lfn, tign, fuel_time, time_start + dt, fuel_frac_end, fire_area)

               call print_2d_stats(ifts, ifte, jfts, jfte, ifts, ifte, jfts, jfte, fuel_frac_end, 'model: fuel_frac end')

               do j = jfts, jfte
                  do i = ifts, ifte
                     t = min(fuel_frac(i, j), fuel_frac_end(i, j))
                     fuel_frac_burnt(i, j) = fuel_frac(i, j) - t
                     fuel_frac(i, j) = t
                  end do
               end do
            elseif (fire_update_fuel_frac .eq. 2) then
               do j = jfts, jfte
                  do i = ifts, ifte
                     if (lfn(i, j) < 0.) then
                        fuel_frac_burnt(i, j) = dt/fuel_time(i, j)
                     else
                        fuel_frac_burnt(i, j) = 0.
                     end if
                  end do
               end do
            else
               call crash('fire_update_fuel_frac value not supported')
            end if

            call print_2d_stats(ifts, ifte, jfts, jfte, ifms, ifme, jfms, jfme, fuel_frac_burnt, 'model: fuel_frac burned')
            !print *, "estou dentro ifun=6 na rotina sfire_model e vou entrar na heat_fluxes"
            !call flush (6)
            call heat_fluxes(dt, fp, &
                             ifms, ifme, jfms, jfme, &
                             ifts, ifte, jfts, jfte, &
                             ifms, ifme, jfms, jfme, &
                             fp%fgip, &
                             fuel_frac_burnt, &
                             fgrnhfx, fgrnqfx)
               !LFR-DBG beg
               do j = jfts, jfte
                  do i = ifts, ifte            
                     WRITE(91,*) i,j,fp%fgip(i,j),fuel_frac_burnt(i,j),fgrnhfx(i,j),fgrnqfx(i,j)
                  enddo
               enddo
               !LFR-DBG end
         case (1, 2)
!$OMP CRITICAL(SFIRE_MODEL_CRIT)
            write (msg, *) "model: expecting fire_hfx to be set in WRF, from wrfinput or wrfrst files"
            call message(msg)
!$OMP END CRITICAL(SFIRE_MODEL_CRIT)

            do j = jfts, jfte
               do i = ifts, ifte
                  fgrnhfx(i, j) = (1.-fire_hfx_latent_part)*fire_hfx(i, j)
                  fgrnqfx(i, j) = fire_hfx_latent_part*fire_hfx(i, j)
               end do
            end do

         case (3)

            call message('artificial heat flux from parameters given in namelist.input')
            !print *, "estou dentro ifun=6 na rotina sfire_model e vou entrar na param_hfx"
            !call flush (6)
            call param_hfx(time_start, &
                           ifms, ifme, jfms, jfme, &
                           ifts, ifte, jfts, jfte, &
                           coord_xf, coord_yf, &
                           hfx, &
                           fire_area, fgrnhfx, fgrnqfx)

         case default
            call crash('bad fire_hfx_given')
         end select

         if (fire_print_msg .ge. stat_lev) then
            tfa = fun_real(REAL_SUM, &
                           ifms, ifme, 1, 1, jfms, jfme, &
                           ifds, ifde, 1, 1, jfds, jfde, &
                           ifts, ifte, 1, 1, jfts, jfte, &
                           0, 0, 0, &
                           fire_area, fire_area)*fdx*fdy
            thf = fun_real(REAL_SUM, &
                           ifms, ifme, 1, 1, jfms, jfme, &
                           ifds, ifde, 1, 1, jfds, jfde, &
                           ifts, ifte, 1, 1, jfts, jfte, &
                           0, 0, 0, &
                           fgrnhfx, fgrnhfx)*fdx*fdy
            mhf = fun_real(REAL_MAX, &
                           ifms, ifme, 1, 1, jfms, jfme, &
                           ifds, ifde, 1, 1, jfds, jfde, &
                           ifts, ifte, 1, 1, jfts, jfte, &
                           0, 0, 0, &
                           fgrnhfx, fgrnhfx)
            tqf = fun_real(REAL_SUM, &
                           ifms, ifme, 1, 1, jfms, jfme, &
                           ifds, ifde, 1, 1, jfds, jfde, &
                           ifts, ifte, 1, 1, jfts, jfte, &
                           0, 0, 0, &
                           fgrnqfx, fgrnqfx)*fdx*fdy
            mqf = fun_real(REAL_MAX, &
                           ifms, ifme, 1, 1, jfms, jfme, &
                           ifds, ifde, 1, 1, jfds, jfde, &
                           ifts, ifte, 1, 1, jfts, jfte, &
                           0, 0, 0, &
                           fgrnqfx, fgrnqfx)
!$OMP MASTER
            write (msg, 91) time_start, 'Fire area           ', tfa, 'm^2'
            call message(msg, stat_lev)
            write (msg, 91) time_start, 'Heat output         ', thf, 'W'
            call message(msg, stat_lev)
            write (msg, 91) time_start, 'Max heat flux       ', mhf, 'W/m^2'
            call message(msg, stat_lev)
            write (msg, 91) time_start, 'Latent heat output  ', tqf, 'W'
            call message(msg, stat_lev)
            write (msg, 91) time_start, 'Max latent heat flux', mqf, 'W/m^2'
            call message(msg, stat_lev)
!$OMP END MASTER
91          format('Time ', f11.3, ' s ', a, e12.3, 1x, a)
         end if

         call print_2d_stats(ifts, ifte, jfts, jfte, &
                             ifms, ifme, jfms, jfme, &
                             fgrnhfx, 'model: heat flux(J/m^2/s)')

      else
!$OMP CRITICAL(SFIRE_MODEL_CRIT)
         write (msg, *) 'sfire_model: bad ifun=', ifun
!$OMP END CRITICAL(SFIRE_MODEL_CRIT)
         call crash(msg)
      end if

   end subroutine sfire_model

   subroutine param_hfx(time_now, &
                        ifms, ifme, jfms, jfme, &
                        ifts, ifte, jfts, jfte, &
                        coord_xf, coord_yf, &
                        hfx, &
                        fire_area, fgrnhfx, fgrnqfx)

      real, intent(in)::time_now
      integer, intent(in):: &
         ifms, ifme, jfms, jfme, &
         ifts, ifte, jfts, jfte
      type(lines_type), intent(in)::hfx
      real, dimension(ifms:ifme, jfms:jfme), intent(in)::coord_xf, coord_yf
      real, dimension(ifms:ifme, jfms:jfme), intent(out)::fire_area, fgrnhfx, fgrnqfx
      character(len=128):: msg

      integer::i, j, k, nfa, ncells
      real:: d2, ax, ay, hfrac, fa, thfx, t, r, radius
      real, parameter:: sigma_mult = 3.
      real:: maxspeed = 100

      do j = jfts, jfte
         do i = ifts, ifte
            fire_area(i, j) = 0
            fgrnhfx(i, j) = 0.
            fgrnqfx(i, j) = 0.
         end do
      end do

      do k = 1, hfx%num_lines
         if (hfx%line(k)%radius > 0.) then

            t = max(hfx%line(k)%start_time - time_now, time_now - hfx%line(k)%end_time)
            if (t > 0.) then
               r = t/hfx%line(k)%trans_time
               hfrac = exp(-(sigma_mult*r)**2/2.)
            else
               hfrac = 1.0
            end if

            ax = safe_prop(time_now, &
                           hfx%line(k)%start_time, &
                           hfx%line(k)%end_time, &
                           hfx%line(k)%start_x, &
                           hfx%line(k)%end_x, &
                           hfx%unit_fxlong*maxspeed)
            ay = safe_prop(time_now, &
                           hfx%line(k)%start_time, &
                           hfx%line(k)%end_time, &
                           hfx%line(k)%start_y, &
                           hfx%line(k)%end_y, &
                           hfx%unit_fxlat*maxspeed)

            radius = hfx%line(k)%radius

!$OMP CRITICAL(SFIRE_MODEL_CRIT)
            write (msg, *) 'hfx line ', i, ' at ', time_now, 's ', hfrac, ' of max ', hfx%line(k)%hfx_value
            call message(msg)
            write (msg, *) 'center ', ax, ay, ' radius ', radius
            call message(msg)
!$OMP END CRITICAL(SFIRE_MODEL_CRIT)

            nfa = 0
            ncells = 0
            do j = jfts, jfte
               do i = ifts, ifte

                  d2 = (hfx%unit_fxlong*(ax - coord_xf(i, j)))**2 + (hfx%unit_fxlat*(ay - coord_yf(i, j)))**2
                  if (d2 < radius**2) then
                     fa = 1.
                  else
                     fa = 0.
                  end if

                  thfx = hfx%line(k)%hfx_value*hfrac*fa
                  fgrnhfx(i, j) = fgrnhfx(i, j) + (1.-fire_hfx_latent_part)*thfx
                  fgrnqfx(i, j) = fgrnqfx(i, j) + fire_hfx_latent_part*thfx

                  fire_area(i, j) = max(fire_area(i, j), fa)

                  nfa = nfa + fa; 
                  ncells = ncells + 1
               end do
            end do

!$OMP CRITICAL(SFIRE_MODEL_CRIT)
            write (msg, *) 'Number of cells in fire area in this tile ', nfa, ' ', (100.*nfa)/ncells, ' %'
            call message(msg)
!$OMP END CRITICAL(SFIRE_MODEL_CRIT)

         end if
      end do
   end subroutine param_hfx

   real function safe_prop(t, t1, t2, x1, x2, ms)

      real, intent(in)::t, t1, t2, x1, x2, ms
      real:: p, x
      if (t2 < t1) call crash('safe_prop: must have t2>t1')
      if (.not. ms > 0.) call crash('safe_prop: must have ms>0')
      if (t1 .eq. t2) then
         if (x1 .eq. x2) then
            x = x1
         else
            call crash('safe_prop: infinite speed')
         end if
      else
         p = (t - t1)/(t2 - t1)
         x = x1*(1.-p) + x2*p
      end if
      safe_prop = x
   end function safe_prop

end module module_fr_sfire_model
