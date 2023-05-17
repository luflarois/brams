module module_fr_sfire_phys

   use module_model_constants, only: cp, xlv
   use module_fr_sfire_util

!!!! INTRODUZIDO POR ISILDA CUNHA MENEZES
!use module_fr_sfire_crown
   use ModNamelistsfireFile, only: grid_config_rec_type
!!!!!

   implicit none
   private

   type fire_params
      real, pointer, dimension(:, :):: vx, vy
      real, pointer, dimension(:, :):: zsf
      real, pointer, dimension(:, :):: dzdxf, dzdyf
      real, pointer, dimension(:, :):: bbb, phisc, phiwc, r_0
      real, pointer, dimension(:, :):: fgip
      real, pointer, dimension(:, :):: ischap
      real, pointer, dimension(:, :):: fuel_time
      real, pointer, dimension(:, :):: fmc_g
      real, pointer, dimension(:, :):: nfuel_cat

!!!!!!!INTRODUZIDO POR ISILDA CUNHA MENEZES
      real, pointer, dimension(:, :):: CFB
      real, pointer, dimension(:, :):: HPA
      real, pointer, dimension(:, :):: bbb_FM10, phisc_FM10, phiwc_FM10, r_0_FM10
!!!!!!!!

   end type fire_params

!!! ALTERADO POR ISILDA CUNHA MENEZES
   public:: init_fuel_cats, heat_fluxes, set_nfuel_cat, set_fire_params, &
            write_fuels_m, fire_risk, fire_intensity, fuel_moisture, advance_moisture, fuel_name, &
            fire_rate_of_spread, set_fire_crown_params, fire_ros, fire_total_intensity
!!!!!
!PUBLIC:: init_fuel_cats,fire_ros,heat_fluxes,set_nfuel_cat,set_fire_params, &
!write_fuels_m,fire_risk,fire_intensity,fuel_moisture,advance_moisture,fuel_name,&
!fire_rate_of_spread

   public:: fire_params

   public:: fire_wind_height, fcz0, fcwh, have_fuel_cats, nfuelcats, no_fuel_cat, no_fuel_cat2, windrf, moisture_classes &
            , itree, fueldepthm
   public:: mfuelcats

   integer, parameter :: mfuelcats = 100
   integer, parameter::max_moisture_classes = 5

   integer, parameter::zm = max_moisture_classes - 3
   integer:: moisture_classes = 3
   real, dimension(max_moisture_classes):: drying_lag, wetting_lag, saturation_moisture, saturation_rain, &
                                           rain_threshold, rec_drying_lag_sec, rec_wetting_lag_sec
   integer, dimension(max_moisture_classes):: drying_model, wetting_model, fmc_gc_initialization

   integer::itmp
   character(len=80), dimension(max_moisture_classes), save :: moisture_class_name
   real, dimension(mfuelcats):: &
      fmc_gw01 = (/(1.0, itmp=1, mfuelcats)/), &
      fmc_gw02 = (/(0.0, itmp=1, mfuelcats)/), &
      fmc_gw03 = (/(0.0, itmp=1, mfuelcats)/), &
      fmc_gw04 = (/(0.0, itmp=1, mfuelcats)/), &
      fmc_gw05 = (/(0.0, itmp=1, mfuelcats)/)

   data moisture_class_name/'1-hour fuel', '10-hour fuel', '100-hour fuel', zm*'NOT USED'/
   data drying_lag/1., 10., 100., zm*0./
   data wetting_lag/14, 140., 1400., zm*0./
   data saturation_moisture/2.5, 2.5, 2.5, zm*0./
   data saturation_rain/8.0, 8.0, 8.0, zm*0./
   data rain_threshold/0.05, 0.05, 0.05, zm*0/
   data drying_model/1, 1, 1, zm*1/
   data wetting_model/1, 1, 1, zm*1/
   data fmc_gc_initialization/2, 2, 2, zm*2/
   real, dimension(7)::eq_p
   data eq_p/1.035e-09, &
      -2.62e-07, &
      2.507e-05, &
      -0.001107, &
      0.02245, &
      -0.05901, &
      3.043/

   real, save:: cmbcnst, hfgl, fuelmc_g, fuelmc_c, fire_wind_height

   real, save:: fuelheat

   data cmbcnst/17.433e+06/
   data hfgl/17.e4/
   data fuelmc_g/0.08/
   data fuelmc_c/1.00/
   data fire_wind_height/6.096/

   integer, parameter :: nf = 14
   integer, save      :: nfuelcats = 13
   integer, parameter :: zf = mfuelcats - nf
   integer, save      :: no_fuel_cat = 14
   integer, save      :: no_fuel_cat2 = 99999999
   integer, save      :: ibeh = 1
   character(len=80), dimension(mfuelcats), save :: fuel_name
   integer, dimension(mfuelcats), save :: ichap, itree
   real, dimension(mfuelcats), save :: windrf, weight, fgi, fci, fci_d, fct, fcbr, &
                                       fueldepthm, fueldens, fuelmce, &
                                       fcwh, fcz0, ffw, &
                                       savr, st, se, adjr0, adjrw, adjrs, &
                                       fmc_gl_stdev, fmc_gl_ndwi_0, fmc_gl_ndwi_rate, fmc_gl_ndwi_stdev, &
                                       CBH, FMC_crown !!!! colocados pelo autor Isilda da Cunha Menezes

   real, dimension(mfuelcats, max_moisture_classes), save :: fmc_gw

   data fuel_name/ &
      'FUEL MODEL 1: Short grass (1 ft)', &
      'FUEL MODEL 2: Timber (grass and understory)', &
      'FUEL MODEL 3: Tall grass (2.5 ft)', &
      'FUEL MODEL 4: Chaparral (6 ft)', &
      'FUEL MODEL 5: Brush (2 ft) ', &
      'FUEL MODEL 6: Dormant brush, hardwood slash', &
      'FUEL MODEL 7: Southern rough', &
      'FUEL MODEL 8: Closed timber litter', &
      'FUEL MODEL 9: Hardwood litter', &
      'FUEL MODEL 10: Timber (litter + understory)', &
      'FUEL MODEL 11: Light logging slash', &
      'FUEL MODEL 12: Medium logging slash', &
      'FUEL MODEL 13: Heavy logging slash', &
      'FUEL MODEL 14: no fuel', &
      zf*' '/
   data windrf/0.36, 0.36, 0.44, 0.55, 0.42, 0.44, 0.44, &
      0.36, 0.36, 0.36, 0.36, 0.43, 0.46, 1e-7, zf*0/
   data fgi/0.166, 0.897, 0.675, 2.468, 0.785, 1.345, 1.092, &
      1.121, 0.780, 2.694, 2.582, 7.749, 13.024, 1.e-7, zf*0./
   data fueldepthm/0.305, 0.305, 0.762, 1.829, 0.61, 0.762, 0.762, &
      0.0610, 0.0610, 0.305, 0.305, 0.701, 0.914, 0.305, zf*0./
   data savr/3500., 2784., 1500., 1739., 1683., 1564., 1562., &
      1889., 2484., 1764., 1182., 1145., 1159., 3500., zf*0./
   data fuelmce/0.12, 0.15, 0.25, 0.20, 0.20, 0.25, 0.40, &
      0.30, 0.25, 0.25, 0.15, 0.20, 0.25, 0.12, zf*0./
   data fueldens/nf*32., zf*0./
   data st/nf*0.0555, zf*0./
   data se/nf*0.010, zf*0./

   data weight/7., 7., 7., 180., 100., 100., 100., &
      900., 900., 900., 900., 900., 900., 7., zf*0./

   data fci_d/0., 0., 0., 1.123, 0., 0., 0., &
      1.121, 1.121, 1.121, 1.121, 1.121, 1.121, 0., zf*0./
   data fct/60., 60., 60., 60., 60., 60., 60., &
      60., 120., 180., 180., 180., 180., 60., zf*0./
   data ichap/0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, zf*0/
   data itree/0, 5, 0, 0, 0, 5, 5, 5, 5, 5, 0, 0, 0, 0, zf*0/
   data fcwh/6.096, 6.096, 6.096, 6.096, 6.096, 6.096, 6.096, &
      6.096, 6.096, 6.096, 6.096, 6.096, 6.096, 6.096, zf*0./

   data fcz0/0.0396, 0.0396, 0.1000, 0.2378, 0.0793, 0.0991, 0.0991, &
      0.0079, 0.0079, 0.0396, 0.0396, 0.0911, 0.1188, 0.0396, zf*0./

   data ffw/nf*0.9, zf*0/
   data fmc_gl_ndwi_0/nf*0.1, zf*0./
   data fmc_gl_ndwi_rate/nf*0.6, zf*0./
   data fmc_gl_ndwi_stdev/nf*0.2, zf*0./
   data fmc_gl_stdev/nf*0.2, zf*0./
   data adjr0/mfuelcats*1./
   data adjrw/mfuelcats*1./
   data adjrs/mfuelcats*1./

!!!!! AUTOR ISILDA DA CUNHA MENEZES
   data CBH/0.305, 0.305, 0.762, 1.829, 0.61, 0.762, 0.762, &
      0.0610, 0.0610, 0.305, 0.305, 0.701, 0.914, 0.305, zf*0./
   data FMC_crown/0.12, 0.15, 0.25, 0.20, 0.20, 0.25, 0.40, &
      0.30, 0.25, 0.25, 0.15, 0.20, 0.25, 0.12, zf*0./

!!!!!

   logical, save :: have_fuel_cats = .false.

contains

   subroutine fuel_moisture( &
      id, &
      nfmc, &
      ids, ide, jds, jde, &
      ims, ime, jms, jme, &
      ips, ipe, jps, jpe, &
      its, ite, jts, jte, &
      ifds, ifde, jfds, jfde, &
      ifms, ifme, jfms, jfme, &
      ifts, ifte, jfts, jfte, &
      ir, jr, &
      nfuel_cat, &
      fndwi, &
      fmc_gc, &
      fmc_g &
      )

      implicit none

      integer, intent(in):: &
         id, nfmc, &
         ids, ide, jds, jde, &
         ims, ime, jms, jme, &
         ips, ipe, jps, jpe, &
         its, ite, jts, jte, &
         ifds, ifde, jfds, jfde, &
         ifms, ifme, jfms, jfme, &
         ifts, ifte, jfts, jfte, &
         ir, jr

      real, intent(in), dimension(ifms:ifme, jfms:jfme):: nfuel_cat, &
                                                          fndwi
      real, intent(inout), dimension(ims:ime, nfmc, jms:jme):: fmc_gc
      real, intent(out), dimension(ifms:ifme, jfms:jfme):: fmc_g

      real, dimension(its - 1:ite + 1, jts - 1:jte + 1):: fmc_k
      real, dimension(ifts:ifte, jfts:jfte):: fmc_f, &
                                              nwdi_f
      integer::i, j, k, n
      integer::ibs, ibe, jbs, jbe
      real::f1, w1, w2, f2, fa, fc

      character(len=128)::msg

      !print *, "estou no modulo_fr_drive_phys e estou na rotina fuel_moisture"; call flush (6)
      !call flush (6)
      call check_mesh_2dim(ifts, ifte, jfts, jfte, ifds, ifde, jfds, jfde)
      call check_mesh_2dim(ifts, ifte, jfts, jfte, ifms, ifme, jfms, jfme)

      do j = jfts, jfte
         do i = ifts, ifte
            fmc_g(i, j) = 0.
         end do
      end do

      ibs = max(ids, its - 1)
      ibe = min(ide, ite + 1)
      jbs = max(jds, jts - 1)
      jbe = min(jde, jte + 1)

      call check_mesh_2dim(ibs, ibe, jbs, jbe, ims, ime, jms, jme)

      do k = 1, moisture_classes

         do j = jbs, jbe
            do i = ibs, ibe
               fmc_k(i, j) = fmc_gc(i, k, j)
            end do
         end do

         !LFR call print_2d_stats(ibs, ibe, jbs, jbe, its - 1, ite + 1, jts - 1, jte + 1, fmc_k, 'fuel_moisture: fmc_k')

         !print *, "estou no modulo_fr_drive_phys e vou chamar interpolate z2fire"; call flush (6)
         !call flush (6)
         call interpolate_z2fire(id, 0, &
                                 ids, ide, jds, jde, &
                                 its - 1, ite + 1, jts - 1, jte + 1, &
                                 ips, ipe, jps, jpe, &
                                 its, ite, jts, jte, &
                                 ifds, ifde, jfds, jfde, &
                                 ifts, ifte, jfts, jfte, &
                                 ifts, ifte, jfts, jfte, &
                                 ir, jr, &
                                 fmc_k, &
                                 fmc_f)
         print *, 'LFR-DBG: Sai de interpolate_z2fire'; call flush (6)
         !LFR call print_2d_stats(ifts, ifte, jfts, jfte, ifts, ifte, jfts, jfte, fmc_f, 'fuel_moisture: fmc_f')

         if (k .eq. kfmc_ndwi) then
            call print_2d_stats(ifts, ifte, jfts, jfte, ifms, ifme, jfms, jfme, fndwi, 'fuel_moisture: fndwi')
            write (msg, '(a,i4)') 'Assimilating NDWI in fuel moisture class ', k
            call message(msg)
         end if

         !print *, "estou no modulo_fr_drive_phys e vou fazer calculo dos w e f e fmc_g"
         !call flush (6)
         do j = jfts, jfte
            do i = ifts, ifte
               n = nfuel_cat(i, j)
               if (n > 0) then
                  if (k .ne. kfmc_ndwi) then
                     fmc_g(i, j) = fmc_g(i, j) + fmc_gw(n, k)*fmc_f(i, j)
                  else
                     f1 = fmc_f(i, j)
                     w1 = fmc_gl_stdev(n)
                     w1 = 1./(w1*w1)
                     w2 = fmc_gl_ndwi_stdev(n)
                     w2 = 1./(w2*w2)
                     f2 = fmc_gl_ndwi_0(n) + fmc_gl_ndwi_rate(n)*fndwi(i, j)
                     fa = (w1*f1 + w2*f2)/(w1 + w2)
                     fc = fmc_gw(n, k)*fa
                     fmc_g(i, j) = fmc_g(i, j) + fc

                  end if
               end if
            end do
         end do

      end do

      !LFR call print_2d_stats(ifts, ifte, jfts, jfte, ifms, ifme, jfms, jfme, fmc_g, 'fuel_moisture: fmc_g')

   end subroutine fuel_moisture

   subroutine advance_moisture( &
      initialize, &
      ims, ime, jms, jme, &
      its, ite, jts, jte, &
      nfmc, &
      moisture_dt, &
      fmep_decay_tlag, &
      rainc, rainnc, &
      t2, q2, psfc, &
      rain_old, &
      t2_old, q2_old, psfc_old, &
      rh_fire, &
      fmc_gc, &
      fmep, &
      fmc_equi, &
      fmc_lag &
      )

      implicit none

      logical, intent(in):: initialize
      integer, intent(in):: &
         ims, ime, jms, jme, &
         its, ite, jts, jte, &
         nfmc
      real, intent(in):: moisture_dt, fmep_decay_tlag
      real, intent(in), dimension(ims:ime, jms:jme):: t2, q2, psfc, rainc, rainnc
      real, intent(inout), dimension(ims:ime, jms:jme):: t2_old, q2_old, psfc_old, rain_old
      real, intent(inout), dimension(ims:ime, nfmc, jms:jme):: fmc_gc
      real, intent(inout), dimension(ims:ime, 2, jms:jme):: fmep
      real, intent(out), dimension(ims:ime, nfmc, jms:jme):: fmc_equi, fmc_lag
      real, intent(out), dimension(ims:ime, jms:jme)::rh_fire

      integer:: i, j, k
      real::rain_int, T, P, Q, QRS, ES, RH, tend, EMC_d, EMC_w, EMC, R, rain_diff, fmc, rlag, equi, &
             d, w, rhmax, rhmin, change, rainmax, rainmin, fmc_old, H, deltaS, deltaE
      real, parameter::tol = 1e-2
      character(len=256)::msg
      logical::bad_wrf
      integer::msglevel = 2
      logical, parameter::check_rh = .false.
      integer::check_data = 2
      real::epsilon, Pws, Pw, t2_min, q2_min, psfc_min
      real::t2_floor = 200.
      real::q2_floor = 1e-8
      real::psfc_floor = 1000.
      integer :: ii,jj !LFR

! !LFR
!                         do ii=its,ite
!                            do jj=jts,jte
!                               write(60,*) ii,jj, fmep(ii, 1, jj)
!                            end do
!                         end do
! !LFR
      !print *, "estou no modulo_fr_drive_phys e estou na rotina advance_moisture"
      !call flush (6)
      if (msglevel > 1) then
!$OMP CRITICAL(SFIRE_PHYS_CRIT)
         write (msg, '(a,f10.2,a,i4,a,i4)') 'advance_moisture dt=', moisture_dt, 's using ', moisture_classes &
            , ' classes from possible ', nfmc
!$OMP END CRITICAL(SFIRE_PHYS_CRIT)
         call message(msg, level=2)
      end if
!print*, "estou no modulo_fr_drive_phys e estou na rotina advance_moisture A1"
!call flush(6)
      if (moisture_classes > nfmc .or. moisture_classes > max_moisture_classes) then
!$OMP CRITICAL(SFIRE_PHYS_CRIT)
         write (msg, *) 'advance_moisture: moisture_classes=', moisture_classes, &
            ' > nfmc=', nfmc, ' or >  max_moisture_classes=', max_moisture_classes
!$OMP END CRITICAL(SFIRE_PHYS_CRIT)
         call crash(msg)
      end if
!print*, "estou no modulo_fr_drive_phys e estou na rotina advance_moisture A2"
!call flush(6)
      call print_2d_stats(its, ite, jts, jte, ims, ime, jms, jme, t2, 'T2')
      call print_2d_stats(its, ite, jts, jte, ims, ime, jms, jme, q2, 'Q2')
      call print_2d_stats(its, ite, jts, jte, ims, ime, jms, jme, psfc, 'PSFC')
!print*, "estou no modulo_fr_drive_phys e estou na rotina advance_moisture A3"
!call flush(6)
      if (initialize) then
         call message('advance_moisture: initializing, copying surface variables to old')
         call copy2old
      else
         !   print*,"estou no modulo_fr_drive_phys vou para print_3d_stats_by_slice A4"
         !   call flush(6)
         call print_3d_stats_by_slice(its, ite, 1, moisture_classes, jts, jte, ims, ime, 1, nfmc, jms, jme, fmc_gc &
                                      , 'before advance fmc_gc')
         !   print*,"estou no modulo_fr_drive_phys  sai do print_3d_stats_by_slice A5"
         !   call flush(6)
      end if

!print*, "estou no modulo_fr_drive_phys e estou na rotina advance_moisture A6"
!call flush(6)
      if (check_data .ge. 2 .or. msglevel .ge. 2) then
         t2_min = huge(t2_min)
         q2_min = huge(q2_min)
         psfc_min = huge(psfc_min)
         do j = jts, jte
            do i = its, ite
               t2_min = min(t2(i, j), t2_min)
               q2_min = min(q2(i, j), q2_min)
               psfc_min = min(psfc(i, j), psfc_min)
            end do
         end do

!print*, "estou no modulo_fr_drive_phys e estou na rotina advance_moisture A7"
!call flush(6)
         bad_wrf = (t2_min < t2_floor .or. psfc_min < psfc_floor .or. q2_min < q2_floor)
         if (bad_wrf .or. msglevel .ge. 2) then
!$OMP CRITICAL(SFIRE_PHYS_CRIT)
91          format(a, 3(2x, a, e11.3))
            write (msg, 91) 'minimal  ', 't2', t2_min, 'q2', q2_min, 'psfc', psfc_min
            call message(msg, level=0)
            write (msg, 91) 'floor    ', 't2', t2_floor, 'q2', q2_floor, 'psfc', psfc_floor
            call message(msg, level=0)
!$OMP END CRITICAL(SFIRE_PHYS_CRIT)
         end if

!print*, "estou no modulo_fr_drive_phys e estou na rotina advance_moisture A8"
!call flush(6)
         if (bad_wrf) then
            if (check_data .ge. 3) then
               call crash('advance_moisture: invalid data passed from WRF')
            else
               call message('WARNING: advance_moisture: nonphysical input values replaced by floor', level=0)
            end if
         end if
      end if

!print*, "estou no modulo_fr_drive_phys e estou na rotina advance_moisture A9"
!call flush(6)
      rhmax = -huge(rhmax)
      rhmin = huge(rhmin)
      rainmax = -huge(rainmax)
      rainmin = huge(rainmin)
!print*, "estou no modulo_fr_drive_phys e estou na rotina advance_moisture A10"
!call flush(6)
      do j = jts, jte
         do k = 1, moisture_classes
            do i = its, ite

               rain_diff = ((rainc(i, j) + rainnc(i, j)) - rain_old(i, j))
               !  print*,'rotina advance_moisture rain_diff',rain_diff
               !  call flush(6)
               !  print*,'rotina advance_moisture moisture_dt',moisture_dt
               !  call flush(6)
               if (moisture_dt > 0.) then
                  rain_int = 3600.*rain_diff/moisture_dt
               else
                  rain_int = 0.
               end if
               rainmax = max(rainmax, rain_int)
               rainmin = min(rainmin, rain_int)
               R = rain_int - rain_threshold(k)

               T = 0.5*(t2_old(i, j) + t2(i, j))
               P = 0.5*(psfc_old(i, j) + psfc(i, j))
               Q = 0.5*(q2_old(i, j) + q2(i, j))
               ! print*,'rotina advance_moisture T',T
               ! print*,'rotina advance_moisture Q',Q
               ! print*,'rotina advance_moisture P',P
               if (check_data .ge. 1) then
                  T = max(T, t2_floor)
                  P = max(P, psfc_floor)
                  Q = max(Q, q2_floor)
               end if

               epsilon = 0.622
               Pw = q*P/(epsilon + (1 - epsilon)*q); 
               Pws = exp(54.842763 - 6763.22/T - 4.210*log(T) + 0.000367*T + &
                         tanh(0.0415*(T - 218.8))*(53.878 - 1331.22/T - 9.44523*log(T) + 0.014025*T))
               RH = Pw/Pws
               rh_fire(i, j) = RH
               rhmax = max(RH, rhmax)
               rhmin = min(RH, rhmin)
               ! print*,"deltaE"
               deltaE = fmep(i, 1, j)
               !if (deltaE .gt. 0.0) then
                  ! print*,"deltaE=",deltaE,fmep(i,1,j)
                  ! call flush(6)
               !end if
               deltaS = fmep(i, 2, j)
               if (deltaS .gt. 0.0) then
                  ! print*,"deltaS=",deltaS,fmep(i,2,j)
                  ! call flush(6)
               end if
               if (.not. check_rh) then
                  RH = min(RH, 1.0)
                  !    print*,'rotina advance_moisture RH', RH
                  !    call flush(6)
               else
                  if (RH < 0.0 .or. RH > 1.0 .or. RH .ne. RH) then
!$OMP CRITICAL(SFIRE_PHYS_CRIT)
                     write (msg, '(a,2i6,5(a,f10.2))') 'At i,j ', i, j, ' RH=', RH, &
                        ' from T=', T, ' P=', P, ' Q=', Q
                     call message(msg)
                     call crash('Relative humidity must be between 0 and 1, saturated water contents must be >0')
!$OMP END CRITICAL(SFIRE_PHYS_CRIT)
                  end if
               end if
               ! print*,"ADVANCE A1"
               ! call flush(6)
               if (R > 0.) then
                  select case (wetting_model(k))
                  case (1)
                     EMC_w = saturation_moisture(k) + deltaS
                     EMC_d = saturation_moisture(k) + deltaS
                     rlag = rec_wetting_lag_sec(k)*(1.-exp(-R/saturation_rain(k)))
                  end select
               else
                  select case (drying_model(k))
                  case (1)
                     H = RH*100.
                     d = 0.942*H**0.679 + 0.4994e-4*exp(0.1*H) + 0.18*(21.1 + 273.15 - T)*(1 - exp(-0.115*H))
                     w = 0.618*H**0.753 + 0.4540e-4*exp(0.1*H) + 0.18*(21.1 + 273.15 - T)*(1 - exp(-0.115*H))
                     if (d .ne. d .or. w .ne. w) call crash('equilibrium moisture calculation failed, result is NaN')
                     d = d*0.01
                     w = w*0.01
                     ! if(isnan(deltaE)) Then
                     !    print *, 'LFR-DBG: i,j,k,deltaE=',i,j,k,deltaE; call flush(6)
                     !    do ii=its,ite
                     !       do jj=jts,jte
                     !          write(55,*) ii,jj, fmep(ii, 1, jj)
                     !       end do
                     !    end do
                     ! endif
                        
                     EMC_d = max(max(d, w) + deltaE, 0.0)
                     EMC_w = max(min(d, w) + deltaE, 0.0)
                     rlag = rec_drying_lag_sec(k)
                  end select
               end if

               ! print*,"ADVANCE A2"
               ! call flush(6)
               if (rlag > 0.0) then

                  if (.not. initialize .or. fmc_gc_initialization(k) .eq. 0) then
                     fmc_old = fmc_gc(i, k, j)
                  elseif (fmc_gc_initialization(k) .eq. 1) then
                     fmc_old = fuelmc_g
                  elseif (fmc_gc_initialization(k) .eq. 2) then
                     fmc_old = 0.5*(EMC_d + EMC_w)
                  else
                     call crash('bad value of fmc_gc_initialization(k), must be between 0 and 2')
                  end if
                  equi = max(min(fmc_old, EMC_d), EMC_w)

                  change = moisture_dt*rlag

                  if (change < tol) then
                     if (fire_print_msg .ge. 3) call message('midpoint method')
                     fmc = fmc_old + (equi - fmc_old)*change*(1.0 - 0.5*change)
                  else
                     if (fire_print_msg .ge. 3) call message('exponential method')
                     fmc = fmc_old + (equi - fmc_old)*(1 - exp(-change))
                  end if
                  fmc_gc(i, k, j) = fmc
                  !    print*,'ADVANCE fmc_gc',  fmc_gc(i,k,j)
                  !    call flush(6)

                  fmc_equi(i, k, j) = equi
                  !   print*,'ADVANCE fmc_equi',  fmc_equi(i,k,j)
                  !   call flush(6)
                  fmc_lag(i, k, j) = 1.0/(3600.0*rlag)
                  !   print*,'ADVANCE fmc_lag',  fmc_lag(i,k,j)
                  !   call flush(6)
                  if (fire_print_msg .ge. 3) then
!$OMP CRITICAL(SFIRE_PHYS_CRIT)
                     write (msg, *) 'i=', i, ' j=', j, 'EMC_w=', EMC_w, ' EMC_d=', EMC_d
                     call message(msg)
                     write (msg, *) 'fmc_old=', fmc, ' equi=', equi, ' change=', change, ' fmc=', fmc
                     call message(msg)
!$OMP END CRITICAL(SFIRE_PHYS_CRIT)
                  end if

               end if
            end do
         end do
      end do

!print*,'ADVANCE A6'
!call flush(6)
      change = moisture_dt/(fmep_decay_tlag*3600.)
      print *,'LFR-DBG: change,tol:',change,tol
      if (change < tol) then
         do j = jts, jte
            do k = 1, 2
               do i = its, ite
                  fmep(i, k, j) = fmep(i, k, j)*(1.0 - change*(1.0 - 0.5*change))
               end do 
            end do
         end do
      else
         do j = jts, jte
            do k = 1, 2
               do i = its, ite
                  fmep(i, k, j) = fmep(i, k, j)*exp(-change)
               end do 
            end do
         end do
      end if



!LFR      do j = jts, jte
!LFR         do k = 1, 2
!LFR            do i = its, ite
!LFR               change = moisture_dt/(fmep_decay_tlag*3600.)
!LFR               if (change < tol) then
!LFR                  fmep(i, k, j) = fmep(i, k, j)*(1.0 - change*(1.0 - 0.5*change))
!LFR                  !       print*,'ADVANCE fmep',  fmep(i,k,j)
!LFR                  !       call flush(6)
!LFR               else
!LFR                  fmep(i, k, j) = fmep(i, k, j)*exp(-change)
!LFR                  !      print*,'ADVANCE fmep 2',  fmep(i,k,j)
!LFR                  !      call flush(6)
!LFR               end if
!LFR            end do
!LFR         end do
!LFR      end do
!print*,"ADVANCE A4"
!call flush(6)
      if (fire_print_msg .ge. 2) then
!$OMP CRITICAL(SFIRE_PHYS_CRIT)
         ! print*,"ADVANCE A4-a"
         ! call flush(6)
         write (msg, 2) 'Rain intensity    min', rainmin, ' max', rainmax, ' mm/h'
         ! print*,"ADVANCE A4-a1"
         ! call flush(6)
         call message(msg)
         ! print*,"ADVANCE A4-a2"
         ! call flush(6)
         if (rainmin < 0.) then
            call message('WARNING rain accumulation must increase')
         end if
         ! print*,"ADVANCE A4-b"
         ! call flush(6)
         write (msg, 2) 'Relative humidity min', 100*rhmin, ' max', 100*rhmax, '%'
         ! print*,"ADVANCE A4-b1"
         ! call flush(6)
         call message(msg)
         ! print*,"ADVANCE A4-c"
         ! call flush(6)
         if (.not. (rhmax <= 1.0 .and. rhmin >= 0)) then
            call message('WARNING Relative humidity must be between 0 and 100%')
         end if
         ! print*,"ADVANCE A4-c1"
         ! call flush(6)
2        format(2(a, f10.2), a)
!$OMP END CRITICAL(SFIRE_PHYS_CRIT)
      end if
!print*,"vou para o print_3d_stats_by_slice A2"
!call flush(6)
      call print_3d_stats_by_slice(its, ite, 1, moisture_classes, jts, jte, ims, ime, 1, nfmc, jms, jme, fmc_equi &
                                   , 'equilibrium fmc_equi')
!print*,"sai print_3d_stats_by_slice A2"
!call flush(6)
!print*,"vou para o print_3d_stats_by_slice A3"
!call flush(6)
      call print_3d_stats_by_slice(its, ite, 1, moisture_classes, jts, jte, ims, ime, 1, nfmc, jms, jme, fmc_lag, 'time lag')
!print*,"sai print_3d_stats_by_slice A3"
!call flush(6)
!print*,"vou para o print_3d_stats_by_slice A4"
!call flush(6)
      call print_3d_stats_by_slice(its, ite, 1, moisture_classes, jts, jte, ims, ime, 1, nfmc, jms, jme, fmc_gc &
                                   , 'after advance fmc_gc')
!print*,"sai print_3d_stats_by_slice A4"
!call flush(6)
      call copy2old
      print *, "advance moiture vou sair da rotina"
      call flush (6)
      return

   contains

      subroutine copy2old

         do j = jts, jte
            do i = its, ite
               rain_old(i, j) = rainc(i, j) + rainnc(i, j)
               t2_old(i, j) = t2(i, j)
               q2_old(i, j) = q2(i, j)
               psfc_old(i, j) = psfc(i, j)
            end do
         end do

      end subroutine copy2old

      subroutine get_equi_moist
      end subroutine get_equi_moist

   end subroutine advance_moisture

   subroutine init_fuel_cats(ifun, init_fuel_moisture, config_flags) !config_flags introduzido por Isilda CM

      implicit none

      logical, intent(in)::init_fuel_moisture
      logical, external:: wrf_dm_on_monitor

      integer:: i, j, k, ii, iounit, ierr, kk
      character(len=128):: msg
      real, dimension(mfuelcats) :: fwh, fz0

! #### INTRODUZIDO POR ISILDA CUNHA MENEZES
      type(grid_config_rec_type), intent(IN)  :: config_flags
      logical::crown
      integer, intent(in):: ifun
!####

      namelist /fuel_scalars/ cmbcnst, hfgl, fuelmc_g, fuelmc_c, nfuelcats, no_fuel_cat, no_fuel_cat2 &
         , fire_wind_height, ibeh
      namelist /fuel_categories/ fuel_name, windrf, fgi, fueldepthm, savr, &
         fuelmce, fueldens, st, se, weight, fci_d, fct, ichap, itree, fwh, fz0, ffw, &
         fmc_gl_ndwi_0, fmc_gl_ndwi_rate, fmc_gl_ndwi_stdev, fmc_gl_stdev, &
         adjr0, adjrw, adjrs, fmc_gw01, fmc_gw02, fmc_gw03, fmc_gw04, fmc_gw05

      namelist /moisture/ moisture_classes, drying_lag, wetting_lag, saturation_moisture &
         , saturation_rain, rain_threshold, &
         drying_model, wetting_model, moisture_class_name, fmc_gc_initialization

!!!!AUTOR ISILDA CUNHA MENEZES
      namelist /fuel_crown/ CBH, FMC_crown

      !print *, "estou no modulo_fr_drive_phys e estou na rotina init_fuel_cats"
      !call flush (6)

      crown = config_flags%crown
!!!!!
!IF ( wrf_dm_on_monitor() ) THEN

      fwh = fcwh
      fz0 = fcz0

      iounit = open_text_file('namelist.fire', 'read')
      read (iounit, fuel_scalars, iostat=ierr)
      if (ierr .ne. 0) then
         call message('init_fuel_cats: error reading namelist fuel_scalars in file namelist.fire', 0)
         call error_namelist(iounit)
      end if
      read (iounit, fuel_categories, iostat=ierr)
      if (ierr .ne. 0) then
         call message('init_fuel_cats: error reading namelist fuel_categories in file namelist.fire')
         call error_namelist(iounit)
      end if
      if (init_fuel_moisture) then
         read (iounit, moisture, iostat=ierr)
         if (ierr .ne. 0) then
            call message('init_fuel_cats: error reading namelist moisture in file namelist.fire')
            call error_namelist(iounit)
         end if
      end if
!!!!AUTOR ISILDA CUNHA MENEZES
      read (iounit, fuel_crown, iostat=ierr)
!    print*,"LI FUEL_CROWN",iounit,ierr
!    call flush(6)
      if (ierr .ne. 0) then
         call message('init_fuel_cats: error reading namelist fuel_crown in file namelist.fire')
         call error_namelist(iounit)
      end if
!!!!!!!

      fmc_gw(1:mfuelcats, 1) = fmc_gw01
      fmc_gw(1:mfuelcats, 2) = fmc_gw02
      fmc_gw(1:mfuelcats, 3) = fmc_gw03
      fmc_gw(1:mfuelcats, 4) = fmc_gw04
      fmc_gw(1:mfuelcats, 5) = fmc_gw05

      close (iounit)

      fcwh = fwh
      fcz0 = fz0

      if (nfuelcats > mfuelcats) then
         write (msg, *) 'nfuelcats=', nfuelcats, ' too large, increase mfuelcats'
         call crash(msg)
      end if
      if (no_fuel_cat >= 1 .and. no_fuel_cat <= nfuelcats) then
         write (msg, *) 'no_fuel_cat=', no_fuel_cat, ' may not be between 1 and nfuelcats=', nfuelcats
         call crash(msg)
      end if
      if (no_fuel_cat > no_fuel_cat2) then
         write (msg, *) 'no_fuel_cat=', no_fuel_cat, ' must not be larger than no_fuel_cat2=', no_fuel_cat2
         call crash(msg)
      end if
!ENDIF

!call wrf_dm_bcast_real(cmbcnst,1)
!call wrf_dm_bcast_real(hfgl,1)
!call wrf_dm_bcast_real(fuelmc_g,1)
!call wrf_dm_bcast_real(fuelmc_c,1)
!call wrf_dm_bcast_real(fire_wind_height,1)
!call wrf_dm_bcast_integer(nfuelcats,1)
!call wrf_dm_bcast_integer(no_fuel_cat,1)
!call wrf_dm_bcast_integer(no_fuel_cat2,1)
!call wrf_dm_bcast_integer(ibeh,1)
!call wrf_dm_bcast_real(windrf,    nfuelcats)
!call wrf_dm_bcast_real(fgi,       nfuelcats)
!call wrf_dm_bcast_real(fueldepthm,nfuelcats)
!call wrf_dm_bcast_real(savr,      nfuelcats)
!call wrf_dm_bcast_real(fuelmce,   nfuelcats)
!call wrf_dm_bcast_real(fueldens,  nfuelcats)
!call wrf_dm_bcast_real(st,        nfuelcats)
!call wrf_dm_bcast_real(se,        nfuelcats)
!call wrf_dm_bcast_real(weight,    nfuelcats)
!call wrf_dm_bcast_real(fci_d,     nfuelcats)
!call wrf_dm_bcast_real(fct,       nfuelcats)
!call wrf_dm_bcast_integer(ichap,  nfuelcats)
!call wrf_dm_bcast_integer(itree,  nfuelcats)
!call wrf_dm_bcast_real(fcwh,      nfuelcats)
!call wrf_dm_bcast_real(fcz0,      nfuelcats)
!call wrf_dm_bcast_real(ffw,       nfuelcats)
!call wrf_dm_bcast_real(adjr0,     nfuelcats)
!call wrf_dm_bcast_real(adjrw,     nfuelcats)
!call wrf_dm_bcast_real(adjrs,     nfuelcats)
!call wrf_dm_bcast_real(fmc_gl_ndwi_0,    nfuelcats)
!call wrf_dm_bcast_real(fmc_gl_ndwi_rate, nfuelcats)
!call wrf_dm_bcast_real(fmc_gl_ndwi_stdev,nfuelcats)
!call wrf_dm_bcast_real(fmc_gl_stdev,     nfuelcats)

!call wrf_dm_bcast_integer(moisture_classes,1)
!call wrf_dm_bcast_real(drying_lag,     max_moisture_classes)
!call wrf_dm_bcast_real(wetting_lag,     max_moisture_classes)
!call wrf_dm_bcast_real(saturation_moisture,     max_moisture_classes)
!call wrf_dm_bcast_real(saturation_rain,     max_moisture_classes)
!call wrf_dm_bcast_real(rain_threshold,     max_moisture_classes)
!call wrf_dm_bcast_integer(drying_model,     max_moisture_classes)
!call wrf_dm_bcast_integer(wetting_model,     max_moisture_classes)
!call wrf_dm_bcast_integer(fmc_gc_initialization,     max_moisture_classes)
!call wrf_dm_bcast_real(fmc_gw,     mfuelcats*max_moisture_classes)

      do i = 1, moisture_classes
         rec_drying_lag_sec(i) = 1.0/(3600.0*drying_lag(i))
         rec_wetting_lag_sec(i) = 1.0/(3600.0*wetting_lag(i))
      end do

      fuelheat = cmbcnst*4.30e-04

      do i = 1, nfuelcats
         fci(i) = (1.+fuelmc_c)*fci_d(i)
         if (fct(i) .ne. 0.) then
            fcbr(i) = fci_d(i)/fct(i)
         else
            fcbr(i) = 0
         end if
      end do

      call message('**********************************************************')
      call message('FUEL COEFFICIENTS')
      write (msg, 8) 'cmbcnst    ', cmbcnst
      call message(msg)
      write (msg, 8) 'hfgl       ', hfgl
      call message(msg)
      write (msg, 8) 'fuelmc_g   ', fuelmc_g
      call message(msg)
      write (msg, 8) 'fuelmc_c   ', fuelmc_c
      call message(msg)
      write (msg, 8) 'fuelheat   ', fuelheat
      call message(msg)
      write (msg, 7) 'nfuelcats  ', nfuelcats
      call message(msg)
      write (msg, 7) 'no_fuel_cat', no_fuel_cat
      call message(msg)
      write (msg, 7) 'no_fuel_cat2', no_fuel_cat2
      call message(msg)
      if (init_fuel_moisture) then
         write (msg, 7) 'moisture_classes', moisture_classes
         call message(msg)
      end if

      j = 1
7     format(a, 5(i9, 4x))
8     format(a, 5(1x, g12.5e2))
9     format(a, 5(1x, a))
10    format(a, i2.2, 2x, 5(1x, g12.5e2))
      do i = 1, nfuelcats, j
         k = min(i + j - 1, nfuelcats)
         call message(' ')
         write (msg, 7) 'CATEGORY  ', (ii, ii=i, k)
         call message(msg)
         write (msg, 9) 'fuel name ', (trim(fuel_name(ii)), ii=i, k)
         call message(msg)
         write (msg, 8) 'fwh       ', (fcwh(ii), ii=i, k)
         call message(msg)
         write (msg, 8) 'fz0       ', (fcz0(ii), ii=i, k)
         call message(msg)
         write (msg, 8) 'windrf    ', (windrf(ii), ii=i, k)
         call message(msg)
         write (msg, 8) 'fgi       ', (fgi(ii), ii=i, k)

         !   print*,"FGI_A",(fgi(ii),ii=i,k)
         !   call flush(6)
         call message(msg)
         write (msg, 8) 'fueldepthm', (fueldepthm(ii), ii=i, k)
         call message(msg)
         write (msg, 8) 'savr      ', (savr(ii), ii=i, k)
         call message(msg)
         write (msg, 8) 'fuelmce   ', (fuelmce(ii), ii=i, k)
         call message(msg)
         write (msg, 8) 'fueldens  ', (fueldens(ii), ii=i, k)
         call message(msg)
         write (msg, 8) 'st        ', (st(ii), ii=i, k)
         call message(msg)
         write (msg, 8) 'se        ', (se(ii), ii=i, k)
         call message(msg)
         write (msg, 8) 'weight    ', (weight(ii), ii=i, k)
         call message(msg)
         write (msg, 8) 'fci_d     ', (fci_d(ii), ii=i, k)
         call message(msg)
         write (msg, 8) 'fct       ', (fct(ii), ii=i, k)
         call message(msg)
         write (msg, 7) 'ichap     ', (ichap(ii), ii=i, k)
         call message(msg)
         write (msg, 7) 'itree     ', (itree(ii), ii=i, k)
         call message(msg)
         write (msg, 8) 'fci       ', (fci(ii), ii=i, k)
         call message(msg)
         write (msg, 8) 'fcbr      ', (fcbr(ii), ii=i, k)
         call message(msg)
         write (msg, 8) 'ffw       ', (ffw(ii), ii=i, k)
         call message(msg)
         write (msg, 8) 'adjr0     ', (adjr0(ii), ii=i, k)
         call message(msg)
         write (msg, 8) 'adjrw     ', (adjrw(ii), ii=i, k)
         call message(msg)
         write (msg, 8) 'adjrs     ', (adjrs(ii), ii=i, k)
         call message(msg)
         if (init_fuel_moisture) then
            do kk = 1, moisture_classes
               write (msg, 10) 'fmc_gw', kk, (fmc_gw(ii, kk), ii=i, k)
               call message(msg)
            end do
         end if
         if (kfmc_ndwi > 0) then
            write (msg, 8) 'fmc_gl_stdev     ', (fmc_gl_stdev(ii), ii=i, k)
            call message(msg)
            write (msg, 8) 'fmc_gl_ndwi_0    ', (fmc_gl_ndwi_0(ii), ii=i, k)
            call message(msg)
            write (msg, 8) 'fmc_gl_ndwi_rate ', (fmc_gl_ndwi_rate(ii), ii=i, k)
            call message(msg)
            write (msg, 8) 'fmc_gl_ndwi_stdev', (fmc_gl_ndwi_stdev(ii), ii=i, k)
            call message(msg)
         end if
      end do
      call message(' ')
      call message('**********************************************************')

      if (init_fuel_moisture) then
         j = 1
         do i = 1, moisture_classes, j
            k = min(i + j - 1, nfuelcats)
            call message(' ')
            write (msg, 7) 'FUEL MOISTURE CLASS', (ii, ii=i, k)
            call message(msg)
            write (msg, 9) 'moisture class name    ', (trim(moisture_class_name(ii)), ii=i, k)
            call message(msg)
            write (msg, 7) 'drying_model           ', (drying_model(ii), ii=i, k)
            call message(msg)
            write (msg, 8) 'drying_lag (h)         ', (drying_lag(ii), ii=i, k)
            call message(msg)
            write (msg, 7) 'wetting_model          ', (wetting_model(ii), ii=i, k)
            call message(msg)
            write (msg, 7) 'fmc_gc_initialization  ', (fmc_gc_initialization(ii), ii=i, k)
            call message(msg)
            write (msg, 8) 'wetting_lag (h)        ', (wetting_lag(ii), ii=i, k)
            call message(msg)
            write (msg, 8) 'saturation_moisture (1)', (saturation_moisture(ii), ii=i, k)
            call message(msg)
            write (msg, 8) 'saturation_rain (mm/h) ', (saturation_rain(ii), ii=i, k)
            call message(msg)
            write (msg, 8) 'rain_threshold (mm/h)  ', (rain_threshold(ii), ii=i, k)
            call message(msg)
         end do
         call message(' ')
         call message('**********************************************************')
         call message(' ')
      end if

!!!!! AUTOR ISILDA CUNHA MENEZES
      if (crown) then
      do i = 1, nfuelcats, j
         k = min(i + j - 1, nfuelcats)
         call message(' ')
         write (msg, 8) 'CBH       ', (CBH(ii), ii=i, k)
         !    print*,"CBH_A",(CBH(ii),ii=i,k)
         !    call flush(6)
         call message(msg)
         write (msg, 8) 'FMC_crown       ', (FMC_crown(ii), ii=i, k)
         !    print*,"FMC_CROWN_A",(FMC_crown(ii),ii=i,k)
         !    call flush(6)
         call message(msg)
      end do
      end if
!!!!!!!!!

      have_fuel_cats = .true.

!IF ( wrf_dm_on_monitor() ) THEN
      call write_fuels_m(ifun, 61, 30., 1., config_flags) !config_flags introduzido por ISILDA CM
      !print *, "vim do driver e estou no modulo_drive_phys e estou na rotina init_fuel_cats e chamei a rotina write_fuels_m"
!ENDIF

   end subroutine init_fuel_cats

   subroutine write_fuels_m(ifun, nsteps, maxwind, maxslope, config_flags) !config_flags e ifun introduzido por Isilda CM

      implicit none
      integer, intent(in):: nsteps
      real, intent(in):: maxwind, maxslope

      integer:: iounit, k, j, i, isave
      type(fire_params)::fp
      real, dimension(1:3, 1:nsteps), target::vx, vy, zsf, dzdxf, dzdyf, bbb, phisc, phiwc, r_0, fgip, ischap, fmc_g &
                                               , wind, nfuel_cat, bbb_FM10, phisc_FM10, phiwc_FM10, r_0_FM10, CFB, HPA ! INTRODUZIDO POR ISILDA CM
      real, dimension(1:3, 1:nsteps)::fuel_time, ros, fwh, fz0
      real::ros_back, ros_wind, ros_slope, propx, propy, r
      real:: ros_total !!!! INTRODUZIDO POR ISILDA CUNHA MENEZES
      integer, intent(in)::ifun !!!INTRODUZIDO POR ISILDA DA CUNHA MENEZES
      integer::ierrx
      character(len=128)::msg

! ## Introduzido por ISILDA CUNHA MENEZES
      type(grid_config_rec_type), intent(IN)  :: config_flags
      logical::crown

      !print *, "estou no modulo_fr_drive_phys e estou na rotina write_fuels_m"
      !call flush (6)
      crown = config_flags%crown

!!!!####

      if (.not. have_fuel_cats) call crash('write_fuels_m: fuel categories not yet set')

      fp%vx => vx
      fp%vy => vy
      fp%dzdxf => dzdxf
      fp%dzdyf => dzdyf
      fp%bbb => bbb
      fp%phisc => phisc
      fp%phiwc => phiwc
      fp%r_0 => r_0
      fp%fgip => fgip
      fp%ischap => ischap
      fp%fmc_g => fmc_g
      fp%nfuel_cat => nfuel_cat
!!!INTRODUZIDO POR ISILDA CUNHA MENEZES
      fp%bbb_FM10 => bbb_FM10
      fp%phisc_FM10 => phisc_FM10
      fp%phiwc_FM10 => phiwc_FM10
      fp%r_0_FM10 => r_0_FM10
      fp%CFB => CFB
      fp%HPA => HPA
!!!!

      iounit = open_text_file('fuels.m', 'write')

10    format('fuel(', i3, ').', a, '=', "'", a, "'", ';% ', a)
      do k = 1, nfuelcats
         write (iounit, 10) k, 'fuel_name', trim(fuel_name(k)), 'FUEL MODEL NAME'
         call write_var(k, 'windrf', windrf(k), 'WIND REDUCTION FACTOR FROM FCWH TO MIDFLAME HEIGHT')
         call write_var(k, 'fwh', fcwh(k), 'WIND HEIGHT TO INTERPOLATE VERTICALLY TO (M)')
         call write_var(k, 'fz0', fcz0(k), 'ROUGHNESS LENGTH FOR VERTICAL WIND LOG INTERPOLATION (M)')
         call write_var(k, 'fgi', fgi(k), 'INITIAL TOTAL MASS OF SURFACE FUEL (KG/M**2)')
         call write_var(k, 'fueldepthm', fueldepthm(k), 'FUEL DEPTH (M)')
         call write_var(k, 'savr', savr(k), 'FUEL PARTICLE SURFACE-AREA-TO-VOLUME RATIO, 1/FT')
         call write_var(k, 'fuelmce', fuelmce(k), 'MOISTURE CONTENT OF EXTINCTION')
         call write_var(k, 'fueldens', fueldens(k), 'OVENDRY PARTICLE DENSITY, LB/FT^3')
         call write_var(k, 'st', st(k), 'FUEL PARTICLE TOTAL MINERAL CONTENT')
         call write_var(k, 'se', se(k), 'FUEL PARTICLE EFFECTIVE MINERAL CONTENT')
         call write_var(k, 'weight', weight(k), 'WEIGHTING PARAMETER THAT DETERMINES THE SLOPE OF THE MASS LOSS CURVE')
         call write_var(k, 'fci_d', fci_d(k), 'INITIAL DRY MASS OF CANOPY FUEL')
         call write_var(k, 'fct', fct(k), 'BURN OUT TIME FOR CANOPY FUEL, AFTER DRY (S)')
         call write_var(k, 'ichap', float(ichap(k)), '1 if chaparral, 0 if not')
         call write_var(k, 'fci', fci(k), 'INITIAL TOTAL MASS OF CANOPY FUEL')
         call write_var(k, 'fcbr', fcbr(k), 'FUEL CANOPY BURN RATE (KG/M**2/S)')
         call write_var(k, 'adjr0', adjr0(k), 'MULTIPLICATIVE ADJUSTMENT OF BACKING SPREAD RATE')
         call write_var(k, 'adjrw', adjrw(k), 'MULTIPLICATIVE ADJUSTMENT OF WIND CONTRIBUTION TO SPREAD RATE')
         call write_var(k, 'adjrs', adjrs(k), 'MULTIPLICATIVE ADJUSTMENT OF SLOPE CONTRIBUTION TO SPREAD RATE')
         call write_var(k, 'ffw', ffw(k), 'FUEL FRACTION CONSUMED IN THE FLAMING ZONE')
         call write_var(k, 'hfgl', hfgl, 'SURFACE FIRE HEAT FLUX THRESHOLD TO IGNITE CANOPY (W/m^2)')
         call write_var(k, 'cmbcnst', cmbcnst, 'JOULES PER KG OF DRY FUEL')
         call write_var(k, 'fuelheat', fuelheat, 'FUEL PARTICLE LOW HEAT CONTENT, BTU/LB')
         call write_vaR(k, 'fuelmc_g', fuelmc_g, 'FUEL PARTICLE (SURFACE) MOISTURE CONTENT')
         call write_var(k, 'fuelmc_c', fuelmc_c, 'FUEL PARTICLE (CANOPY) MOISTURE CONTENT')

 !!!!!INTRODUZIDO POR ISILDA CUNHA MENEZES
         if (crown) then
            call write_var(k, 'CBH', CBH(k), 'Canopy base heigt (m)')
            call write_var(k, 'FMC_crown', FMC_crown(k), 'CROWN FOLIAR MOISTURE CONTENT, (%)')

         end if
!!!!!!

         nfuel_cat = k
         do j = 1, nsteps
            fmc_g(1, j) = fuelmc_g
            fmc_g(2, j) = fuelmc_g
            fmc_g(3, j) = (fuelmce(k)*(j - 1))/(nsteps - 2)
         end do
         isave = fire_fmc_read
         fire_fmc_read = 0
         call set_fire_params( &
            1, 3, 1, nsteps, &
            1, 3, 1, nsteps, &
            1, 3, 1, nsteps, &
            0., 0., k, &
            nfuel_cat, fuel_time, &
            fp)
         fire_fmc_read = isave

         propx = 1.
         propy = 0.
         do j = 1, nsteps
            r = float(j - 1)/float(nsteps - 1)

            wind(1, j) = maxwind*r
            vx(1, j) = wind(1, j)*windrf(k)
            vy(1, j) = 0.
            dzdxf(1, j) = 0.
            dzdyf(1, j) = 0.

            vx(2, j) = 0.
            vy(2, j) = 0.
            dzdxf(2, j) = maxslope*r
            dzdyf(2, j) = 0.

            vx(3, j) = 0.
            vy(3, j) = 0.
            dzdxf(3, j) = 0.
            dzdyf(3, j) = 0.
         end do
         do j = 1, nsteps
            do i = 1, 3
              !!!!! INTRODUZIDO POR ISILDA CUNHA MENEZES
               call fire_ros(ifun, ros_back, ros_wind, ros_slope, &
                             propx, propy, i, j, fp, ierrx, msg, config_flags, ros_total) !config_flags e ifun introduzido por ISILDA CM
               ros(i, j) = ros_total
               !          print*, "estou no modulo_fr_phys na rotina  write_fuels_m e chamei a rotina fire_ros para escrever a var ros_total"
               !          call flush(6)
                  !!!!!
               !  call fire_ros(ros_back,ros_wind,ros_slope, &
               !     propx,propy,i,j,fp,ierrx,msg)
               !  ros(i,j)=ros_back+ros_wind+ros_slope
            end do
            write (iounit, 13) k, 'wind', j, wind(1, j), 'wind speed at 6.1m'
            write (iounit, 13) k, 'ros_wind', j, ros(1, j), 'rate of spread for the wind speed at 6.1m'
            write (iounit, 13) k, 'slope', j, dzdxf(2, j), 'slope'
            write (iounit, 13) k, 'ros_slope', j, ros(2, j), 'rate of spread for the slope'
            write (iounit, 13) k, 'fmc_g', j, fmc_g(3, j), 'fuel moisture content'
            write (iounit, 13) k, 'ros_fmc_g', j, ros(3, j), 'rate of spread for the fuel moisture content'
         end do
      end do
13    format('fuel(', i3, ').', a, '(', i3, ')=', g12.5e2, ';% ', a)

      close (iounit)

   contains

      subroutine write_var(k, name, value, descr)

         integer, intent(in)::k
         character(len=*), intent(in)::name, descr
         real, intent(in)::value
         write (iounit, 11) k, name, value
         write (iounit, 12) k, name, descr
11       format('fuel(', i3, ').', a, '=', g12.5e2, ';')
12       format('fuel(', i3, ').', a, "_descr='", a, "';")
      end subroutine write_var

   end subroutine write_fuels_m

   subroutine set_fire_params( &
      ifds, ifde, jfds, jfde, &
      ifms, ifme, jfms, jfme, &
      ifts, ifte, jfts, jfte, &
      fdx, fdy, nfuel_cat0, &
      nfuel_cat, fuel_time, &
      fp)

      implicit none

      integer, intent(in)::ifds, ifde, jfds, jfde
      integer, intent(in)::ifts, ifte, jfts, jfte
      integer, intent(in)::ifms, ifme, jfms, jfme
      real, intent(in):: fdx, fdy
      integer, intent(in)::nfuel_cat0
      real, intent(inout), dimension(ifms:ifme, jfms:jfme)::nfuel_cat
      real, intent(out), dimension(ifms:ifme, jfms:jfme)::fuel_time
      type(fire_params), intent(inout)::fp

      real::  fuelload, fueldepth, rtemp1, rtemp2, &
             qig, epsilon, rhob, wn, betaop, e, c, &
             xifr, etas, etam, a, gammax, gamma, ratio, ir, &
             fuelloadm, fdxinv, fdyinv, betafl, bmst
      integer:: i, j, k
      integer::nerr
      character(len=128)::msg

      !print *, "estou no modulo_fr_drive_phys e estou na rotina set_fire_params"
      !call flush (6)
      if (.not. have_fuel_cats) call crash('set_fire_params: fuel categories not yet set')

      call print_2d_stats(ifts, ifte, jfts, jfte, ifms, ifme, jfms, jfme, fp%fmc_g, 'set_fire_params: fmc_g')

      nerr = 0
      do j = jfts, jfte
         do i = ifts, ifte

            k = int(nfuel_cat(i, j))
            if (k .ge. no_fuel_cat .and. k .le. no_fuel_cat2) then
               fp%fgip(i, j) = 0.
               fp%ischap(i, j) = 0.
               fp%phisc(i, j) = 0.
               fp%bbb(i, j) = 1.
               fuel_time(i, j) = 7./0.85
               fp%phiwc(i, j) = 0.
               fp%r_0(i, j) = 0.
            else

               if (k .lt. 1 .or. k .gt. nfuelcats) then
!$OMP CRITICAL(SFIRE_PHYS_CRIT)
                  write (msg, '(3(a,i5))') 'nfuel_cat(', i, ',', j, ')=', k
                  call message(msg, level=0)
                  if (k .eq. 0) then
                     call message('Possibly nfuel_cat is uninitialized on input', level=0)
                  end if
                  write (msg, '(a,i5)') 'nfuel_cat must between 1 and nfuel_cat=', nfuelcats
                  call message(msg, level=0)
                  write (msg, '(a,i5,a,i10,a)') ' or ', no_fuel_cat, ' and', no_fuel_cat2, ' for no fuel'
                  call message(msg, level=0)
!$OMP END CRITICAL(SFIRE_PHYS_CRIT)
                  call crash('set_fire_params: fuel category out of bounds')
               end if

               fuel_time(i, j) = weight(k)/0.85

               fp%ischap(i, j) = ichap(k)
               fp%fgip(i, j) = fgi(k)
               if (fire_fmc_read .eq. 1) then
                  fp%fmc_g(i, j) = fuelmc_g
               end if

               bmst = fp%fmc_g(i, j)/(1.+fp%fmc_g(i, j))
               fuelloadm = (1.-bmst)*fgi(k)
               fuelload = fuelloadm*(.3048)**2*2.205
               fueldepth = fueldepthm(k)/0.3048
               betafl = fuelload/(fueldepth*fueldens(k))
               betaop = 3.348*savr(k)**(-0.8189)
               qig = 250.+1116.*fp%fmc_g(i, j)
               epsilon = exp(-138./savr(k))
               rhob = fuelload/fueldepth

               c = 7.47*exp(-0.133*savr(k)**0.55)
               fp%bbb(i, j) = 0.02526*savr(k)**0.54

               e = 0.715*exp(-3.59e-4*savr(k))
               fp%phiwc(i, j) = c*(betafl/betaop)**(-e)

               fp%phisc(i, j) = 5.275*(betafl)**(-0.3)

               rtemp2 = savr(k)**1.5
               gammax = rtemp2/(495.+0.0594*rtemp2)
               a = 1./(4.774*savr(k)**0.1 - 7.27)
               ratio = betafl/betaop
               gamma = gammax*(ratio**a)*exp(a*(1.-ratio))

               wn = fuelload/(1 + st(k))
               rtemp1 = fp%fmc_g(i, j)/fuelmce(k)
               etam = 1.-2.59*rtemp1 + 5.11*rtemp1**2 - 3.52*rtemp1**3
               etam = max(etam, 0.0)
               etas = 0.174*se(k)**(-0.19)
               ir = gamma*wn*fuelheat*etam*etas

               xifr = exp((0.792 + 0.681*savr(k)**0.5) &
                          *(betafl + 0.1))/(192.+0.2595*savr(k))

               fp%r_0(i, j) = ir*xifr/(rhob*epsilon*qig)
               fp%r_0(i, j) = fp%r_0(i, j)*.00508
               fp%phiwc(i, j) = fp%phiwc(i, j)*fp%r_0(i, j)
               fp%phisc(i, j) = fp%phisc(i, j)*fp%r_0(i, j)

               fp%r_0(i, j) = fp%r_0(i, j)*adjr0(k)
               fp%phiwc(i, j) = fp%phiwc(i, j)*adjrw(k)
               fp%phisc(i, j) = fp%phisc(i, j)*adjrs(k)

               if (fp%r_0(i, j) > 1e-6 .and. fp%fmc_g(i, j) > fuelmce(k)) then
!$OMP CRITICAL(SFIRE_PHYS_CRIT)
                  write (msg, '(a,2i6,3(a,e14.6))') 'set_fire_params: at ', i, j, ' base rate of spread', fp%r_0(i, j) &
                     , ' moisture ', fp%fmc_g(i, j), '> extinction =', fuelmce(k)
                  call message(msg, level=0)
                  write (msg, '(5(a,e14.6))') 'rtemp1=', rtemp1, ' etam=', etam
                  call message(msg, level=0)
!$OMP END CRITICAL(SFIRE_PHYS_CRIT)
                  call warning('propagation above extinction moisture', level=0)
               end if
            end if
         end do
      end do

      if (nerr .gt. 1) then
!$OMP CRITICAL(SFIRE_PHYS_CRIT)
         write (msg, '(a,i6,a)') 'set_fire_params: WARNING: fuel category 0 replaced in', nerr, ' cells'
!$OMP END CRITICAL(SFIRE_PHYS_CRIT)
         call message(msg)
      end if

   end subroutine set_fire_params

!!!!!! INTRODUZIDO POR ISILDA CUNHA MENEZES

   subroutine set_fire_crown_params( &
      ifds, ifde, jfds, jfde, &
      ifms, ifme, jfms, jfme, &
      ifts, ifte, jfts, jfte, &
      fdx, fdy, nfuel_cat0, &
      nfuel_cat, fuel_time, &
      fp)

      implicit none

      integer, intent(in)::ifds, ifde, jfds, jfde
      integer, intent(in)::ifts, ifte, jfts, jfte
      integer, intent(in)::ifms, ifme, jfms, jfme
      real, intent(in):: fdx, fdy
      integer, intent(in)::nfuel_cat0
      real, intent(inout), dimension(ifms:ifme, jfms:jfme)::nfuel_cat
      real, intent(out), dimension(ifms:ifme, jfms:jfme)::fuel_time
      type(fire_params), intent(inout)::fp

      real::  fuelload, fueldepth, rtemp1, rtemp2, &
             qig, epsilon, rhob, wn, betaop, e, c, &
             xifr, etas, etam, a, gammax, gamma, ratio, ir, &
             fuelloadm, fdxinv, fdyinv, betafl, bmst
      integer:: i, j, k
      integer::nerr
      character(len=128)::msg

      !print*,"estou no modulo_fr_drive_phys e estou na rotina set_fire_crown_params",have_fuel_cats
      !call flush(6)
      if (.not. have_fuel_cats) call crash('set_fire_params: fuel categories &
 &      notyet set')

      call print_2d_stats(ifts, ifte, jfts, jfte, ifms, ifme, jfms, jfme, fp%fmc_g, &
                          'set_fire_params:fmc_g')

      nerr = 0
      do j = jfts, jfte
         do i = ifts, ifte

            k = 10

            fuel_time(i, j) = weight(k)/0.85

            fp%ischap(i, j) = ichap(k)
            fp%fgip(i, j) = fgi(k)
            if (fire_fmc_read .eq. 1) then
               fp%fmc_g(i, j) = fuelmc_g
            end if

            bmst = fp%fmc_g(i, j)/(1.+fp%fmc_g(i, j))
            fuelloadm = (1.-bmst)*fgi(k)
            fuelload = fuelloadm*(.3048)**2*2.205
            fueldepth = fueldepthm(k)/0.3048
            betafl = fuelload/(fueldepth*fueldens(k))
            betaop = 3.348*savr(k)**(-0.8189)
            qig = 250.+1116.*fp%fmc_g(i, j)
            epsilon = exp(-138./savr(k))
            rhob = fuelload/fueldepth

            c = 7.47*exp(-0.133*savr(k)**0.55)
            fp%bbb_FM10(i, j) = 0.02526*savr(k)**0.54

            e = 0.715*exp(-3.59e-4*savr(k))
            fp%phiwc_FM10(i, j) = c*(betafl/betaop)**(-e)

            fp%phisc_FM10(i, j) = 5.275*(betafl)**(-0.3)

            rtemp2 = savr(k)**1.5
            gammax = rtemp2/(495.+0.0594*rtemp2)
            a = 1./(4.774*savr(k)**0.1 - 7.27)
            ratio = betafl/betaop
            gamma = gammax*(ratio**a)*exp(a*(1.-ratio))

            wn = fuelload/(1 + st(k))
            rtemp1 = fp%fmc_g(i, j)/fuelmce(k)
            etam = 1.-2.59*rtemp1 + 5.11*rtemp1**2 - 3.52*rtemp1**3
            etam = max(etam, 0.0)
            etas = 0.174*se(k)**(-0.19)
            ir = gamma*wn*fuelheat*etam*etas

            xifr = exp((0.792 + 0.681*savr(k)**0.5) &
                       *(betafl + 0.1))/(192.+0.2595*savr(k))

            fp%r_0_FM10(i, j) = ir*xifr/(rhob*epsilon*qig)
            !print *,'LFR-DBG: fp%r_0_FM10(',i,',',',',j,'), ir,xifr,resto =',fp%r_0_FM10(i, j),ir,xifr,(rhob*epsilon*qig)
            fp%r_0_FM10(i, j) = fp%r_0_FM10(i, j)*.00508
            fp%phiwc_FM10(i, j) = fp%phiwc_FM10(i, j)*fp%r_0_FM10(i, j)
            fp%phisc_FM10(i, j) = fp%phisc_FM10(i, j)*fp%r_0_FM10(i, j)

            fp%r_0_FM10(i, j) = fp%r_0_FM10(i, j)*adjr0(k)
            fp%phiwc_FM10(i, j) = fp%phiwc_FM10(i, j)*adjrw(k)
            fp%phisc_FM10(i, j) = fp%phisc_FM10(i, j)*adjrs(k)

            if (fp%r_0_FM10(i, j) > 1e-6 .and. fp%fmc_g(i, j) > fuelmce(k)) then
!$OMP CRITICAL(SFIRE_PHYS_CRIT)
               write (msg, '(a,2i6,3(a,e14.6))') 'set_fire_params: at &
  &             ', i, j, ' base rate of spread', fp%r_0(i, j), ' moisture ', &
                fp%fmc_g(i, j), '>extinction =', fuelmce(k)
               call message(msg, level=0)
               write (msg, '(5(a,e14.6))') 'rtemp1=', rtemp1, ' etam=', etam
               call message(msg, level=0)
!$OMP END CRITICAL(SFIRE_PHYS_CRIT)
               call warning('propagation above extinction &
  &              moisture', level=0)
            end if
         end do
      end do
      if (nerr .gt. 1) then
!$OMP CRITICAL(SFIRE_PHYS_CRIT)
         write (msg, '(a,i6,a)') 'set_fire_params: WARNING: fuel category 0 &
    &     replaced in', nerr, ' cells'
!$OMP END CRITICAL(SFIRE_PHYS_CRIT)
         call message(msg)
!print*,"estou na rotina consume_crown"
      end if
      !  print*,"fp%bbb_FM10",fp%bbb_FM10(19,19),"fp%phiwc_FM10",fp%phiwc_FM10(19,19),&
      !  "fp%phisc_FM10",fp%phisc_FM10(19,19),"fp%r_0_FM10",fp%r_0_FM10(19,19)
      !  call flush(6)
   end subroutine set_fire_crown_params

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!! INTRODUZIDO POR ISILDA CUNHA MENEZES
   subroutine consume_crown(i, j, rr_FM10, I_initiation, rr_active, fp)
      !subroutine consume_crown(i,j,rr_FM10,fp%r_0_FM10,&
      !I_initiation,rr_active,fp%CFB,fp%HPA)

      implicit none

      real, intent(out)::rr_active
      real, intent(in)::I_initiation, rr_FM10
      integer, intent(in)::i, j
      integer :: k, l
      real :: R_initiation, t_R
      real :: FME
      type(fire_params), intent(inout)::fp
      real, parameter :: FME_0 = 0.0007383

!print*,"estou no modulo_fr_drive_phys e estou na rotina consume_crown"
!call flush(6)
      l = 10
      t_R = 12.595/savr(l)
      k = int(fp%nfuel_cat(i, j))
      FME = ((1.5 - 0.00275*FMC_crown(l))**4.)/(460.+25.9*FMC_crown(l))
!print*,"FME",FME,"FMC_crown",FMC_crown(l)
!call flush(6)
      fp%HPA(i, j) = fp%r_0_FM10(i, j)*t_R
      !print *, "LFR-DBG fp%r_0_FM10(",i,",",j,")=",fp%r_0_FM10(i, j)
      !call flush (6)
      R_initiation = (60*I_initiation)/fp%HPA(i, j)
!print*,"R_initiation",R_initiation,"I_initiation",I_initiation
!call flush(6)
      rr_active = 3.34*(FME/FME_0)*rr_FM10
!print*,"rr_active 0.166666667 =",rr_active
!call flush(6)
      if ((rr_active > R_initiation) .and. (rr_active > 0.166666667)) then
         fp%CFB(i, j) = 1 - exp(-0.23*(rr_active - R_initiation))
!print*,"CFB",fp%CFB(i,j)
!call flush(6)
      end if

!if((i .eq. 19) .and. (j .eq. 19)) then
!print*,"fp%HPA",fp%HPA(19,19),"fp%r_0_FM10",fp%r_0_FM10(19,19),R_initiation,R_initiation,&
! "fp%CFB",fp%CFB(19,19),"rr_active",rr_active,"FME",FME,"rr_FM10",rr_FM10
!call flush(6)
!endif
   end subroutine consume_crown

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!! INTRODUZIDO POR ISILDA CUNHA MENEZES
!!! CALCULA A INTENSIDADE DE BYRAM (INTENSIDADE DE SUPERFICIE) E A INTENSIDADE DE INICIACAO NA COPA
   subroutine fire_suplementary_intensity(dt, fp, i, j, &
                                          ros, Ibyram, I_initiation)

      implicit none
      type(fire_params), intent(in)::fp
      real, intent(in) :: ros
      real, intent(out):: Ibyram
      real, intent(out):: I_initiation

      integer:: i, j, k
      real, intent(in)::dt

      real:: dmass, bmst

!print*,"estou no modulo_fr_drive_phys e estou na rotina fire_suplementary_intensity"
!call flush(6)
!print*,"fgip",fp%fgip(i,j)
!call flush(6)
!print*,"ros",ros
!call flush(6)
!print*,"fmc_g",fp%fmc_g(i,j)
!call flush(6)
!print*,"cmbcn",cmbcnst
!call flush(6)
      dmass = fp%fgip(i, j)*ros
      bmst = fp%fmc_g(i, j)/(1.+fp%fmc_g(i, j))
      Ibyram = (dmass/dt)*(1.-bmst)*cmbcnst

      k = int(fp%nfuel_cat(i, j))
      !  print*,"fp nfuel_cat =",i,j,k
      !  call flush(6)
      !  print*,"ffw",ffw(k)
      !  call flush(6)
      !  print*,"CBH",CBH(k)
      !  call flush(6)
      !  print*,"FMC_crown",FMC_crown(k)
      !  call flush(6)
      Ibyram = Ibyram*ffw(k)
      I_initiation = (CBH(k)*(460.+25.9*FMC_crown(k))/100.)**(3./2.)
      !  print*,"I_initiation",I_initiation
      !  call flush(6)
      !  if((i .eq. 19) .and. (j .eq. 19)) then
      !  print*,"CBH",CBH(k),"FMC_crown",FMC_crown(k)
      !  call flush(6)
      ! endif

   end subroutine fire_suplementary_intensity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine heat_fluxes(dt, fp, &
                          ifms, ifme, jfms, jfme, &
                          ifts, ifte, jfts, jfte, &
                          iffs, iffe, jffs, jffe, &
                          fgip, fuel_frac_burnt, &
                          grnhft, grnqft)
      implicit none

      type(fire_params), intent(in)::fp
      real, intent(in)::dt
      integer, intent(in)::ifts, ifte, jfts, jfte, ifms, ifme, jfms, jfme, iffs, iffe, jffs, jffe
      real, intent(in), dimension(ifms:ifme, jfms:jfme):: fgip
      real, intent(in), dimension(iffs:iffe, jffs:jffe):: fuel_frac_burnt
      real, intent(out), dimension(ifms:ifme, jfms:jfme):: grnhft
      real, intent(out), dimension(ifms:ifme, jfms:jfme), optional:: grnqft

      integer::i, j
      real:: dmass, bmst
      logical::latent

      !print *, "estou no modulo_fr_drive_phys e estou na rotina heat_fluxes"
      !call flush (6)
      latent = present(grnqft)
      do j = jfts, jfte
         do i = ifts, ifte
            dmass = &
               fgip(i, j) &
               *fuel_frac_burnt(i, j)
            bmst = fp%fmc_g(i, j)/(1.+fp%fmc_g(i, j))
            grnhft(i, j) = (dmass/dt)*(1.-bmst)*cmbcnst
            if (latent) grnqft(i, j) = (bmst + (1.-bmst)*.56)*(dmass/dt)*xlv
            !if(grnhft(i,j)>0) write(90,*) i,j,grnhft(i,j), dmass,dt,bmst,cmbcnst,fp%fmc_g(i, j),fgip(i, j),fuel_frac_burnt(i, j)!LFR-DBG
         end do
      end do

   end subroutine heat_fluxes

   subroutine set_nfuel_cat( &
      ifms, ifme, jfms, jfme, &
      ifts, ifte, jfts, jfte, &
      ifuelread, nfuel_cat0, zsf, nfuel_cat)

      implicit none

      integer, intent(in)::   ifts, ifte, jfts, jfte, &
                            ifms, ifme, jfms, jfme

      integer, intent(in)::ifuelread, nfuel_cat0
      real, intent(in), dimension(ifms:ifme, jfms:jfme)::zsf
      real, intent(out), dimension(ifms:ifme, jfms:jfme)::nfuel_cat

      integer:: i, j, iu1
      real:: t1
      character(len=128) msg

      !print *, "estou no modulo_fr_drive_phys e estou na rotina set_nfuel_cat"
      !call flush (6)
!$OMP CRITICAL(SFIRE_PHYS_CRIT)
      write (msg, '(a,i3)') 'set_nfuel_cat: ifuelread=', ifuelread
!$OMP END CRITICAL(SFIRE_PHYS_CRIT)
      call message(msg)

      if (ifuelread .eq. -1 .or. ifuelread .eq. 2) then
!$OMP CRITICAL(SFIRE_PHYS_CRIT)
         call message('set_nfuel_cat: assuming nfuel_cat initialized already')
         call message(msg)
!$OMP END CRITICAL(SFIRE_PHYS_CRIT)
      else if (ifuelread .eq. 0) then

         do j = jfts, jfte
            do i = ifts, ifte
               nfuel_cat(i, j) = real(nfuel_cat0)
            end do
         end do
!$OMP CRITICAL(SFIRE_PHYS_CRIT)
         write (msg, '(a,i3)') 'set_nfuel_cat: fuel initialized with category', nfuel_cat0
!$OMP END CRITICAL(SFIRE_PHYS_CRIT)
         call message(msg)

      else if (ifuelread .eq. 1) then

         do j = jfts, jfte
            do i = ifts, ifte

               t1 = zsf(i, j)
               if (t1 .le. 1524.) then
                  nfuel_cat(i, j) = 3
               else if (t1 .ge. 1524. .and. t1 .le. 2073.) then
                  nfuel_cat(i, j) = 2
               else if (t1 .ge. 2073. .and. t1 .le. 2438.) then
                  nfuel_cat(i, j) = 8
               else if (t1 .gt. 2438. .and. t1 .le. 3354.) then

                  nfuel_cat(i, j) = 10
               else if (t1 .gt. 3354. .and. t1 .le. 3658.) then
                  nfuel_cat(i, j) = 1
               else if (t1 .gt. 3658.) then
                  nfuel_cat(i, j) = 14
               end if
            end do
         end do

         call message('set_nfuel_cat: fuel initialized by altitude')
      else

         call crash('set_nfuel_cat: bad ifuelread')
      end if

   end subroutine set_nfuel_cat

   real function fire_rate_of_spread(ifun, propx, propy, i, j, fp, config_flags)
!config_flags e ifun introduzido por ISILDA CUNHA MENEZES

      implicit none

      real, intent(in)::propx, propy
      integer, intent(in)::i, j
!COMENTADO POR ISILDA
!type(fire_params),intent(in)::fp
!!!!! INTRODUZIDO POR ISILDA CUNHA MENEZES
      type(fire_params), intent(inout)::fp
      type(grid_config_rec_type), intent(IN)  :: config_flags
      integer, intent(in)::ifun
!!!!!!

      real:: ros_back, ros_wind, ros_slope, nvx, nvy, scale, rr
      real:: ros_total !!!!! INTRODUZIDO POR ISILDA CUNHA MENEZES
      integer::ierrx
      character(len=128)::msg

!print*,"estou no modulo_fr_drive_phys e estou na rotina fire_rate_of_spread"
!call flush(6)
      scale = sqrt(propx*propx + propy*propy)
      if (.not. scale > 0.) scale = 1.
      nvx = propx/scale
      nvy = propy/scale

!!!! INTRODUZIDO POR ISILDA CUNHA MENEZES
      call fire_ros(ifun, ros_back, ros_wind, ros_slope, nvx, nvy, i, j, fp, &
                    ierrx, msg, config_flags, ros_total) !config_flags e ifun  introduzido por ISILDA CUNHA MENEZES
      rr = ros_total
!print*,"estou no modulo_fr_drive_phys e estou na rotina fire_rate_of_spread rr=",rr
!!!!!

!call fire_ros(ros_back,ros_wind,ros_slope, nvx,nvy,i,j,fp,ierrx,msg)
!rr = ros_back+ros_wind+ros_slope
      if (fire_grows_only .gt. 0) rr = max(rr, 0.)
!print*,"fire_rate_of_spread =", rr
      fire_rate_of_spread = rr
!print*,"vou sair da funcao fire_rate_of_spread"
   end function fire_rate_of_spread

!!!!! COMENTADO POR ISILDA CUNHA MENEZES
!subroutine fire_ros(ros_back,ros_wind,ros_slope, &
!propx,propy,i,j,fp,ierrx,msg)

!implicit none

!real, intent(out)::ros_back,ros_wind,ros_slope
!real, intent(in)::propx,propy
!integer, intent(in)::i,j
!type(fire_params),intent(in)::fp
!integer, intent(out)::ierrx
!character(len=*), intent(out)::msg

!real:: speed, tanphi
!real::cor_wind,cor_slope,nvx,nvy,scale

!scale=1.
!nvx=propx/scale
!nvy=propy/scale
!if (fire_advection.ne.0) then

!    speed =  sqrt(fp%vx(i,j)*fp%vx(i,j)+ fp%vy(i,j)*fp%vy(i,j))+tiny(speed)

!    tanphi = sqrt(fp%dzdxf(i,j)*fp%dzdxf(i,j) + fp%dzdyf(i,j)*fp%dzdyf(i,j))+tiny(tanphi)

!    cor_wind =  max(0.,(fp%vx(i,j)*nvx + fp%vy(i,j)*nvy)/speed)

!    cor_slope = max(0., (fp%dzdxf(i,j)*nvx + fp%dzdyf(i,j)*nvy)/tanphi)
!else

!    speed =  fp%vx(i,j)*nvx + fp%vy(i,j)*nvy

!    tanphi = fp%dzdxf(i,j)*nvx + fp%dzdyf(i,j)*nvy
!    cor_wind=1.
!    cor_slope=1.
!endif

!call fire_ros_cawfe(ros_back,ros_wind,ros_slope, &
!speed,tanphi,cor_wind,cor_slope,i,j,fp,ierrx,msg)

!end subroutine fire_ros
!!!!!!!!

!!!!!#######################################################################
!!! ROTINA ALTERADA PARA O MODELO CROWN E INTRODUZIDO POR ISILDA CUNHA MENEZES
!########## ROTINA QUE CALCULA A VELOCIDADE DE PROPAGACAO DE ROTHERMEL_FOGO_SUPERFI�CIE (rr_surface) E DE FOGO TOTAL (rr_TOTAL) - Estou a manter o mesmo nome da originaL
   subroutine fire_ros(ifun, ros_back, ros_wind, ros_slope, &
                       propx, propy, i, j, fp, ierrx, msg, config_flags, total_rr) !config_flags e ifun introduzido por ISILDA CM

      implicit none

      real, intent(out):: ros_back, ros_wind, ros_slope
      real, intent(out), optional:: total_rr
      real, intent(in)::propx, propy
      integer, intent(in)::i, j
      type(fire_params), intent(inout)::fp
      integer, intent(out)::ierrx
      character(len=*), intent(out)::msg

      real:: speed, tanphi
      real::cor_wind, cor_slope, nvx, nvy, scale
      real:: rr_FM10, FM10_ros_back, FM10_ros_wind, FM10_ros_slope
      real:: rr_active, surface_rr
      real:: Ibyram, I_initiation

!!## INTRODUZIDO POR ISILDA CUNHA MENEZES
      type(grid_config_rec_type), intent(IN)  :: config_flags
      logical:: crown
      integer, intent(in)::ifun

!print*,"estou no modulo_fr_drive_phys e estou na rotina fire_ros"
!call flush(6)
      crown = config_flags%crown

!!!!!####

      scale = 1.
      nvx = propx/scale
      nvy = propy/scale
      if (fire_advection .ne. 0) then

         speed = sqrt(fp%vx(i, j)*fp%vx(i, j) + fp%vy(i, j)*fp%vy(i, j)) + tiny(speed)

         tanphi = sqrt(fp%dzdxf(i, j)*fp%dzdxf(i, j) + fp%dzdyf(i, j)*fp%dzdyf(i, j)) + tiny(tanphi)

         cor_wind = max(0., (fp%vx(i, j)*nvx + fp%vy(i, j)*nvy)/speed)

         cor_slope = max(0., (fp%dzdxf(i, j)*nvx + fp%dzdyf(i, j)*nvy)/tanphi)
      else

         speed = fp%vx(i, j)*nvx + fp%vy(i, j)*nvy

         tanphi = fp%dzdxf(i, j)*nvx + fp%dzdyf(i, j)*nvy
         cor_wind = 1.
         cor_slope = 1.
      end if

!####### chamada para rr_suface

      call fire_ros_cawfe(ros_back, ros_wind, ros_slope, &
                          speed, tanphi, cor_wind, cor_slope, i, j, fp, ierrx, msg)

!####### chamada para rr_total

      surface_rr = ros_back + ros_wind + ros_slope
!print*,"velocidade surface =",crown,surface_rr
!call flush(6)
      total_rr = surface_rr
      if (ifun .ge. 2) then
      if (crown) then
      if (int(fp%nfuel_cat(i, j)) .eq. 14) then
         total_rr = surface_rr
      else
!print*,"estou no modulo_fr_drive_phys e estou na rotina fire_ros - vou entrar fire_suplementary_intensity"
!call flush(6)
         call fire_suplementary_intensity(1., fp, i, j, &
                                          surface_rr, Ibyram, I_initiation)

         ! print*,"estou no modulo_fr_drive_phys e estou na rotina fire_ros - vou comparar Ibyram com I_initiation"
         ! print*,"Ibyram .gt. I_initiation",Ibyram, I_initiation
         ! call flush(6)
         if (Ibyram .gt. I_initiation) then !(fogo de copa passivo e activo)

            !  print*,"estou no modulo_fr_drive_phys e estou na rotina fire_ros - vou entrar fire_ros_cawfe_crown"
            !  call flush(6)
            call fire_ros_cawfe_crown(FM10_ros_back, FM10_ros_wind, FM10_ros_slope, &
                                      speed, tanphi, cor_wind, cor_slope, i, j, fp, ierrx, msg)

            rr_FM10 = FM10_ros_back + FM10_ros_wind + FM10_ros_slope
            !print*,"estou na rotina fire_ros - vou entrar em consume_crown"
            call consume_crown(i, j, rr_FM10, I_initiation, rr_active, fp)

            total_rr = surface_rr + fp%CFB(i, j)*(rr_active - surface_rr)
            !    print*,"velocidade crown =",total_rr

         end if
      end if
      end if
      end if
      ! print*, "estou no modulo_fr_drive_phys e estou na rotina fire_ros - vou sair"
      ! call flush(6)

   end subroutine fire_ros

!!!!!#######################################################################

!!! ROTINA ALTERADA PARA O MODELO CROWN E INTRODUZIDA POR ISILDA CUNHA MENEZES
!!!!!!####################### calculo de R_FM10    - VELOCIDADE PROPAGACAO de Rothermel PARA FM10

   subroutine fire_ros_cawfe_crown(ros_back_FM10, ros_wind_FM10, ros_slope_FM10, &
                                   speed, tanphi, cor_wind, cor_slope, i, j, fp, ierrx, msg)

      implicit none

      real, intent(out)::ros_back_FM10, ros_wind_FM10, ros_slope_FM10
      real, intent(in)::speed, tanphi, cor_wind, cor_slope
      integer, intent(in)::i, j
      type(fire_params), intent(in)::fp
      integer, intent(out)::ierrx
      character(len=*), intent(out)::msg

      real:: umid, spdms, umidm, excess, tanphim
      real, parameter::ros_max_FM10 = 6.
      integer::k
      real:: ros_FM10, excess_FM10

!print*,"estou no modulo_fr_drive_phys e estou na rotina fire_ros_cawfe_crown"
!call flush(6)

      if (ibeh .eq. 1) then

         spdms = max(speed, 0.)
         umidm = min(spdms, 30.)
         umid = umidm*196.850
         !   print*,"phiwc_FM10",fp%phiwc_FM10(i,j)
         !   call flush(6)
         !   print*,"bbb_FM10",fp%bbb_FM10(i,j)
         !   call flush(6)
         !   print*,"umid",umid
         !   call flush(6)
         ros_wind_FM10 = fp%phiwc_FM10(i, j)*(umid**fp%bbb_FM10(i, j))
         !   print*,"ros_wind_FM10",ros_wind_FM10
         !   call flush(6)
         !   print*,"tanphi",tanphi
         !   call flush(6)
         tanphim = max(tanphi, 0.0)
         tanphim = min(tanphim, 5.0)

         ros_slope_FM10 = fp%phisc_FM10(i, j)*tanphim**2
         !   print*,"r_0_FM10",fp%r_0_FM10(i,j)
         !   call flush(6)
         ros_back_FM10 = fp%r_0_FM10(i, j)

         !   if((i .eq. 19) .and. (j .eq. 19)) then
         !   print*,"fp%phisc_FM10",fp%phisc_FM10(19,19),"fp%r_0_FM10",fp%r_0_FM10(19,19),"ros_slope_FM10",ros_slope_FM10,"ros_back_FM10",ros_back_FM10
         !  call flush(6)
         !  endif

      end if
      k = int(fp%nfuel_cat(i, j))
      ros_FM10 = ros_back_FM10 + ros_wind_FM10 + ros_slope_FM10
      if (ros_FM10 > 1e-6 .and. fp%fmc_g(i, j) > fuelmce(k)) then
!$OMP CRITICAL(SFIRE_PHYS_CRIT)
         write (msg, '(a,2i6,3(a,e13.5))') 'fire_ros_cawfe: at ', i, j, ' rate of spread', ros_FM10, ' moisture ' &
            , fp%fmc_g(i, j), '> extinction =', fuelmce(k)
!$OMP END CRITICAL(SFIRE_PHYS_CRIT)

         ierrx = 1
      end if

      ros_wind_FM10 = ros_wind_FM10*cor_wind
      ros_slope_FM10 = ros_slope_FM10*cor_slope
!print*,"ros_wind_FM10",ros_wind_FM10,"ros_slope_FM10",ros_slope_FM10,"ros_back_FM10",ros_back_FM10,"ros_max_FM10",ros_max_FM10
!call flush(6)
      excess_FM10 = ros_back_FM10 + ros_wind_FM10 + ros_slope_FM10 - ros_max_FM10
!print*,"excess_FM10",excess_FM10
!call flush(6)
      if (excess_FM10 > 0.) then

         ros_wind_FM10 = ros_wind_FM10 - excess_FM10*ros_wind_FM10/(ros_wind_FM10 + ros_slope_FM10)
         ros_slope_FM10 = ros_slope_FM10 - excess_FM10*ros_slope_FM10/(ros_wind_FM10 + ros_slope_FM10)
      end if
!print*,"vou sair da rotina fire_ros_cawfe_crown"
!call flush(6)
      return

   end subroutine fire_ros_cawfe_crown

!!!!!#######################################################################

!!!

   subroutine fire_ros_cawfe(ros_back, ros_wind, ros_slope, &
                             speed, tanphi, cor_wind, cor_slope, i, j, fp, ierrx, msg)

      implicit none

      real, intent(out)::ros_back, ros_wind, ros_slope
      real, intent(in)::speed, tanphi, cor_wind, cor_slope
      integer, intent(in)::i, j
      type(fire_params), intent(in)::fp
      integer, intent(out)::ierrx
      character(len=*), intent(out)::msg

      real:: umid, phis, phiw, spdms, umidm, excess, tanphim, ros
      real, parameter::ros_max = 6.
      integer::k

!print*,"estou no modulo_fr_drive_phys e entrei na rotina fire_ros_cawfe"
!call flush(6)
      ierrx = 0

      if (.not. fp%ischap(i, j) > 0.) then
         if (ibeh .eq. 1) then

            spdms = max(speed, 0.)
            umidm = min(spdms, 30.)
            umid = umidm*196.850

            ros_wind = fp%phiwc(i, j)*(umid**fp%bbb(i, j))
            tanphim = max(tanphi, 0.0)
            tanphim = min(tanphim, 5.0)

            ros_slope = fp%phisc(i, j)*tanphim**2

            ros_back = fp%r_0(i, j)
         elseif (ibeh .eq. 2) then
            ros_back = 0.
            ros_wind = max(speed, 0.)
            ros_slope = 0.
         elseif (ibeh .eq. 3) then
            ros_back = 0.
            ros_wind = speed
            ros_slope = 0.
         elseif (ibeh .eq. 0) then

            ros_back = 0.18
            ros_wind = 0.18*exp(0.8424*max(speed, 0.)) - ros_back
            ros_slope = 0.

         else

            ros_back = -999.
            ros_wind = -999.
            ros_slope = -999.
         end if
         k = int(fp%nfuel_cat(i, j))
         ros = ros_back + ros_wind + ros_slope
         !   print*,"SURF fire_ros_cawfe: at ",i,j,"rate of spread",ros,"moisture",fp%fmc_g(i,j),"extinction =",fuelmce(k)
         !   call flush(6)
         !   print*,"SURF k=",k,"fp%phiwc",fp%phiwc(i,j),"umid",umid,"fp%bbb",fp%bbb(i,j), "tanphim",tanphim,"fp%r_0",fp%r_0(i,j),"speed",speed
         !   call flush(6)
         if (ros > 1e-6 .and. fp%fmc_g(i, j) > fuelmce(k)) then
!$OMP CRITICAL(SFIRE_PHYS_CRIT)
            write (msg, '(a,2i6,3(a,e13.5))') 'fire_ros_cawfe: at ', i, j, ' rate of spread', ros, ' moisture ', fp%fmc_g(i, j) &
               , '> extinction =', fuelmce(k)
!$OMP END CRITICAL(SFIRE_PHYS_CRIT)
            ierrx = 1
         end if

      else

         spdms = max(speed, 0.)

         ros_back = .03333
         ros_wind = 1.2974*spdms**1.41
         ros_wind = max(ros_wind, ros_back) - ros_back
         ros_slope = 0.

      end if

      ros_wind = ros_wind*cor_wind
      ros_slope = ros_slope*cor_slope

      excess = ros_back + ros_wind + ros_slope - ros_max
!print*,"SURF excess",excess
!call flush(6)
      if (excess > 0.) then

         ros_wind = ros_wind - excess*ros_wind/(ros_wind + ros_slope)
         ros_slope = ros_slope - excess*ros_slope/(ros_wind + ros_slope)
      end if

!print*,"SURF excess*ros_slope/(ros_wind+ros_slope)","ros_slope",ros_slope,"ros_wind",ros_wind,"ros_slope",ros_slope
!call flush(6)

      return

   contains
      real function nrm2(u, v)
         real, intent(in)::u, v
         nrm2 = sqrt(u*u + v*v)
      end function nrm2

   end subroutine fire_ros_cawfe

   subroutine fire_risk(fp, &
                        ifms, ifme, jfms, jfme, &
                        ifts, ifte, jfts, jfte, &
                        nfuel_cat, &
                        f_ros0, f_rosx, f_rosy, f_ros, &
                        f_int, f_lineint, f_lineint2)

      type(fire_params), intent(in)::fp
      integer, intent(in):: &
         ifms, ifme, jfms, jfme, &
         ifts, ifte, jfts, jfte
      real, intent(in), dimension(ifms:ifme, jfms:jfme) :: nfuel_cat
      real, intent(out), dimension(ifms:ifme, jfms:jfme) :: &
         f_ros0, f_rosx, f_rosy, f_ros, &
         f_int, f_lineint, f_lineint2

      integer:: i, j, k, ierrx
      real:: cor_wind = 1., cor_slope = 1., dt_fake = 1.
      real:: ros_back, ros_wind, ros_slope, speed, tanphi, front_speed, ros_x, ros_y
      character(len=128)::msg

!print*,"estou no modulo_fr_drive_phys e entrei na rotina fire_risk"
!call flush(6)
      do j = jfts, jfte
         do i = ifts, ifte

            speed = sqrt(fp%vx(i, j)*fp%vx(i, j) + fp%vy(i, j)*fp%vy(i, j)) + tiny(speed)

            tanphi = sqrt(fp%dzdxf(i, j)*fp%dzdxf(i, j) + fp%dzdyf(i, j)*fp%dzdyf(i, j)) + tiny(tanphi)

            call fire_ros_cawfe(ros_back, ros_wind, ros_slope, &
                                speed, tanphi, cor_wind, cor_slope, i, j, fp, ierrx, msg)

            ros_x = ros_wind*fp%vx(i, j)/speed + ros_slope*fp%dzdxf(i, j)/tanphi
            ros_y = ros_wind*fp%vy(i, j)/speed + ros_slope*fp%dzdyf(i, j)/tanphi

            f_ros0(i, j) = ros_back
            f_rosx(i, j) = ros_x
            f_rosy(i, j) = ros_y

            f_ros(i, j) = ros_back + sqrt(ros_x*ros_x + ros_y*ros_y)

         end do
      end do

      call fire_intensity(fp, &
                          ifms, ifme, jfms, jfme, &
                          ifts, ifte, jfts, jfte, &
                          ifms, ifme, jfms, jfme, &
                          f_ros, nfuel_cat, &
                          f_lineint, f_lineint2, f_int)

   end subroutine fire_risk

   subroutine fire_intensity(fp, &
                             ifms, ifme, jfms, jfme, &
                             ifts, ifte, jfts, jfte, &
                             irms, irme, jrms, jrme, &
                             ros, nfuel_cat, &
                             fibyram, filimit, f_int)

      type(fire_params), intent(in)::fp
      integer, intent(in):: &
         ifms, ifme, jfms, jfme, &
         ifts, ifte, jfts, jfte, &
         irms, irme, jrms, jrme
      real, intent(in), dimension(irms:irme, jrms:jrme) :: ros
      real, intent(in), dimension(ifms:ifme, jfms:jfme) :: nfuel_cat
      real, intent(out), dimension(ifms:ifme, jfms:jfme) :: &
         fibyram, filimit
      real, intent(out), dimension(ifms:ifme, jfms:jfme), optional :: f_int

      integer:: i, j, k
      real, dimension(ifts:ifte, jfts:jfte):: rate_frac
      real:: dt_fake = 1.

!print*,"estou no modulo_fr_drive_phys e entrei na rotina fire_intensity"
!call flush(6)
      call heat_fluxes(dt_fake, fp, &
                       ifms, ifme, jfms, jfme, &
                       ifts, ifte, jfts, jfte, &
                       irms, irme, jrms, jrme, &
                       fp%fgip, ros, &
                       fibyram)

      do j = jfts, jfte
         do i = ifts, ifte
            k = int(nfuel_cat(i, j))
            fibyram(i, j) = fibyram(i, j)*ffw(k)
         end do
      end do

      do j = jfts, jfte
         do i = ifts, ifte
            rate_frac(i, j) = 0.5*ros(i, j)/fp%fuel_time(i, j)
         end do
      end do

      call heat_fluxes(dt_fake, fp, &
                       ifms, ifme, jfms, jfme, &
                       ifts, ifte, jfts, jfte, &
                       ifts, ifte, jfts, jfte, &
                       fp%fgip, rate_frac, &
                       filimit)

      if (present(f_int)) then
         do j = jfts, jfte
            do i = ifts, ifte
               k = int(nfuel_cat(i, j))

               rate_frac(i, j) = ffw(k)/(fp%fuel_time(i, j)*(-log(1.-ffw(k))))
            end do
         end do

         call heat_fluxes(dt_fake, fp, &
                          ifms, ifme, jfms, jfme, &
                          ifts, ifte, jfts, jfte, &
                          ifts, ifte, jfts, jfte, &
                          fp%fgip, rate_frac, &
                          f_int)

      end if

   end subroutine fire_intensity

!!!!!!INTRODUZIDO POR ISILDA CUNHA MENEZES

   subroutine fire_total_intensity(fp, &
                                   ifms, ifme, jfms, jfme, &
                                   ifts, ifte, jfts, jfte, &
                                   irms, irme, jrms, jrme, &
                                   fgip, ros, &
                                   I_final)

      implicit none

      type(fire_params), intent(in)::fp

      integer, intent(in)::ifms, ifme, jfms, jfme, ifts, ifte, jfts, jfte, &
                            irms, irme, jrms, jrme
      real, intent(in), dimension(ifms:ifme, jfms:jfme):: fgip
      real, intent(in), dimension(irms:irme, jrms:jrme)::  ros
!real, dimension(ifms,ifme,jfms,jfme):: grnhft
      real, intent(out), dimension(ifms:ifme, jfms:jfme):: I_final
      integer:: i, j
      real :: dmass, bmst, grnhft

      real, parameter :: dt = 1.0

!print*,"estou no modulo_fr_drive_phys e entrei na rotina fire_total_intensity"
!call flush(6)
      do j = jfts, jfte
         do i = ifts, ifte
            dmass = fgip(i, j)*ros(i, j)
            bmst = fp%fmc_g(i, j)/(1.+fp%fmc_g(i, j))
            grnhft = (dmass/dt)*(1.-bmst)*cmbcnst
            I_final(i, j) = ((fp%HPA(i, j) + (grnhft*fp%CFB(i, j)))*ros(i, j))/60.0
         end do
      end do
!if ((i .eq. 19) .and. (j .eq. 19)) then
!print*,"I_final",I_final(19,19),"fp%HPA",fp%HPA(19,19),"grnhft",grnhft,&
! "fp%CFB",fp%CFB(19,19),"ros",ros(19,19)
!call flush(6)
!endif

   end subroutine fire_total_intensity

!!!!!!!!!!!

end module module_fr_sfire_phys
