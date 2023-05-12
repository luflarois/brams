
module module_fr_sfire_atm

   use module_model_constants, only: cp, xlv, g
   use module_fr_sfire_phys, only: fire_wind_height, have_fuel_cats, fcz0, fcwh, nfuelcats, no_fuel_cat, no_fuel_cat2, windrf &
                                   , itree, fueldepthm, mfuelcats
   use module_fr_sfire_util
!USE module_dm , only : wrf_dm_sum_reals

   use module_fr_sfire_phys, only: mfuelcats, nfuelcats
   use module_state_description, only: num_tracer
   use module_fr_sfire_phys, only: fuel_name
   use module_state_description, only: num_chem

   use module_state_description, only: &
      p_tr17_1, &
      p_tr17_2, &
      p_tr17_3, &
      p_tr17_4, &
      p_tr17_5, &
      p_tr17_6, &
      p_tr17_7, &
      p_tr17_8

   use module_configure, only: &
      p_co, &
      p_ch4, &
      p_h2, &
      p_no, &
      p_no2, &
      p_so2, &
      p_nh3, &
      p_p25, &
      p_p25i, &
      p_p25j, &
      p_oc1, &
      p_oc2, &
      p_bc1, &
      p_bc2, &
      p_ald, &
      p_csl, &
      p_eth, &
      p_hc3, &
      p_hc5, &
      p_hcho, &
      p_iso, &
      p_ket, &
      p_mgly, &
      p_ol2, &
      p_olt, &
      p_oli, &
      p_ora2, &
      p_tol, &
      p_xyl, &
      p_bigalk, &
      p_bigene, &
      p_c10h16, &
      p_c2h4, &
      p_c2h5oh, &
      p_c2h6, &
      p_c3h6, &
      p_c3h8, &
      p_ch3cooh, &
      p_ch3oh, &
      p_cres, &
      p_glyald, &
      p_isopr, &
      p_macr, &
      p_mek, &
      p_mvk, &
      !p_smoke, &
      p_sulf, &
      p_dms, &
      p_msa, &
      p_dust_1, &
      p_dust_2, &
      p_dust_3, &
      p_dust_4, &
      p_dust_5, &
      p_seas_1, &
      p_seas_2, &
      p_seas_3, &
      p_seas_4, &
      p_p10

   implicit none
! COMENTADO POR ISILDA CUNHA MENEZES
!integer, parameter, private:: num_chem=0

   logical, save :: have_wind_log_interpolation = .false.

   real, dimension(mfuelcats), save:: &
      co = 0., &
      ch4 = 0., &
      h2 = 0., &
      no = 0., &
      no2 = 0., &
      so2 = 0., &
      nh3 = 0., &
      oc1 = 0., &
      oc2 = 0., &
      bc1 = 0., &
      bc2 = 0., &
      ald = 0., &
      csl = 0., &
      eth = 0., &
      hc3 = 0., &
      p25 = 0., &
      p25i = 0., &
      p25j = 0., &
      hc5 = 0., &
      hcho = 0., &
      iso = 0., &
      ket = 0., &
      mgly = 0., &
      ol2 = 0., &
      olt = 0., &
      oli = 0., &
      ora2 = 0., &
      tol = 0., &
      xyl = 0., &
      bigalk = 0., &
      bigene = 0., &
      c10h16 = 0., &
      c2h4 = 0., &
      c2h5oh = 0., &
      c2h6 = 0., &
      c3h6 = 0., &
      c3h8 = 0., &
      ch3cooh = 0., &
      ch3oh = 0., &
      cres = 0., &
      glyald = 0., &
      isopr = 0., &
      macr = 0., &
      mek = 0., &
      mvk = 0., &
      smoke = 0., &
      sulf = 0., &
      dms = 0., &
      msa = 0., &
      dust_1 = 0., &
      dust_2 = 0., &
      dust_3 = 0., &
      dust_4 = 0., &
      dust_5 = 0., &
      seas_1 = 0., &
      seas_2 = 0., &
      seas_3 = 0., &
      seas_4 = 0., &
      p10 = 0., &
      tr17_1 = 0., &
      tr17_2 = 0., &
      tr17_3 = 0., &
      tr17_4 = 0., &
      tr17_5 = 0., &
      tr17_6 = 0., &
      tr17_7 = 0., &
      tr17_8 = 0.

   real, parameter:: &
      imw_co = 1./28.010, &
      imw_ch4 = 1./16.04, &
      imw_h2 = 1./2.016, &
      imw_no = 1./30.006, &
      imw_no2 = 1./46.006, &
      imw_so2 = 1./64.066, &
      imw_nh3 = 1./17.031

   real, parameter:: mw_air = 28.97

   integer, parameter:: max_tracer = 10, max_chem = 64

   real, save, dimension(max_chem)::c_chem
   real, save, dimension(mfuelcats)::c_fuel
   real, save, dimension(max_tracer)::c_tracer

   logical, save:: emis_read = .false.
   integer, save:: msglevel = 1, printsums = 0
   integer, parameter:: line = 5

contains

   subroutine read_emissions_table(chem_opt, tracer_opt)
      implicit none
      integer, intent(in)::chem_opt, tracer_opt
      logical, external:: wrf_dm_on_monitor
      external::wrf_dm_bcast_integer, wrf_dm_bcast_real
      integer, dimension(10)::compatible_chem_opt
      integer:: iounit, ierr, i
      character(len=128)::msg
      namelist /emissions/ compatible_chem_opt, printsums, &
         co, &
         ch4, &
         h2, &
         no, &
         no2, &
         so2, &
         nh3, &
         p25, &
         p25i, &
         p25j, &
         oc1, &
         oc2, &
         bc1, &
         bc2, &
         ald, &
         csl, &
         eth, &
         hc3, &
         hc5, &
         hcho, &
         iso, &
         ket, &
         mgly, &
         ol2, &
         olt, &
         oli, &
         ora2, &
         tol, &
         xyl, &
         bigalk, &
         bigene, &
         c10h16, &
         c2h4, &
         c2h5oh, &
         c2h6, &
         c3h6, &
         c3h8, &
         ch3cooh, &
         ch3oh, &
         cres, &
         glyald, &
         isopr, &
         macr, &
         mek, &
         mvk, &
         smoke, &
         sulf, &
         dms, &
         msa, &
         dust_1, &
         dust_2, &
         dust_3, &
         dust_4, &
         dust_5, &
         seas_1, &
         seas_2, &
         seas_3, &
         seas_4, &
         p10, &
         tr17_1, &
         tr17_2, &
         tr17_3, &
         tr17_4, &
         tr17_5, &
         tr17_6, &
         tr17_7, &
         tr17_8

!IF ( wrf_dm_on_monitor() ) THEN

      iounit = open_text_file('namelist.fire_emissions', 'read')
      compatible_chem_opt = 0
      read (iounit, emissions, iostat=ierr)
      if (ierr .ne. 0) call crash('read_emissions_table: error reading namelist emissions in file namelist.fire_emissions')
      close (iounit)
      write (msg, '(a,i3,a)') 'reading emissions table for', nfuelcats, ' fuel categories'
      call message(msg, level=0)
      if (.not. any(compatible_chem_opt .eq. chem_opt)) then
         write (msg, '(a,i4,a)') 'read_emissions_table: chem_opt=', chem_opt &
            , ' not between given compatible_chem_opt in namelist.fire_emissions'
         call message(msg, level=0)
         write (msg, '(a,10i4)') 'compatible_chem_opt=', compatible_chem_opt
         call message(msg, level=0)
         call crash('chem_opt in namelist.input is not consistent with namelist.fire_emissions')
      end if
!ENDIF
!call wrf_dm_bcast_integer(printsums, 1)
!call wrf_dm_bcast_real(co,  nfuelcats)
!call wrf_dm_bcast_real(ch4,  nfuelcats)
!call wrf_dm_bcast_real(h2,  nfuelcats)
!call wrf_dm_bcast_real(no,  nfuelcats)
!call wrf_dm_bcast_real(no2,  nfuelcats)
!call wrf_dm_bcast_real(so2,  nfuelcats)
!call wrf_dm_bcast_real(nh3,  nfuelcats)
!call wrf_dm_bcast_real(p25,  nfuelcats)
!call wrf_dm_bcast_real(p25i,  nfuelcats)
!call wrf_dm_bcast_real(p25j,  nfuelcats)
!call wrf_dm_bcast_real(oc1,  nfuelcats)
!call wrf_dm_bcast_real(oc2,  nfuelcats)
!call wrf_dm_bcast_real(bc1,  nfuelcats)
!call wrf_dm_bcast_real(bc2,  nfuelcats)
!call wrf_dm_bcast_real(ald,  nfuelcats)
!call wrf_dm_bcast_real(csl,  nfuelcats)
!call wrf_dm_bcast_real(eth,  nfuelcats)
!call wrf_dm_bcast_real(hc3,  nfuelcats)
!call wrf_dm_bcast_real(hc5,  nfuelcats)
!call wrf_dm_bcast_real(hcho,  nfuelcats)
!call wrf_dm_bcast_real(iso,  nfuelcats)
!call wrf_dm_bcast_real(ket,  nfuelcats)
!call wrf_dm_bcast_real(mgly,  nfuelcats)
!call wrf_dm_bcast_real(ol2, nfuelcats)
!call wrf_dm_bcast_real(olt, nfuelcats)
!call wrf_dm_bcast_real(oli, nfuelcats)
!call wrf_dm_bcast_real(ora2,nfuelcats)
!call wrf_dm_bcast_real(tol, nfuelcats)
!call wrf_dm_bcast_real(xyl, nfuelcats)
!call wrf_dm_bcast_real(bigalk, nfuelcats)
!call wrf_dm_bcast_real(bigene, nfuelcats)
!call wrf_dm_bcast_real(c10h16, nfuelcats)
!call wrf_dm_bcast_real(c2h4, nfuelcats)
!call wrf_dm_bcast_real(c2h5oh, nfuelcats)
!call wrf_dm_bcast_real(c2h6, nfuelcats)
!call wrf_dm_bcast_real(c3h6, nfuelcats)
!call wrf_dm_bcast_real(c3h8, nfuelcats)
!call wrf_dm_bcast_real(ch3cooh, nfuelcats)
!call wrf_dm_bcast_real(ch3oh, nfuelcats)
!call wrf_dm_bcast_real(cres, nfuelcats)
!call wrf_dm_bcast_real(glyald, nfuelcats)

!call wrf_dm_bcast_real(isopr, nfuelcats)
!call wrf_dm_bcast_real(macr, nfuelcats)
!call wrf_dm_bcast_real(mek, nfuelcats)
!call wrf_dm_bcast_real(mvk, nfuelcats)
!call wrf_dm_bcast_real(smoke, nfuelcats)
!call wrf_dm_bcast_real(sulf, nfuelcats)
!call wrf_dm_bcast_real(dms, nfuelcats)
!call wrf_dm_bcast_real(msa, nfuelcats)
!call wrf_dm_bcast_real(dust_1, nfuelcats)
!call wrf_dm_bcast_real(dust_2, nfuelcats)
!call wrf_dm_bcast_real(dust_3, nfuelcats)
!call wrf_dm_bcast_real(dust_4, nfuelcats)
!call wrf_dm_bcast_real(dust_5, nfuelcats)
!call wrf_dm_bcast_real(seas_1, nfuelcats)
!call wrf_dm_bcast_real(seas_2, nfuelcats)
!call wrf_dm_bcast_real(seas_3, nfuelcats)
!call wrf_dm_bcast_real(seas_4, nfuelcats)
!call wrf_dm_bcast_real(p10, nfuelcats)
!call wrf_dm_bcast_real(tr17_1, nfuelcats)
!call wrf_dm_bcast_real(tr17_2, nfuelcats)
!call wrf_dm_bcast_real(tr17_3, nfuelcats)
!call wrf_dm_bcast_real(tr17_4, nfuelcats)
!call wrf_dm_bcast_real(tr17_5, nfuelcats)
!call wrf_dm_bcast_real(tr17_6, nfuelcats)
!call wrf_dm_bcast_real(tr17_7, nfuelcats)
!call wrf_dm_bcast_real(tr17_8, nfuelcats)

      if (fire_print_msg .ge. msglevel .and. printsums .gt. 0) then

         write (msg, '(3(a,i3,1x))') 'c_chem size', num_chem, 'c_tracer size', num_tracer, 'c_fuel size', nfuelcats
         call message(msg, level=2)
         if (num_chem > max_chem) then
            call crash('read_emissions_table: increase max_chem')
         end if
         c_chem = 0.
         if (num_tracer > max_tracer) then
            call crash('read_emissions_table: increase max_tracer')
         end if
         c_tracer = 0.
         c_fuel = 0.
      end if

      emis_read = .true.

   end subroutine read_emissions_table

   subroutine add_fire_emissions(tracer_opt, dt, dx, dy, &
                                 ifms, ifme, jfms, jfme, &
                                 ifts, ifte, jtfs, jfte, &
                                 ids, ide, kds, kde, jds, jde, &
                                 ims, ime, kms, kme, jms, jme, &
                                 its, ite, kts, kte, jts, jte, &
                                 rho, dz8w, &
                                 fgip, fuel_frac_burnt, nfuel_cat, &
                                 chem, tracer)

!ISILDA retirei desta rotina o chem_opt que estava em primeiro lugados parametros

      implicit none

      integer, intent(in)::tracer_opt !,chem_opt
      real, intent(in):: dt, dx, dy
      integer, intent(in)::its, ite, kts, kte, jts, jte, ims, ime, kms, kme, jms, jme &
                            , ids, ide, kds, kde, jds, jde
      integer, intent(in)::ifts, ifte, jtfs, jfte, ifms, ifme, jfms, jfme
      real, intent(in)::rho(ims:ime, kms:kme, jms:jme), &
                         dz8w(ims:ime, kms:kme, jms:jme)
      real, intent(in), dimension(ifms:ifme, jfms:jfme):: fgip, &
                                                          fuel_frac_burnt, &
                                                          nfuel_cat

      real, intent(inout)::chem(ims:ime, kms:kme, jms:jme, num_chem), &
                            tracer(ims:ime, kms:kme, jms:jme, num_tracer)

      integer:: i, i_f, j, j_f, ir, jr, isz1, isz2, jsz1, jsz2, ioff, joff, ibase, jbase, cat, k1, areaw, m, k, k_p, errors
      real::fuel_burnt, vol, air, conv, avgw, emis
      character(len=128) msg

      real, dimension(mfuelcats)::s_fuel, t_fuel, r_fuel
      integer, parameter:: chem_np = 58 !!!!!!!!!!!!!!!!!!!!!!!!!!! MUDADO POR ISILDA CUNHA MENEZES - ESTAVA ZERO
      real, dimension(num_chem) ::s_chem, t_chem, r_chem
      real, dimension(num_chem) ::a_chem, g_chem
      integer:: chem_pointers(chem_np)
      character(len=8)::chem_names(chem_np)
      integer, parameter:: tracer_np_max = 8
      integer:: tracer_pointers(tracer_np_max)

      real, dimension(num_tracer) ::s_tracer, t_tracer, r_tracer
      character(len=8)::tracer_names(tracer_np_max)
      integer::tracer_np

      call check_mesh_2dim(its, ite, jts, jte, ims, ime, jms, jme)
      call check_mesh_2dim(ifts, ifte, jtfs, jfte, ifms, ifme, jfms, jfme)

      if (.not. emis_read) call crash('add_fire_emissions: read_emissions_table must be called first')

!!!!!!! INTRODUZIDO POR ISILDA CUNHA MENEZES

      chem_pointers = (/ &
                      p_co, &
                      p_ch4, &
                      p_h2, &
                      p_no, &
                      p_no2, &
                      p_so2, &
                      p_nh3, &
                      p_p25, &
                      p_p25i, &
                      p_p25j, &
                      p_oc1, &
                      p_oc2, &
                      p_bc1, &
                      p_bc2, &
                      p_ald, &
                      p_csl, &
                      p_eth, &
                      p_hc3, &
                      p_hc5, &
                      p_hcho, &
                      p_iso, &
                      p_ket, &
                      p_mgly, &
                      p_ol2, &
                      p_olt, &
                      p_oli, &
                      p_ora2, &
                      p_tol, &
                      p_xyl, &
                      p_bigalk, &
                      p_bigene, &
                      p_c10h16, &
                      p_c2h4, &
                      p_c2h5oh, &
                      p_c2h6, &
                      p_c3h6, &
                      p_c3h8, &
                      p_ch3cooh, &
                      p_ch3oh, &
                      p_cres, &
                      p_glyald, &
                      p_isopr, &
                      p_macr, &
                      p_mek, &
                      p_mvk, &
                      p_sulf, &
                      p_dms, &
                      p_msa, &
                      p_dust_1, &
                      p_dust_2, &
                      p_dust_3, &
                      p_dust_4, &
                      p_dust_5, &
                      p_seas_1, &
                      p_seas_2, &
                      p_seas_3, &
                      p_seas_4, &
                      p_p10/)

      chem_names = (/ &
                   'co      ', &
                   'ch4     ', &
                   'h2      ', &
                   'no      ', &
                   'no2     ', &
                   'so2     ', &
                   'nh3     ', &
                   'p25     ', &
                   'p25i    ', &
                   'p25j    ', &
                   'oc1     ', &
                   'oc2     ', &
                   'bc1     ', &
                   'bc2     ', &
                   'ald     ', &
                   'csl     ', &
                   'eth     ', &
                   'hc3     ', &
                   'hc5     ', &
                   'hcho    ', &
                   'iso     ', &
                   'ket     ', &
                   'mgly    ', &
                   'ol2     ', &
                   'olt     ', &
                   'oli     ', &
                   'ora2    ', &
                   'tol     ', &
                   'xyl     ', &
                   'bigalk  ', &
                   'bigene  ', &
                   'c10h16  ', &
                   'c2h4    ', &
                   'c2h5oh  ', &
                   'c2h6    ', &
                   'c3h6    ', &
                   'c3h8    ', &
                   'ch3cooh ', &
                   'ch3oh   ', &
                   'cres    ', &
                   'glyald  ', &
                   'isopr   ', &
                   'macr    ', &
                   'mek     ', &
                   'mvk     ', &
                   'sulf    ', &
                   'dms     ', &
                   'msa     ', &
                   'dust_1  ', &
                   'dust_2  ', &
                   'dust_3  ', &
                   'dust_4  ', &
                   'dust_5  ', &
                   'seas_1  ', &
                   'seas_2  ', &
                   'seas_3  ', &
                   'seas_4  ', &
                   'p10     '/)

!!!!!!!!!!!!!!!!!!!

      if (tracer_opt == 0) then
         call message('add_fire_emissions: no tracers')
         tracer_np = 0
      elseif (tracer_opt == 2) then
         tracer_pointers(1:8) = (/ &
                                p_tr17_1, &
                                p_tr17_2, &
                                p_tr17_3, &
                                p_tr17_4, &
                                p_tr17_5, &
                                p_tr17_6, &
                                p_tr17_7, &
                                p_tr17_8/)
         tracer_names(1:8) = (/ &
                             'tr17_1  ', &
                             'tr17_2  ', &
                             'tr17_3  ', &
                             'tr17_4  ', &
                             'tr17_5  ', &
                             'tr17_6  ', &
                             'tr17_7  ', &
                             'tr17_8  '/)
         tracer_np = 8
      else
         call crash('add_fire_emissions: tracer_opt value not supported')
      end if

      if (fire_print_msg .ge. msglevel) then
!$OMP CRITICAL(SFIRE_ATM_CRIT)
900      format(6(a, "=", i3, 1x))
         !   write(msg,900)'add_fire_emissions: chem_opt',chem_opt,'tracer_opt',tracer_opt, &
         !      'num_chem',num_chem,'num_tracer',num_tracer, &
         !       'max tracer',maxval(tracer_pointers)
         call message(msg, level=msglevel)
         call print_dims_4d('tracer', tracer)
!    call print_dims_4d('chem',chem)
!$OMP END CRITICAL(SFIRE_ATM_CRIT)
      end if

!!!!!! INTRODUZIDO POR ISILDA CUNHA MENEZES

      call check_pointers('chem', chem, chem_np, chem_names, chem_pointers)
      if (tracer_opt > 0) then
         call check_pointers('tracer', tracer, tracer_np, tracer_names, tracer_pointers)
      end if

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

      avgw = 1.0/(ir*jr)

      t_chem = 0.
      t_fuel = 0.
      do i = 1, num_tracer
         t_tracer(i) = 0.
      end do

      k1 = kts
!$OMP CRITICAL(SFIRE_ATM_CRIT)
      write (msg, '(a,i3)') 'Fire tracers into atmosphere level', k1
!$OMP END CRITICAL(SFIRE_ATM_CRIT)
      call message(msg, level=msglevel)

      if (fire_print_msg .ge. msglevel .and. printsums .gt. 0) then

         a_chem = 0.0
         g_chem = 0.0
         errors = 0
         do j = jts, jte
            do k = 1, chem_np
               k_p = chem_pointers(k)
               do i = its, ite
                  if (chem(i, k1, j, k_p) .ne. chem(i, k1, j, k_p)) errors = errors + 1
                  a_chem(k_p) = a_chem(k_p) + chem(i, k1, j, k_p)
               end do
            end do
         end do
         if (errors > 0) call crash('NaN before chem update')
         !   call wrf_dm_sum_reals(a_chem,g_chem)
         call message('Layer1 raw sums before adding fire emissions', level=msglevel)
         do i = 1, chem_np, line
            m = min(i + line - 1, chem_np)
            write (msg, 80) 'Emissions ', (trim(chem_names(j)), j=i, m)
            call message(msg, level=msglevel)

                !!!!!! INTRODUZIDO POR ISILDA APAGAR DEPOIS

            write (msg, 81) 'Layer1 end', (a_chem(chem_pointers(j)), j=i, m)

                !!!!!!!
            ! write(msg,81)'Layer1 beg',(g_chem(chem_pointers(j)),j=i,m)
            call message(msg, level=msglevel)
            call message(' ', level=msglevel)
         end do
      end if

!!!!!!!!!!!!!!!!!!!

      call message("Inserting chemical species", level=msglevel)
      do j = max(jds + 1, jts), min(jte, jde - 1)
         jbase = jtfs + jr*(j - jts)
         do i = max(ids + 1, its), min(ite, ide - 1)
            ibase = ifts + ir*(i - its)

            do joff = 0, jr - 1
               j_f = joff + jbase
               do ioff = 0, ir - 1
                  i_f = ioff + ibase

                  fuel_burnt = fgip(i_f, j_f)*fuel_frac_burnt(i_f, j_f)
                  cat = nfuel_cat(i_f, j_f)

                  if (cat .lt. no_fuel_cat) t_fuel(cat) = t_fuel(cat) + fuel_burnt

!!!!! INTRODUZIDO POR ISILDA CUNHA MENEZES

                  conv = avgw*1e3*mw_air/(rho(i, k1, j)*dz8w(i, k1, j))

                  emis = co(cat)*fuel_burnt
                  t_chem(p_co) = t_chem(p_co) + emis

                  chem(i, k1, j, p_co) = chem(i, k1, j, p_co) + emis*conv*imw_co

                  emis = ch4(cat)*fuel_burnt
                  t_chem(p_ch4) = t_chem(p_ch4) + emis
                  chem(i, k1, j, p_ch4) = chem(i, k1, j, p_ch4) + emis*conv*imw_ch4

                  emis = h2(cat)*fuel_burnt
                  t_chem(p_h2) = t_chem(p_h2) + emis
                  chem(i, k1, j, p_h2) = chem(i, k1, j, p_h2) + emis*conv*imw_h2

                  emis = no(cat)*fuel_burnt
                  t_chem(p_no) = t_chem(p_no) + emis
                  chem(i, k1, j, p_no) = chem(i, k1, j, p_no) + emis*conv*imw_no

                  emis = no2(cat)*fuel_burnt
                  t_chem(p_no2) = t_chem(p_no2) + emis
                  chem(i, k1, j, p_no2) = chem(i, k1, j, p_no2) + emis*conv*imw_no2

                  emis = so2(cat)*fuel_burnt
                  t_chem(p_so2) = t_chem(p_so2) + emis
                  chem(i, k1, j, p_so2) = chem(i, k1, j, p_so2) + emis*conv*imw_so2

                  emis = nh3(cat)*fuel_burnt
                  t_chem(p_nh3) = t_chem(p_nh3) + emis
                  chem(i, k1, j, p_nh3) = chem(i, k1, j, p_nh3) + emis*conv*imw_nh3

                  emis = ald(cat)*fuel_burnt
                  t_chem(p_ald) = t_chem(p_ald) + emis
                  chem(i, k1, j, p_ald) = chem(i, k1, j, p_ald) + emis*conv

                  emis = csl(cat)*fuel_burnt
                  t_chem(p_csl) = t_chem(p_csl) + emis
                  chem(i, k1, j, p_csl) = chem(i, k1, j, p_csl) + emis*conv

                  emis = eth(cat)*fuel_burnt
                  t_chem(p_eth) = t_chem(p_eth) + emis
                  chem(i, k1, j, p_eth) = chem(i, k1, j, p_eth) + emis*conv

                  emis = hc3(cat)*fuel_burnt
                  t_chem(p_hc3) = t_chem(p_hc3) + emis
                  chem(i, k1, j, p_hc3) = chem(i, k1, j, p_hc3) + emis*conv

                  emis = hc5(cat)*fuel_burnt
                  t_chem(p_hc5) = t_chem(p_hc5) + emis
                  chem(i, k1, j, p_hc5) = chem(i, k1, j, p_hc5) + emis*conv

                  emis = hcho(cat)*fuel_burnt
                  t_chem(p_hcho) = t_chem(p_hcho) + emis
                  chem(i, k1, j, p_hcho) = chem(i, k1, j, p_hcho) + emis*conv

                  emis = iso(cat)*fuel_burnt
                  t_chem(p_iso) = t_chem(p_iso) + emis
                  chem(i, k1, j, p_iso) = chem(i, k1, j, p_iso) + emis*conv

                  emis = ket(cat)*fuel_burnt
                  t_chem(p_ket) = t_chem(p_ket) + emis
                  chem(i, k1, j, p_ket) = chem(i, k1, j, p_ket) + emis*conv

                  emis = mgly(cat)*fuel_burnt
                  t_chem(p_mgly) = t_chem(p_mgly) + emis
                  chem(i, k1, j, p_mgly) = chem(i, k1, j, p_mgly) + emis*conv

                  emis = ol2(cat)*fuel_burnt
                  t_chem(p_ol2) = t_chem(p_ol2) + emis
                  chem(i, k1, j, p_ol2) = chem(i, k1, j, p_ol2) + emis*conv

                  emis = olt(cat)*fuel_burnt
                  t_chem(p_olt) = t_chem(p_olt) + emis
                  chem(i, k1, j, p_olt) = chem(i, k1, j, p_olt) + emis*conv

                  emis = oli(cat)*fuel_burnt
                  t_chem(p_oli) = t_chem(p_oli) + emis
                  chem(i, k1, j, p_oli) = chem(i, k1, j, p_oli) + emis*conv

                  emis = ora2(cat)*fuel_burnt
                  t_chem(p_ora2) = t_chem(p_ora2) + emis
                  chem(i, k1, j, p_ora2) = chem(i, k1, j, p_ora2) + emis*conv

                  emis = tol(cat)*fuel_burnt
                  t_chem(p_tol) = t_chem(p_tol) + emis
                  chem(i, k1, j, p_tol) = chem(i, k1, j, p_tol) + emis*conv

                  emis = xyl(cat)*fuel_burnt
                  t_chem(p_xyl) = t_chem(p_xyl) + emis
                  chem(i, k1, j, p_xyl) = chem(i, k1, j, p_xyl) + emis*conv

                  emis = bigalk(cat)*fuel_burnt
                  t_chem(p_bigalk) = t_chem(p_bigalk) + emis
                  chem(i, k1, j, p_bigalk) = chem(i, k1, j, p_bigalk) + emis*conv

                  emis = bigene(cat)*fuel_burnt
                  t_chem(p_bigene) = t_chem(p_bigene) + emis
                  chem(i, k1, j, p_bigene) = chem(i, k1, j, p_bigene) + emis*conv

                  emis = c10h16(cat)*fuel_burnt
                  t_chem(p_c10h16) = t_chem(p_c10h16) + emis
                  chem(i, k1, j, p_c10h16) = chem(i, k1, j, p_c10h16) + emis*conv

                  emis = c2h4(cat)*fuel_burnt
                  t_chem(p_c2h4) = t_chem(p_c2h4) + emis
                  chem(i, k1, j, p_c2h4) = chem(i, k1, j, p_c2h4) + emis*conv

                  emis = c2h5oh(cat)*fuel_burnt
                  t_chem(p_c2h5oh) = t_chem(p_c2h5oh) + emis
                  chem(i, k1, j, p_c2h5oh) = chem(i, k1, j, p_c2h5oh) + emis*conv

                  emis = c2h6(cat)*fuel_burnt
                  t_chem(p_c2h6) = t_chem(p_c2h6) + emis
                  chem(i, k1, j, p_c2h6) = chem(i, k1, j, p_c2h6) + emis*conv

                  emis = c3h6(cat)*fuel_burnt
                  t_chem(p_c3h6) = t_chem(p_c3h6) + emis
                  chem(i, k1, j, p_c3h6) = chem(i, k1, j, p_c3h6) + emis*conv

                  emis = c3h8(cat)*fuel_burnt
                  t_chem(p_c3h8) = t_chem(p_c3h8) + emis
                  chem(i, k1, j, p_c3h8) = chem(i, k1, j, p_c3h8) + emis*conv

                  emis = ch3cooh(cat)*fuel_burnt
                  t_chem(p_ch3cooh) = t_chem(p_ch3cooh) + emis
                  chem(i, k1, j, p_ch3cooh) = chem(i, k1, j, p_ch3cooh) + emis*conv

                  emis = ch3oh(cat)*fuel_burnt
                  t_chem(p_ch3oh) = t_chem(p_ch3oh) + emis
                  chem(i, k1, j, p_ch3oh) = chem(i, k1, j, p_ch3oh) + emis*conv

                  emis = cres(cat)*fuel_burnt
                  t_chem(p_cres) = t_chem(p_cres) + emis
                  chem(i, k1, j, p_cres) = chem(i, k1, j, p_cres) + emis*conv

                  emis = glyald(cat)*fuel_burnt
                  t_chem(p_glyald) = t_chem(p_glyald) + emis
                  chem(i, k1, j, p_glyald) = chem(i, k1, j, p_glyald) + emis*conv

                  emis = isopr(cat)*fuel_burnt
                  t_chem(p_isopr) = t_chem(p_isopr) + emis
                  chem(i, k1, j, p_isopr) = chem(i, k1, j, p_isopr) + emis*conv

                  emis = macr(cat)*fuel_burnt
                  t_chem(p_macr) = t_chem(p_macr) + emis
                  chem(i, k1, j, p_macr) = chem(i, k1, j, p_macr) + emis*conv

                  emis = mek(cat)*fuel_burnt
                  t_chem(p_mek) = t_chem(p_mek) + emis
                  chem(i, k1, j, p_mek) = chem(i, k1, j, p_mek) + emis*conv

                  emis = mvk(cat)*fuel_burnt
                  t_chem(p_mvk) = t_chem(p_mvk) + emis
                  chem(i, k1, j, p_mvk) = chem(i, k1, j, p_mvk) + emis*conv

                  conv = avgw*1e6/(rho(i, k1, j)*dz8w(i, k1, j))

                  emis = p25(cat)*fuel_burnt
                  t_chem(p_p25) = t_chem(p_p25) + emis
                  chem(i, k1, j, p_p25) = chem(i, k1, j, p_p25) + emis*conv

                  emis = p25i(cat)*fuel_burnt
                  t_chem(p_p25i) = t_chem(p_p25i) + emis
                  chem(i, k1, j, p_p25i) = chem(i, k1, j, p_p25i) + emis*conv

                  emis = p25j(cat)*fuel_burnt
                  t_chem(p_p25j) = t_chem(p_p25j) + emis
                  chem(i, k1, j, p_p25j) = chem(i, k1, j, p_p25j) + emis*conv

                  emis = oc1(cat)*fuel_burnt
                  t_chem(p_oc1) = t_chem(p_oc1) + emis
                  chem(i, k1, j, p_oc1) = chem(i, k1, j, p_oc1) + emis*conv

                  emis = oc2(cat)*fuel_burnt
                  t_chem(p_oc2) = t_chem(p_oc2) + emis
                  chem(i, k1, j, p_oc2) = chem(i, k1, j, p_oc2) + emis*conv

                  emis = bc1(cat)*fuel_burnt
                  t_chem(p_bc1) = t_chem(p_bc1) + emis
                  chem(i, k1, j, p_bc1) = chem(i, k1, j, p_bc1) + emis*conv

                  emis = bc2(cat)*fuel_burnt
                  t_chem(p_bc2) = t_chem(p_bc2) + emis
                  chem(i, k1, j, p_bc2) = chem(i, k1, j, p_bc2) + emis*conv

                  emis = sulf(cat)*fuel_burnt
                  t_chem(p_sulf) = t_chem(p_sulf) + emis
                  chem(i, k1, j, p_sulf) = chem(i, k1, j, p_sulf) + emis*conv

                  emis = dms(cat)*fuel_burnt
                  t_chem(p_dms) = t_chem(p_dms) + emis
                  chem(i, k1, j, p_dms) = chem(i, k1, j, p_dms) + emis*conv

                  emis = msa(cat)*fuel_burnt
                  t_chem(p_msa) = t_chem(p_msa) + emis
                  chem(i, k1, j, p_msa) = chem(i, k1, j, p_msa) + emis*conv

                  emis = dust_1(cat)*fuel_burnt
                  t_chem(p_dust_1) = t_chem(p_dust_1) + emis
                  chem(i, k1, j, p_dust_1) = chem(i, k1, j, p_dust_1) + emis*conv

                  emis = dust_2(cat)*fuel_burnt
                  t_chem(p_dust_2) = t_chem(p_dust_2) + emis
                  chem(i, k1, j, p_dust_2) = chem(i, k1, j, p_dust_2) + emis*conv

                  emis = dust_3(cat)*fuel_burnt
                  t_chem(p_dust_3) = t_chem(p_dust_3) + emis
                  chem(i, k1, j, p_dust_3) = chem(i, k1, j, p_dust_3) + emis*conv

                  emis = dust_4(cat)*fuel_burnt
                  t_chem(p_dust_4) = t_chem(p_dust_4) + emis
                  chem(i, k1, j, p_dust_4) = chem(i, k1, j, p_dust_4) + emis*conv

                  emis = dust_5(cat)*fuel_burnt
                  t_chem(p_dust_5) = t_chem(p_dust_5) + emis
                  chem(i, k1, j, p_dust_5) = chem(i, k1, j, p_dust_5) + emis*conv

                  emis = seas_1(cat)*fuel_burnt
                  t_chem(p_seas_1) = t_chem(p_seas_1) + emis
                  chem(i, k1, j, p_seas_1) = chem(i, k1, j, p_seas_1) + emis*conv

                  emis = seas_2(cat)*fuel_burnt
                  t_chem(p_seas_2) = t_chem(p_seas_2) + emis
                  chem(i, k1, j, p_seas_2) = chem(i, k1, j, p_seas_2) + emis*conv

                  emis = seas_3(cat)*fuel_burnt
                  t_chem(p_seas_3) = t_chem(p_seas_3) + emis
                  chem(i, k1, j, p_seas_3) = chem(i, k1, j, p_seas_3) + emis*conv

                  emis = seas_4(cat)*fuel_burnt
                  t_chem(p_seas_4) = t_chem(p_seas_4) + emis
                  chem(i, k1, j, p_seas_4) = chem(i, k1, j, p_seas_4) + emis*conv

                  emis = p10(cat)*fuel_burnt
                  t_chem(p_p10) = t_chem(p_p10) + emis
                  chem(i, k1, j, p_p10) = chem(i, k1, j, p_p10) + emis*conv

!!!!!!!!!!!!!!!!!!!!

                  conv = avgw*1e3*mw_air/(rho(i, k1, j)*dz8w(i, k1, j))
                  if (tracer_opt > 0) then

                     conv = avgw*1e6/(rho(i, k1, j)*dz8w(i, k1, j))
                  end if

                  if (tracer_opt == 2) then
                     emis = tr17_1(cat)*fuel_burnt
                     t_tracer(p_tr17_1) = t_tracer(p_tr17_1) + emis
                     tracer(i, k1, j, p_tr17_1) = tracer(i, k1, j, p_tr17_1) + emis*conv

                     emis = tr17_2(cat)*fuel_burnt
                     t_tracer(p_tr17_2) = t_tracer(p_tr17_2) + emis
                     tracer(i, k1, j, p_tr17_2) = tracer(i, k1, j, p_tr17_2) + emis*conv

                     emis = tr17_3(cat)*fuel_burnt
                     t_tracer(p_tr17_3) = t_tracer(p_tr17_3) + emis
                     tracer(i, k1, j, p_tr17_3) = tracer(i, k1, j, p_tr17_3) + emis*conv

                     emis = tr17_4(cat)*fuel_burnt
                     t_tracer(p_tr17_4) = t_tracer(p_tr17_4) + emis
                     tracer(i, k1, j, p_tr17_4) = tracer(i, k1, j, p_tr17_4) + emis*conv

                     emis = tr17_5(cat)*fuel_burnt
                     t_tracer(p_tr17_5) = t_tracer(p_tr17_5) + emis
                     tracer(i, k1, j, p_tr17_5) = tracer(i, k1, j, p_tr17_5) + emis*conv

                     emis = tr17_6(cat)*fuel_burnt
                     t_tracer(p_tr17_6) = t_tracer(p_tr17_6) + emis
                     tracer(i, k1, j, p_tr17_6) = tracer(i, k1, j, p_tr17_6) + emis*conv

                     emis = tr17_7(cat)*fuel_burnt
                     t_tracer(p_tr17_7) = t_tracer(p_tr17_7) + emis
                     tracer(i, k1, j, p_tr17_7) = tracer(i, k1, j, p_tr17_7) + emis*conv

                     emis = tr17_8(cat)*fuel_burnt
                     t_tracer(p_tr17_8) = t_tracer(p_tr17_8) + emis
                     tracer(i, k1, j, p_tr17_8) = tracer(i, k1, j, p_tr17_8) + emis*conv
                  end if
               end do
            end do
         end do
      end do

!!!!!!!!!!!!!! INTRODUZIDO POR ISILDA CUNHA MENEZES

      if (fire_print_msg .ge. msglevel .and. printsums .gt. 0) then

         a_chem = 0.0
         g_chem = 0.0
         errors = 0
         do j = jts, jte
            do k = 1, chem_np
               k_p = chem_pointers(k)
               do i = its, ite
                  if (chem(i, k1, j, k_p) .ne. chem(i, k1, j, k_p)) errors = errors + 1
                  a_chem(k_p) = a_chem(k_p) + chem(i, k1, j, k_p)
               end do
            end do
         end do
         if (errors > 0) call crash('NaN after chem update')
         !  call wrf_dm_sum_reals(a_chem,g_chem)
         call message('Layer1 raw sums after adding fire emissions', level=msglevel)
         do i = 1, chem_np, line
            m = min(i + line - 1, chem_np)
            write (msg, 80) 'Emissions ', (trim(chem_names(j)), j=i, m)
            call message(msg, level=msglevel)

                  !!!!! INTRODUZIDO POR ISILDA APAGAR DEPOIS

            write (msg, 81) 'Layer1 end', (a_chem(chem_pointers(j)), j=i, m)
         !!!!!!!!
            ! write(msg,81)'Layer1 end',(g_chem(chem_pointers(j)),j=i,m)
            call message(msg, level=msglevel)
            call message(' ', level=msglevel)
         end do
      end if

      if (fire_print_msg .ge. msglevel .and. printsums .gt. 0) then
         call message("Computing totals over all processes", level=msglevel)

         !  call wrf_dm_sum_reals(t_fuel,s_fuel)

  !!!!! INTRODUZIDO POR ISILDA APAGAR DEPOIS

         s_fuel = t_fuel

  !!!!!!!!
         s_fuel = s_fuel*dx*dy

         r_fuel = s_fuel/dt

         if (mfuelcats < nfuelcats) call crash('add_fire_emissions: bad size c_fuel')
         do i = 1, nfuelcats
            c_fuel(i) = c_fuel(i) + s_fuel(i)
         end do

         ! call wrf_dm_sum_reals(a_chem,g_chem)

         ! call wrf_dm_sum_reals(t_chem,s_chem)

   !!!!! INTRODUZIDO POR ISILDA APAGAR DEPOIS
         s_chem = t_chem

        !!!!!!
         s_chem = s_chem*dx*dy
         r_chem = s_chem/dt
         if (max_chem < num_chem) call crash('add_fire_emissions: bad size c_chem')
         do i = 1, num_chem
            c_chem(i) = c_chem(i) + s_chem(i)
         end do

         call message('Total emissions in g or mol per the table', level=msglevel)
         do i = 1, chem_np, line
            m = min(i + line - 1, chem_np)
            write (msg, 80) 'Emissions ', (trim(chem_names(j)), j=i, m)
            call message(msg, level=msglevel)
            write (msg, 81) 'Cumulative', (c_chem(chem_pointers(j)), j=i, m)
            call message(msg, level=msglevel)
            write (msg, 81) 'Rate per s', (r_chem(chem_pointers(j)), j=i, m)
            call message(msg, level=msglevel)
            call message(' ', level=msglevel)
         end do

!!!!!!!!!!!!!

         if (num_tracer > 0) then
            call message("Computing totals of tracers over all processes", level=msglevel)

            !   call wrf_dm_sum_reals(t_tracer,s_tracer)
                !!!!! INTRODUZIDO POR ISILDA APAGAR DEPOIS
            s_tracer = t_tracer

                !!!!!

            s_tracer = s_tracer*dx*dy
            r_tracer = s_tracer/dt
            if (max_tracer < num_tracer) call crash('add_fire_emissions: bad size c_tracer')
            do i = 1, num_tracer
               c_tracer(i) = c_tracer(i) + s_tracer(i)
            end do

            call message('Total emissions in g', level=msglevel)
            do i = 1, tracer_np, line
               m = min(i + line - 1, tracer_np)
               write (msg, 80) 'Emissions ', (trim(tracer_names(j)), j=i, m)
               call message(msg, level=msglevel)
               write (msg, 81) 'Cumulative', (c_tracer(tracer_pointers(j)), j=i, m)
               call message(msg, level=msglevel)
               write (msg, 81) 'Rate per s', (r_tracer(tracer_pointers(j)), j=i, m)
               call message(msg, level=msglevel)
               call message(' ', level=msglevel)
            end do
         end if

         do cat = 1, nfuelcats
            if (c_fuel(cat) > 0.) then
               write (msg, 83) fuel_name(cat), ' burned', c_fuel(cat), 'kg, rate', r_fuel(cat), ' kg/s'
               call message(msg, level=msglevel)
            end if
         end do
         write (msg, 83) 'Total fuel', ' burned', sum(c_fuel), 'kg, rate', sum(r_fuel), ' kg/s'
         call message(msg, level=msglevel)

      end if
80    format(a, 8a11)
81    format(a, 8e11.3)
83    format(a30, a, g14.4, a, g14.4, a, a)

      return

9     continue
!$OMP CRITICAL(SFIRE_ATM_CRIT)
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
!$OMP END CRITICAL(SFIRE_ATM_CRIT)
      call crash('add_fire_emissions: bad mesh sizes')

   end subroutine add_fire_emissions

   subroutine check_pointers(array_name, array, np, pointer_names, pointers)
      implicit none

      character(len=*), intent(in)::array_name
      real, dimension(:, :, :, :), intent(in)::array
      integer, intent(in)::np
      character(len=*), dimension(:), intent(in)::pointer_names
      integer, dimension(:), intent(in)::pointers

      integer::i, m, j, ps, pe
      character(len=256)::msg

!$OMP CRITICAL(SFIRE_ATM_CRIT)
      do i = 1, np, line
         m = min(i + line - 1, np)
         write (msg, '(a7,8(1x,a8))') array_name, (trim(pointer_names(j)), j=i, m)
         call message(msg, level=msglevel)
         write (msg, '(a,8i9)') 'Pointer', (pointers(j), j=i, m)
         call message(msg, level=msglevel)
         call message(' ')
      end do

      ps = lbound(array, 4)
      pe = ubound(array, 4)
      do i = 1, np
         write (msg, 995) array_name, 'species', trim(pointer_names(i)), 'pointer', pointers(i)
995      format(a, 1x, a, 1x, a, 1x, a, i3)
         call message(msg, level=2)
         if (pointers(i) < ps .or. pointers(i) > pe) then
            call print_dims_4d(array_name, array)
            write (msg, 995) array_name, 'species', trim(pointer_names(i)), 'pointer', pointers(i)
            call message(msg, level=-1)
            call crash('add_fire_emissions: pointer out of bounds')
         end if
      end do
!$OMP END CRITICAL(SFIRE_ATM_CRIT)
   end subroutine check_pointers

   subroutine print_dims_4d(array_name, array)

      character(len=*), intent(in) ::array_name
      real, dimension(:, :, :, :), intent(in) ::array

      character(len=256)::msg

993   format(3a, 4(1x, i3, ':', i3))
      write (msg, 993) 'array ', array_name, ' has dimensions ', &
         lbound(array, 1), ubound(array, 1), &
         lbound(array, 2), ubound(array, 2), &
         lbound(array, 3), ubound(array, 3), &
         lbound(array, 4), ubound(array, 4)
      call message(msg, level=-1)
   end subroutine print_dims_4d

   subroutine fire_tendency( &
      ids, ide, kds, kde, jds, jde, &
      ims, ime, kms, kme, jms, jme, &
      its, ite, kts, kte, jts, jte, &
      grnhfx, grnqfx, canhfx, canqfx, &
      alfg, alfc, z1can, &
      zs, z_at_w, dz8w, mu, rho, &
      rthfrten, rqvfrten)

      implicit none

      integer, intent(in) :: ids, ide, kds, kde, jds, jde, &
                             ims, ime, kms, kme, jms, jme, &
                             its, ite, kts, kte, jts, jte

      real, intent(in), dimension(ims:ime, jms:jme) :: grnhfx, grnqfx
      real, intent(in), dimension(ims:ime, jms:jme) :: canhfx, canqfx
      real, intent(in), dimension(ims:ime, jms:jme) :: zs
      real, intent(in), dimension(ims:ime, jms:jme) :: mu

      real, intent(in), dimension(ims:ime, kms:kme, jms:jme) :: z_at_w
      real, intent(in), dimension(ims:ime, kms:kme, jms:jme) :: dz8w
      real, intent(in), dimension(ims:ime, kms:kme, jms:jme) :: rho

      real, intent(in) :: alfg
      real, intent(in) :: alfc
      real, intent(in) :: z1can

      real, intent(out), dimension(ims:ime, kms:kme, jms:jme) :: &
         rthfrten, &
         rqvfrten

      integer :: i, j, k
      integer :: i_st, i_en, j_st, j_en, k_st, k_en

      real :: cp_i
      real :: rho_i
      real :: xlv_i
      real :: z_w
      real :: fact_g, fact_c
      real :: alfg_i, alfc_i

      real, dimension(its:ite, kts:kte, jts:jte) :: hfx, qfx

      do j = jts, jte
         do k = kts, min(kte + 1, kde)
            do i = its, ite
               rthfrten(i, k, j) = 0.
               rqvfrten(i, k, j) = 0.
            end do
         end do
      end do

      cp_i = 1./cp
      xlv_i = 1./xlv
      alfg_i = 1./alfg
      alfc_i = 1./alfc

      call print_2d_stats(its, ite, jts, jte, ims, ime, jms, jme, grnhfx, 'fire_tendency:grnhfx')
      call print_2d_stats(its, ite, jts, jte, ims, ime, jms, jme, grnqfx, 'fire_tendency:grnqfx')

      i_st = max(its, ids + 1)
      i_en = min(ite, ide - 1)
      k_st = kts
      k_en = min(kte, kde - 1)
      j_st = max(jts, jds + 1)
      j_en = min(jte, jde - 1)

      do j = j_st, j_en
         do k = k_st, k_en
            do i = i_st, i_en

               z_w = z_at_w(i, k, j) - zs(i, j)

               fact_g = cp_i*exp(-alfg_i*z_w)
               if (z_w < z1can) then
                  fact_c = cp_i
               else
                  fact_c = cp_i*exp(-alfc_i*(z_w - z1can))
               end if
               hfx(i, k, j) = fact_g*grnhfx(i, j) + fact_c*canhfx(i, j)

               fact_g = xlv_i*exp(-alfg_i*z_w)
               if (z_w < z1can) then
                  fact_c = xlv_i
               else
                  fact_c = xlv_i*exp(-alfc_i*(z_w - z1can))
               end if
               qfx(i, k, j) = fact_g*grnqfx(i, j) + fact_c*canqfx(i, j)

            end do
         end do
      end do

      do j = j_st, j_en
         do k = k_st, k_en - 1
            do i = i_st, i_en

               rho_i = 1./rho(i, k, j)

               rthfrten(i, k, j) = -mu(i, j)*rho_i*(hfx(i, k + 1, j) - hfx(i, k, j))/dz8w(i, k, j)
               rqvfrten(i, k, j) = -mu(i, j)*rho_i*(qfx(i, k + 1, j) - qfx(i, k, j))/dz8w(i, k, j)

            end do
         end do
      end do

      call print_3d_stats(its, ite, kts, kte, jts, jte, ims, ime, kms, kme, jms, jme, rthfrten, 'fire_tendency:rthfrten')
      call print_3d_stats(its, ite, kts, kte, jts, jte, ims, ime, kms, kme, jms, jme, rqvfrten, 'fire_tendency:rqvfrten')

      return

   end subroutine fire_tendency

   subroutine interpolate_atm2fire(id, &
                                   ids, ide, kds, kde, jds, jde, &
                                   ims, ime, kms, kme, jms, jme, &
                                   ips, ipe, jps, jpe, &
                                   its, ite, jts, jte, &
                                   ifds, ifde, jfds, jfde, &
                                   ifms, ifme, jfms, jfme, &
                                   ifps, ifpe, jfps, jfpe, &
                                   ifts, ifte, jfts, jfte, &
                                   ir, jr, &
                                   u_frame, v_frame, &
                                   u, v, &
                                   ph, phb, &
                                   z0, zs, &
                                   uah, vah, &
                                   uf, vf)

      implicit none

      integer, intent(in)::id
      integer, intent(in):: &
         ids, ide, kds, kde, jds, jde, &
         ims, ime, kms, kme, jms, jme, &
         ips, ipe, jps, jpe, &
         its, ite, jts, jte, &
         ifds, ifde, jfds, jfde, &
         ifms, ifme, jfms, jfme, &
         ifps, ifpe, jfps, jfpe, &
         ifts, ifte, jfts, jfte, &
         ir, jr
      real, intent(in):: u_frame, v_frame
      real, intent(in), dimension(ims:ime, kms:kme, jms:jme):: &
         u, v, &
         ph, phb
      real, intent(in), dimension(ims:ime, jms:jme):: &
         z0, &
         zs
      real, intent(out), dimension(ims:ime, jms:jme):: &
         uah, &
         vah
      real, intent(out), dimension(ifms:ifme, jfms:jfme):: &
         uf, vf

      character(len=256)::msg
      real, dimension(its - 2:ite + 2, jts - 2:jte + 2):: ua, va
      real, dimension(its - 2:ite + 2, kds:kde, jts - 2:jte + 2):: altw, altub, altvb, hgtu, hgtv
      integer:: i, j, k, ifts1, ifte1, jfts1, jfte1, ite1, jte1
      integer::itst, itet, jtst, jtet, itsu, iteu, jtsu, jteu, itsv, itev, jtsv, jtev
      integer::kdmax, its1, jts1, ips1, jps1
      integer::itsou, iteou, jtsou, jteou, itsov, iteov, jtsov, jteov
      real:: ground, loght, loglast, logz0, logfwh, ht, zr
      real::r_nan
      integer::i_nan
      equivalence(i_nan, r_nan)

      i_nan = 2147483647
      ua = r_nan
      va = r_nan
      altw = r_nan
      altub = r_nan
      hgtu = r_nan
      hgtv = r_nan

      if (kds .ne. 1) then
!$OMP CRITICAL(SFIRE_DRIVER_CRIT)
         write (msg, *) 'WARNING: bottom index kds=', kds, ' should be 1?'
         call message(msg)
!$OMP END CRITICAL(SFIRE_DRIVER_CRIT)
      end if

      ite1 = snode(ite, ide, 1)
      jte1 = snode(jte, jde, 1)

      call print_3d_stats(its, ite1, kds, kde, jts, jte, ims, ime, kms, kme, jms, jme, u, 'wind U in')
      call print_3d_stats(its, ite, kds, kde, jts, jte1, ims, ime, kms, kme, jms, jme, v, 'wind V in')

      if (fire_print_msg .gt. 0) then
!$OMP CRITICAL(SFIRE_DRIVER_CRIT)
         write (msg, '(a,f7.2,a)') 'interpolate_atm2fire: log-interpolation of wind to', fire_wind_height, ' m'
         call message(msg)
!$OMP END CRITICAL(SFIRE_DRIVER_CRIT)
      end if

      itst = ifval(ids .eq. its, its, its - 1)
      itet = ifval(ide .eq. ite, ite, ite + 1)
      jtst = ifval(jds .ge. jts, jts, jts - 1)
      jtet = ifval(jde .eq. jte, jte, jte + 1)

      itsu = ifval(ids .eq. its, its + 1, its)
      iteu = ifval(ide .eq. ite, ite, ite + 1)
      jtsu = ifval(jds .ge. jts, jts, jts - 1)
      jteu = ifval(jde .eq. jte, jte, jte + 1)

      jtsv = ifval(jds .eq. jts, jts + 1, jts)
      jtev = ifval(jde .eq. jte, jte, jte + 1)
      itsv = ifval(ids .ge. its, its, its - 1)
      itev = ifval(ide .eq. ite, ite, ite + 1)

      if (fire_print_msg .ge. 1) then
!$OMP CRITICAL(SFIRE_DRIVER_CRIT)
         write (msg, 7001) 'atm input  ', 'tile', its, ite, jts, jte
         call message(msg)
         write (msg, 7001) 'altw       ', 'tile', itst, itet, jtst, jtet
         call message(msg)
!$OMP END CRITICAL(SFIRE_DRIVER_CRIT)
      end if
7001  format(a, ' dimensions ', a4, ':', i6, ' to ', i6, ' by ', i6, ' to ', i6)

      kdmax = kde - 1

      do j = jtst, jtet
         do k = kds, kdmax + 1
            do i = itst, itet
         !!!! ALTERADO POR ISILDA CUNHA MENEZES
               ! altw(i,k,j) = (ph(i,k,j)+phb(i,k,j))/g
               ! SUBSTITUIDO POR
               altw(i, k, j) = phb(i, k, j)
         !!!!!!
            end do
         end do
      end do

      if (fire_print_msg .ge. 1) then
!$OMP CRITICAL(SFIRE_DRIVER_CRIT)
         write (msg, 7001) 'u interp at', 'tile', itsu, iteu, jtsu, jteu
         call message(msg)
!$OMP END CRITICAL(SFIRE_DRIVER_CRIT)
      end if

      do j = jtsu, jteu
         do k = kds, kdmax + 1
            do i = itsu, iteu
               altub(i, k, j) = 0.5*(altw(i - 1, k, j) + altw(i, k, j))
            end do
         end do
         do k = kds, kdmax
            do i = itsu, iteu
               hgtu(i, k, j) = 0.5*(altub(i, k, j) + altub(i, k + 1, j)) - altub(i, kds, j)
            end do
         end do
      end do

      if (fire_print_msg .ge. 1) then
!$OMP CRITICAL(SFIRE_DRIVER_CRIT)
         write (msg, 7001) 'v interp at', 'tile', itsv, itev, jtsv, jtev
         call message(msg)
!$OMP END CRITICAL(SFIRE_DRIVER_CRIT)
      end if

      do j = jtsv, jtev
         do k = kds, kdmax + 1
            do i = itsv, itev
               altvb(i, k, j) = 0.5*(altw(i, k, j - 1) + altw(i, k, j))
            end do
         end do
         do k = kds, kdmax
            do i = itsv, itev
               hgtv(i, k, j) = 0.5*(altvb(i, k, j) + altvb(i, k + 1, j)) - altvb(i, kds, j)
            end do
         end do
      end do

      logfwh = log(fire_wind_height)

      do j = jtsu, jteu
         do i = itsu, iteu
            zr = 0.5*(z0(i, j) + z0(i - 1, j))
            if (fire_wind_height > zr) then
               do k = kds, kdmax
                  ht = hgtu(i, k, j)
                  if (.not. ht < fire_wind_height) then
                     loght = log(ht)
                     if (k .eq. kds) then
                        logz0 = log(zr)
                        ua(i, j) = u(i, k, j)*(logfwh - logz0)/(loght - logz0)
                     else
                        loglast = log(hgtu(i, k - 1, j))
                        ua(i, j) = u(i, k - 1, j) + (u(i, k, j) - u(i, k - 1, j))*(logfwh - loglast)/(loght - loglast)
                     end if
                     goto 10
                  end if
                  if (k .eq. kdmax) then
                     ua(i, j) = u(i, k, j)
                  end if
               end do
10             continue
            else
               ua(i, j) = 0.
            end if
         end do
      end do

      do j = jtsv, jtev
         do i = itsv, itev
            zr = 0.5*(z0(i, j - 1) + z0(i, j))
            if (fire_wind_height > zr) then
               do k = kds, kdmax
                  ht = hgtv(i, k, j)
                  if (.not. ht < fire_wind_height) then
                     loght = log(ht)
                     if (k .eq. kds) then
                        logz0 = log(zr)
                        va(i, j) = v(i, k, j)*(logfwh - logz0)/(loght - logz0)
                     else
                        loglast = log(hgtv(i, k - 1, j))
                        va(i, j) = v(i, k - 1, j) + (v(i, k, j) - v(i, k - 1, j))*(logfwh - loglast)/(loght - loglast)
                     end if
                     goto 11
                  end if
                  if (k .eq. kdmax) then
                     va(i, j) = v(i, k, j)
                  end if
               end do
11             continue
            else
               va(i, j) = 0.
            end if
         end do
      end do

      ips1 = ifval(ips .eq. ids, ips + 1, ips)
      call continue_at_boundary(1, 1, 0., &
                                its - 2, ite + 2, jts - 2, jte + 2, &
                                ids + 1, ide, jds, jde, &
                                ips1, ipe, jps, jpe, &
                                itsu, iteu, jtsu, jteu, &
                                itsou, iteou, jtsou, jteou, &
                                ua)

      jps1 = ifval(jps .eq. jds, jps + 1, jps)
      call continue_at_boundary(1, 1, 0., &
                                its - 2, ite + 2, jts - 2, jte + 2, &
                                ids, ide, jds + 1, jde, &
                                ips, ipe, jps1, jpe, &
                                itsv, itev, jtsv, jtev, &
                                itsov, iteov, jtsov, jteov, &
                                va)

      do j = jts, jte1
         do i = its, ite1
            uah(i, j) = ua(i, j)
            vah(i, j) = va(i, j)
         end do
      end do

!!!!!$OMP CRITICAL(SFIRE_DRIVER_CRIT)

      !write (msg, 12) 'atm mesh wind U at', fire_wind_height, ' m'
      !print *,'LFR-DBG ua: ',size(ua,1),size(ua,2)
      !print *,'LFR-DBG ua: ',maxval(ua),minval(ua)
      !LFR call print_2d_stats(itsou, iteou, jtsou, jteou, its - 2, ite + 2, jts - 2, jte + 2, ua, msg)
      !write (msg, 12) 'atm mesh wind V at', fire_wind_height, ' m'
      !LFR call print_2d_stats(itsov, iteov, jtsov, jteov, its - 2, ite + 2, jts - 2, jte + 2, va, msg)
!12    format(a, f6.2, a)
      !LFR call print_2d_stats(its, ite1, jts, jte, ims, ime, jms, jme, uah, 'UAH')
      !LFR call print_2d_stats(its, ite, jts, jte1, ims, ime, jms, jme, vah, 'VAH')

!!!!!!$OMP END CRITICAL(SFIRE_DRIVER_CRIT)

      call interpolate_2d( &
         its - 2, ite + 2, jts - 2, jte + 2, &
         itsou, iteou, jtsou, jteou, &
         ifms, ifme, jfms, jfme, &
         ifts, ifte, jfts, jfte, &
         ir, jr, &
         real(ids), real(jds), ifds - 0.5, jfds + (jr - 1)*0.5, &
         ua, &
         uf)

      call interpolate_2d( &
         its - 2, ite + 2, jts - 2, jte + 2, &
         itsov, iteov, jtsov, jteov, &
         ifms, ifme, jfms, jfme, &
         ifts, ifte, jfts, jfte, &
         ir, jr, &
         real(ids), real(jds), ifds + (ir - 1)*0.5, jfds - 0.5, &
         va, &
         vf)

      !LFR call print_2d_stats_vec(ifts, ifte, jfts, jfte, ifms, ifme, jfms, jfme, uf, vf, 'fire wind (m/s)')

      return

   end subroutine interpolate_atm2fire

   subroutine apply_windrf( &
      ifms, ifme, jfms, jfme, &
      ifts, ifte, jfts, jfte, &
      nfuel_cat, uf, vf)
      integer:: &
         ifms, ifme, jfms, jfme, &
         ifts, ifte, jfts, jfte
      real, intent(in), dimension(ifms:ifme, jfms:jfme)::nfuel_cat
      real, intent(inout), dimension(ifms:ifme, jfms:jfme)::uf, vf

      integer::i, j, k
!print*,"Estou na rotina windrf"

      do j = jfts, jfte
         do i = ifts, ifte
            k = int(nfuel_cat(i, j))
            if (k .lt. no_fuel_cat) then
               uf(i, j) = uf(i, j)*windrf(k)
               vf(i, j) = vf(i, j)*windrf(k)
               !  print*,"windrf k ",k
               !  print*,"windrf(k) ",windrf(k)
               !  print*,"windrf uf(i,j) ", uf(i,j)
            else
               !  print*,"windrf ij ",i,j
               uf(i, j) = 0.
               vf(i, j) = 0.
               !  print*,"windrf ELSE uf(i,j) ", uf(i,j)
               !  print*,"windrf ELSE vf(i,j) ", vf(i,j)
            end if
         end do
      end do
!print*,"sai da rotina windrf"
   end subroutine apply_windrf

   subroutine setup_wind_log_interpolation( &
      ids, ide, jds, jde, &
      ims, ime, jms, jme, &
      ips, ipe, jps, jpe, &
      its, ite, jts, jte, &
      ifds, ifde, jfds, jfde, &
      ifms, ifme, jfms, jfme, &
      ifts, ifte, jfts, jfte, &
      ir, jr, &
      z0, &
      nfuel_cat, &
      fz0, fwh, wz0)

      integer, intent(in):: &
         ids, ide, jds, jde, &
         ims, ime, jms, jme, &
         ips, ipe, jps, jpe, &
         its, ite, jts, jte, &
         ifds, ifde, jfds, jfde, &
         ifms, ifme, jfms, jfme, &
         ifts, ifte, jfts, jfte, &
         ir, jr
      real, intent(in), dimension(ims:ime, jms:jme)::z0
      real, intent(in), dimension(ifms:ifme, jfms:jfme)::nfuel_cat
      real, intent(out), dimension(ifms:ifme, jfms:jfme):: &
         fz0, &
         fwh, &
         wz0

      integer::i, j, ii, jj, k, id = 0
      character(len=128)::msg
      real::r

      if (.not. have_fuel_cats) call crash('setup_wind_log_interpolation: fuel categories not yet set')

      select case (fire_wind_log_interp)

      case (1)
         call message('fire_wind_log_interp=1: log interpolation on fire mesh, roughness and wind height from fuel categories')
         do j = jfts, jfte
            do i = ifts, ifte
               k = int(nfuel_cat(i, j))
               if (k .ge. no_fuel_cat .and. k .le. no_fuel_cat2) then
                  fz0(i, j) = -1.
                  fwh(i, j) = -1.
               elseif (k < 1 .or. k > nfuelcats) then
!$OMP CRITICAL(SFIRE_ATM_CRIT)
                  write (msg, *) 'i,j,nfuel_cat,nfuelcats=', i, j, k, nfuelcats
!$OMP END CRITICAL(SFIRE_ATM_CRIT)
                  call message(msg)
                  call crash('setup_wind_log_interpolation: fuel category out of bounds')
               else
                  fz0(i, j) = fcz0(k)
                  fwh(i, j) = fcwh(k)
               end if
            end do
         end do
         do j = jts, jte
            do i = its, ite
               do jj = (j - 1)*jr + 1, (j - 1)*jr + jr
                  do ii = (i - 1)*ir + 1, (i - 1)*ir + ir
                     wz0(ii, jj) = z0(i, j)
                  end do
               end do
            end do
         end do
      case (2)
         call message('fire_wind_log_interp=2: log interpolation on fire mesh'// &
                      'piecewise constant roughness from landuse, constant fire_wind_height')
         do j = jts, jte
            do i = its, ite
               do jj = (j - 1)*jr + 1, (j - 1)*jr + jr
                  do ii = (i - 1)*ir + 1, (i - 1)*ir + ir
                     fz0(ii, jj) = z0(i, j)
                     wz0(ii, jj) = z0(i, j)
                  end do
               end do
            end do
         end do
         do j = jfts, jfte
            do i = ifts, ifte
               k = int(nfuel_cat(i, j))
               if (k .lt. no_fuel_cat) then
                  fwh(i, j) = fcwh(k)
               else
                  fz0(i, j) = -1.
                  fwh(i, j) = -1.
               end if
            end do
         end do

      case (3)
         call message('fire_wind_log_interp=3: log interpolation on fire mesh,'// &
                      ' interpolated roughness from landuse, constant fire_wind_height')
         call interpolate_z2fire(id, 1, &
                                 ids, ide, jds, jde, &
                                 ims, ime, jms, jme, &
                                 ips, ipe, jps, jpe, &
                                 its, ite, jts, jte, &
                                 ifds, ifde, jfds, jfde, &
                                 ifms, ifme, jfms, jfme, &
                                 ifts, ifte, jfts, jfte, &
                                 ir, jr, &
                                 z0, &
                                 fz0)
         do j = jfts, jfte
            do i = ifts, ifte
               k = int(nfuel_cat(i, j))
               if (k .ne. no_fuel_cat) then
                  fwh(i, j) = fcwh(k)
               else
                  fz0(i, j) = -1.
                  fwh(i, j) = -1.
               end if
               wz0(i, j) = 0
            end do
         end do

      case (4)
         call message('fire_wind_log_interp=4: log interpolation on atmospheric'// &
                      ' mesh, roughness from landuse, constant fire_wind_height')
         return

      case default
         !$OMP CRITICAL(SFIRE_ATM_CRIT)
         write (msg, *) 'setup_wind_log_interpolation: invalid fire_wind_log_interp=', fire_wind_log_interp
         !$OMP END CRITICAL(SFIRE_ATM_CRIT)
         call crash('msg')

      end select

      select case (fire_use_windrf)

      case (0)
         call message('setup_wind_log_interpolation: not using wind reduction factors')

      case (1)
         call message('setup_wind_log_interpolation: multiplying wind by reduction factors')

      case (2)
         call message('setup_wind_log_interpolation: resetting wind interpolation height from wind reduction factors')
         do j = jfts, jfte
            do i = ifts, ifte
               k = int(nfuel_cat(i, j))
               if (k .ne. no_fuel_cat) then
                  fwh(i, j) = fz0(i, j)**(1.-windrf(k))*fire_wind_height**windrf(k)

                  if (.not. fz0(i, j) > 0. .or. .not. fwh(i, j) > fz0(i, j)) then
!$OMP CRITICAL(SFIRE_ATM_CRIT)
                     write (msg, *) 'category ', k, 'windrf=', windrf(k), ' fire_wind_height=', fire_wind_height
                     call message(msg, level=-1)
                     write (msg, *) 'i=', i, ' j=', j, ' fz0(i,j)=', fz0(i, j), ' fwh(i,j)=', fwh(i, j)
                     call message(msg, level=-1)
!$OMP END CRITICAL(SFIRE_ATM_CRIT)
                     call crash('setup_wind_log_interpolation: must have fwh > fz0 > 0')
                  end if

               end if
            end do
         end do

      case (3)
         if (fire_wind_log_interp .eq. 2 .or. fire_wind_log_interp .eq. 3) then
            call message('setup_wind_log_interpolation: adjusting wind interpolation height for LANDUSE roughness height')
            do j = jfts, jfte
               do i = ifts, ifte
                  k = int(nfuel_cat(i, j))
                  if (k .lt. no_fuel_cat) then
                     r = log(fcwh(k)/fcz0(k))/log(fire_wind_height/fcz0(k))
                     fwh(i, j) = fz0(i, j)**(1.-r)*fire_wind_height**r

                     if (.not. fz0(i, j) > 0. .or. .not. fwh(i, j) > fz0(i, j)) then
!$OMP CRITICAL(SFIRE_ATM_CRIT)
                        write (msg, *) 'category ', k, 'roughness ', fcz0(k), ' midflame height ', fcwh(k), ' fire_wind_height=' &
                           , fire_wind_height
                        call message(msg, level=-1)
                        write (msg, *) 'computed wind reduction factor ', r
                        call message(msg, level=-1)
                        write (msg, *) 'i=', i, ' j=', j, ' fz0(i,j)=', fz0(i, j), ' fwh(i,j)=', fwh(i, j)
                        call message(msg, level=-1)
!$OMP END CRITICAL(SFIRE_ATM_CRIT)
                        call crash('setup_wind_log_interpolation: must have fwh > fz0 > 0')
                     end if

                  end if
               end do
            end do
         else
            call message('setup_wind_log_interpolation: not using wind reduction factors')
         end if

      case (4)
         call message('setup_wind_log_interpolation: resetting wind interpolation height from wind reduction factors (DVM NEW)')
         do j = jfts, jfte
            do i = ifts, ifte
               k = int(nfuel_cat(i, j))
               if (k .ne. no_fuel_cat) then
                  fwh(i, j) = fz0(i, j)**(1.-windrf(k))*fire_wind_height**windrf(k)

                  if (.not. fz0(i, j) > 0. .or. .not. fwh(i, j) > fz0(i, j)) then
!$OMP CRITICAL(SFIRE_ATM_CRIT)
                     write (msg, *) 'category ', k, 'windrf=', windrf(k), ' fire_wind_height=', fire_wind_height
                     call message(msg, level=-1)
                     write (msg, *) 'i=', i, ' j=', j, ' fz0(i,j)=', fz0(i, j), ' fwh(i,j)=', fwh(i, j)
                     call message(msg, level=-1)
!$OMP END CRITICAL(SFIRE_ATM_CRIT)
                     call crash('setup_wind_log_interpolation: must have fwh > fz0 > 0')
                  end if

               end if
            end do
         end do

      case default
         !$OMP CRITICAL(SFIRE_ATM_CRIT)
         write (msg, *) 'setup_wind_log_interpolation: invalid fire_use_windrf=', fire_use_windrf
         !$OMP END CRITICAL(SFIRE_ATM_CRIT)
         call crash('msg')

      end select

      do j = jfts, jfte
         do i = ifts, ifte
            k = int(nfuel_cat(i, j))
            if (k .lt. no_fuel_cat) then
               if (.not. fz0(i, j) > 0. .or. .not. fwh(i, j) > fz0(i, j)) then
!$OMP CRITICAL(SFIRE_ATM_CRIT)
                  write (msg, *) 'i=', i, ' j=', j, ' fz0(i,j)=', fz0(i, j), ' fwh(i,j)=', fwh(i, j)
                  call message(msg, level=-1)
!$OMP END CRITICAL(SFIRE_ATM_CRIT)
                  call crash('setup_wind_log_interpolation: must have fwh > fz0 > 0')
               end if
            else
               if (.not. fwh(i, j) < 0.) then
!$OMP CRITICAL(SFIRE_ATM_CRIT)
                  write (msg, *) 'i=', i, ' j=', j, ' fz0(i,j)=', fz0(i, j), ' fwh(i,j)=', fwh(i, j)
                  call message(msg, level=-1)
!$OMP END CRITICAL(SFIRE_ATM_CRIT)
                  call crash('setup_wind_log_interpolation: no fuel must be signalled by fwh<0')
               end if
            end if
         end do
      end do

      have_wind_log_interpolation = .true.

      call print_2d_stats(ifts, ifte, jfts, jfte, ifms, ifme, jfms, jfme, fz0, 'setup_wind_log:fz0')
      call print_2d_stats(ifts, ifte, jfts, jfte, ifms, ifme, jfms, jfme, fz0, 'setup_wind_log:fwh')

   end subroutine setup_wind_log_interpolation

   subroutine interpolate_wind2fire_height(id, &
                                           ids, ide, kds, kde, jds, jde, &
                                           ims, ime, kms, kme, jms, jme, &
                                           ips, ipe, jps, jpe, &
                                           its, ite, jts, jte, &
                                           ifds, ifde, jfds, jfde, &
                                           ifms, ifme, jfms, jfme, &
                                           ifps, ifpe, jfps, jfpe, &
                                           ifts, ifte, jfts, jfte, &
                                           ir, jr, &
                                           u_frame, v_frame, &
                                           u, v, ph, phb, &
                                           fz0, fwh, &
                                           uf, vf)

      implicit none

      integer, intent(in):: id, &
                            ids, ide, kds, kde, jds, jde, &
                            ims, ime, kms, kme, jms, jme, &
                            ips, ipe, jps, jpe, &
                            its, ite, jts, jte, &
                            ifds, ifde, jfds, jfde, &
                            ifms, ifme, jfms, jfme, &
                            ifps, ifpe, jfps, jfpe, &
                            ifts, ifte, jfts, jfte, &
                            ir, jr
      real, intent(in):: u_frame, v_frame
      real, intent(in), dimension(ims:ime, kms:kme, jms:jme):: &
         u, v, &
         ph, phb
      real, intent(in), dimension(ifms:ifme, jfms:jfme):: &
         fz0, &
         fwh
      real, intent(out), dimension(ifms:ifme, jfms:jfme):: &
         uf, &
         vf

      integer:: i, j, k, jcb, jcm, icb, icm, kdmax, kmin, kmax
      integer::itst, itet, jtst, jtet
      integer::iftst, iftet, jftst, jftet
      real:: wjcb, wjcm, wicb, wicm, ht, i_g, loght, zr, ht_last, logwh, wh, loght_last, uk, vk, uk1, vk1, z0, logz0
      real, dimension(its - 1:ite + 1, kds:kde, jts - 1:jte + 1)::z, zw
      character(len=128)::msg

      if (.not. have_wind_log_interpolation) call crash('interpolate_wind2fire_height: wind_log_interpolation must be set up first')

      kdmax = kde - 1

      itst = ifval(ids .eq. its, its, its - 1)
      itet = ifval(ide .eq. ite, ite, ite + 1)
      jtst = ifval(jds .ge. jts, jts, jts - 1)
      jtet = ifval(jde .eq. jte, jte, jte + 1)

      i_g = 1./g
      do j = jtst, jtet
         do k = kds, kdmax + 1
            do i = itst, itet
        !!!!! ALTERADO POR ISILDA CUNHA MENEZES
               zw(i, k, j) = (ph(i, k, j) + phb(i, k, j))*i_g
               ! SUBSTITUIDO POR
               zw(i, k, j) = phb(i, k, j)
        !!!!!!!
            end do
         end do
         do k = kds, kdmax
            do i = itst, itet
               z(i, k, j) = (zw(i, k, j) + zw(i, k + 1, j))*0.5 - zw(i, kds, j)
            end do
         end do
      end do

      iftst = ifval(ifds .eq. ifts, ifts + ir/2, ifts)
      iftet = ifval(ifde .eq. ifte, ifte - ir/2, ifte)
      jftst = ifval(jfds .ge. jfts, jfts + jr/2, jfts)
      jftet = ifval(jfde .eq. jfte, jfte - jr/2, jfte)

      do j = jfts, jfte
         do i = ifts, ifte
            uf(i, j) = 0.
            vf(i, j) = 0.
         end do
      end do

      kmin = kde
      kmax = kds

      loop_j: do j = jftst, jftet
         call coarse(j, jr, -2, jcb, wjcb)
         call coarse(j, jr, ir, jcm, wjcm)
         loop_i: do i = iftst, iftet
            call coarse(i, ir, -2, icb, wicb)
            call coarse(i, ir, ir, icm, wicm)
            z0 = fz0(i, j)
            wh = fwh(i, j)

            if (wh > z0 .and. z0 > 0) then

               ht_last = z0
               loop_k: do k = kds, kdmax

                  ht = interpolate_h(its - 1, ite + 1, kds, kde, jts - 1, jte + 1, icm, k, jcm, wicm, wjcm, z)
                  if (.not. ht < wh) exit loop_k
                  ht_last = ht
               end do loop_k

               if (k .gt. kdmax) then
                  goto 91
               end if
               kmin = min(k, kmin)
               kmax = max(k, kmax)

               logz0 = log(z0)
               logwh = log(wh)
               loght_last = log(ht_last)
               loght = log(ht)

               uk = interpolate_h(ims, ime, kms, kme, jms, jme, icb, k, jcm, wicb, wjcm, u) - u_frame

               vk = interpolate_h(ims, ime, kms, kme, jms, jme, icm, k, jcb, wicm, wjcb, v) - v_frame

               if (k .gt. kds) then

                  uk1 = interpolate_h(ims, ime, kms, kme, jms, jme, icb, k - 1, jcm, wicb, wjcm, u)

                  vk1 = interpolate_h(ims, ime, kms, kme, jms, jme, icm, k - 1, jcb, wicm, wjcb, v)
               else
                  uk1 = 0.
                  vk1 = 0.
               end if

               uf(i, j) = uk1 + (uk - uk1)*(logwh - loght_last)/(loght - loght_last)
               vf(i, j) = vk1 + (vk - vk1)*(logwh - loght_last)/(loght - loght_last)

            else

               uf(i, j) = 0.
               vf(i, j) = 0.
            end if

         end do loop_i
      end do loop_j

!$OMP CRITICAL(SFIRE_ATM_CRIT)
      write (msg, *) 'wind interpolated from layers', kmin, ' to ', kmax
      call message(msg)
!$OMP END CRITICAL(SFIRE_ATM_CRIT)

      return

91    call crash('interpolate_wind2fire_height: fire wind height too large, increase kdmax or atm height')
92    continue
!$OMP CRITICAL(SFIRE_ATM_CRIT)
      write (msg, *) 'fz0(', i, j, ')=', fz0(i, j), 'fwh(', i, j, ')=', fwh(i, j)
      call message(msg)
!$OMP END CRITICAL(SFIRE_ATM_CRIT)
      call crash('interpolate_wind2fire_height: must have fire wind height > roughness height > 0')

   contains

      real function interpolate_h(ims, ime, kms, kme, jms, jme, ic, kc, jc, wic, wjc, a)

         implicit none

         integer, intent(in)::ims, ime, jms, kms, kme, jme, ic, kc, jc
         real, intent(in)::wic, wjc, a(ims:ime, kms:kme, jms:jme)

         interpolate_h = &
            a(ic, kc, jc)*wic*wjc + &
            a(ic, kc, jc + 1)*wic*(1.-wjc) + &
            a(ic + 1, kc, jc)*(1.-wic)*wjc + &
            a(ic + 1, kc, jc + 1)*(1.-wic)*(1.-wjc)
      end function interpolate_h

      subroutine coarse(ix, nr, ia, ic, w)

         implicit none
         integer, intent(in)::ix, nr, ia
         integer, intent(out)::ic
         real, intent(out)::w

         real:: c, a

         a = (ia + 1)*0.5
         c = 1 + (ix - a)/nr
         ic = floor(c)
         w = (1 + ic) - c

      end subroutine coarse

   end subroutine interpolate_wind2fire_height

   subroutine find_trees_fmesh(id, &
                               ids, ide, kds, kde, jds, jde, &
                               ims, ime, kms, kme, jms, jme, &
                               ips, ipe, jps, jpe, &
                               its, ite, jts, jte, &
                               ifds, ifde, jfds, jfde, &
                               ifms, ifme, jfms, jfme, &
                               ifps, ifpe, jfps, jfpe, &
                               ifts, ifte, jfts, jfte, &
                               ir, jr, &
                               can_top, nfuel_cat)

      implicit none

      integer, intent(in):: id, &
                            ids, ide, kds, kde, jds, jde, &
                            ims, ime, kms, kme, jms, jme, &
                            ips, ipe, jps, jpe, &
                            its, ite, jts, jte, &
                            ifds, ifde, jfds, jfde, &
                            ifms, ifme, jfms, jfme, &
                            ifps, ifpe, jfps, jfpe, &
                            ifts, ifte, jfts, jfte, &
                            ir, jr

      real, intent(in), dimension(ifms:ifme, jfms:jfme)::nfuel_cat
      real, intent(inout), dimension(ifms:ifme, jfms:jfme)::can_top

      integer :: i, j, iveg, tree_type
      real, dimension(8) :: z_max

      z_max = (/10., 20., 15., 24., 22.5, 10., 2.2, 0.73/)

      do j = jfts, jfte
         do i = ifts, ifte

            iveg = nfuel_cat(i, j)
            tree_type = itree(iveg)

            if (tree_type /= 0) then

               can_top(i, j) = z_max(tree_type)

            else

               can_top(i, j) = 0

            end if

         end do
      end do

   end subroutine find_trees_fmesh

   subroutine massman_fwh(id, &
                          ids, ide, kds, kde, jds, jde, &
                          ims, ime, kms, kme, jms, jme, &
                          ips, ipe, jps, jpe, &
                          its, ite, jts, jte, &
                          ifds, ifde, jfds, jfde, &
                          ifms, ifme, jfms, jfme, &
                          ifps, ifpe, jfps, jfpe, &
                          ifts, ifte, jfts, jfte, &
                          ir, jr, &
                          cuf, cvf, &
                          nfuel_cat, &
                          can_top, &
                          uf, vf)

      implicit none

      integer, intent(in):: id, &
                            ids, ide, kds, kde, jds, jde, &
                            ims, ime, kms, kme, jms, jme, &
                            ips, ipe, jps, jpe, &
                            its, ite, jts, jte, &
                            ifds, ifde, jfds, jfde, &
                            ifms, ifme, jfms, jfme, &
                            ifps, ifpe, jfps, jfpe, &
                            ifts, ifte, jfts, jfte, &
                            ir, jr

      real, intent(in), dimension(ifms:ifme, jfms:jfme):: &
         cuf, cvf, &
         nfuel_cat, &
         can_top

      real, intent(inout), dimension(ifms:ifme, jfms:jfme)::uf, vf

      real :: xi_og, k, c1, c3, c2, mid_fh
      real :: ha, z_og, drag_index, cantop_frac, C_surf, N, cd, zeta_sum, u_scale, z_top
      real, dimension(8) :: xi_max, sigma_u, sigma_l, PAI, big_cd, z_max, d1
      real, dimension(20) :: z_agl, xi, fa, zeta, Ub, Ut, u
      integer :: iveg, i, j, z, znum

      xi_og = 0.0025
      k = 0.40
      c1 = 0.38
      c3 = 15
      c2 = c1 + k/log(xi_og)

      znum = 20

      xi_max = (/0.36, 0.60, 0.58, 0.60, 0.84, 0.60, 0.94, 0.62/)
      sigma_u = (/0.60, 0.30, 0.20, 0.10, 0.13, 0.38, 0.03, 0.50/)
      sigma_l = (/0.20, 0.10, 0.20, 0.27, 0.30, 0.16, 0.60, 0.45/)
      PAI = (/3.28, 2.41, 2.14, 3.78, 4.93, 5.73, 2.94, 3.10/)
      big_cd = (/0.20, 0.25, 0.27, 0.20, 0.09, 0.20, 0.30, 0.30/)
      z_max = (/10., 20., 15., 24., 22.5, 10., 2.2, 0.73/)
      d1 = (/0., 0., -0.7, 0., 2., 0., 0., 0./)

      do j = jfts, jfte
         do i = ifts, ifte

            z_top = can_top(i, j)

            if (z_top == 0) then
               cycle
            end if

            iveg = itree(nfuel_cat(i, j))

            do z = 1, znum
               if (z > 1) then
                  z_agl(z) = z_agl(z - 1) + z_top/(znum - 1)
               else
                  z_agl(z) = 0
               end if
            end do

            xi = z_agl/z_top

            do z = 1, znum
               if (xi(z) >= xi_max(iveg)) then
                  fa(z) = exp((-(xi(z) - xi_max(iveg))**2)/(sigma_u(iveg)**2))
               else
                  fa(z) = exp((-(xi_max(iveg) - xi(z))**2)/(sigma_l(iveg)**2))
               end if
            end do

            ha = sum((fa/sum(fa))*PAI(iveg))

            z_og = xi_og*z_top

            drag_index = big_cd(iveg)*PAI(iveg)

            cantop_frac = c1 - c2*exp(-c3*drag_index)

            C_surf = (cantop_frac**2)*2

            N = drag_index/C_surf

            cd = big_cd(iveg)

            zeta_sum = 0

            do z = 1, znum
               zeta(z) = cd*(fa(z)/sum(fa))*PAI(iveg) + zeta_sum
               zeta_sum = zeta(z)
            end do

            do z = 1, znum
               if (xi(z) >= xi_og) then
                  Ub(z) = log(xi(z)/xi_og)/log(1/xi_og)
               else
                  Ub(z) = 0
               end if
            end do

            Ut = cosh((N*zeta)/drag_index)/cosh(N)

            u = Ut*Ub

            mid_fh = fueldepthm(nfuel_cat(i, j))

            do z = 1, znum
               if (z_agl(z) >= mid_fh) then
                  u_scale = (((z_agl(z - 1) - mid_fh)*(u(z - 1) - u(z)))/(z_agl(z) - z_agl(z - 1))) + u(z - 1)
                  exit
               end if
            end do

            uf(i, j) = cuf(i, j)*u_scale
            vf(i, j) = cvf(i, j)*u_scale

         end do
      end do

   end subroutine massman_fwh

end module module_fr_sfire_atm
