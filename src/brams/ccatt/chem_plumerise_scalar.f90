
module mod_chem_plumerise_scalar

   implicit none

   private

   public :: plumerise_driver   ! Subroutine

contains

   !----------------------------------------------------------------------------
   ! Plume rise model for vegetation fires (CPTEC/INPE 2005-2006)
   ! Refs.:
   ! Freitas, S. R., K. M. Longo, R. Chatfield, D. Latham, M. A. F. Silva Dias,
   !  M. O. Andreae,
   ! E. Prins, J. C. Santos, R. Gielow and J. A. Carvalho Jr.: Including the
   ! sub-grid scale
   ! plume rise of vegetation fires in low resolution atmospheric transport
   ! models.
   !  Atmospheric Chemistry and Physics and discussion, 6, 2006.
   !-
   ! Freitas, S. R.; Longo, K. M.; M. Andreae. Impact of including the plume
   ! rise of vegetation
   ! fires in numerical simulations of associated atmospheric pollutants.
   ! Geophys. Res. Lett.,
   ! 33, L17808, doi:10.1029/2006GL026608, 2006.
   !----------------------------------------------------------------------------

   subroutine plumerise_driver(mzp, mxp, myp, ia, iz, ja, jz, &
                               srctime1, chem_nspecies, spc_chem_alloc, src, &
                               on, off, nmodes, aer_nspecies, spc_aer_alloc, &
                               aer_bburn, nsrc, bburn, nveg_agreg, tropical_forest, &
                               boreal_forest, savannah, grassland, nzpmax, dtlt, &
                               time, zt_rams, zm_rams, dzt_rams, g, cp, cpor, p00, &
                               rgas, theta, pp, pi0, rv, up, vp, rtgt, lpw, pr_time, &
                               aerosol, chem1_src_g, aer1_g, plume_mean_g, spc_aer_name, &
                               plume_fre_g, &
                               plumerise_flag)
      use mem_chem1, only: &
         chem1_src_vars      ! Type
      use mem_aer1, only: &
         aer1_vars           ! Type
      use mem_plume_chem1, only: &
         plume_mean_vars, &     ! Type
         plume_fre_vars, &     ! Type
         iflam_frac, &
         imean_frp, &
         istd_frp, &
         imean_size, &
         istd_size
      use node_mod, only: &
         mynum,          & ! INTENT(IN)
         mchnum,        &
         master_num
      use mem_sfire, only: sfire

      integer, intent(IN)    :: mzp
      integer, intent(IN)    :: mxp
      integer, intent(IN)    :: myp
      integer, intent(IN)    :: ia
      integer, intent(IN)    :: iz
      integer, intent(IN)    :: ja
      integer, intent(IN)    :: jz
      integer, intent(IN)    :: plumerise_flag
      ! grid_dims
      integer, intent(IN)    :: nzpmax
      ! mem_basic
      real, intent(IN)    :: theta(:, :, :)
      real, intent(IN)    :: pp(:, :, :)
      real, intent(IN)    :: pi0(:, :, :)
      real, intent(IN)    :: rv(:, :, :)
      real, intent(IN)    :: up(:, :, :)!srf-AWE
      real, intent(IN)    :: vp(:, :, :)!srf-AWE
      ! mem_grid
      real, intent(IN)    :: rtgt(:, :)
      integer, intent(IN)    :: lpw(:, :)
      real, intent(IN)    :: dtlt
      real, intent(IN)    :: time

      real, intent(IN)    :: zt_rams(:)
      real, intent(IN)    :: zm_rams(:)
      real, intent(IN)    :: dzt_rams(:)
      ! chem_sources
      real, intent(IN)    :: pr_time
      real, intent(IN)    :: srctime1
      ! chem1_list
      integer, intent(IN)    :: chem_nspecies
      integer, intent(IN)    :: spc_chem_alloc(:, :)
      integer, intent(IN)    :: src
      integer, intent(IN)    :: on
      integer, intent(IN)    :: off
      ! aer1_list
      integer, intent(IN)    :: nmodes
      integer, intent(IN)    :: aer_nspecies
      integer, intent(IN)    :: spc_aer_alloc(:, :, :)
      integer, intent(IN)    :: aer_bburn
      ! mem_chem1
      integer, intent(IN)    :: nsrc
      type(chem1_src_vars), intent(INOUT) :: chem1_src_g(:, :, :)
      integer, intent(IN)    :: bburn
      ! mem_aer1
      type(aer1_vars), intent(INOUT) :: aer1_g(:, :)
      integer, intent(IN)    :: aerosol
      ! mem_plume_chem1
      type(plume_mean_vars), intent(IN)    :: plume_mean_g(:)
      type(plume_fre_vars), intent(INOUT) :: plume_fre_g(5)
      integer, intent(IN)    :: nveg_agreg
      integer, intent(IN)    :: tropical_forest
      integer, intent(IN)    :: boreal_forest
      integer, intent(IN)    :: savannah
      integer, intent(IN)    :: grassland
      ! rconstants
      real, intent(IN)    :: g
      real, intent(IN)    :: cp
      real, intent(IN)    :: cpor
      real, intent(IN)    :: p00
      real, intent(IN)    :: rgas
      character(LEN=8), intent(IN) :: spc_aer_name(nmodes, aer_nspecies)

      logical, parameter :: IS2PRINT = .false. !.TRUE.
      character(LEN=8)   :: ctime
      character(LEN=4)   :: cmynum
      integer            :: nrec

      character(LEN=*), parameter :: h = "**(plumerise_driver)**"
      character(LEN=8) :: c0
      !TO
      integer, parameter :: nkp = 200
      integer, parameter :: ntime = 200
      integer            :: n_setgrid = 0

      real :: dz
      real :: zt(nkp)
      real :: zm(nkp)
      real :: dzt(nkp)
      real :: dzm(nkp)
      integer :: i,j
      character(len=*), parameter :: fm = '(2(I3.3,1X),3(E15.5,1X))'

      if (mod(time + .001, pr_time) .le. (dtlt) .or. time .lt. .01 .or. &
          abs(time - srctime1) .lt. 1.e-5 .or. sfire>0) then
!if (mchnum==master_num) print *,'LFR-DBG - Chamando plumerise = timestep ',time
!!$       WRITE(c0,"(f8.1)") time
!!$       CALL MsgOne(h,' -----------------------------------------------------')
!!$       CALL MsgOne(h,' plumerise tendencies updated')
!!$       CALL MsgOne(h,' time = '//TRIM(ADJUSTL(c0))//' seconds')
!!$       CALL MsgOne(h,' -----------------------------------------------------')


         call plumerise(mzp, mxp, myp, ia, iz, ja, jz, &
                        theta, pp, pi0, rv, up, vp, rtgt, lpw, zt_rams, zm_rams, &
                        nzpmax, dzt_rams, chem_nspecies, spc_chem_alloc, &
                        src, on, off, nmodes, aer_nspecies, spc_aer_alloc, &
                        aer_bburn, nsrc, chem1_src_g, bburn, aer1_g, aerosol, &
                        plume_mean_g, nveg_agreg, tropical_forest, &
                        boreal_forest, savannah, grassland, g, cp, cpor, &
                        p00, rgas, spc_aer_name, &
                        plume_fre_g, &
                        plumerise_flag, nkp, ntime, n_setgrid, dz, dzm, dzt, zm, zt)
      end if

   end subroutine plumerise_driver
   !-------------------------------------------------------------------------

   subroutine plumerise(m1, m2, m3, ia, iz, ja, jz, &
                        theta, pp, pi0, rv, up, vp, rtgt, lpw, zt_rams, zm_rams, &
                        nzpmax, dzt_rams, chem_nspecies, spc_chem_alloc, src, &
                        on, off, nmodes, aer_nspecies, spc_aer_alloc, aer_bburn, &
                        nsrc, chem1_src_g, bburn, aer1_g, aerosol, plume_mean_g, &
                        nveg_agreg, tropical_forest, boreal_forest, savannah, &
                        grassland, g, cp, cpor, p00, rgas, spc_aer_name, &
                        plume_fre_g, &
                        plumerise_flag, nkp, ntime, n_setgrid, dz, dzm, dzt, zm, zt)
      use mem_chem1, only: &
         chem1_src_vars      ! Type
      use mem_aer1, only: &
         aer1_vars           ! Type
      use mem_plume_chem1, only: &
         plume_mean_vars, &     ! Type
         plume_fre_vars, &     ! Type
         iflam_frac, &
         imean_frp, &
         istd_frp, &
         imean_size, &
         istd_size
      use node_mod, only: &
         mynum,          & ! INTENT(IN)
         mchnum,        &
         master_num

      integer, intent(IN)      :: nkp
      integer, intent(IN)      :: ntime
      integer, intent(INOUT)      :: n_setgrid
      real, intent(INOUT) :: dz
      real, intent(INOUT) :: dzm(nkp)
      real, intent(INOUT) :: dzt(nkp)
      real, intent(INOUT) :: zt(nkp)
      real, intent(INOUT) :: zm(nkp)
      integer, intent(IN)     :: m1
      integer, intent(IN)     :: m2
      integer, intent(IN)     :: m3
      integer, intent(IN)     :: ia
      integer, intent(IN)     :: iz
      integer, intent(IN)     :: ja
      integer, intent(IN)     :: jz
      integer, intent(IN)     :: plumerise_flag
      real, intent(IN)     :: theta(:, :, :)
      real, intent(IN)     :: pp(:, :, :)
      real, intent(IN)     :: pi0(:, :, :)
      real, intent(IN)     :: rv(:, :, :)
      real, intent(IN)     :: up(:, :, :)!srf-AWE
      real, intent(IN)     :: vp(:, :, :)!srf-AWE
      real, intent(IN)     :: rtgt(:, :)
      integer, intent(IN)     :: lpw(:, :)
      real, intent(IN)     :: zt_rams(:)
      real, intent(IN)     :: zm_rams(:)
      ! grid_dims
      integer, intent(IN)     :: nzpmax
      ! mem_grid
      real, intent(IN)     :: dzt_rams(:)
      ! chem1_list
      integer, intent(IN)     :: chem_nspecies
      integer, intent(IN)     :: spc_chem_alloc(:, :)
      integer, intent(IN)     :: src
      integer, intent(IN)     :: on
      integer, intent(IN)     :: off
      ! aer1_list
      integer, intent(IN)     :: nmodes
      integer, intent(IN)     :: aer_nspecies
      integer, intent(IN)     :: spc_aer_alloc(:, :, :)
      integer, intent(IN)     :: aer_bburn
      ! mem_chem1
      integer, intent(IN)     :: nsrc
      type(chem1_src_vars), intent(INOUT)  :: chem1_src_g(:, :, :)
      integer, intent(IN)     :: bburn

      ! mem_aer1
      type(aer1_vars), intent(INOUT) :: aer1_g(:, :)
      integer, intent(IN)    :: aerosol
      ! mem_plume_chem1
      type(plume_mean_vars), intent(IN)    :: plume_mean_g(:)
      type(plume_fre_vars), intent(INOUT) :: plume_fre_g(5)
      integer, intent(IN)    :: nveg_agreg
      integer, intent(IN)    :: tropical_forest
      integer, intent(IN)    :: boreal_forest
      integer, intent(IN)    :: savannah
      integer, intent(IN)    :: grassland
      ! rconstants
      real, intent(IN)    :: g
      real, intent(IN)    :: cp
      real, intent(IN)    :: cpor
      real, intent(IN)    :: p00
      real, intent(IN)    :: rgas

      character(LEN=8), intent(IN) :: spc_aer_name(nmodes, aer_nspecies)

      integer            :: i, j, k, ixx, iveg_ag, imm, k1, k2, ispc, imode, iloop
      integer            :: kmt
      real               :: burnt_area, STD_burnt_area, dz_flam, dzi, FRP, convert_smold_to_flam
      real               :: ztopmax(2)
      real               :: W_VMD(nkp, 2), VMD(nkp, 2)
      real               :: q_smold_kgm2
      integer            :: it1, it2
      integer, parameter :: use_last = 0
      integer, parameter :: izprint = 0 ! if = 0 => no printout
      integer, parameter :: wind_eff = 1

      real    :: area
      integer :: maxtime
      real    :: alpha
      real    :: rsurf
      real    :: fmoist
      real    :: tdur
      real    :: zsurf
      real    :: heating(ntime)
      real    :: thtcon(nkp)
      real    :: rvcon(nkp)
      real    :: picon(nkp)
      real    :: zcon(nkp)
      real    :: zzcon(nkp)
      real    :: qv(nkp)
      real    :: qh(nkp)
      real    :: qi(nkp)
      real    :: qc(nkp)
      real    :: txs(nkp)
      real    :: cvi(nkp)
      real    :: vth(nkp)
      real    :: w(nkp)
      real    :: t(nkp)
      real    :: qsat(nkp)
      real    :: rho(nkp)
      real    :: radius(nkp)
      real    :: visc(nkp)
      real    :: wc(nkp)
      real    :: wt(nkp)
      real    :: est(nkp)
      real    :: pe(nkp)
      real    :: te(nkp)
      real    :: tt(nkp)
      real    :: qvenv(nkp)
      real    :: cvh(nkp)
      real    :: vti(nkp)
      real    :: pke(nkp)
      real    :: the(nkp)
      real    :: thve(nkp)
      real    :: dne(nkp)
      real    :: sc(nkp)
      real    :: sct(nkp)
      real    :: ucon(nkp)
      real    :: vcon(nkp)
      real    :: upe(nkp)
      real    :: vpe(nkp)
      real    :: vel_e(nkp)
      real    :: vel_p(nkp)
      real    :: rad_p(nkp)
      real    :: vel_t(nkp)
      real    :: rad_t(nkp)
      integer omp_get_num_threads, omp_get_thread_num

      !- for biomass burn, only one memory allocation (only one time level)
      it1 = 1

      !- zera o termo fonte associado `as emissoes com plumerise (k>2)
      !- chemistry section
      do ispc = 1, chem_nspecies
         if (spc_chem_alloc(src, ispc) == off) cycle
         chem1_src_g(it1, bburn, ispc)%sc_src(3:m1, :, :) = 0.
      end do

      !- aerosol section
      if (aerosol == 1) then
         do imode = 1, nmodes
            !- only for biomass burning aerosols
            if (spc_aer_alloc(src, imode, aer_bburn) == off) cycle
            aer1_g(imode, aer_bburn)%sc_src(3:m1, :, :) = 0.
         end do
      elseif (aerosol == 2) then ! for MATRIX and only for boc_ocar and boc_bcar species

         do ispc = 1, aer_nspecies
            do imode = 1, nmodes
               if (spc_aer_alloc(src, imode, ispc) == off) cycle

               !-only for bburn aerosols)
               if (spc_aer_name(imode, ispc) == "boc_bcar" .or. &
                   spc_aer_name(imode, ispc) == "boc_ocar") then
                  aer1_g(imode, ispc)%sc_src(3:m1, :, :) = 0.
               end if
            end do
         end do
      end if

      do j = ja, jz
         do i = ia, iz
            convert_smold_to_flam = 0.0

            !- if the max value of flaming is close to zero => there is not
            !  emission with plume rise => cycle
            if (PLUMERISE_flag == 1) then
               if (plume_mean_g(tropical_forest)%flam_frac(i, j) + &
                   plume_mean_g(boreal_forest)%flam_frac(i, j) + &
                   plume_mean_g(savannah)%flam_frac(i, j) + &
                   plume_mean_g(grassland)%flam_frac(i, j) < 1.e-6) cycle
               !- loop over the four types of aggregate biomes with fires
               iloop = nveg_agreg

            elseif (PLUMERISE_flag == 2) then
               if (plume_fre_g(iflam_frac)%pvar(i, j) < 1.e-6) cycle
               iloop = 1

            end if

            do k = 1, m1
               ucon(k) = up(k, i, j)                ! u wind (m/s)
               vcon(k) = vp(k, i, j)                ! v wind (m/s)
               !wcon  (k)=wp(k,i,j)                ! w wind (m/s)
               thtcon(k) = theta(k, i, j)             ! pot temperature (K)
               picon(k) = (pp(k, i, j) + pi0(k, i, j))   ! exner function
               !tmpcon(k)=thtcon(k)*picon(k)/cp    ! temperature (K)
               !dncon (k)=dn0(k,i,j)                ! dry air density (basic state) (kg/m3)
               !prcon (k)=(picon(k)/cp)**cpor*p00  ! pressure (Pa)
               rvcon(k) = rv(k, i, j)                ! water vapor mixing ratio (kg/kg)

               zcon(k) = zt_rams(k)*rtgt(i, j)    ! termod-point height (m)
               zzcon(k) = zm_rams(k)*rtgt(i, j)    ! W-point height      (m)

               zsurf = 0.
            end do ! end do k

            !- get envinronmental state (temp, water vapor mix ratio, ...)
            call get_env_condition(lpw(i, j), m1 - 1, kmt, qvenv, pke, the, te, pe, thve, dne, upe, vpe, vel_e, &
                                   zcon, zt, thtcon, rvcon, picon, ucon, vcon, zm, dzm, dzt, dz, zsurf, &
                                   wind_eff, g, cp, cpor, p00, rgas, nkp, n_setgrid)

            !- only for testing (using an external sounding)
            !CALL get_env_condition_sound(lpw(i,j),m1-1,kmt,qvenv,pke,the,te,pe,thve,dne,upe,vpe,vel_e &
            !                         ,zcon,zt,thtcon,rvcon,picon,ucon,vcon,zm,dzm,dzt,dz,zsurf,n_setgrid,wind_eff)

            !- loop over the four types of aggregate biomes with fires
            do iveg_ag = 1, iloop

               !- verifica a existencia de emissao flaming para um bioma especifico
               if (PLUMERISE_flag == 1) then
                  if (plume_mean_g(iveg_ag)%flam_frac(i, j) < 1.e-6) cycle

                  !-burnt area and standard deviation
                  !burnt_area   = 50.*1e4           !m^2, only for testing with a sounding
                  burnt_area = plume_mean_g(iveg_ag)%fire_size(i, j)

                  STD_burnt_area = 0.!not em use

                  !-number to calculate the flaming emission from the amount emitted
                  !-during the smoldering phase
                  convert_smold_to_flam = plume_mean_g(iveg_ag)%flam_frac(i, j)

               elseif (PLUMERISE_flag == 2) then
                  !-number to calculate the emission during the flaming pahse
                  !-from the amount emitted during the smoldering phase

                  convert_smold_to_flam = plume_fre_g(iflam_frac)%pvar(i, j)
                  !
                  !- check if there is only one fire in a given grid box (=> std =0.)
                  if (plume_fre_g(istd_frp)%pvar(i, j) < 1.0e-6) then

                     !- if yes, we will set it as a 20% of the mean frp as a gross estimation
                     !- of the retrieval uncertainty by the sensors.
                     !- (we are not taking care about the fire size retrieval)
                     !print*,"xx=",I,J,plume_fre_g(istd_frp )%pvar(i,j),0.2*plume_fre_g(imean_frp )%pvar(i,j)
                     plume_fre_g(istd_frp)%pvar(i, j) = 0.2*plume_fre_g(imean_frp)%pvar(i, j)
                  end if
               end if

               !- loop over the minimum and maximum heat fluxes/FRP
               do imm = 1, 2
                  !--------------------
                  !ixx=iveg_ag*10 + imm
                  ixx = 0 + imm
                  !--------------------

                  if (PLUMERISE_flag == 2) then
                     if (imm == 1) then
                        !for imm = 1 => lower injection height
                        burnt_area = max(1.0e4, plume_fre_g(imean_size)%pvar(i, j) - 0.5*plume_fre_g(istd_size)%pvar(i, j))
                        FRP = max(1000., plume_fre_g(imean_frp)%pvar(i, j) - 0.5*plume_fre_g(istd_frp)%pvar(i, j))

                     elseif (imm == 2) then
                        !for imm = 2 => higher injection height
                        burnt_area = max(1.0e4, plume_fre_g(imean_size)%pvar(i, j) + 0.5*plume_fre_g(istd_size)%pvar(i, j))
                        FRP = max(1000., plume_fre_g(imean_frp)%pvar(i, j) + 0.5*plume_fre_g(istd_frp)%pvar(i, j))
                     end if
                  end if
                  ! print*,"imm",imm,plume_fre_g(imean_size)%pvar(i,j), plume_fre_g(istd_size)%pvar(i,j),&
                  !                    plume_fre_g(imean_frp )%pvar(i,j), plume_fre_g(istd_frp )%pvar(i,j)
                  !
                  !-just for testing...
                  !burnt_area =  burnt_area*30.
                  !FRP            =  FRP*600.
                  !-----------
                  !print*,"size(ha)/frp=",i,j,imm, burnt_area/1.e+4, 0.88*1000.*FRP/(burnt_area) ; call flush(6)
                  !write(*,100)"size(ha)/frp(kW/m2)",i,j,imm, burnt_area/1.e+4, 0.88*1.e-3*FRP/(burnt_area)
                  !100 format (1x,A20,1x,3I4,1x,2F15.4)

                  !- get fire properties (burned area, plume radius, heating rates ...)
                  call get_fire_properties(imm, iveg_ag, burnt_area, use_last, &
                                           heating, area, maxtime, alpha, &
                                           rsurf, fmoist, tdur, &
                                           plumerise_flag, FRP, ntime)

                  !------  generates the plume rise    ------

                  if (PLUMERISE_flag == 1) then
                     !-- only one value for eflux of GRASSLAND
                     if (iveg_ag == GRASSLAND .and. imm == 2) then
                        ztopmax(2) = ztopmax(1)
                        ztopmax(1) = zzcon(1)
                        cycle
                     end if
                  end if

                  call makeplume(i, j, kmt, ztopmax(imm), ixx, imm, use_last, w, &
                                 t, pe, qsat, est, rho, zt, te, txs, wc, dz, maxtime, zm, &
                                 wt, dzm, tt, qv, qh, qi, qc, qvenv, radius, alpha, visc, &
                                 cvi, vth, cvh, vti, dzt, sc, area, rsurf, tdur, heating, &
                                 fmoist, vel_e, vel_p, rad_p, vel_t, rad_t, W_VMD, &
                                 izprint, nkp, ntime)

               end do ! enddo of the loop imm

               !- defines the vertical layer where the flaming phase emission will
               !  be released in the 3d atmospheric model
               call set_flam_vert(ztopmax, k1, k2, zzcon, W_VMD, VMD, nkp)

               !- thickness of the vertical layer
               dz_flam = zzcon(k2) - zzcon(k1 - 1)

               !print*,"   z1/z2=",k1,k2,zzcon(k1-1),zzcon(k2)!,convert_smold_to_flam ; call flush(6)
               !print*,"=================================================" ; call flush(6)

               !- emission during flamming phase is evenly distributed between levels k1 and k2
               do k = k1, k2
                  !use this in case the emission src is already in mixing ratio
                  !rhodzi= 1./(dn0(k,i,j) * dz_flam)

                  !use this in case the emission src is tracer density
                  dzi = 1./(dz_flam)

                  !- chemistry section
                  do ispc = 1, chem_nspecies
                     if (spc_chem_alloc(src, ispc) == off) cycle

                     !- get back the smoldering emission in kg/m2  (actually in
                     !   1e-9 kg/m2)

                     !use this in case the emission src is already in mixing
                     ! ratio
                     !q_smold_kgm2 = (rtgt(i,j)/dzt_rams(2) *  dn0(2,i,j)        )*   &
                     !                  chem1_src_g(it1,bburn,ispc,ng)%sc_src(2,i,j)

                     !use this in case the emission src is tracer density
                     q_smold_kgm2 = (rtgt(i, j)/dzt_rams(2))* &
                                    chem1_src_g(it1, bburn, ispc)%sc_src(2, i, j)

                     ! units = already in ppbm,  don't need "fcu" factor
                     chem1_src_g(it1, bburn, ispc)%sc_src(k, i, j) = chem1_src_g(it1, bburn, ispc)%sc_src(k, i, j) + &
                                                                     convert_smold_to_flam* &
                                                                     !plume_mean_g(iveg_ag)%flam_frac(i,j)      *&
                                                                     q_smold_kgm2* &
                                                                     dzi    !use this in case the emission src is tracer density
                     !rhodzi !use this in case the emission src is already in mixing ratio

                  end do
                  !- aerosol section ( only for biomass burning aerosols )
                  if (aerosol == 1) then

                     do imode = 1, nmodes

                        if (spc_aer_alloc(src, imode, aer_bburn) == off) cycle

                        !- get back the smoldering emission in kg/m2
                        !  (actually in 1e-9 kg/m2)

                        !use this in case the emission src is already in
                        !mixing ratio
                        !q_smold_kgm2 = (rtgt(i,j)/dzt_rams(2) * dn0(2,i,j) )*  &
                        !                       aer1_g(imode,aer_bburn)%sc_src(2,i,j)

                        !use this in case the emission src is tracer density

                        q_smold_kgm2 = (rtgt(i, j)/dzt_rams(2))* &
                                       aer1_g(imode, aer_bburn)%sc_src(2, i, j)

                        ! units = already in ppbm,  don't need "fcu" factor

                        aer1_g(imode, aer_bburn)%sc_src(k, i, j) = aer1_g(imode, aer_bburn)%sc_src(k, i, j) + &
                                                                   convert_smold_to_flam* &
                                                                   !plume_mean_g(iveg_ag)%flam_frac(i,j)  *&
                                                                   q_smold_kgm2* &
                                                                   dzi    !use this in case the emission src is tracer density
                        !rhodzi !use this in case the emission src is already in mixing ratio

                     end do

                  elseif (aerosol == 2) then ! for MATRIX and only for boc_ocar and boc_bcar species

                     do ispc = 1, aer_nspecies
                        do imode = 1, nmodes
                           if (spc_aer_alloc(src, imode, ispc) == off) cycle

                           !-only for bburn aerosols - MATRIX MECH=1
                           if (spc_aer_name(imode, ispc) == "boc_bcar" .or. &
                               spc_aer_name(imode, ispc) == "boc_ocar") then

                              q_smold_kgm2 = (rtgt(i, j)/dzt_rams(2))* &
                                             aer1_g(imode, ispc)%sc_src(2, i, j)

                              ! units = already in ppbm,  don't need "fcu" factor
                              aer1_g(imode, ispc)%sc_src(k, i, j) = aer1_g(imode, ispc)%sc_src(k, i, j) + &
                                                                    !plume_mean_g(iveg_ag)%flam_frac(i,j) *&
                                                                    convert_smold_to_flam* &
                                                                    q_smold_kgm2* &
                                                                    dzi     !use this in case the emission src is tracer density
                              !rhodzi !use this in case the emission src is already in mixing ratio

                           end if
                        end do
                     end do
                  end if

               end do

            end do ! enddo do loop em iveg_ag
         end do  ! loop em i
      end do   ! loop em j
   end subroutine plumerise
   !-------------------------------------------------------------------------

   subroutine get_env_condition(k1, k2, kmt, qvenv, pke, the, te, pe, thve, dne, &
                                upe, vpe, vel_e, zcon, zt, thtcon, rvcon, picon, &
                                ucon, vcon, zm, dzm, dzt, dz, zsurf, &! n_setgrid,    &
                                wind_eff, g, cp, cpor, p00, rgas, nkp, n_setgrid)

      integer, intent(IN)    :: nkp
      integer, intent(INOUT)    :: n_setgrid
      integer, intent(IN)    :: k1
      integer, intent(IN)    :: k2
      integer, intent(INOUT) :: kmt
      integer, intent(IN)    :: wind_eff
      ! plumegen_coms
      real, intent(INOUT) :: dz
      real, intent(IN)    :: zsurf
!   INTEGER, INTENT(INOUT) :: n_setgrid
!!$    REAL,    INTENT(INOUT) :: qvenv (nkp)
!!$    REAL,    INTENT(INOUT) :: pke   (nkp)
!!$    REAL,    INTENT(INOUT) :: the   (nkp)
!!$    REAL,    INTENT(INOUT) :: te    (nkp)
!!$    REAL,    INTENT(INOUT) :: pe    (nkp)
!!$    REAL,    INTENT(INOUT) :: thve  (nkp)
!!$    REAL,    INTENT(INOUT) :: dne   (nkp)
!!$    REAL,    INTENT(INOUT) :: upe   (nkp)
!!$    REAL,    INTENT(INOUT) :: vpe   (nkp)
!!$    REAL,    INTENT(INOUT) :: vel_e (nkp)
      real, intent(INOUT) :: zcon(nkp)
!!$    REAL,    INTENT(INOUT) :: zt    (nkp)
!!$    REAL,    INTENT(INOUT) :: thtcon(nkp)
!!$    REAL,    INTENT(INOUT) :: rvcon (nkp)
!!$    REAL,    INTENT(INOUT) :: picon (nkp)
!!$    REAL,    INTENT(INOUT) :: ucon  (nkp)
!!$    REAL,    INTENT(INOUT) :: vcon  (nkp)
!!$    REAL,    INTENT(INOUT) :: zm    (nkp)
!!$    REAL,    INTENT(INOUT) :: dzm   (nkp)
!!$    REAL,    INTENT(INOUT) :: dzt   (nkp)
      real, intent(INOUT) :: qvenv(:)
      real, intent(INOUT) :: pke(:)
      real, intent(INOUT) :: the(:)
      real, intent(INOUT) :: te(:)
      real, intent(INOUT) :: pe(:)
      real, intent(INOUT) :: thve(:)
      real, intent(INOUT) :: dne(:)
      real, intent(INOUT) :: upe(:)
      real, intent(INOUT) :: vpe(:)
      real, intent(INOUT) :: vel_e(:)
!!$    REAL,    INTENT(INOUT) :: zcon  (:)
      real, intent(INOUT) :: zt(:)
      real, intent(INOUT) :: thtcon(:)
      real, intent(INOUT) :: rvcon(:)
      real, intent(INOUT) :: picon(:)
      real, intent(INOUT) :: ucon(:)
      real, intent(INOUT) :: vcon(:)
      real, intent(INOUT) :: zm(:)
      real, intent(INOUT) :: dzm(:)
      real, intent(INOUT) :: dzt(:)
      ! rconstants
      real, intent(IN)    :: g
      real, intent(IN)    :: cp
      real, intent(IN)    :: cpor
      real, intent(IN)    :: p00
      real, intent(IN)    :: rgas

      character(LEN=*), parameter :: h = "**(get_env_condition)**"

      integer :: k, nk, kl
      real    :: znz
      logical :: found

      !!LFR:  As linhas abaixo (IF) foram comentadas pois set_grid calcula os valores de zt
      !!LFR:  zt é usado nas chamada de plumerise em intervalos de tempo regulares.
      !!LFR:  Entretanto na rotina de plumerise_driver, cabeça da estrutura, zt é uma
      !!LFR:  variável não preservada entre chamadas. Com isso, na primeira vez td 
      !!LFR:  mas nas próximas zt contém lixo e afeta a interpolação.

      !LFR IF( n_setgrid == 0) THEN
      !LFR   n_setgrid = 1
      ! define vertical grid of plume model
      call set_grid(zt, zm, dzm, dzt, dz, zsurf, nkp)
      ! zt(k) =  thermo and water levels
      ! zm(k) =  dynamical levels
      !LFR  ENDIF
      znz = zcon(k2)

      found = .false.
      do k = nkp, 1, -1
         if (zt(k) .lt. znz) then
            found = .true.
            exit
         end if
      end do

      if (.not. found) then
         stop ' ERROR: Envir stop 12 - chem_plumerise_scalar'
!!$       CALL FatalError(h//" Envir stop 12")
      end if

!--(DMK-BRAMS-5-INI)------------------------------------------------------
      kmt = min(k, nkp - 1)
!--(DMK-BRAMS-5-OLD)------------------------------------------------------
!    kmt=k
!--(DMK-BRAMS-5-FIM)------------------------------------------------------

      nk = k2 - k1 + 1

      if (nk > nkp) then
         stop ' ERROR: nk > nkp - chem_plumerise_scalar'
      end if

!!$    CALL htint(nk,  ucon,zcon(k1:nk+k1-1),kmt,upe  ,zt)
!!$    CALL htint(nk,  vcon,zcon(k1:nk+k1-1),kmt,vpe  ,zt)
!!$    CALL htint(nk,thtcon,zcon(k1:nk+k1-1),kmt,the  ,zt)
!!$    CALL htint(nk, rvcon,zcon(k1:nk+k1-1),kmt,qvenv,zt)
      !print *,'P1 - chem_plumerise_scalar ucon'
      !print *,"nk  =",nk,"kmt=",kmt
      !print *,"ucon=",(ucon(kl),kl=1,nkp)
      !print *,"zcon=",(zcon(kl),kl=1,nkp)
      !print *,"upe =",(upe(kl),kl=1,kmt)
      !print *,"zt  =",(zt(kl),kl=1,kmt)
      !print *,'------------------------------------------------------------------'
      call htint(nk, ucon, zcon, kmt, upe, zt)
      !print *,'P1 - chem_plumerise_scalar vcon'
      call htint(nk, vcon, zcon, kmt, vpe, zt)
      !print *,'P1 - chem_plumerise_scalar thtcon'
      call htint(nk, thtcon, zcon, kmt, the, zt)
      !print *,'P1 - chem_plumerise_scalar rvcon'
      call htint(nk, rvcon, zcon, kmt, qvenv, zt)
      do k = 1, kmt
         qvenv(k) = max(qvenv(k), 1e-8)
      end do

      pke(1) = picon(1)
      do k = 1, kmt
         thve(k) = the(k)*(1.+.61*qvenv(k)) ! virtual pot temperature
      end do
      do k = 2, kmt
         pke(k) = pke(k - 1) - g*2.*(zt(k) - zt(k - 1)) & ! exner function
                  /(thve(k) + thve(k - 1))
      end do
      do k = 1, kmt
         te(k) = the(k)*pke(k)/cp         ! temperature (K)
         pe(k) = (pke(k)/cp)**cpor*p00    ! pressure (Pa)
         dne(k) = pe(k)/(rgas*te(k)*(1.+.61*qvenv(k))) !  dry air density (kg/m3)

         vel_e(k) = sqrt(upe(k)**2 + vpe(k)**2)         !-env wind (m/s)
         !print*,'k,vel_e(k),te(k)=',vel_e(k),te(k)
      end do

      !-ewe - env wind effect
      if (wind_eff < 1) vel_e(1:kmt) = 0.

      !-use este para gerar o RAMS.out
      ! ------- print environment state
      !print*,'k,zt(k),pe(k),te(k)-273.15,qvenv(k)*1000'
      !do k=1,kmt
      ! write(*,100)  k,zt(k),pe(k),te(k)-273.15,qvenv(k)*1000.
      !enddo
      !stop 333

      !--------- nao eh necessario este calculo
      !do k=1,kmt
      !  call thetae(pe(k),te(k),qvenv(k),thee(k))
      !enddo

      !--------- converte press de Pa para kPa para uso modelo de plumerise

      do k = 1, kmt
         pe(k) = pe(k)*1.e-3
      end do

      return

      !para testes    ----------------------
      !
      !      open(13,file='prn.out')
      !      do 112, i=1,nkp
    !!  112 read (13,1000) zt(i),pe(i),te(i),rhe(i),dne(i),qvenv(i)
      !  112 read (13,1000) dummy,pe(i),te(i),rhe(i),dne(i),qvenv(i)
      !1000 format(1x,6F15.4)
      !      close(13)
      !
      !te(:)=te(:)+273.15
      !--------------------------------------

   end subroutine get_env_condition
   !-------------------------------------------------------------------------

   subroutine set_grid(zt, zm, dzm, dzt, dz, zsurf, nkp)

      ! plumegen_coms
      real, intent(INOUT) :: dz ! (DMK) alterado (OUT) para (INOUT)
      real, intent(IN)    :: zsurf
      real, intent(INOUT) :: zt(:)
      real, intent(INOUT) :: zm(:)
      real, intent(INOUT) :: dzm(:)
      real, intent(INOUT) :: dzt(:)
      integer, intent(IN) :: nkp
      integer :: k, mzp, kk

      dz = 100. ! set constant grid spacing of plume grid model(meters)

      mzp = nkp
      zt(1) = zsurf
      zm(1) = zsurf
      zt(2) = zt(1) + 0.5*dz
      zm(2) = zm(1) + dz
      do k = 3, mzp
         zt(k) = zt(k - 1) + dz ! thermo and water levels
         zm(k) = zm(k - 1) + dz ! dynamical levels
      end do
      !print*,zsurf
      !Print*,zt(:)
      do k = 1, mzp - 1
         dzm(k) = 1./(zt(k + 1) - zt(k))
      end do
      dzm(mzp) = dzm(mzp - 1)

      do k = 2, mzp
         dzt(k) = 1./(zm(k) - zm(k - 1))
      end do
      dzt(1) = dzt(2)*dzt(2)/dzt(3)

      !   dzm(1) = 0.5/dz
      !   dzm(2:mzp) = 1./dz
      print *, 'Zt output=', (zt(kk), kk=1, mzp)
      return
   end subroutine set_grid
   !-------------------------------------------------------------------------

   subroutine set_flam_vert(ztopmax, k1, k2, zzcon, W_VMD, VMD, nkp)

!!$    REAL,    INTENT(IN)  :: ztopmax(2)
      real, intent(IN)  :: ztopmax(:)
      integer, intent(OUT) :: k1
      integer, intent(OUT) :: k2
      ! plumegen_coms
!!$    REAL,    INTENT(IN)  :: zzcon(nkp)
      real, intent(IN)  :: zzcon(:)
      !- version 2
!!$    REAL,    INTENT(IN)  :: W_VMD(nkp,2)
!!$    REAL,    INTENT(OUT) :: VMD(nkp,2)
      real, intent(IN)  :: W_VMD(:, :)
      real, intent(OUT) :: VMD(:, :)
      integer, intent(IN)  :: nkp

      integer :: imm, k
      integer :: k_lim(2)
      real    :: w_thresold, xxx
      integer :: k_initial, k_final, ko, kk4, kl

      !- version 1
      do imm = 1, 2
         ! checar
         !    do k=1,m1-1
         do k = 1, nkp - 1
            if (zzcon(k) > ztopmax(imm)) exit
         end do
         k_lim(imm) = k
      end do
      k1 = max(3, k_lim(1))
      k2 = max(3, k_lim(2))

      if (k2 < k1) then
         !print*,'1: ztopmax k=',ztopmax(1), k1
         !print*,'2: ztopmax k=',ztopmax(2), k2
         k2 = k1
         !stop 1234
      end if

      !- version 2
      !- vertical mass distribution
      !-
      w_thresold = 1.
      do imm = 1, 2

         VMD(1:nkp, imm) = 0.
         xxx = 0.
         k_initial = 0
         k_final = 0

         !- define range of the upper detrainemnt layer
         do ko = nkp - 10, 2, -1

            if (w_vmd(ko, imm) < w_thresold) cycle

            if (k_final == 0) k_final = ko

            if (w_vmd(ko, imm) - 1. > w_vmd(ko - 1, imm)) then
               k_initial = ko
               exit
            end if

         end do
         !- if there is a non zero depth layer, make the mass vertical
         !  distribution
         if (k_final > 0 .and. k_initial > 0) then

            k_initial = int((k_final + k_initial)*0.5)

            !- parabolic vertical distribution between k_initial and k_final
            kk4 = k_final - k_initial + 2
            do ko = 1, kk4 - 1
               kl = ko + k_initial - 1
               VMD(kl, imm) = 6.*float(ko)/float(kk4)**2*(1.-float(ko)/float(kk4))
            end do
            if (sum(VMD(1:NKP, imm)) .ne. 1.) then
               xxx = (1.-sum(VMD(1:NKP, imm)))/float(k_final - k_initial + 1)
               do ko = k_initial, k_final
                  VMD(ko, imm) = VMD(ko, imm) + xxx !- values between 0 and 1.
               end do
               ! print*,'new mass=',sum(mass)*100.,xxx
               !pause
            end if
         end if !k_final > 0 .and. k_initial >

      end do

   end subroutine set_flam_vert

   !-------------------------------------------------------------------------

   subroutine get_fire_properties(imm, iveg_ag, burnt_area, use_last, &
                                  heating, area, maxtime, alpha, rsurf, fmoist, tdur, &
                                  plumerise_flag, FRP, ntime)

      integer, intent(IN)    :: ntime
      integer, intent(IN)    :: imm
      integer, intent(IN)    :: iveg_ag
      real, intent(IN)    :: burnt_area, FRP
      integer, intent(IN)    :: use_last
      integer, intent(IN)    :: plumerise_flag
      ! plumegen_coms
      real, intent(INOUT) :: heating(:)
      real, intent(INOUT) :: area
      integer, intent(INOUT) :: maxtime
      real, intent(OUT)   :: alpha
      real, intent(OUT)   :: rsurf
      real, intent(OUT)   :: fmoist
      real, intent(OUT)   :: tdur

      character(LEN=*), parameter :: h = "**(get_fire_properties)**"

      integer :: moist, i, icount
      real    :: bfract, effload, heat, hinc, heat_fluxW
      real    :: heat_flux(2, 4)

      integer :: mdur           ! (DMK) scratch, not used
      real    :: bload          ! (DMK) scratch, not used

      data heat_flux/ &
         !---------------------------------------------------------------------
         !  heat flux      !IGBP Land Cover            !
         ! min  ! max      !Legend and                    ! reference
         !    kW/m^2       !description              !
         !--------------------------------------------------------------------
         30.0, 80.0, &! Tropical Forest         ! igbp 2 & 4
         30.0, 80.0, &! Boreal forest           ! igbp 1 & 3
         4.4, 23.0, &! cerrado/woody savanna   | igbp  5 thru 9
         3.3, 3.3/! Grassland/cropland      ! igbp 10 thru 17
      !--------------------------------------------------------------------
      real, parameter :: beta = 0.88  !ref.: Paugam et al., 2015
      !real, parameter :: beta = 5.0   !ref.: Wooster et al., 2005

      !-- fire at the surface
      !area = 20.e+4   ! area of burn, m^2
      area = burnt_area! area of burn, m^2

      if (PLUMERISE_flag == 1) then
         !fluxo de calor para o bioma
         heat_fluxW = heat_flux(imm, iveg_ag)*1000. ! converte para W/m^2

      elseif (PLUMERISE_flag == 2) then
         ! "beta" factor converts FRP to convective energy
         heat_fluxW = beta*(FRP/area)/0.55 ! in W/m^2

      end if

      mdur = 53        ! duration of burn, minutes
      bload = 10.      ! total loading, kg/m**2
      moist = 10       ! fuel moisture, %. average fuel moisture,percent dry
      maxtime = mdur + 2  ! model time, min
      !maxtime =mdur-1  ! model time, min

      !heat = 21.e6    !- joules per kg of fuel consumed
      !heat = 15.5e6   !joules/kg - cerrado
      heat = 19.3e6    !joules/kg - floresta em alta floresta (mt)
      !alpha = 0.1      !- entrainment constant
      alpha = 0.05      !- entrainment constant

      !-------------------- printout ----------------------------------------

    !!WRITE ( * ,  * ) ' SURFACE =', ZSURF, 'M', '  LCL =', ZBASE, 'M'
      !
      !PRINT*,'======================================================='
      !print * , ' FIRE BOUNDARY CONDITION   :'
      !print * , ' DURATION OF BURN, MINUTES =',MDUR
      !print * , ' AREA OF BURN, HA              =',AREA*1.e-4
      !print * , ' HEAT FLUX, kW/m^2              =',heat_fluxW*1.e-3
      !print * , ' TOTAL LOADING, KG/M**2    =',BLOAD
      !print * , ' FUEL MOISTURE, %              =',MOIST !average fuel moisture,percent dry
      !print * , ' MODEL TIME, MIN.              =',MAXTIME
      !
      !
      !
      ! ******************** fix up inputs *********************************
      !

      !IF (MOD (MAXTIME, 2) .NE.0) MAXTIME = MAXTIME+1  !make maxtime even

      MAXTIME = MAXTIME*60  ! and put in seconds
      !
      RSURF = sqrt(AREA/3.14159) !- entrainment surface radius (m)

      FMOIST = MOIST/100.       !- fuel moisture fraction
      !
      !
      ! calculate the energy flux and water content at lboundary.
      ! fills heating() on a minute basis. could ask for a file at this po
      ! in the program. whatever is input has to be adjusted to a one
      ! minute timescale.
      !

      do I = 1, ntime         !- make sure of energy release
         HEATING(I) = 0.0001  !- avoid possible divide by 0
      end do
      !
      TDUR = MDUR*60.       !- number of seconds in the burn

      bfract = 1.             !- combustion factor

      EFFLOAD = BLOAD*BFRACT  !- patchy burning

      !     spread the burning evenly over the interval
      !     except for the first few minutes for stability
      ICOUNT = 1
      !
      if (MDUR > NTIME) &
         stop 'Increase time duration (ntime) in min - see file "plumerise_mod.f90"'
!!$         CALL FatalError(h//" Increase time duration (ntime) in min - see file chem_plumerise_scalar.f90")

      do while (ICOUNT .le. MDUR)
         !  HEATING (ICOUNT) = HEAT * EFFLOAD / TDUR  ! W/m**2
         !  HEATING (ICOUNT) = 80000.  * 0.55         ! W/m**2

         HEATING(ICOUNT) = heat_fluxW*0.55     ! W/m**2 (0.55 converte para energia convectiva)
         ICOUNT = ICOUNT + 1
      end do
      !     ramp for 5 minutes
      if (use_last /= 1) then

         HINC = HEATING(1)/4.
         HEATING(1) = 0.1
         HEATING(2) = HINC
         HEATING(3) = 2.*HINC
         HEATING(4) = 3.*HINC
      else
         if (imm == 1) then
            HINC = HEATING(1)/4.
            HEATING(1) = 0.1
            HEATING(2) = HINC
            HEATING(3) = 2.*HINC
            HEATING(4) = 3.*HINC
         else
            stop "check units and use heat_fluxW"
            HINC = (HEATING(1) - heat_flux(imm - 1, iveg_ag)*1000.*0.55)/4.
            HEATING(1) = heat_flux(imm - 1, iveg_ag)*1000.*0.55 + 0.1
            HEATING(2) = HEATING(1) + HINC
            HEATING(3) = HEATING(2) + HINC
            HEATING(4) = HEATING(3) + HINC
         end if
      end if
      !srf-24jan2007 - end!<<<<<<<<<<<<<<<<<<<<02082007

      return
   end subroutine get_fire_properties

   !----------------------------------------------------------------------------
   !
   subroutine MAKEPLUME(i, j, kmt, ztopmax, ixx, imm, use_last, w, t, pe, qsat, &
                        est, rho, zt, te, txs, wc, dz, maxtime, zm, wt, dzm, tt, qv, qh, &
                        qi, qc, qvenv, radius, alpha, visc, cvi, vth, cvh, vti, dzt, &
                        sc, area, rsurf, tdur, heating, fmoist, vel_e, vel_p, rad_p, &
                        vel_t, rad_t, W_VMD, izprint, nkp, ntime)
      !
      ! *********************************************************************
      !
      !    EQUATION SOURCE--Kessler Met.Monograph No. 32 V.10 (K)
      !    Alan Weinstein, JAS V.27 pp 246-255. (W),
      !    Ogura and Takahashi, Monthly Weather Review V.99,pp895-911 (OT)
      !    Roger Pielke,Mesoscale Meteorological Modeling,Academic Press,1984
      !    Originally developed by: Don Latham (USFS)
      !
      !
      ! ************************ VARIABLE ID ********************************
      !
      !     DT=COMPUTING TIME INCREMENT (SEC)
      !     DZ=VERTICAL INCREMENT (M)
      !     LBASE=LEVEL ,CLOUD BASE
      !
      !     CONSTANTS:
      !       G = GRAVITATIONAL ACCELERATION 9.80796 (M/SEC/SEC).
      !       R = DRY AIR GAS CONSTANT (287.04E6 JOULE/KG/DEG K)
      !       CP = SPECIFIC HT. (1004 JOULE/KG/DEG K)
      !       HEATCOND = HEAT OF CONDENSATION (2.5E6 JOULE/KG)
      !       HEATFUS = HEAT OF FUSION (3.336E5 JOULE/KG)
      !       HEATSUBL = HEAT OF SUBLIMATION (2.83396E6 JOULE/KG)
      !       EPS = RATIO OF MOL.WT. OF WATER VAPOR TO THAT OF DRY AIR (0.622)
      !       DES = DIFFERENCE BETWEEN VAPOR PRESSURE OVER WATER AND ICE (MB)
      !       TFREEZE = FREEZING TEMPERATURE (K)
      !
      !
      !     PARCEL VALUES:
      !       T = TEMPERATURE (K)
      !       TXS = TEMPERATURE EXCESS (K)
      !       QH = HYDROMETEOR WATER CONTENT (G/G DRY AIR)
      !       QHI = HYDROMETEOR ICE CONTENT (G/G DRY AIR)
      !       QC = WATER CONTENT (G/G DRY AIR)
      !       QVAP = WATER VAPOR MIXING RATIO (G/G DRY AIR)
      !       QSAT = SATURATION MIXING RATIO (G/G DRY AIR)
      !       RHO = DRY AIR DENSITY (G/M**3) MASSES = RHO*Q'S IN G/M**3
      !       ES = SATURATION VAPOR PRESSURE (kPa)
      !
      !     ENVIRONMENT VALUES:
      !       TE = TEMPERATURE (K)
      !       PE = PRESSURE (kPa)
      !       QVENV = WATER VAPOR (G/G)
      !       RHE = RELATIVE HUMIDITY FRACTION (e/esat)
      !       DNE = dry air density (kg/m^3)
      !
      !     HEAT VALUES:
      !       HEATING = HEAT OUTPUT OF FIRE (WATTS/M**2)
      !       MDUR = DURATION OF BURN, MINUTES
      !
      !       W = VERTICAL VELOCITY (M/S)
      !       RADIUS=ENTRAINMENT RADIUS (FCN OF Z)
      !        RSURF = ENTRAINMENT RADIUS AT GROUND (SIMPLE PLUME, TURNER)
      !        ALPHA = ENTRAINMENT CONSTANT
      !       MAXTIME = TERMINATION TIME (MIN)
      !
      !
      !**********************************************************************
      !**********************************************************************

      use node_mod, only: &
         mynum

      integer, intent(IN) :: nkp
      integer, intent(IN) :: ntime
      integer, intent(IN)    :: kmt, i, j
      real, intent(INOUT) :: ztopmax  ! (DMK) alterado (OUT) para (INOUT)
      integer, intent(IN)    :: ixx
      integer, intent(IN)    :: imm
      integer, intent(IN)    :: use_last
      integer, intent(IN)    :: izprint
      ! plumegen_coms
      real, intent(IN)    :: dz
      integer, intent(IN)    :: maxtime
      real, intent(IN)    :: alpha
      real, intent(IN)    :: area
      real, intent(IN)    :: rsurf
      real, intent(IN)    :: tdur
      real, intent(IN)    :: fmoist
!!$    REAL,    INTENT(IN)    :: heating(ntime)
!!$    REAL,    INTENT(INOUT) :: w     (nkp)
!!$    REAL,    INTENT(INOUT) :: t     (nkp)
!!$    REAL,    INTENT(IN)    :: pe    (nkp)
!!$    REAL,    INTENT(INOUT) :: qsat  (nkp) ! (DMK) alterado (OUT) para (INOUT)
!!$    REAL,    INTENT(INOUT) :: est   (nkp) ! (DMK) alterado (OUT) para (INOUT)
!!$    REAL,    INTENT(INOUT) :: rho   (nkp) ! (DMK) alterado (OUT) para (INOUT)
!!$    REAL,    INTENT(IN)    :: zt    (nkp)
!!$    REAL,    INTENT(IN)    :: te    (nkp)
!!$    REAL,    INTENT(INOUT) :: txs   (nkp) ! (DMK) alterado (OUT) para (INOUT)
!!$    REAL,    INTENT(INOUT) :: wc    (nkp)
!!$    REAL,    INTENT(IN)    :: zm    (nkp)
!!$    REAL,    INTENT(INOUT) :: wt    (nkp)
!!$    REAL,    INTENT(IN)    :: dzm   (nkp)
!!$    REAL,    INTENT(INOUT) :: tt    (nkp)
!!$    REAL,    INTENT(INOUT) :: qv    (nkp)
!!$    REAL,    INTENT(INOUT) :: qh    (nkp)
!!$    REAL,    INTENT(INOUT) :: qi    (nkp)
!!$    REAL,    INTENT(INOUT) :: qc    (nkp)
!!$    REAL,    INTENT(IN)    :: qvenv (nkp)
!!$    REAL,    INTENT(INOUT) :: radius(nkp)
!!$    REAL,    INTENT(INOUT) :: visc  (nkp)
!!$    REAL,    INTENT(INOUT) :: cvi   (nkp)
!!$    REAL,    INTENT(INOUT) :: vth   (nkp)
!!$    REAL,    INTENT(INOUT) :: cvh   (nkp)
!!$    REAL,    INTENT(INOUT) :: vti   (nkp)
!!$    REAL,    INTENT(IN)    :: dzt   (nkp)
!!$    REAL,    INTENT(IN)    :: sc    (nkp)
!!$!srf-awe
!!$    REAL,    INTENT(IN)    :: vel_e (nkp)
!!$    REAL,    INTENT(INOUT) :: vel_p (nkp)
!!$    REAL,    INTENT(INOUT) :: rad_p (nkp)
!!$    REAL,    INTENT(INOUT) :: vel_t (nkp)
!!$    REAL,    INTENT(INOUT) :: rad_t (nkp)
!!$    REAL,    INTENT(OUT)   :: W_VMD (nkp,2)
!!$!srf-awe
      real, intent(IN)    :: heating(:)
      real, intent(INOUT) :: w(:)
      real, intent(INOUT) :: t(:)
      real, intent(IN)    :: pe(:)
      real, intent(INOUT) :: qsat(:) ! (DMK) alterado (OUT) para (INOUT)
      real, intent(INOUT) :: est(:) ! (DMK) alterado (OUT) para (INOUT)
      real, intent(INOUT) :: rho(:) ! (DMK) alterado (OUT) para (INOUT)
      real, intent(IN)    :: zt(:)
      real, intent(IN)    :: te(:)
      real, intent(INOUT) :: txs(:) ! (DMK) alterado (OUT) para (INOUT)
      real, intent(INOUT) :: wc(:)
      real, intent(IN)    :: zm(:)
      real, intent(INOUT) :: wt(:)
      real, intent(IN)    :: dzm(:)
      real, intent(INOUT) :: tt(:)
      real, intent(INOUT) :: qv(:)
      real, intent(INOUT) :: qh(:)
      real, intent(INOUT) :: qi(:)
      real, intent(INOUT) :: qc(:)
      real, intent(IN)    :: qvenv(:)
      real, intent(INOUT) :: radius(:)
      real, intent(INOUT) :: visc(:)
      real, intent(INOUT) :: cvi(:)
      real, intent(INOUT) :: vth(:)
      real, intent(INOUT) :: cvh(:)
      real, intent(INOUT) :: vti(:)
      real, intent(IN)    :: dzt(:)
      real, intent(IN)    :: sc(:)
!srf-awe
      real, intent(IN)    :: vel_e(:)
      real, intent(INOUT) :: vel_p(:)
      real, intent(INOUT) :: rad_p(:)
      real, intent(INOUT) :: vel_t(:)
      real, intent(INOUT) :: rad_t(:)
      real, intent(OUT)   :: W_VMD(:, :)
!srf-awe

      !logical :: endspace
      integer :: k, kk, kkmax, deltak, ilastprint, &
                 nrectotal, i_micro, n_sub_step, nrec

      real    :: vc, g, r, cp, eps, &
                 tmelt, heatsubl, heatfus, heatcond, tfreeze, &
                 wmax, rmaxtime, es, dt_save

      real, external :: esat

      real    :: DELZ_THRESOLD

      character(len=2) :: cixx

      !
      !
      ! ******************* SOME CONSTANTS **********************************
      !
      !      XNO=10.0E06 median volume diameter raindrop (K table 4)
      !      VC = 38.3/(XNO**.125) mean volume fallspeed eqn. (K)
      !
      parameter(vc=5.107387)
      parameter(g=9.80796, r=287.04, cp=1004., eps=0.622, tmelt=273.3)
      parameter(heatsubl=2.834e6, heatfus=3.34e5, heatcond=2.501e6)
      parameter(tfreeze=269.3)
      !

      real    :: viscosity     ! (DMK) scratch
      real    :: tstpf         ! (DMK) scratch
      integer :: mintime       ! (DMK) scratch
      real    :: ztop          ! (DMK) scratch
      real    :: dqsdz         ! (DMK) scratch
      real    :: time          ! (DMK) scratch
      integer :: nm1           ! (DMK) scratch
      integer :: l             ! (DMK) scratch
      real    :: adiabat       ! (DMK) scratch
      real    :: vhrel         ! (DMK) scratch
      real    :: wbar          ! (DMK) scratch
      real    :: dt            ! (DMK) scratch
      real    :: ztop_(ntime)  ! (DMK) scratch
      real    :: scr1(nkp)     ! (DMK) scratch
      real    :: qvt(nkp)     ! (DMK) scratch
      real    :: qct(nkp)     ! (DMK) scratch
      real    :: qht(nkp)     ! (DMK) scratch
      real    :: qit(nkp)     ! (DMK) scratch

      tstpf = 2.0     !- timestep factor
      viscosity = 100.!- viscosity constant (original value: 0.001)

      nrectotal = 150
      nrec = 0
      !
      !*************** PROBLEM SETUP AND INITIAL CONDITIONS *****************
      mintime = 1
      ztopmax = 0.
      ztop = 0.
      time = 0.
      dt = 1.
      wmax = 1.
      kkmax = 10
      deltaK = 20
      ilastprint = 0
      L = 1   ! L initialization
      wbar = 0.
      W_VMD(:, :) = 0.

      !--- initialization
      ! CALL INITIAL(kmt)
      DELZ_THRESOLD = 1.*dz
      if (imm == 1 .and. use_last == 1) &
         call INITIAL(kmt, txs, w, t, wc, wt, qv, vth, vti, qh, qi, qc, est, qsat, rho, &
                      radius, visc, te, qvenv, rsurf, alpha, zt, viscosity, pe, &
                      area, dt, l, wbar, dqsdz, cvi, tdur, heating, &
                      mintime, fmoist, time, VEL_P, rad_p, nkp)
      if (use_last /= 1) &
         call INITIAL(kmt, txs, w, t, wc, wt, qv, vth, vti, qh, qi, qc, est, qsat, rho, &
                      radius, visc, te, qvenv, rsurf, alpha, zt, viscosity, pe, &
                      area, dt, l, wbar, dqsdz, cvi, tdur, heating, &
                      mintime, fmoist, time, VEL_P, rad_p, nkp)
      !
      !--- initial print fields:
      if (izprint .ne. 0) then
         write (cixx(1:2), '(i2.2)') ixx
         open (2, file='debug.'//cixx//'.dat')
         open (19, file='plumegen.'//cixx//'.gra', &
               form='unformatted', access='direct', status='unknown', &
               recl=4*nrectotal)  !PC
         !     recl=1*nrectotal) !sx6 e tupay
         call printout(izprint, nrectotal, mintime, dt, time, ztop, &
                       pe, t, te, qv, qsat, qc, qh, &
                       qi, zt, w, vth, sc, qvenv, nrec)
         ilastprint = 2
      end if

      ! ******************* model evolution ******************************
      rmaxtime = float(maxtime)
      !
      do while (TIME .le. RMAXTIME)  !beginning of time loop

         !   do itime=1,120

         !-- set model top integration
         nm1 = min(kmt, kkmax + deltak)

         !-- set timestep
         !dt = (zm(2)-zm(1)) / (tstpf * wmax)
         dt = min(5., (zm(2) - zm(1))/(tstpf*wmax))
         if (dt < 0.25 .or. maxval(T(:)) < 0.) then
            ztop = 0.
            ztop_(mintime) = ztop
            ztopmax = 0.
            kkmax = 1
            print *, "PRM_UNST=", mynum, i, j, dt
            call flush (6)
            exit
         end if

         !-- elapsed time, sec
         time = time + dt
         !-- elapsed time, minutes
         mintime = 1 + int(time)/60
         wmax = 1.  !no zeroes allowed.
         !************************** BEGIN SPACE LOOP **************************

         !-- zerout all model tendencies
         call tend0_plumerise(nm1, wt, tt, qvt, qct, qht, qit, vel_t, rad_t)

         !-- bounday conditions (k=1)
         L = 1
         wbar = W(1)
         call lbound(qh, qi, qc, w, t, wc, vth, vti, txs, visc, rho, qv, est, qsat, pe, &
                     alpha, area, rsurf, te, viscosity, dt, qvenv, l, wbar, dqsdz, cvi, &
                     tdur, heating, mintime, fmoist, time, VEL_P, rad_p)

         !-- dynamics for the level k>1
         !-- W advection
         !   call vel_advectc_plumerise(NM1,WC,WT,DNE,DZM)
         call vel_advectc_plumerise(NM1, WC, WT, RHO, DZM)

         !-- scalars advection 1
         call scl_advectc_plumerise(NM1, scr1, dt, w, wc, rho, dzm, zt, zm, dzt, &
                                    t, tt, qv, qvt, qc, qct, qi, qit, qh, qht, vel_p, &
                                    vel_t, rad_p, rad_t, nkp)

         !-- scalars advection 2
         !call scl_advectc_plumerise2('SC',NM1)

         !-- scalars entrainment, adiabatic
         call scl_misc(NM1, wbar, w, adiabat, alpha, radius, tt, t, te, qvt, &
                       qv, qvenv, qct, qc, qht, qh, qit, qi, vel_e, vel_p, vel_t, &
                       rad_p, rad_t)

         !-- scalars dinamic entrainment
         call scl_dyn_entrain(NM1, wbar, w, adiabat, alpha, radius, tt, t, &
                              te, qvt, qv, qvenv, qct, qc, qht, qh, qit, qi, vel_e, &
                              vel_p, vel_t, rad_p, rad_t)

         !-- gravity wave damping using Rayleigh friction layer fot T
         call damp_grav_wave(1, nm1, deltak, dt, zt, zm, w, t, tt, te)

         !-- microphysics
         !   goto 101 ! bypass microphysics
         dt_save = dt
         n_sub_step = 3
         dt = dt/float(n_sub_step)

         do i_micro = 1, n_sub_step
            !-- sedim ?
            call fallpart(NM1, rho, vth, vhrel, w, cvh, vti, cvi, qh, qi, zm, qht, qit)
            !-- microphysics
            do L = 2, nm1 - 1
               WBAR = 0.5*(W(L) + W(L - 1))
               ES = 0.1*ESAT(T(L))            !BLOB SATURATION VAPOR PRESSURE, EM KPA
               QSAT(L) = (EPS*ES)/(PE(L) - ES)  !BLOB SATURATION LWC G/G DRY AIR
               EST(L) = ES
               RHO(L) = 3483.8*PE(L)/T(L) ! AIR PARCEL DENSITY , G/M**3
               !srf18jun2005
               !        IF (W(L) .ge. 0.) DQSDZ = (QSAT(L  ) - QSAT(L-1)) / (ZT(L  ) -ZT(L-1))
               !        IF (W(L) .lt. 0.) DQSDZ = (QSAT(L+1) - QSAT(L  )) / (ZT(L+1) -ZT(L  ))
               if (W(L) .ge. 0.) then
                  DQSDZ = (QSAT(L + 1) - QSAT(L - 1))/(ZT(L + 1) - ZT(L - 1))
               else
                  DQSDZ = (QSAT(L + 1) - QSAT(L - 1))/(ZT(L + 1) - ZT(L - 1))
               end if

               call waterbal(qc, l, qh, qi, qv, t, qsat, wbar, dqsdz, dt, rho, est, cvi)
            end do
         end do
         dt = dt_save
         !
101      continue
         !
         !-- W-viscosity for stability
         call visc_W(nm1, kmt, zt, visc, zm, w, t, qv, qh, qc, qi, wt, tt, qvt, qct, qht, &
                     qit, vel_p, vel_t, rad_p, rad_t)

         !-- update scalars
         call update_plumerise(nm1, 'S', wt, dt, tt, qvt, qct, qht, w, t, qv, qc, qh, &
                               qit, qi, vel_p, vel_t, rad_p, rad_t)

         call hadvance_plumerise(1, nm1, WC, W, mintime)

         !-- Buoyancy
         call buoyancy_plumerise(nm1, t, te, qv, qvenv, qh, qi, qc, wt, scr1)

         !-- Entrainment
         call entrainment(nm1, w, wt, radius, alpha, vel_p, vel_e)

         !-- update W
         call update_plumerise(nm1, 'W', wt, dt, tt, qvt, qct, qht, w, t, qv, qc, &
                               qh, qit, qi, vel_p, vel_t, rad_p, rad_t)

         call hadvance_plumerise(2, nm1, WC, W, mintime)

         !-- misc
         do k = 2, nm1
            !    pe esta em kpa  - esat do rams esta em mbar = 100 Pa = 0.1 kpa
            es = 0.1*esat(t(k)) !blob saturation vapor pressure, em kPa
            !    rotina do plumegen calcula em kPa
            !    es    = esat_pr (t(k))  !blob saturation vapor pressure, em kPa
            qsat(k) = (eps*es)/(pe(k) - es)  !blob saturation lwc g/g dry air
            est(k) = es
            txs(k) = t(k) - te(k)
            rho(k) = 3483.8*pe(k)/t(k) ! air parcel density , g/m**3
            ! no pressure diff with radius

            if ((abs(wc(k))) .gt. wmax) wmax = abs(wc(k)) ! keep wmax largest w
         end do

         ! Gravity wave damping using Rayleigh friction layer for W
         call damp_grav_wave(2, nm1, deltak, dt, zt, zm, w, t, tt, te)
         !---

         !- update radius
         do k = 2, nm1
            radius(k) = rad_p(k)
         end do

         !-- try to find the plume top (above surface height)
         kk = 1
         do while (w(kk) .gt. 1.)
            kk = kk + 1
            ztop = zm(kk)
            !print*,'W=',w (kk)
         end do
         !
         ztop_(mintime) = ztop
         ztopmax = max(ztop, ztopmax)
         kkmax = max(kk, kkmax)
         !print * ,'ztopmax=', mintime,'mn ',ztop_(mintime), ztopmax

         !
         ! if the solution is going to a stationary phase, exit
         if (mintime > 10) then
            !   if(mintime > 20) then
            !    if( abs(ztop_(mintime)-ztop_(mintime-10)) < DZ ) exit
            if (abs(ztop_(mintime) - ztop_(mintime - 10)) < DELZ_THRESOLD) then

               !- determine W parameter to determine the VMD
               do k = 2, nm1
                  W_VMD(k, imm) = w(k)
               end do
               exit ! finish the integration
            end if
         end if

         if (ilastprint == mintime .and. izprint .ne. 0) then
            call printout(izprint, nrectotal, mintime, dt, time, ztop, pe, t, te, &
                          qv, qsat, qc, qh, qi, zt, w, vth, sc, qvenv, nrec)
            ilastprint = mintime + 1
         end if

      end do   !do next timestep

      !print * ,' ztopmax=',ztopmax,'m',mintime,'mn '
      !print*,'======================================================='
      !
      !the last printout
      if (izprint .ne. 0) then
         call printout(izprint, nrectotal, mintime, dt, time, ztop, &
                       pe, t, te, qv, qsat, qc, qh, qi, zt, w, vth, sc, qvenv, &
                       nrec)
         close (2)
         close (19)
      end if

      return
   end subroutine MAKEPLUME

   !----------------------------------------------------------------------------
   !
   subroutine BURN(EFLUX, WATER, tdur, time, heating, mintime, dt, fmoist)
      !
      !- calculates the energy flux and water content at lboundary

      real, intent(OUT) :: EFLUX
      real, intent(OUT) :: WATER
      ! plumegen_coms
      real, intent(IN)  :: tdur
      real, intent(IN)  :: time
      integer, intent(IN)  :: mintime
      real, intent(IN)  :: dt
      real, intent(IN)  :: fmoist
!!$    REAL,    INTENT(IN)  :: heating(ntime)
      real, intent(IN)  :: heating(:)

      !real, parameter :: HEAT = 21.E6 !Joules/kg
      !real, parameter :: HEAT = 15.5E6 !Joules/kg - cerrado
      real, parameter :: HEAT = 19.3e6 !Joules/kg - floresta em Alta Floresta (MT)
      !
      ! The emission factor for water is 0.5. The water produced, in kg,
      ! is then  fuel mass*0.5 + (moist/100)*mass per square meter.
      ! The fire burns for DT out of TDUR seconds, the total amount of
      ! fuel burned is AREA*BLOAD*(DT/TDUR) kg. this amount of fuel is
      ! considered to be spread over area AREA and so the mass burned per
      ! unit area is BLOAD*(DT/TDUR), and the rate is BLOAD/TDUR.
      !
      if (TIME .gt. TDUR) then !is the burn over?
         EFLUX = 0.000001    !prevent a potential divide by zero
         WATER = 0.
         return
      else
         !
         EFLUX = HEATING(MINTIME)                          ! Watts/m**2
         !  WATER = EFLUX * (DT / HEAT) * (0.5 + FMOIST)       ! kg/m**2
         WATER = EFLUX*(DT/HEAT)*(0.5 + FMOIST)/0.55 ! kg/m**2
         WATER = WATER*1000.                              ! g/m**2
         !
         !        print*,'BURN:',time,EFLUX/1.e+9
      end if
      !
      return
   end subroutine BURN
   !----------------------------------------------------------------------------
   !
   subroutine lbound(qh, qi, qc, w, t, wc, vth, vti, txs, visc, rho, qv, est, qsat, pe, &
                     alpha, area, rsurf, te, viscosity, dt, qvenv, l, wbar, dqsdz, cvi, &
                     tdur, heating, mintime, fmoist, time, VEL_P, rad_p)

      !
      ! ********** BOUNDARY CONDITIONS AT ZSURF FOR PLUME AND CLOUD ********
      !
      ! source of equations: J.S. Turner Buoyancy Effects in Fluids
      !                      Cambridge U.P. 1973 p.172,
      !                      G.A. Briggs Plume Rise, USAtomic Energy Commissio
      !                      TID-25075, 1969, P.28
      !
      ! fundamentally a point source below ground. at surface, this produces
      ! a velocity w(1) and temperature T(1) which vary with time. There is
      ! also a water load which will first saturate, then remainder go into
      ! QC(1).
      ! EFLUX = energy flux at ground,watt/m**2 for the last DT
      !

      ! plumegen_coms
      real, intent(IN)    :: alpha
      real, intent(IN)    :: area
      real, intent(IN)    :: rsurf
      real, intent(IN)    :: viscosity
      real, intent(IN)    :: dt
      integer, intent(IN)    :: l
      real, intent(IN)    :: wbar
      real, intent(IN)    :: dqsdz
      real, intent(IN)    :: tdur
      integer, intent(IN)    :: mintime
      real, intent(IN)    :: fmoist
      real, intent(IN)    :: time
!!$    REAL,    INTENT(IN)    :: heating(ntime)
!!$    REAL,    INTENT(INOUT) :: qh   (nkp)
!!$    REAL,    INTENT(INOUT) :: qi   (nkp)
!!$    REAL,    INTENT(INOUT) :: qc   (nkp)
!!$    REAL,    INTENT(INOUT) :: w    (nkp) ! (DMK) alterado (OUT) para (INOUT)
!!$    REAL,    INTENT(INOUT) :: t    (nkp)
!!$    REAL,    INTENT(INOUT) :: wc   (nkp) ! (DMK) alterado (OUT) para (INOUT)
!!$    REAL,    INTENT(INOUT) :: vth  (nkp) ! (DMK) alterado (OUT) para (INOUT)
!!$    REAL,    INTENT(INOUT) :: vti  (nkp) ! (DMK) alterado (OUT) para (INOUT)
!!$    REAL,    INTENT(INOUT) :: txs  (nkp) ! (DMK) alterado (OUT) para (INOUT)
!!$    REAL,    INTENT(INOUT) :: visc (nkp) ! (DMK) alterado (OUT) para (INOUT)
!!$    REAL,    INTENT(INOUT) :: rho  (nkp) ! (DMK) alterado (OUT) para (INOUT)
!!$    REAL,    INTENT(INOUT) :: qv   (nkp)
!!$    REAL,    INTENT(INOUT) :: est  (nkp) ! (DMK) alterado (OUT) para (INOUT)
!!$    REAL,    INTENT(INOUT) :: qsat (nkp) ! (DMK) alterado (OUT) para (INOUT)
!!$    REAL,    INTENT(IN)    :: pe   (nkp)
!!$    REAL,    INTENT(IN)    :: te   (nkp)
!!$    REAL,    INTENT(IN)    :: qvenv(nkp)
!!$    REAL,    INTENT(IN)    :: cvi  (nkp)
!!$    REAL,    INTENT(INOUT) :: VEL_P(nkp)
!!$    REAL,    INTENT(INOUT) :: rad_p(nkp)
      real, intent(IN)    :: heating(:)
      real, intent(INOUT) :: qh(:)
      real, intent(INOUT) :: qi(:)
      real, intent(INOUT) :: qc(:)
      real, intent(INOUT) :: w(:) ! (DMK) alterado (OUT) para (INOUT)
      real, intent(INOUT) :: t(:)
      real, intent(INOUT) :: wc(:) ! (DMK) alterado (OUT) para (INOUT)
      real, intent(INOUT) :: vth(:) ! (DMK) alterado (OUT) para (INOUT)
      real, intent(INOUT) :: vti(:) ! (DMK) alterado (OUT) para (INOUT)
      real, intent(INOUT) :: txs(:) ! (DMK) alterado (OUT) para (INOUT)
      real, intent(INOUT) :: visc(:) ! (DMK) alterado (OUT) para (INOUT)
      real, intent(INOUT) :: rho(:) ! (DMK) alterado (OUT) para (INOUT)
      real, intent(INOUT) :: qv(:)
      real, intent(INOUT) :: est(:) ! (DMK) alterado (OUT) para (INOUT)
      real, intent(INOUT) :: qsat(:) ! (DMK) alterado (OUT) para (INOUT)
      real, intent(IN)    :: pe(:)
      real, intent(IN)    :: te(:)
      real, intent(IN)    :: qvenv(:)
      real, intent(IN)    :: cvi(:)
      real, intent(INOUT) :: VEL_P(:)
      real, intent(INOUT) :: rad_p(:)

      real, parameter :: g = 9.80796, r = 287.04, cp = 1004.6, eps = 0.622, tmelt = 273.3
      real, parameter :: tfreeze = 269.3, pi = 3.14159, e1 = 1./3., e2 = 5./3.
      real            :: es, eflux, water, pres, c1, c2, f, zv, denscor, xwater

      real, external :: esat

      !
      QH(1) = QH(2)   !soak up hydrometeors
      QI(1) = QI(2)
      QC(1) = 0.       !no cloud here
      !
      !
      call BURN(EFLUX, WATER, tdur, time, heating, mintime, dt, fmoist)
      !
      !  calculate parameters at boundary from a virtual buoyancy point source
      !
      PRES = PE(1)*1000.   !need pressure in N/m**2

      C1 = 5./(6.*ALPHA)  !alpha is entrainment constant

      C2 = 0.9*ALPHA

      F = EFLUX/(PRES*CP*PI)

      F = G*R*F*AREA  !buoyancy flux

      ZV = C1*RSURF  !virtual boundary height

      W(1) = C1*((C2*F)**E1)/ZV**E1  !boundary velocity

      DENSCOR = C1*F/G/(C2*F)**E1/ZV**E2   !density correction

      T(1) = TE(1)/(1.-DENSCOR)    !temperature of virtual plume at zsurf

      !
      WC(1) = W(1)

      VEL_P(1) = 0.
      rad_p(1) = rsurf
      !SC(1) = SCE(1)+F/1000.*dt  ! gas/particle (g/g)

      ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      !     match dw/dz,dt/dz at the boundary. F is conserved.
      !
      !WBAR = W (1) * (1. - 1. / (6. * ZV) )
      !ADVW = WBAR * W (1) / (3. * ZV)
      !ADVT = WBAR * (5. / (3. * ZV) ) * (DENSCOR / (1. - DENSCOR) )
      !ADVC = 0.
      !ADVH = 0.
      !ADVI = 0.
      !ADIABAT = - WBAR * G / CP
      VTH(1) = -4.
      VTI(1) = -3.
      TXS(1) = T(1) - TE(1)

      VISC(1) = VISCOSITY

      RHO(1) = 3483.8*PE(1)/T(1)   !air density at level 1, g/m**3

      XWATER = WATER/(W(1)*DT*RHO(1))   !firewater mixing ratio

      QV(1) = XWATER + QVENV(1)  !plus what's already there

      !  PE esta em kPa  - ESAT do RAMS esta em mbar = 100 Pa = 0.1 kPa
      ES = 0.1*ESAT(T(1)) !blob saturation vapor pressure, em kPa
      !  rotina do plumegen ja calcula em kPa
      !  ES       = ESAT_PR (T(1))  !blob saturation vapor pressure, em kPa

      EST(1) = ES
      QSAT(1) = (EPS*ES)/(PE(1) - ES)   !blob saturation lwc g/g dry air

      if (QV(1) .gt. QSAT(1)) then
         QC(1) = QV(1) - QSAT(1) + QC(1)  !remainder goes into cloud drops
         QV(1) = QSAT(1)
      end if
      !
      call WATERBAL(qc, l, qh, qi, qv, t, qsat, wbar, dqsdz, dt, rho, est, cvi)
      !
      return
   end subroutine LBOUND

   !----------------------------------------------------------------------------
   !
   subroutine INITIAL(kmt, txs, w, t, wc, wt, qv, vth, vti, qh, qi, qc, est, qsat, rho, &
                      radius, visc, te, qvenv, rsurf, alpha, zt, viscosity, pe, &
                      area, dt, l, wbar, dqsdz, cvi, tdur, heating, mintime, &
                      fmoist, time, VEL_P, rad_p, nkp)

      !
      ! ************* SETS UP INITIAL CONDITIONS FOR THE PROBLEM ************
      integer, intent(IN) :: nkp
      integer, intent(IN) :: kmt

      ! plumegen_coms
      real, intent(IN)    :: rsurf
      real, intent(IN)    :: alpha
      real, intent(IN)    :: viscosity
      real, intent(IN)    :: area
      real, intent(IN)    :: dt
      integer, intent(IN)    :: l
      real, intent(IN)    :: wbar
      real, intent(IN)    :: dqsdz
      real, intent(IN)    :: tdur
      integer, intent(IN)    :: mintime
      real, intent(IN)    :: fmoist
      real, intent(IN)    :: time
!!$    REAL,    INTENT(IN)    :: heating(ntime)
!!$    REAL,    INTENT(INOUT) :: txs   (nkp) ! (DMK) alterado (OUT) para (INOUT)
!!$    REAL,    INTENT(INOUT) :: w     (nkp) ! (DMK) alterado (OUT) para (INOUT)
!!$    REAL,    INTENT(INOUT) :: t     (nkp) ! (DMK) alterado (OUT) para (INOUT)
!!$    REAL,    INTENT(INOUT) :: wc    (nkp) ! (DMK) alterado (OUT) para (INOUT)
!!$    REAL,    INTENT(INOUT) :: wt    (nkp) ! (DMK) alterado (OUT) para (INOUT)
!!$    REAL,    INTENT(INOUT) :: qv    (nkp) ! (DMK) alterado (OUT) para (INOUT)
!!$    REAL,    INTENT(INOUT) :: vth   (nkp) ! (DMK) alterado (OUT) para (INOUT)
!!$    REAL,    INTENT(INOUT) :: vti   (nkp) ! (DMK) alterado (OUT) para (INOUT)
!!$    REAL,    INTENT(INOUT) :: qh    (nkp) ! (DMK) alterado (OUT) para (INOUT)
!!$    REAL,    INTENT(INOUT) :: qi    (nkp) ! (DMK) alterado (OUT) para (INOUT)
!!$    REAL,    INTENT(INOUT) :: qc    (nkp) ! (DMK) alterado (OUT) para (INOUT)
!!$    REAL,    INTENT(INOUT) :: est   (nkp) ! (DMK) alterado (OUT) para (INOUT)
!!$    REAL,    INTENT(INOUT) :: qsat  (nkp) ! (DMK) alterado (OUT) para (INOUT)
!!$    REAL,    INTENT(INOUT) :: rho   (nkp) ! (DMK) alterado (OUT) para (INOUT)
!!$    REAL,    INTENT(INOUT) :: radius(nkp)
!!$    REAL,    INTENT(INOUT) :: visc  (nkp) ! (DMK) alterado (OUT) para (INOUT)
!!$    REAL,    INTENT(IN)    :: te    (nkp)
!!$    REAL,    INTENT(IN)    :: qvenv (nkp)
!!$    REAL,    INTENT(IN)    :: zt    (nkp)
!!$    REAL,    INTENT(IN)    :: pe    (nkp)
!!$    REAL,    INTENT(IN)    :: cvi   (nkp)
!!$    REAL,    INTENT(INOUT) :: VEL_P (nkp)
!!$    REAL,    INTENT(INOUT) :: rad_p (nkp)
      real, intent(IN)    :: heating(:)
      real, intent(INOUT) :: txs(:) ! (DMK) alterado (OUT) para (INOUT)
      real, intent(INOUT) :: w(:) ! (DMK) alterado (OUT) para (INOUT)
      real, intent(INOUT) :: t(:) ! (DMK) alterado (OUT) para (INOUT)
      real, intent(INOUT) :: wc(:) ! (DMK) alterado (OUT) para (INOUT)
      real, intent(INOUT) :: wt(:) ! (DMK) alterado (OUT) para (INOUT)
      real, intent(INOUT) :: qv(:) ! (DMK) alterado (OUT) para (INOUT)
      real, intent(INOUT) :: vth(:) ! (DMK) alterado (OUT) para (INOUT)
      real, intent(INOUT) :: vti(:) ! (DMK) alterado (OUT) para (INOUT)
      real, intent(INOUT) :: qh(:) ! (DMK) alterado (OUT) para (INOUT)
      real, intent(INOUT) :: qi(:) ! (DMK) alterado (OUT) para (INOUT)
      real, intent(INOUT) :: qc(:) ! (DMK) alterado (OUT) para (INOUT)
      real, intent(INOUT) :: est(:) ! (DMK) alterado (OUT) para (INOUT)
      real, intent(INOUT) :: qsat(:) ! (DMK) alterado (OUT) para (INOUT)
      real, intent(INOUT) :: rho(:) ! (DMK) alterado (OUT) para (INOUT)
      real, intent(INOUT) :: radius(:)
      real, intent(INOUT) :: visc(:) ! (DMK) alterado (OUT) para (INOUT)
      real, intent(IN)    :: te(:)
      real, intent(IN)    :: qvenv(:)
      real, intent(IN)    :: zt(:)
      real, intent(IN)    :: pe(:)
      real, intent(IN)    :: cvi(:)
      real, intent(INOUT) :: VEL_P(:)
      real, intent(INOUT) :: rad_p(:)

      real, parameter :: tfreeze = 269.3
      integer         :: k
      real            :: es

      integer         :: n ! (DMK) scratch

      real, external :: esat

      !

      N = kmt
      ! initialize temperature structure,to the end of equal spaced sounding,
      do k = 1, N
         TXS(k) = 0.0
         W(k) = 0.0
         T(k) = TE(k)       !blob set to environment
         WC(k) = 0.0
         WT(k) = 0.0
         QV(k) = QVENV(k)   !blob set to environment
         VTH(k) = 0.                !initial rain velocity = 0
         VTI(k) = 0.                !initial ice  velocity = 0
         QH(k) = 0.                !no rain
         QI(k) = 0.                !no ice
         QC(k) = 0.                !no cloud drops
         !  PE esta em kPa  - ESAT do RAMS esta em mbar = 100 Pa = 0.1 kPa
         ES = 0.1*ESAT(T(k)) !blob saturation vapor pressure, em kPa
         !  rotina do plumegen calcula em kPa
         !  ES       = ESAT_PR (T(k))  !blob saturation vapor pressure, em kPa
         EST(k) = ES
         QSAT(k) = (.622*ES)/(PE(k) - ES) !saturation lwc g/g
         RHO(k) = 3483.8*PE(k)/T(k)         !dry air density g/m**3

         VEL_P(k) = 0.
         rad_p(k) = 0.

      end do

      ! Initialize the entrainment radius, Turner-style plume
      radius(1) = rsurf
      rad_p(1) = rsurf
      do k = 2, N
         radius(k) = radius(k - 1) + (6./5.)*alpha*(zt(k) - zt(k - 1))
         rad_p(k) = radius(k)
      end do

      !  Initialize the viscosity
      VISC(1) = VISCOSITY
      do k = 2, N
         !  VISC (k) = VISCOSITY!max(1.e-3,visc(k-1) - 1.* VISCOSITY/float(nkp))
         VISC(k) = max(1.e-3, visc(k - 1) - 1.*VISCOSITY/float(nkp))
      end do
      !--   Initialize gas/concentration
      !DO k =10,20
      !   SC(k) = 20.
      !ENDDO
      !stop 333

      call lbound(qh, qi, qc, w, t, wc, vth, vti, txs, visc, rho, qv, est, qsat, pe, &
                  alpha, area, rsurf, te, viscosity, dt, qvenv, l, wbar, dqsdz, cvi, &
                  tdur, heating, mintime, fmoist, time, VEL_P, rad_p)
      return
   end subroutine INITIAL

   !----------------------------------------------------------------------------
   !
   subroutine damp_grav_wave(ifrom, nm1, deltak, dt, zt, zm, w, t, tt, te)

      integer, intent(IN)    :: ifrom
      integer, intent(IN)    :: nm1
      integer, intent(IN)    :: deltak
      real, intent(IN)    :: dt
!!$    REAL,    INTENT(IN)    :: zt(nm1)
!!$    REAL,    INTENT(IN)    :: zm(nm1)
!!$    REAL,    INTENT(INOUT) :: w (nm1)
!!$    REAL,    INTENT(INOUT) :: t (nm1)
!!$    REAL,    INTENT(INOUT) :: tt(nm1)
!!$    REAL,    INTENT(IN)    :: te(nm1)
      real, intent(IN)    :: zt(:)
      real, intent(IN)    :: zm(:)
      real, intent(INOUT) :: w(:)
      real, intent(INOUT) :: t(:)
      real, intent(INOUT) :: tt(:)
      real, intent(IN)    :: te(:)

      real :: dummy(nm1)

      if (ifrom == 1) then
         call friction(ifrom, nm1, deltak, dt, zt, zm, t, tt, te)
         !call friction(ifrom,nm1,dt,zt,zm,qv,qvt,qvenv)
         ! call friction(ifrom,nm1,deltak,dt,zt,zm,vel_p,vel_t,vel_e)
         return
      end if

      dummy(:) = 0.
      if (ifrom == 2) call friction(ifrom, nm1, deltak, dt, zt, zm, w, dummy, dummy)
      !call friction(ifrom,nm1,dt,zt,zm,qi,qit ,dummy)
      !call friction(ifrom,nm1,dt,zt,zm,qh,qht ,dummy)
      !call friction(ifrom,nm1,dt,zt,zm,qc,qct ,dummy)
      return
   end subroutine damp_grav_wave

   !----------------------------------------------------------------------------
   !
   subroutine friction(ifrom, nm1, deltak, dt, zt, zm, var1, vart, var2)

      integer, intent(IN)    :: ifrom
      integer, intent(IN)    :: nm1
      integer, intent(IN)    :: deltak
      real, intent(IN)    :: dt
      real, intent(IN)    :: zt(nm1)
      real, intent(IN)    :: zm(nm1)
      real, intent(INOUT) :: var1(nm1)
      real, intent(INOUT) :: vart(nm1)
      real, intent(IN)    :: var2(nm1)

      integer :: k, kf
      real    :: zmkf, ztop, distim, c1, c2

      !nfpt=50
      !kf = nm1 - nfpt
      !kf = nm1 - INT(deltak/2)
      kf = nm1 - int(deltak)

      zmkf = zm(kf) !old: float(kf )*dz
      ztop = zm(nm1)
      !distim = min(4.*dt,200.)
      !distim = 60. ! orig
      distim = min(3.*dt, 60.)

      c1 = 1./(distim*(ztop - zmkf))
      c2 = dt*c1

      if (ifrom == 1) then
         do k = nm1, 2, -1
            if (zt(k) .le. zmkf) cycle
            vart(k) = vart(k) + c1*(zt(k) - zmkf)*(var2(k) - var1(k))
         end do
      elseif (ifrom == 2) then
         do k = nm1, 2, -1
            if (zt(k) .le. zmkf) cycle
            var1(k) = var1(k) + c2*(zt(k) - zmkf)*(var2(k) - var1(k))
         end do
      end if
      return
   end subroutine friction

   !----------------------------------------------------------------------------
   !
   subroutine vel_advectc_plumerise(m1, wc, wt, rho, dzm)

      integer, intent(IN)    :: m1
!!$    REAL,    INTENT(IN)    :: wc (m1)
!!$    REAL,    INTENT(INOUT) :: wt (m1)
!!$    REAL,    INTENT(IN)    :: rho(m1)
!!$    REAL,    INTENT(IN)    :: dzm(m1)
      real, intent(IN)    :: wc(:)
      real, intent(INOUT) :: wt(:)
      real, intent(IN)    :: rho(:)
      real, intent(IN)    :: dzm(:)

      integer :: k
      real    :: flxw(m1)
      real    :: dn0(m1) ! var local
      real    :: c1z

      !dzm(:)= 1./dz

      dn0(1:m1) = rho(1:m1)*1.e-3 ! converte de cgs para mks

      flxw(1) = wc(1)*dn0(1)

      do k = 2, m1 - 1
         flxw(k) = wc(k)*.5*(dn0(k) + dn0(k + 1))
      end do

      ! Compute advection contribution to W tendency

      c1z = .5

      do k = 2, m1 - 2

         wt(k) = wt(k) &
                 + c1z*dzm(k)/(dn0(k) + dn0(k + 1))*( &
                 (flxw(k) + flxw(k - 1))*(wc(k) + wc(k - 1)) &
                 - (flxw(k) + flxw(k + 1))*(wc(k) + wc(k + 1)) &
                 + (flxw(k + 1) - flxw(k - 1))*2.*wc(k))

      end do

      return
   end subroutine vel_advectc_plumerise

   !----------------------------------------------------------------------------
   !
   subroutine hadvance_plumerise(iac, m1, wc, wp, mintime)

      integer, intent(IN)    :: iac
      integer, intent(IN)    :: m1
!!$    REAL,    INTENT(INOUT) :: wc(m1)
!!$    REAL,    INTENT(INOUT) :: wp(m1)
      real, intent(INOUT) :: wc(:)
      real, intent(INOUT) :: wp(:)
      integer, intent(IN)    :: mintime

      real :: dummy(m1)
      real :: eps

      !     It is here that the Asselin filter is applied.  For the velocities
      !     and pressure, this must be done in two stages, the first when
      !     IAC=1 and the second when IAC=2.

      eps = .2
      if (mintime == 1) eps = 0.5

      !     For both IAC=1 and IAC=2, call PREDICT for U, V, W, and P.
      !
      call predict_plumerise(m1, wc, wp, dummy, iac, eps)
      !print*,'mintime',mintime,eps
      !do k=1,m1
      !   print*,'W-HAD',k,wc(k),wp(k),wt(k)
      !enddo
      return
   end subroutine hadvance_plumerise

   !----------------------------------------------------------------------------
   !
   subroutine predict_plumerise(npts, ac, ap, af, iac, epsu)

      integer, intent(IN)    :: npts
!!$    REAL,    INTENT(INOUT) :: ac(npts)
!!$    REAL,    INTENT(INOUT) :: ap(npts)
!!$    REAL,    INTENT(INOUT) :: af(npts) ! (DMK) alterado (OUT) para (INOUT)
      real, intent(INOUT) :: ac(:)
      real, intent(INOUT) :: ap(:)
      real, intent(INOUT) :: af(:) ! (DMK) alterado (OUT) para (INOUT)
      integer, intent(IN)    :: iac
      real, intent(IN)    :: epsu

      integer :: m

      !     For IAC=3, this routine moves the arrays AC and AP forward by
      !     1 time level by adding in the prescribed tendency. It also
      !     applies the Asselin filter given by:

      !              {AC} = AC + EPS * (AP - 2 * AC + AF)

      !     where AP,AC,AF are the past, current and future time levels of A.
      !     All IAC=1 does is to perform the {AC} calculation without the AF
      !     term present.  IAC=2 completes the calculation of {AC} by adding
      !     the AF term only, and advances AC by filling it with input AP
      !     values which were already updated in ACOUSTC.
      !

      if (iac .eq. 1) then
         do m = 1, npts
            ac(m) = ac(m) + epsu*(ap(m) - 2.*ac(m))
         end do
         return
      elseif (iac .eq. 2) then
         do m = 1, npts
            af(m) = ap(m)
            ap(m) = ac(m) + epsu*af(m)
         end do
         !elseif (iac .eq. 3) then
         !   do m = 1,npts
         !      af(m) = ap(m) + dtlp * fa(m)
         !   enddo
         !   do m = 1,npts
         !      ap(m) = ac(m) + epsu * (ap(m) - 2. * ac(m) + af(m))
         !   enddo
      end if

      do m = 1, npts
         ac(m) = af(m)
      end do
      return
   end subroutine predict_plumerise

   !----------------------------------------------------------------------------
   !
   subroutine buoyancy_plumerise(m1, T, TE, QV, QVENV, QH, QI, QC, WT, scr1)

      integer, intent(IN)    :: m1
!!$    REAL,    INTENT(IN)    :: T    (m1)
!!$    REAL,    INTENT(IN)    :: TE   (m1)
!!$    REAL,    INTENT(IN)    :: QV   (m1)
!!$    REAL,    INTENT(IN)    :: QVENV(m1)
!!$    REAL,    INTENT(IN)    :: QH   (m1)
!!$    REAL,    INTENT(IN)    :: QI   (m1)
!!$    REAL,    INTENT(IN)    :: QC   (m1)
!!$    REAL,    INTENT(INOUT) :: wt   (m1)
!!$    REAL,    INTENT(INOUT) :: scr1 (m1) ! (DMK) alterado (OUT) para (INOUT)
      real, intent(IN)    :: T(:)
      real, intent(IN)    :: TE(:)
      real, intent(IN)    :: QV(:)
      real, intent(IN)    :: QVENV(:)
      real, intent(IN)    :: QH(:)
      real, intent(IN)    :: QI(:)
      real, intent(IN)    :: QC(:)
      real, intent(INOUT) :: wt(:)
      real, intent(INOUT) :: scr1(:) ! (DMK) alterado (OUT) para (INOUT)

      real, parameter :: g = 9.8, eps = 0.622, gama = 0.5 ! mass virtual coeff.
      real, parameter :: mu = 0.15

      integer         :: k
      real            :: TV, TVE, QWTOTL, umgamai

      !- orig
      umgamai = 1./(1.+gama) ! compensa a falta do termo de aceleracao associado `as
      ! pertubacoes nao-hidrostaticas no campo de pressao

      !- new                 ! Siesbema et al, 2004
      !umgamai = 1./(1.-2.*mu)

      do k = 2, m1 - 1

         TV = T(k)*(1.+(QV(k)/EPS))/(1.+QV(k))  !blob virtual temp.
         TVE = TE(k)*(1.+(QVENV(k)/EPS))/(1.+QVENV(k))  !and environment

         QWTOTL = QH(k) + QI(k) + QC(k)                       ! QWTOTL*G is drag
         !- orig
         !scr1(k)= G*( umgamai*(  TV - TVE) / TVE   - QWTOTL)
         scr1(k) = G*umgamai*((TV - TVE)/TVE - QWTOTL)

         !if(k .lt. 10)print*,'BT',k,TV,TVE,TVE,QWTOTL
      end do

      do k = 2, m1 - 2
         wt(k) = wt(k) + 0.5*(scr1(k) + scr1(k + 1))
         !   print*,'W-BUO',k,wt(k),scr1(k),scr1(k+1)
      end do

   end subroutine buoyancy_plumerise

   !----------------------------------------------------------------------------
   !
   subroutine ENTRAINMENT(m1, w, wt, radius, ALPHA, vel_p, vel_e)

      integer, intent(IN)    :: m1
      real, intent(IN)    :: alpha
!!$    REAL,    INTENT(IN)    :: w     (m1)
!!$    REAL,    INTENT(INOUT) :: wt    (m1)
!!$    REAL,    INTENT(IN)    :: radius(m1)
!!$    REAL,    INTENT(IN)    :: vel_p (m1)
!!$    REAL,    INTENT(IN)    :: vel_e (m1)
      real, intent(IN)    :: w(:)
      real, intent(INOUT) :: wt(:)
      real, intent(IN)    :: radius(:)
      real, intent(IN)    :: vel_p(:)
      real, intent(IN)    :: vel_e(:)

      real, parameter :: mu = 0.15, gama = 0.5 ! mass virtual coeff.

      integer         :: k
      real            :: DMDTM, WBAR, RADIUS_BAR, umgamai, DYN_ENTR

      !- new - Siesbema et al, 2004
      !umgamai = 1./(1.-2.*mu)

      !- orig
      !umgamai = 1
      umgamai = 1./(1.+gama) ! compensa a falta do termo de aceleracao associado `as
      ! pertubacoes nao-hidrostaticas no campo de pressao

      !
      !-- ALPHA/RADIUS(L) = (1/M)DM/DZ  (W 14a)
      do k = 2, m1 - 1

         !-- for W: WBAR is only W(k)
         WBAR = W(k)
         RADIUS_BAR = 0.5*(RADIUS(k) + RADIUS(k - 1))
         ! orig
         !DMDTM =           2. * ALPHA * ABS (WBAR) / RADIUS_BAR  != (1/M)DM/DT
         DMDTM = umgamai*2.*ALPHA*abs(WBAR)/RADIUS_BAR  != (1/M)DM/DT

         !--  DMDTM*W(L) entrainment,
         wt(k) = wt(k) - DMDTM*abs(WBAR)

         !if(VEL_P (k) - VEL_E (k) > 0.) cycle

         !-   dynamic entrainment
         DYN_ENTR = (2./3.1416)*0.5*abs(VEL_P(k) - VEL_E(k) + VEL_P(k - 1) - &
                                        VEL_E(k - 1))/RADIUS_BAR

         wt(k) = wt(k) - DYN_ENTR*abs(WBAR)

         !- entraiment acceleration for output only
         !dwdt_entr(k) =  - DMDTM*ABS (WBAR)- DYN_ENTR*ABS (WBAR)

      end do
   end subroutine ENTRAINMENT

   !----------------------------------------------------------------------------
   !
   subroutine scl_advectc_plumerise(mzp, scr1, dt, w, wc, rho, dzm, zt, zm, dzt, &
                                    t, tt, qv, qvt, qc, qct, qi, qit, qh, qht, vel_p, &
                                    vel_t, rad_p, rad_t, nkp)

      integer, intent(IN) :: nkp
      integer, intent(IN) :: mzp
      ! plumegen_coms
      real, intent(IN)    :: dt
!!$    REAL, INTENT(INOUT) :: scr1 (nkp)
!!$    REAL, INTENT(IN)    :: w    (nkp)
!!$    REAL, INTENT(IN)    :: wc   (nkp)
!!$    REAL, INTENT(IN)    :: rho  (nkp)
!!$    REAL, INTENT(IN)    :: dzm  (nkp)
!!$    REAL, INTENT(IN)    :: zt   (nkp)
!!$    REAL, INTENT(IN)    :: zm   (nkp)
!!$    REAL, INTENT(IN)    :: dzt  (nkp)
!!$    REAL, INTENT(IN)    :: t    (nkp)
!!$    REAL, INTENT(INOUT) :: tt   (nkp)
!!$    REAL, INTENT(IN)    :: qv   (nkp)
!!$    REAL, INTENT(INOUT) :: qvt  (nkp)
!!$    REAL, INTENT(IN)    :: qc   (nkp)
!!$    REAL, INTENT(INOUT) :: qct  (nkp)
!!$    REAL, INTENT(IN)    :: qi   (nkp)
!!$    REAL, INTENT(INOUT) :: qit  (nkp)
!!$    REAL, INTENT(IN)    :: qh   (nkp)
!!$    REAL, INTENT(INOUT) :: qht  (nkp)
!!$    REAL, INTENT(IN)    :: vel_p(nkp)
!!$    REAL, INTENT(INOUT) :: vel_t(nkp)
!!$    REAL, INTENT(IN)    :: rad_p(nkp)
!!$    REAL, INTENT(INOUT) :: rad_t(nkp)
      real, intent(INOUT) :: scr1(:)
      real, intent(IN)    :: w(:)
      real, intent(IN)    :: wc(:)
      real, intent(IN)    :: rho(:)
      real, intent(IN)    :: dzm(:)
      real, intent(IN)    :: zt(:)
      real, intent(IN)    :: zm(:)
      real, intent(IN)    :: dzt(:)
      real, intent(IN)    :: t(:)
      real, intent(INOUT) :: tt(:)
      real, intent(IN)    :: qv(:)
      real, intent(INOUT) :: qvt(:)
      real, intent(IN)    :: qc(:)
      real, intent(INOUT) :: qct(:)
      real, intent(IN)    :: qi(:)
      real, intent(INOUT) :: qit(:)
      real, intent(IN)    :: qh(:)
      real, intent(INOUT) :: qht(:)
      real, intent(IN)    :: vel_p(:)
      real, intent(INOUT) :: vel_t(:)
      real, intent(IN)    :: rad_p(:)
      real, intent(INOUT) :: rad_t(:)

      real    :: dtlto2
      integer :: k
      real    :: vt3dc(nkp) ! (DMK) scratch
      real    :: vt3df(nkp) ! (DMK) scratch
      real    :: vt3dg(nkp) ! (DMK) scratch
      real    :: vt3dk(nkp) ! (DMK) scratch
      real    :: vctr1(nkp) ! (DMK) scratch
      real    :: vctr2(nkp) ! (DMK) scratch

      !  wp => w
      !- Advect  scalars
      dtlto2 = .5*dt
      !  vt3dc(1) =      (w(1) + wc(1)) * dtlto2 * dne(1)
      vt3dc(1) = (w(1) + wc(1))*dtlto2*rho(1)*1.e-3!converte de CGS p/ MKS
      vt3df(1) = .5*(w(1) + wc(1))*dtlto2*dzm(1)

      do k = 2, mzp
         !     vt3dc(k) =  (w(k) + wc(k)) * dtlto2 *.5 * (dne(k) + dne(k+1))
         vt3dc(k) = (w(k) + wc(k))*dtlto2*.5*(rho(k) + rho(k + 1))*1.e-3
         vt3df(k) = (w(k) + wc(k))*dtlto2*.5*dzm(k)
         !print*,'vt3df-vt3dc',k,vt3dc(k),vt3df(k)
      end do

      !-srf-24082005
      !  do k = 1,mzp-1
      do k = 1, mzp
         vctr1(k) = (zt(k + 1) - zm(k))*dzm(k)
         vctr2(k) = (zm(k) - zt(k))*dzm(k)
         !    vt3dk(k) = dzt(k) / dne(k)
         vt3dk(k) = dzt(k)/(rho(k)*1.e-3)
         !print*,'VT3dk',k,dzt(k) , dne(k)
      end do

      !      scalarp => scalar_tab(n,ngrid)%var_p
      !      scalart => scalar_tab(n,ngrid)%var_t

      !- temp advection tendency (TT)
      scr1 = T
!!$    CALL fa_zc_plumerise(mzp,T,scr1(1),vt3dc(1),vt3df(1),vt3dg(1), &
!!$                         vt3dk(1),vctr1,vctr2)
      call fa_zc_plumerise(mzp, T, scr1, vt3dc, vt3df, vt3dg, &
                           vt3dk, vctr1, vctr2)

!!$    CALL advtndc_plumerise(mzp,T,scr1(1),TT,dt)
      call advtndc_plumerise(mzp, T, scr1, TT, dt)

      !- water vapor advection tendency (QVT)
      scr1 = QV
!!$    CALL fa_zc_plumerise(mzp,QV,scr1(1),vt3dc(1),vt3df(1),vt3dg(1), &
!!$                         vt3dk(1),vctr1,vctr2)
      call fa_zc_plumerise(mzp, QV, scr1, vt3dc, vt3df, vt3dg, &
                           vt3dk, vctr1, vctr2)

!!$    CALL advtndc_plumerise(mzp,QV,scr1(1),QVT,dt)
      call advtndc_plumerise(mzp, QV, scr1, QVT, dt)

      !- liquid advection tendency (QCT)
      scr1 = QC
!!$    CALL fa_zc_plumerise(mzp,QC,scr1(1),vt3dc(1),vt3df(1),vt3dg(1), &
!!$                         vt3dk(1),vctr1,vctr2)
      call fa_zc_plumerise(mzp, QC, scr1, vt3dc, vt3df, vt3dg, &
                           vt3dk, vctr1, vctr2)

!!$    CALL advtndc_plumerise(mzp,QC,scr1(1),QCT,dt)
      call advtndc_plumerise(mzp, QC, scr1, QCT, dt)

      !- ice advection tendency (QIT)
      scr1 = QI
!!$    CALL fa_zc_plumerise(mzp,QI,scr1(1),vt3dc(1),vt3df(1),vt3dg(1), &
!!$                         vt3dk(1),vctr1,vctr2)
      call fa_zc_plumerise(mzp, QI, scr1, vt3dc, vt3df, vt3dg, &
                           vt3dk, vctr1, vctr2)

!!$    CALL advtndc_plumerise(mzp,QI,scr1(1),QIT,dt)
      call advtndc_plumerise(mzp, QI, scr1, QIT, dt)

      !- hail/rain advection tendency (QHT)
      !   if(ak1 > 0. .or. ak2 > 0.) then

      scr1 = QH
!!$    CALL fa_zc_plumerise(mzp,QH,scr1(1),vt3dc(1),vt3df(1),vt3dg(1), &
!!$                         vt3dk(1),vctr1,vctr2)
      call fa_zc_plumerise(mzp, QH, scr1, vt3dc, vt3df, vt3dg, &
                           vt3dk, vctr1, vctr2)

!!$    CALL advtndc_plumerise(mzp,QH,scr1(1),QHT,dt)
      call advtndc_plumerise(mzp, QH, scr1, QHT, dt)
      !   endif

      !- horizontal wind advection tendency (VEL_T)
      scr1 = VEL_P
!!$    CALL fa_zc_plumerise(mzp                       &
!!$                            ,VEL_P     ,scr1  (1)  &
!!$                            ,vt3dc (1) ,vt3df (1)  &
!!$                            ,vt3dg (1) ,vt3dk (1)  &
!!$                            ,vctr1,vctr2             )
      call fa_zc_plumerise(mzp &
                           , VEL_P, scr1 &
                           , vt3dc, vt3df &
                           , vt3dg, vt3dk &
                           , vctr1, vctr2)

!!$    CALL advtndc_plumerise(mzp,VEL_P,scr1(1),VEL_T,dt)
      call advtndc_plumerise(mzp, VEL_P, scr1, VEL_T, dt)

      !- vertical radius transport

      scr1 = rad_p
!!$    CALL fa_zc_plumerise(mzp                   &
!!$                             ,rad_p     ,scr1  (1)  &
!!$                             ,vt3dc (1) ,vt3df (1)  &
!!$                             ,vt3dg (1) ,vt3dk (1)  &
!!$                             ,vctr1,vctr2               )
      call fa_zc_plumerise(mzp &
                           , rad_p, scr1 &
                           , vt3dc, vt3df &
                           , vt3dg, vt3dk &
                           , vctr1, vctr2)

!!$    CALL advtndc_plumerise(mzp,rad_p,scr1(1),rad_t,dt)
      call advtndc_plumerise(mzp, rad_p, scr1, rad_t, dt)

      return

   end subroutine scl_advectc_plumerise

   !----------------------------------------------------------------------------
   !
   subroutine fa_zc_plumerise(m1, scp, scr1, vt3dc, vt3df, vt3dg, vt3dk, vctr1, vctr2)

      integer, intent(IN)    :: m1
!!$    REAL,    INTENT(IN)    :: scp  (m1)
!!$    REAL,    INTENT(INOUT) :: scr1 (m1)
!!$    REAL,    INTENT(IN)    :: vt3dc(m1)
!!$    REAL,    INTENT(IN)    :: vt3df(m1)
!!$    REAL,    INTENT(INOUT) :: vt3dg(m1) ! (DMK) alterado (OUT) para (INOUT)
!!$    REAL,    INTENT(IN)    :: vt3dk(m1)
!!$    REAL,    INTENT(IN)    :: vctr1(m1)
!!$    REAL,    INTENT(IN)    :: vctr2(m1)
      real, intent(IN)    :: scp(:)
      real, intent(INOUT) :: scr1(:)
      real, intent(IN)    :: vt3dc(:)
      real, intent(IN)    :: vt3df(:)
      real, intent(INOUT) :: vt3dg(:) ! (DMK) alterado (OUT) para (INOUT)
      real, intent(IN)    :: vt3dk(:)
      real, intent(IN)    :: vctr1(:)
      real, intent(IN)    :: vctr2(:)

      integer :: k
      real    :: dfact

      dfact = .5

      ! Compute scalar flux VT3DG
      do k = 1, m1 - 1
         vt3dg(k) = vt3dc(k) &
                    *(vctr1(k)*scr1(k) &
                      + vctr2(k)*scr1(k + 1) &
                      + vt3df(k)*(scr1(k) - scr1(k + 1)))
      end do

      ! Modify fluxes to retain positive-definiteness on scalar quantities.
      !    If a flux will remove 1/2 quantity during a timestep,
      !    reduce to first order flux. This will remain positive-definite
      !    under the assumption that ABS(CFL(i)) + ABS(CFL(i-1)) < 1.0 if
      !    both fluxes are evacuating the box.

      do k = 1, m1 - 1
         if (vt3dc(k) .gt. 0.) then
            if (vt3dg(k)*vt3dk(k) .gt. dfact*scr1(k)) then
               vt3dg(k) = vt3dc(k)*scr1(k)
            end if
         elseif (vt3dc(k) .lt. 0.) then
            if (-vt3dg(k)*vt3dk(k + 1) .gt. dfact*scr1(k + 1)) then
               vt3dg(k) = vt3dc(k)*scr1(k + 1)
            end if
         end if

      end do

      ! Compute flux divergence

      do k = 2, m1 - 1
         scr1(k) = scr1(k) &
                   + vt3dk(k)*(vt3dg(k - 1) - vt3dg(k) &
                               + scp(k)*(vt3dc(k) - vt3dc(k - 1)))
      end do
      return
   end subroutine fa_zc_plumerise

   !----------------------------------------------------------------------------
   !
   subroutine advtndc_plumerise(m1, scp, sca, sct, dtl)

      integer, intent(IN)    :: m1
      real, intent(IN)    :: dtl
!!$    REAL,    INTENT(IN)    :: scp(m1)
!!$    REAL,    INTENT(IN)    :: sca(m1)
!!$    REAL,    INTENT(INOUT) :: sct(m1)
      real, intent(IN)    :: scp(:)
      real, intent(IN)    :: sca(:)
      real, intent(INOUT) :: sct(:)

      integer :: k
      real    :: dtli

      dtli = 1./dtl
      do k = 2, m1 - 1
         sct(k) = sct(k) + (sca(k) - scp(k))*dtli
      end do
      return
   end subroutine advtndc_plumerise

   !----------------------------------------------------------------------------
   !
   subroutine tend0_plumerise(nm1, wt, tt, qvt, qct, qht, qit, vel_t, rad_t)

      ! plumegen_coms
      integer, intent(IN)    :: nm1
!!$    REAL,    INTENT(INOUT) :: wt   (nkp)
!!$    REAL,    INTENT(INOUT) :: tt   (nkp)
!!$    REAL,    INTENT(INOUT) :: qvt  (nkp)
!!$    REAL,    INTENT(INOUT) :: qct  (nkp)
!!$    REAL,    INTENT(INOUT) :: qht  (nkp)
!!$    REAL,    INTENT(INOUT) :: qit  (nkp)
!!$    REAL,    INTENT(INOUT) :: vel_t(nkp)
!!$    REAL,    INTENT(INOUT) :: rad_t(nkp)
      real, intent(INOUT) :: wt(:)
      real, intent(INOUT) :: tt(:)
      real, intent(INOUT) :: qvt(:)
      real, intent(INOUT) :: qct(:)
      real, intent(INOUT) :: qht(:)
      real, intent(INOUT) :: qit(:)
      real, intent(INOUT) :: vel_t(:)
      real, intent(INOUT) :: rad_t(:)

      wt(1:nm1) = 0.
      tt(1:nm1) = 0.
      qvt(1:nm1) = 0.
      qct(1:nm1) = 0.
      qht(1:nm1) = 0.
      qit(1:nm1) = 0.
      vel_t(1:nm1) = 0.
      rad_t(1:nm1) = 0.
      !sct (1:nm1)  = 0.
   end subroutine tend0_plumerise

   !     ****************************************************************
   subroutine scl_misc(m1, wbar, w, adiabat, alpha, radius, tt, t, te, qvt, qv, &
                       qvenv, qct, qc, qht, qh, qit, qi, vel_e, vel_p, vel_t, rad_p, &
                       rad_t)

      integer, intent(IN)    :: m1
      ! plumegen_coms
      real, intent(INOUT) :: wbar    ! (DMK) alterado (OUT) para (INOUT)
      real, intent(INOUT) :: adiabat ! (DMK) alterado (OUT) para (INOUT)
      real, intent(IN)    :: alpha
!!$    REAL,    INTENT(IN)    :: w     (nkp)
!!$    REAL,    INTENT(IN)    :: radius(nkp)
!!$    REAL,    INTENT(INOUT) :: tt    (nkp)
!!$    REAL,    INTENT(IN)    :: t     (nkp)
!!$    REAL,    INTENT(IN)    :: te    (nkp)
!!$    REAL,    INTENT(INOUT) :: qvt   (nkp)
!!$    REAL,    INTENT(IN)    :: qv    (nkp)
!!$    REAL,    INTENT(IN)    :: qvenv (nkp)
!!$    REAL,    INTENT(INOUT) :: qct   (nkp)
!!$    REAL,    INTENT(IN)    :: qc    (nkp)
!!$    REAL,    INTENT(INOUT) :: qht   (nkp)
!!$    REAL,    INTENT(IN)    :: qh    (nkp)
!!$    REAL,    INTENT(INOUT) :: qit   (nkp)
!!$    REAL,    INTENT(IN)    :: qi    (nkp)
!!$    REAL,    INTENT(IN)    :: vel_e (nkp)
!!$    REAL,    INTENT(IN)    :: vel_p (nkp)
!!$    REAL,    INTENT(INOUT) :: vel_t (nkp)
!!$    REAL,    INTENT(INOUT) :: rad_T (nkp)
!!$    REAL,    INTENT(IN)    :: rad_p (nkp)
      real, intent(IN)    :: w(:)
      real, intent(IN)    :: radius(:)
      real, intent(INOUT) :: tt(:)
      real, intent(IN)    :: t(:)
      real, intent(IN)    :: te(:)
      real, intent(INOUT) :: qvt(:)
      real, intent(IN)    :: qv(:)
      real, intent(IN)    :: qvenv(:)
      real, intent(INOUT) :: qct(:)
      real, intent(IN)    :: qc(:)
      real, intent(INOUT) :: qht(:)
      real, intent(IN)    :: qh(:)
      real, intent(INOUT) :: qit(:)
      real, intent(IN)    :: qi(:)
      real, intent(IN)    :: vel_e(:)
      real, intent(IN)    :: vel_p(:)
      real, intent(INOUT) :: vel_t(:)
      real, intent(INOUT) :: rad_T(:)
      real, intent(IN)    :: rad_p(:)

      real, parameter :: g = 9.81, cp = 1004.
      integer         :: k
      real            :: dmdtm

      do k = 2, m1 - 1
         WBAR = 0.5*(W(k) + W(k - 1))
         !-- dry adiabat
         ADIABAT = -WBAR*G/CP
         !
         !-- entrainment
         DMDTM = 2.*ALPHA*abs(WBAR)/RADIUS(k)  != (1/M)DM/DT

         !-- tendency temperature = adv + adiab + entrainment
         TT(k) = TT(K) + ADIABAT - DMDTM*(T(k) - TE(k))

         !-- tendency water vapor = adv  + entrainment
         QVT(K) = QVT(K) - DMDTM*(QV(k) - QVENV(k))

         QCT(K) = QCT(K) - DMDTM*(QC(k))
         QHT(K) = QHT(K) - DMDTM*(QH(k))
         QIT(K) = QIT(K) - DMDTM*(QI(k))

         !-- tendency horizontal speed = adv  + entrainment
         VEL_T(K) = VEL_T(K) - DMDTM*(VEL_P(k) - VEL_E(k))

         !-- tendency horizontal speed = adv  + entrainment
         rad_t(K) = rad_t(K) + 0.5*DMDTM*(6./5.)*RADIUS(k)

         !-- tendency gas/particle = adv  + entrainment
         !      SCT(K) = SCT(K)         - DMDTM * ( SC (k) -   SCE (k) )

      end do
   end subroutine scl_misc

!     ****************************************************************
   subroutine scl_dyn_entrain(m1, wbar, w, adiabat, alpha, radius, tt, t, &
                              te, qvt, qv, qvenv, qct, qc, qht, qh, qit, qi, &
                              vel_e, vel_p, vel_t, rad_p, rad_t)

      integer, intent(IN)    :: m1
      ! plumegen_coms
      real, intent(INOUT) :: wbar
      real, intent(INOUT) :: adiabat
      real, intent(IN)    :: alpha
!!$    REAL,    INTENT(IN)    :: w     (nkp)
!!$    REAL,    INTENT(IN)    :: radius(nkp)
!!$    REAL,    INTENT(INOUT) :: tt    (nkp)
!!$    REAL,    INTENT(IN)    :: t     (nkp)
!!$    REAL,    INTENT(IN)    :: te    (nkp)
!!$    REAL,    INTENT(INOUT) :: qvt   (nkp)
!!$    REAL,    INTENT(IN)    :: qv    (nkp)
!!$    REAL,    INTENT(IN)    :: qvenv (nkp)
!!$    REAL,    INTENT(INOUT) :: qct   (nkp)
!!$    REAL,    INTENT(IN)    :: qc    (nkp)
!!$    REAL,    INTENT(INOUT) :: qht   (nkp)
!!$    REAL,    INTENT(IN)    :: qh    (nkp)
!!$    REAL,    INTENT(INOUT) :: qit   (nkp)
!!$    REAL,    INTENT(IN)    :: qi    (nkp)
!!$    REAL,    INTENT(IN)    :: vel_e (nkp)
!!$    REAL,    INTENT(IN)    :: vel_p (nkp)
!!$    REAL,    INTENT(INOUT) :: vel_t (nkp)
!!$    REAL,    INTENT(INOUT) :: rad_T (nkp)
!!$    REAL,    INTENT(IN)    :: rad_p (nkp)
      real, intent(IN)    :: w(:)
      real, intent(IN)    :: radius(:)
      real, intent(INOUT) :: tt(:)
      real, intent(IN)    :: t(:)
      real, intent(IN)    :: te(:)
      real, intent(INOUT) :: qvt(:)
      real, intent(IN)    :: qv(:)
      real, intent(IN)    :: qvenv(:)
      real, intent(INOUT) :: qct(:)
      real, intent(IN)    :: qc(:)
      real, intent(INOUT) :: qht(:)
      real, intent(IN)    :: qh(:)
      real, intent(INOUT) :: qit(:)
      real, intent(IN)    :: qi(:)
      real, intent(IN)    :: vel_e(:)
      real, intent(IN)    :: vel_p(:)
      real, intent(INOUT) :: vel_t(:)
      real, intent(INOUT) :: rad_T(:)
      real, intent(IN)    :: rad_p(:)

      real, parameter :: g = 9.81, cp = 1004., pi = 3.1416
      integer         :: k
      real            :: dmdtm

      do k = 2, m1 - 1
         !
         !-- tendency horizontal radius from dyn entrainment
         !rad_t(K) = rad_t(K)   +        (vel_e(k)-vel_p(k)) /pi
         rad_t(K) = rad_t(K) + abs((vel_e(k) - vel_p(k)))/pi

         !-- entrainment
         !DMDTM = (2./3.1416)  *     (VEL_E (k) - VEL_P (k)) / RADIUS (k)
         DMDTM = (2./3.1416)*abs(VEL_E(k) - VEL_P(k))/RADIUS(k)

         !-- tendency horizontal speed  from dyn entrainment
         VEL_T(K) = VEL_T(K) - DMDTM*(VEL_P(k) - VEL_E(k))

         !     if(VEL_P (k) - VEL_E (k) > 0.) cycle

         !-- tendency temperature  from dyn entrainment
         TT(k) = TT(K) - DMDTM*(T(k) - TE(k))

         !-- tendency water vapor  from dyn entrainment
         QVT(K) = QVT(K) - DMDTM*(QV(k) - QVENV(k))

         QCT(K) = QCT(K) - DMDTM*(QC(k))
         QHT(K) = QHT(K) - DMDTM*(QH(k))
         QIT(K) = QIT(K) - DMDTM*(QI(k))

         !-- tendency gas/particle  from dyn entrainment
         !         SCT(K) = SCT(K)         - DMDTM * ( SC (k) - SCE (k) )

      end do
   end subroutine scl_dyn_entrain
   !     ****************************************************************

   subroutine visc_W(m1, kmt, zt, visc, zm, w, t, qv, qh, qc, qi, wt, tt, qvt, &
                     qct, qht, qit, vel_p, vel_t, rad_p, rad_t)

      integer, intent(IN)    :: m1
      integer, intent(IN)    :: kmt
      ! plumegen_coms
!!$    REAL,    INTENT(IN)    :: zt   (nkp)
!!$    REAL,    INTENT(IN)    :: visc (nkp)
!!$    REAL,    INTENT(IN)    :: zm   (nkp)
!!$    REAL,    INTENT(IN)    :: w    (nkp)
!!$    REAL,    INTENT(IN)    :: t    (nkp)
!!$    REAL,    INTENT(IN)    :: qv   (nkp)
!!$    REAL,    INTENT(IN)    :: qh   (nkp)
!!$    REAL,    INTENT(IN)    :: qc   (nkp)
!!$    REAL,    INTENT(IN)    :: qi   (nkp)
!!$    REAL,    INTENT(INOUT) :: wt   (nkp)
!!$    REAL,    INTENT(INOUT) :: tt   (nkp)
!!$    REAL,    INTENT(INOUT) :: qvt  (nkp)
!!$    REAL,    INTENT(INOUT) :: qct  (nkp)
!!$    REAL,    INTENT(INOUT) :: qht  (nkp)
!!$    REAL,    INTENT(INOUT) :: qit  (nkp)
!!$    REAL,    INTENT(IN)    :: vel_p(nkp)
!!$    REAL,    INTENT(INOUT) :: vel_t(nkp)
!!$    REAL,    INTENT(INOUT) :: rad_T(nkp)
!!$    REAL,    INTENT(IN)    :: rad_p(nkp)
      real, intent(IN)    :: zt(:)
      real, intent(IN)    :: visc(:)
      real, intent(IN)    :: zm(:)
      real, intent(IN)    :: w(:)
      real, intent(IN)    :: t(:)
      real, intent(IN)    :: qv(:)
      real, intent(IN)    :: qh(:)
      real, intent(IN)    :: qc(:)
      real, intent(IN)    :: qi(:)
      real, intent(INOUT) :: wt(:)
      real, intent(INOUT) :: tt(:)
      real, intent(INOUT) :: qvt(:)
      real, intent(INOUT) :: qct(:)
      real, intent(INOUT) :: qht(:)
      real, intent(INOUT) :: qit(:)
      real, intent(IN)    :: vel_p(:)
      real, intent(INOUT) :: vel_t(:)
      real, intent(INOUT) :: rad_T(:)
      real, intent(IN)    :: rad_p(:)

      integer :: k, m2
      real    :: dz1t, dz1m, dz2t, dz2m, d2wdz, d2tdz, &
                 d2qvdz, d2qhdz, d2qcdz, d2qidz, &
                 d2vel_pdz, d2rad_dz

      !srf--- 17/08/2005
      !m2=min(m1+deltak,kmt)
      !m2=MIN(m1,kmt)

      !do k=2,m1-1
      do k = 2, min(m1, kmt - 1)!m2-1
         DZ1T = 0.5*(ZT(K + 1) - ZT(K - 1))
         DZ2T = VISC(k)/(DZ1T*DZ1T)
         DZ1M = 0.5*(ZM(K + 1) - ZM(K - 1))
         DZ2M = VISC(k)/(DZ1M*DZ1M)
         D2WDZ = (W(k + 1) - 2*W(k) + W(k - 1))*DZ2M
         D2TDZ = (T(k + 1) - 2*T(k) + T(k - 1))*DZ2T
         D2QVDZ = (QV(k + 1) - 2*QV(k) + QV(k - 1))*DZ2T
         D2QHDZ = (QH(k + 1) - 2*QH(k) + QH(k - 1))*DZ2T
         D2QCDZ = (QC(k + 1) - 2*QC(k) + QC(k - 1))*DZ2T
         D2QIDZ = (QI(k + 1) - 2*QI(k) + QI(k - 1))*DZ2T
         !D2SCDZ = (SC (k + 1) - 2 * SC (k) + SC (k - 1) ) * DZ2T

         d2vel_pdz = (vel_P(k + 1) - 2*vel_P(k) + vel_P(k - 1))*DZ2T
         d2rad_dz = (rad_p(k + 1) - 2*rad_p(k) + rad_p(k - 1))*DZ2T

         WT(k) = WT(k) + D2WDZ
         TT(k) = TT(k) + D2TDZ
         QVT(k) = QVT(k) + D2QVDZ
         QCT(k) = QCT(k) + D2QCDZ
         QHT(k) = QHT(k) + D2QHDZ
         QIT(k) = QIT(k) + D2QIDZ

         vel_t(k) = vel_t(k) + d2vel_pdz
         rad_t(k) = rad_t(k) + d2rad_dz

         !SCT(k) =  SCT(k) + D2SCDZ
         !print*,'W-VISC=',k,D2WDZ
      end do

   end subroutine visc_W

   !     ****************************************************************

   subroutine update_plumerise(m1, varn, wt, dt, tt, qvt, qct, qht, w, t, qv, qc, &
                               qh, qit, qi, vel_p, vel_t, rad_p, rad_t)

      integer, intent(IN)    :: m1
      character(len=*), intent(IN)    :: varn
      ! plumegen_coms
      real, intent(IN)    :: dt
!!$    REAL,             INTENT(IN)    :: wt   (nkp)
!!$    REAL,             INTENT(IN)    :: tt   (nkp)
!!$    REAL,             INTENT(IN)    :: qvt  (nkp)
!!$    REAL,             INTENT(IN)    :: qct  (nkp)
!!$    REAL,             INTENT(IN)    :: qht  (nkp)
!!$    REAL,             INTENT(INOUT) :: w    (nkp)
!!$    REAL,             INTENT(INOUT) :: t    (nkp)
!!$    REAL,             INTENT(INOUT) :: qv   (nkp)
!!$    REAL,             INTENT(INOUT) :: qc   (nkp)
!!$    REAL,             INTENT(INOUT) :: qh   (nkp)
!!$    REAL,             INTENT(IN)    :: qit  (nkp) ! (DMK) alterado (INOUT) (IN)
!!$    REAL,             INTENT(INOUT) :: qi   (nkp)
!!$    REAL,             INTENT(INOUT) :: vel_p(nkp)
!!$    REAL,             INTENT(IN)    :: vel_t(nkp)
!!$    REAL,             INTENT(IN)    :: rad_T(nkp)
!!$    REAL,             INTENT(INOUT) :: rad_p(nkp)
      real, intent(IN)    :: wt(:)
      real, intent(IN)    :: tt(:)
      real, intent(IN)    :: qvt(:)
      real, intent(IN)    :: qct(:)
      real, intent(IN)    :: qht(:)
      real, intent(INOUT) :: w(:)
      real, intent(INOUT) :: t(:)
      real, intent(INOUT) :: qv(:)
      real, intent(INOUT) :: qc(:)
      real, intent(INOUT) :: qh(:)
      real, intent(IN)    :: qit(:) ! (DMK) alterado (INOUT) (IN)
      real, intent(INOUT) :: qi(:)
      real, intent(INOUT) :: vel_p(:)
      real, intent(IN)    :: vel_t(:)
      real, intent(IN)    :: rad_T(:)
      real, intent(INOUT) :: rad_p(:)

      integer :: k

      if (varn == 'W') then

         do k = 2, m1 - 1
            W(k) = W(k) + WT(k)*DT
         end do
         return

      else
         do k = 2, m1 - 1
            T(k) = T(k) + TT(k)*DT

            QV(k) = QV(k) + QVT(k)*DT

            QC(k) = QC(k) + QCT(k)*DT !cloud drops travel with air
            QH(k) = QH(k) + QHT(k)*DT
            QI(k) = QI(k) + QIT(k)*DT
            ! SC(k) = SC(k) + SCT(k) * DT

            !srf---18jun2005
            QV(k) = max(0., QV(k))
            QC(k) = max(0., QC(k))
            QH(k) = max(0., QH(k))
            QI(k) = max(0., QI(k))

            VEL_P(k) = VEL_P(k) + VEL_T(k)*DT

            rad_p(k) = rad_p(k) + rad_t(k)*DT

            ! SC(k) = max(0., SC(k))

         end do
      end if
   end subroutine update_plumerise

   !----------------------------------------------------------------------------
   !
   subroutine fallpart(m1, rho, vth, vhrel, w, cvh, vti, cvi, qh, qi, zm, qht, qit)

      integer, intent(IN)    :: m1
      ! plumegen_coms
      real, intent(INOUT) :: vhrel   ! (DMK) alterado (OUT) para (INOUT)
!!$    REAL,    INTENT(IN)    :: rho(nkp)
!!$    REAL,    INTENT(INOUT) :: vth(nkp)
!!$    REAL,    INTENT(IN)    :: w  (nkp)
!!$    REAL,    INTENT(INOUT) :: cvh(nkp) ! (DMK) alterado (OUT) para (INOUT)
!!$    REAL,    INTENT(INOUT) :: vti(nkp)
!!$    REAL,    INTENT(INOUT) :: cvi(nkp) ! (DMK) alterado (OUT) para (INOUT)
!!$    REAL,    INTENT(IN)    :: qh (nkp)
!!$    REAL,    INTENT(IN)    :: qi (nkp)
!!$    REAL,    INTENT(IN)    :: zm (nkp)
!!$    REAL,    INTENT(INOUT) :: qht(nkp)
!!$    REAL,    INTENT(INOUT) :: qit(nkp)
      real, intent(IN)    :: rho(:)
      real, intent(INOUT) :: vth(:)
      real, intent(IN)    :: w(:)
      real, intent(INOUT) :: cvh(:) ! (DMK) alterado (OUT) para (INOUT)
      real, intent(INOUT) :: vti(:)
      real, intent(INOUT) :: cvi(:) ! (DMK) alterado (OUT) para (INOUT)
      real, intent(IN)    :: qh(:)
      real, intent(IN)    :: qi(:)
      real, intent(IN)    :: zm(:)
      real, intent(INOUT) :: qht(:)
      real, intent(INOUT) :: qit(:)

      !srf==================================
      !   verificar se o gradiente esta correto
      !
      !srf==================================
      !
      !     XNO=1.E7  [m**-4] median volume diameter raindrop,Kessler
      !     VC = 38.3/(XNO**.125), median volume fallspeed eqn., Kessler
      !     for ice, see (OT18), use F0=0.75 per argument there. rho*q
      !     values are in g/m**3, velocities in m/s

      real, parameter :: VCONST = 5.107387, EPS = 0.622, F0 = 0.75
      real, parameter :: G = 9.81, CP = 1004.

      integer         :: k
      real            :: vtc, dfhz, dfiz, dz1
      real            :: virel ! (DMK) scratch

      !
      do k = 2, m1 - 1

         VTC = VCONST*RHO(k)**.125   ! median volume fallspeed (KTable4)

         !  hydrometeor assembly velocity calculations (K Table4)
         !  VTH(k)=-VTC*QH(k)**.125  !median volume fallspeed, water
         VTH(k) = -4.            !small variation with qh

         VHREL = W(k) + VTH(k)  !relative to surrounding cloud

         !  rain ventilation coefficient for evaporation
         CVH(k) = 1.6 + 0.57e-3*(abs(VHREL))**1.5
         !
         !  VTI(k)=-VTC*F0*QI(k)**.125    !median volume fallspeed,ice
         VTI(k) = -3.                !small variation with qi

         VIREL = W(k) + VTI(k)       !relative to surrounding cloud
         !
         !  ice ventilation coefficient for sublimation
         CVI(k) = 1.6 + 0.57e-3*(abs(VIREL))**1.5/F0
         !
         !
         if (VHREL .ge. 0.0) then
            DFHZ = QH(k)*(RHO(k)*VTH(k) - RHO(k - 1)*VTH(k - 1))/RHO(k - 1)
         else
            DFHZ = QH(k)*(RHO(k + 1)*VTH(k + 1) - RHO(k)*VTH(k))/RHO(k)
         end if
         !
         !
         if (VIREL .ge. 0.0) then
            DFIZ = QI(k)*(RHO(k)*VTI(k) - RHO(k - 1)*VTI(k - 1))/RHO(k - 1)
         else
            DFIZ = QI(k)*(RHO(k + 1)*VTI(k + 1) - RHO(k)*VTI(k))/RHO(k)
         end if

         DZ1 = ZM(K) - ZM(K - 1)

         qht(k) = qht(k) - DFHZ/DZ1 !hydrometeors don't

         qit(k) = qit(k) - DFIZ/DZ1  !nor does ice? hail, what about

      end do

   end subroutine fallpart

   !----------------------------------------------------------------------------
   !
   subroutine printout(izprint, nrectotal, mintime, dt, time, ztop, pe, t, te, qv, &
                       qsat, qc, qh, qi, zt, w, vth, sc, qvenv, nrec)

      integer, intent(IN)    :: izprint
      integer, intent(IN)    :: nrectotal
      ! plumegen_coms
      integer, intent(IN)    :: mintime
      real, intent(IN)    :: dt
      real, intent(IN)    :: time
      real, intent(IN)    :: ztop
!!$    REAL,    INTENT(IN)    :: pe   (nkp)
!!$    REAL,    INTENT(IN)    :: t    (nkp)
!!$    REAL,    INTENT(IN)    :: te   (nkp)
!!$    REAL,    INTENT(IN)    :: qv   (nkp)
!!$    REAL,    INTENT(IN)    :: qsat (nkp)
!!$    REAL,    INTENT(IN)    :: qc   (nkp)
!!$    REAL,    INTENT(IN)    :: qh   (nkp)
!!$    REAL,    INTENT(IN)    :: qi   (nkp)
!!$    REAL,    INTENT(IN)    :: zt   (nkp)
!!$    REAL,    INTENT(IN)    :: w    (nkp)
!!$    REAL,    INTENT(IN)    :: vth  (nkp)
!!$    REAL,    INTENT(IN)    :: sc   (nkp)
!!$    REAL,    INTENT(IN)    :: qvenv(nkp)
      real, intent(IN)    :: pe(:)
      real, intent(IN)    :: t(:)
      real, intent(IN)    :: te(:)
      real, intent(IN)    :: qv(:)
      real, intent(IN)    :: qsat(:)
      real, intent(IN)    :: qc(:)
      real, intent(IN)    :: qh(:)
      real, intent(IN)    :: qi(:)
      real, intent(IN)    :: zt(:)
      real, intent(IN)    :: w(:)
      real, intent(IN)    :: vth(:)
      real, intent(IN)    :: sc(:)
      real, intent(IN)    :: qvenv(:)
      ! save attribute
      integer, intent(INOUT) :: nrec

      real, parameter :: tmelt = 273.3
      integer         :: ko, interval
      real            :: pea, btmp, etmp, vap1, vap2, gpkc, gpkh, gpki, deficit

      interval = 1              !debug time interval,min

      !
      if (IZPRINT .eq. 0) return

      if (MINTIME == 1) nrec = 0
      !
      write (2, 430) MINTIME, DT, TIME
      write (2, 431) ZTOP
      write (2, 380)
      !
      ! do the print
      !
      do 390 KO = 1, nrectotal, interval

         PEA = PE(KO)*10.       !pressure is stored in decibars(kPa),print in mb;
         BTMP = T(KO) - TMELT     !temps in Celsius
         ETMP = T(KO) - TE(KO)   !temperature excess
         VAP1 = QV(KO)*1000.  !printout in g/kg for all water,
         VAP2 = QSAT(KO)*1000.  !vapor (internal storage is in g/g)
         GPKC = QC(KO)*1000.  !cloud water
         GPKH = QH(KO)*1000.  !raindrops
         GPKI = QI(KO)*1000.  !ice particles
         DEFICIT = VAP2 - VAP1     !vapor deficit
         !
         write (2, 400) zt(KO)/1000., PEA, W(KO), BTMP, ETMP, VAP1, &
            VAP2, GPKC, GPKH, GPKI, VTH(KO), SC(KO)
         !
         !
         !                                    !end of printout

390      continue

         nrec = nrec + 1
         write (19, rec=nrec) (W(KO), KO=1, nrectotal)
         nrec = nrec + 1
         write (19, rec=nrec) (T(KO), KO=1, nrectotal)
         nrec = nrec + 1
         write (19, rec=nrec) (TE(KO), KO=1, nrectotal)
         nrec = nrec + 1
         write (19, rec=nrec) (QV(KO)*1000., KO=1, nrectotal)
         nrec = nrec + 1
         write (19, rec=nrec) (QC(KO)*1000., KO=1, nrectotal)
         nrec = nrec + 1
         write (19, rec=nrec) (QH(KO)*1000., KO=1, nrectotal)
         nrec = nrec + 1
         write (19, rec=nrec) (QI(KO)*1000., KO=1, nrectotal)
         nrec = nrec + 1
         !   write (19,rec=nrec) (SC(KO), KO=1,nrectotal)
         write (19, rec=nrec) (QSAT(KO)*1000., KO=1, nrectotal)
         nrec = nrec + 1
         write (19, rec=nrec) (QVENV(KO)*1000., KO=1, nrectotal)

         print *, 'ntimes=', nrec/(9)
         !
         return
         !
         ! ************** FORMATS *********************************************
         !
380      format(/, ' Z(KM) P(MB) W(MPS) T(C)  T-TE   VAP   SAT   QC    QH' &
                 '     QI    VTH(MPS) SCAL'/)
         !
400      format(1h, F4.1, F8.2, F8.2, F7.1, 6f6.2, F7.2, 1x, F6.2)
         !
430      format(1h, //I5, ' MINUTES       DT= ', F6.2, ' SECONDS   TIME= ' &
                , F8.2, ' SECONDS')
431      format(' ZTOP= ', F10.2)
         !
         end subroutine printout
         !
         ! *********************************************************************

         subroutine WATERBAL(qc, l, qh, qi, qv, t, qsat, wbar, dqsdz, dt, rho, est, cvi)

            ! plumegen_coms
            integer, intent(IN)    :: l
            real, intent(IN)    :: wbar
            real, intent(IN)    :: dqsdz
            real, intent(IN)    :: dt
!!$       REAL,    INTENT(INOUT) :: qc  (nkp)
!!$       REAL,    INTENT(INOUT) :: qh  (nkp)
!!$       REAL,    INTENT(INOUT) :: qi  (nkp)
!!$       REAL,    INTENT(INOUT) :: qv  (nkp)
!!$       REAL,    INTENT(INOUT) :: t   (nkp)
!!$       REAL,    INTENT(IN)    :: qsat(nkp)
!!$       REAL,    INTENT(IN)    :: rho (nkp)
!!$       REAL,    INTENT(IN)    :: est (nkp)
!!$       REAL,    INTENT(IN)    :: cvi (nkp)
            real, intent(INOUT) :: qc(:)
            real, intent(INOUT) :: qh(:)
            real, intent(INOUT) :: qi(:)
            real, intent(INOUT) :: qv(:)
            real, intent(INOUT) :: t(:)
            real, intent(IN)    :: qsat(:)
            real, intent(IN)    :: rho(:)
            real, intent(IN)    :: est(:)
            real, intent(IN)    :: cvi(:)
            !

            if (QC(L) .le. 1.0e-10) QC(L) = 0.  !DEFEAT UNDERFLOW PROBLEM
            if (QH(L) .le. 1.0e-10) QH(L) = 0.
            if (QI(L) .le. 1.0e-10) QI(L) = 0.
            !
            call EVAPORATE(qc, qh, qi, qv, t, qsat, l, wbar, dqsdz, dt, rho, est, cvi)  !vapor to cloud,cloud to vapor
            !
            call SUBLIMATE(t, l, qv, qsat, rho, qi, est, dt)  !vapor to ice
            !
            call GLACIATE(qh, l, qv, qsat, t, dt, qi) !rain to ice

            call MELT(qi, l, t, dt, rho, cvi, qh)  !ice to rain
            !
            !if(ak1 > 0. .or. ak2 > 0.) &
            call CONVERT(t, l, qc, rho, qh, dt) !(auto)conversion and accretion
            !CALL CONVERT2 () !(auto)conversion and accretion
            !

            return
         end subroutine WATERBAL

         ! *********************************************************************
         subroutine EVAPORATE(qc, qh, qi, qv, t, qsat, l, wbar, dqsdz, dt, rho, est, cvi)
            !
            !- evaporates cloud,rain and ice to saturation
            !

            ! plumegen_coms
            integer, intent(IN)    :: l
            real, intent(IN)    :: wbar
            real, intent(IN)    :: dqsdz
            real, intent(IN)    :: dt
!!$       REAL,    INTENT(INOUT) :: qc  (nkp)
!!$       REAL,    INTENT(INOUT) :: qh  (nkp)
!!$       REAL,    INTENT(INOUT) :: qi  (nkp)
!!$       REAL,    INTENT(INOUT) :: qv  (nkp)
!!$       REAL,    INTENT(INOUT) :: t   (nkp)
!!$       REAL,    INTENT(IN)    :: qsat(nkp)
!!$       REAL,    INTENT(IN)    :: rho (nkp)
!!$       REAL,    INTENT(IN)    :: est (nkp)
!!$       REAL,    INTENT(IN)    :: cvi (nkp)
            real, intent(INOUT) :: qc(:)
            real, intent(INOUT) :: qh(:)
            real, intent(INOUT) :: qi(:)
            real, intent(INOUT) :: qv(:)
            real, intent(INOUT) :: t(:)
            real, intent(IN)    :: qsat(:)
            real, intent(IN)    :: rho(:)
            real, intent(IN)    :: est(:)
            real, intent(IN)    :: cvi(:)

            !
            !     XNO=10.0E06
            !     HERC = 1.93*1.E-6*XN035        !evaporation constant
            !
            real, parameter :: HERC = 5.44e-4, CP = 1.004, HEATCOND = 2.5e3
            real, parameter :: HEATSUBL = 2834., TMELT = 273., TFREEZE = 269.3
            real, parameter :: FRC = HEATCOND/CP, SRC = HEATSUBL/CP

            real            :: evhdt, evidt, evrate, evap, sd, &
                               quant, dividend, divisor, devidt

            !
            !
            SD = QSAT(L) - QV(L)  !vapor deficit
            if (SD .eq. 0.0) return
            !IF (abs(SD).lt.1.e-7)  RETURN

            EVHDT = 0.
            EVIDT = 0.
            !evrate =0.; evap=0.; sd=0.0; quant=0.0; dividend=0.0; divisor=0.0; devidt=0.0

            EVRATE = abs(WBAR*DQSDZ)   !evaporation rate (Kessler 8.32)
            EVAP = EVRATE*DT            !what we can get in DT

            if (SD .le. 0.0) then  !     condense. SD is negative

               if (EVAP .ge. abs(SD)) then    !we get it all

                  QC(L) = QC(L) - SD  !deficit,remember?
                  QV(L) = QSAT(L)       !set the vapor to saturation
                  T(L) = T(L) - SD*FRC  !heat gained through condensation
                  !per gram of dry air
                  return

               else

                  QC(L) = QC(L) + EVAP         !get what we can in DT
                  QV(L) = QV(L) - EVAP         !remove it from the vapor
                  T(L) = T(L) + EVAP*FRC   !get some heat

                  return

               end if
               !
            else                                !SD is positive, need some water
               !
               ! not saturated. saturate if possible. use everything in order
               ! cloud, rain, ice. SD is positive

               if (EVAP .le. QC(L)) then        !enough cloud to last DT
                  !

                  if (SD .le. EVAP) then          !enough time to saturate

                     QC(L) = QC(L) - SD       !remove cloud
                     QV(L) = QSAT(L)          !saturate
                     T(L) = T(L) - SD*FRC   !cool the parcel
                     return  !done
                     !

                  else   !not enough time

                     SD = SD - EVAP               !use what there is
                     QV(L) = QV(L) + EVAP     !add vapor
                     T(L) = T(L) - EVAP*FRC !lose heat
                     QC(L) = QC(L) - EVAP     !lose cloud
                     !go on to rain.
                  end if
                  !
               else                !not enough cloud to last DT
                  !
                  if (SD .le. QC(L)) then   !but there is enough to sat

                     QV(L) = QSAT(L)  !use it
                     QC(L) = QC(L) - SD
                     T(L) = T(L) - SD*FRC
                     return

                  else            !not enough to sat
                     SD = SD - QC(L)
                     QV(L) = QV(L) + QC(L)
                     T(L) = T(L) - QC(L)*FRC
                     QC(L) = 0.0  !all gone

                  end if       !on to rain
               end if          !finished with cloud
               !
               !  but still not saturated, so try to use some rain
               !  this is tricky, because we only have time DT to evaporate.
               !  if there is enough rain, we can evaporate it for dt. ice can
               !  also sublimate at the same time. there is a compromise
               !  here.....use rain first, then
               !  ice. saturation may not be possible in one DT time.
               !  rain evaporation rate (W12),(OT25),(K Table 4).
               !  evaporate rain first
               !  sd is still positive or we wouldn't be here.

               if (QH(L) .le. 1.e-10) goto 33

               !srf-25082005
               !  QUANT = ( QC (L)  + QV (L) - QSAT (L) ) * RHO (L)   !g/m**3
               QUANT = (QSAT(L) - QC(L) - QV(L))*RHO(L)   !g/m**3
               !
               EVHDT = (DT*HERC*(QUANT)*(QH(L)*RHO(L))**.65)/RHO(L)
               !             rain evaporation in time DT

               if (EVHDT .le. QH(L)) then           !enough rain to last DT

                  if (SD .le. EVHDT) then                  !enough time to saturate
                     QH(L) = QH(L) - SD           !remove rain
                     QV(L) = QSAT(L)                  !saturate
                     T(L) = T(L) - SD*FRC          !cool the parcel

                     return                          !done
                     !
                  else                               !not enough time
                     SD = SD - EVHDT                   !use what there is
                     QV(L) = QV(L) + EVHDT           !add vapor
                     T(L) = T(L) - EVHDT*FRC           !lose heat
                     QH(L) = QH(L) - EVHDT           !lose rain

                  end if                                    !go on to ice.
                  !
               else  !not enough rain to last DT
                  !
                  if (SD .le. QH(L)) then             !but there is enough to sat
                     QV(L) = QSAT(L)                !use it
                     QH(L) = QH(L) - SD
                     T(L) = T(L) - SD*FRC
                     return
                     !
                  else                              !not enough to sat
                     SD = SD - QH(L)
                     QV(L) = QV(L) + QH(L)
                     T(L) = T(L) - QH(L)*FRC
                     QH(L) = 0.0                   !all gone

                  end if                             !on to ice
                  !

               end if                                !finished with rain
               !
               !
               !  now for ice
               !  equation from (OT); correction factors for units applied
               !
33             continue
               if (QI(L) .le. 1.e-10) return            !no ice there
               !
               DIVIDEND = ((1.e6/RHO(L))**0.475)*(SD/QSAT(L) &
                                                  - 1)*(QI(L)**0.525)*1.13
               DIVISOR = 7.e5 + 4.1e6/(10.*EST(L))

               DEVIDT = -CVI(L)*DIVIDEND/DIVISOR   !rate of change

               EVIDT = DEVIDT*DT                      !what we could get
               !
               ! logic here is identical to rain. could get fancy and make
               ! subroutine but duplication of code is easier. God bless the
               ! screen editor.
               !

               if (EVIDT .le. QI(L)) then             !enough ice to last DT
                  !

                  if (SD .le. EVIDT) then                    !enough time to saturate
                     QI(L) = QI(L) - SD             !remove ice
                     QV(L) = QSAT(L)                    !saturate
                     T(L) = T(L) - SD*SRC            !cool the parcel

                     return                            !done
                     !

                  else                                !not enough time

                     SD = SD - EVIDT                    !use what there is
                     QV(L) = QV(L) + EVIDT            !add vapor
                     T(L) = T(L) - EVIDT*SRC            !lose heat
                     QI(L) = QI(L) - EVIDT            !lose ice

                  end if                                    !go on,unsatisfied
                  !
               else                                   !not enough ice to last DT
                  !
                  if (SD .le. QI(L)) then             !but there is enough to sat

                     QV(L) = QSAT(L)                !use it
                     QI(L) = QI(L) - SD
                     T(L) = T(L) - SD*SRC

                     return
                     !
                  else                                 !not enough to sat
                     SD = SD - QI(L)
                     QV(L) = QV(L) + QI(L)
                     T(L) = T(L) - QI(L)*SRC
                     QI(L) = 0.0                      !all gone

                  end if                                !on to better things
                  !finished with ice
               end if
               !
            end if                                   !finished with the SD decision
            !
            return
            !
         end subroutine EVAPORATE

         !
         ! *********************************************************************
         subroutine CONVERT(t, l, qc, rho, qh, dt)
            !
            !- ACCRETION AND AUTOCONVERSION
            !

            ! plumegen_coms
            integer, intent(IN)    :: l
            real, intent(IN)    :: dt
!!$       REAL,    INTENT(IN)    :: t  (nkp)
!!$       REAL,    INTENT(INOUT) :: qc (nkp)
!!$       REAL,    INTENT(IN)    :: rho(nkp)
!!$       REAL,    INTENT(INOUT) :: qh (nkp)
            real, intent(IN)    :: t(:)
            real, intent(INOUT) :: qc(:)
            real, intent(IN)    :: rho(:)
            real, intent(INOUT) :: qh(:)

            real, parameter ::  AK1 = 0.001    !conversion rate constant
            real, parameter ::  AK2 = 0.0052   !collection (accretion) rate
            real, parameter ::  TH = 0.5      !Kessler threshold
            integer, parameter ::  iconv = 1      !- Kessler conversion (=0)
            !real,   parameter :: ANBASE =  50.!*1.e+6   !Berry-number at cloud base #/m^3(maritime)
            real, parameter :: ANBASE = 100000.!*1.e+6 !Berry-number at cloud base #/m^3(continental)
            !real,   parameter :: BDISP = 0.366          !Berry--size dispersion (maritime)
            real, parameter :: BDISP = 0.146          !Berry--size dispersion (continental)
            real, parameter :: TFREEZE = 269.3        !ice formation temperature
            !
            real               :: accrete, con, q, h, total

            if (T(L) .le. TFREEZE) return  !process not allowed above ice
            !
            if (QC(L) .eq. 0.) return

            ACCRETE = 0.
            CON = 0.
            Q = RHO(L)*QC(L)
            H = RHO(L)*QH(L)
            !
            !     selection rules
            !
            !
            if (QH(L) .gt. 0.) ACCRETE = AK2*Q*(H**.875)  !accretion, Kessler
            !
            if (ICONV .ne. 0) then   !select Berry or Kessler
               !
               !old   BC1 = 120.
               !old   BC2 = .0266 * ANBASE * 60.
               !old   CON = BDISP * Q * Q * Q / (BC1 * Q * BDISP + BC2)

               CON = Q*Q*Q*BDISP/(60.*(5.*Q*BDISP + 0.0366*ANBASE))
               !
            else
               !
               !   CON = AK1 * (Q - TH)   !Kessler autoconversion rate
               !
               !   IF (CON.LT.0.0) CON = 0.0   !havent reached threshold

               CON = max(0., AK1*(Q - TH)) ! versao otimizada
               !
            end if
            !
            !
            TOTAL = (CON + ACCRETE)*DT/RHO(L)

            !
            if (TOTAL .lt. QC(L)) then
               !
               QC(L) = QC(L) - TOTAL
               QH(L) = QH(L) + TOTAL    !no phase change involved
               return
               !
            else
               !
               QH(L) = QH(L) + QC(L)    !uses all there is
               QC(L) = 0.0
               !
            end if
            !
            return
            !
         end subroutine CONVERT

         ! ice - effect on temperature
         !      TTD = 0.0
         !      TTE = 0.0
         !       CALL ICE(QSATW,QSATE,Y(1),Y(2),Y(3), &
         !               TTA,TTB,TTC,DZ,ROH,D,C,TTD,TTE)
         !       DYDX(1) = DYDX(1) + TTD  + TTE ! DT/DZ on Temp
         !
         !**********************************************************************
         !
         subroutine SUBLIMATE(t, l, qv, qsat, rho, qi, est, dt)
            !
            ! ********************* VAPOR TO ICE (USE EQUATION OT22)***************

            ! plumegen_coms
            integer, intent(IN)    :: l
            real, intent(IN)    :: dt
!!$       REAL,    INTENT(INOUT) :: t   (nkp)
!!$       REAL,    INTENT(INOUT) :: qv  (nkp)
!!$       REAL,    INTENT(IN)    :: qsat(nkp)
!!$       REAL,    INTENT(IN)    :: rho (nkp)
!!$       REAL,    INTENT(INOUT) :: qi  (nkp)
!!$       REAL,    INTENT(IN)    :: est (nkp)
            real, intent(INOUT) :: t(:)
            real, intent(INOUT) :: qv(:)
            real, intent(IN)    :: qsat(:)
            real, intent(IN)    :: rho(:)
            real, intent(INOUT) :: qi(:)
            real, intent(IN)    :: est(:)

            !
            real, parameter :: EPS = 0.622, HEATFUS = 334., HEATSUBL = 2834., CP = 1.004
            real, parameter :: SRC = HEATSUBL/CP, FRC = HEATFUS/CP, TMELT = 273.3
            real, parameter :: TFREEZE = 269.3

            real            :: dtsubh, dividend, divisor, subl
            !
            DTSUBH = 0.
            !
            !selection criteria for sublimation
            if (T(L) .gt. TFREEZE) return
            if (QV(L) .le. QSAT(L)) return
            !
            !     from (OT); correction factors for units applied
            !
            DIVIDEND = ((1.e6/RHO(L))**0.475)*(QV(L)/QSAT(L) &
                                               - 1)*(QI(L)**0.525)*1.13
            DIVISOR = 7.e5 + 4.1e6/(10.*EST(L))
            !

            DTSUBH = abs(DIVIDEND/DIVISOR)   !sublimation rate
            SUBL = DTSUBH*DT                  !and amount possible
            !
            !     again check the possibilities
            !
            if (SUBL .lt. QV(L)) then
               !
               QV(L) = QV(L) - SUBL             !lose vapor
               QI(L) = QI(L) + SUBL                !gain ice
               T(L) = T(L) + SUBL*SRC         !energy change, warms air

               !print*,'5',l,qi(l),SUBL

               return
               !
            else
               !
               QI(L) = QV(L)                    !use what there is
               T(L) = T(L) + QV(L)*SRC      !warm the air
               QV(L) = 0.0
               !print*,'6',l,qi(l)
               !
            end if
            !
            return
         end subroutine SUBLIMATE

         !
         ! *********************************************************************
         !
         subroutine GLACIATE(qh, l, qv, qsat, t, dt, qi)
            !
            ! *********************** CONVERSION OF RAIN TO ICE *******************
            !     uses equation OT 16, simplest. correction from W not applied, but
            !     vapor pressure differences are supplied.
            !

            ! plumegen_coms
            integer, intent(IN)    :: l
            real, intent(IN)    :: dt
!!$       REAL   , INTENT(INOUT) :: qh  (nkp)
!!$       REAL   , INTENT(IN)    :: qv  (nkp)
!!$       REAL   , INTENT(IN)    :: qsat(nkp)
!!$       REAL   , INTENT(INOUT) :: t   (nkp)
!!$       REAL   , INTENT(INOUT) :: qi  (nkp)
            real, intent(INOUT) :: qh(:)
            real, intent(IN)    :: qv(:)
            real, intent(IN)    :: qsat(:)
            real, intent(INOUT) :: t(:)
            real, intent(INOUT) :: qi(:)

            !
            real, parameter :: HEATFUS = 334., CP = 1.004, EPS = 0.622, HEATSUBL = 2834.
            real, parameter :: FRC = HEATFUS/CP, FRS = HEATSUBL/CP, TFREEZE = 269.3
            real, parameter :: GLCONST = 0.025   !glaciation time constant, 1/sec
            real            :: dfrzh
            !

            DFRZH = 0.    !rate of mass gain in ice
            !
            !selection rules for glaciation
            if (QH(L) .le. 0.) return
            if (QV(L) .lt. QSAT(L)) return
            if (T(L) .gt. TFREEZE) return
            !
            !      NT=TMELT-T(L)
            !      IF (NT.GT.50) NT=50
            !

            DFRZH = DT*GLCONST*QH(L)    ! from OT(16)
            !
            if (DFRZH .lt. QH(L)) then
               !
               QI(L) = QI(L) + DFRZH
               QH(L) = QH(L) - DFRZH
               T(L) = T(L) + FRC*DFRZH  !warms air

               !print*,'7',l,qi(l),DFRZH

               return
               !
            else
               !
               QI(L) = QI(L) + QH(L)
               T(L) = T(L) + FRC*QH(L)
               QH(L) = 0.0

               !print*,'8',l,qi(l), QH (L)
               !
            end if
            !
            return
            !
         end subroutine GLACIATE

         !
         !
         ! *********************************************************************
         subroutine MELT(qi, l, t, dt, rho, cvi, qh)
            !
            ! ******************* MAKES WATER OUT OF ICE **************************

            ! plumegen_coms
            integer, intent(IN)    :: l
            real, intent(IN)    :: dt
!!$       REAL   , INTENT(INOUT) :: qi (nkp)
!!$       REAL   , INTENT(INOUT) :: t  (nkp)
!!$       REAL   , INTENT(IN)    :: rho(nkp)
!!$       REAL   , INTENT(IN)    :: cvi(nkp)
!!$       REAL   , INTENT(INOUT) :: qh (nkp)
            real, intent(INOUT) :: qi(:)
            real, intent(INOUT) :: t(:)
            real, intent(IN)    :: rho(:)
            real, intent(IN)    :: cvi(:)
            real, intent(INOUT) :: qh(:)

            !
            real, parameter :: FRC = 332.27, TMELT = 273., F0 = 0.75   !ice velocity factor
            real            :: DTMELT
            !
            DTMELT = 0.   !conversion,ice to rain
            !
            !selection rules
            if (QI(L) .le. 0.0) return
            if (T(L) .lt. TMELT) return
            !
            !OT(23,24)
            DTMELT = DT*(2.27/RHO(L))*CVI(L)*(T(L) - TMELT)*((RHO(L) &
                                                              *QI(L)*1.e-6)**0.525)*(F0**(-0.42))
            !after Mason,1956
            !
            !     check the possibilities
            !
            if (DTMELT .lt. QI(L)) then
               !
               QH(L) = QH(L) + DTMELT
               QI(L) = QI(L) - DTMELT
               T(L) = T(L) - FRC*DTMELT     !cools air
               !print*,'9',l,qi(l),DTMELT

               return
               !
            else
               !
               QH(L) = QH(L) + QI(L)   !get all there is to get
               T(L) = T(L) - FRC*QI(L)
               QI(L) = 0.0
               !print*,'10',l,qi(l)
               !
            end if
            !
            return
            !
         end subroutine MELT

         !-------------------------------------------------------------------------
         function ESAT_PR(TEM)

            real, intent(IN) :: TEM

            !
            ! ******* Vapor Pressure  A.L. Buck JAM V.20 p.1527. (1981) ***********
            !
            real, parameter :: CI1 = 6.1115, CI2 = 22.542, CI3 = 273.48
            real, parameter :: CW1 = 6.1121, CW2 = 18.729, CW3 = 257.87, CW4 = 227.3
            real, parameter :: TMELT = 273.3

            real            :: ESAT_PR, temc, esatm
            !
            !     formulae from Buck, A.L., JAM 20,1527-1532
            !     custom takes esat wrt water always. formula for h2o only
            !     good to -40C so:
            !
            !
            TEMC = TEM - TMELT
            if (TEMC .gt. -40.0) then
               ESATM = CW1*exp(((CW2 - (TEMC/CW4))*TEMC)/(TEMC + CW3))
               ESAT_PR = ESATM/10.        !kPa
            else
               ESATM = CI1*exp(CI2*TEMC/(TEMC + CI3))  !ice, millibars
               ESAT_PR = ESATM/10.        !kPa
            end if

            return

         end function ESAT_PR

         end module mod_chem_plumerise_scalar
