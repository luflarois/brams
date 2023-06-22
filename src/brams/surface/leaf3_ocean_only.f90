!---- this routine is prepared only for 1 grid -----
!
!----------------------------------------------------------------------
!- This routine calculates surface fluxes from the ocean using the
!  approach of the LEAF3 suface scheme. It intends to overwrite the fluxes
!  calculated by the JULES scheme, if desired.
!----------------------------------------------------------------------

module mod_leaf3_ocean_only

  use mem_grid, only:   ngrid,ngrids
  use cuparm_grell3,    only : g3d_g
  use modConvParGF, only : use_gustiness

  real, allocatable, dimension(:,:,:) :: can_temp, can_rvap, ustar  ,tstar  ,rstar
  real, allocatable, dimension(:,:)   :: sflux_t ,sflux_r

  logical :: firsttime = .true.  


contains

!*****************************************************************************

subroutine alloc_ocean_only(m2,m3)

        integer, intent(in) :: m2
        !# number of points in x direction
        integer, intent(in) :: m3
        !# number of points in y direction
        
	if(ngrids > 1) stop "alloc_ocean_only routine is prepared only for 1 grid" 
	
	allocate(can_temp(m2,m3,1), can_rvap(m2,m3,1), ustar  (m2,m3,1)  &
                ,tstar   (m2,m3,1), rstar   (m2,m3,1), sflux_t(m2,m3)    &
		,sflux_r (m2,m3),  STAT=ierr)

        if (ierr/=0) call fatal_error("Allocating ocean_only")
        can_temp=0.0; can_rvap=0.0; ustar=0.0 ;tstar=0.0 ;rstar=0.0
        sflux_t=0.0;  sflux_r=0.0

end subroutine alloc_ocean_only
!*****************************************************************************

subroutine sfclyr_ocean_only(mzp,mxp,myp,ia,iz,ja,jz,ibcon)

  use mem_all
  implicit none

  !Arguments:
  integer, intent(in) :: mzp,mxp,myp,ia,iz,ja,jz,ibcon

  !Local Variables
  real :: rslif
  integer :: ng

  real, dimension(mxp,myp) :: l_ths2, l_rvs2, l_pis2, l_dens2, &
       l_ups2, l_vps2, l_zts2

  ! Interface necessary to use pointer as argument - TEB
  interface
     subroutine leaf3_ocean_only(m1,m2,m3,mzg,mzs,np,ia,iz,ja,jz    &
          ,leaf,basic,turb,radiate,grid,cuparm,micro     &
          ,ths2,rvs2,pis2,dens2,ups2,vps2,zts2           &
          !
          )
       ! for interface
       use mem_leaf, only: leaf_vars  ! type
       use mem_basic, only: basic_vars! type
       use mem_turb, only: turb_vars! type
       use mem_radiate, only: radiate_vars! type
       use mem_grid, only: grid_vars! type
       use mem_cuparm, only: cuparm_vars! type
       use mem_micro, only: micro_vars! type
       !
       integer, intent(in) :: m1,m2,m3,mzg,mzs,np,ia,iz,ja,jz
       type (leaf_vars)    :: leaf
       type (basic_vars)   :: basic
       type (turb_vars)    :: turb
       type (radiate_vars) :: radiate
       type (grid_vars)    :: grid
       type (cuparm_vars)  :: cuparm
       type (micro_vars)   :: micro
       real, dimension(m2,m3), intent(out) :: ths2,rvs2,pis2,dens2,ups2,vps2,zts2
     end subroutine leaf3_ocean_only
  end interface

  if (nstbot == 0) return

  ng=ngrid


  call sub_leaf3_ocean_only(mzp,mxp,myp,nzg,nzs,npatch,ia,iz,ja,jz             &
       ,leaf_g (ng), basic_g (ng), turb_g (ng), radiate_g(ng)   &
       ,grid_g (ng), cuparm_g(ng), micro_g(ng)                  &
       ,l_ths2(1,1), l_rvs2(1,1), l_pis2(1,1)                   &
       ,l_dens2(1,1),l_ups2(1,1), l_vps2(1,1)                   &
       ,l_zts2(1,1)                                             &
       !
       )

  ! Apply lateral boundary conditions to leaf3 arrays

  call leaf_bcond(mxp,myp,nzg,nzs,1,jdim     &
       ,leaf_g(ng)%soil_water (1,1,1,1) ,leaf_g(ng)%sfcwater_mass  (1,1,1,1)  &
       ,leaf_g(ng)%soil_energy(1,1,1,1) ,leaf_g(ng)%sfcwater_energy(1,1,1,1)  &
       ,leaf_g(ng)%soil_text  (1,1,1,1) ,leaf_g(ng)%sfcwater_depth (1,1,1,1)  &
       ,leaf_g(ng)%ustar        (1,1,1) ,leaf_g(ng)%tstar            (1,1,1)  &
       ,leaf_g(ng)%rstar        (1,1,1) ,leaf_g(ng)%veg_albedo       (1,1,1)  &
       ,leaf_g(ng)%veg_fracarea (1,1,1) ,leaf_g(ng)%veg_lai          (1,1,1)  &
       ,leaf_g(ng)%veg_tai      (1,1,1)                                       &
       ,leaf_g(ng)%veg_rough    (1,1,1) ,leaf_g(ng)%veg_height       (1,1,1)  &
       ,leaf_g(ng)%patch_area   (1,1,1) ,leaf_g(ng)%patch_rough      (1,1,1)  &
       ,leaf_g(ng)%patch_wetind (1,1,1) ,leaf_g(ng)%leaf_class       (1,1,1)  &
       ,leaf_g(ng)%soil_rough   (1,1,1) ,leaf_g(ng)%sfcwater_nlev    (1,1,1)  &
       ,leaf_g(ng)%stom_resist  (1,1,1) ,leaf_g(ng)%ground_rsat      (1,1,1)  &
       ,leaf_g(ng)%ground_rvap  (1,1,1) ,leaf_g(ng)%veg_water        (1,1,1)  &
       ,leaf_g(ng)%veg_temp     (1,1,1) ,leaf_g(ng)%can_rvap         (1,1,1)  &
       ,leaf_g(ng)%can_temp     (1,1,1) ,leaf_g(ng)%veg_ndvip        (1,1,1)  &
       ,leaf_g(ng)%veg_ndvic    (1,1,1) ,leaf_g(ng)%veg_ndvif        (1,1,1)  )

  return
end subroutine sfclyr_ocean_only

!*****************************************************************************

subroutine sub_leaf3_ocean_only(m1,m2,m3,mzg,mzs,np,ia,iz,ja,jz  &
     ,leaf,basic,turb,radiate,grid,cuparm,micro     &
     ,ths2,rvs2,pis2,dens2,ups2,vps2,zts2           &
     )

  use mem_leaf, only: leaf_vars, isfcl, dthcon, drtcon, pctlcon
  use mem_basic, only: basic_vars
  use mem_turb, only: turb_vars
  use mem_radiate, only: radiate_vars, ilwrtyp, iswrtyp
  use mem_grid, only: grid_vars, zt, zm, ngrid, time, dtlt, if_adap, jdim, dzt, istp, itime1
  use mem_micro, only: micro_vars
  use mem_cuparm, only: cuparm_vars, nnqparm
  use micphys, only: level
  use io_params, only: iupdsst, ssttime1, ssttime2
  use leaf_coms, only: timefac_sst, niter_leaf, niter_can, dtll, dtll_factor, &
       dtlc_factor, dtlc, hcapcan, wcapcan, hcapveg, dtllohcc, dtllowcc, &
       dtlcohcc, dtlcowcc, dtlcohcv, z0fac_water, snowrough, ups, vps, ths, rvs, &
       zts, pis, dens, prss, vels, gzotheta, tempk, snowfac, thetacan, &
       vels_pat, rdi, pcpgl,fracliq
  use rconstants, only: g, cpor, p00, vonk, cp, cpi, rocp, alvl,stefan
  use node_mod, only : mynum  ! INTENT(IN)

  use ccatt_start, only: &
       ccatt


  implicit none

  ! Arguments:
  integer, intent(in) :: m1,m2,m3,mzg,mzs,np,ia,iz,ja,jz
  type (leaf_vars)    :: leaf
  type (basic_vars)   :: basic
  type (turb_vars)    :: turb
  type (radiate_vars) :: radiate
  type (grid_vars)    :: grid
  type (cuparm_vars)  :: cuparm
  type (micro_vars)   :: micro
  real, dimension(m2,m3), intent(out) :: ths2,rvs2,pis2,dens2,ups2,vps2,zts2

  ! Local variables:
  integer :: i,j,ip,iter_leaf
  real :: rslif
  integer :: k2
  real :: alb,dvelu,dvelv,velnew,sflux_uv,cosine1,sine1,cosz,patch_area,albedt,l_area &
         ,gust,zws,zkhvfl,pahfs,pqhfl,temps,pgeoh,hpbl,vels2

  real, dimension(m2,m3)    :: sflux_u, sflux_v, sflux_w & !sflux_t, sflux_r, O_albedt &
                             , O_rlongup, O_albedt
  real, dimension(m2,m3,1) :: patch_rough, ground_rsat

  real, parameter :: beta= 1. ,min_ocean=0.1

  INTERFACE

     subroutine sfclmcv(i,j,ustar, tstar, rstar, vels, vels_pat, ups, vps, &
          gzotheta, patch_area, sflux_u, sflux_v, sflux_w, sflux_t, sflux_r)
       real, intent(in)    :: ustar, tstar, rstar, vels,vels_pat, ups, vps, &
            gzotheta, patch_area
       real, intent(inout) :: sflux_u, sflux_v, sflux_w, sflux_t, sflux_r
       integer, intent(in) :: i,j
     end subroutine sfclmcv

  END INTERFACE


  if(firsttime) then 
    firsttime = .FALSE. 
    !allocate(can_temp(m2,m3,1), can_rvap(m2,m3,1), ustar      (m2,m3,1)  &
    !           ,tstar(m2,m3,1), rstar   (m2,m3,1))
    call alloc_ocean_only(m2,m3)
    can_temp(:,:,1) = leaf%can_temp(:,:,1) 
    can_rvap(:,:,1) = leaf%can_rvap(:,:,1)
    ustar   (:,:,1) = leaf%ustar   (:,:,1)
    tstar   (:,:,1) = leaf%tstar   (:,:,1) 	          	       
    rstar   (:,:,1) = leaf%rstar   (:,:,1)
    
    sflux_u  = 0.
    sflux_v  = 0.
    sflux_w  = 0.
    O_albedt  = 0.
    O_rlongup = 0.
  endif


  ! Time interpolation factor for updating SST

  if (iupdsst == 0) then
     timefac_sst = 0.
  else
     timefac_sst = (time - ssttime1(ngrid)) / (ssttime2(ngrid) - ssttime1(ngrid))
  endif

  ! Define leaf3 and canopy time-split timesteps here.  This ensures that leaf3
  ! will not use a timestep longer than about 40 seconds, and canopy will not
  ! use a timestep longer than about 15 seconds.  This allows values of
  ! hcapcan = 2.e4, wcapcan = 2.e1, and hcapveg = 3.e4 as are now defined below.

  niter_leaf  = max(1,nint(dtlt/40.+.4))
  niter_can   = max(1,nint(dtll/15.+.4))
  dtll_factor = 1. / float(niter_leaf)
  dtll        = dtlt * dtll_factor
  dtlc_factor = 1. / float(niter_can)
  dtlc        = dtll * dtlc_factor

  hcapcan = 2.0e4
  wcapcan = 2.0e1
  hcapveg = 3.e4

  dtllohcc = dtll / hcapcan
  dtllowcc = dtll / wcapcan
  dtlcohcc = dtlc / hcapcan
  dtlcowcc = dtlc / wcapcan
  dtlcohcv = dtlc / hcapveg

  z0fac_water = .016 / g
  snowrough = .001
  
  hpbl = 500. ! typical PBL height over the oceans (m)

  ! Copy surface atmospheric variables into 2d arrays for input to leaf
  call sfc_fields(m1,m2,m3,ia,iz,ja,jz,jdim                       &
          ,basic%theta(1,1,1) ,basic%rv (1,1,1) ,basic%up(1,1,1)  &
          ,basic%vp   (1,1,1) ,basic%dn0(1,1,1) ,basic%pp(1,1,1)  &
          ,basic%pi0  (1,1,1) ,grid%rtgt(1,1)   ,zt               &
          ,ths2,rvs2,ups2,vps2,pis2,dens2,zts2                    )

  do j = ja,jz
     do i = ia,iz

        if(leaf%patch_area(i,j,1) < min_ocean) cycle 
       
        ! Copy surface variables to single-column values

        ups = ups2(i,j)
        vps = vps2(i,j)
        ths = ths2(i,j)
        rvs = rvs2(i,j)
        zts = zts2(i,j)
        pis = pis2(i,j)
        dens = dens2(i,j)

	gzotheta = g * zts / ths
        prss = pis ** cpor * p00
        vels = sqrt(ups ** 2 + vps ** 2)
        temps = ths * pis 
	
        !--- downdraft mass flux for the gustiness parameterization
	!--- based on Redelsperger et al (2000).
	gust = g3d_g(ngrid)%xmb_deep_dd(i,j)
        gust = min( 0.1, max(0.0, gust ) )
	gust = log(1. + 600.4*gust -4375.*gust**2) 
	
	!--- add the gustiness to the grid scale wind
	vels2 = vels**2 + use_gustiness*(gust**2)

        !--- get the convective-scale velocity w*
        !--- local le and h fluxes for W*
        pahfs=-sflux_t(i,j) *dens*1004.64  !W/m^2
        pqhfl=-sflux_r(i,j)                !kg/m^2/s

        !--- buoyancy flux (h+le)
        zkhvfl= (pahfs/1004.64+0.608*temps*pqhfl)/dens ! K m s-1
        pgeoh = hpbl*g
      
        !--- convective-scale velocity W* (m/s)
        zws = max(0.,0.001-1.5*0.41*zkhvfl*pgeoh/temps)! m+3 s-3
        zws = 1.2*zws**.3333 ! m/s

        !---add enhancement by boundary layer convection
	!--- based on Redelsperger et al (2000) with beta ~1.
        vels2 = vels2 + (beta*zws)**2 ! m^2/s^2
	
	!--- get the final vels2 for the fluxes
	vels = sqrt(vels2)
	
        ! Update water internal energy from time-dependent SST
        leaf%soil_energy(mzg,i,j,1) = 334000.  &
             + 4186. * (leaf%seatp(i,j) + (leaf%seatf(i,j) - leaf%seatp(i,j))  &
             * timefac_sst - 273.15)

        !---reset the fluxes 
	sflux_u  (i,j) = 0.
        sflux_v  (i,j) = 0.
        sflux_w  (i,j) = 0.
        sflux_t  (i,j) = 0.
        sflux_r  (i,j) = 0.
        O_albedt (i,j) = 0.
        O_rlongup(i,j) = 0.
	
        ! Begin patch loop
        do ip = 1,1

           !if(mynum == 1) print*,"1ocean only",patch_area(i,j,ip), sflux_t(i,j)

           ! Begin leaf small timestep here.
           do iter_leaf = 1,niter_leaf

              ! Calculate radiative fluxes between atmosphere, vegetation, and ground/snow
              ! based on already-computed downward shortwave and longwave fluxes from
              ! the atmosphere.  Fill tempk array with soil and snow temperature (C) and
              ! fracliq array with liquid fraction of water content in soil and snow.
              ! Other snowcover properties are also computed here.
              cosz       = radiate%cosz  (i,j)  
	      patch_area = leaf%patch_area(i,j,ip) 
              l_area     = 1.-patch_area
	      alb        = 0.
	      
	      if (cosz>0.03) then
                  alb    = min(max(-.0139 + .0467*tan(acos(cosz)), 0.03), 0.999)
              endif

	      O_albedt(i,j)= O_albedt(i,j) + patch_area*alb
              
	      call qtk(leaf%soil_energy(mzg,i,j,ip), tempk(mzg), fracliq(mzg))
              
	      O_rlongup(i,j) = O_rlongup(i,j) + stefan*tempk(mzg)**4


              ! For water surface (patch 1), compute surface saturation mixing ratio
              ! and roughness length based on previous ustar.
              ! For soil patches, compute roughness length based on vegetation and snow.

               ground_rsat(i,j,ip) = rslif(prss,tempk(mzg))
               patch_rough(i,j,ip) = max(z0fac_water * ustar(i,j,ip) ** 2,.0001)


              ! Calculate turbulent fluxes between atmosphere and canopy (or "canopy")

              thetacan = can_temp(i,j,ip) / pis

              call stars(ustar(i,j,ip),tstar(i,j,ip)                 &
                   ,rstar(i,j,ip),ths,rvs,thetacan,can_rvap(i,j,ip)  &
                   ,zts,patch_rough(i,j,ip),leaf%patch_area(i,j,ip)       &
                   ,vels,vels_pat,vonk,dtllohcc,dens,dtll                      &
                   ,leaf%R_aer(i,j,ip))

              call sfclmcv(i,j,ustar(i,j,ip), tstar(i,j,ip),        &
                      rstar(i,j,ip), vels, vels_pat, ups, vps, gzotheta, &
                      leaf%patch_area(i,j,ip), sflux_u(i,j),		      &
                      sflux_v(i,j), sflux_w(i,j),		              &        
                      sflux_t(i,j), sflux_r(i,j)		              &
                      )

              ! For water patches, update temperature and moisture of "canopy" from
              ! divergence of fluxes with water surface and atmosphere.  rdi = ustar/5
              ! is the viscous sublayer conductivity from Garratt (1992).

              rdi = .2 * ustar(i,j,1)

              can_temp(i,j,1) = can_temp(i,j,1)        &
                      + dtllohcc * dens * cp                          &
                      * ((tempk(mzg) - can_temp(i,j,1)) * rdi    &
                      + ustar(i,j,1) * tstar(i,j,1) * pis)

              can_rvap(i,j,1) = can_rvap(i,j,1) + dtllowcc * dens       &
                      * ((ground_rsat(i,j,1) - can_rvap(i,j,1)) * rdi  &
                      + ustar(i,j,1) * rstar(i,j,1))

               if (CCATT==1) then
                    leaf%R_aer(i,j,ip) = rdi
               endif

           enddo
        enddo

     enddo
  enddo

  ! Normalize accumulated fluxes and albedo seen by atmosphere over model
  ! timestep dtlt.
  do j = ja,jz
     do i = ia,iz
        sflux_u(i,j) = sflux_u(i,j) * dtll_factor * dens2(i,j)
        sflux_v(i,j) = sflux_v(i,j) * dtll_factor * dens2(i,j)
        sflux_w(i,j) = sflux_w(i,j) * dtll_factor * dens2(i,j)
        sflux_t(i,j) = sflux_t(i,j) * dtll_factor * dens2(i,j)
        sflux_r(i,j) = sflux_r(i,j) * dtll_factor * dens2(i,j)
     enddo
  enddo
!--- combine the fluxes over the land from JULES + over the ocean calculated here 
  do j = ja,jz
     do i = ia,iz  
        if(leaf%patch_area(i,j,1) < min_ocean) cycle 
        patch_area = leaf%patch_area(i,j,1)
	l_area     = 1.-patch_area
        turb%sflux_u(i,j) = l_area*turb%sflux_u(i,j) + patch_area * sflux_u(i,j) 
        turb%sflux_v(i,j) = l_area*turb%sflux_v(i,j) + patch_area * sflux_v(i,j) 
        turb%sflux_w(i,j) = l_area*turb%sflux_w(i,j) + patch_area * sflux_w(i,j) 
        turb%sflux_t(i,j) = l_area*turb%sflux_t(i,j) + patch_area * sflux_t(i,j) 
        turb%sflux_r(i,j) = l_area*turb%sflux_r(i,j) + patch_area * sflux_r(i,j) 
     enddo
  enddo

  if (ilwrtyp > 0 .or. iswrtyp > 0) then
     do j = ja,jz
        do i = ia,iz
           O_albedt (i,j) = O_albedt (i,j) * dtll_factor
           O_rlongup(i,j) = O_rlongup(i,j) * dtll_factor
        enddo
     enddo

     do j = ja,jz
        do i = ia,iz
           if(leaf%patch_area(i,j,1) < min_ocean) cycle 
           patch_area = leaf%patch_area(i,j,1)
	   l_area     = 1.-patch_area
           radiate%albedt (i,j) = l_area*radiate%albedt (i,j)  + patch_area * O_albedt (i,j) 
           radiate%rlongup(i,j) = l_area*radiate%rlongup(i,j)  + patch_area * O_rlongup(i,j)
        enddo
     enddo
  endif

!if(mynum == 1) print*,"2ocean only"
end subroutine sub_leaf3_ocean_only

!###########################################################################
end module mod_leaf3_ocean_only
