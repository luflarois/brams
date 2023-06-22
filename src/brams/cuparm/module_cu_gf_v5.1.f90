MODULE module_cu_gf2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     This convective parameterization is build to attempt     !
!     a smooth transition to cloud resolving scales as proposed!
!     by Arakawa et al (2011, ACP). It currently does not use  !
!     subsidencespreading as in G3. Difference and details     !
!     will be described in a forthcoming paper by              !
!     Grell and Freitas (2013). The parameterization also      !
!     offers options to couple with aerosols. Both, the smooth !
!     transition part as well as the aerosol coupling are      !
!     experimental. While the smooth transition part is turned !
!     on, nd has been tested dow to a resolution of about 3km  !
!     the aerosol coupling is turned off.                      !
!                                                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   USE module_gate

   IMPLICIT NONE
   REAL   , PARAMETER,PRIVATE:: c1=.001

CONTAINS

   SUBROUTINE GFDRV2(mgmxp,mgmyp,mgmzp,ngrid,ngrids_cp,iens     &
              ,mynum,i0,j0,time,mzp,mxp,myp                     &
              ,DT                                               &
              ,DX                                               &
              ,autoconv                                         & !
              ,aeroevap                                         & !
              ,rho                                              &
              ,PRATEC                                           &
              ,U                                                &
              ,V                                                &
              ,theta                                            &
              ,thetail                                          &
              ,pp                                               &
              ,pi0                                              &
              ,W                                                &
              ,rv                                               &
              ,rtp                                              &
              ,rtgt                                             &
              ,pt                                               &
              ,XLV                                              &
              ,CP                                               &
              ,G                                                &
              ,r_v                                              &
              ,p00                                              &
              ,cpor                                             &
              ,rgas                                             &
              ,zm                                               &
              ,zt                                               &
              ,APR_GR                                           &
              ,APR_W                                            &
              ,APR_MC                                           &
              ,APR_ST                                           &
              ,APR_AS                                           &
              ,MASS_FLUX                                        &
              ,err_deep                                         &
              ,xmb_shallow                                      &
              ,weight_GR                                        &
              ,weight_W                                         &
              ,weight_MC                                        &
              ,weight_ST                                        &
              ,weight_AS                                        &
              ,training                                         &
              ,topt                                             &
              ,patch_area                                       &
              ,npat                                             &
              ,rshort                                           &
              ,cugd_avedx                                       &
              ,imomentum                                        &
              ,ensdim,maxiens,maxens,maxens2,maxens3,ichoice    &
              ,ishallow_gf                                      &
              ,ids,ide, jds,jde, kds,kde                        &
              ,ims,ime, jms,jme, kms,kme                        &
              ,ips,ipe, jps,jpe, kps,kpe                        &
              ,its,ite, jts,jte, kts,kte                        &
              ,RTHCUTEN                                         &
              ,RQVCUTEN                                         &
              ,RQCCUTEN                                         &
              ,RUCUTEN                                          &
              ,RVCUTEN                                          &
              ,MUP                                              &
              ,RTHFTEN                                          &
              ,RQVFTEN                                          &
              ,rthblten                                         &
              ,rqvblten                                         &
              ,level                                            &
              ,rcp                                              &
              ,aot500                                           &
              ,temp2m                                           &
              ,sflux_r        &
              ,sflux_t        &
              ,tke            &
              ,tkmin          &
              ,akmin          &
	      ,COUPL_MPHYSICS &
!- for convective transport-start
              ,ierr4d                    &
              ,jmin4d                    &
              ,kdet4d                    &
              ,k224d                     &
              ,kbcon4d                   &
              ,ktop4d                    &
              ,kpbl4d                    &
              ,kstabi4d                  &
              ,kstabm4d                  &
              ,xmb4d                     &
              ,edt4d                     &
              ,pwav4d                    &
              ,pcup5d                    &
              ,up_massentr5d             &
              ,up_massdetr5d             &
              ,dd_massentr5d             &
              ,dd_massdetr5d             &
              ,zup5d                     &
              ,zdn5d                     &
              ,prup5d                    &
              ,prdn5d                    &
              ,clwup5d                   &
              ,tup5d                     &
!- for convective transport- end
!
              ,F_QV    ,F_QC    ,F_QR    ,F_QI    ,F_QS    &
!
                                     )
!-----------------------------------------------------------------------------

   IMPLICIT NONE
      ! autoconv, 1=kessler c0, 2=berry c0
      ! aeroevap, 1=orig      , 2=ccn dependence, 3=average
      integer, parameter :: use_excess   =1
      integer, parameter :: use_excess_sh=1
      real   , parameter :: zustart_dp   =0.3
      real   , parameter :: zufinal_dp   =1.0
      real   , parameter :: zutop_dp     =0.2
      real   , parameter :: zustart_sh   =0.1
      real   , parameter :: zufinal_sh   =1.0
      real   , parameter :: zutop_sh     =0.1
      integer, parameter :: imid_gf      =0
      integer, parameter :: ideep_gf     =1
      !-- ichoice_s: 0 ensemble, 1 Wstar, 4 heat-engine or 7 BLQE
      integer, parameter :: ichoice_s    =0
      real   , parameter :: ccnclean     =250.
      integer, parameter :: ens4_spread = 3 ! max(3,cugd_avedx)
      integer, parameter :: ens4=ens4_spread*ens4_spread
      real   , parameter :: tcrit=258.
      integer, parameter :: ipr_shal =0, jpr_shal=0 &
                           ,ipr_deep =0, jpr_deep=0

      real   ::  entr_rate_sh   =1.2e-3 !1.2e-3     ! 5.0e-3
      real   ::  entr_rate_deep =7.0e-5 !1.e-4      ! 7.0e-5

      logical,parameter :: MASSCON=.false. !- new mass conserv form flag
      logical,parameter :: DICYCLE=.true.  !- diurnal cycle flag
      logical,parameter :: ENTRNEW=.TRUE. !- new entr formulation
      logical,intent(IN):: COUPL_MPHYSICS  !- coupling with cloud
                                           !- microphysics param
!-------------------------------------------------------------
   INTEGER,      INTENT(IN   ) ::                               &
                                  ids,ide, jds,jde, kds,kde,    &
                                  ims,ime, jms,jme, kms,kme,    &
                                  ips,ipe, jps,jpe, kps,kpe,    &
                                  its,ite, jts,jte, kts,kte,    &
!
                                  mgmxp,mgmyp,mgmzp,ngrid,ngrids_cp,  &
                                  iens,mynum,i0,j0,mzp,mxp,myp

   INTEGER,      INTENT(IN   ) :: cugd_avedx, &
                                  ishallow_gf,imomentum,&
                                  autoconv,             & !
                                  aeroevap
   INTEGER, INTENT (in   )              ::                      &
                       ensdim,maxiens,maxens,maxens2,maxens3,ichoice,NPAT,level&
                       ,training

   REAL,         INTENT(IN   ) :: XLV, R_v,tkmin,akmin
   REAL,         INTENT(IN   ) :: CP,G, cpor, p00,rgas
   REAL,         INTENT(IN   ) :: DT, DX, time

   INTEGER ::  ITIMESTEP=0

   REAL,  DIMENSION(kms:kme) , INTENT(IN   ) ::  zm,zt

   REAL,  DIMENSION(kms:kme ,  ims:ime , jms:jme ), INTENT(IN   ) :: &
                                                          U,    &
                                                          V,    &
                                                          W,    &
                                                          rv,   &
                                                          rtp,  &
                                                          theta   ,&
                                                          thetail ,&
                                                          pp      ,&
                                                          pi0     ,&
                                                          rho     ,&
                                                          pt      ,&
                                                          rcp,tke

   REAL, DIMENSION(ims:ime , jms:jme,npat), intent(in)  :: patch_area


   REAL, DIMENSION( ims:ime , jms:jme ),INTENT(IN) :: rshort,topt,rtgt,aot500,temp2m &
                                                      ,sflux_r,sflux_t


   REAL, DIMENSION( ims:ime , jms:jme ),  INTENT(INOUT) ::      &
                          pratec,mass_flux,                     &
                          apr_gr,apr_w,apr_mc,apr_st,apr_as,    &
                          xmb_shallow,err_deep
   REAL, DIMENSION( ims:ime , jms:jme ),   INTENT(IN) ::        &
                  weight_GR,weight_W,weight_MC,weight_ST,weight_AS

   REAL, DIMENSION( its:ite , jts:jte ) ::RAINCV,               &
                                          APR_CAPMA,APR_CAPME,APR_CAPMI

!
! Optionals
!
   REAL, DIMENSION(kms:kme , ims:ime ,  jms:jme ),              &
         OPTIONAL,                                              &
         INTENT(IN) ::      RTHFTEN,  RQVFTEN,                  &
                            RTHBLTEN,RQVBLTEN


   REAL, DIMENSION(kms:kme , ims:ime ,  jms:jme ), OPTIONAL,    &
         INTENT(INOUT) ::                           RTHCUTEN    &
                                                   ,RQVCUTEN    &
                                                   ,RQCCUTEN    &
                                                   ,RUCUTEN     &
                                                   ,RVCUTEN     &
                                                   ,MUP

!
! Flags relating to the optional tendency arrays declared above
! Models that carry the optional tendencies will provdide the
! optional arguments at compile time; these flags all the model
! to determine at run-time whether a particular tracer is in
! use or not.
!
   LOGICAL, OPTIONAL ::     F_QV      &
                           ,F_QC      &
                           ,F_QR      &
                           ,F_QI      &
                           ,F_QS
!- for convective transport-start
  integer, dimension(mgmxp,mgmyp,maxiens,ngrids_cp) ::         &
               ierr4d                       &
              ,jmin4d                       &
              ,kdet4d                       &
              ,k224d                        &
              ,kbcon4d                      &
              ,ktop4d                       &
              ,kpbl4d                       &
              ,kstabi4d                     &
              ,kstabm4d

   real,dimension(mgmxp,mgmyp,maxiens,ngrids_cp) :: &
               xmb4d                     &
              ,edt4d                     &
              ,pwav4d
   real,dimension(mgmzp,mgmxp,mgmyp,maxiens,ngrids_cp) ::&
               pcup5d                    &
              ,up_massentr5d             &
              ,up_massdetr5d             &
              ,dd_massentr5d             &
              ,dd_massdetr5d             &
              ,zup5d                     &
              ,zdn5d                     &
              ,prup5d                    &
              ,prdn5d                    &
              ,clwup5d                   &
              ,tup5d
!- for convective transport-end

!----------------------------------------------------------------------
! LOCAL VARS
!
! basic environmental input includes moisture convergence (mconv)
! omega (omeg), windspeed (us,vs), and a flag (aaeq) to turn off
! convection for this call only and at that particular gridpoint
!
     real,    dimension(its:ite,jts:jte,1:ensdim) ::  xf_ens,pr_ens
     real,    dimension(its:ite,jts:jte,1:ensdim) ::  massflni,xfi_ens,pri_ens

     real,    dimension(its:ite,jts:jte )         ::  massi_flx    &
                         ,apri_gr,apri_w,apri_mc,apri_st,apri_as   &
                         ,apri_capma,apri_capme,apri_capmi         &
                         ,gswi                                     &
                         ,iweight_gr,iweight_w,iweight_mc,iweight_st,iweight_as

     real,    dimension (its:ite,kms:kme) ::phh,phhl,phil
     real,    dimension (its:ite,kts:kte) ::                                  &
                            qcheck,zo,temp_old,qv_old,PO,P2d,US,VS,rhoi       &
			   ,temp_new_dp,qv_new_dp,temp_new_sh,qv_new_sh,z2d   &
                           ,tkeg,rcpg                                         &
                           ,subt,subq,outt,outq,outqc,subm,cupclw,cupclws,dhdt&
                           ,outts,outqs,outqcs,outu,outv

    real,    dimension (its:ite)         ::  &
         pret, ter11, xlandi,pbl             &
        ,zws,ccn,psur,aaeq,cuten,umean,vmean,pmean,xmbs,ztexec,zqexec &
        ,xmb,h_sfc_flux,le_sfc_flux,tsur

    real,    dimension (its:ite,kts:kte,1:ens4) ::  omeg

    real,    dimension (its:ite,1:ens4)         ::  mconv

    integer, dimension (its:ite)         :: kbcon, ktop,kpbli,k22s,kbcons,ktops&
                                           ,ierr,ierrs
    integer :: i,j,k,kr,n,itf,jtf,ktf

    character*50 :: ierrc(its:ite), ierrcs(its:ite)

    real :: tscl_kf,dp,dq,exner, cpdtdt, beta &
           ,PTEN,PQEN,PAPH,ZRHO,PAHFS,PQHFL,ZKHVFL,PGEOH

    !-don�t change this
    itf=ite
    ktf=kte-1
    jtf=jte
!----------------------------------------------------------------------


!-- big loop over j dimension
   DO 100 j = jts,jtf

   !-- initialization
     !- deep convection section
     DO I= its,itf
        APRi_GR(i,j)    =apr_gr(i,j)
        APRi_w(i,j)     =apr_w(i,j)
        APRi_mc(i,j)    =apr_mc(i,j)
        APRi_st(i,j)    =apr_st(i,j)
        APRi_as(i,j)    =apr_as(i,j)
        APRi_capma(i,j) =apr_capma(i,j)
        APRi_capme(i,j) =apr_capme(i,j)
        APRi_capmi(i,j) =apr_capmi(i,j)
     ENDDO
     DO I= its,itf
        ierrc(i)=" "
        xfi_ens(i,j,1:ensdim)=0.0
        pri_ens(i,j,1:ensdim)=0.0

        ierr (i) =0
        kbcon(i) =0
        ktop (i) =0

        aaeq  (i) =0.0
        pret  (i) =0.0
        umean (i) =0.0
        vmean (i) =0.0
        pmean (i) =0.0
        cuten (i) =0.0
        xmb   (i) =0.0
        ztexec(i) =0.0
        zqexec(i) =0.0
        !
        mass_flux(i,j)  =0.0
        massi_flx(i,j)  =0.0
        raincv   (i,j)  =0.0
        pratec   (i,j)  =0.0
        !- tendencies
        cupclw(i,:)=0.0
        subt  (i,:)=0.0
        subq  (i,:)=0.0
        outt  (i,:)=0.0
        outu  (i,:)=0.0
        outv  (i,:)=0.0
        outq  (i,:)=0.0
        outqc (i,:)=0.0
     ENDDO
     !- shallow convection section
     DO I= its,itf
        ierrs (i) =0
        k22s  (i) =0
        kbcons(i) =0
        ktops (i) =0
        !-mass flux at cloud base
        xmbs  (i)   =0.0
	xmb_shallow(i,j)=0.0
        !- tendencies
        outqcs (i,:)=0.0
        outts  (i,:)=0.0
        outqs  (i,:)=0.0
        cupclws(i,:)=0.0
     ENDDO
!
     DO n=1,ens4
       do i= its,itf
         mconv(i,n)=0.0
       enddo
       do k=kts,kte
        do i= its,itf
         omeg(i,k,n)=0.0
        enddo
       enddo
     ENDDO
     !
     DO k=1,ensdim
      do i= its,itf
         massflni(i,j,k)=0.0
      enddo
     ENDDO

     IF(autoconv == 2) THEN
      DO I= its,itf
             ccn(i) = max( 100., ( 370.37*(0.01+MAX(0.,aot500(i,j))))**1.555 )
      ENDDO
     ELSE
      DO I= its,itf
             ccn(i) = 100.
      ENDDO
     ENDIF
     !  put hydrostatic pressure on half levels
     !DO k=kts,ktf
     ! do i=its,itf
     !    phh(i,k) =((pp(k,i,j)+pi0(k,i,j))/cp)**cpor*p00
     ! enddo
     !ENDDO

     DO i=its,itf
         gswi(i,j) = rshort(i,j)

	 xlandi(i) = patch_area(i,j,1) !flag < 1 para land
                                       !flag  =1 para water

         psur (i)  = 0.5*( ((pp(1,i,j)+pi0(1,i,j))/cp)**cpor*p00 +  &
                           ((pp(2,i,j)+pi0(2,i,j))/cp)**cpor*p00 )*1.e-2
         tsur (i)  = temp2m(i,j)
         ter11(i)  = max(0.,topt(i,j))
     ENDDO

     DO k=kts,ktf
       DO i=its,itf
         kr=k+1
!
         !- heigths, current pressure, temp and water vapor mix ratio
	 zo      (i,k)  = zt(kr)*rtgt(i,j)+topt(i,j)
         po      (i,k)  = ((pp(kr,i,j)+pi0(kr,i,j))/cp)**cpor*p00*.01
         p2d     (i,k)  = po(i,k)
         temp_old(i,k)  = theta(kr,i,j)*(pp(kr,i,j)+pi0(kr,i,j))/cp
         qv_old  (i,k)  = max(1.e-08, rv(kr,i,j))

	 qcheck  (i,k)  = qv_old(i,k)

         !- air density, TKE and cloud liq water mixing ratio
         rhoi    (i,k)  = rho(kr,i,j)
         tkeg    (i,k)  = tke(kr,i,j)
         rcpg    (i,k)  = rcp(kr,i,j)
         !- wind velocities
         us      (i,k)  = 0.5*( u(kr,i,j) + u(kr,max(1,i-1),j) )
         vs      (i,k)  = 0.5*( v(kr,i,j) + v(kr,i,max(1,j-1)) )
         omeg    (i,k,:)= -g*rho(kr,i,j)*w(kr,i,j)

       ENDDO
     ENDDO

     !- pbl  (i) = depth of pbl layer (m)
     !- kpbli(i) = index of zo(i,k)
     call get_zi_gf(j,its,ite,kts,kte,its,itf,ktf,ierrs,kpbli,pbl,&
                  tkeg,rcpg,zo,ter11,tkmin)

     !- begin: for GATE soundings-------------------------------------------
     !- this section is intended for model developments only and must
     !- not be used for normal runs.
     if(use_gate) then
       do i=its,itf
         do k=kts,ktf
           p2d      (i,k) = ppres(jl,k)
           temp_old (i,k) = ptemp(jl,k)+273.15
           qv_old   (i,k) = pq(jl,k)/1000.
           us  (i,k) = pu(jl,k)
           vs  (i,k) = pv(jl,k)
           omeg(i,k,:)=pvervel(jl,k)
           phil(i,k) = pgeo(jl,k)*g   !geo
           po  (i,k) = p2d (i,k)
           zo  (i,k) = phil(i,k)/g    !meters
           rhoi(i,k) =1.e2*p2d(i,k)/(rgas*temp_old(i,k))
          enddo
          ter11(i)  = phil(i,1)/g  ! phil is given in g*h.
          psur(i)   = p2d(i,1)
          tsur(i)   = temp2m(i,j) !temp_old(i,1)
          gswi(i,1) = 0.0
          kpbli(i)  = 4
          zws  (i)  = 1.0 ! wstar
          do k=kts,ktf
           temp_new_dp(i,k) = temp_old(i,k) + dt *(zadvt(jl,k)+zqr(jl,k))/86400.
           qv_new_dp  (i,k) = qv_old  (i,k) + dt * zadvq(jl,k)
          enddo
          pbl (i)  = zo(i,kpbli(i))
       enddo
     endif
     !- end:   for GATE soundings-------------------------------------------

     DO I=ITS,ITF
       pten = temp_old(i,1)
       pqen = qv_old(i,1)
       paph = 100.*psur(i)
       zrho = paph/(287.04*(temp_old(i,1)*(1.+0.608*qv_old(i,1))))
       !- sensible and latent sfc flux for new heat-engine closure
       h_sfc_flux (i)=zrho*cp *sflux_t(i,j)
       le_sfc_flux(i)=zrho*xlv*sflux_r(i,j)
       !
       !- le and h fluxes
       pahfs=-sflux_t(i,j) *zrho*1004.64!w/m2
       pqhfl=-sflux_r(i,j)!(kg/m^2/s)
       !- buoyancy flux (h+le)
       zkhvfl= (pahfs/1004.64+0.608*pten*pqhfl)/zrho
       !- depth of 1st model layer
       !pgeoh = dz8w(i,1,j)*g
       !pgeoh = zt(2)*rtgt(i,j)*g
       pgeoh =  ( zo(i,1)-topt(i,j) )*g !zo(i,1)=zt(2)*rtgt(i,j)+topt(i,j)
       !-convective-scale velocity w*
       !- in the future, change 0.001 by ustar^3
       zws(i) = max(0.,0.001-1.5*0.41*zkhvfl*pgeoh/pten)

       if(zws(i) > tiny(pgeoh)) then
         !-convective-scale velocity w*
         zws(i) = 1.2*zws(i)**.3333
         !- temperature excess
         ztexec(i)     = max(0.,-1.5*pahfs/(zrho*zws(i)*1004.64))
         !- moisture  excess
         zqexec(i)     = max(0.,-1.5*pqhfl/(zrho*zws(i)))
       ! zws(i)=zws(i)*zrho

       endif   ! zws > 0
       !- zws for shallow convection closure (grant 2001)
       !- height of the pbl
       pgeoh = pbl(i)*g
       !-convective-scale velocity w*
       zws(i) = max(0.,0.001-1.5*0.41*zkhvfl*pgeoh/pten)
       zws(i) = 1.2*zws(i)**.3333
       zws(i) = zws(i)*zrho !check if zrho is correct
       !print*,"surf=",pbl (i), tsur(i)-273.,h_sfc_flux (i),le_sfc_flux(i),zws(i)
     enddo

     DO k=  kts+1,ktf-1
      DO i = its,itf
         if((p2d(i,1)-p2d(i,k)).gt.150.and.p2d(i,k).gt.300)then
            dp=-.5*(p2d(i,k+1)-p2d(i,k-1))
            umean(i)=umean(i)+us(i,k)*dp
            vmean(i)=vmean(i)+vs(i,k)*dp
            pmean(i)=pmean(i)+dp
         endif
      ENDDO
     ENDDO
     DO n=1,ens4
      DO k=kts,ktf-1
       DO I = its,itf
        dq        =  qv_old(i,k+1)-qv_old(i,k)
        mconv(i,n)=mconv(i,n)+omeg(i,k,n)*dq/g
       ENDDO
      ENDDO
     ENDDO
     DO n=1,ens4
      DO i = its,itf
        if(mconv(i,n).lt.0.)mconv(i,n)=0.
      ENDDO
     ENDDO
!
!------ CALL CUMULUS PARAMETERIZATION
!
!-- shallow convection
    IF(ishallow_gf.eq.1)THEN

       tscl_kf=dx/25.

       DO i=its,itf
          DO k=kts,ktf
           kr=k+1

           exner= pp(kr,i,j)+pi0(kr,i,j)
           cpdtdt= exner*(rthblten(kr,i,j)+rthften(kr,i,j)) + theta(kr,i,j)*pt(kr,i,j)

           temp_new_sh(i,k)=temp_old(i,k) + cpdtdt/cp*dt
           qv_new_sh  (i,k)=  qv_old(i,k) + (rqvblten(kr,i,j)+rqvften(kr,i,j))*dt

           !- all forcings change moist static energy
           !dhdt(i,k)=cpdtdt + xlv*(rqvblten(kr,i,j) + rqvften(kr,i,j)

           !- only pbl forcing changes moist static energy
           dhdt(i,k)= exner*(rthblten(kr,i,j)) + theta(kr,i,j)*pt(kr,i,j) + &
                      xlv  *(rqvblten(kr,i,j))

          ENDDO
        ENDDO
        DO i=its,ite
           aaeq (i)=0.0
           ierrs(i)=0
        ENDDO
        call CUP_gf_sh(J,COUPL_MPHYSICS,TCRIT  ,ccnclean,tscl_kf &
                      ,xlv,r_v,cp,g,ichoice_s,ipr_shal,jpr_shal &
                      ,autoconv,itf,jtf,ktf,itimestep		&
                      ,use_excess_sh,its,ite, jts,jte, kts,kte  &
	              !input data
	              ,dt,dx				     &
	              ,entr_rate_sh			     &
                      ,zustart_sh,zufinal_sh,zutop_sh	     &
	              ,h_sfc_flux,le_sfc_flux		     &
                      ,tsur				     &
	              ,psur				     &
                      ,ter11 ,xlandi,gswi		     &
	              ,ztexec				     &
	              ,zqexec				     &
	              ,ccn				     &
	              ,rhoi				     &
	              ,dhdt				     &
	              ,zws				     &
	              ,zo				     &
	              ,aaeq				     &
	              ,temp_old 			     &
	              ,qv_old				     &
                      ,temp_new_sh			     &
	              ,qv_new_sh			     &
	              ,p2d				     &
	              ,p2d				     &
	              ,us				     &
	              ,vs				     &
	              !output data			     &
	              ,ierrs,ierrcs			     &
	              ,xmbs				     &
	              ,outts				     &
	              ,outqs				     &
	              ,outqcs				     &
                      ,cupclws  			     &
	              ,kpbli,kbcons,ktops,k22s  	     &
                       !- for shallow convective transport
                       !- here iens is assumed equal to "2"
                      ,mgmxp,  mgmzp ,mynum	          &
                      ,ierr4d	     (1,j,2,ngrid)	  &
                      ,jmin4d	     (1,j,2,ngrid)	  &
                      ,kdet4d	     (1,j,2,ngrid)	  &
                      ,k224d	     (1,j,2,ngrid)	  &
                      ,kbcon4d       (1,j,2,ngrid)	  &
                      ,ktop4d	     (1,j,2,ngrid)	  &
                      ,kpbl4d	     (1,j,2,ngrid)	  &
                      ,kstabi4d      (1,j,2,ngrid)	  &
                      ,kstabm4d      (1,j,2,ngrid)	  &
                      ,xmb4d	     (1,j,2,ngrid)	  &
                      ,edt4d	     (1,j,2,ngrid)	  &
                      ,pwav4d	     (1,j,2,ngrid)	  &
                      ,pcup5d	     (1,1,j,2,ngrid)	  &
                      ,up_massentr5d (1,1,j,2,ngrid)	  &
                      ,up_massdetr5d (1,1,j,2,ngrid)	  &
                      ,dd_massentr5d (1,1,j,2,ngrid)	  &
                      ,dd_massdetr5d (1,1,j,2,ngrid)	  &
                      ,zup5d	     (1,1,j,2,ngrid)	  &
                      ,zdn5d	     (1,1,j,2,ngrid)	  &
                      ,prup5d	     (1,1,j,2,ngrid)	  &
                      ,prdn5d	     (1,1,j,2,ngrid)	  &
                      ,clwup5d       (1,1,j,2,ngrid)	  &
                      ,tup5d	     (1,1,j,2,ngrid)	  &
                                                          )
                      !- for convective transport-end

     ENDIF

     !--- deep convection
     IF(ideep_gf ==1 )then

      if(.not. use_gate) then ! GATE
       do k=kts,ktf
         do i=its,itf
         kr=k+1
         exner  = pp(kr,i,j)+pi0(kr,i,j)
         cpdtdt =  theta(kr,i,j)*pt(kr,i,j) + exner*(rthblten(kr,i,j) + rthften(kr,i,j))

         temp_new_dp(i,k)= temp_old(i,k)  +  cpdtdt/cp                           *dt !! + dt*outts(i,k)
         qv_new_dp  (i,k)= qv_old  (i,k)  +  (rqvften(kr,i,j) + rqvblten(kr,i,j))*dt !! + dt*outqs(i,k)

         temp_new_dp(i,k)=max(200.  ,temp_old  (i,k))
         qv_new_dp  (i,k)=max(1.e-08,qv_new_dp (i,k))
       enddo
      enddo
     endif
     CALL CUP_gf(j,coupl_mphysics,maxiens,maxens,maxens2,maxens3,ensdim      &
          	,its,ite, jts,jte, kts,kte, itf,jtf,ktf       &
          	,xlv,r_v,cp,g,ichoice,ipr_deep,jpr_deep,ens4  &
          	,beta,autoconv,aeroevap,training	      &
          	,imid_gf     ,itimestep,masscon   &
          	,zustart_dp,zufinal_dp,zutop_dp   &
          	,dicycle,entr_rate_deep,entrnew   &
          	,tcrit  			  &
          	,iens				  &
          	! input data
          	,dx				  &
          	,dt	      &
          	,kpbli        &
          	,use_excess   &
          	,ztexec       &
          	,zqexec       &
          	,ccn	      &
          	,ccnclean     &
          	,rhoi	      &
          	,mconv        &
          	,omeg	      &
          	,aaeq	      &
          	,temp_old     &
          	,qv_old	      &
          	,ter11,xlandi,gswi    &
          	,temp_new_dp  &
          	,qv_new_dp    &
          	,zo	      &
          	,po	      &
          	,p2d	      &
          	,psur	      &
          	,us	      &
          	,vs	      &
          	,zws	      &
          	,dhdt	      &
          	!output data
          	,ierr,ierrc   &
          	,xmb	      &
          	,subm,subt,subq &
          	,pret	      &
          	,outt	      &
          	,outq	      &
          	,outqc        &
          	,outu	      &
          	,outv	      &
          	,cupclw       &
          	,APRi_GR,APRi_W,APRi_MC,APRi_ST,APRi_AS     &
          	,APRi_CAPMA,APRi_CAPME,APRi_CAPMI	    &
          	,kbcon,ktop				    &
          	,xfi_ens,pri_ens			    &
          	!- for training
          	,iweight_GR    &
          	,iweight_W     &
          	,iweight_MC    &
          	,iweight_ST    &
          	,iweight_AS    &
          	!- for convective transport-start
          	!- here iens is assumed equal to "1"
          	,mgmxp,  mgmzp ,mynum	&
          	,ierr4d   (1,j,1,ngrid) 	&
          	,jmin4d   (1,j,1,ngrid) 	&
          	,kdet4d   (1,j,1,ngrid) 	&
          	,k224d    (1,j,1,ngrid) 	&
          	,kbcon4d  (1,j,1,ngrid) 	&
          	,ktop4d   (1,j,1,ngrid) 	&
          	,kpbl4d   (1,j,1,ngrid) 	&
          	,kstabi4d (1,j,1,ngrid) 	&
          	,kstabm4d (1,j,1,ngrid) 	&
          	,xmb4d    (1,j,1,ngrid) 	&
          	,edt4d    (1,j,1,ngrid) 	&
          	,pwav4d   (1,j,1,ngrid) 	&
          	,pcup5d        (1,1,j,1,ngrid)  &
          	,up_massentr5d (1,1,j,1,ngrid)  &
          	,up_massdetr5d (1,1,j,1,ngrid)  &
          	,dd_massentr5d (1,1,j,1,ngrid)  &
          	,dd_massdetr5d (1,1,j,1,ngrid)  &
          	,zup5d         (1,1,j,1,ngrid)  &
          	,zdn5d         (1,1,j,1,ngrid)  &
          	,prup5d        (1,1,j,1,ngrid)  &
          	,prdn5d        (1,1,j,1,ngrid)  &
          	,clwup5d       (1,1,j,1,ngrid)  &
          	,tup5d         (1,1,j,1,ngrid)  &
          	)
               !- for convective transport-end

           DO i=its,itf
             DO k=kts,ktf
              qcheck(i,k)=qv_old(i,k)+outqs(i,k)*dt
             ENDDO
           ENDDO
      CALL neg_check(j,subt,subq,dt,qcheck,outq,outt,outqc,pret,its,ite,kts,kte,itf,ktf)
     ENDIF


     !------------ feedback

     !-- deep convection
     DO I=its,itf
       cuten(i)=0.
       if(pret(i).gt.0.)then
          !--   rainfall output
          pratec(i,j)=pratec(i,j)+pret(i)
          !print*,"GF prec=",pratec(i,j)
               !
          cuten(i)=1.
       endif
     ENDDO
     !-- shallow convection
     DO I=its,itf
        if(xmbs(i).gt.0.)then
          xmb_shallow(i,j)=xmbs(i)
        endif
     ENDDO

     !-- deep + shallow convection
     DO i = its,itf
      DO k = kts,ktf-1
             kr=k+1

             !- feedback the tendencies for convection
             !- adding         shallow (s) +        deep
             RTHCUTEN(kr,i,j)= outts (i,k)+(subt(i,k)+outt(i,k))*cuten(i)
             RQVCUTEN(kr,i,j)= outqs (i,k)+(subq(i,k)+outq(i,k))*cuten(i)
             RQCCUTEN(kr,i,j)= outqcs(i,k)+          outqc(i,k) *cuten(i)

             !- use this formulation to add total water to rt tendency
             !RQVCUTEN(kr,i,j)= RQVCUTEN(kr,i,j)+RQCCUTEN(kr,i,J)
             !RQCCUTEN(kr,i,j)= 0.0

             if(imomentum == 1) then ! - only deep convection
               RUCUTEN(kr,i,j)=outu(i,k)*cuten(i)
               RVCUTEN(kr,i,j)=outv(i,k)*cuten(i)
             endif
      ENDDO
      !
      !- converting Dtemp/Dt to Dtheta/ Dt
      RTHCUTEN(2:kme,i,j)=RTHCUTEN(2:kme,i,j)*cp/(pp(2:kme,i,j) + pi0(2:kme,i,j))

      !- setting tendencies at k=1
      RTHCUTEN(1,i,j)=  RTHCUTEN(2,i,j)
      RQVCUTEN(1,i,j)=  RQVCUTEN(2,i,j)
      RQCCUTEN(1,i,J)=  RQCCUTEN(2,i,J)
      if(imomentum == 1) then
        RUCUTEN(1,i,j)=RUCUTEN(2,i,j)
        RVCUTEN(1,i,j)=RVCUTEN(2,i,j)
      endif
     ENDDO

     DO n= 1,ensdim
            DO I= its,itf
              xf_ens(i,j,n)=xfi_ens(i,j,n)
              pr_ens(i,j,n)=pri_ens(i,j,n)
            ENDDO
     ENDDO
     DO I= its,itf
               APR_GR(i,j)=apri_gr(i,j)
               APR_w(i,j) =apri_w(i,j)
               APR_mc(i,j)=apri_mc(i,j)
               APR_st(i,j)=apri_st(i,j)
               APR_as(i,j)=apri_as(i,j)
               APR_capma(i,j)=apri_capma(i,j)
               APR_capme(i,j)=apri_capme(i,j)
               APR_capmi(i,j)=apri_capmi(i,j)
               mass_flux(i,j)=xmb(i)
     ENDDO
     ! added output por updraft mass flux
     DO I=its,itf
       mup(:,i,j)=0.0
       if(pret(i).gt.0.)then
          DO k = kts,ktf-1
                 kr=k+1
             mup(kr,i,j)=xmb(i)*zup5d(k,i,j,1,ngrid) !deep convection
          ENDDO
       endif
       if(xmbs(i).gt.0.)then
          DO k = kts,ktf-1
                 kr=k+1
             mup(kr,i,j)=mup(kr,i,j)+xmbs(i)*zup5d(k,i,j,2,ngrid) !deep + shallow convection
          ENDDO
       endif
     ENDDO

 100    continue

 END SUBROUTINE GFDRV2

!========================================================================================
   SUBROUTINE CUP_gf(J,coupl_mphysics,maxiens,maxens,maxens2,maxens3,ensdim   &
     		    ,its,ite, jts,jte, kts,kte ,itf,jtf,ktf	   &
     		    ,xl,rv,cp,g,ichoice,ipr,jpr,ens4		   &
     		    ,beta,autoconv,aeroevap,training		   &
     		    ,imid   ,ktau,masscon			   &
     		    ,zustart_dp,zufinal_dp,zutop_dp		   &
     		    ,dicycle,entr_rate, entrnew 		   &
     		    ,TCRIT	       &
     		    ,iens	       &
     		    !input data
     		    ,dx 	       &
     		    ,DTIME	       &
     		    ,kpbl	       &
     		    ,use_excess        &
     		    ,ztexec	       &
     		    ,zqexec	       &
     		    ,ccn	       &
     		    ,ccnclean	       &
     		    ,rho	       &
     		    ,mconv	       &
     		    ,omeg	       &
     		    ,AAEQ	       &
     		    ,T  	       &
     		    ,Q  	       &
     		    ,Z1  ,xland,gsw    &
     		    ,TN 	       &
     		    ,QO 	       &
     		    ,zo 	       &
     		    ,PO 	       &
     		    ,P  	       &
     		    ,PSUR	       &
     		    ,US 	       &
     		    ,VS 	       &
     		    ,zws	       &
     		    ,dhdt	       &
     		    !output data
     		    ,ierr,ierrc        &
     		    ,xmb_out	       &
     		    ,sub_mas ,subt,subq&
     		    ,PRE	       &
     		    ,OUTT	       &
     		    ,OUTQ	       &
     		    ,OUTQC	       &
     		    ,outu	       &
     		    ,outv	       &
     		    ,cupclw	       &
     		    ,APR_GR,APR_W,APR_MC,APR_ST,APR_AS     &
     		    ,APR_CAPMA,APR_CAPME,APR_CAPMI	   &
     		    ,kbcon,ktop 			   &
     		    ,xf_ens,pr_ens			   &
     		    !--srf - for training
     		    ,weight_GR,weight_W,weight_MC,weight_ST,weight_AS	&
     		    !- for convective transport-start
     		    ,mgmxp,  mgmzp,mynum &
     		    ,ierr4d		 &
     		    ,jmin4d		 &
     		    ,kdet4d		 &
     		    ,k224d		 &
     		    ,kbcon4d		 &
     		    ,ktop4d		 &
     		    ,kpbl4d		 &
     		    ,kstabi4d		 &
     		    ,kstabm4d		 &
     		    ,xmb4d		 &
     		    ,edt4d		 &
     		    ,pwav4d		 &
     		    ,pcup5d		 &
     		    ,up_massentr5d	 &
     		    ,up_massdetr5d	 &
     		    ,dd_massentr5d	 &
     		    ,dd_massdetr5d	 &
     		    ,zup5d		 &
     		    ,zdn5d		 &
     		    ,prup5d		 &
     		    ,prdn5d		 &
     		    ,clwup5d		 &
     		    ,tup5d		 )
                    !- for convective transport-end

     IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
        autoconv,aeroevap,itf,jtf,ktf,ktau,training,use_excess,        &
        its,ite, jts,jte, kts,kte,ipr,jpr,ens4,imid
     integer, intent (in   )              ::                           &
        j,ensdim,maxiens,maxens,maxens2,maxens3,ichoice,iens
  !
  !
  !
     real,    dimension (its:ite,jts:jte,1:ensdim)                     &
        ,intent (inout)                   ::                           &
        xf_ens,pr_ens
     real,    dimension (its:ite,jts:jte)                              &
        ,intent (inout )                  ::                           &
               APR_GR,APR_W,APR_MC,APR_ST,APR_AS,APR_CAPMA,     &
               APR_CAPME,APR_CAPMI
!    real, dimension( its:ite , jts:jte )                               &
!          :: weight_GR,weight_W,weight_MC,weight_ST,weight_AS
     real,    dimension (its:ite,jts:jte)                              &
        ,intent (in   )                   ::                           &
               gsw
  ! outtem = output temp tendency (per s)
  ! outq   = output q tendency (per s)
  ! outqc  = output qc tendency (per s)
  ! pre    = output precip
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (inout  )                   ::                           &
        outu,outv,OUTT,OUTQ,OUTQC,subt,subq,sub_mas,cupclw
     real,    dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
        pre,xmb_out
     integer,    dimension (its:ite)                                   &
        ,intent (out  )                   ::                           &
        kbcon,ktop
     integer,    dimension (its:ite)                                   &
        ,intent (in  )                   ::                           &
        kpbl
  !
  ! basic environmental input includes moisture convergence (mconv)
  ! omega (omeg), windspeed (us,vs), and a flag (aaeq) to turn off
  ! convection for this call only and at that particular gridpoint
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        dhdt,rho,T,PO,P,US,VS,tn
     real,    dimension (its:ite,kts:kte,1:ens4)                       &
        ,intent (inout   )                   ::                           &
        omeg
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (inout)                   ::                           &
         Q,QO
     real, dimension (its:ite)                                         &
        ,intent (in   )                   ::                           &
        ccn,Z1,PSUR,xland
     real, dimension (its:ite)                                         &
        ,intent (inout   )                   ::                           &
        AAEQ
     real, dimension (its:ite)                                         &
        ,intent (inout   )                   ::                           &
        zws,ztexec,zqexec
     real, dimension (its:ite,1:ens4)                                         &
        ,intent (in   )                   ::                           &
        mconv


       real  ,intent (in   )                   ::                          &
        dx,ccnclean,dtime,tcrit,xl,cp,rv,g,zustart_dp,zufinal_dp,zutop_dp

       real  ,intent (inout   )                   ::                       &
        beta,entr_rate


!
!  local ensemble dependent variables in this routine
!
     real,    dimension (its:ite,1:maxens)  ::                         &
        xaa0_ens
     real,    dimension (1:maxens)  ::                                 &
        mbdt_ens
     real,    dimension (1:maxens2) ::                                 &
        edt_ens
     real,    dimension (its:ite,1:maxens2) ::                         &
        edtc
     real,    dimension (its:ite,kts:kte,1:maxens2) ::                 &
        dellat_ens,dellaqc_ens,dellaq_ens,pwo_ens,subt_ens,subq_ens
!
!
!
!***************** the following are your basic environmental
!                  variables. They carry a "_cup" if they are
!                  on model cloud levels (staggered). They carry
!                  an "o"-ending (z becomes zo), if they are the forced
!                  variables. They are preceded by x (z becomes xz)
!                  to indicate modification by some typ of cloud
!
  ! z           = heights of model levels
  ! q           = environmental mixing ratio
  ! qes         = environmental saturation mixing ratio
  ! t           = environmental temp
  ! p           = environmental pressure
  ! he          = environmental moist static energy
  ! hes         = environmental saturation moist static energy
  ! z_cup       = heights of model cloud levels
  ! q_cup       = environmental q on model cloud levels
  ! qes_cup     = saturation q on model cloud levels
  ! t_cup       = temperature (Kelvin) on model cloud levels
  ! p_cup       = environmental pressure
  ! he_cup = moist static energy on model cloud levels
  ! hes_cup = saturation moist static energy on model cloud levels
  ! gamma_cup = gamma on model cloud levels
!
!
  ! hcd = moist static energy in downdraft
  ! zd normalized downdraft mass flux
  ! dby = buoancy term
  ! entr = entrainment rate
  ! zd   = downdraft normalized mass flux
  ! entr= entrainment rate
  ! hcd = h in model cloud
  ! bu = buoancy term
  ! zd = normalized downdraft mass flux
  ! gamma_cup = gamma on model cloud levels
  ! qcd = cloud q (including liquid water) after entrainment
  ! qrch = saturation q in cloud
  ! pwd = evaporate at that level
  ! pwev = total normalized integrated evaoprate (I2)
  ! entr= entrainment rate
  ! z1 = terrain elevation
  ! entr = downdraft entrainment rate
  ! jmin = downdraft originating level
  ! kdet = level above ground where downdraft start detraining
  ! psur        = surface pressure
  ! z1          = terrain elevation
  ! pr_ens = precipitation ensemble
  ! xf_ens = mass flux ensembles
  ! massfln = downdraft mass flux ensembles used in next timestep
  ! omeg = omega from large scale model
  ! mconv = moisture convergence from large scale model
  ! zd      = downdraft normalized mass flux
  ! zu      = updraft normalized mass flux
  ! dir     = "storm motion"
  ! mbdt    = arbitrary numerical parameter
  ! dtime   = dt over which forcing is applied
  ! iact_gr_old = flag to tell where convection was active
  ! kbcon       = LFC of parcel from k22
  ! k22         = updraft originating level
  ! icoic       = flag if only want one closure (usually set to zero!)
  ! dby = buoancy term
  ! ktop = cloud top (output)
  ! xmb    = total base mass flux
  ! hc = cloud moist static energy
  ! hkb = moist static energy at originating level

     real,    dimension (its:ite,kts:kte) ::                           &
        entr_rate_2d,mentrd_rate_2d,he,hes,qes,z,                      &
        heo,heso,qeso,zo,                                              &
        xhe,xhes,xqes,xz,xt,xq,                                        &

        qes_cup,q_cup,he_cup,hes_cup,z_cup,p_cup,gamma_cup,t_cup,      &
        qeso_cup,qo_cup,heo_cup,heso_cup,zo_cup,po_cup,gammao_cup,     &
        tn_cup,                                                        &
        xqes_cup,xq_cup,xhe_cup,xhes_cup,xz_cup,xp_cup,xgamma_cup,     &
        xt_cup,hcot,                                                   &

        xlamue,dby,qc,qrcd,pwd,pw,hcd,qcd,dbyd,hc,qrc,zu,zd,clw_all,   &
        dbyo,qco,qrcdo,pwdo,pwo,hcdo,qcdo,dbydo,hco,qrco,zuo,zdo,      &
        xdby,xqc,xqrcd,xpwd,xpw,xhcd,xqcd,xhc,xqrc,xzu,xzd,            &

  ! cd  = detrainment function for updraft
  ! cdd = detrainment function for downdraft
  ! dellat = change of temperature per unit mass flux of cloud ensemble
  ! dellaq = change of q per unit mass flux of cloud ensemble
  ! dellaqc = change of qc per unit mass flux of cloud ensemble

        cd,cdd,DELLAH,DELLAQ,DELLAT,DELLAQC,dsubt,dsubh,dsubq,          &
        u_cup,v_cup,uc,vc,ucd,vcd,dellu,dellv

  ! aa0 cloud work function for downdraft
  ! edt = epsilon
  ! aa0     = cloud work function without forcing effects
  ! aa1     = cloud work function with forcing effects
  ! xaa0    = cloud work function with cloud effects (ensemble dependent)
  ! edt     = epsilon

     real,    dimension (its:ite) ::                                   &
       edt,edto,edtx,AA1,AA0,XAA0,HKB,                          &
       HKBO,XHKB,QKB,QKBO,                                    &
       XMB,XPWAV,XPWEV,PWAV,PWEV,PWAVO,                                &
       PWEVO,BU,BUD,BUO,cap_max,xland1,                                    &
       cap_max_increment,closure_n,psum,psumh,sig,sigd,zuhe
     real,    dimension (its:ite,1:ens4) ::                                   &
        axx
     integer,    dimension (its:ite) ::                                &
       kzdown,KDET,K22,KB,JMIN,kstabi,kstabm,K22x,        &   !-lxz
       KBCONx,KBx,KTOPx,ierr,ierr2,ierr3,KBMAX

     integer                              ::                           &
       iloop,nall,iedt,nens,nens3,ki,I,K,KK,iresult,nvar
     real                                 ::                           &
      day,dz,dzo,mbdt,radius,entrd_rate,mentrd_rate,  &
      zcutdown,edtmax,edtmin,depth_min,zkbmax,z_detr,zktop,      &
      massfld,dh,cap_maxs,trash,frh,xlamdd,radiusd,frhd
      real detdo1,detdo2,entdo,dp,subin,detdo,entup,                &
      detup,subdown,entdoj,entupk,detupk,totmas
      real :: zusafe,xxx,xx1,xx2,zutop,power_entr,zustart,zufinal,dzm1,dzp1


     integer :: jprnt,k1,k2,kbegzu,kdefi,kfinalzu,kstart,jmini,levadj
     logical :: keep_going
     real tot_time_hr,xff_shal(9),blqe,xkshal
     character*50 :: ierrc(its:ite)
     real,    dimension (its:ite,kts:kte) ::                           &
       up_massentr,up_massdetr,dd_massentr,dd_massdetr                 &
      ,up_massentro,up_massdetro,dd_massentro,dd_massdetro
     real dts,fp,fpi,lambau,pgcon,up_massent,up_massdet
     real,    dimension (kts:kte) :: smth
     real :: xff_mid(its:ite,2)

!-srf- begin
!-srf- for training
    real, dimension( its:ite , jts:jte )   , intent (in)               &
          :: weight_GR,weight_W,weight_MC,weight_ST,weight_AS
!-srf- for convective transport-start
     integer, intent (in   )  ::    mgmxp,  mgmzp ,mynum
     integer, intent (inout)  ::     &
               ierr4d   (mgmxp)      &
              ,jmin4d          (mgmxp)             &
              ,kdet4d          (mgmxp)             &
              ,k224d        (mgmxp)      &
              ,kbcon4d         (mgmxp)             &
              ,ktop4d          (mgmxp)             &
              ,kpbl4d          (mgmxp)             &
              ,kstabi4d        (mgmxp)             &
              ,kstabm4d        (mgmxp)

     real   , intent (inout)  ::     &
               xmb4d        (mgmxp)             &
              ,edt4d        (mgmxp)             &
              ,pwav4d        (mgmxp)

     real   , intent (inout)  ::     &
               pcup5d               (mgmzp,mgmxp)      &
              ,up_massentr5d (mgmzp,mgmxp)      &
              ,up_massdetr5d (mgmzp,mgmxp)      &
              ,dd_massentr5d (mgmzp,mgmxp)      &
              ,dd_massdetr5d (mgmzp,mgmxp)      &
              ,zup5d             (mgmzp,mgmxp)      &
              ,zdn5d                (mgmzp,mgmxp)      &
              ,prup5d               (mgmzp,mgmxp)      &
              ,prdn5d               (mgmzp,mgmxp)      &
              ,clwup5d              (mgmzp,mgmxp)      &
              ,tup5d                (mgmzp,mgmxp)
!-srf- for convective transport-end
!-srf- for diurnal cycle
!----------------------------------------------------------------------
!
      integer :: iversion=1
      real :: umean,T_star
      logical, intent(IN) :: DICYCLE,ENTRNEW,MASSCON,COUPL_MPHYSICS

      real, dimension (its:ite)         :: aa1_bl,hkbo_bl,tau_bl,tau_ecmwf,wmean
      real, dimension (its:ite,kts:kte) :: tn_bl, qo_bl, qeso_bl, heo_bl, heso_bl &
                                          ,qeso_cup_bl,qo_cup_bl, heo_cup_bl,heso_cup_bl&
                                          ,gammao_cup_bl,tn_cup_bl,hco_bl,DBYo_bl
      real, dimension(its:ite) :: xf_dicycle
      real :: C_up, E_dn,G_rain,trash2

!srf- end
!
!proportionality constant to estimate pressure gradient of updraft (Zhang and Wu, 2003, JAS
!
!- ecmwf formulation
     lambau=2.
     pgcon=0.
!    if(imid.eq.1)then
!- SAS formulation
!     lambau=0.
!     pgcon=-.55
!    endif

      zustart=zustart_dp
      zufinal=zufinal_dp
      zutop = zutop_dp
      if(imid.eq.1)zutop=.1

      levadj=5
!      power_entr=2. ! 1.2
!      day=86400.


!     cap_maxs=225.
!     if(imid.eq.1)cap_maxs=150.
      cap_maxs=150.
      if(imid.eq.1)cap_maxs=150.
      do i=its,itf
        edto(i) = 0.0
        closure_n(i)=16.
        xland1(i)=xland(i) ! 1.
        xmb_out(i)=0.
        cap_max(i)=cap_maxs
        cap_max_increment(i)=20.
        if(imid.eq.1)cap_max_increment(i)=10.
        if(xland(i).gt.1.5 .or. xland(i).lt.0.5)then
            xland1(i)=0.
!           ztexec(i)=0.
!           zqexec(i)=0.
            if(imid.eq.0)cap_max(i)=cap_maxs-50.
            if(imid.eq.1)cap_max(i)=cap_maxs-50.
        endif
        ierrc(i)=" "
!       cap_max_increment(i)=1.
      enddo
!
!--- initial entrainment rate (these may be changed later on in the


      if(imid.eq.1)entr_rate=3.e-4

      if(ENTRNEW) then
         entr_rate  = 1.e-3
         mentrd_rate= 0.3*entr_rate
         !-this is valid only for deltax=5000 meters
	     !sig (:)=0.2
	     !sig (:)=0.1 !5km
	      sig (:)=0.4 !7.5km
         !-this is valid only for deltax=20000 meters
	      !sig (:)=0.9
	      !sig (:)=0.7 !03032016
	      sigd(:)=1.
      else
        !entr_rate  = 1.e-4 !-21/10/2015 7.0e-5
         entr_rate  = 7.0e-5
         mentrd_rate= 2.0*entr_rate

         radius =.2/entr_rate
         radiusd=.2/mentrd_rate
         frh =3.14*radius*radius/dx/dx
         frhd=3.14*radiusd*radiusd/dx/dx
         !orig      if(frh .gt. 0.55)then
         !orig         frh=.55
         if(frh .gt. 0.7)then
           frh=.7
           radius=sqrt(frh*dx*dx/3.14)
           entr_rate=.2/radius
         endif

         if(frhd .gt. 0.7)then
           frhd=.7
           radiusd=sqrt(frhd*dx*dx/3.14)
           mentrd_rate=.2/radiusd
         endif
         do i=its,itf
           sig(i)=(1.-frh)**2
           sigd(i)=1.
           !        sigd(i)=sig(i)/(1.-frhd)**2
           !        if(imid.eq.0)sig(i)=sig(i)
         enddo
      endif

!--- entrainment of mass
!
      xlamdd=mentrd_rate
!
!--- initial detrainmentrates
!
      do k=kts,ktf
       do i=its,itf
        hcot(i,k)=0.
        z(i,k)=zo(i,k)
        xz(i,k)=zo(i,k)
        cupclw(i,k)=0.
        cd(i,k)=1.*entr_rate
        cdd(i,k)=xlamdd
        hcdo(i,k)=0.
        qrcdo(i,k)=0.
        dellaqc(i,k)=0.
       enddo
      enddo
!
!--- max/min allowed value for epsilon (ratio downdraft base mass flux/updraft
!    base mass flux
!
!      edtmax=1. !srf 29112015
      edtmax=0.3
      if(imid.eq.1)edtmax=.3
      edtmin=.1
!
!--- minimum depth (m), clouds must have
!
      depth_min=1000.
      if(imid.eq.1)depth_min=500.
!
!--- maximum depth (mb) of capping
!--- inversion (larger cap = no convection)
!
      DO i=its,itf
        kbmax(i)=1
        aa0(i)=0.
        aa1(i)=0.
        edt(i)=0.
        kstabm(i)=ktf-1
        IERR2(i)=0
        IERR3(i)=0
      enddo
!     do i=its,itf
!         cap_max(i)=cap_maxs
!         cap_max3(i)=25.
!         if(gsw(i,j).lt.1.)cap_max(i)=25.
        iresult=0

!     enddo
!      write(11,*)'ipr,ips = ',ipr,its,cap_max(its)
!
!--- max height(m) above ground where updraft air can originate
!
      zkbmax=4000.
      if(imid.eq.1)zkbmax=3000.
!
!--- height(m) above which no downdrafts are allowed to originate
!
      zcutdown=3000.
!
!--- depth(m) over which downdraft detrains all its mass
!
      z_detr=1000.
!     if(imid.eq.1)z_detr=800.
!
      do nens=1,maxens
         mbdt_ens(nens)=(float(nens)-3.)*dtime*1.e-3+dtime*5.E-03
      enddo
      mbdt_ens(1)=10. !dtime*1.E-02
      do nens=1,maxens2
         edt_ens(nens)=.95-float(nens)*.01
      enddo
!
!--- environmental conditions, FIRST HEIGHTS
!
      do i=its,itf
         if(ierr(i).ne.20)then
            do k=1,maxens*maxens2*maxens3
               xf_ens(i,j,(iens-1)*maxens*maxens2*maxens3+k)=0.
               pr_ens(i,j,(iens-1)*maxens*maxens2*maxens3+k)=0.
            enddo
         endif
      enddo
!
!--- calculate moist static energy, heights, qes
!
      call cup_env(z,qes,he,hes,t,q,p,z1, &
           psur,ierr,tcrit,-1,xl,cp,   &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
      call cup_env(zo,qeso,heo,heso,tn,qo,po,z1, &
           psur,ierr,tcrit,-1,xl,cp,   &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)

!
!--- environmental values on cloud levels
!
      call cup_env_clev(t,qes,q,he,hes,z,p,qes_cup,q_cup,he_cup, &
           hes_cup,z_cup,p_cup,gamma_cup,t_cup,psur, &
           ierr,z1,xl,rv,cp,          &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
      call cup_env_clev(tn,qeso,qo,heo,heso,zo,po,qeso_cup,qo_cup, &
           heo_cup,heso_cup,zo_cup,po_cup,gammao_cup,tn_cup,psur,  &
           ierr,z1,xl,rv,cp,          &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
      do i=its,itf
        if(ierr(i).eq.0)then
          u_cup(i,kts)=us(i,kts)
          v_cup(i,kts)=vs(i,kts)
          do k=kts+1,ktf
           u_cup(i,k)=.5*(us(i,k-1)+us(i,k))
           v_cup(i,k)=.5*(vs(i,k-1)+vs(i,k))
          enddo
        endif
      enddo
      do i=its,itf
        if(ierr(i).eq.0)then
        if(aaeq(i).lt.-0.1)then
           ierr(i)=20
        endif
!     if(ierr(i).eq.0)then
!
      do k=kts,ktf
        if(zo_cup(i,k).gt.zkbmax+z1(i))then
          kbmax(i)=k
          go to 25
        endif
      enddo
 25   continue
!
!--- level where detrainment for downdraft starts
!
      do k=kts,ktf
        if(zo_cup(i,k).gt.z_detr+z1(i))then
          kdet(i)=k
          go to 26
        endif
      enddo
 26   continue
!
      endif
      enddo
      if(j.eq.jpr)then
      i=ipr
        write(0,*)i,j,'k,z(i,k),he(i,k),hes(i,k)'
      do k=kts,ktf
        write(0,*)k,z(i,k),he(i,k),hes(i,k)
      enddo
      do k=kts,ktf
        write(0,*)k,zo(i,k),heo(i,k),heso(i,k)
      enddo
      endif

!
!
!
!------- DETERMINE LEVEL WITH HIGHEST MOIST STATIC ENERGY CONTENT - K22
!
      !-srf to increase are coverage change "2" below to "1", the start point
      CALL cup_MAXIMI(HEO_CUP,2,KBMAX,K22,ierr, &
!srf21112015     CALL cup_MAXIMI(HEO_CUP,1,KBMAX,K22,ierr, &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
       DO 36 i=its,itf
         IF(ierr(I).eq.0)THEN
         IF(K22(I).GE.KBMAX(i))then
           ierr(i)=2
           ierrc(i)="could not find k22"
           ktop(i)=0
           k22(i)=0
           kbcon(i)=0
         endif
         endif
 36   CONTINUE
!
!--- DETERMINE THE LEVEL OF CONVECTIVE CLOUD BASE  - KBCON
!

      do i=its,itf
       IF(ierr(I).eq.0)THEN
         if(use_excess == 2) then
             k1=max(1,k22(i)-1)
             k2=k22(i)+1
             hkb(i)=he_cup(i,k22(i))
             hkbo(i)=sum(heo_cup(i,k1:k2))/float(k2-k1+1)!+(xl*zqexec(i)+cp*ztexec(i))/float(k2-k1+1)
         else if(use_excess <= 1)then
             hkb (i)=he_cup (i,k22(i))  +float(use_excess)*(xl*zqexec(i)+cp*ztexec(i))
             hkbo(i)=heo_cup(i,k22(i))  +float(use_excess)*(xl*zqexec(i)+cp*ztexec(i))
         endif  ! excess

       endif ! ierr
      enddo
      if(j.eq.jpr)write(0,*)'k22,kbmax,cap_inc,cap =',k22(ipr),kbmax(ipr),cap_max_increment,cap_max(ipr)
      jprnt=0
      if(j.eq.jpr)jprnt=1
      iloop=1
!     if(imid.eq.1)iloop=5
      call cup_kbcon(ierrc,cap_max_increment,iloop,k22,kbcon,heo_cup,heso_cup, &
           hkbo,ierr,kbmax,po_cup,cap_max, &
           xl,cp,ztexec,zqexec,use_excess,       &
           jprnt,itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte, &
           z_cup,entr_rate,heo)
!
!--- increase detrainment in stable layers
!
      CALL cup_minimi(HEso_cup,Kbcon,kstabm,kstabi,ierr,  &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
      DO i=its,itf
         IF(ierr(I).eq.0)THEN
            kdefi=0
            do k=kts,ktf

               frh = min(qo_cup(i,k)/qeso_cup(i,k),1.)
                   !
               !-------------------------------------------
               if(ENTRNEW) then
               !- v 2
                  if(k >= kbcon(i)) then
                     entr_rate_2d(i,k)=entr_rate*(1.3-frh)*(qeso_cup(i,k)/qeso_cup(i,kbcon(i)))**3
                  else
                     entr_rate_2d(i,k)=entr_rate*(1.3-frh)
                  endif
                  cd(i,k)=0.75e-4*(1.6-frh)

               ELSE
                  !- v 1
                  entr_rate_2d(i,k)=entr_rate*(1.3-frh)
                  cd(i,k)=1.*entr_rate

               ENDIF
               !
               !if(i.eq.ipr.and.j.eq.jpr)write(0,*)'entr_rate',entr_rate,qo_cup(i,k),qeso_cup(i,k)
            enddo
          endif
       enddo
!
! get entrainment and detrainmentrates for updraft
!
!srf
     !call rates_up_curvefitting('UP',ktop,ierr,po_cup,entr_rate_2d,hkbo,heo,heso_cup,zo_cup, &
     !               imid,kstabi,k22,kbcon,its,ite,itf,kts,kte,ktf,zuo,zustart_dp,zufinal_dp,zutop_dp)
     call rates_up_pdf('deep',ktop,ierr,po_cup,entr_rate_2d,hkbo,heo,heso_cup,zo_cup, &
                    kstabi,k22,kbcon,its,ite,itf,kts,kte,ktf,zuo,kpbl)
!
! calculate mass entrainment and detrainment
!
      do k=kts,ktf
      do i=its,itf
         uc(i,k)=0.
         vc(i,k)=0.
         hc(i,k)=0.
         dby (i,k)=0.
         hco (i,k)=0.
         dbyo(i,k)=0.
         up_massentro(i,k)=0.
         up_massdetro(i,k)=0.
         up_massentr (i,k)=0.
         up_massdetr (i,k)=0.
      enddo
      enddo
      do i=its,itf
       IF(ierr(I).eq.0)THEN
         do k=1,kbcon(i)-1
            hc(i,k)=hkb(i)
            hco(i,k)=hkbo(i)
            uc(i,k)=u_cup(i,k22(i))
            vc(i,k)=v_cup(i,k22(i))
         enddo
         k=kbcon(i)
         hc(i,k)=hkb(i)
         uc(i,k)=u_cup(i,k22(i))
         vc(i,k)=v_cup(i,k22(i))
!            hkbo(i)=sum(heo_cup(i,k1:k2))/float(k2-k1+1)+(xl*zqexec(i)+cp*ztexec(i))/float(k2-k1+1)
         DBY(I,Kbcon(i))=Hkb(I)-HES_cup(I,K)
         hco(i,k)=hkbo(i)
         DBYo(I,Kbcon(i))=Hkbo(I)-HESo_cup(I,K)
       endif ! ierr
      enddo
!
!
      do i=its,itf
         if(ierr(i).eq.0)then

         do k=1,ktop(i)+1
          xzu(i,k)= zuo(i,k)
          zu (i,k)= zuo(i,k)
         enddo

         !- mass entrainment and detrinament is defined on model levels

         do k=2,maxloc(zuo(i,:),1)
         !=> below location of maximum value zu -> change entrainment
           dz=zo_cup(i,k)-zo_cup(i,k-1)

           up_massdetro(i,k-1)=cd(i,k-1)*dz*zuo(i,k-1)
           up_massentro(i,k-1)=zuo(i,k)-zuo(i,k-1)+up_massdetro(i,k-1)
           if(up_massentro(i,k-1).lt.0.)then
              up_massentro(i,k-1)=0.
              up_massdetro(i,k-1)=zuo(i,k-1)-zuo(i,k)
              if(zuo(i,k-1).gt.0.) &
                cd(i,k-1)=up_massdetro(i,k-1)/(dz*zuo(i,k-1))
           endif
           if(zuo(i,k-1).gt.0.) &
             entr_rate_2d(i,k-1)=(up_massentro(i,k-1))/(dz*zuo(i,k-1))
         enddo
         do k=maxloc(zuo(i,:),1)+1,ktop(i)
         !=> above location of maximum value zu -> change detrainment
           dz=zo_cup(i,k)-zo_cup(i,k-1)
           up_massentro(i,k-1)=entr_rate_2d(i,k-1)*dz*zuo(i,k-1)
           up_massdetro(i,k-1)=zuo(i,k-1)+up_massentro(i,k-1)-zuo(i,k)
           if(up_massdetro(i,k-1).lt.0.)then
              up_massdetro(i,k-1)=0.
              up_massentro(i,k-1)=zuo(i,k)-zuo(i,k-1)
              if(zuo(i,k-1).gt.0.) &
                entr_rate_2d(i,k-1)=(up_massentro(i,k-1))/(dz*zuo(i,k-1))
           endif

           if(zuo(i,k-1).gt.0.)cd(i,k-1)=up_massdetro(i,k-1)/(dz*zuo(i,k-1))
         enddo

         do k=2,ktf-1
          up_massentr(i,k-1)=up_massentro(i,k-1)
          up_massdetr(i,k-1)=up_massdetro(i,k-1)
         enddo

!---mass con
!        do k=kbcon(i)+1,ktop(i)  ! original
         do k=k22(i)  +1,ktop(i)  ! mass cons option
!---mass con

          hc(i,k)=(hc(i,k-1)*zu(i,k-1)-.5*up_massdetr(i,k-1)*hc(i,k-1)+ &
                         up_massentr(i,k-1)*he(i,k-1))   /            &
                         (zu(i,k-1)-.5*up_massdetr(i,k-1)+up_massentr(i,k-1))
          uc(i,k)=(uc(i,k-1)*zu(i,k-1)-.5*up_massdetr(i,k-1)*uc(i,k-1) &
                         -lambau*up_massdetr(i,k-1)*uc(i,k-1) +         &
                         (up_massentr(i,k-1)+lambau*up_massdetr(i,k-1))*us(i,k-1) &
                         -pgcon*.5*(zu(i,k)+zu(i,k-1))*(u_cup(i,k)-u_cup(i,k-1)))  /            &
                         (zu(i,k-1)-.5*up_massdetr(i,k-1)+up_massentr(i,k-1))
          vc(i,k)=(vc(i,k-1)*zu(i,k-1)-.5*up_massdetr(i,k-1)*vc(i,k-1)                &
                         -lambau*up_massdetr(i,k-1)*vc(i,k-1) +                        &
                         (up_massentr(i,k-1)+lambau*up_massdetr(i,k-1))*vs(i,k-1)      &
                         -pgcon*.5*(zu(i,k)+zu(i,k-1))*(u_cup(i,k)-u_cup(i,k-1))) /    &
                         (zu(i,k-1)-.5*up_massdetr(i,k-1)+up_massentr(i,k-1))
          dby(i,k)=hc(i,k)-hes_cup(i,k)
          hco(i,k)=(hco(i,k-1)*zuo(i,k-1)-.5*up_massdetro(i,k-1)*hco(i,k-1)+ &
                         up_massentro(i,k-1)*heo(i,k-1))   /            &
                         (zuo(i,k-1)-.5*up_massdetro(i,k-1)+up_massentro(i,k-1))
          dbyo(i,k)=hco(i,k)-heso_cup(i,k)
         enddo

         if(ktop(i).lt.kbcon(i)+2)then
            ierr(i)=5
            ierrc(i)='ktop too small'
         endif
         !
         do k=ktop(i)+1,ktf
           HC(i,K)=hes_cup(i,k)
           UC(i,K)=u_cup(i,k)
           VC(i,K)=v_cup(i,k)
           HCo(i,K)=heso_cup(i,k)
           DBY(I,K)=0.
           DBYo(I,K)=0.
           zu(i,k)=0.
           zuo(i,k)=0.
           cd(i,k)=0.
           entr_rate_2d(i,k)=0.
           up_massentr(i,k)=0.
           up_massdetr(i,k)=0.
           up_massentro(i,k)=0.
           up_massdetro(i,k)=0.
         enddo
         if(i.eq.ipr.and.j.eq.jpr)then
            write(0,*)'hcnew = '
            do k=1,ktf
              write(0,*)k,hco(i,k),dbyo(i,k)
            enddo
         endif
      endif
      enddo
!
      DO 37 i=its,itf
         kzdown(i)=0
         if(ierr(i).eq.0)then
            zktop=(zo_cup(i,ktop(i))-z1(i))*.6
            zktop=min(zktop+z1(i),zcutdown+z1(i))
            do k=kts,ktf
              if(zo_cup(i,k).gt.zktop)then
                 kzdown(i)=k
                 go to 37
              endif
              enddo
         endif
 37   CONTINUE
!
!--- DOWNDRAFT ORIGINATING LEVEL - JMIN
!
!     call cup_minimi(HEso_cup,K22,kstabi,JMIN,ierr, &
      call cup_minimi(HEso_cup,K22,kzdown,JMIN,ierr, &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
      DO 100 i=its,ite
         if(i.eq.ipr.and.j.eq.jpr)then
              write(0,*)i,j,'k,p_cup(i,k),heo_cup(i,k),heso_cup(i,k),hco(i,k)'
            do k=kts,ktf
              write(0,*)k,p_cup(i,k),heo_cup(i,k),heso_cup(i,k),hco(i,k)
            enddo
         endif
         IF(ierr(I).eq.0)THEN
!
!--- check whether it would have buoyancy, if there where
!--- no entrainment/detrainment
!
         jmini = jmin(i)
         keep_going = .TRUE.
         do while ( keep_going )
           keep_going = .FALSE.
           if ( jmini - 1 .lt. kdet(i)   ) kdet(i) = jmini-1
           if ( jmini     .ge. ktop(i)-1 ) jmini = ktop(i) - 2
           ki = jmini
           hcdo(i,ki)=heso_cup(i,ki)
           DZ=Zo_cup(i,Ki+1)-Zo_cup(i,Ki)
           dh=0.
           do k=ki-1,1,-1
             hcdo(i,k)=heso_cup(i,jmini)
             DZ=Zo_cup(i,K+1)-Zo_cup(i,K)
             dh=dh+dz*(HCDo(i,K)-heso_cup(i,k))
             if(dh.gt.0.)then
               jmini=jmini-1
               if ( jmini .gt. 5 ) then
                 keep_going = .TRUE.
               else
                 ierr(i) = 9
                 ierrc(i) = "could not find jmini9"
                 exit
               endif
             endif
           enddo
         enddo
         jmin(i) = jmini
         if ( jmini .le. 5 ) then
           ierr(i)=4
           ierrc(i) = "could not find jmini4"
         endif
       ENDIF
100   continue
!
! - Must have at least depth_min m between cloud convective base
!     and cloud top.
!
      do i=its,itf
         IF(ierr(I).eq.0)THEN
            if ( jmin(i) - 1 .lt. kdet(i)   ) kdet(i) = jmin(i)-1
            IF(-zo_cup(I,KBCON(I))+zo_cup(I,KTOP(I)).LT.depth_min)then
               ierr(i)=6
               ierrc(i)="cloud depth very shallow"
            endif
         endif
      enddo

!
!--- normalized downdraft mass flux profile,also work on bottom detrainment
!--- in this routine
!
      do k=kts,ktf
      do i=its,itf
       zd(i,k)=0.
       zdo(i,k)=0.
       cdd(i,k)=0.
       dd_massentr(i,k)=0.
       dd_massdetr(i,k)=0.
       dd_massentro(i,k)=0.
       dd_massdetro(i,k)=0.
       hcdo(i,k)=heso_cup(i,k)
       ucd(i,k)=u_cup(i,k)
       vcd(i,k)=v_cup(i,k)
       dbydo(i,k)=0.
      enddo
      enddo
      do i=its,itf
          beta=.05
          if(imid.eq.1)beta=.0
!         if(xland1(i).gt.0.1)beta=.1
          bud(i)=0.

          IF(ierr(I).eq.0)then

        !- this calls routine to get downdrafts normalized mass flux
        !
        zutop  = 0.2 ! zd at jmin
        zustart= beta  ! zd at surface
        zufinal= 1.  ! zd at jmin-lev_adjust
        mentrd_rate_2d(i,:)=mentrd_rate
        cdd(i,1:jmin(i))=xlamdd
        cdd(i,jmin(i))=0.
        dd_massdetro(i,:)=0.
        dd_massentro(i,:)=0.

        call get_zu_zd_pdf("DOWN",ierr(i),kdet(i),jmin(i),zdo(i,:),kts,kte,ktf,kpbl(i))
!
!        if(jmin(i)-levadj.gt.kdet(i))then
!           call get_zu_zd("DOWN",ierr(i),jmin(i)-levadj,jmin(i),zustart,zufinal,zutop,zdo(i,:),kts,kte,ktf)
!        else
!           call get_zu_zd("DOWNM",ierr(i),jmin(i)-levadj,jmin(i),zustart,zufinal,zutop,zdo(i,:),kts,kte,ktf)
!        endif

        xzd(i,jmin(i))= zdo(i,jmin(i))
        zd (i,jmin(i))= zdo(i,jmin(i))
!        write(92,*)'k,zdo(i,k)'
!       do k=jmin(i),1,-1
!        write(92,*)k,zdo(i,k)
!       enddo

!-srf-21jul2015 mass cons
       !do ki=jmin(i)-1,maxloc(zdo(i,:),1),-1
        do ki=jmin(i)  ,maxloc(zdo(i,:),1),-1
!-srf mass cons
          !=> from jmin to maximum value zd -> change entrainment
          dzo=zo_cup(i,ki+1)-zo_cup(i,ki)
          dd_massdetro(i,ki)=cdd(i,ki)*dzo*zdo(i,ki+1)
          dd_massentro(i,ki)=zdo(i,ki)-zdo(i,ki+1)+dd_massdetro(i,ki)
          if(dd_massentro(i,ki).lt.0.)then
             dd_massentro(i,ki)=0.
             dd_massdetro(i,ki)=zdo(i,ki+1)-zdo(i,ki)
             if(zdo(i,ki+1) > 0.0)&
               cdd(i,ki)=dd_massdetro(i,ki)/(dzo*zdo(i,ki+1))
          endif
          if(zdo(i,ki+1) > 0.0)&
            mentrd_rate_2d(i,ki)=dd_massentro(i,ki)/(dzo*zdo(i,ki+1))
        enddo
        mentrd_rate_2d(i,1)=0.
        do ki=maxloc(zdo(i,:),1)-1,1,-1
          !=> from maximum value zd to surface -> change detrainment
          dzo=zo_cup(i,ki+1)-zo_cup(i,ki)
          dd_massentro(i,ki)=mentrd_rate_2d(i,ki)*dzo*zdo(i,ki+1)
          dd_massdetro(i,ki) = zdo(i,ki+1)+dd_massentro(i,ki)-zdo(i,ki)
          if(dd_massdetro(i,ki).lt.0.)then
            dd_massdetro(i,ki)=0.
            dd_massentro(i,ki)=zdo(i,ki)-zdo(i,ki+1)
            if(zdo(i,ki+1) > 0.0)&
              mentrd_rate_2d(i,ki)=dd_massentro(i,ki)/(dzo*zdo(i,ki+1))
          endif
          if(zdo(i,ki+1) > 0.0)&
            cdd(i,ki)= dd_massdetro(i,ki)/(dzo*zdo(i,ki+1))
        enddo

!         write(92,*)'k,zdo(i,k),dd_massentro(i,k),dd_massdetro(i,k)'
        do k=jmin(i),1,-1
          xzd(i,k)= zdo(i,k)
          zd (i,k)= zdo(i,k)
          dd_massentr(i,k)=dd_massentro(i,k)
          dd_massdetr(i,k)=dd_massdetro(i,k)
!         write(92,*)k,zdo(i,k),dd_massentro(i,k),dd_massdetro(i,k)
        enddo

! downdraft moist static energy + moisture budget
            dbydo(i,jmin(i))=hcdo(i,jmin(i))-heso_cup(i,jmin(i))
            bud(i)=dbydo(i,jmin(i))*(zo_cup(i,jmin(i)+1)-zo_cup(i,jmin(i)))
!-srf-21jul2015 mass cons
           !do ki=jmin(i)-1,1,-1
            do ki=jmin(i)  ,1,-1
!-srf mass cons
             dzo=zo_cup(i,ki+1)-zo_cup(i,ki)
             ucd(i,ki)=(ucd(i,ki+1)*zdo(i,ki+1)                       &
                         -.5*dd_massdetro(i,ki)*ucd(i,ki+1)+ &
                        dd_massentro(i,ki)*us(i,ki)          &
                        -pgcon*zdo(i,ki+1)*(us(i,ki+1)-us(i,ki)))   /     &
                        (zdo(i,ki+1)-.5*dd_massdetro(i,ki)+dd_massentro(i,ki))
             vcd(i,ki)=(vcd(i,ki+1)*zdo(i,ki+1)                       &
                         -.5*dd_massdetro(i,ki)*vcd(i,ki+1)+ &
                        dd_massentro(i,ki)*vs(i,ki)            &
                        -pgcon*zdo(i,ki+1)*(vs(i,ki+1)-vs(i,ki)))   /  &
                        (zdo(i,ki+1)-.5*dd_massdetro(i,ki)+dd_massentro(i,ki))
             hcdo(i,ki)=(hcdo(i,ki+1)*zdo(i,ki+1)                       &
                         -.5*dd_massdetro(i,ki)*hcdo(i,ki+1)+ &
                        dd_massentro(i,ki)*heo(i,ki))   /            &
                        (zdo(i,ki+1)-.5*dd_massdetro(i,ki)+dd_massentro(i,ki))
             dbydo(i,ki)=hcdo(i,ki)-heso_cup(i,ki)
             if(i.eq.ipr.and.j.eq.jpr)write(0,*)'ki,bud = ',ki,bud(i),hcdo(i,ki)
             bud(i)=bud(i)+dbydo(i,ki)*dzo
            enddo
          endif

        if(bud(i).gt.0)then
          ierr(i)=7
          ierrc(i)='downdraft is not negatively buoyant '
        endif
      enddo
!
!--- calculate moisture properties of downdraft
!
      call cup_dd_moisture_new(ierrc,zdo,hcdo,heso_cup,qcdo,qeso_cup, &
           pwdo,qo_cup,zo_cup,dd_massentro,dd_massdetro,jmin,ierr,gammao_cup, &
           pwevo,bu,qrcdo,qo,heo,tn_cup,1,xl, &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
!
!--- calculate moisture properties of updraft
!
      if(imid.eq.0)then
      call cup_up_moisture('deep',ierr,zo_cup,qco,qrco,pwo,pwavo, &
           ccnclean,p_cup,kbcon,ktop,cd,dbyo,clw_all, &
           t_cup,qo,GAMMAo_cup,zuo,qeso_cup,k22,qo_cup,xl,        &
           ZQEXEC,use_excess,ccn,rho,up_massentr,up_massdetr,psum,psumh,&
           autoconv,aeroevap,1,itf,jtf,ktf,j,ipr,jpr, &
           its,ite, jts,jte, kts,kte)
      else if(imid.eq.1)then
      call cup_up_moisture('shallow',ierr,zo_cup,qco,qrco,pwo,pwavo, &
           ccnclean,p_cup,kbcon,ktop,cd,dbyo,clw_all, &
           t_cup,qo,GAMMAo_cup,zuo,qeso_cup,k22,qo_cup,xl,        &
           ZQEXEC,use_excess,ccn,rho,up_massentr,up_massdetr,psum,psumh,&
           autoconv,aeroevap,1,itf,jtf,ktf,j,ipr,jpr, &
           its,ite, jts,jte, kts,kte)
      endif
      do i=its,itf
      if(ierr(i).eq.0)then
      do k=kts,ktop(i)
         cupclw(i,k)=qrco(i,k)        ! my mod
      enddo
      endif
      enddo
!
!--- calculate workfunctions for updrafts
!
      call cup_up_aa0(aa0,z,zu,dby,GAMMA_CUP,t_cup, &
           kbcon,ktop,ierr,           &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
      call cup_up_aa0(aa1,zo,zuo,dbyo,GAMMAo_CUP,tn_cup, &
           kbcon,ktop,ierr,           &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
      do i=its,itf
         if(ierr(i).eq.0)then
           if(aa1(i).eq.0.)then
               ierr(i)=17
               ierrc(i)="cloud work function zero"
           endif
           aaeq(i)=aa1(i)
         endif
      enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!=================================================================
!-srf- begin
      !--- AA1 from boundary layer (bl) processes only
      aa1_bl    (:) = 0.0
      xf_dicycle(:) = 0.0
      !- way to calculate the fraction of cape consumed by shallow convection
      iversion=1 ! ecmwf
      !iversion=0 ! orig
      !
      ! Bechtold et al 2008 time-scale of cape removal
      DO i=its,itf
            if(ierr(i).eq.0)then
                !- mean vertical velocity
!               wmean(i) = 5.0 ! m/s ! in the future change for Wmean == integral( W dz) / cloud_depth
                wmean(i) = 10.0 ! m/s ! in the future change for Wmean == integral( W dz) / cloud_depth
!               wmean(i) = 7.0 ! m/s ! in the future change for Wmean == integral( W dz) / cloud_depth
                !- time-scale cape removal from  Betchold et al. 2008
                tau_ecmwf(i)=( zo_cup(i,ktop(i))- zo_cup(i,kbcon(i)) ) / wmean(i)
                tau_ecmwf(i)= tau_ecmwf(i) * (1.0061 + 1.23E-2 * (dx/1000.))! dx must be in meters
            endif
      enddo
      !
      IF(dicycle) then

        DO i=its,itf

            if(ierr(i).eq.0)then
                if(xland(i) > 0.9 ) then
                  !- over water
                  umean= 2.0+sqrt(2.0*(US(i,1)**2+VS(i,1)**2+US(i,kbcon(i))**2+VS(i,kbcon(i))**2))
                  tau_bl(i) = (zo_cup(i,kbcon(i))- z1(i)) /umean
                else
                  !- over land
                  tau_bl(i) =( zo_cup(i,ktop(i))- zo_cup(i,kbcon(i)) ) / wmean(i)
                endif

            endif
        ENDDO

        if(iversion == 1) then
        !-- version ecmwf

           !- T_star = temp scale in original paper = 1 K
           ! T_star = 1.0
           T_star = 3.0

           !-- calculate pcape from BL forcing only
            call cup_up_aa1bl(aa1_bl,t,tn,q,qo,dtime, &
                              zo_cup,zuo,dbyo_bl,GAMMAo_CUP_bl,tn_cup_bl, &
                              kbcon,ktop,ierr,               &
                              itf,jtf,ktf,its,ite, jts,jte, kts,kte)

            DO i=its,itf

                  if(ierr(i).eq.0)then

                      !- only for convection rooting in the PBL
                      !if(zo_cup(i,kbcon(i))-z1(i) > 500.0) then !- instead 500 -> zo_cup(kpbl(i))
                      !        aa1_bl(i) = 0.0
                      !else
                      !        !- multiply aa1_bl/T_* by the " time-scale" - tau_bl
                         aa1_bl(i) = (aa1_bl(i)/T_star) * tau_bl(i)
                      !endif
                      !print*,'aa0,aa1bl=',aa0(i),aa1_bl(i),aa0(i)-aa1_bl(i),tau_bl(i)!,dtime,xland(i)
                  endif
            ENDDO

        else

          !- version for real cloud-work function

          !-get the profiles modified only by bl tendencies
          DO i=its,itf
           tn_bl(i,:)=0.;qo_bl(i,:)=0.
           if ( ierr(i) == 0 )then
            !below kbcon -> modify profiles
            tn_bl(i,1:kbcon(i)) = tn(i,1:kbcon(i))
            qo_bl(i,1:kbcon(i)) = qo(i,1:kbcon(i))
            !above kbcon -> keep environment profiles
            tn_bl(i,kbcon(i)+1:kte) = t(i,kbcon(i)+1:kte)
            qo_bl(i,kbcon(i)+1:kte) = q(i,kbcon(i)+1:kte)
           endif
          ENDDO
          !--- calculate moist static energy, heights, qes, ... only by bl tendencies
          call cup_env(zo,qeso_bl,heo_bl,heso_bl,tn_bl,qo_bl,po,z1, &
                     psur,ierr,tcrit,-1,xl,cp,   &
                     itf,jtf,ktf, its,ite, jts,jte, kts,kte)
          !--- environmental values on cloud levels only by bl tendencies
          call cup_env_clev(tn_bl,qeso_bl,qo_bl,heo_bl,heso_bl,zo,po,qeso_cup_bl,qo_cup_bl, &
                              heo_cup_bl,heso_cup_bl,zo_cup,po_cup,gammao_cup_bl,tn_cup_bl,psur,  &
                              ierr,z1,xl,rv,cp,               &
                              itf,jtf,ktf,its,ite, jts,jte, kts,kte)
          DO i=its,itf
            IF(ierr(I).eq.0)THEN
             if(use_excess == 2) then
               k1=max(1,k22(i)-1)
               k2=k22(i)+1
               hkbo_bl(i)=sum(heo_cup_bl(i,k1:k2))/float(k2-k1+1) !+(xl*zqexec(i)+cp*ztexec(i))/float(k2-k1+1)
             else if(use_excess <= 1)then
               hkbo_bl(i)=heo_cup_bl(i,k22(i))  +float(use_excess)*(xl*zqexec(i)+cp*ztexec(i))
             endif  ! excess
            endif ! ierr
          ENDDO
          DO k=kts,ktf
           do i=its,itf
             hco_bl (i,k)=0.
             DBYo_bl(i,k)=0.
           enddo
          ENDDO
          DO i=its,itf
            IF(ierr(I).eq.0)THEN
             do k=1,kbcon(i)-1
              hco_bl(i,k)=hkbo_bl(i)
             enddo
             k=kbcon(i)
             hco_bl (i,k)=hkbo_bl(i)
             DBYo_bl(i,k)=Hkbo_bl(i) - HESo_cup_bl(i,k)
            ENDIF
          ENDDO
!
!
          DO i=its,itf
            if(ierr(i).eq.0)then
               do k=kbcon(i)+1,ktop(i)
                    hco_bl(i,k)=(hco_bl(i,k-1)*zuo(i,k-1)-.5*up_massdetro(i,k-1)*hco_bl(i,k-1)+ &
                               up_massentro(i,k-1)*heo_bl(i,k-1))   /               &
                               (zuo(i,k-1)-.5*up_massdetro(i,k-1)+up_massentro(i,k-1))
                    dbyo_bl(i,k)=hco_bl(i,k)-heso_cup_bl(i,k)
               enddo
               do k=ktop(i)+1,ktf
                  hco_bl (i,k)=heso_cup_bl(i,k)
                  dbyo_bl(i,k)=0.0
               enddo
            endif
          ENDDO

          !--- calculate workfunctions for updrafts
          call cup_up_aa0(aa1_bl,zo,zuo,dbyo_bl,GAMMAo_CUP_bl,tn_cup_bl, &
                        kbcon,ktop,ierr,           &
                        itf,jtf,ktf,its,ite, jts,jte, kts,kte)

          DO i=its,itf

            if(ierr(i).eq.0)then
                !- get the increment on AA0 due the BL processes
                aa1_bl(i) = aa1_bl(i) - aa0(i)
                !- only for convection rooting in the PBL
                !if(zo_cup(i,kbcon(i))-z1(i) > 500.0) then !- instead 500 -> zo_cup(kpbl(i))
                !   aa1_bl(i) = 0.0
                !else
                !   !- multiply aa1_bl the "normalized time-scale" - tau_bl/ model_timestep
                   aa1_bl(i) = aa1_bl(i)* tau_bl(i)/ dtime
                !endif
                print*,'aa0,aa1bl=',aa0(i),aa1_bl(i),aa0(i)-aa1_bl(i),tau_bl(i)!,dtime,xland(i)
            endif
           ENDDO
        ENDIF
     ENDIF  ! version of implementation
        !i=maxloc(aa1_bl,1)
        !if(aa1_bl(i)>0. .and. ierr(i) == 0) print*,'aa0,aa1bl=',j,i,aa0(i),aa1_bl(i), xland(i) &
        !                    ,tau_bl(i) ; call flush(6)!,tau_bl(i),dtime,xland(i)
        !aa1_bl(:)=0.0 !tmp
!-srf -end
!=================================================================


       do i=1,ens4
       axx(:,i)=aa1(:)
       enddo

!
!--- DETERMINE DOWNDRAFT STRENGTH IN TERMS OF WINDSHEAR
!
      call cup_dd_edt(ierr,us,vs,zo,ktop,kbcon,edt,po,pwavo, &
           pwo,ccn,pwevo,edtmax,edtmin,maxens2,edtc,psum,psumh, &
           ccnclean,rho,aeroevap,itf,jtf,ktf,j,ipr,jpr, &
           its,ite, jts,jte, kts,kte)
      do 250 iedt=1,maxens2
        do i=its,itf
         if(ierr(i).eq.0)then
         edt(i)=sigd(i)*edtc(i,iedt)
         edto(i)=sigd(i)*edtc(i,iedt)
         edtx(i)=sigd(i)*edtc(i,iedt)
         if(maxens2.eq.3)then
            edt(i)=sigd(i)*edtc(i,3)
            edto(i)=sigd(i)*edtc(i,3)
            edtx(i)=sigd(i)*edtc(i,3)
         endif
         endif
        enddo
        do k=kts,ktf
        do i=its,itf
           subt_ens(i,k,iedt)=0.
           subq_ens(i,k,iedt)=0.
           dellat_ens(i,k,iedt)=0.
           dellaq_ens(i,k,iedt)=0.
           dellaqc_ens(i,k,iedt)=0.
           pwo_ens(i,k,iedt)=0.
        enddo
        enddo
!
!     do i=its,itf
!       if(j.eq.jpr.and.i.eq.ipr)then
!         do k=1,ktop(i)+1
!           write(0,*)zu(i,k),zd(i,k),pwo(i,k),pwdo(i,k)
!         enddo
!       endif
!     enddo
!
!--- change per unit mass that a model cloud would modify the environment
!
!--- 1. in bottom layer
!
        dellu=0.
        dellv=0.
        dellah=0.
        dsubt=0.
        dsubh=0.
        dellat=0.
        dellaq=0.
        dellaqc=0.
        dsubq=0.
!
!----------------------------------------------  cloud level ktop
!
!- - - - - - - - - - - - - - - - - - - - - - - - model level ktop-1
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!
!----------------------------------------------  cloud level k+2
!
!- - - - - - - - - - - - - - - - - - - - - - - - model level k+1
!
!----------------------------------------------  cloud level k+1
!
!- - - - - - - - - - - - - - - - - - - - - - - - model level k
!
!----------------------------------------------  cloud level k
!
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!
!----------------------------------------------  cloud level 3  _cup
!
!- - - - - - - - - - - - - - - - - - - - - - - - model level 2
!
!----------------------------------------------  cloud level 2  _cup
!
!- - - - - - - - - - - - - - - - - - - - - - - - model level 1

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

if(.not. masscon) then
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

!-------------------------------------------------- form orig
      do i=its,itf
        if(ierr(i).eq.0)then
         dp=100.*(po_cup(i,1)-po_cup(i,2))
         dellu(i,1)=(edto(i)*zdo(i,2)*ucd(i,2)   &
                     -edto(i)*zdo(i,2)*u_cup(i,2))*g/dp
         dellv(i,1)=(edto(i)*zdo(i,2)*vcd(i,2)   &
                     -edto(i)*zdo(i,2)*v_cup(i,2))*g/dp
         dellah(i,1)=(edto(i)*zdo(i,2)*hcdo(i,2)   &
                     -edto(i)*zdo(i,2)*heo_cup(i,2))*g/dp
         dellaq(i,1)=(edto(i)*zdo(i,2)*qrcdo(i,2)   &
                     -edto(i)*zdo(i,2)*qo_cup(i,2))*g/dp
         dsubh(i,1)=0.
         dsubq(i,1)=0.
         if(j.eq.jpr.and.i.eq.ipr)then
             write(0,*)'dellas(1) '
             write(0,*)edto(i),zdo(i,2),dp
             write(0,*)heo_cup(i,2),hcdo(i,2),qo_cup(i,2),qrcdo(i,2)
         endif
         do k=kts+1,ktop(i)
! these three are only used at or near mass detrainment and/or entrainment levels
            entupk=0.
            detupk=0.
            entdoj=0.
! detrainment and entrainment for fowndrafts
            detdo=edto(i)*dd_massdetro(i,k)
            entdo=edto(i)*dd_massentro(i,k)
! entrainment/detrainment for updraft
            entup=up_massentro(i,k)
            detup=up_massdetro(i,k)
! subsidence by downdrafts only
            subin=-zdo(i,k+1)*edto(i)
            subdown=-zdo(i,k)*edto(i)
!
!         SPECIAL LEVELS
!
! updraft originates at k22, only updraft term at k22-1 is (zu-entupk)*he
!           if(k.eq.k22(i)-1)then
!              entupk=zuo(i,k+1)
!           endif
! downdraft originating level, similiar to k22-1 for updraft
!-srf- 21jul2015 mass con
            !if(k.eq.jmin(i))then
            !   entdoj=edto(i)*zdo(i,k)
            !endif
!-srf- mass con
            if(k.eq.ktop(i))then
               detupk=zuo(i,ktop(i))
               subin=0.
               subdown=0.
               detdo=0.
               entdo=0.
               entup=0.
               detup=0.
            endif
            totmas=subin-subdown+detup-entup-entdo+ &
             detdo-entupk-entdoj+detupk+zuo(i,k+1)-zuo(i,k)
!               print *,'*********************',k,totmas
!              write(0,123)k,subin+zuo(i,k+1),subdown-zuo(i,k),detup,entup, &
!                          detdo,entdo,entupk,detupk
!             write(8,*)'totmas = ',k,totmas
            if(abs(totmas).gt.1.e-6)then
               write(0,*)'*********************',i,j,k,totmas
!              print *,jmin(i),k22(i),kbcon(i),ktop(i)
               write(0,123)k,subin,subdown,detup,entup, &
                           detdo,entdo,entupk,detupk
123     formAT(1X,i2,8E12.4)
!        call wrf_error_fatal ( 'totmas .gt.1.e-6' )
            endif
            dp=100.*(po_cup(i,k)-po_cup(i,k+1))
            dellah(i,k)=(detup*.5*(HCo(i,K+1)+HCo(i,K)) &
                    +detdo*.5*(HCDo(i,K+1)+HCDo(i,K)) &
                    -entup*heo(i,k) &
                    -entdo*heo(i,k) &
                    +subin*heo_cup(i,k+1) &
                    -subdown*heo_cup(i,k) &
                    +detupk*(hco(i,ktop(i))-heo_cup(i,ktop(i)))    &
                    -entupk*heo_cup(i,k22(i)) &
                    -entdoj*heo_cup(i,jmin(i)) &
                     )*g/dp
            dellu(i,k)=((detup+lambau*detup)*.5*(uC(i,K+1)+uC(i,K)) &
                    +detdo*.5*(UCD(i,K+1)+UCD(i,K)) &
                    -(entup+lambau*detup)*us(i,k) &
                    -entdo*us(i,k) &
                    +subin*u_cup(i,k+1) &
                    -subdown*u_cup(i,k) &
                    +detupk*(uc(i,ktop(i))-u_cup(i,ktop(i)))    &
                    -entupk*u_cup(i,k22(i)) &
                    -entdoj*u_cup(i,jmin(i)) &
                     )*g/dp
            dellv(i,k)=((detup+lambau*detup)*.5*(vC(i,K+1)+vC(i,K)) &
                    +detdo*.5*(vCD(i,K+1)+vCD(i,K)) &
                    -(entup+lambau*detup)*vs(i,k) &
                    -entdo*vs(i,k) &
                    +subin*v_cup(i,k+1) &
                    -subdown*v_cup(i,k) &
                    +detupk*(vc(i,ktop(i))-v_cup(i,ktop(i)))    &
                    -entupk*v_cup(i,k22(i)) &
                    -entdoj*v_cup(i,jmin(i)) &
                     )*g/dp
            dellaq(i,k)=(detup*.5*(qco(i,K+1)+qco(i,K)-qrco(i,k+1)-qrco(i,k)) &
                    +detdo*.5*(qrcdo(i,K+1)+qrcdo(i,K)) &
                    -entup*qo(i,k) &
                    -entdo*qo(i,k) &
                    +subin*qo_cup(i,k+1) &
                    -subdown*qo_cup(i,k) &
                    +detupk*(qco(i,ktop(i))-qrco(i,ktop(i))-qo_cup(i,ktop(i)))    &
                    -entupk*qo_cup(i,k22(i)) &
                    -entdoj*qo_cup(i,jmin(i)) &
                     )*g/dp
!
! updraft subsidence only
!
           if(k.lt.ktop(i))then
             dsubh(i,k)=(zuo(i,k+1)*heo_cup(i,k+1) &
                    -zuo(i,k)*heo_cup(i,k))*g/dp
             dsubq(i,k)=(zuo(i,k+1)*qo_cup(i,k+1) &
                    -zuo(i,k)*qo_cup(i,k))*g/dp
             dellu(i,k)=dellu(i,k)+(zuo(i,k+1)*u_cup(i,k+1) &
                    -zuo(i,k)*u_cup(i,k) &
                    +pgcon*.5*(zuo(i,k)+zuo(i,k+1))*(u_cup(i,k+1)-u_cup(i,k)))*g/dp
             dellv(i,k)=dellv(i,k)+(zuo(i,k+1)*v_cup(i,k+1) &
                    -zuo(i,k)*v_cup(i,k) &
                    +pgcon*.5*(zuo(i,k)+zuo(i,k+1))*(v_cup(i,k+1)-v_cup(i,k)))*g/dp
!            dellv(i,k)=dellv(i,k)+(zuo(i,k+1)*v_cup(i,k+1) &
!                   -zuo(i,k)*v_cup(i,k)   &
!                   +pgcon*.5*(zuo(i,k)+zuo(i,k+1))*(v_cup(i,k+1)-v_cup(i,k)))*g/dp
           if(i.eq.ipr.and.j.eq.jpr)then
            write(0,*)'dq3',k,zuo(i,k+1)*qo_cup(i,k+1),zuo(i,k)*qo_cup(i,k)
           endif
!          elseif (k.eq.ktop(i))then
!            dsubt(i,k)=(-zuo(i,k)*heo_cup(i,k))*g/dp
!            dsubq(i,k)=(-zuo(i,k)*qo_cup(i,k))*g/dp
           endif
           if(i.eq.ipr .and. j.eq.jpr .and. k.eq.ktop(i))then
            write(0,*)'dq4',k,dellaq(i,k),detup*.5*(qco(i,K+1)+qco(i,K)-qrco(i,k+1)-qrco(i,k)), &
                       -entup*qo(i,k),detupk*(qco(i,ktop(i))-qrco(i,ktop(i))-qo_cup(i,ktop(i)))
            write(0,*)'dq4.2',k,dellaq(i,k),detup*(qco(i,K)-qrco(i,k)), &
                       -entup*qo(i,k),detupk*(qco(i,ktop(i))-qrco(i,ktop(i))-qo_cup(i,ktop(i)))
           endif
!
! in igh res case, subsidence terms are for meighbouring points only. This has to be
! done mass consistent with the della term
       enddo   ! k

        endif
      enddo
!
!-- take out cloud liquid water for detrainment
!
      do k=kts,ktf-1
      do i=its,itf
       dellaqc(i,k)=0.
       if(ierr(i).eq.0)then
         if(k.eq.ktop(i)-0)dellaqc(i,k)= &
                      .01*zuo(i,ktop(i))*qrco(i,ktop(i))* &
                      9.81/(po_cup(i,k)-po_cup(i,k+1))
         if(k.lt.ktop(i).and.k.gt.kbcon(i))then
           dz=zo_cup(i,k+1)-zo_cup(i,k)
           dellaqc(i,k)=.01*9.81*up_massdetro(i,k)*.5*(qrco(i,k)+qrco(i,k+1))/ &
                        (po_cup(i,k)-po_cup(i,k+1))
         endif
         if(dellaqc(i,k).lt.0)write(0,*)'neg della',i,j,k,ktop(i),qrco(i,k), &
              qrco(i,k+1),up_massdetro(i,k),zuo(i,ktop(i))
         dellaqc(i,k)=max(0.,dellaqc(i,k))
         if(j.eq.jpr.and.i.eq.ipr)write(0,*)'dellaqc,qrco = ',k,dellaqc(i,k),qrco(i,k)
       endif
      enddo
      enddo
      dellah(:,:)=dellah(:,:)+dsubh(:,:)
      dellaq(:,:)=dellaq(:,:)+dsubq(:,:)
      dsubh(:,:)=0.0
      dsubq(:,:)=0.0

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
else

      do i=its,itf
        if(ierr(i).eq.0)then
         dp=100.*(po_cup(i,1)-po_cup(i,2))
         dellu(i,1)=(edto(i)*zdo(i,2)*ucd(i,2)   &
                     -edto(i)*zdo(i,2)*u_cup(i,2))*g/dp
         dellv(i,1)=(edto(i)*zdo(i,2)*vcd(i,2)   &
                     -edto(i)*zdo(i,2)*v_cup(i,2))*g/dp

         do k=kts+1,ktop(i)
            ! these three are only used at or near mass detrainment and/or entrainment levels
            entupk=0.
            detupk=0.
            entdoj=0.
            ! detrainment and entrainment for fowndrafts
            detdo=edto(i)*dd_massdetro(i,k)
            entdo=edto(i)*dd_massentro(i,k)
            ! entrainment/detrainment for updraft
            entup=up_massentro(i,k)
            detup=up_massdetro(i,k)
            ! subsidence by downdrafts only
            subin=-zdo(i,k+1)*edto(i)
            subdown=-zdo(i,k)*edto(i)
!
            !         SPECIAL LEVELS
            ! downdraft originating level, similiar to k22-1 for updraft
!-srf- 21jul2015 mass con
!            if(k.eq.jmin(i))then
!               entdoj=edto(i)*zdo(i,k)
!            endif
!-srf- mass con
            if(k.eq.ktop(i))then
               detupk=zuo(i,ktop(i))
               subin=0.
               subdown=0.
               detdo=0.
               entdo=0.
               entup=0.
               detup=0.
            endif
            totmas=subin-subdown+detup-entup-entdo+ &
                   detdo-entupk-entdoj+detupk+zuo(i,k+1)-zuo(i,k)
!              print *,'*********************',k,totmas
!              write(0,123)k,subin+zuo(i,k+1),subdown-zuo(i,k),detup,entup, &
!                          detdo,entdo,entupk,detupk
!             write(8,*)'totmas = ',k,totmas
            if(abs(totmas).gt.1.e-6)then
               write(0,*)'*********************',i,j,k,totmas
               !print *,jmin(i),k22(i),kbcon(i),ktop(i)
               write(0,123)k,subin,subdown,detup,entup, &
                           detdo,entdo,entupk,detupk
               !123     formAT(1X,i2,8E12.4)
               ! call wrf_error_fatal ( 'totmas .gt.1.e-6' )
            endif
            dp=100.*(po_cup(i,k)-po_cup(i,k+1))

            dellu(i,k)=((detup+lambau*detup)*.5*(uC(i,K+1)+uC(i,K)) &
                    +detdo*.5*(UCD(i,K+1)+UCD(i,K)) &
                    -(entup+lambau*detup)*us(i,k) &
                    -entdo*us(i,k) &
                    +subin*u_cup(i,k+1) &
                    -subdown*u_cup(i,k) &
                    +detupk*(uc(i,ktop(i))-u_cup(i,ktop(i)))    &
                    -entupk*u_cup(i,k22(i)) &
                    -entdoj*u_cup(i,jmin(i)) &
                     )*g/dp
            dellv(i,k)=((detup+lambau*detup)*.5*(vC(i,K+1)+vC(i,K)) &
                    +detdo*.5*(vCD(i,K+1)+vCD(i,K)) &
                    -(entup+lambau*detup)*vs(i,k) &
                    -entdo*vs(i,k) &
                    +subin*v_cup(i,k+1) &
                    -subdown*v_cup(i,k) &
                    +detupk*(vc(i,ktop(i))-v_cup(i,ktop(i)))    &
                    -entupk*v_cup(i,k22(i)) &
                    -entdoj*v_cup(i,jmin(i)) &
                     )*g/dp


!
! updraft subsidence only
!
           if(k.lt.ktop(i))then

             dellu(i,k)=dellu(i,k)+(zuo(i,k+1)*u_cup(i,k+1) &
                    -zuo(i,k)*u_cup(i,k) &
                    +pgcon*.5*(zuo(i,k)+zuo(i,k+1))*(u_cup(i,k+1)-u_cup(i,k)))*g/dp
             dellv(i,k)=dellv(i,k)+(zuo(i,k+1)*v_cup(i,k+1) &
                    -zuo(i,k)*v_cup(i,k) &
                    +pgcon*.5*(zuo(i,k)+zuo(i,k+1))*(v_cup(i,k+1)-v_cup(i,k)))*g/dp
           endif
       enddo   ! k

        endif
      enddo


      do i=its,itf
        trash  = 0.0
        trash2 = 0.0
        if(ierr(i).eq.0)then
         dsubq(i,:)  =0.0
         dsubh(i,:)  =0.0

         dp=100.*(po_cup(i,1)-po_cup(i,2))

         dellah(i,1)=(edto(i)*zdo(i,2)*hcdo(i,2)   &
                     -edto(i)*zdo(i,2)*heo_cup(i,2))*g/dp

         dellaqc(i,1)=0.0
         dellaq (i,1)=(edto(i)*zdo(i,2)*qcdo(i,2)   &
                      -edto(i)*zdo(i,2)*qo_cup(i,2))*g/dp

         G_rain=  0.5*(pwo (i,1)+pwo (i,2))*g/dp
         E_dn  = -0.5*(pwdo(i,1)+pwdo(i,2))*g/dp*edto(i)  ! pwdo < 0 and E_dn must > 0
         dellaq(i,1) = dellaq(i,1)+ E_dn-G_rain

         !--- conservation check
         !- water mass balance
         trash = trash  + (dellaq(i,1)+dellaqc(i,1)+G_rain-E_dn)*dp/g
         !- H  budget
         trash2 = trash2+ (dellah(i,1))*dp/g

         !write(0,*) "delH1=",1,dellah(i,k),dellaq(i,k),E_dn
         !write(3,*)'=>H k22 kbcon ktop= ',k22(i),kbcon(i),ktop(i)
         !write(3,*)'=>H= ',1,real(trash2,4),real(dellah(i,1),4)

         !write(4,*)'=>W k22 kbcon ktop= ',k22(i),kbcon(i),ktop(i)
         !write(4,*)'=>W= ',1,real(trash,4),real(dellaq(i,1),4)

         do k=kts+1,ktop(i)
            ! these three are only used at or near mass detrainment and/or entrainment levels
            entupk=0.
            detupk=0.
            entdoj=0.
            ! detrainment and entrainment for downdrafts
            detdo=edto(i)*dd_massdetro(i,k)
            entdo=edto(i)*dd_massentro(i,k)
            ! entrainment/detrainment for updraft
            entup=up_massentro(i,k)
            detup=up_massdetro(i,k)
            ! subsidence by downdrafts only
            subin=-zdo(i,k+1)*edto(i)
            subdown=-zdo(i,k)*edto(i)
            ! downdraft originating level, similiar to k22-1 for updraft

            !write(0,*)"down",k,edto(i),detdo,entdo,subin,subdown

            !if(k.eq.jmin(i))then
            !   entdoj=edto(i)*zdo(i,k)
            !endif
            !if(k.eq.ktop(i))then
            !   detupk=zuo(i,ktop(i))
            !   subin=0.
            !   subdown=0.
            !   detdo=0.
            !   entdo=0.
            !   entup=0.
            !  detup=0.
            !endif
            !totmas=subin-subdown+detup-entup-entdo+ &
            !       detdo-entupk-entdoj+detupk+zuo(i,k+1)-zuo(i,k)

            dp=100.*(po_cup(i,k)-po_cup(i,k+1))

            dellah(i,k) =-(zuo(i,k+1)*(hco (i,k+1)-heo_cup(i,k+1) ) -                 &
                           zuo(i,k  )*(hco (i,k  )-heo_cup(i,k  ) ) )*g/dp            &
                         +(zdo(i,k+1)*(hcdo(i,k+1)-heo_cup(i,k+1) ) -                 &
                           zdo(i,k  )*(hcdo(i,k  )-heo_cup(i,k  ) ) )*g/dp*edto(i)


            !- check H conservation
             trash2 = trash2+ (dellah(i,k))*dp/g


            !===tmp
            !entdoj=0.0
            !if(k.eq.jmin(i))then
            !   entdoj=edto(i)*zdo(i,k)
            !   !write(0,*)"entdoj=",zdo(i,k),zdo(i,k+1),dd_massentro(i,k),k!;stop 43
            !endif
            !===tmp

            !-- take out cloud liquid water for detrainment
            detup=up_massdetro(i,k)
            dellaqc(i,k) = detup*0.5*(qrco(i,k+1)+qrco(i,k)) *g/dp

            !--10jul latest form of Georg - to be tested (see cup_up_moisture)
            !dellaqc(i,k) = zuo(i,k)*c1*qrco(i,k)*dz/dp*g !detup*0.5*(qrco(i,k+1)+qrco(i,k)) *g/dp
            !if(k.eq.ktop(i))dellaqc(i,k)= detup*0.5*(qrco(i,k+1)+qrco(i,k)) *g/dp
            !

            !---
            G_rain=  0.5*(pwo (i,k)+pwo (i,k+1))*g/dp
            E_dn  = -0.5*(pwdo(i,k)+pwdo(i,k+1))*g/dp*edto(i) ! pwdo < 0 and E_dn must > 0
            !
            !write(0,*) "eva=",k,pwdo(i,k),E_dn,zdo(i,k  )
            !
            !-- condensation source term = detrained + flux divergence of
            !-- cloud liquid water (qrco) + converted to rain

            C_up = dellaqc(i,k)+(zuo(i,k+1)* qrco(i,k+1) -       &
                                 zuo(i,k  )* qrco(i,k  )  )*g/dp + G_rain

            !-- water vapor budget
            !-- = flux divergence z*(Q_c - Q_env)_up_and_down &
            !--   - condensation term + evaporation
            dellaq(i,k) =-(zuo(i,k+1)*(qco (i,k+1)-qo_cup(i,k+1) ) -                 &
                           zuo(i,k  )*(qco (i,k  )-qo_cup(i,k  ) ) )*g/dp            &
                         +(zdo(i,k+1)*(qcdo(i,k+1)-qo_cup(i,k+1) ) -                 &
                           zdo(i,k  )*(qcdo(i,k  )-qo_cup(i,k  ) ) )*g/dp*edto(i)    &
                         - C_up + E_dn
            !- check water conservation liq+condensed (including rainfall)
             trash= trash+ (dellaq(i,k)+dellaqc(i,k)+ G_rain-E_dn)*dp/g

         !write(3,*)'=>H= ',k,real(trash2,4),real(dellah(i,k),4)
         !write(4,*)'=>W= ',k,real(trash,4),real(dellaq(i,k),4)

         enddo   ! k
         !write(0,*)'=>H/W-FINAL= ',real(trash2,4),real(trash,4),k22(i),kbcon(i),ktop(i)
         !if(abs(trash)>1.e-6 .or. abs(trash2) > 1.e-6) then
         !  write(0,*)'=> not water mass or H cons for deep= ',i,trash,trash2
         !  !stop 33
         !endif

        endif

      enddo
!!!=====srf tmp
!      do k=kts,ktf-1
!      do i=its,itf
!         entup=up_massentro(i,k)
!         detup=up_massdetro(i,k)
!         dp=100.*(po_cup(i,k)-po_cup(i,k+1))
!
!         !erro
!         trash2 = ( -  (zuo(i,k+1)*hco(i,k+1) - zuo(i,k)*hco(i,k)     )&
!                    -  (detup*.5*(HCo(i,K+1)+HCo(i,K)) -entup*heo(i,k)) &
!                      )*g/dp
!
!       !dsubh(i,k)=0.0
!       if(ierr(i).eq.0)then
!        !write(0,*) "delH_N=",k,dellah(i,k),dsubh(i,k),trash2!,dellaq(i,k)+dellaqc(i,k)
!        write(0,*) "delH_N=",k,trash2,(zuo(i,k+1)-zuo(i,k))-entup+detup,kbcon(i)!,dellaq(i,k)+dellaqc(i,k)
!        dsubh(i,k)=0.0
!        endif
!      enddo
!      enddo
!!!=====srf tmp
endif
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
!
!--- using dellas, calculate changed environmental profiles
!
!     do 200 nens=1,maxens
      mbdt=mbdt_ens(1)
      do i=its,itf
      xaa0_ens(i,:)=0.
      enddo

!      if(j.eq.jpr)then
!               write(0,*)'xt',xl,'DELLAH(I,K),DELLAQ(I,K),dsubq(I,K),dsubt(i,k)'
!      endif
!     if(j.eq.jpr) write(0,*)'DELLAS'
      do k=kts,ktf
      do i=its,itf
         dellat(i,k)=0.
         if(ierr(i).eq.0)then
            !adding dellaqc to dellaq:
            !dellaq (i,k)= dellaq(i,k)+dellaqc(i,k)
            !dellaqc(i,k)=0.0
            !
            XHE(I,K)=(dsubh(i,k)+DELLAH(I,K))*MBDT+HEO(I,K)
            XQ (I,K)=(dsubq(i,k)+DELLAQ(I,K)+DELLAQC(i,k))*MBDT+QO(I,K)

            !- don�t feed dellat with dellaqc if
            !- the detrainment of liquid water will be used as
            !- a source for cloud microphysics
            if(COUPL_MPHYSICS) then
              DELLAT(I,K)=(1./cp)*(DELLAH(I,K)-xl*DELLAQ(I,K))
            else
	      DELLAT(I,K)=(1./cp)*(DELLAH(I,K)-xl*(DELLAQ(I,K)+DELLAQC(i,k)))
            endif
	    dSUBT(I,K)=(1./cp)*(dsubh(i,k)-xl*dsubq(i,k))
            XT(I,K)= (DELLAT(I,K)+dsubt(i,k)-xl/cp*dellaqc(i,k))*MBDT+TN(I,K)
            IF(XQ(I,K).LE.0.)XQ(I,K)=1.E-08
         ENDIF
      enddo
      enddo
      do i=its,itf
      if(ierr(i).eq.0)then
      XHKB(I)=(dsubh(i,k22(i))+DELLAH(I,K22(i)))*MBDT+HKBO(I)
      XHE(I,ktf)=HEO(I,ktf)
      XQ(I,ktf)=QO(I,ktf)
      XT(I,ktf)=TN(I,ktf)
      IF(XQ(I,ktf).LE.0.)XQ(I,ktf)=1.E-08
      endif
      enddo
!
!--- calculate moist static energy, heights, qes
!
      call cup_env(xz,xqes,xhe,xhes,xt,xq,po,z1, &
           psur,ierr,tcrit,-1,xl,cp,   &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
!
!--- environmental values on cloud levels
!
      call cup_env_clev(xt,xqes,xq,xhe,xhes,xz,po,xqes_cup,xq_cup, &
           xhe_cup,xhes_cup,xz_cup,po_cup,gamma_cup,xt_cup,psur,   &
           ierr,z1,xl,rv,cp,          &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
!
!
!**************************** static control
!
!--- moist static energy inside cloud
!
!     do i=its,itf
!       if(ierr(i).eq.0)then
!         xhkb(i)=xhe(i,k22(i))
!       endif
!     enddo
      do k=kts,ktf
      do i=its,itf
         xhc(i,k)=0.
         xDBY(I,K)=0.
      enddo
      enddo
      do i=its,itf
        if(ierr(i).eq.0)then
!        if(use_excess == 2) then
!            k1=max(1,k22(i)-1)
!            k2=max(1,min(kbcon(i)-1,k22(i)+1))
!            k1=1
!            k2=k22(i)+1
!            xhkb(i) =sum(xhe_cup(i,k1:k2))/float(k2-k1+1)+xl*zqexec(i)+cp*ztexec(i)
!        else if(use_excess <= 1) then
!        xhkb(i)=xhe_cup(i,k22(i))+float(use_excess)*(xl*zqexec(i)+cp*ztexec(i))

!       endif

!        do k=1,k22(i)
!           xhc(i,k)=xhe_cup(i,k)
!        enddo
         do k=1,kbcon(i)-1
            xhc(i,k)=xhkb(i)
         enddo
          k=kbcon(i)
          xhc(i,k)=xhkb(i)
          xDBY(I,Kbcon(i))=xHkb(I)-xHES_cup(I,K)
        endif !ierr
      enddo
!
!
      do i=its,itf
      if(ierr(i).eq.0)then
!ss     xzu(i,:)=zuo(i,:)
!-- mass con
!     do k=kbcon(i)+1,ktop(i)  ! orig
      do k=k22(i)  +1,ktop(i)  ! mass cons option
!-- mass con
       xhc(i,k)=(xhc(i,k-1)*xzu(i,k-1)-.5*up_massdetro(i,k-1)*xhc(i,k-1)+ &
                         up_massentro(i,k-1)*xhe(i,k-1))   /            &
                         (xzu(i,k-1)-.5*up_massdetro(i,k-1)+up_massentro(i,k-1))
       xdby(i,k)=xhc(i,k)-xhes_cup(i,k)
      enddo
      do k=ktop(i)+1,ktf
           xHC(i,K)=xhes_cup(i,k)
           xDBY(I,K)=0.
           xzu(i,k)=0.
      enddo
      endif
      enddo

!
!--- workfunctions for updraft
!
      call cup_up_aa0(xaa0,xz,xzu,xdby,GAMMA_CUP,xt_cup, &
           kbcon,ktop,ierr,           &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
      if(j.eq.jpr)write(0,*)'xaa0,aa0,aa1 = ',xaa0(1),aa0(1),aa1(1)
      do 200 nens=1,maxens
      do i=its,itf
         if(ierr(i).eq.0)then
           xaa0_ens(i,nens)=xaa0(i)
           nall=(iens-1)*maxens3*maxens*maxens2 &
                +(iedt-1)*maxens*maxens3 &
                +(nens-1)*maxens3
           do k=kts,ktf
              if(k.le.ktop(i))then
                 do nens3=1,maxens3
                 if(nens3.eq.7)then
!--- b=0
                 pr_ens(i,j,nall+nens3)=pr_ens(i,j,nall+nens3)  &
!                                +edto(i)*pwdo(i,k)             &
                                    +pwo(i,k)
!--- b=beta
                 else if(nens3.eq.8)then
                 pr_ens(i,j,nall+nens3)=pr_ens(i,j,nall+nens3)+ &
                                    pwo(i,k)
!--- b=beta/2
                 else if(nens3.eq.9)then
                 pr_ens(i,j,nall+nens3)=pr_ens(i,j,nall+nens3)  &
!                                +.5*edto(i)*pwdo(i,k)          &
                                 +  pwo(i,k)
                 else
                 pr_ens(i,j,nall+nens3)=pr_ens(i,j,nall+nens3)+ &
                                    pwo(i,k) ! +edto(i)*pwdo(i,k)
                 endif
                 enddo
              endif
           enddo
         if(pr_ens(i,j,nall+7).lt.1.e-6)then
            ierr(i)=18
            ierrc(i)="total normalized condensate too small"
            if(i.eq.ipr.and.j.eq.jpr)write(0,*)ierr(i),ierrc(i)
            do nens3=1,maxens3
               pr_ens(i,j,nall+nens3)=0.
            enddo
         endif
         do nens3=1,maxens3
           if(pr_ens(i,j,nall+nens3).lt.1.e-4)then
            pr_ens(i,j,nall+nens3)=0.
           endif
         enddo
         endif
      if(i.eq.ipr.and.j.eq.jpr)write(0,*)'ierrc = ',ierr(i),ierrc(i)
      enddo
 200  continue
!
!--- LARGE SCALE FORCING
!
!
!------- CHECK wether aa0 should have been zero, assuming this
!        ensemble is chosen
!
!
      do i=its,itf
         ierr2(i)=ierr(i)
         ierr3(i)=ierr(i)
         k22x(i)=k22(i)
      enddo
      IF(maxens.gt.0)then
       call cup_MAXIMI(HEO_CUP,3,KBMAX,K22x,ierr, &
!srf21112015       call cup_MAXIMI(HEO_CUP,1,KBMAX,K22x,ierr, &
            itf,jtf,ktf, &
            its,ite, jts,jte, kts,kte)
       call cup_kbcon(ierrc,cap_max_increment,2,k22x,kbconx,heo_cup, &
            heso_cup,hkbo,ierr2,kbmax,po_cup,cap_max, &
            xl,cp,ztexec,zqexec,use_excess,       &
            0,itf,jtf,ktf, &
            its,ite, jts,jte, kts,kte, &
            z_cup,entr_rate,heo)
       call cup_kbcon(ierrc,cap_max_increment,3,k22x,kbconx,heo_cup, &
            heso_cup,hkbo,ierr3,kbmax,po_cup,cap_max, &
            xl,cp,ztexec,zqexec,use_excess,       &
            0,itf,jtf,ktf, &
            its,ite, jts,jte, kts,kte, &
            z_cup,entr_rate,heo)
      ENDIF
!
!--- calculate cloud base mass flux
!

      call cup_forcing_ens_3d(closure_n,xland1,aa0,aa1,xaa0_ens,mbdt_ens,dtime,   &
           ierr,ierr2,ierr3,xf_ens,j,'deeps',axx,                 &
           maxens,iens,iedt,maxens2,maxens3,mconv,            &
           po_cup,ktop,omeg,zdo,k22,zuo,pr_ens,edto,kbcon,    &
           ensdim,ichoice,     &
           imid,ipr,jpr,itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte,ens4,ktau, &
           dicycle,tau_ecmwf,aa1_bl,xf_dicycle)
!     call cup_forcing_ens_3d(closure_n,xland1,aa0,aa1,xaa0_ens,mbdt_ens,dtime,   &
!          ierr,ierr2,ierr3,xf_ens,j,'deeps',axx,                 &
!          tot_time_hr,maxens,iens,iedt,maxens2,maxens3,mconv,            &
!          po_cup,ktop,omeg,zdo,k22,zuo,pr_ens,edto,kbcon,    &
!          ensdim,ichoice,     &
!          ipr,jpr,itf,jtf,ktf, &
!          its,ite, jts,jte, kts,kte,ens4,ktau)
!
      do k=kts,ktf
      do i=its,itf
        if(ierr(i).eq.0)then
           subt_ens(i,k,iedt)=dsubt(i,k)
           subq_ens(i,k,iedt)=dsubq(i,k)
           dellat_ens(i,k,iedt)=dellat(i,k)
           dellaq_ens(i,k,iedt)=dellaq(i,k)
           dellaqc_ens(i,k,iedt)=dellaqc(i,k)
           pwo_ens(i,k,iedt)=pwo(i,k)+edto(i)*pwdo(i,k)
        else
           subt_ens(i,k,iedt)=0.
           subq_ens(i,k,iedt)=0.
           dellat_ens(i,k,iedt)=0.
           dellaq_ens(i,k,iedt)=0.
           dellaqc_ens(i,k,iedt)=0.
           pwo_ens(i,k,iedt)=0.
        endif
        if(i.eq.ipr.and.j.eq.jpr)then
          write(0,*)'2',k,subt_ens(i,k,iedt),subq_ens(i,k,iedt)
        endif
      enddo
      enddo
 250  continue
!
!--- FEEDBACK
!
       IF(imid.eq.1)THEN
         do i=its,itf
          xff_mid(i,1)=0.
          xff_mid(i,2)=0.
          if(ierr(i).eq.0)then
          blqe=0.
          trash=0.
          if(k22(i).lt.kpbl(i)+1)then
             do k=1,kpbl(i)
                blqe=blqe+100.*dhdt(i,k)*(po_cup(i,k)-po_cup(i,k+1))/g
             enddo
             trash=max((hco(i,kbcon(i))-heo_cup(i,kbcon(i))),1.e1)
             xff_mid(i,1)=max(0.,blqe/trash)
             xff_mid(i,1)=min(0.1,xff_mid(i,1))
          endif
          xff_mid(i,2)=.03*zws(i)
          endif
         enddo
       ENDIF
       call cup_output_ens_3d(xff_mid,xf_ens,ierr,dellat_ens,dellaq_ens, &
            dellaqc_ens,subt_ens,subq_ens,subt,subq,outt,     &
            outq,outqc,zuo,sub_mas,pre,pwo_ens,xmb,ktop,      &
            j,'deep',maxens2,maxens,iens,ierr2,ierr3,         &
            pr_ens,maxens3,ensdim,                    &
            sig,APR_GR,APR_W,APR_MC,APR_ST,APR_AS,                &
            APR_CAPMA,APR_CAPME,APR_CAPMI,closure_n,xland1,   &
            weight_GR,weight_W,weight_MC,weight_ST,weight_AS,training, &
            ichoice,imid,ipr,jpr,itf,jtf,ktf,                        &
            its,ite, jts,jte, kts,kte, &
            dicycle,xf_dicycle )
      k=1
      do i=its,itf
          if(ierr(i).eq.0) then
!            outu(i,ktf)=edto(i)
!            outv(i,ktf)=xmb(i)
             PRE(I)=MAX(PRE(I),0.)
             xmb_out(i)=xmb(i)
             !print*,"precip deep=",pre(i),xmb(i)
             do k=kts,ktop(i)
               outu(i,k)=dellu(i,k)*xmb(i)
               outv(i,k)=dellv(i,k)*xmb(i)
               sub_mas(i,k)=zuo(i,k)*xmb(i)
             enddo
          else
             ktop(i)=0
          endif
      enddo
!     do i=its,itf
!         if(ierr(i).eq.0) then
!            dts=0.
!            fpi=0.
!            do k=kts+1,ktop(i)
!               dp=p_cup(i,k)-p_cup(i,k-1)
!               dts= dts +(outu(i,k)*us(i,k)+outv(i,k)*vs(i,k))*dp/g
!               fpi = fpi  +sqrt(outu(i,k)*outu(i,k) + outv(i,k)*outv(i,k))*dp
!            enddo
!            if(fpi.gt.0.)then
!            do k=kts+1,ktop(i)
!               fp= sqrt((outu(i,k)*outu(i,k)+outv(i,k)*outv(i,k)))/fpi/cp*g*dts
!               outt(i,k)=outt(i,k)+fp
!            enddo
!            endif
!         endif
!     enddo

!------- temporary use of apr_* arrays to save several cloud work
!------- functions for debug purposes
    IF(training == 0) then
     do i=its,itf
            apr_w (i,j)=aa0(i)
            apr_st(i,j)=aa1(i)
            apr_gr(i,j)=xaa0(i)
            apr_mc(i,j)=(XAA0(I)-AA1(I))/MBDT
            apr_as(i,j)=aa1_bl(i)
         !if(ierr(i).eq.0) then
         !endif
     enddo
    ENDIF
!------------------------- temporary


!
!- for tracer convective transport-start
 !-sub-grid scale ice/liquid cloud water for radiation
 do i=its,itf
   ierr4d(i)   = ierr(i)
   xmb4d(i)    =  xmb(i)
   do k=1,ktf
        zup5d  (k,i) = zuo (i,k)
        clwup5d(k,i) = qrco(i,k) ! ice/liquid water
 enddo;enddo
 !---data saved for the tracer convective transport:
    do i=its,itf
       jmin4d(i)   = jmin(i)
       kdet4d(i)   = kdet(i)
       k224d(i)    = k22(i)
       kbcon4d(i)  = kbcon(i)
       ktop4d(i)   = ktop(i)
       kstabi4d(i) = kstabi(i)
       kstabm4d(i) = kstabm(i)
       !kpbl4d(i)   = kpbl(i)
       kpbl4d(i)   =  k22(i)
       edt4d(i)    =  edto(i)
       pwav4d(i)   = pwavo(i)
     do k=1,ktf
        !zcup5d(k,i) = zo_cup(i,k)
        pcup5d(k,i) = po_cup(i,k)
        up_massentr5d(k,i) = up_massentro(i,k)
        up_massdetr5d(k,i) = up_massdetro(i,k)
        dd_massentr5d(k,i) = dd_massentro(i,k)
        dd_massdetr5d(k,i) = dd_massdetro(i,k)
        zdn5d(k,i)  = zdo(i,k)

        prup5d (k,i) =   pwo(i,k) !for updraft
        prdn5d (k,i) =   pwdo(i,k)!for downdraft
        tup5d  (k,i) = t_cup(i,k)
        !>>> em verdade deveria ser a temperatura da parcela
        !>>> de ar no updraft e _NAO_ a temperatura ambiente
     enddo
    enddo
!- for convective transport-end
!---------------------------done------------------------------
!
!- begin: for GATE soundings-------------------------------------------
  if(use_gate) then
    do i=its,itf
     if(ierr(i).eq.0) then
      !- 2-d section
      do k=kts,ktop(i)

       nvar=1
       call setgradsvar(jl,k,nvar,xmb(i)*zuo(i,k),"mup",'m flux up (kg/s/m^2)')

       nvar=2
       call setgradsvar(jl,k,nvar,-edto(i)*xmb(i)*zdo(i,k),"mdn" ,'m flux down (kg/s/m^2)')

       nvar=3
       call setgradsvar(jl,k,nvar,dellah(i,k),"delh" ,'dellah')

       nvar=4
       call setgradsvar(jl,k,nvar,dellaq(i,k)*86400.*xl/cp*xmb(i), "dellq" ,'dellaq')

       nvar=5
       call setgradsvar(jl,k,nvar,dellaqc(i,k)*86400.*xl/cp*xmb(i),"dellqc" ,'dellaqc')

       nvar=6
       call setgradsvar(jl,k,nvar,up_massentro(i,k),"upent" ,'up_massentr')

       nvar=7
       call setgradsvar(jl,k,nvar,up_massdetro(i,k),"updet" ,'up_massdetr')

       nvar=8
       call setgradsvar(jl,k,nvar,outt(i,k)*86400.,"outt" ,'outt K/s')

       nvar=9
       call setgradsvar(jl,k,nvar,outq(i,k)*86400.*xl/cp,"outq" ,'outq K/s')

       nvar=10
       call setgradsvar(jl,k,nvar,outqc(i,k)*86400.*xl/cp,"outqc" ,'outqc K/s')

       nvar=11
       call setgradsvar(jl,k,nvar,pre(i)*3600.,"precip" ,'precip mm')

       nvar=12
       call setgradsvar(jl,k,nvar,outu(i,k),"outu" ,'out_U m/s/s')

       nvar=13
       call setgradsvar(jl,k,nvar,outv(i,k),"outqv" ,'out_V m/s/s')

       !outvars(i,k,7)=up_massentro(i,k)
       !outvars(i,k,8)=up_massdetro(i,k)
      enddo
     endif
    enddo
  endif
!- end  : for GATE soundings-------------------------------------------


   END SUBROUTINE CUP_gf
!========================================================================================
   SUBROUTINE SETGRADSVAR(i,k,nvar,f,name1,name2)
     implicit none
     integer, intent(in) :: nvar,i,k
     real, intent(in) :: f
     character*(*) :: name1,name2

     cupout(nvar)%varp(i,k)= f
     cupout(nvar)%varn(1)=name1
     cupout(nvar)%varn(2)=name2

   END SUBROUTINE SETGRADSVAR
!========================================================================================
   SUBROUTINE cup_dd_edt(ierr,us,vs,z,ktop,kbcon,edt,p,pwav, &
              pw,ccn,pwev,edtmax,edtmin,maxens2,edtc,psum2,psumh, &
              ccnclean,rho,aeroevap,itf,jtf,ktf,j,ipr,jpr,          &
              its,ite, jts,jte, kts,kte                     )

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
        j,ipr,jpr,aeroevap,itf,jtf,ktf,           &
        its,ite, jts,jte, kts,kte
     integer, intent (in   )              ::                           &
        maxens2
  !
  ! ierr error value, maybe modified in this routine
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        rho,us,vs,z,p,pw
     real,    dimension (its:ite,1:maxens2)                            &
        ,intent (out  )                   ::                           &
        edtc
     real,    dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
        edt
     real,    dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        pwav,pwev,ccn,psum2,psumh
     real                                                              &
        ,intent (in   )                   ::                           &
        ccnclean,edtmax,edtmin
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        ktop,kbcon
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
!
!  local variables in this routine
!

     integer i,k,kk
     real    einc,pef,pefb,prezk,zkbc
     real,    dimension (its:ite)         ::                           &
      vshear,sdp,vws
     real :: prop_c,pefc,aeroadd,alpha3,beta3,rhoc
     prop_c=8. !10.386
     alpha3 = 1.9
     beta3  = -1.13
     pefc=0.

!
!--- DETERMINE DOWNDRAFT STRENGTH IN TERMS OF WINDSHEAR
!
! */ calculate an average wind shear over the depth of the cloud
!
       do i=its,itf
        edt(i)=0.
        vws(i)=0.
        sdp(i)=0.
        vshear(i)=0.
       enddo
       do k=1,maxens2
       do i=its,itf
        edtc(i,k)=0.
       enddo
       enddo
       do kk = kts,ktf-1
         do 62 i=its,itf
          IF(ierr(i).ne.0)GO TO 62
          if (kk .le. min0(ktop(i),ktf) .and. kk .ge. kbcon(i)) then
             vws(i) = vws(i)+ &
              (abs((us(i,kk+1)-us(i,kk))/(z(i,kk+1)-z(i,kk))) &
          +   abs((vs(i,kk+1)-vs(i,kk))/(z(i,kk+1)-z(i,kk)))) * &
              (p(i,kk) - p(i,kk+1))
            sdp(i) = sdp(i) + p(i,kk) - p(i,kk+1)
          endif
          if (kk .eq. ktf-1)vshear(i) = 1.e3 * vws(i) / sdp(i)
   62   continue
       end do
      do i=its,itf
         IF(ierr(i).eq.0)then
            pef=(1.591-.639*VSHEAR(I)+.0953*(VSHEAR(I)**2) &
               -.00496*(VSHEAR(I)**3))
            if(pef.gt.0.9)pef=0.9
            if(pef.lt.0.1)pef=0.1
!
!--- cloud base precip efficiency
!
            zkbc=z(i,kbcon(i))*3.281e-3
            prezk=.02
            if(zkbc.gt.3.)then
               prezk=.96729352+zkbc*(-.70034167+zkbc*(.162179896+zkbc &
               *(- 1.2569798E-2+zkbc*(4.2772E-4-zkbc*5.44E-6))))
            endif
            if(zkbc.gt.25)then
               prezk=2.4
            endif
            pefb=1./(1.+prezk)
!           write(11,*)pefb,prezk,zkbc
            if(pefb.gt.0.9)pefb=0.9
            if(pefb.lt.0.1)pefb=0.1
            EDT(I)=1.-.5*(pefb+pef)
            if(aeroevap.gt.1)then
               aeroadd=(ccnclean**beta3)*((psumh(i))**(alpha3-1)) !*1.e6
               if(i.eq.ipr.and.j.eq.jpr)write(0,*)'edt',ccnclean,psumh(i),aeroadd
!              prop_c=.9/aeroadd
               prop_c=.5*(pefb+pef)/aeroadd
               aeroadd=(ccn(i)**beta3)*((psum2(i))**(alpha3-1)) !*1.e6
               if(i.eq.ipr.and.j.eq.jpr)write(0,*)'edt',ccn(i),psum2(i),aeroadd,prop_c
               aeroadd=prop_c*aeroadd
               pefc=aeroadd
               if(pefc.gt.0.9)pefc=0.9
               if(pefc.lt.0.1)pefc=0.1
               EDT(I)=1.-pefc
               if(aeroevap.eq.2)EDT(I)=1.-.25*(pefb+pef+2.*pefc)
            endif


!--- edt here is 1-precipeff!
            einc=.2*edt(i)
            do k=1,maxens2
                edtc(i,k)=edt(i)+float(k-2)*einc
            enddo
         endif
      enddo
      do i=its,itf
         IF(ierr(i).eq.0)then
            do k=1,maxens2
               EDTC(I,K)=-EDTC(I,K)*pwav(I)/PWEV(I)
               IF(EDTC(I,K).GT.edtmax)EDTC(I,K)=edtmax
               IF(EDTC(I,K).LT.edtmin)EDTC(I,K)=edtmin
            enddo
         endif
      enddo

   END SUBROUTINE cup_dd_edt


   SUBROUTINE cup_dd_moisture_new(ierrc,zd,hcd,hes_cup,qcd,qes_cup,    &
              pwd,q_cup,z_cup,dd_massentr,dd_massdetr,jmin,ierr,            &
              gamma_cup,pwev,bu,qrcd,                        &
              q,he,t_cup,iloop,xl,           &
              itf,jtf,ktf,                     &
              its,ite, jts,jte, kts,kte                     )

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
                                  itf,jtf,ktf,           &
                                  its,ite, jts,jte, kts,kte
  ! cdd= detrainment function
  ! q = environmental q on model levels
  ! q_cup = environmental q on model cloud levels
  ! qes_cup = saturation q on model cloud levels
  ! hes_cup = saturation h on model cloud levels
  ! hcd = h in model cloud
  ! bu = buoancy term
  ! zd = normalized downdraft mass flux
  ! gamma_cup = gamma on model cloud levels
  ! mentr_rate = entrainment rate
  ! qcd = cloud q (including liquid water) after entrainment
  ! qrch = saturation q in cloud
  ! pwd = evaporate at that level
  ! pwev = total normalized integrated evaoprate (I2)
  ! entr= entrainment rate
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        zd,t_cup,hes_cup,hcd,qes_cup,q_cup,z_cup,                      &
        dd_massentr,dd_massdetr,gamma_cup,q,he
     real                                                              &
        ,intent (in   )                   ::                           &
        xl
     integer                                                           &
        ,intent (in   )                   ::                           &
        iloop
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        jmin
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (out  )                   ::                           &
        qcd,qrcd,pwd
     real,    dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
        pwev,bu
     character*50 :: ierrc(its:ite)
!
!  local variables in this routine
!

     integer                              ::                           &
        i,k,ki
     real                                 ::                           &
        dh,dz,dqeva

      do i=its,itf
         bu(i)=0.
         pwev(i)=0.
      enddo
      do k=kts,ktf
      do i=its,itf
         qcd(i,k)=0.
         qrcd(i,k)=0.
         pwd(i,k)=0.
      enddo
      enddo
!
!
!
      do 100 i=its,itf
      IF(ierr(I).eq.0)then
      k=jmin(i)
      DZ=Z_cup(i,K+1)-Z_cup(i,K)
      qcd(i,k)=q_cup(i,k)
      DH=HCD(I,k)-HES_cup(I,K)
      if(dh.lt.0)then
        QRCD(I,K)=(qes_cup(i,k)+(1./XL)*(GAMMA_cup(i,k) &
                  /(1.+GAMMA_cup(i,k)))*DH)
      else
          qrcd(i,k)=qes_cup(i,k)
      endif
      pwd(i,jmin(i))=zd(i,jmin(i))*min(0.,qcd(i,k)-qrcd(i,k))
      qcd(i,k)=qrcd(i,k)
      pwev(i)=pwev(i)+pwd(i,jmin(i)) ! *dz
!
      bu(i)=dz*dh
      do ki=jmin(i)-1,1,-1
         DZ=Z_cup(i,Ki+1)-Z_cup(i,Ki)
!        QCD(i,Ki)=(qCD(i,Ki+1)*(1.-.5*CDD(i,Ki+1)*DZ) &
!                 +entr*DZ*q(i,Ki) &
!                )/(1.+entr*DZ-.5*CDD(i,Ki+1)*DZ)
!        dz=qcd(i,ki)
         qcd(i,ki)=(qcd(i,ki+1)*zd(i,ki+1)                          &
                  -.5*dd_massdetr(i,ki)*qcd(i,ki+1)+ &
                  dd_massentr(i,ki)*q(i,ki))   /            &
                  (zd(i,ki+1)-.5*dd_massdetr(i,ki)+dd_massentr(i,ki))
!        write(0,*)'qcd in dd_moi = ',qcd(i,ki)

!
!--- to be negatively buoyant, hcd should be smaller than hes!
!--- ideally, dh should be negative till dd hits ground, but that is not always
!--- the case
!
         DH=HCD(I,ki)-HES_cup(I,Ki)
         bu(i)=bu(i)+dz*dh
         QRCD(I,Ki)=qes_cup(i,ki)+(1./XL)*(GAMMA_cup(i,ki) &
                  /(1.+GAMMA_cup(i,ki)))*DH
         dqeva=qcd(i,ki)-qrcd(i,ki)
         if(dqeva.gt.0.)then
          dqeva=0.
          qrcd(i,ki)=qcd(i,ki)
         endif
         pwd(i,ki)=zd(i,ki)*dqeva
         qcd(i,ki)=qrcd(i,ki)
         pwev(i)=pwev(i)+pwd(i,ki) ! *dz
!        if(iloop.eq.1.and.i.eq.102.and.j.eq.62)then
!         print *,'in cup_dd_moi ', hcd(i,ki),HES_cup(I,Ki),dh,dqeva
!        endif
      enddo
!
!--- end loop over i
       if(pwev(I).eq.0.and.iloop.eq.1)then
!        print *,'problem with buoy in cup_dd_moisture',i
         ierr(i)=7
         ierrc(i)="problem with buoy in cup_dd_moisture"
       endif
       if(BU(I).GE.0.and.iloop.eq.1)then
!        print *,'problem with buoy in cup_dd_moisture',i
         ierr(i)=7
         ierrc(i)="problem2 with buoy in cup_dd_moisture"
       endif
      endif
100    continue

   END SUBROUTINE cup_dd_moisture_new

   SUBROUTINE cup_env(z,qes,he,hes,t,q,p,z1,                 &
              psur,ierr,tcrit,itest,xl,cp,                   &
              itf,jtf,ktf,                     &
              its,ite, jts,jte, kts,kte                     )

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
        itf,jtf,ktf,           &
        its,ite, jts,jte, kts,kte
  !
  ! ierr error value, maybe modified in this routine
  ! q           = environmental mixing ratio
  ! qes         = environmental saturation mixing ratio
  ! t           = environmental temp
  ! tv          = environmental virtual temp
  ! p           = environmental pressure
  ! z           = environmental heights
  ! he          = environmental moist static energy
  ! hes         = environmental saturation moist static energy
  ! psur        = surface pressure
  ! z1          = terrain elevation
  !
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        p,t,q
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (out  )                   ::                           &
        he,hes,qes
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (inout)                   ::                           &
        z
     real,    dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        psur,z1
     real                                                              &
        ,intent (in   )                   ::                           &
        xl,cp
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
     integer                                                           &
        ,intent (in   )                   ::                           &
        itest
!
!  local variables in this routine
!

     integer                              ::                           &
       i,k,iph
      real, dimension (1:2) :: AE,BE,HT
      real, dimension (its:ite,kts:kte) :: tv
      real :: tcrit,e,tvbar
!      real, external :: satvap
!      real :: satvap


      HT(1)=XL/CP
      HT(2)=2.834E6/CP
      BE(1)=.622*HT(1)/.286
      AE(1)=BE(1)/273.+ALOG(610.71)
      BE(2)=.622*HT(2)/.286
      AE(2)=BE(2)/273.+ALOG(610.71)
!      print *, 'TCRIT = ', tcrit,its,ite
      DO k=kts,ktf
      do i=its,itf
        if(ierr(i).eq.0)then
!Csgb - IPH is for phase, dependent on TCRIT (water or ice)
        IPH=1
        IF(T(I,K).LE.TCRIT)IPH=2
!       print *, 'AE(IPH),BE(IPH) = ',AE(IPH),BE(IPH),AE(IPH)-BE(IPH),T(i,k),i,k
!       E=EXP(AE(IPH)-BE(IPH)/T(I,K))
!       print *, 'P, E = ', P(I,K), E
!       QES(I,K)=.622*E/(100.*P(I,K)-E)
        e=satvap(t(i,k))
        qes(i,k)=0.622*e/max(1.e-8,(p(i,k)-e))
        IF(QES(I,K).LE.1.E-08)QES(I,K)=1.E-08
        IF(QES(I,K).LT.Q(I,K))QES(I,K)=Q(I,K)
!       IF(Q(I,K).GT.QES(I,K))Q(I,K)=QES(I,K)
        TV(I,K)=T(I,K)+.608*Q(I,K)*T(I,K)
        endif
      enddo
      enddo
!
!--- z's are calculated with changed h's and q's and t's
!--- if itest=2
!
      if(itest.eq.1 .or. itest.eq.0)then
         do i=its,itf
           if(ierr(i).eq.0)then
             Z(I,1)=max(0.,Z1(I))-(ALOG(P(I,1))- &
                 ALOG(PSUR(I)))*287.*TV(I,1)/9.81
           endif
         enddo

! --- calculate heights
         DO K=kts+1,ktf
         do i=its,itf
           if(ierr(i).eq.0)then
              TVBAR=.5*TV(I,K)+.5*TV(I,K-1)
              Z(I,K)=Z(I,K-1)-(ALOG(P(I,K))- &
               ALOG(P(I,K-1)))*287.*TVBAR/9.81
           endif
         enddo
         enddo
      else if(itest.eq.2)then
         do k=kts,ktf
         do i=its,itf
           if(ierr(i).eq.0)then
             z(i,k)=(he(i,k)-1004.*t(i,k)-2.5e6*q(i,k))/9.81
             z(i,k)=max(1.e-3,z(i,k))
           endif
         enddo
         enddo
      else if(itest.eq.-1)then
      endif
!
!--- calculate moist static energy - HE
!    saturated moist static energy - HES
!
       DO k=kts,ktf
       do i=its,itf
         if(ierr(i).eq.0)then
         if(itest.le.0)HE(I,K)=9.81*Z(I,K)+1004.*T(I,K)+2.5E06*Q(I,K)
         HES(I,K)=9.81*Z(I,K)+1004.*T(I,K)+2.5E06*QES(I,K)
         IF(HE(I,K).GE.HES(I,K))HE(I,K)=HES(I,K)
         endif
      enddo
      enddo

   END SUBROUTINE cup_env


   SUBROUTINE cup_env_clev(t,qes,q,he,hes,z,p,qes_cup,q_cup,   &
              he_cup,hes_cup,z_cup,p_cup,gamma_cup,t_cup,psur, &
              ierr,z1,xl,rv,cp,                                &
              itf,jtf,ktf,                       &
              its,ite, jts,jte, kts,kte                       )

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
        itf,jtf,ktf,           &
        its,ite, jts,jte, kts,kte
  !
  ! ierr error value, maybe modified in this routine
  ! q           = environmental mixing ratio
  ! q_cup       = environmental mixing ratio on cloud levels
  ! qes         = environmental saturation mixing ratio
  ! qes_cup     = environmental saturation mixing ratio on cloud levels
  ! t           = environmental temp
  ! t_cup       = environmental temp on cloud levels
  ! p           = environmental pressure
  ! p_cup       = environmental pressure on cloud levels
  ! z           = environmental heights
  ! z_cup       = environmental heights on cloud levels
  ! he          = environmental moist static energy
  ! he_cup      = environmental moist static energy on cloud levels
  ! hes         = environmental saturation moist static energy
  ! hes_cup     = environmental saturation moist static energy on cloud levels
  ! gamma_cup   = gamma on cloud levels
  ! psur        = surface pressure
  ! z1          = terrain elevation
  !
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        qes,q,he,hes,z,p,t
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (out  )                   ::                           &
        qes_cup,q_cup,he_cup,hes_cup,z_cup,p_cup,gamma_cup,t_cup
     real,    dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        psur,z1
     real                                                              &
        ,intent (in   )                   ::                           &
        xl,rv,cp
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
!
!  local variables in this routine
!

     integer                              ::                           &
       i,k


      do k=kts,ktf
      do i=its,itf
        qes_cup(i,k)=0.
        q_cup(i,k)=0.
        hes_cup(i,k)=0.
        he_cup(i,k)=0.
        z_cup(i,k)=0.
        p_cup(i,k)=0.
        t_cup(i,k)=0.
        gamma_cup(i,k)=0.
      enddo
      enddo
      do k=kts+1,ktf
      do i=its,itf
        if(ierr(i).eq.0)then
        qes_cup(i,k)=.5*(qes(i,k-1)+qes(i,k))
        q_cup(i,k)=.5*(q(i,k-1)+q(i,k))
        hes_cup(i,k)=.5*(hes(i,k-1)+hes(i,k))
        he_cup(i,k)=.5*(he(i,k-1)+he(i,k))
        if(he_cup(i,k).gt.hes_cup(i,k))he_cup(i,k)=hes_cup(i,k)
        z_cup(i,k)=.5*(z(i,k-1)+z(i,k))
        p_cup(i,k)=.5*(p(i,k-1)+p(i,k))
        t_cup(i,k)=.5*(t(i,k-1)+t(i,k))
        gamma_cup(i,k)=(xl/cp)*(xl/(rv*t_cup(i,k) &
                       *t_cup(i,k)))*qes_cup(i,k)
        endif
      enddo
      enddo
      do i=its,itf
        if(ierr(i).eq.0)then
        qes_cup(i,1)=qes(i,1)
        q_cup(i,1)=q(i,1)
!       hes_cup(i,1)=hes(i,1)
!       he_cup(i,1)=he(i,1)
        hes_cup(i,1)=9.81*z1(i)+1004.*t(i,1)+2.5e6*qes(i,1)
        he_cup(i,1)=9.81*z1(i)+1004.*t(i,1)+2.5e6*q(i,1)
        z_cup(i,1)=.5*(z(i,1)+z1(i))
        p_cup(i,1)=.5*(p(i,1)+psur(i))
        z_cup(i,1)=z1(i)
        p_cup(i,1)=psur(i)
        t_cup(i,1)=t(i,1)
        gamma_cup(i,1)=xl/cp*(xl/(rv*t_cup(i,1) &
                       *t_cup(i,1)))*qes_cup(i,1)
        endif
      enddo

   END SUBROUTINE cup_env_clev

   SUBROUTINE cup_forcing_ens_3d(closure_n,xland,aa0,aa1,xaa0,mbdt,dtime,ierr,ierr2,ierr3,&
              xf_ens,j,name,axx,maxens,iens,iedt,maxens2,maxens3,mconv,    &
              p_cup,ktop,omeg,zd,k22,zu,pr_ens,edt,kbcon,      &
              ensdim,icoic,            &
              imid,ipr,jpr,itf,jtf,ktf,               &
              its,ite, jts,jte, kts,kte,ens4,ktau, &
              dicycle,tau_ecmwf,aa1_bl,xf_dicycle  )

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
        imid,ipr,jpr,itf,jtf,ktf,           &
        its,ite, jts,jte, kts,kte,ens4,ktau
     integer, intent (in   )              ::                           &
        j,ensdim,maxens,iens,iedt,maxens2,maxens3
  !
  ! ierr error value, maybe modified in this routine
  ! pr_ens = precipitation ensemble
  ! xf_ens = mass flux ensembles
  ! massfln = downdraft mass flux ensembles used in next timestep
  ! omeg = omega from large scale model
  ! mconv = moisture convergence from large scale model
  ! zd      = downdraft normalized mass flux
  ! zu      = updraft normalized mass flux
  ! aa0     = cloud work function without forcing effects
  ! aa1     = cloud work function with forcing effects
  ! xaa0    = cloud work function with cloud effects (ensemble dependent)
  ! edt     = epsilon
  ! dir     = "storm motion"
  ! mbdt    = arbitrary numerical parameter
  ! dtime   = dt over which forcing is applied
  ! iact_gr_old = flag to tell where convection was active
  ! kbcon       = LFC of parcel from k22
  ! k22         = updraft originating level
  ! icoic       = flag if only want one closure (usually set to zero!)
  ! name        = deep or shallow convection flag
  !
     real,    dimension (its:ite,jts:jte,1:ensdim)                     &
        ,intent (inout)                   ::                           &
        pr_ens
     real,    dimension (its:ite,jts:jte,1:ensdim)                     &
        ,intent (out  )                   ::                           &
        xf_ens
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        zd,zu,p_cup
     real,    dimension (its:ite,kts:kte,1:ens4)                              &
        ,intent (in   )                   ::                           &
        omeg
     real,    dimension (its:ite,1:maxens)                             &
        ,intent (in   )                   ::                           &
        xaa0
     real,    dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        aa1,edt,xland
     real,    dimension (its:ite,1:ens4)                                      &
        ,intent (in   )                   ::                           &
        mconv,axx
     real,    dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        aa0,closure_n
     real,    dimension (1:maxens)                                     &
        ,intent (in   )                   ::                           &
        mbdt
     real                                                              &
        ,intent (in   )                   ::                           &
        dtime
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        k22,kbcon,ktop
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr,ierr2,ierr3
     integer                                                           &
        ,intent (in   )                   ::                           &
        icoic
      character *(*), intent (in)         ::                           &
       name
!-srf begin
      logical, intent(IN) :: DICYCLE
      real,    intent(IN)   , dimension (its:ite) :: aa1_bl,tau_ecmwf
      real,    intent(INOUT), dimension (its:ite) :: xf_dicycle
      !- local var
      real  :: xff_dicycle
!-srf end
!
!  local variables in this routine
!

     real,    dimension (1:maxens3)       ::                           &
       xff_ens3
     real,    dimension (1:maxens)        ::                           &
       xk
     integer                              ::                           &
       i,k,nall,n,ne,nens,nens3,iresult,iresultd,iresulte,mkxcrt,kclim
     parameter (mkxcrt=15)
     real                                 ::                           &
       fens4,a1,massfld,a_ave,xff0,xff00,xxx,xomg,aclim1,aclim2,aclim3,aclim4
     real,    dimension(1:mkxcrt)         ::                           &
       pcrit,acrit,acritt

     integer :: nall2,ixxx,irandom
     integer,  dimension (8) :: seed
     real, dimension (its:ite) :: ens_adj


      DATA PCRIT/850.,800.,750.,700.,650.,600.,550.,500.,450.,400.,    &
                 350.,300.,250.,200.,150./
      DATA ACRIT/.0633,.0445,.0553,.0664,.075,.1082,.1521,.2216,       &
                 .3151,.3677,.41,.5255,.7663,1.1686,1.6851/
!  GDAS DERIVED ACRIT
      DATA ACRITT/.203,.515,.521,.566,.625,.665,.659,.688,             &
                  .743,.813,.886,.947,1.138,1.377,1.896/

!
       ens_adj=1.
       seed=0
       do i=its,itf
        if(ierr(i).eq.0)then
          seed(1)=int(aa0(i))
          seed(2)=int(aa1(i))
          exit
        endif
       enddo

       nens=0
       irandom=0
       fens4=float(ens4)

!--- LARGE SCALE FORCING
!
       DO 100 i=its,itf
          if(name.eq.'deeps'.and.ierr(i).gt.995)then
           aa0(i)=0.
           ierr(i)=0
          endif
          IF(ierr(i).eq.0)then
             ens_adj(i)=1.
             if(ierr2(i).gt.0.and.ierr3(i).eq.0)ens_adj(i)=0. ! 2./3.
             if(ierr2(i).gt.0.and.ierr3(i).gt.0)ens_adj(i)=0.
!
!---
!
             if(name.eq.'deeps')then
!
                a_ave=0.
                do ne=1,ens4
                  a_ave=a_ave+axx(i,ne)
!               if(i.eq.ipr.and.j.eq.jpr)write(0,*)'in forcing, a_ave,axx(i,ne)
!               = ',a_ave,axx(i,ne),maxens,xland(i)
                enddo
                a_ave=max(0.,a_ave/fens4)
                a_ave=min(a_ave,aa1(i))
                a_ave=max(0.,a_ave)
                do ne=1,16
                  xff_ens3(ne)=0.
                enddo
                xff0= (AA1(I)-AA0(I))/DTIME
                xff_ens3(1)=max(0.,(AA1(I)-AA0(I))/dtime)
                xff_ens3(2)=max(0.,(AA1(I)-AA0(I))/dtime)
!               xff_ens3(2)=max(0.,(a_ave-AA0(I))/dtime)

!               if(i.eq.ipr.and.j.eq.jpr)write(0,*)AA1(I),AA0(I),xff_ens3(1),xff_ens3(2),dtime
                if(irandom.eq.1)then
                   call random_number (xxx)
                   ixxx=min(ens4,max(1,int(fens4*xxx+1.e-8)))
                   xff_ens3(3)=max(0.,(axx(i,ixxx)-AA0(I))/dtime)
                   call random_number (xxx)
                   ixxx=min(ens4,max(1,int(fens4*xxx+1.e-8)))
                   xff_ens3(13)=max(0.,(axx(i,ixxx)-AA0(I))/dtime)
                else
                   xff_ens3(3)=max(0.,(AA1(I)-AA0(I))/dtime)
                   xff_ens3(13)=max(0.,(AA1(I)-AA0(I))/dtime)
                endif
!
!--- omeg is in bar/s, mconv done with omeg in Pa/s
!     more like Brown (1979), or Frank-Cohen (199?)
!
                xff_ens3(14)=0.
                do ne=1,ens4
                  xff_ens3(14)=xff_ens3(14)-omeg(i,k22(i),ne)/(fens4*9.81)
                enddo
                if(xff_ens3(14).lt.0.)xff_ens3(14)=0.
                xff_ens3(5)=0.
                do ne=1,ens4
                  xff_ens3(5)=xff_ens3(5)-omeg(i,kbcon(i),ne)/(fens4*9.81)
                enddo
                if(xff_ens3(5).lt.0.)xff_ens3(5)=0.
!
! minimum below kbcon
!
                xff_ens3(4)=-omeg(i,2,1)/9.81
                do k=2,kbcon(i)-1
                   do ne=1,ens4
                     xomg=-omeg(i,k,ne)/9.81
                     if(xomg.lt.xff_ens3(4))xff_ens3(4)=xomg
                   enddo
                enddo
                if(xff_ens3(4).lt.0.)xff_ens3(4)=0.
!
! max below kbcon
                xff_ens3(6)=-omeg(i,2,1)/9.81
                do k=2,kbcon(i)-1
                   do ne=1,ens4
                     xomg=-omeg(i,k,ne)/9.81
                     if(xomg.gt.xff_ens3(6))xff_ens3(6)=xomg
                   enddo
                enddo
                if(xff_ens3(6).lt.0.)xff_ens3(6)=0.
                xff_ens3(5)=xff_ens3(6)
                xff_ens3(4)=xff_ens3(6)
!               if(i.eq.ipr.and.j.eq.jpr)write(0,*)xff_ens3(4),xff_ens3(5)
!
!--- more like Krishnamurti et al.; pick max and average values
!
                xff_ens3(7)=mconv(i,1)
                xff_ens3(8)=mconv(i,1)
                xff_ens3(9)=mconv(i,1)
                if(ens4.gt.1)then
                   do ne=2,ens4
                      if (mconv(i,ne).gt.xff_ens3(7))xff_ens3(7)=mconv(i,ne)
                   enddo
                   do ne=2,ens4
                      if (mconv(i,ne).lt.xff_ens3(8))xff_ens3(8)=mconv(i,ne)
                   enddo
                   do ne=2,ens4
                      xff_ens3(9)=xff_ens3(9)+mconv(i,ne)
                   enddo
                   xff_ens3(9)=xff_ens3(9)/fens4
                endif
!               if(i.eq.ipr.and.j.eq.jpr)write(0,*)xff_ens3(7),xff_ens3(8)
!
                if(irandom.eq.1)then
                   call random_number (xxx)
                   ixxx=min(ens4,max(1,int(fens4*xxx+1.e-8)))
                   xff_ens3(15)=mconv(i,ixxx)
                else
                   xff_ens3(15)=mconv(i,1)
                endif
!
!--- more like Fritsch Chappel or Kain Fritsch (plus triggers)
!
                xff_ens3(10)=AA0(i)/(60.*20.)
                xff_ens3(11)=AA0(I)/(60.*20.)
                xff_ens3(16)=AA0(I)/(60.*20.)
                if(irandom.eq.1)then
                call random_number (xxx)
                   ixxx=min(ens4,max(1,int(fens4*xxx+1.e-8)))
                   xff_ens3(12)=AA0(I)/(60.*20.)
                else
                   xff_ens3(12)=AA0(I)/(60.*20.)
                endif
!
!-srf begin
!- more like Bechtold et al. (JAS 2014)
                !IF(DICYCLE) then
                ! xff_ens3(13)=(AA0(i)-AA1_BL(i))/tau_ecmwf(i)
                ! xff_ens3(14)=(AA0(i)-AA1_BL(i))/tau_ecmwf(i)
                ! xff_ens3(15)=(AA0(i)-AA1_BL(i))/tau_ecmwf(i)
                !ENDIF
!-srf implementation of 2014
!                xff_ens3(13)=(AA0(i))/tau_ecmwf(i)
!                xff_ens3(14)=(AA0(i))/tau_ecmwf(i)
!                xff_ens3(15)=(AA0(i))/tau_ecmwf(i)
!-srf implementation of 2015 - AA1 already appplied the forcing
                xff_ens3(13)=(AA1(i))/tau_ecmwf(i)
                xff_ens3(14)=(AA1(i))/tau_ecmwf(i)
                xff_ens3(15)=(AA1(i))/tau_ecmwf(i)

                xff_dicycle = AA1_BL(i)/tau_ecmwf(i)
!-srf end
!
!gtest
                if(icoic.eq.0)then
                if(xff0.lt.0.)then
                     xff_ens3(1)=0.
                     xff_ens3(2)=0.
                     xff_ens3(3)=0.
                     !srf if(.not. dicycle) xff_ens3(13)=0.
                     xff_ens3(10)=0.
                     xff_ens3(11)=0.
                     xff_ens3(12)=0.
                     !-srf begin
                     xff_ens3(13)= 0.
                     xff_ens3(14)= 0.
                     xff_ens3(15)= 0.
                     xff_dicycle = 0.
                     !-srf end
                endif
                  if(xff0.lt.0 .and. xland(i).lt.0.1 .and. imid.eq.0)then
                     xff_ens3(:)=0.
                  endif
                endif

!                  if(i.eq.ipr.and.j.eq.jpr)write(0,*)'xff_ens =
!                  ',i,j,ipr,jpr,xff_ens3


                do nens=1,maxens
                   XK(nens)=(XAA0(I,nens)-AA1(I))/MBDT(1)
!                  if(i.eq.ipr.and.j.eq.jpr)write(0,*)'xks =
!                  ',xk(nens),XAA0(I,nens),AA1(I),mbdt
                   !if(xk(nens).le.0.and.xk(nens).gt.-1.e-2) &
                   !        xk(nens)=-1.e-2
                   !lixo if(xk(nens).le.0.and.xk(nens).gt.-10.*mbdt(1)) &
                   !lixo        xk(nens)=-10.*mbdt(1)
                   if(xk(nens).le.0.and.xk(nens).gt.-.1*mbdt(1)) &
                           xk(nens)=-.1*mbdt(1)

                   if(xk(nens).gt.0.and.xk(nens).lt.1.e-2) &
                           xk(nens)=1.e-2
                enddo
!
!--- add up all ensembles
!
                do 350 ne=1,maxens
!
!--- for every xk, we have maxens3 xffs
!--- iens is from outermost ensemble (most expensive!
!
!--- iedt (maxens2 belongs to it)
!--- is from second, next outermost, not so expensive
!
!--- so, for every outermost loop, we have maxens*maxens2*3
!--- ensembles!!! nall would be 0, if everything is on first
!--- loop index, then ne would start counting, then iedt, then iens....
!
                   iresult=0
                   iresultd=0
                   iresulte=0
                   nall=(iens-1)*maxens3*maxens*maxens2 &
                        +(iedt-1)*maxens*maxens3 &
                        +(ne-1)*maxens3
!                 if(i.eq.ipr.and.j.eq.jpr)write(0,*)'maxens',ne,nall,iens,maxens3,maxens,maxens2,iedt
!
! over water, enfor!e small cap for some of the closures
!
                if(maxens.gt.0 .and. xland(i).lt.0.1)then
                 if(ierr2(i).gt.0.or.ierr3(i).gt.0)then
                      xff_ens3(1) =ens_adj(i)*xff_ens3(1)
                      xff_ens3(2) =ens_adj(i)*xff_ens3(2)
                      xff_ens3(3) =ens_adj(i)*xff_ens3(3)
                      xff_ens3(13) =ens_adj(i)*xff_ens3(13)
                      xff_ens3(10) =ens_adj(i)*xff_ens3(10)
                      xff_ens3(11) =ens_adj(i)*xff_ens3(11)
                      xff_ens3(12) =ens_adj(i)*xff_ens3(12)
                      xff_ens3(16) =ens_adj(i)*xff_ens3(16)
                      xff_ens3(7) =ens_adj(i)*xff_ens3(7)
                      xff_ens3(8) =ens_adj(i)*xff_ens3(8)
                      xff_ens3(9) =ens_adj(i)*xff_ens3(9)
                      xff_ens3(15) =ens_adj(i)*xff_ens3(15)
                      !srf
                      xff_dicycle = ens_adj(i)*xff_dicycle
                      !srf end
!                     xff_ens3(7) =0.
!                     xff_ens3(8) =0.
!                     xff_ens3(9) =0.
                 endif
                endif
!
! end water treatment
!
!
!--- check for upwind convection
!                  iresult=0
                   massfld=0.

                   IF(XK(ne).lt.0.and.xff0.gt.0.)iresultd=1
                   iresulte=max(iresult,iresultd)
                   iresulte=1
                   if(iresulte.eq.1)then
!
!--- special treatment for stability closures
!
!                      if(i.eq.ipr.and.j.eq.jpr)write(0,*)'xffs =
!                      ',xff_ens3(1:16)

                      if(xff0.gt.0.)then
                         if(xff_ens3(1).gt.0)xf_ens(i,j,nall+1)=max(0.,-xff_ens3(1)/xk(ne))
                         if(xff_ens3(2).gt.0)xf_ens(i,j,nall+2)=max(0.,-xff_ens3(2)/xk(ne))
                         if(xff_ens3(3).gt.0)xf_ens(i,j,nall+3)=max(0.,-xff_ens3(3)/xk(ne))
                         !srf if(.not. dicycle) then
                         !srf  if(xff_ens3(13).gt.0)xf_ens(i,j,nall+13)=max(0.,-xff_ens3(13)/xk(ne))
                         !srf endif

                      else
                         xff_ens3(1)=0
                         xff_ens3(2)=0
                         xff_ens3(3)=0
                         !srf if(.not. dicycle)xff_ens3(13)=0
                      endif
!
!--- if iresult.eq.1, following independent of xff0
!
                         xf_ens(i,j,nall+4)=max(0.,xff_ens3(4))
                         xf_ens(i,j,nall+5)=max(0.,xff_ens3(5))
                         xf_ens(i,j,nall+6)=max(0.,xff_ens3(6))
                         !srf if(.not. dicycle) xf_ens(i,j,nall+14)=max(0.,xff_ens3(14))

                         a1=max(1.e-3,pr_ens(i,j,nall+7))
                         xf_ens(i,j,nall+7)=max(0.,xff_ens3(7)/a1)
                         a1=max(1.e-3,pr_ens(i,j,nall+8))
                         xf_ens(i,j,nall+8)=max(0.,xff_ens3(8)/a1)
                         a1=max(1.e-3,pr_ens(i,j,nall+9))
                         xf_ens(i,j,nall+9)=max(0.,xff_ens3(9)/a1)
                         !srf a1=max(1.e-3,pr_ens(i,j,nall+15))
                         !srf if(.not. dicycle) xf_ens(i,j,nall+15)=max(0.,xff_ens3(15)/a1)

                         if(XK(ne).lt.0.)then
                            xf_ens(i,j,nall+10)=max(0.,-xff_ens3(10)/xk(ne))
                            xf_ens(i,j,nall+11)=max(0.,-xff_ens3(11)/xk(ne))
                            xf_ens(i,j,nall+12)=max(0.,-xff_ens3(12)/xk(ne))
                            xf_ens(i,j,nall+16)=max(0.,-xff_ens3(16)/xk(ne))
                         endif
!srf-begin
                         if(XK(ne).lt.0.)then
                            xf_ens(i,j,nall+13)=max(0.,-xff_ens3(13)/xk(ne))
                            xf_ens(i,j,nall+14)=max(0.,-xff_ens3(14)/xk(ne))
                            xf_ens(i,j,nall+15)=max(0.,-xff_ens3(15)/xk(ne))
                            !ask georg if the line below should be here -OR- out of this if statement
                            !xf_dicycle(i)      =  max(0.,-xff_dicycle /xk(ne))
                         endif

                         !testing out
                         xf_dicycle(i)     =    max(0.,-xff_dicycle /xk(ne))

!srf-end
                      if(icoic.ge.1)then
                      closure_n(i)=0.
                      xf_ens(i,j,nall+1)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+2)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+3)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+4)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+5)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+6)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+7)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+8)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+9)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+10)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+11)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+12)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+13)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+14)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+15)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+16)=xf_ens(i,j,nall+icoic)
                      endif
                     !print*,"xf_dic2=",xf_dicycle(i),xff_dicycle,maxval(xf_ens(i,j,:));call flush(6)
!
! 16 is a randon pick from the oher 15
!
                if(irandom.eq.1)then
                   call random_number (xxx)
                   ixxx=min(15,max(1,int(15.*xxx+1.e-8)))
                   xf_ens(i,j,nall+16)=xf_ens(i,j,nall+ixxx)
!               else
!                  xf_ens(i,j,nall+16)=xf_ens(i,j,nall+1)
                endif
!
!
!--- do some more on the caps!!! ne=1 for 175, ne=2 for 100,....
!
!     do not care for caps here for closure groups 1 and 5,
!     they are fine, do not turn them off here
!
!!!!    NOT USED FOR "NORMAL" APPLICATION (maxens=1)
!
                if(maxens.gt.1)then
                if(ne.eq.2.and.ierr2(i).gt.0)then
                      xf_ens(i,j,nall+1) =0.
                      xf_ens(i,j,nall+2) =0.
                      xf_ens(i,j,nall+3) =0.
                      xf_ens(i,j,nall+4) =0.
                      xf_ens(i,j,nall+5) =0.
                      xf_ens(i,j,nall+6) =0.
                      xf_ens(i,j,nall+7) =0.
                      xf_ens(i,j,nall+8) =0.
                      xf_ens(i,j,nall+9) =0.
                      xf_ens(i,j,nall+10)=0.
                      xf_ens(i,j,nall+11)=0.
                      xf_ens(i,j,nall+12)=0.
                      xf_ens(i,j,nall+13)=0.
                      xf_ens(i,j,nall+14)=0.
                      xf_ens(i,j,nall+15)=0.
                      xf_ens(i,j,nall+16)=0.
                endif
                if(ne.eq.3.and.ierr3(i).gt.0)then
                      xf_ens(i,j,nall+1) =0.
                      xf_ens(i,j,nall+2) =0.
                      xf_ens(i,j,nall+3) =0.
                      xf_ens(i,j,nall+4) =0.
                      xf_ens(i,j,nall+5) =0.
                      xf_ens(i,j,nall+6) =0.
                      xf_ens(i,j,nall+7) =0.
                      xf_ens(i,j,nall+8) =0.
                      xf_ens(i,j,nall+9) =0.
                      xf_ens(i,j,nall+10)=0.
                      xf_ens(i,j,nall+11)=0.
                      xf_ens(i,j,nall+12)=0.
                      xf_ens(i,j,nall+13)=0.
                      xf_ens(i,j,nall+14)=0.
                      xf_ens(i,j,nall+15)=0.
                      xf_ens(i,j,nall+16)=0.
                endif
                endif

                   endif
 350            continue
                if(maxens.gt.1)then
! ne=1, cap=175
!
                   nall=(iens-1)*maxens3*maxens*maxens2 &
                        +(iedt-1)*maxens*maxens3
! ne=2, cap=100
!
                   nall2=(iens-1)*maxens3*maxens*maxens2 &
                        +(iedt-1)*maxens*maxens3 &
                        +(2-1)*maxens3
                      xf_ens(i,j,nall+4) = xf_ens(i,j,nall2+4)
                      xf_ens(i,j,nall+5) =xf_ens(i,j,nall2+5)
                      xf_ens(i,j,nall+6) =xf_ens(i,j,nall2+6)
                      xf_ens(i,j,nall+14) =xf_ens(i,j,nall2+14)
                      xf_ens(i,j,nall+7) =xf_ens(i,j,nall2+7)
                      xf_ens(i,j,nall+8) =xf_ens(i,j,nall2+8)
                      xf_ens(i,j,nall+9) =xf_ens(i,j,nall2+9)
                      xf_ens(i,j,nall+15) =xf_ens(i,j,nall2+15)
                      xf_ens(i,j,nall+10)=xf_ens(i,j,nall2+10)
                      xf_ens(i,j,nall+11)=xf_ens(i,j,nall2+11)
                      xf_ens(i,j,nall+12)=xf_ens(i,j,nall2+12)
!                     if(i.eq.ipr.and.j.eq.jpr)write(0,*)'should not be here'
                   endif
!         if(i.eq.ipr.and.j.eq.jpr)write(0,*)'xff_ens3=',xff_ens3
                go to 100
             endif
          elseif(ierr(i).ne.20.and.ierr(i).ne.0)then
             do n=1,ensdim
               xf_ens(i,j,n)=0.
               !srf begin
               xf_dicycle(i) = 0.
               !srf end
             enddo
          endif
!         if(i.eq.ipr.and.j.eq.jpr)write(0,*)'xff_ens3=',xff_ens3
 100   continue

   END SUBROUTINE cup_forcing_ens_3d

   SUBROUTINE cup_kbcon_old(ierrc,cap_inc,iloop,k22,kbcon,he_cup,hes_cup, &
              hkb,ierr,kbmax,p_cup,cap_max,                         &
              xl,cp,ztexec,zqexec,use_excess,       &
              jprnt,itf,jtf,ktf,                        &
              its,ite, jts,jte, kts,kte                        )

   IMPLICIT NONE
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
        use_excess,jprnt,itf,jtf,ktf,           &
        its,ite, jts,jte, kts,kte
  !
  !
  !
  ! ierr error value, maybe modified in this routine
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        he_cup,hes_cup,p_cup
     real,    dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        ztexec,zqexec,cap_max,cap_inc
     real,intent (in   )                  ::                           &
        xl,cp
     real,    dimension (its:ite)                                      &
        ,intent (inout   )                   ::                           &
        hkb
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        kbmax
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        kbcon,k22,ierr
     integer                                                           &
        ,intent (in   )                   ::                           &
        iloop
     character*50 :: ierrc(its:ite)
!
!  local variables in this routine
!

     integer                              ::                           &
        i,k,k1,k2
     real                                 ::                           &
        pbcdif,plus,hetest
!
!--- DETERMINE THE LEVEL OF CONVECTIVE CLOUD BASE  - KBCON
!
       DO 27 i=its,itf
      kbcon(i)=1
      IF(ierr(I).ne.0)GO TO 27
      KBCON(I)=K22(I)+1
      if(iloop.eq.5)KBCON(I)=K22(I)
      GO TO 32
 31   CONTINUE
      KBCON(I)=KBCON(I)+1
      IF(KBCON(I).GT.KBMAX(i)+2)THEN
         if(iloop.ne.4)then
                ierr(i)=3
                ierrc(i)="could not find reasonable kbcon in cup_kbcon"
         endif
        GO TO 27
      ENDIF
 32   CONTINUE
      hetest=hkb(i) ! HE_cup(I,K22(I))
      if(iloop.eq.5)then
       hetest=HKB(I)
!      do k=1,k22(i)
!        hetest=max(hetest,he_cup(i,k))
!      enddo
      endif
      IF(HETEST.LT.HES_cup(I,KBCON(I)))then
        if(jprnt.eq.1)write(0,*)'htest',k22(i),kbcon(i),HETEST,-P_cup(I,KBCON(I))+P_cup(I,K22(I))
        GO TO 31
      endif

!     cloud base pressure and max moist static energy pressure
!     i.e., the depth (in mb) of the layer of negative buoyancy
      if(KBCON(I)-K22(I).eq.1)go to 27
      if(iloop.eq.5 .and. (KBCON(I)-K22(I)).le.2)go to 27
      PBCDIF=-P_cup(I,KBCON(I))+P_cup(I,K22(I))
      plus=max(25.,cap_max(i)-float(iloop-1)*cap_inc(i))
      if(iloop.eq.4)plus=cap_max(i)
!
! for shallow convection, if cap_max is greater than 25, it is the pressure at pbltop
      if(iloop.eq.5)plus=150.
      if(iloop.eq.5.and.cap_max(i).gt.25)pbcdif=-P_cup(I,KBCON(I))+cap_max(i)
      IF(PBCDIF.GT.plus)THEN
       if(jprnt.eq.1)write(0,*)'htest',k22(i),kbcon(i),plus,-P_cup(I,KBCON(I))+P_cup(I,K22(I))
        K22(I)=K22(I)+1
        KBCON(I)=K22(I)+1
         if(use_excess == 2) then
             k1=max(1,k22(i)-1)
             k2=max(1,min(kbcon(i)-1,k22(i)+1))  !kbcon(i)-1
             k2=k22(i)+1
             hkb(i) =sum(he_cup(i,k1:k2))/float(k2-k1+1)!+(xl*zqexec(i)+cp*ztexec(i))/float(k2-k1+1)
        else if(use_excess <= 1)then
         hkb(i)=he_cup(i,k22(i)) !+float(use_excess)*(xl*zqexec(i)+cp*ztexec(i))
        endif  ! excess

        if(iloop.eq.5)KBCON(I)=K22(I)
        IF(KBCON(I).GT.KBMAX(i)+2)THEN
            if(iloop.ne.4)then
                ierr(i)=3
                ierrc(i)="could not find reasonable kbcon in cup_kbcon"
            endif
            GO TO 27
        ENDIF
        GO TO 32
      ENDIF
 27   CONTINUE

   END SUBROUTINE cup_kbcon_old


   SUBROUTINE cup_ktop(ierrc,ilo,dby,kbcon,ktop,ierr,              &
              itf,jtf,ktf,                     &
              its,ite, jts,jte, kts,kte                     )

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
        itf,jtf,ktf,           &
        its,ite, jts,jte, kts,kte
  ! dby = buoancy term
  ! ktop = cloud top (output)
  ! ilo  = flag
  ! ierr error value, maybe modified in this routine
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (inout)                   ::                           &
        dby
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        kbcon
     integer                                                           &
        ,intent (in   )                   ::                           &
        ilo
     integer, dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
        ktop
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
     character*50 :: ierrc(its:ite)
!
!  local variables in this routine
!

     integer                              ::                           &
        i,k
!
        DO 42 i=its,itf
        ktop(i)=1
         IF(ierr(I).EQ.0)then
          DO 40 K=KBCON(I)+1,ktf-1
            IF(DBY(I,K).LE.0.)THEN
                KTOP(I)=K-1
                GO TO 41
             ENDIF
  40      CONTINUE
          if(ilo.eq.1)ierr(i)=5
          if(ilo.eq.1)ierrc(i)="problem with defining ktop"
!         if(ilo.eq.2)ierr(i)=998
          GO TO 42
  41     CONTINUE
         do k=ktop(i)+1,ktf
           dby(i,k)=0.
         enddo
         if(kbcon(i).eq.ktop(i))then
            ierr(i)=55
            ierrc(i)="kbcon == ktop "
         endif
         endif
  42     CONTINUE

   END SUBROUTINE cup_ktop


   SUBROUTINE cup_MAXIMI(ARRAY,KS,KE,MAXX,ierr,              &
              itf,jtf,ktf,                     &
              its,ite, jts,jte, kts,kte                     )

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
         itf,jtf,ktf,                                    &
         its,ite, jts,jte, kts,kte
  ! array input array
  ! x output array with return values
  ! kt output array of levels
  ! ks,kend  check-range
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
         array
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
         ierr,ke
     integer                                                           &
        ,intent (in   )                   ::                           &
         ks
     integer, dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
         maxx
     real,    dimension (its:ite)         ::                           &
         x
     real                                 ::                           &
         xar
     integer                              ::                           &
         i,k

       DO 200 i=its,itf
       MAXX(I)=KS
       if(ierr(i).eq.0)then
      X(I)=ARRAY(I,KS)
!
       DO 100 K=KS,KE(i)
         XAR=ARRAY(I,K)
         IF(XAR.GE.X(I)) THEN
            X(I)=XAR
            MAXX(I)=K
         ENDIF
 100  CONTINUE
      endif
 200  CONTINUE

   END SUBROUTINE cup_MAXIMI


   SUBROUTINE cup_minimi(ARRAY,KS,KEND,KT,ierr,              &
              itf,jtf,ktf,                     &
              its,ite, jts,jte, kts,kte                     )

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
         itf,jtf,ktf,                                    &
         its,ite, jts,jte, kts,kte
  ! array input array
  ! x output array with return values
  ! kt output array of levels
  ! ks,kend  check-range
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
         array
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
         ierr,ks,kend
     integer, dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
         kt
     real,    dimension (its:ite)         ::                           &
         x
     integer                              ::                           &
         i,k,kstop

       DO 200 i=its,itf
      KT(I)=KS(I)
      if(ierr(i).eq.0)then
      X(I)=ARRAY(I,KS(I))
       KSTOP=MAX(KS(I)+1,KEND(I))
!
       DO 100 K=KS(I)+1,KSTOP
         IF(ARRAY(I,K).LT.X(I)) THEN
              X(I)=ARRAY(I,K)
              KT(I)=K
         ENDIF
 100  CONTINUE
      endif
 200  CONTINUE

   END SUBROUTINE cup_MINIMI


   SUBROUTINE cup_up_aa0(aa0,z,zu,dby,GAMMA_CUP,t_cup,       &
              kbcon,ktop,ierr,                               &
              itf,jtf,ktf,                     &
              its,ite, jts,jte, kts,kte                     )

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
        itf,jtf,ktf,                                     &
        its,ite, jts,jte, kts,kte
  ! aa0 cloud work function
  ! gamma_cup = gamma on model cloud levels
  ! t_cup = temperature (Kelvin) on model cloud levels
  ! dby = buoancy term
  ! zu= normalized updraft mass flux
  ! z = heights of model levels
  ! ierr error value, maybe modified in this routine
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        z,zu,gamma_cup,t_cup,dby
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        kbcon,ktop
!
! input and output
!


     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
     real,    dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
        aa0
!
!  local variables in this routine
!

     integer                              ::                           &
        i,k
     real                                 ::                           &
        dz,da
!
        do i=its,itf
         aa0(i)=0.
        enddo
        DO k=kts+1,ktf
            DO i=its,itf
               IF(ierr(i).ne.0) cycle !GO TO 100
               IF(K.LE.KBCON(I)) cycle !GO TO 100
               IF(K.Gt.KTOP(I)) cycle !GO TO 100
               DZ=Z(I,K)-Z(I,K-1)
               da=zu(i,k)*DZ*(9.81/(1004.*( &
               (T_cup(I,K)))))*DBY(I,K-1)/ &
               (1.+GAMMA_CUP(I,K))
               IF(K.eq.KTOP(I).and.da.le.0.) cycle!go to 100
               AA0(I)=AA0(I)+da
               if(aa0(i).lt.0.)aa0(i)=0.
            end do
         end do
!100     continue

   END SUBROUTINE cup_up_aa0

!====================================================================
   SUBROUTINE g3init(RTHCUTEN,RQVCUTEN,RQCCUTEN,RQICUTEN,           &
                        MASS_FLUX,cp,restart,                       &
                        P_QC,P_QI,P_FIRST_SCALAR,                   &
                        RTHFTEN, RQVFTEN,                           &
                        APR_GR,APR_W,APR_MC,APR_ST,APR_AS,          &
                        APR_CAPMA,APR_CAPME,APR_CAPMI,              &
                        cugd_tten,cugd_ttens,cugd_qvten,            &
                        cugd_qvtens,cugd_qcten,                     &
                        allowed_to_read,                            &
                        ids, ide, jds, jde, kds, kde,               &
                        ims, ime, jms, jme, kms, kme,               &
                        its, ite, jts, jte, kts, kte               )
!--------------------------------------------------------------------
   IMPLICIT NONE
!--------------------------------------------------------------------
   LOGICAL , INTENT(IN)           ::  restart,allowed_to_read
   INTEGER , INTENT(IN)           ::  ids, ide, jds, jde, kds, kde, &
                                      ims, ime, jms, jme, kms, kme, &
                                      its, ite, jts, jte, kts, kte
   INTEGER , INTENT(IN)           ::  P_FIRST_SCALAR, P_QI, P_QC
   REAL,     INTENT(IN)           ::  cp

   REAL,     DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(OUT) ::       &
                                                          CUGD_TTEN,         &
                                                          CUGD_TTENS,        &
                                                          CUGD_QVTEN,        &
                                                          CUGD_QVTENS,       &
                                                          CUGD_QCTEN
   REAL,     DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(OUT) ::       &
                                                          RTHCUTEN, &
                                                          RQVCUTEN, &
                                                          RQCCUTEN, &
                                                          RQICUTEN

   REAL,     DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(OUT) ::       &
                                                          RTHFTEN,  &
                                                          RQVFTEN

   REAL,     DIMENSION( ims:ime , jms:jme ) , INTENT(OUT) ::        &
                                APR_GR,APR_W,APR_MC,APR_ST,APR_AS,  &
                                APR_CAPMA,APR_CAPME,APR_CAPMI,      &
                                MASS_FLUX

   INTEGER :: i, j, k, itf, jtf, ktf

   jtf=min0(jte,jde-1)
   ktf=min0(kte,kde-1)
   itf=min0(ite,ide-1)

   IF(.not.restart)THEN
     DO j=jts,jte
     DO k=kts,kte
     DO i=its,ite
        RTHCUTEN(k,i,j)=0.
        RQVCUTEN(k,i,j)=0.
     ENDDO
     ENDDO
     ENDDO
     DO j=jts,jte
     DO k=kts,kte
     DO i=its,ite
       cugd_tten(k,i,j)=0.
       cugd_ttens(k,i,j)=0.
       cugd_qvten(k,i,j)=0.
       cugd_qvtens(k,i,j)=0.
     ENDDO
     ENDDO
     ENDDO

     DO j=jts,jtf
     DO k=kts,ktf
     DO i=its,itf
        RTHFTEN(k,i,j)=0.
        RQVFTEN(k,i,j)=0.
     ENDDO
     ENDDO
     ENDDO

     IF (P_QC .ge. P_FIRST_SCALAR) THEN
        DO j=jts,jtf
        DO k=kts,ktf
        DO i=its,itf
           RQCCUTEN(k,i,j)=0.
           cugd_qcten(k,i,j)=0.
        ENDDO
        ENDDO
        ENDDO
     ENDIF

     IF (P_QI .ge. P_FIRST_SCALAR) THEN
        DO j=jts,jtf
        DO k=kts,ktf
        DO i=its,itf
           RQICUTEN(k,i,j)=0.
        ENDDO
        ENDDO
        ENDDO
     ENDIF

     DO j=jts,jtf
     DO i=its,itf
        mass_flux(i,j)=0.
     ENDDO
     ENDDO

     DO j=jts,jtf
     DO i=its,itf
        APR_GR(i,j)=0.
        APR_ST(i,j)=0.
        APR_W(i,j)=0.
        APR_MC(i,j)=0.
        APR_AS(i,j)=0.
        APR_CAPMA(i,j)=0.
        APR_CAPME(i,j)=0.
        APR_CAPMI(i,j)=0.
     ENDDO
     ENDDO

   ENDIF

   END SUBROUTINE g3init

   SUBROUTINE neg_check(j,subt,subq,dt,q,outq,outt,outqc,pret,its,ite,kts,kte,itf,ktf)

   INTEGER,      INTENT(IN   ) ::            j,its,ite,kts,kte,itf,ktf

     real, dimension (its:ite,kts:kte  )                    ,                 &
      intent(inout   ) ::                                                     &
       outq,outt,outqc,subt,subq
     real, dimension (its:ite,kts:kte  )                    ,                 &
      intent(inout   ) ::                                                     &
       q
     real, dimension (its:ite  )                            ,                 &
      intent(inout   ) ::                                                     &
       pret
     real                                                                     &
        ,intent (in  )                   ::                                   &
        dt
     real :: thresh,qmem,qmemf,qmem2,qtest,qmem1
     integer :: i,k
!
! first do check on vertical heating rate
!
!      thresh=300.01
     thresh=100.01 !srf 29fev2016
      do i=its,itf
      qmemf=1.
      qmem=0.
      do k=kts,ktf
         qmem=(subt(i,k)+outt(i,k))*86400.
         if(qmem.gt.2.*thresh)then
           qmem2=2.*thresh/qmem
           qmemf=min(qmemf,qmem2)
!
!
!          print *,'1',' adjusted massflux by factor ',i,j,k,qmem,qmem2,qmemf,dt
         endif
         if(qmem.lt.-thresh)then
           qmem2=-thresh/qmem
           qmemf=min(qmemf,qmem2)
!
!
!          print *,'2',' adjusted massflux by factor ',i,j,k,qmem,qmem2,qmemf,dt
         endif
      enddo
!     if(qmemf.lt.1)then
!          write(0,*)'1',' adjusted massflux by factor ',i,j,qmemf
!     endif
      do k=kts,ktf
         subq(i,k)=subq(i,k)*qmemf
         subt(i,k)=subt(i,k)*qmemf
         outq(i,k)=outq(i,k)*qmemf
         outt(i,k)=outt(i,k)*qmemf
         outqc(i,k)=outqc(i,k)*qmemf
      enddo
      pret(i)=pret(i)*qmemf
      enddo
!
! check whether routine produces negative q's. This can happen, since
! tendencies are calculated based on forced q's. This should have no
! influence on conservation properties, it scales linear through all
! tendencies
!
      thresh=1.e-10
      do i=its,itf
      qmemf=1.
      do k=kts,ktf-1
         qmem=subq(i,k)+outq(i,k)
         if(abs(qmem).gt.0.)then
         qtest=q(i,k)+(subq(i,k)+outq(i,k))*dt
         if(qtest.lt.thresh)then
!
! qmem2 would be the maximum allowable tendency
!
           qmem1=outq(i,k)+subq(i,k)
           qmem2=(thresh-q(i,k))/dt
           qmemf=min(qmemf,qmem2/qmem1)
           qmemf=max(0.,qmemf)
!          write(0,*)'4 adjusted tendencies ',i,k,qmem,qmem2,qmemf
!          write(0,*)'4 adjusted tendencies ',i,j,k,q(i,k),qmem1,qmemf
         endif
         endif
      enddo
!     if(qmemf.lt.1.)write(0,*)'4 adjusted tendencies ',i,j,qmemf
      do k=kts,ktf
         subq(i,k)=subq(i,k)*qmemf
         subt(i,k)=subt(i,k)*qmemf
         outq(i,k)=outq(i,k)*qmemf
         outt(i,k)=outt(i,k)*qmemf
         outqc(i,k)=outqc(i,k)*qmemf
      enddo
      pret(i)=pret(i)*qmemf
      enddo

   END SUBROUTINE neg_check


   SUBROUTINE cup_output_ens_3d(xff_mid,xf_ens,ierr,dellat,dellaq,dellaqc,  &
              subt_ens,subq_ens,subt,subq,outtem,outq,outqc,     &
              zu,sub_mas,pre,pw,xmb,ktop,                 &
              j,name,nx,nx2,iens,ierr2,ierr3,pr_ens,             &
              maxens3,ensdim,                            &
              sig,APR_GR,APR_W,APR_MC,APR_ST,APR_AS,                 &
              APR_CAPMA,APR_CAPME,APR_CAPMI,closure_n,xland1,    &
              weight_GR,weight_W,weight_MC,weight_ST,weight_AS,training,  &
              ichoice,imid,ipr,jpr,itf,jtf,ktf, &
              its,ite, jts,jte, kts,kte, &
              dicycle,xf_dicycle )

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
        ichoice,imid,ipr,jpr,itf,jtf,ktf,     &
        its,ite, jts,jte, kts,kte
     integer, intent (in   )              ::                           &
        j,ensdim,nx,nx2,iens,maxens3,training
  ! xf_ens = ensemble mass fluxes
  ! pr_ens = precipitation ensembles
  ! dellat = change of temperature per unit mass flux of cloud ensemble
  ! dellaq = change of q per unit mass flux of cloud ensemble
  ! dellaqc = change of qc per unit mass flux of cloud ensemble
  ! outtem = output temp tendency (per s)
  ! outq   = output q tendency (per s)
  ! outqc  = output qc tendency (per s)
  ! pre    = output precip
  ! xmb    = total base mass flux
  ! xfac1  = correction factor
  ! pw = pw -epsilon*pd (ensemble dependent)
  ! ierr error value, maybe modified in this routine
  !
     real,    dimension (its:ite,jts:jte,1:ensdim)                     &
        ,intent (inout)                   ::                           &
       xf_ens,pr_ens
!srf ------
!    real,    dimension (its:ite,jts:jte)                              &
     real,    dimension (its:ite,jts:jte)                              &
         ,intent (inout)                   ::                          &
               APR_GR,APR_W,APR_MC,APR_ST,APR_AS,APR_CAPMA,            &
               APR_CAPME,APR_CAPMI
     real, dimension( its:ite , jts:jte )                      &
         ,intent(in) :: weight_gr,weight_w,weight_mc,weight_st,weight_as
!-srf---
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (out  )                   ::                           &
        outtem,outq,outqc,subt,subq,sub_mas
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in  )                   ::                           &
        zu
     real,   dimension (its:ite)                                      &
         ,intent (in  )                   ::                           &
        sig
     real,   dimension (its:ite,2)                                      &
         ,intent (in  )                   ::                           &
        xff_mid
     real,    dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
        pre,xmb
     real,    dimension (its:ite)                                      &
        ,intent (inout  )                   ::                           &
        closure_n,xland1
     real,    dimension (its:ite,kts:kte,1:nx)                     &
        ,intent (in   )                   ::                           &
       subt_ens,subq_ens,dellat,dellaqc,dellaq,pw
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        ktop
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr,ierr2,ierr3
!-srf begin
      logical, intent(IN) :: DICYCLE
      real,    intent(IN), dimension (its:ite) :: xf_dicycle
!-srf end
!
!  local variables in this routine
!

     integer                              ::                           &
        i,k,n,ncount
     real                                 ::                           &
        outtes,ddtes,dtt,dtq,dtqc,dtpw,prerate,clos_wei,xmbhelp
     real                                 ::                           &
        dtts,dtqs
     real,    dimension (its:ite)         ::                           &
       xfac1,xfac2
     real,    dimension (its:ite)::                           &
       xmb_ske,xmb_ave,xmb_std,xmb_cur,xmbweight
     real,    dimension (its:ite)::                           &
       pr_ske,pr_ave,pr_std,pr_cur
     real,    dimension (its:ite,jts:jte)::                           &
               pr_gr,pr_w,pr_mc,pr_st,pr_as,pr_capma,     &
               pr_capme,pr_capmi
     real, dimension (5) :: weight,wm,wm1,wm2,wm3
     real, dimension (its:ite,5) :: xmb_w

!
      character *(*), intent (in)        ::                           &
       name

!
     weight(1) = -999.  !this will turn off weights
     wm(1)=-999.

!
!
      DO k=kts,ktf
      do i=its,itf
        outtem(i,k)=0.
        outq(i,k)=0.
        outqc(i,k)=0.
        subt(i,k)=0.
        subq(i,k)=0.
        sub_mas(i,k)=0.
      enddo
      enddo
      do i=its,itf
        pre(i)=0.
        xmb(i)=0.
         xfac1(i)=0.
         xfac2(i)=0.
        xmbweight(i)=1.
      enddo
      do i=its,itf
        IF(ierr(i).eq.0)then
        do n=(iens-1)*nx*nx2*maxens3+1,iens*nx*nx2*maxens3
           if(pr_ens(i,j,n).le.0.)then
             if(i.eq.ipr.and.j.eq.jpr)write(0,*)'pr_ens',n,pr_ens(i,j,n),xf_ens(i,j,n)
             xf_ens(i,j,n)=0.
           endif
        enddo
        endif
      enddo
!
!--- calculate ensemble average mass fluxes
!

       xmb_w=0.
!
!-- now do feedback
!
      ddtes=100.
      if(imid.eq.0)then
      do i=its,itf
        if(ierr(i).eq.0)then
         k=0
         xmb_ave(i)=0.
         do n=(iens-1)*nx*nx2*maxens3+1,iens*nx*nx2*maxens3
          k=k+1
          xmb_ave(i)=xmb_ave(i)+xf_ens(i,j,n)
          if(i.eq.ipr.and.j.eq.jpr)write(0,*)'xmb_ave = ',k,xmb_ave(i),xf_ens(i,j,n)
         enddo
         xmb_ave(i)=xmb_ave(i)/float(k)
         !srf begin
         !print*,"XMB=",xmb_ave(i),xf_dicycle(i),xmb_ave(i) - xf_dicycle(i);call  flush(6)
         if(DICYCLE) then
            xmb_ave(i)=xmb_ave(i) - xf_dicycle(i)
         endif
         !srf end
         if(xmb_ave(i).le.0.)then
              ierr(i)=13
              xmb_ave(i)=0.
         endif
         xmb(i)=sig(i)*xmb_ave(i)
! --- Now use proper count of how many closures were actually
!       used in cup_forcing_ens (including screening of some
!       closures over water) to properly normalize xmb
           clos_wei=16./max(1.,closure_n(i))

           if(xmb(i).eq.0.)then
              ierr(i)=19
           endif
           if(xmb(i).gt.100.)then
              ierr(i)=19
           endif
           xfac1(i)=xmb(i)
           xfac2(i)=xmb(i)

        endif
      ENDDO
      else
         do i=its,itf
         IF(ierr(i).eq.0)then
           xmb(i)=.5*sig(i)*(xff_mid(i,1)+xff_mid(i,2))
           if(ichoice.gt.0)xmb(i)=sig(i)*xff_mid(i,ichoice)
         endif
         enddo
      endif

      DO i=its,itf
       IF(ierr(i).eq.0)then
         DO k=kts,ktop(i)
           dtt =0.
           dtts=0.
           dtq =0.
           dtqs=0.
           dtqc=0.
           dtpw=0.
           do n=1,nx
              dtt =dtt  + dellat  (i,k,n)
              dtts=dtts + subt_ens(i,k,n)
              dtq =dtq  + dellaq  (i,k,n)
              dtqs=dtqs + subq_ens(i,k,n)
              dtqc=dtqc + dellaqc (i,k,n)
              dtpw=dtpw + pw      (i,k,n)
           enddo
           OUTTEM(I,K)= XMB(I)* dtt /float(nx)
           SUBT  (I,K)= XMB(I)* dtts/float(nx)
           OUTQ  (I,K)= XMB(I)* dtq /float(nx)
           SUBQ  (I,K)= XMB(I)* dtqs/float(nx)
           OUTQC (I,K)= XMB(I)* dtqc/float(nx)
           PRE(I)=PRE(I)+XMB(I)*dtpw/float(nx)
           xf_ens(i,j,:)=sig(i)*xf_ens(i,j,:)*dtpw/float(nx)
           sub_mas(i,k)=zu(i,k)*xmb(i)
         ENDDO
       ENDIF
      ENDDO

      do i=its,itf
        if(ierr(i).eq.0)then
        do k=(iens-1)*nx*nx2*maxens3+1,iens*nx*nx2*maxens3
          xf_ens(i,j,k)=xf_ens(i,j,k)*xfac1(i)
        enddo
        endif
      ENDDO
      !do i=its,itf
      !  if(ierr(i).ne. 0)then
      !     !-check max/min values for outt
      !     if(minval(outtem(i,:)*86400.) < -100.) &
!             !print*,"==> low deep outt=",i, minval(outtem(i,:)*86400.) &
!             !;call flush(6)
!           !-check max/min values for outt
!           if(maxval(outtem(i,:)*86400.) > +100.) &
!             print*,"==> high deep outt=",i, maxval(outtem(i,:)*86400.)&
!             ;call flush(6)
!        endif
!      ENDDO
!srf-fix for preci
      do i=its,itf
        if(ierr(i).ne. 0)then
            apr_w (i,j)=0.0
            apr_st(i,j)=0.0
            apr_gr(i,j)=0.0
            apr_mc(i,j)=0.0
            apr_as(i,j)=0.0
        endif
      ENDDO
!srf
   END SUBROUTINE cup_output_ens_3d
!-------------------------------------------------------
   SUBROUTINE cup_up_moisture(name,ierr,z_cup,qc,qrc,pw,pwav,     &
              ccnclean,p_cup,kbcon,ktop,cd,dby,clw_all,&
              t_cup,q,GAMMA_cup,zu,qes_cup,k22,qe_cup,xl,         &
              ZQEXEC,use_excess,ccn,rho, &
              up_massentr,up_massdetr,psum,psumh,                 &
              autoconv,aeroevap,itest,itf,jtf,ktf,j,ipr,jpr,                &
              its,ite, jts,jte, kts,kte                     )

   IMPLICIT NONE
  real, parameter :: BDISPM = 0.366       !Berry--size dispersion (maritime)
  REAL, PARAMETER :: BDISPC = 0.146       !Berry--size dispersion (continental)
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
                                  use_excess,itest,autoconv,aeroevap,itf,jtf,ktf,           &
                                  its,ite, jts,jte,j,ipr,jpr, kts,kte
  ! cd= detrainment function
  ! q = environmental q on model levels
  ! qe_cup = environmental q on model cloud levels
  ! qes_cup = saturation q on model cloud levels
  ! dby = buoancy term
  ! cd= detrainment function
  ! zu = normalized updraft mass flux
  ! gamma_cup = gamma on model cloud levels
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        t_cup,p_cup,rho,q,zu,gamma_cup,qe_cup,                         &
        up_massentr,up_massdetr,dby,qes_cup,z_cup,cd
     real,    dimension (its:ite)                              &
        ,intent (in   )                   ::                           &
        zqexec
  ! entr= entrainment rate
     real                                                              &
        ,intent (in   )                   ::                           &
        ccnclean,xl
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        kbcon,ktop,k22
!
! input and output
!

   ! ierr error value, maybe modified in this routine

     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
      character *(*), intent (in)        ::                           &
       name
   ! qc = cloud q (including liquid water) after entrainment
   ! qrch = saturation q in cloud
   ! qrc = liquid water content in cloud after rainout
   ! pw = condensate that will fall out at that level
   ! pwav = totan normalized integrated condensate (I1)
   ! c0 = conversion rate (cloud to rain)

     real,    dimension (its:ite,kts:kte)                              &
        ,intent (out  )                   ::                           &
        qc,qrc,pw,clw_all
     real,    dimension (its:ite,kts:kte) ::                           &
        qch,qrcb,pwh,clw_allh
     real,    dimension (its:ite)         ::                           &
        pwavh
     real,    dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
        pwav,psum,psumh
     real,    dimension (its:ite)                                      &
        ,intent (in  )                   ::                           &
        ccn
!
!  local variables in this routine
!

     integer                              ::                           &
        iounit,iprop,iall,i,k,k1,k2
     real                                 ::                           &
        prop_ave,qrcb_h,bdsp,dp,g,rhoc,dh,qrch,c0,dz,radius,berryc0,q1,berryc
     real,    dimension (kts:kte)         ::                           &
        prop_b
!
        prop_b(kts:kte)=0
        iall=0
        c0=.002
        g=9.81
        bdsp=BDISPM
!
!--- no precip for small clouds
!
        if(name.eq.'shallow')then
            c0=0.
            iall=1
        endif
        do i=its,itf
          pwav(i)=0.
          pwavh(i)=0.
          psum(i)=0.
          psumh(i)=0.
        enddo
        do k=kts,ktf
        do i=its,itf
          pw(i,k)=0.
          pwh(i,k)=0.
          qc(i,k)=0.
          if(ierr(i).eq.0)qc(i,k)=qe_cup(i,k)
          if(ierr(i).eq.0)qch(i,k)=qe_cup(i,k)
          clw_all(i,k)=0.
          clw_allh(i,k)=0.
          qrc(i,k)=0.
          qrcb(i,k)=0.
        enddo
        enddo
!     if(use_excess < 2 ) then
      do i=its,itf
      if(ierr(I).eq.0)then

!-- mass con
!     do k=2,kbcon(i)-1  ! orignal
      do k=2,k22(i)      ! mass cons option
!-- mass con

        DZ=Z_cup(i,K)-Z_cup(i,K-1)
        qc (i,k)=qe_cup(i,k22(i))+float(use_excess)*zqexec(i)
        qch(i,k)=qe_cup(i,k22(i))+float(use_excess)*zqexec(i)
!
!       if(qc(i,k).gt.qes_cup(i,kbcon(i)-1))then
!           pw(i,k)=zu(i,k)*(qc(i,k)-qes_cup(i,kbcon(i)-1))
!           qc(i,k)=qes_cup(i,kbcon(i)-1)
!           qch(i,k)=qes_cup(i,kbcon(i)-1)
!           PWAV(I)=PWAV(I)+PW(I,K)
!           Psum(I)=Psum(I)+pw(I,K)*dz
!       endif
        if(i.eq.ipr.and.j.eq.jpr)write(0,*)'cup_mois ',k,pw(i,k)
      enddo
      endif
      enddo
!     else if(use_excess == 2) then
!        write(0,*)'in cup_mois'
!       do i=its,itf

!        if(ierr(I).eq.0)then
!            k1=max(1,k22(i)-1)
!            k2=k22(i)+1
!         do k=2,kbcon(i)-1
!            DZ=Z_cup(i,K)-Z_cup(i,K-1)
!            write(0,*)'cup_mois',k1,k2,zqexec(i),qc(i,k),sum(qe_cup(i,k1:k2))/float(k2-k1+1)
!            qc (i,k)=sum(qe_cup(i,k1:k2))/float(k2-k1+1) !+zqexec(i)
!            qch(i,k)=sum(qe_cup(i,k1:k2))/float(k2-k1+1) !+zqexec(i)
!            if(qc(i,k).gt.qes_cup(i,kbcon(i)-1))then
!                pw(i,k)=zu(i,k)*(qc(i,k)-qes_cup(i,kbcon(i)-1))
!                qc(i,k)=qes_cup(i,kbcon(i)-1)
!                qch(i,k)=qes_cup(i,kbcon(i)-1)
!                PWAV(I)=PWAV(I)+PW(I,K)
!                Psum(I)=Psum(I)+pw(I,K)*dz
!            endif

!         enddo !k
!        endif  !ierr
!       enddo !i
!     endif  ! use_excess

        DO k=kts+1,ktf
            DO i=its,itf
         IF(ierr(i).ne.0) cycle !GO TO 100

!-- mass con
!         IF(K.Lt.KBCON(I))GO TO 100     !  original
         IF (K.Lt.k22(I)+1) cycle!GO TO 100    !  mass cons option
!-- mass con


         IF(K.Gt.KTOP(I)) cycle!GO TO 100
         rhoc=.5*(rho(i,k)+rho(i,k-1))
         DZ=Z_cup(i,K)-Z_cup(i,K-1)
         DP=p_cup(i,K)-p_cup(i,K-1)
!
!--- saturation  in cloud, this is what is allowed to be in it
!
         QRCH=QES_cup(I,K)+(1./XL)*(GAMMA_cup(i,k) &
              /(1.+GAMMA_cup(i,k)))*DBY(I,K)
!
!------    1. steady state plume equation, for what could
!------       be in cloud without condensation
!
!
       qc(i,k)=   (qc(i,k-1)*zu(i,k-1)-.5*up_massdetr(i,k-1)* qc(i,k-1)+ &
                         up_massentr(i,k-1)*q(i,k-1))   /            &
                         (zu(i,k-1)-.5*up_massdetr(i,k-1)+up_massentr(i,k-1))
       qch(i,k)= (qch(i,k-1)*zu(i,k-1)-.5*up_massdetr(i,k-1)*qch(i,k-1)+ &
                         up_massentr(i,k-1)*q(i,k-1))   /            &
                         (zu(i,k-1)-.5*up_massdetr(i,k-1)+up_massentr(i,k-1))

        if(qc(i,k).le.qrch)then
          qc(i,k)=qrch
!         ierr(i)=232
!         go to 100
        endif
        if(qch(i,k).le.qrch)then
          qch(i,k)=qrch
!         ierr(i)=232
!         go to 100
        endif
!
!------- Total condensed water before rainout
!
        clw_all(i,k)=QC(I,K)-QRCH
        QRC(I,K)=(QC(I,K)-QRCH) ! /(1.+C0*DZ*zu(i,k))
        clw_allh(i,k)=QCH(I,K)-QRCH
        QRCB(I,K)=(QCH(I,K)-QRCH) ! /(1.+C0*DZ*zu(i,k))
    IF(autoconv.eq.2) then


!
! normalized berry
!
! first calculate for average conditions, used in cup_dd_edt!
! this will also determine proportionality constant prop_b, which, if applied,
! would give the same results as c0 under these conditions
!
         q1=1.e3*rhoc*qrcb(i,k)  ! g/m^3 ! g[h2o]/cm^3
         berryc0=q1*q1/(60.0*(5.0 + 0.0366*CCNclean/ &
                ( q1 * BDSP)  ) ) !/(
      !if(i.eq.ipr.and.j.eq.jpr)write(0,*)'cupm',k,rhoc,rho(i,k)
!        qrcb_h=qrcb(i,k)/(1.+c0*dz)
         qrcb_h=((QCH(I,K)-QRCH)*zu(i,k)-qrcb(i,k-1)*(.5*up_massdetr(i,k-1)))/ &
                   (zu(i,k)+.5*up_massdetr(i,k-1)+c0*dz*zu(i,k))
         prop_b(k)=c0*qrcb_h*zu(i,k)/(1.e-3*berryc0+1.e-8)
      !if(i.eq.ipr.and.j.eq.jpr)write(0,*)'cupm',berryc0,prop_b(k),qrcb_h
         pwh(i,k)=zu(i,k)*1.e-3*berryc0*dz*prop_b(k) ! 2.
         berryc=qrcb(i,k)
         qrcb(i,k)=((QCh(I,K)-QRCH)*zu(i,k)-pwh(i,k)-qrcb(i,k-1)*(.5*up_massdetr(i,k-1)))/ &
                   (zu(i,k)+.5*up_massdetr(i,k-1))
         if(qrcb(i,k).lt.0.)then
           berryc0=(qrcb(i,k-1)*(.5*up_massdetr(i,k-1))-(QCh(I,K)-QRCH)*zu(i,k))/zu(i,k)*1.e-3*dz*prop_b(k)
           pwh(i,k)=zu(i,k)*1.e-3*berryc0*dz*prop_b(k)
           qrcb(i,k)=0.
         endif
      !if(i.eq.ipr.and.j.eq.jpr)write(0,*)'cupm',zu(i,k),pwh(i,k),dz,qrch,qrcb(i,k),clw_allh(i,k)
      QCh(I,K)=QRCb(I,K)+qrch
      PWAVH(I)=PWAVH(I)+pwh(I,K)
      Psumh(I)=Psumh(I)+clw_allh(I,K)*zu(i,k) *dz
!
! then the real berry
!
          q1=1.e3*rhoc*qrc(i,k)  ! g/m^3 ! g[h2o]/cm^3
          berryc0=q1*q1/(60.0*(5.0 + 0.0366*CCN(i)/ &
                ( q1 * BDSP)  ) ) !/(
          berryc0=1.e-3*berryc0*dz*prop_b(k) ! 2.
          berryc=qrc(i,k)
         qrc(i,k)=((QC(I,K)-QRCH)*zu(i,k)-zu(i,k)*berryc0-qrc(i,k-1)*(.5*up_massdetr(i,k-1)))/ &
                   (zu(i,k)+.5*up_massdetr(i,k-1))
          if(qrc(i,k).lt.0.)then
            berryc0=((QC(I,K)-QRCH)*zu(i,k)-qrc(i,k-1)*(.5*up_massdetr(i,k-1)))/zu(i,k)
            qrc(i,k)=0.
          endif
          pw(i,k)=max(0.,berryc0*zu(i,k))
          QC(I,K)=QRC(I,K)+qrch
!
!  if not running with berry at all, do the following
!
       ELSE       !c0=.002
        if(iall.eq.1)then
          qrc(i,k)=0.
          pw(i,k)=(QC(I,K)-QRCH)*zu(i,k)
          if(pw(i,k).lt.0.)pw(i,k)=0.
        else

!---gg 10jul
         QRC(I,K)=(QC(I,K)-QRCH)/(1.+C0*DZ)      ! original form
!        QRC(I,K)=(QC(I,K)-QRCH)/(1.+(c1+C0)*DZ) ! new to be tested
!--
         PW(i,k)=c0*dz*QRC(I,K)*zu(i,k)
!        qrc(i,k)=((QC(I,K)-QRCH)*zu(i,k)-qrc(i,k-1)*(.5*up_massdetr(i,k-1)))/ &
!                  (zu(i,k)+.5*up_massdetr(i,k-1)+c0*dz*zu(i,k))
!        PW(i,k)=c0*dz*qrc(i,k)*zu(i,k)
         if(qrc(i,k).lt.0)then
           qrc(i,k)=0.
           pw(i,k)=0.
         endif
        endif
        QC(I,K)=QRC(I,K)+qrch
      endif !autoconv
!
!--- integrated normalized ondensate
!
         PWAV(I)=PWAV(I)+PW(I,K)
         Psum(I)=Psum(I)+clw_all(I,K)*zu(i,k) *dz
      end do
      end do
!   100     CONTINUE
       prop_ave=0.
       iprop=0
       do k=kts,kte
        prop_ave=prop_ave+prop_b(k)
        if(prop_b(k).gt.0)iprop=iprop+1
       enddo
       iprop=max(iprop,1)
!      write(11,*)'prop_ave = ',prop_ave/float(iprop)
!      print *,'pwav = ',pwav(1)

   END SUBROUTINE cup_up_moisture
!====================================================================
   SUBROUTINE gdinit(RTHCUTEN,RQVCUTEN,RQCCUTEN,RQICUTEN,           &
                        MASS_FLUX,cp,restart,                       &
                        P_QC,P_QI,P_FIRST_SCALAR,                   &
                        RTHFTEN, RQVFTEN,                           &
                        APR_GR,APR_W,APR_MC,APR_ST,APR_AS,          &
                        APR_CAPMA,APR_CAPME,APR_CAPMI,              &
                        allowed_to_read,                            &
                        ids, ide, jds, jde, kds, kde,               &
                        ims, ime, jms, jme, kms, kme,               &
                        its, ite, jts, jte, kts, kte               )
!--------------------------------------------------------------------
   IMPLICIT NONE
!--------------------------------------------------------------------
   LOGICAL , INTENT(IN)           ::  restart,allowed_to_read
   INTEGER , INTENT(IN)           ::  ids, ide, jds, jde, kds, kde, &
                                      ims, ime, jms, jme, kms, kme, &
                                      its, ite, jts, jte, kts, kte
   INTEGER , INTENT(IN)           ::  P_FIRST_SCALAR, P_QI, P_QC
   REAL,     INTENT(IN)           ::  cp

   REAL,     DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(OUT) ::       &
                                                          RTHCUTEN, &
                                                          RQVCUTEN, &
                                                          RQCCUTEN, &
                                                          RQICUTEN

   REAL,     DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(OUT) ::       &
                                                          RTHFTEN,  &
                                                          RQVFTEN

   REAL,     DIMENSION( ims:ime , jms:jme ) , INTENT(OUT) ::        &
                                APR_GR,APR_W,APR_MC,APR_ST,APR_AS,  &
                                APR_CAPMA,APR_CAPME,APR_CAPMI,      &
                                MASS_FLUX

   IF(.not.restart)THEN
        RTHCUTEN=0.
        RQVCUTEN=0.
        RTHFTEN=0.
        RQVFTEN=0.

     IF (P_QC .ge. P_FIRST_SCALAR) THEN
           RQCCUTEN=0.
     ENDIF

     IF (P_QI .ge. P_FIRST_SCALAR) THEN
           RQICUTEN=0.
     ENDIF

        mass_flux=0.

   ENDIF
        APR_GR=0.
        APR_ST=0.
        APR_W=0.
        APR_MC=0.
        APR_AS=0.
        APR_CAPMA=0.
        APR_CAPME=0.
        APR_CAPMI=0.

   END SUBROUTINE gdinit

!--------------------------------------------------------------------

      real function satvap(temp2)
      implicit none
      real :: temp2, temp, toot, toto, eilog, tsot,  &
     &        ewlog, ewlog2, ewlog3, ewlog4
      temp = temp2-273.155
      if (temp.lt.-20.) then   !!!! ice saturation
        toot = 273.16 / temp2
        toto = 1 / toot
        eilog = -9.09718 * (toot - 1) - 3.56654 * (log(toot) / &
     &    log(10.)) + .876793 * (1 - toto) + (log(6.1071) / log(10.))
        satvap = 10 ** eilog
      else
        tsot = 373.16 / temp2
        ewlog = -7.90298 * (tsot - 1) + 5.02808 * &
     &             (log(tsot) / log(10.))
        ewlog2 = ewlog - 1.3816e-07 * &
     &             (10 ** (11.344 * (1 - (1 / tsot))) - 1)
        ewlog3 = ewlog2 + .0081328 * &
     &             (10 ** (-3.49149 * (tsot - 1)) - 1)
        ewlog4 = ewlog3 + (log(1013.246) / log(10.))
        satvap = 10 ** ewlog4
      end if
      return
      end function
!--------------------------------------------------------------------
   SUBROUTINE CUP_gf_sh(J,COUPL_MPHYSICS,TCRIT  ,ccnclean,tscl_kf  &
             ,xl,rv,cp,g,ichoice,ipr,jpr            &
             ,autoconv,itf,jtf,ktf	     ,ktau  &
             ,use_excess,its,ite, jts,jte, kts,kte  &
	     !input data                         &
	     ,DTIME,dx                           &
	     ,entr_rate  			 &
	     ,zustart_sh,zufinal_sh,zutop_sh	 &
             ,h_sfc_flux,le_sfc_flux		 &
	     ,tsur				 &
	     ,PSUR				 &
             ,Z1,xland,gsw			 &
             ,ztexec				 &
	     ,zqexec				 &
	     ,ccn				 &
	     ,rho				 &
	     ,dhdt        			 &
             ,zws				 &
	     ,zo				 &
	     ,AAEQ				 &
	     ,T 				 &
	     ,Q 				 &
             ,TN				 &
	     ,QO				 &
	     ,PO				 &
	     ,P 				 &
	     ,US				 &
	     ,VS    				 &
	     !output data			 &
	     ,ierr,ierrc   			 &
             ,xmb_out				 &
	     ,OUTT				 &
	     ,OUTQ				 &
	     ,OUTQC				 &
             ,cupclw				 &
	     ,kpbl,kbcon,ktop,k22                &
!- for convective transport-start
             ,mgmxp,  mgmzp,mynum          &
             ,ierr4d	      &
             ,jmin4d	      &
             ,kdet4d	      &
             ,k224d	      &
             ,kbcon4d	      &
             ,ktop4d	      &
             ,kpbl4d	      &
             ,kstabi4d        &
             ,kstabm4d        &
             ,xmb4d	      &
             ,edt4d	      &
             ,pwav4d	      &
             ,pcup5d	      &
             ,up_massentr5d   &
             ,up_massdetr5d   &
             ,dd_massentr5d   &
             ,dd_massdetr5d   &
             ,zup5d	      &
             ,zdn5d	      &
             ,prup5d	      &
             ,prdn5d	      &
             ,clwup5d	      &
             ,tup5d	      )
!- for convective transport-end

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
        autoconv,itf,jtf,ktf,ktau,use_excess,        &
        its,ite, jts,jte, kts,kte,ipr,jpr
     integer, intent (in   )              ::                           &
        j,ichoice
  !
  !
  !
     real,    dimension (its:ite,jts:jte)                              &
        ,intent (in   )                   ::                           &
               gsw
  ! outtem = output temp tendency (per s)
  ! outq   = output q tendency (per s)
  ! outqc  = output qc tendency (per s)
  ! pre    = output precip
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (inout  )                   ::                           &
        OUTT,OUTQ,OUTQC,cupclw
     real,    dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
        xmb_out
     integer,    dimension (its:ite)                                   &
        ,intent (out  )                   ::                           &
        kbcon,ktop,k22
     integer,    dimension (its:ite)                                   &
        ,intent (in  )                   ::                           &
        kpbl
  !
  ! basic environmental input includes moisture convergence (mconv)
  ! omega (omeg), windspeed (us,vs), and a flag (aaeq) to turn off
  ! convection for this call only and at that particular gridpoint
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        rho,T,PO,P,US,VS,tn,dhdt
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (inout)                   ::                           &
         Q,QO
     real, dimension (its:ite)                                         &
        ,intent (in   )                   ::                           &
        zws,ztexec,zqexec,ccn,Z1,PSUR,AAEQ,xland,h_sfc_flux,le_sfc_flux&
        ,tsur

       real                                                            &
        ,intent (in   )                   ::                           &
        tscl_kf,dx,ccnclean,dtime,tcrit,xl,cp,rv,g,zustart_sh,zufinal_sh,&
        zutop_sh,entr_rate


!
!
!***************** the following are your basic environmental
!                  variables. They carry a "_cup" if they are
!                  on model cloud levels (staggered). They carry
!                  an "o"-ending (z becomes zo), if they are the forced
!                  variables. They are preceded by x (z becomes xz)
!                  to indicate modification by some typ of cloud
!
  ! z           = heights of model levels
  ! q           = environmental mixing ratio
  ! qes         = environmental saturation mixing ratio
  ! t           = environmental temp
  ! p           = environmental pressure
  ! he          = environmental moist static energy
  ! hes         = environmental saturation moist static energy
  ! z_cup       = heights of model cloud levels
  ! q_cup       = environmental q on model cloud levels
  ! qes_cup     = saturation q on model cloud levels
  ! t_cup       = temperature (Kelvin) on model cloud levels
  ! p_cup       = environmental pressure
  ! he_cup = moist static energy on model cloud levels
  ! hes_cup = saturation moist static energy on model cloud levels
  ! gamma_cup = gamma on model cloud levels
!
!
  ! hcd = moist static energy in downdraft
  ! zd normalized downdraft mass flux
  ! dby = buoancy term
  ! entr = entrainment rate
  ! zd   = downdraft normalized mass flux
  ! entr= entrainment rate
  ! hcd = h in model cloud
  ! bu = buoancy term
  ! zd = normalized downdraft mass flux
  ! gamma_cup = gamma on model cloud levels
  ! qcd = cloud q (including liquid water) after entrainment
  ! qrch = saturation q in cloud
  ! pwd = evaporate at that level
  ! pwev = total normalized integrated evaoprate (I2)
  ! entr= entrainment rate
  ! z1 = terrain elevation
  ! entr = downdraft entrainment rate
  ! jmin = downdraft originating level
  ! kdet = level above ground where downdraft start detraining
  ! psur        = surface pressure
  ! z1          = terrain elevation
  ! pr_ens = precipitation ensemble
  ! xf_ens = mass flux ensembles
  ! massfln = downdraft mass flux ensembles used in next timestep
  ! omeg = omega from large scale model
  ! mconv = moisture convergence from large scale model
  ! zd      = downdraft normalized mass flux
  ! zu      = updraft normalized mass flux
  ! dir     = "storm motion"
  ! mbdt    = arbitrary numerical parameter
  ! dtime   = dt over which forcing is applied
  ! iact_gr_old = flag to tell where convection was active
  ! kbcon       = LFC of parcel from k22
  ! k22         = updraft originating level
  ! icoic       = flag if only want one closure (usually set to zero!)
  ! dby = buoancy term
  ! ktop = cloud top (output)
  ! xmb    = total base mass flux
  ! hc = cloud moist static energy
  ! hkb = moist static energy at originating level

     real,    dimension (its:ite,kts:kte) ::                           &
        entr_rate_2d,mentrd_rate_2d,he,hes,qes,z,                      &
        heo,heso,qeso,zo,                                              &
        xhe,xhes,xqes,xz,xt,xq,                                        &

        qes_cup,q_cup,he_cup,hes_cup,z_cup,p_cup,gamma_cup,t_cup,      &
        qeso_cup,qo_cup,heo_cup,heso_cup,zo_cup,po_cup,gammao_cup,     &
        tn_cup,                                                        &
        xqes_cup,xq_cup,xhe_cup,xhes_cup,xz_cup,xp_cup,xgamma_cup,     &
        xt_cup,                                                        &

        xlamue,dby,qc,qrcd,pwd,pw,hcd,qcd,dbyd,hc,qrc,zu,zd,clw_all,   &
        dbyo,qco,qrcdo,pwdo,pwo,hcdo,qcdo,dbydo,hco,qrco,zuo,zdo,      &
        xdby,xqc,xqrcd,xpwd,xpw,xhcd,xqcd,xhc,xqrc,xzu,xzd, tempco,    &

  ! cd  = detrainment function for updraft
  ! cdd = detrainment function for downdraft
  ! dellat = change of temperature per unit mass flux of cloud ensemble
  ! dellaq = change of q per unit mass flux of cloud ensemble
  ! dellaqc = change of qc per unit mass flux of cloud ensemble

        cd,cdd,DELLAH,DELLAQ,DELLAT,DELLAQC,dsubt,dsubh,dsubq,subt,subq

  ! aa0 cloud work function for downdraft
  ! edt = epsilon
  ! aa0     = cloud work function without forcing effects
  ! aa1     = cloud work function with forcing effects
  ! xaa0    = cloud work function with cloud effects (ensemble dependent)
  ! edt     = epsilon

     real,    dimension (its:ite) ::                                   &
       edt,edto,edtx,AA1,AA0,XAA0,HKB,                          &
       HKBO,XHKB,QKB,QKBO,                                    &
       xmbmax,XMB,XPWAV,XPWEV,PWAV,PWEV,PWAVO,                                &
       PWEVO,BU,BUD,BUO,cap_max,xland1,                                    &
       cap_max_increment,closure_n,psum,psumh,sig,zuhe,cape
     integer,    dimension (its:ite) ::                                &
       kzdown,KDET,KB,JMIN,kstabi,kstabm,K22x,        &   !-lxz
       KBCONx,KBx,KTOPx,ierr,ierr2,ierr3,KBMAX

     integer                              ::                           &
       nall,iedt,nens,nens3,ki,I,K,KK,iresult
     real                                 ::                           &
      day,dz,dzo,mbdt,radius,entrd_rate,mentrd_rate,  &
      zcutdown,edtmax,edtmin,depth_min,zkbmax,z_detr,zktop,      &
      massfld,dh,cap_maxs,trash,frh,xlamdd,fsum

      real detdo1,detdo2,entdo,dp,subin,detdo,entup,                &
      detup,subdown,entdoj,entupk,detupk,totmas
      real :: zufinals,power_entr,zustart,zufinal,dzm1,dzp1


     integer :: icount,tun_lim,jprnt,k1,k2,kbegzu,kfinalzu,kstart,jmini,levadj,nvar
     logical :: keep_going
     real xff_shal(9),blqe,xkshal
     character*50 :: ierrc(its:ite)
     real,    dimension (its:ite,kts:kte) ::                           &
       up_massentr,up_massdetr,dd_massentr,dd_massdetr                 &
      ,up_massentro,up_massdetro,dd_massentro,dd_massdetro
     real,    dimension (kts:kte) :: smth
     real :: C_up, tcold,thot,efic,fin,trash2
     logical, intent(IN) :: COUPL_MPHYSICS

!-srf- for convective transport-start
     integer, intent (in   )              ::    mgmxp,  mgmzp ,mynum
     integer, intent (inout)  ::     &
               ierr4d   (mgmxp)      &
              ,jmin4d          (mgmxp)             &
              ,kdet4d          (mgmxp)             &
              ,k224d           (mgmxp)             &
              ,kbcon4d         (mgmxp)             &
              ,ktop4d          (mgmxp)             &
              ,kpbl4d          (mgmxp)             &
              ,kstabi4d        (mgmxp)             &
              ,kstabm4d        (mgmxp)

     real   , intent (inout)  ::     &
               xmb4d        (mgmxp)             &
              ,edt4d        (mgmxp)             &
              ,pwav4d        (mgmxp)

     real   , intent (inout)  ::     &
               pcup5d               (mgmzp,mgmxp)      &
              ,up_massentr5d (mgmzp,mgmxp)      &
              ,up_massdetr5d (mgmzp,mgmxp)      &
              ,dd_massentr5d (mgmzp,mgmxp)      &
              ,dd_massdetr5d (mgmzp,mgmxp)      &
              ,zup5d             (mgmzp,mgmxp)      &
              ,zdn5d                (mgmzp,mgmxp)      &
              ,prup5d               (mgmzp,mgmxp)      &
              ,prdn5d               (mgmzp,mgmxp)      &
              ,clwup5d              (mgmzp,mgmxp)      &
              ,tup5d                (mgmzp,mgmxp)
!-srf- for convective transport-end

      zustart=zustart_sh
      zufinal=zufinal_sh


      levadj=4
      power_entr=2.
      icount=0
      day=86400.
      do i=its,itf
        xmb_out(i)=0.
        xland1(i)=xland(i) ! 1.
        if(xland(i).gt.1.5)xland1(i)=0.
        cap_max_increment(i)=25.
        ierrc(i)=" "
      enddo
!
!--- initial entrainment rate
!--- initial detrainment rates
!
      do k=kts,ktf
      do i=its,itf
        up_massentro(i,k)=0.
        up_massdetro(i,k)=0.
        z(i,k)=zo(i,k)
        xz(i,k)=zo(i,k)
        qrco(i,k)=0.
        cd(i,k)=1.*entr_rate
        dellaqc(i,k)=0.
        cupclw(i,k)=0.
      enddo
      enddo
!
!--- max/min allowed value for epsilon (ratio downdraft base mass flux/updraft
!
!--- minimum depth (m), clouds must have
!
      depth_min=50.
!
!--- maximum depth (mb) of capping
!--- inversion (larger cap = no convection)
!
      cap_maxs=75.
      DO i=its,itf
        kbmax(i)=1
        aa0(i)=0.
        aa1(i)=0.
      enddo
      do i=its,itf
        cap_max(i)=cap_maxs
        iresult=0
      enddo
!
!--- max height(m) above ground where updraft air can originate
!
      zkbmax=3000.
!
!--- calculate moist static energy, heights, qes
!
      call cup_env(z,qes,he,hes,t,q,p,z1, &
           psur,ierr,tcrit,-1,xl,cp,   &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
      call cup_env(zo,qeso,heo,heso,tn,qo,po,z1, &
           psur,ierr,tcrit,-1,xl,cp,   &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)

!
!--- environmental values on cloud levels
!
      call cup_env_clev(t,qes,q,he,hes,z,p,qes_cup,q_cup,he_cup, &
           hes_cup,z_cup,p_cup,gamma_cup,t_cup,psur, &
           ierr,z1,xl,rv,cp,          &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
      call cup_env_clev(tn,qeso,qo,heo,heso,zo,po,qeso_cup,qo_cup, &
           heo_cup,heso_cup,zo_cup,po_cup,gammao_cup,tn_cup,psur,  &
           ierr,z1,xl,rv,cp,          &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
      do i=its,itf
        if(ierr(i).eq.0)then
        if(aaeq(i).lt.-0.1)then
           ierr(i)=20
        endif
!
      do k=kts,ktf
        if(zo_cup(i,k).gt.zkbmax+z1(i))then
          kbmax(i)=k
          go to 25
        endif
      enddo
 25   continue
!
      kbmax(i)=min(kbmax(i),ktf-4)
      endif
      enddo

!
!
!
!------- DETERMINE LEVEL WITH HIGHEST MOIST STATIC ENERGY CONTENT - K22
!
      CALL cup_MAXIMI(HEO_CUP,1,KBMAX,K22,ierr, &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
       DO 36 i=its,itf
!srf     k22(i)=kpbl(i)
         k22(i)=max(2,kpbl(i)-2)
!         k22(i)=2
         if(kpbl(i).gt.3)cap_max(i)=po_cup(i,kpbl(i))
         IF(ierr(I).eq.0)THEN
         IF(K22(I).GT.KBMAX(i))then
           ierr(i)=2
           ierrc(i)="could not find k22"
         endif
!           if(kpbl(i).gt.3)then
!              k22(i)=kpbl(i)
!              ierr(i)=0
!              ierrc(i)=" to zero becausof kpbl"
!            endif
!        else
!            ierrc(i)="why here? "
         endif
       if(j.eq.jpr .and. i.eq.ipr)write(0,*)'initial k22 = ',k22(ipr),kpbl(i)
 36   CONTINUE
!
!--- DETERMINE THE LEVEL OF CONVECTIVE CLOUD BASE  - KBCON
!

      do i=its,itf
       IF(ierr(I).eq.0)THEN
        ! if(use_excess == 2) then
        !     k1=max(1,k22(i)-1)
        !     k2=k22(i)+1
        !     hkb(i) =sum(he_cup(i,k1:k2))/float(k2-k1+1)
        !     hkbo(i)=sum(heo_cup(i,k1:k2))/float(k2-k1+1)
        !     qkbo(i)=sum(qo_cup(i,k1:k2))/float(k2-k1+1)
!       !     write(0,*)sum(heo_cup(i,k1:k2))/float(k2-k1+1),heo_cup(i,k1),heo(i,k1:k2)
        !else if(use_excess <= 1) then
        !     hkb(i)=he_cup(i,k22(i))
        !     hkbo(i)=heo_cup(i,k22(i))
        !     qkbo(i)=qo_cup(i,k22(i))
        !endif  ! excess
        if(use_excess == 2) then
             k1=max(1,k22(i)-1)
             k2=k22(i)
             hkb(i) =sum(he_cup (i,k1:k2))/float(k2-k1+1)+xl*zqexec(i)+cp*ztexec(i)
             hkbo(i)=sum(heo_cup(i,k1:k2))/float(k2-k1+1)+xl*zqexec(i)+cp*ztexec(i)
             qkbo(i)=sum(qo_cup (i,k1:k2))/float(k2-k1+1)+ zqexec(i)
!            write(0,*)sum(heo_cup(i,k1:k2))/float(k2-k1+1),heo_cup(i,k1),heo(i,k1:k2)
        else if(use_excess <= 1) then
!<<<<<<<<<<<<<<<<<<<<<<<
!             k22(i)=min(1,kpbl(i))!<<<<<<<<<<<<<<<<<<<<<<<
!<<<<<<<<<<<<<<<<<<<<<<<
             hkb (i)=he_cup (i,1)+float(use_excess)*(xl*zqexec(i)+cp*ztexec(i))
             hkbo(i)=heo_cup(i,1)+float(use_excess)*(xl*zqexec(i)+cp*ztexec(i))
             qkbo(i)=qo_cup (i,1)+float(use_excess)*zqexec(i)
!             hkb (i)=he_cup (i,k22(i))+float(use_excess)*(xl*zqexec(i)+cp*ztexec(i))
!             hkbo(i)=heo_cup(i,k22(i))+float(use_excess)*(xl*zqexec(i)+cp*ztexec(i))
!             qkbo(i)=qo_cup (i,k22(i))+float(use_excess)*zqexec(i)
        endif  ! excess
         !do k=1,k22(i)
         !   hkb(i) =max(hkb (i),he_cup(i,k))
         !   hkbo(i)=max(hkbo(i),heo_cup(i,k))
         !   qkbo(i)=max(qkbo(i),qo_cup(i,k))
         !enddo
!         if(use_excess >= 0) then
!             hkbo(i)=hkbo(i) ! +float(use_excess)*(xl*zqexec(i)+cp*ztexec(i))
! do not cause artifical cooling
!            do k=1,k22(i)
!              he_cup(i,k)=min(hes_cup(i,k),hkb(i))
!              heo_cup(i,k)=min(heso_cup(i,k),hkbo(i))
!              qo_cup(i,k)=min(qeso_cup(i,k),qkbo(i))
! do not go past saturation
!              if(hkb(i).gt.hes_cup(i,k))hkb(i)=hes_cup(i,k)
!              if(hkbo(i).gt.heso_cup(i,k))hkbo(i)=heso_cup(i,k)
!              if(qkbo(i).gt.qes_cup(i,k))qkbo(i)=qeso_cup(i,k)
!            enddo
!        endif  ! excess
       endif ! ierr
      enddo

      do i=its,itf
      do k=kts,ktf
          dbyo(i,k)=hkbo(i)-heso_cup(i,k)
      enddo
      enddo
      do i=its,itf
      do k=kts,ktf
         if(dbyo(i,k).lt.0)then
            ktop(i)=k-1
            go to 441
         endif
      enddo
 441       continue
      enddo


      call cup_kbcon(ierrc,cap_max_increment,5,k22,kbcon,heo_cup,heso_cup, &
           hkbo,ierr,kbmax,po_cup,cap_max, &
           xl,cp,ztexec,zqexec,use_excess, &
           0,itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte, &
           z_cup,entr_rate,heo)



!
!--- increase detrainment in stable layers
!
      DO i=its,itf
         IF(ierr(I).eq.0)THEN
         !print*,"kbcon/k22/kpbl=",kbcon(i),k22(i),kpbl(i);call flush(6)
            if(kbcon(i).gt.ktf-4)then
                ierr(i)=231
            endif
            do k=kts,ktf
               frh = min(qo_cup(i,k)/qeso_cup(i,k),1.)
               entr_rate_2d(i,k)=entr_rate*(1.3-frh)
               cd(i,k)=entr_rate_2d(i,k)
            enddo
         endif
      enddo

      !call rates_up_curvefitting('shallow',ktop,ierr,po_cup,entr_rate_2d,hkbo,heo,heso_cup,zo_cup, &
      !0,kbcon,k22,kbcon,its,ite,itf,kts,kte,ktf,zuo,zustart_sh,zufinal_sh,zutop_sh)
       call rates_up_pdf('shallow',ktop,ierr,po_cup,entr_rate_2d,hkbo,heo,heso_cup,zo_cup, &
                    kstabi,k22,kbcon,its,ite,itf,kts,kte,ktf,zuo,kpbl)
      do i=its,itf
         if(ierr(i).eq.0)then
         do k=ktop(i)-1,1,-1
             if(zuo(i,k).lt.1.e-6)then
     !          k22(i)=k !<<<<<<<<<<<<<<<<<<<
               exit
             endif
         enddo

         do k=ktop(i)+1,ktf
           zuo(i,k)=0.
         enddo
      !   k22(i)=max(2,k22(i))!<<<<<<<<<<<<<<<<<<<<
         endif
      enddo
!
! calculate mass entrainment and detrainment
!
      do k=kts,ktf
      do i=its,itf
         hc  (i,k)=0.
         qco (i,k)=0.
         qrco(i,k)=0.
         DBY (I,K)=0.
         hco (i,k)=0.
         DBYo(I,K)=0.
      enddo
      enddo
      do i=its,itf
       IF(ierr(I).eq.0)THEN
         do k=1,kbcon(i)-1
            hc(i,k)=hkb(i)
            hco(i,k)=hkbo(i)
            qco(i,k)=qkbo(i)
            if(qco(i,k).gt.qeso_cup(i,k))then
               qrco  (i,k)= qco(i,k)-qeso_cup(i,k)
               qco   (i,k)= qeso_cup(i,k)
               cupclw(i,k)= qrco(i,k)
            endif
         enddo
         k=kbcon(i)
         hc(i,k)=hkb(i)
         qco(i,k)=qkbo(i)
         DBY(I,Kbcon(i))=Hkb(I)-HES_cup(I,K)
         hco(i,k)=hkbo(i)
         DBYo(I,Kbcon(i))=Hkbo(I)-HESo_cup(I,K)
         !- in-cloud saturation value
         trash=QESo_cup(I,K)+(1./XL)*(GAMMAo_cup(i,k) &
                             /(1.+GAMMAo_cup(i,k)))*DBYo(I,K)

         !- do not allow qco be sub-saturated at cloud base.
         qco(i,k)=max(trash,qco(i,k))

         if(qco(i,k)>=trash) then
           qrco(i,k)= qco(i,k)-trash
           qco (i,k)= trash
         else
           qrco(i,k)= 0.0
         endif
         cupclw(i,k)=qrco(i,k)
       endif ! ierr
      enddo
!
!
      do 42 i=its,itf
         if(ierr(i).eq.0)then
         !srf zu (i,1)=0.
         zuo(i,1)=0.
         xzu(i,1)= zuo(i,1)
         zu (i,1)= zuo(i,1)

         !- mass entrainment and detrinament is defined on model levels

         do k=2,maxloc(zuo(i,:),1)
         !=> below maximum value zu -> change entrainment
           dz=zo_cup(i,k)-zo_cup(i,k-1)

           up_massdetro(i,k-1)=cd(i,k-1)*dz*zuo(i,k-1)
           up_massentro(i,k-1)=zuo(i,k)-zuo(i,k-1)+up_massdetro(i,k-1)
           if(zuo(i,k-1).gt.0.) &
              entr_rate_2d(i,k-1)=(up_massentro(i,k-1))/(dz*zuo(i,k-1))
         enddo
         do k=maxloc(zuo(i,:),1)+1,ktf-1
         !=> above maximum value zu -> change detrainment
           dz=zo_cup(i,k)-zo_cup(i,k-1)
           up_massentro(i,k-1)=entr_rate_2d(i,k-1)*dz*zuo(i,k-1)

           !-special treatment for ktop
           if(k-1==ktop(i))up_massentro(i,k-1)=0.0

           up_massdetro(i,k-1)=zuo(i,k-1)+up_massentro(i,k-1)-zuo(i,k)

           if(zuo(i,k-1).gt.0.) &
             cd(i,k-1)=up_massdetro(i,k-1)/(dz*zuo(i,k-1))
         enddo


         do k=2,ktf-1
          xzu(i,k)= zuo(i,k)
          zu (i,k)= zuo(i,k)
          up_massentr(i,k-1)=up_massentro(i,k-1)
          up_massdetr(i,k-1)=up_massdetro(i,k-1)
         enddo
         do k=kbcon(i)+1,ktop(i)
          hc(i,k)=(hc(i,k-1)*zu(i,k-1)-.5*up_massdetr(i,k-1)*hc(i,k-1)+ &
                         up_massentr(i,k-1)*he(i,k-1))   /            &
                         (1.e-16+zu(i,k-1)-.5*up_massdetr(i,k-1)+up_massentr(i,k-1))
          dby(i,k)=hc(i,k)-hes_cup(i,k)

          hco(i,k)=(hco(i,k-1)*zuo(i,k-1)-.5*up_massdetro(i,k-1)*hco(i,k-1)+ &
                         up_massentro(i,k-1)*heo(i,k-1))   /            &
                         (1.e-16+zuo(i,k-1)-.5*up_massdetro(i,k-1)+up_massentro(i,k-1))

          dbyo(i,k)=hco(i,k)-heso_cup(i,k)

         enddo

         do k=kbcon(i)+1,ktop(i)
          !- in-cloud saturation value
          trash =QESo_cup(I,K)+(1./XL)*(GAMMAo_cup(i,k) &
                             /(1.+GAMMAo_cup(i,k)))*DBYo(I,K)

          !- total water liq+vapour
          trash2  = qco(i,k-1)+qrco(i,k-1)
          qco (i,k)=   (trash2* ( zuo(i,k-1)-0.5*up_massdetr(i,k-1)) + &
                       up_massentr(i,k-1)*qo(i,k-1))   /               &
                       (1.e-16+zuo(i,k-1)-.5*up_massdetr(i,k-1)+up_massentr(i,k-1))

!          qco (i,k)=   (qco(i,k-1)*zuo(i,k-1)-.5*up_massdetr(i,k-1)* qco(i,k-1)+ &
!                       up_massentr(i,k-1)*qo(i,k-1))   /            &
!                       (zuo(i,k-1)-.5*up_massdetr(i,k-1)+up_massentr(i,k-1))
!
!          qrco(i,k)=   (qrco(i,k-1)*zuo(i,k-1)-.5*up_massdetr(i,k-1)* qrco(i,k-1)) /            &
!                       (zuo(i,k-1)-.5*up_massdetr(i,k-1)+up_massentr(i,k-1))
!

          if(qco(i,k)>=trash) then
              ! cloud liquid water
              qrco(i,k)= qco(i,k)-trash
              ! cloud water vapor
              qco (i,k)= trash

          else
              qrco(i,k)= 0.0
          endif
          cupclw(i,k)=qrco(i,k)
           !print*,"tas=",ktop(i),k,trash,DBYo(I,K),qco (i,k),qrco(i,k)
         enddo

         !-
         do k=ktop(i)+1,ktf-1
           HC(i,K)=hes_cup(i,k)
           HCo(i,K)=heso_cup(i,k)
           qco (i,k)=QESo_cup(I,K)
           qrco(i,k)=0.0
           DBY(I,K)=0.
           DBYo(I,K)=0.
           zu(i,k)=0.
           xzu(i,k)=0.
           zuo(i,k)=0.
           cd(i,k)=0.
           entr_rate_2d(i,k)=0.
           up_massentr(i,k)=0.
           up_massdetr(i,k)=0.
           up_massentro(i,k)=0.
           up_massdetro(i,k)=0.
         enddo
         if(i.eq.ipr.and.j.eq.jpr)then
            write(0,*)'hcnew = '
            do k=1,ktf
              write(0,*)k,hco(i,k),dbyo(i,k)
            enddo
         endif
      endif
42    continue
!
!--- calculate workfunctions for updrafts
!
      call cup_up_aa0(aa0,z,zu,dby,GAMMA_CUP,t_cup, &
           kbcon,ktop,ierr,           &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
      call cup_up_aa0(aa1,zo,zuo,dbyo,GAMMAo_CUP,tn_cup, &
           kbcon,ktop,ierr,           &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
      do i=its,itf
         if(ierr(i).eq.0)then
           if(aa1(i).eq.0.)then
               ierr(i)=17
               ierrc(i)="cloud work function zero"
           endif
         endif
      enddo
!
!--- calculate in-cloud air temperature for CAPE
!
      do i=its,itf
         tempco(i,:)=t_cup(i,:)
         if(ierr(i)== 0)then
            !print*,"tempco",ktop(i)
           do k=1,ktop(i)+1
            tempco(i,k) = (1./cp)*(hco(i,k)-g*zo_cup(i,k)-xl*qco(i,k))
            !print*,"tempco",k,tempco(i,k),t_cup(i,k)
           enddo
          endif
      enddo
      call cup_up_cape(cape,z,zu,dby,GAMMA_CUP,t_cup, &
           kbcon,ktop,ierr,tempco,qco,qrco, qo_cup  , &
           itf,jtf,ktf,its,ite, jts,jte, kts,kte)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
!--- change per unit mass that a model cloud would modify the environment
!
!--- 1. in bottom layer
!
      do k=kts,ktf
      do i=its,itf
        dellah(i,k)=0.
        dsubt(i,k)=0.
        dsubh(i,k)=0.
        dellaq(i,k)=0.
        dsubq(i,k)=0.
      enddo
      enddo
!
!----------------------------------------------  cloud level ktop
!
!- - - - - - - - - - - - - - - - - - - - - - - - model level ktop-1
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!
!----------------------------------------------  cloud level k+2
!
!- - - - - - - - - - - - - - - - - - - - - - - - model level k+1
!
!----------------------------------------------  cloud level k+1
!
!- - - - - - - - - - - - - - - - - - - - - - - - model level k
!
!----------------------------------------------  cloud level k
!
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!
!----------------------------------------------  cloud level 3
!
!- - - - - - - - - - - - - - - - - - - - - - - - model level 2
!
!----------------------------------------------  cloud level 2
!
!- - - - - - - - - - - - - - - - - - - - - - - - model level 1
!
      do i=its,itf
        trash=0.
        if(ierr(i).eq.0)then
         do k=kts,ktop(i)+1

           ! entrainment/detrainment for updraft
            entup=up_massentro(i,k)
            detup=up_massdetro(i,k)

            totmas=detup-entup+zuo(i,k+1)-zuo(i,k)
            if(abs(totmas).gt.1.e-6)then
               write(0,*)'*********************',i,j,k,totmas
               write(0,*)k22(i),kbcon(i),ktop(i)
               write(0,123)k,detup,entup,detupk,zuo(i,k+1),zuo(i,k)
               123     format(1X,i2,5E12.4)
            !        call wrf_error_fatal ( 'totmas .gt.1.e-6' )
            endif
            dp=100.*(po_cup(i,k)-po_cup(i,k+1))

            dellah(i,k) =-(zuo(i,k+1)*(hco(i,k+1)-heo_cup(i,k+1) )-     &
                           zuo(i,k  )*(hco(i,k  )-heo_cup(i,k  ) ))*g/dp

            !-- take out cloud liquid water for detrainment
            dellaqc(i,k)=   detup*0.5*(qrco(i,k+1)+qrco(i,k)) *g/dp

            !-- condensation source term = detrained + flux divergence of
            !-- cloud liquid water (qrco)
            C_up = dellaqc(i,k)+(zuo(i,k+1)* qrco(i,k+1) -       &
                                 zuo(i,k  )* qrco(i,k  )  )*g/dp

            !-- water vapor budget (flux divergence of Q_up-Q_env - condensation term)
            dellaq(i,k) =-(zuo(i,k+1)*(qco(i,k+1)-qo_cup(i,k+1) ) -         &
                           zuo(i,k  )*(qco(i,k  )-qo_cup(i,k  ) ) )*g/dp &
                           - C_up

            !- check water conservation liq+condensed
            trash=trash+ (dellaq(i,k)+dellaqc(i,k))*dp/g
         enddo   ! k
         !
         if(abs(trash)>1.e-6) then
           write(6,*)'=> not water mass cons for shallow= ',i,trash
         endif

       endif
      enddo
!

!--- using dellas, calculate changed environmental profiles
!
       mbdt=3.e-4
      !do k=kts,ktf
      !  write(0,*)'zuo,k22 = ',k,zuo(1,k),k22(1)
      !enddo
      do k=kts,ktf
      do i=its,itf
         dellat(i,k)=0.
         if(ierr(i).eq.0)then
            dsubh(i,k)=0.0
            dsubq(i,k)=0.0
            dsubt(i,k)=0.0
            !adding dellaqc to dellaq:
            !dellaq (i,k)= dellaq(i,k)+dellaqc(i,k)
            !dellaqc(i,k)=0.0
            !
            XHE(I,K)=(dsubh(i,k)+DELLAH(I,K))*MBDT+HEO(I,K)
            XQ(I,K) =(dsubq(i,k)+DELLAQ(I,K)+DELLAQC(i,k))*MBDT+QO(I,K)

            !- don�t feed dellat with dellaqc if
            !- the detrainment of liquid water will be used as
            !- a source for cloud microphysics (then microphysics will
            !- evaporate this excess and provide the cooling at the
            !- detrainment region)
            if(COUPL_MPHYSICS) then
	      DELLAT(I,K)=(1./cp)*(DELLAH(I,K)-xl*(DELLAQ(I,K)))
            else
	      DELLAT(I,K)=(1./cp)*(DELLAH(I,K)-xl*(DELLAQ(I,K)+dellaqc(i,k)))
            endif
            DSUBT(I,K)=(1./cp)*(dsubt(i,k)-xl*dsubq(i,k))
            XT(I,K)= (DELLAT(I,K)+DSUBT(I,K)-XL/CP*DELLAQC(I,K))*MBDT+TN(I,K)
            IF(XQ(I,K).LE.0.)XQ(I,K)=1.E-08
         ENDIF
      enddo
      enddo
      do i=its,itf
      if(ierr(i).eq.0)then
      xhkb(i)=hkbo(i)+(dsubh(i,k22(i))+DELLAH(I,K22(i)))*MBDT
      XHE(I,ktf)=HEO(I,ktf)
      XQ(I,ktf)=QO(I,ktf)
      XT(I,ktf)=TN(I,ktf)
      IF(XQ(I,ktf).LE.0.)XQ(I,ktf)=1.E-08
      endif
      enddo
!
!--- calculate moist static energy, heights, qes
!
      call cup_env(xz,xqes,xhe,xhes,xt,xq,po,z1, &
           psur,ierr,tcrit,-1,xl,cp,   &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
!
!--- environmental values on cloud levels
!
      call cup_env_clev(xt,xqes,xq,xhe,xhes,xz,po,xqes_cup,xq_cup, &
           xhe_cup,xhes_cup,xz_cup,po_cup,gamma_cup,xt_cup,psur,   &
           ierr,z1,xl,rv,cp,          &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
!
!
!**************************** static control
!
!--- moist static energy inside cloud
!
!     do i=its,itf
!       if(ierr(i).eq.0)then
!         xhkb(i)=xhe(i,k22(i))
!       endif
!     enddo
      do k=kts,ktf
      do i=its,itf
         xhc(i,k)=0.
         xDBY(I,K)=0.
      enddo
      enddo
      do i=its,itf
        if(ierr(i).eq.0)then
!        if(use_excess == 2) then
!            k1=max(1,k22(i)-1)
!            k2=k22(i)+1
!            xhkb(i) =sum(xhe_cup(i,k1:k2))/float(k2-k1+1)+xl*zqexec(i)+cp*ztexec(i)
!        else if(use_excess <= 1) then
!            xhkb(i)=xhe_cup(i,k22(i))+float(use_excess)*(xl*zqexec(i)+cp*ztexec(i))
!        endif

         do k=1,kbcon(i)-1
            xhc(i,k)=xhkb(i)
         enddo
          k=kbcon(i)
          xhc(i,k)=xhkb(i)
          xDBY(I,Kbcon(i))=xHkb(I)-xHES_cup(I,K)
        endif !ierr
      enddo
!
!
      do i=its,itf
      if(ierr(i).eq.0)then
!     xzu(i,:)=zuo(i,:)
      xzu(i,1:ktf)=zuo(i,1:ktf)        !ss 2/19/14
      do k=kbcon(i)+1,ktop(i)
       xhc(i,k)=(xhc(i,k-1)*xzu(i,k-1)-.5*up_massdetro(i,k-1)*xhc(i,k-1)+ &
                         up_massentro(i,k-1)*xhe(i,k-1))   /            &
                         (xzu(i,k-1)-.5*up_massdetro(i,k-1)+up_massentro(i,k-1)+1.e-8)
       xdby(i,k)=xhc(i,k)-xhes_cup(i,k)
      enddo
      do k=ktop(i)+1,ktf
           xHC(i,K)=xhes_cup(i,k)
           xDBY(I,K)=0.
           xzu(i,k)=0.
      enddo
      endif
      enddo

!
!--- workfunctions for updraft
!
      call cup_up_aa0(xaa0,xz,xzu,xdby,GAMMA_CUP,xt_cup, &
           kbcon,ktop,ierr,           &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
!
! now for shallow forcing
!
       do i=its,itf
        xmb(i)=0.
        xff_shal(1:9)=0.
        if(ierr(i).eq.0)then

!         xmbmax(i)=0.1
          xmbmax(i)=100.*(p(i,kbcon(i))-p(i,kbcon(i)+1))/(g*dtime)
!----
!         write(0,*)'xmbmax = ',xmbmax(i),100.*(p_cup(i,kbcon(i))-p_cup(i,kbcon(i)+1))
!         xkshal=(xaa0(i)-aa1(i))/mbdt
!         if(xkshal.ge.0.)xkshal=+1.e6
!          if(xkshal.gt.-1.e-4 .and. xkshal.lt.0.)xkshal=-1.e-4
!         xff_shal(1)=max(0.,-(aa1(i)-aa0(i))/(xkshal*dtime))
!         xff_shal(1)=min(xmbmax(i),xff_shal(1))
!         xff_shal(2)=max(0.,-(aa1(i)-aa0(i))/(xkshal*dtime))
!         xff_shal(2)=min(xmbmax(i),xff_shal(2))
!         xff_shal(3)=max(0.,-(aa1(i)-aa0(i))/(xkshal*dtime))
!         xff_shal(3)=min(xmbmax(i),xff_shal(3))
!         if(aa1(i).le.0)then
!          xff_shal(1)=0.
!          xff_shal(2)=0.
!          xff_shal(3)=0.
!         endif
!         if(aa1(i)-aa0(i).le.0.)then
!          xff_shal(1)=0.
!          xff_shal(2)=0.
!          xff_shal(3)=0.
!         endif
!----
!
!- closure from Grant (200X)
           xff_shal(1)=.03*zws(i)
           xff_shal(2)=.03*zws(i)
           xff_shal(3)=.03*zws(i)
!----
!- closure from boundary layer QE (Raymond 1995)
            blqe=0.
            trash=0.
            if(k22(i).lt.kpbl(i)+1)then
               do k=1,kbcon(i)-1
                  blqe=blqe+100.*dhdt(i,k)*(po_cup(i,k)-po_cup(i,k+1))/g
               enddo
               trash=max((hc(i,kbcon(i))-he_cup(i,kbcon(i))),1.e1)
               xff_shal(7)=max(0.,blqe/trash)
               !print*,"blqe=", xff_shal(7),blqe,trash

            else
             xff_shal(7)=0.0
            endif
            xff_shal(8)= xff_shal(7)
            xff_shal(9)= xff_shal(7)
!----
!- closure from the heat-engine principle
!- Renn� and Ingersoll(1996), Souza et al (1999)
           !- get the averaged environment temperature between cloud base
           !- and cloud top
           tcold=0.
           do k=kbcon(i),ktop(i)
             dp   = po_cup(i,k)-po_cup(i,k+1)
             tcold= tcold + t_cup(i,k)*dp
           enddo
           tcold=tcold/(po_cup(i,kbcon(i))-po_cup(i,ktop(i)+1))

           !-surface temperature
           thot=tsur(i)  ! + ztexec(i)
           !- thermodynamic eficiency
           !efic = max(0.1, (thot-tcold)/thot )
           efic = max(0.0, (thot-tcold)/thot )

           !- total heat flux from surface
           fin = max(0.0, h_sfc_flux(i)+le_sfc_flux(i))

           !- mass flux at cloud base
           !if(cape(i) > 0.0 .and. h_sfc_flux(i) >0.0 ) then
           if(cape(i) > 0.0  ) then
            xff_shal(4) = efic * fin / cape(i)
            !if(xff_shal(4)>0.0001) then
            ! print*,"xff_shal1=",i,h_sfc_flux(i),le_sfc_flux(i),tsur(i)
            ! print*,"xff_shal2=",xff_shal(4) , efic , fin , cape(i)
            ! call flush(6)
            !endif
           else
            xff_shal(4) = 0.0
           endif
           xff_shal(5)=xff_shal(4)
           xff_shal(6)=xff_shal(4)

           !print*,"mshflx=",xff_shal(4) ,xff_shal(7),xff_shal(1)
!----
!          if(xkshal.lt.-1.1e-04)then ! .and.  &
!            ((aa1(i)-aa0(i).gt.0.) .or. (xff_shal(7).gt.0)))then
!          xff_shal(4)=max(0.,-aa0(i)/(xkshal*tscl_KF))
!          xff_shal(4)=min(xmbmax(i),xff_shal(4))
!          xff_shal(5)=xff_shal(4)
!          xff_shal(6)=xff_shal(4)
!          else
!           xff_shal(4)=0.
!           xff_shal(5)=0.
!           xff_shal(6)=0.
!          endif
!         write(0,888)'i0=',i,j,kpbl(i),blqe,xff_shal(7)
!888       format(a3,3(1x,i3),2e12.4)
!          xff_shal(4)= xff_shal(7)
!          xff_shal(5)= xff_shal(1)
!          xff_shal(6)= (xff_shal(1)+xff_shal(7))*.5
!----


          fsum=0.
          do k=1,9
           !srf- heat engine closure is providing too low values for mass fluxes.
           !srf- until this is checked, the ensemble closure will be calculatd
           !srf- only using closures BLQE and Wstar
           if(k.ge.4 .and. k.le.6) cycle

           xmb(i)=xmb(i)+xff_shal(k)
           fsum=fsum+1.
          enddo
          !- ensemble average of mass flux
          xmb(i)=min(xmbmax(i),xmb(i)/fsum)

          if(ichoice.gt.0)xmb(i)=min(xmbmax(i),xff_shal(ichoice))

          if(i.eq.ipr.and.j.eq.jpr)write(0,*)',ierr,xffs',ierr(i),xff_shal(1:9),xmb(i),xmbmax(i)
          if(xmb(i).eq.0.)ierr(i)=22
          if(xmb(i).eq.0.)ierrc(i)='22'
          if(xmb(i).lt.0.)then
             ierr(i)=21
             ierrc(i)='21'
             write(0,*)'neg xmb,i,j,xmb for shallow = ',i,j,k22(i),ierr(i)
          endif
        endif

        if(ierr(i).ne.0)then
           k22(i)=0
           kbcon(i)=0
           ktop(i)=0
           xmb(i)=0
           do k=kts,ktf
              outt(i,k)=0.
              outq(i,k)=0.
              outqc(i,k)=0.
           enddo
        else if(ierr(i).eq.0)then
!
! got the mass flux, sanity check, first for heating rates
!
          trash=0.
!         kmaxx=0
          do k=2,ktop(i)
           trash=max(trash,86400.*(dsubt(i,k)+dellat(i,k))*xmb(i))
          enddo
          if(trash.gt.100.)then
             xmb(i)=xmb(i)*100./trash
          endif
          trash=0.
          do k=2,ktop(i)
           trash=min(trash,86400.*(dsubt(i,k)+dellat(i,k))*xmb(i))
          enddo
          if(trash.lt.-100.)then
              xmb(i)=-xmb(i)*100./trash
          endif
!
! sanity check on moisture tendencies: do not allow anything that may allow neg
! tendencies
!
          do k=2,ktop(i)
           trash=q(i,k)+(dsubq(i,k)+dellaq(i,k))*xmb(i)*dtime
          if(trash.lt.1.e-12)then
! max allowable tendency over tendency that would lead to too small mix ratios
!
            trash=(1.e-12 -q(i,k))/((dsubq(i,k)+dellaq(i,k))*dtime)
            xmb(i)=(1.e-12 -q(i,k))/((dsubq(i,k)+dellaq(i,k))*dtime)
          endif
          enddo
          xmb_out(i)=xmb(i)
!
! final tendencies
!
          do k=2,ktop(i)
           outt (i,k)=(dsubt(i,k)+dellat(i,k))*xmb(i)
           outq (i,k)=(dsubq(i,k)+dellaq(i,k))*xmb(i)
           outqc(i,k)=(dellaqc(i,k)          )*xmb(i)
          enddo
        endif
       enddo

!- checking max/min values for outt tendencies
!       do i=its,itf
!         if(ierr(i) == 0)then
!           !-check max/min values for outt
!           if(minval(outt(i,:)*86400.) < -100.) &
!             print*,"==> low sh outt=",i, minval(outt(i,:)*86400.) &
!             ;call flush(6)
!           !-check max/min values for outt
!           if(maxval(outt(i,:)*86400.) > +100.) &
!             print*,"==> high sh outt=",i, maxval(outt(i,:)*86400.)&
!             ;call flush(6)
!         endif
!       enddo

!-sub-grid scale ice/liquid cloud water for radiation
 do i=its,itf
   ierr4d(i)   = ierr(i)
   xmb4d(i)    =  xmb(i)
   do k=1,ktf
        zup5d  (k,i) = zuo (i,k)
        clwup5d(k,i) = qrco(i,k) ! ice/liquid water
 enddo;enddo
 !---data saved for the tracer convective transport:
    do i=its,itf
       jmin4d(i)   = 1 ! no downdrafts
       kdet4d(i)   = 1 ! no downdrafts
       k224d(i)    = k22(i)
       kbcon4d(i)  = kbcon(i)
       ktop4d(i)   = ktop(i)
       kstabi4d(i) = kstabi(i)
       kstabm4d(i) = kstabm(i)
       !kpbl4d(i)  = kpbl(i)
       kpbl4d(i)   = k22(i)
       edt4d(i)    = 0.  ! no downdrafts
       pwav4d(i)   = 0.  ! no downdrafts
      do k=1,ktf
        !zcup5d(k,i) = zo_cup(i,k) !not in use
        pcup5d(k,i) = po_cup(i,k)
        up_massentr5d(k,i) = up_massentro(i,k)
        up_massdetr5d(k,i) = up_massdetro(i,k)
        dd_massentr5d(k,i) = 0.
        dd_massdetr5d(k,i) = 0.
        zdn5d  (k,i) = 0.
        prup5d (k,i) =   0.  ! no precip for shallow
        prdn5d (k,i) =   0.  ! no precip for shallow
!srf-30102018
!       tup5d  (k,i) =   tempco(i,k) !temperatura da parcela
!                                   ! de ar no updraft
        tup5d  (k,i) =   t_cup(i,k) !temperatura da parcela
!srf-30102018
      enddo
      !if(ierr(i)==0) print*,"1shallow",j,i,ierr4d(i),xmb4d(i),up_massentr5d(kbcon(i),i)
    enddo
!
!- begin: for GATE soundings-------------------------------------------
  if(use_gate) then
    do i=its,itf
     if(ierr(i).eq.0) then
      !- 2-d section
      do k=kts,ktop(i)

       nvar=60
       call setgradsvar(jl,k,nvar,xmb(i)*zuo(i,k),"mupsh",'m flux up sh (kg/s/m^2)')

       nvar=61
       call setgradsvar(jl,k,nvar,dellah(i,k),"delh_sh" ,'dellah_sh')

       nvar=62
       call setgradsvar(jl,k,nvar,dellaq(i,k)*86400.*xl/cp*xmb(i), "dellq_sh" ,'dellaq_sh')

       nvar=63
       call setgradsvar(jl,k,nvar,dellaqc(i,k)*86400.*xl/cp*xmb(i),"dellqc_sh" ,'dellaqc_sh')

       nvar=64
       call setgradsvar(jl,k,nvar,up_massentro(i,k),"upent_sh" ,'up_massentr_sh')

       nvar=65
       call setgradsvar(jl,k,nvar,up_massdetro(i,k),"updet_sh" ,'up_massdetr_sh')

       nvar=66
       call setgradsvar(jl,k,nvar,outt(i,k)*86400.,"outt_sh" ,'outt_sh K/s')

       nvar=67
       call setgradsvar(jl,k,nvar,outq(i,k)*86400.*xl/cp,"outq_sh" ,'outq_sh K/s')

       nvar=68
       call setgradsvar(jl,k,nvar,outqc(i,k)*86400.*xl/cp,"outqc_sh" ,'outqc_sh K/s')

      enddo
     endif
   enddo
  endif
!- end  : for GATE soundings-------------------------------------------

! done shallow
END SUBROUTINE CUP_gf_sh
!
!--------------------------------------------------------
   SUBROUTINE cup_up_aa1bl(aa0,t,tn,q,qo,dtime, &
              z,zu,dby,GAMMA_CUP,t_cup,         &
              kbcon,ktop,ierr,                  &
              itf,jtf,ktf,                      &
              its,ite, jts,jte, kts,kte         )

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
        itf,jtf,ktf,                                     &
        its,ite, jts,jte, kts,kte
  ! aa0 cloud work function
  ! gamma_cup = gamma on model cloud levels
  ! t_cup = temperature (Kelvin) on model cloud levels
  ! dby = buoancy term
  ! zu= normalized updraft mass flux
  ! z = heights of model levels
  ! ierr error value, maybe modified in this routine
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        z,zu,gamma_cup,t_cup,dby,t,tn,q,qo
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        kbcon,ktop
     real, intent(in) :: dtime
!
! input and output
!


     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
     real,    dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
        aa0
!
!  local variables in this routine
!

     integer                              ::                           &
        i,k
     real                                 ::                           &
        dz,dA
!
        DO i=its,itf
         AA0(I)=0.
        ENDDO
        DO k=kts+1,ktf
        DO i=its,itf
         IF(ierr(i).ne.0 ) cycle!GO TO 100
         IF(k.gt.KBCON(i)) cycle !GO TO 100

         DZ=Z(I,K)-Z(I,K-1)
         !print*,"dz=",i,k,z(i,k),Z(I,K-1),dz
         !da=zu(i,k)*DZ*(9.81/(1004.*( &
         !        (T_cup(I,K)))))*DBY(I,K-1)/ &
         !     (1.+GAMMA_CUP(I,K))
         ! IF(K.eq.KTOP(I).and.da.le.0.)go to 100

         dA=  DZ*9.81*( tn(i,k)-t(i,k) + 0.608*(qo(i,k)-q(i,k)))/dtime
         AA0(I)=AA0(I)+dA
      end do
      end do
!100     CONTINUE

   END SUBROUTINE cup_up_aa1bl
!----------------------------------------------------------------------
  SUBROUTINE get_zi_gf(j,its,ite,kts,kte,istart,iend,ktf,ierr,kzi,pbl,tkeg, &
                  rcpg,z,ztop,tkmin)

  implicit none
  integer,intent(in):: j,its,ite,kts,kte, ktf,istart,iend
  integer :: kzimax,ktke_max,i,k
  real tkmin,tke_tmp
  real,    dimension(its:ite,kts:kte) :: tkeg,rcpg,z
  real,    dimension(its:ite)          :: ztop,pbl
  integer, dimension(its:ite)          :: kzi,ierr

  real, parameter :: rcpmin=1.e-5 , pblhmax=3000.

  do i=istart,iend
     kzi(i)  = 1
!    if(ierr(i).eq.0)then
!         tke_tmp = 0.
         ktke_max= 1
         kzimax  =ktf-1
         !---  max level for kzi
         DO K=kts,ktf
           if(z(i,k).ge. pblhmax+ztop(i)) then
              kzimax = min(k,ktf-1)
              !if(j==8 .and. i==10) print*,"1",z(i,k), pblhmax,ztop(i),kzimax
              exit
           endif
         enddo
         !---

!         !level of max tke  below kzimax and w/out clouds
!         do  k=kts,kzimax
!           !print*,k,tkeg(i,k), tke_tmp,ktke_max,kzimax
!           if(rcpg(i,k) .lt. rcpmin) then
!             if( tkeg(i,k) .ge. tke_tmp) then
!               tke_tmp = tkeg(i,k)
!               cycle
!             else
!               ktke_max= max(1,k-1)
!               exit
!             endif
!           endif
!         enddo
!    201         continue
!             print*,ktke_max

         do k=ktke_max,kzimax+1

!           if(tkeg(i,k) .gt. 1.1*tkmin .and. rcpg(i,k) .lt. rcpmin )  then
            if(tkeg(i,k) .gt. 1.1*tkmin )  then
              kzi(i) = k
              !if(j==8 .and. i==10) print*,"I",k,rcpg(i,k),tkeg(i,k),kzi(i),z(i,k)-ztop(i)
              cycle

            else
               kzi(i) = max(1,k-1)
               !if(j==8 .and. i==10) print*,"F",k,rcpg(i,k),tkeg(i,k),kzi(i),z(i,k)-ztop(i)
               exit
            endif


         enddo
         kzi(i) = max(1     ,kzi(i))
         !print*,"1",kzi(i),i
         kzi(i) = min(kzimax,kzi(i))
         !print*,"2",kzi(i),i;call flush(6)
         pbl(i) = max( z(i,kzi(i))-ztop(i), z(i,1)-ztop(i) )
 enddo

 END SUBROUTINE get_zi_gf
!-------------------------------------------------------------------------
 SUBROUTINE rates_up_curvefitting(draft,ktop,ierr,p_cup,entr_rate_2d,hkbo,heo,heso_cup,z_cup, &
                            imid,kstabi,k22,kbcon,its,ite,itf,kts,kte,ktf,zuo,zustart,zufinal,zutop)
     implicit none
     integer, intent(in) :: imid,its,ite,itf,kts,kte,ktf
     real, dimension (its:ite,kts:kte),intent (inout) :: entr_rate_2d,zuo
     real, dimension (its:ite,kts:kte),intent (in) ::p_cup, heo,heso_cup,z_cup
     real, dimension (its:ite),intent (in) :: hkbo
     integer, dimension (its:ite),intent (in) :: kstabi,k22,kbcon
     integer, dimension (its:ite),intent (inout) :: ierr,ktop
     real, dimension (its:ite,kts:kte) :: hcot
     real, intent (in) ::zustart,zufinal,zutop
     character*(*), intent(in) ::draft

     real :: dz,dh, &
             dbythresh
     integer :: i,k,ipr,kdefi,kstart,kbegzu,kfinalzu
     !integer :: lev_adj,power_entr
     dbythresh=0.
     DO i=its,itf
      if(ierr(i).eq.0)then
        hcot(i,kbcon(i))=hkbo(i)
        do k=kbcon(i)+1,ktf-2
           dz=z_cup(i,k)-z_cup(i,k-1)

           hcot(i,k)=( (1.-0.5*entr_rate_2d(i,k-1)*dz)*hcot(i,k-1) &
                      + entr_rate_2d(i,k-1)*dz*heo(i,k-1))/ &
                      (1.+0.5*entr_rate_2d(i,k-1)*dz)
       enddo
       do k=kbcon(i),ktf-2
          if((hcot(i,k)-heso_cup(i,k)).lt.dbythresh)then
              kfinalzu=k  - 1
              ktop(i)=kfinalzu
              go to 412
          endif
       enddo
       kfinalzu=ktf-2
       ktop(i)=kfinalzu
412    continue
       if(kfinalzu.le.kbcon(i)+2)then
        ierr(i)=41
        ktop(i)= 0
       else
        call get_zu_zd_pdf_orig(draft,ierr(i),kbcon(i),kfinalzu,zustart,zufinal,zutop,zuo(i,kts:kte),kts,kte,ktf)

       endif
      ENDIF
     ENDDO
  END SUBROUTINE rates_up_curvefitting
!-----------------------------------------------------------
 SUBROUTINE rates_up_pdf(name,ktop,ierr,p_cup,entr_rate_2d,hkbo,heo,heso_cup,z_cup, &
                            kstabi,k22,kbcon,its,ite,itf,kts,kte,ktf,zuo,kpbl)
     implicit none
     integer, intent(in) :: its,ite,itf,kts,kte,ktf
     real, dimension (its:ite,kts:kte),intent (inout) :: entr_rate_2d,zuo
     real, dimension (its:ite,kts:kte),intent (in) ::p_cup, heo,heso_cup,z_cup
     real, dimension (its:ite),intent (in) :: hkbo
     integer, dimension (its:ite),intent (in) :: kstabi,k22,kbcon,kpbl
     integer, dimension (its:ite),intent (inout) :: ierr,ktop
     real, dimension (its:ite,kts:kte) :: hcot
      character *(*), intent (in)         ::                           &
       name
     real :: dz,dh, &
             dbythresh
     real :: dby(kts:kte)
     integer :: i,k,ipr,kdefi,kstart,kbegzu,kfinalzu

     dbythresh=1.0  ! the range of this parameter is 0-1, higher => lower
                    ! overshoting (cheque AA0 calculation)
                    ! rainfall is too sensible this parameter
                    ! for now, keep =1.
     if(name == 'shallow')dbythresh=1.0

     DO i=its,itf
      dby(:)=0.0
      if(ierr(i).eq.0)then
        hcot(i,kbcon(i))=hkbo(i)
        dz=z_cup(i,kbcon(i))-z_cup(i,kbcon(i)-1)
        dby(kbcon(i))=(hcot(i,kbcon(i))-heso_cup(i,kbcon(i)))*dz

        do k=kbcon(i)+1,ktf-2
           dz=z_cup(i,k)-z_cup(i,k-1)

           hcot(i,k)=( (1.-0.5*entr_rate_2d(i,k-1)*dz)*hcot(i,k-1) &
                      + entr_rate_2d(i,k-1)*dz*heo(i,k-1))/ &
                      (1.+0.5*entr_rate_2d(i,k-1)*dz)
           dby(k)=dby(k-1)+(hcot(i,k)-heso_cup(i,k))*dz

       enddo
!-orig
!       do k=kbcon(i),ktf-2
!          if((hcot(i,k)-heso_cup(i,k)).lt.dbythresh)then
!              kfinalzu=k  - 1
!              ktop(i)=kfinalzu
!              go to 412
!          endif
!       enddo
!
!-new
       do k=maxloc(dby(:),1),ktf-2
          if(dby(k).lt.dbythresh*maxval(dby))then
              kfinalzu = k - 1
              ktop(i)  = kfinalzu
              go to 412
          endif
       enddo
!-
       kfinalzu=ktf-2
       ktop(i)=kfinalzu
412    continue
       if( name == 'deep' ) then
          if(kfinalzu.le.kbcon(i)+2)then
              ierr(i)=41
              ktop(i)= 0
          else
           call get_zu_zd_pdf("UP",ierr(i),kbcon(i),kfinalzu,zuo(i,kts:kte),kts,kte,ktf,kpbl(i))
          endif
       elseif ( name == 'shallow' ) then
          if(kfinalzu.le.kbcon(i)+1)then
              ierr(i)=41
              ktop(i)= 0
              !write(0,*)'!!!!!shallow kfinalzu,kbcon = ',i,kfinalzu,kbcon(i)
          else
              !write(0,*)'shallow ktop = ',i,ktop(i),kbcon(i)
              call get_zu_zd_pdf("SH1",ierr(i),k22(i),kfinalzu,zuo(i,kts:kte),kts,kte,ktf,kpbl(i))
          endif
         endif
      ENDIF
     ENDDO
  END SUBROUTINE rates_up_pdf
!-------------------------------------------------------------------------
subroutine get_zu_zd_pdf(draft,ierr,kb,kt,zu,kts,kte,ktf,kpbli)

implicit none
integer, intent(in) ::kb,kt,kts,kte,ktf,kpbli
real, intent(inout) :: zu(kts:kte)
real  :: zuh(kts:kte)
integer, intent(inout) :: ierr
character*(*), intent(in) ::draft

!- local var
integer :: kk,add,i,nrec=0,k,kb_adj,kpbli_adj
real ::zumax,ztop_adj
real ::a2,beta, alpha,kratio,tunning,FZU,krmax,dzudk
!- kb cannot be at 1st level

   !-- fill zu with zeros
   zu=0.0
   zuh=0.0

IF(draft == "UP") then
   kb_adj=max(kb,2)
  !beta=4.  !=> must be larger than 1
            !=> higher makes the profile sharper
            !=> around the maximum zu
  !- 2nd approach for beta and alpha parameters
  !- the tunning parameter must be between 0.5 (low  level max zu)
  !-                                   and 1.5 (high level max zu)
  tunning = 1.0
!  tunning = 0.9
  beta    = 2.0/tunning
  alpha   = tunning*beta
  !
  !- this alpha constrains the location of the maximun ZU to be at
  !- "kb_adj" vertical level
  !alpha=1. + (beta-1.0)*(float(kb_adj)/float(kt+1))/(1.0-(float(kb_adj)/float(kt+1)))
  !
  fzu=1.
  do k=kts+1,min(kte,kt+1)
      kratio= float(k)/float(kt+1)
      zu(k) = FZU*kratio**(alpha-1.0) * (1.0-kratio)**(beta-1.0)
   enddo

ELSEIF(draft == "SH1") then
  tunning = 0.6
  beta    =2.2/tunning
  alpha   = tunning*beta
  fzu=1.
  do k=kts,min(kt+1,ktf)
      kratio= float(k)/float(kt+1)
      zuh(k) = FZU*kratio**(alpha-1.0) * (1.0-kratio)**(beta-1.0)
   enddo
   do k=maxloc(zuh(:),1),1,-1
      kk=kb+k-maxloc(zuh(:),1)
      if(kk.gt.1)zu(kk)=zuh(k)
   enddo
   do k=maxloc(zuh(:),1)+1,kt
      kk=kb+k-maxloc(zuh(:),1)
      if(kk.le.kt)zu(kk)=zuh(k)
   enddo

   !do k=kts,min(kt+1,ktf)
   !  print*,"shallow",k,zu(k),kt
   !enddo
ELSEIF(draft == "SH") then
   !alpha= 3.
   !beta = 2.*alpha
   kb_adj=kb ! level where mass flux starts
   kpbli_adj=kpbli
   if(kpbli_adj < kb_adj) then
      !print*,"kpbli < kb", kpbli_adj, kb_adj
      kpbli_adj = kb_adj + 1
   endif
  ! kpbli_adj=max(kpbli_adj,kb_adj + 1)

   !- location of the maximum Zu
   krmax=float(kpbli_adj-kb_adj+1)/float(kt+1)
   !
   beta= 4.
   !- this alpha imposes the maximum zu at kpbli
   alpha=1.+krmax*(beta-1.)/(1.-krmax)

   !- to check if dZu/dk = 0 at k=kpbli_adj
   !kratio=krmax
   !dzudk=(alpha-1.)*(kratio**(alpha-2.)) * (1.-kratio)**(beta-1.) - &
   !          (kratio**(alpha-1.))*((1.-kratio)**(beta-2.))*(beta-1.)

   !- Beta PDF
!srf   do k=kts+kb_adj-1,min(kte,kt+1)
   do k=1+kts+kb_adj-1,min(kte,kt+1)
      kratio=float(k+1-kb_adj)/float(kt+1)  !-kb_adj+1)

      zu(k)=kratio**(alpha-1.0) * (1.0-kratio)**(beta-1.0)
   enddo
   !print*,"========================================="
   !   zu(kts:min(kte,kt+1))= zu(kts:min(kte,kt+1))/ maxval(zu(kts:min(kte,kt+1)))
   !  do k=kts,min(kt+1,ktf)
   !    print*,"shallow2",k,zu(k),kt,kb_adj,kpbli;call flush(6)
   !  enddo
   !print*,"========================================="

ELSEIF(draft == "DOWN" .or. draft == "DOWNM") then

 tunning = 0.8
 beta    = 3.0/tunning
!  tunning = 2.0
!  beta    =4.0/tunning
  alpha   = tunning*beta
  fzu=1.
  zuh(:)=0.
  do k=kts,min(kt+1,ktf)
      kratio= float(k)/float(kt+1)
      zuh(k) = FZU*kratio**(alpha-1.0) * (1.0-kratio)**(beta-1.0)
      !write(0,*)k,zuh(k)
   enddo
   if(maxloc(zuh(:),1).ge.kb)then
      do k=maxloc(zuh(:),1),1,-1
         kk=kb+k-maxloc(zuh(:),1)
         if(kk.gt.1)zu(kk)=zuh(k)
         !if(kk.gt.1)write(0,*)kk,zu(kk)
      enddo
      do k=maxloc(zuh(:),1)+1,kt
         kk=kb+k-maxloc(zuh(:),1)
         if(kk.le.kt)zu(kk)=zuh(k)
         !if(kk.ge.1)write(0,*)kk,zu(kk)
      enddo
   else
      do k=2,kt ! maxloc(zuh(:),1)
        zu(k)=zuh(k-1)
        !write(0,*)k,zu(k)
      enddo
   endif
ENDIF


   !- normalize ZU
   zu(kts:min(kte,kt+1))= zu(kts:min(kte,kt+1))/ maxval(zu(kts:min(kte,kt+1)))

return


end subroutine get_zu_zd_pdf
!-------------------------------------------------------------------------
subroutine get_zu_zd_pdf_orig(draft,ierr,kb,kt,zs,zuf,ztop,zu,kts,kte,ktf)

implicit none
integer, intent(in) ::kb,kt,kts,kte,ktf
real, intent(in) :: Zs,Zuf,Ztop
real, intent(inout) :: zu(kts:kte)
integer, intent(inout) :: ierr
character*(*), intent(in) ::draft

!- local var
integer :: add,i,nrec=0,k,kb_adj
real ::zumax,ztop_adj
real ::beta, alpha,kratio,tunning

!- kb cannot be at 1st level
kb_adj=max(kb,2)

!-- fill zu with zeros
zu=0.0

IF(draft == "UP" .or. draft == "up" ) then
  if(kt<=kb_adj) then
    !stop "ktop must be larger than kbcon"
    ierr=99
    return
  endif
  !beta=4.  !=> must larger than 1
            !=> higher makes the profile sharper
            !=> around the maximum zu
  add=0     !=> additional levels above kbon, where
            !=> the maximum zu will resides
  kb_adj=kb_adj+add
  kb_adj=max(10,kb_adj)

  !- this alpha constrains the location of the maximun ZU to be at
  !- "kb_adj" vertical level
  !alpha=1. + (beta-1.0)*(float(kb_adj)/float(kt+1))/(1.0-(float(kb_adj)/float(kt+1)))

  !- 2nd approach for beta and alpha parameters
  !- the tunning parameter must be between 0.5 (low  level max zu)
  !-                                   and 1.5 (high level max zu)
  tunning = 0.6
  beta    = 2.0/tunning
  alpha   = tunning*beta

   !- Beta PDF
  do k=kts,min(kte,kt+1)
      kratio= float(k)/float(kt+1)

      zu(k) = kratio**(alpha-1.0) * (1.0-kratio)**(beta-1.0)
  enddo

ELSEIF(draft == "DOWN" .or. draft == "DOWNM") then
  add=0    !=> additional levels above kbcon, where
           !=> the maximum zu will resides
  beta=4.  !=> must larger than 1
           !=> higher makes the profile sharper
           !=> around the maximum zu
  alpha= 0.25*beta

  !- for downdrafts kb = jmin(i)-levadj
  kb_adj=kb_adj+add

  !- 2nd approach for beta and alpha parameters
  !- the tunning parameter must be between 0.5 (low  level max zu)
  !-                                   and 1.5 (high level max zu)
  tunning = 1.
  beta    = 2.0/tunning
  alpha   = tunning*beta

  !- Beta PDF
  do k=kts,min(kte,kt)
      kratio= float(k)/float(kt)

      zu(k+1) = kratio**(alpha-1.0) * (1.0-kratio)**(beta-1.0)
  enddo

ELSEIF(draft == "shallow" .or. draft == "SHALLOW") then

 alpha= 3.
 beta = 2.*alpha
 kb_adj=1 ! level where mass flux starts

 !- Beta PDF
 do k=kts+kb_adj-1,min(kte,kt+1)
    kratio=float(k+1-kb_adj)/float(kt+1)  !-kb_adj+1)

    zu(k)=kratio**(alpha-1.0) * (1.0-kratio)**(beta-1.0)
 enddo

ELSE
 print*, "unknown type of flow" ,draft
 stop "routine get_zu_zd"

ENDIF

   !- normalize ZU
   zu(kts:min(kte,kt+1))= zu(kts:min(kte,kt+1))/ maxval(zu(kts:min(kte,kt+1)))

   !--- Sanity checks
   if(beta <= 1) stop "beta must be larger than 1"

   if(minval(zu(:)) < 0.0 ) then
     print*," zu has negative values for ", draft
     stop   " zu < zero"
   endif
   if(maxval(zu(:)) > 1.0 ) then
     print*," zu has values greater than 1 for ", draft
     stop   " zu  >  one"
   endif

return

!OPEN(19,FILE= 'zu.gra', FORM='unformatted',ACCESS='direct'&
!       ,STATUS='unknown',RECL=4)
! DO k = kts,kte
!    nrec=nrec+1
!     WRITE(19,REC=nrec) zu(k)
! END DO
!close (19)

end subroutine get_zu_zd_pdf_orig
!----------------------------------------------------------------------
subroutine get_zu_zd_polynomial(draft,ierr,kb,kt,zs,zuf,ztop,zu,kts,kte,ktf)

implicit none
integer, intent(in) ::kb,kt,kts,kte,ktf
real, intent(in) :: Zs,Zuf,Ztop
real, intent(inout) :: zu(kts:kte)
integer, intent(inout) :: ierr
character*(*), intent(in) ::draft

!- local var
integer :: nmatrix
REAL,ALLOCATABLE, DIMENSION(:  ) :: B ! vector indep and solution
REAL,ALLOCATABLE, DIMENSION(:,:) :: Amatrix
integer :: add,i,nrec=0,k,kint,j,kt2,kb_adj
real ::zumax,zint,fx,fb,ztop_adj
real, parameter :: maxzu=1000.
real, parameter ::  threshold=-1.e-6
LOGICAL :: solution_accepted

!- default settings
solution_accepted = .false.
nmatrix=4

!- kb cannot be at 1st level
kb_adj=max(kb,2)


IF(draft == "UP") then
  if(kt<=kb_adj) then
    !stop "ktop must be larger than kbcon"
    ierr=99
    return
  endif
  fb=0.3 !=> must be in between 0 and 1
         !=> higher makes the maximum zu height lower
  if(kb_adj < nint(0.5*(kt-1))) then
      kint=max(kb_adj+1,nint(fb*kb_adj)+nint((1.-fb)*kt))
      kint=max(kint,nint(0.5*(kt-1)))
      kint=min(kt-3,kint)
  else
      !kint=min(kb_adj-3,nint(fb*kb_adj)+nint((1.-fb)*1))
      !kint=min(kb_adj-3,nint((1.-fb)*kb_adj)+nint(fb*1))
      !kint=max(3,kint)
      kint=max(nint(0.5*kb_adj),3)
  endif
  zint=zuf-0.2
  ztop_adj=(1.-fb)*ztop
  add=0

  if(kint ==  kb_adj  .or. kint ==  kt) nmatrix=3

ELSEIF(draft == "DOWN") then
  add=0
  fb=0.8 !=> must be in between 0 and 1
         !=> higher makes the maximum zu height lower
  !- for downdrafts kb = jmin(i)-levadj
  !- kint is betwen surface (k=1) and kb(= jmin(i)-levadj)
  kint=max(2,int(fb*1)+int((1.-fb)*kb_adj))
  zint=0.1 !zuf-0.7
  ztop_adj=(1.-fb)*ztop
  if(kint ==  kb_adj  .or. kint ==  kt) nmatrix=3

ELSEIF(draft == "DOWNM") then
  nmatrix=3
  add=0
  ztop_adj=ztop
ENDIF

!-- form of the matrix (3rd order)
!           | 1        1       1    1|
!Amatrix=   | kb  **3  kb  **2 kb   1|
!           | kint**3  kint**2 kint 1|
!           | kt  **3  kt  **2 kt   1|
!
!-
will_try_3rd_or_2nd_orders_only:   DO WHILE (solution_accepted .eqv. .false.)

   IF(.NOT. allocated(Amatrix)) then
    allocate (Amatrix(nmatrix,nmatrix)); amatrix=1.
    allocate (B(nmatrix))              ; B=0.
   ENDIF
   IF(size(B) .ne. nmatrix) then
    deallocate (Amatrix,B)
    allocate (Amatrix(nmatrix,nmatrix)); amatrix=1.
    allocate (B(nmatrix))              ; B=0.
   ENDIF
   !
   if(nmatrix-1 == 2) then ! 2nd order polynomial
     do j=1,nmatrix-1
      amatrix(2,j)=float(kb_adj**(nmatrix-j))
      amatrix(3,j)=float(kt    **(nmatrix-j))
     enddo
     ! vector indep
     B(1)=Zs
     B(2)=Zuf
     B(3)=Ztop_adj
   elseif(nmatrix-1 == 3) then ! 3rd order polynomial
     do j=1,nmatrix-1
       amatrix(2,j)=float(kb_adj**(nmatrix-j))
       amatrix(3,j)=float(kint  **(nmatrix-j))
       amatrix(4,j)=float(kt    **(nmatrix-j))
     enddo
     ! vector indep
     B(1)=Zs
     B(2)=Zuf
     B(3)=Zint
     B(4)=Ztop_adj
   else
    stop 333
   endif

   ! in : B = indep vector
   CALL solve(Amatrix,B,nmatrix)
   ! out: B = solution

   !-- fill zu with zeros
   zu=0.0

   !-- build ZU profile
   IF    (nmatrix-1 == 3) then ! 3rd order polynomial

     do k=kts,min(kte,kt+add)
       !zu(k)=MAX(0.0,MIN(maxzu,B(1)*float(k)**3+B(2)*float(k)**2+B(3)*float(k)+B(4)))
       zu(k)=B(1)*float(k)**3+B(2)*float(k)**2+B(3)*float(k)+B(4)
     enddo

   ELSEIF(nmatrix-1 == 2) then ! 2nd order polynomial

     do k=kts,min(kte,kt+add)
      !zu(k)=MAX(0.01,MIN(maxzu,B(1)*float(k)**2+B(2)*float(k)+B(3)))
      zu(k)=B(1)*float(k)**2+B(2)*float(k)+B(3)
     enddo
   ELSE
     stop 444
   ENDIF

   if(minval(zu(1:kt)) < threshold .and. nmatrix==4) then
     !- will try a solution using 2nd order (parabolic profile)
     nmatrix=3
     solution_accepted=.false.
   elseif(minval(zu(1:kt)) < threshold .and. nmatrix==3) then
     !- no acceptable solution was found using
     !- 2nd or 3rd order polynomials. Convection will not be allowed
     !- for this column...
     ierr=100
     exit will_try_3rd_or_2nd_orders_only
   else
     solution_accepted=.true.
   endif

ENDDO will_try_3rd_or_2nd_orders_only


!fx=maxval(zu)
!zu=zu/fx
! print*,"zu=",zu

!--- Sanity checks
if(abs(zuf-zu(int(kb_adj)))> 1.e-4 .or.abs(ztop_adj-zu(int(kt)))>1.e-4)then
 print*,"DRAFT= ",draft, " - vector solution B:"
 print*,"B=",B
 print*,zuf ,zu(int(kb_adj)),int(kb_adj),zuf - zu(int(kb_adj))
 print*,ztop_adj,zu(int(kt)),int(kt),ztop_adj-zu(int(kt))
 if(nmatrix-1 == 3) print*,zint,zu(int(kint)),int(kint)
 stop 331
end if

if(maxval(abs(b)) > 1000.) then
 print*,"non acceptable solution for B"
 print*,"nmatrix,B=",nmatrix,B
 print*,"draft= ",draft
 print*,kb_adj  ,kint,kt
 print*,Zs,Zuf,Zint,Ztop_adj
 print*,zu
 stop 333
ENDIF

!--- deallocate memory
! deallocate (Amatrix,B)
return

!OPEN(19,FILE= 'zu.gra', FORM='unformatted',ACCESS='direct'&
!       ,STATUS='unknown',RECL=4)
! DO k = kts,kte
!    nrec=nrec+1
!     WRITE(19,REC=nrec) zu(k)
! END DO
!close (19)

end subroutine get_zu_zd_polynomial

!==============================================================
SUBROUTINE Solve(coef,indep,n)
  !
  !Surboutine Solve. Solve a system of differential/linear equations
  IMPLICIT NONE

  INTEGER,INTENT(IN)              :: n
  REAL,DIMENSION(n,n),INTENT(IN)  :: coef  !Matrix coeficients to solve
  REAL,DIMENSION(n),INTENT(INOUT) :: indep !Independent terms at input and solutions at output

  REAL,DIMENSION(n,n) :: a
  INTEGER,DIMENSION(n) :: indx
  a=coef !Copy a input coeficients to manipulate array

  CALL ludecomp(a,n,indx)

  CALL lusolver(a,n,indep,indx)

END SUBROUTINE Solve
!==============================================================
SUBROUTINE ludecomp (a, n,indx)
  !
  !Subroutine ludecomp -  LU decomposition from a
  !                         L.U = A
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: n                    !Size of quadratic system
  REAL,DIMENSION(n,n),INTENT(INOUT) :: a      !Coeficients Terms of equation (input)
                                              !LU decomposition of matrix
                                              !(output)
  INTEGER,DIMENSION(N),INTENT(OUT) :: indx    !Vector with index of row permutation
  DOUBLE PRECISION ,PARAMETER :: zero=1.0D-30

  DOUBLE PRECISION ,DIMENSION(n) :: vv
  INTEGER :: i, imax, j, k
  DOUBLE PRECISION :: aamax, dum, temp

  do i = 1, n
    aamax = 0.
    do j = 1, n
      if (abs(a(i,j))>aamax) aamax=abs(a(i,j))
    end do
    if (aamax==0.) THEN
          Write (*,*) 'stopping: singular matrix in ludecomp'
          stop
    END IF
    vv(i) = 1./aamax
  end do

  do j = 1, n
    do i = 1, j - 1
      temp = a (i, j)
      do k = 1, i - 1
        temp = temp - a (i, k) * a (k, j)
      end do
      a (i, j) = temp
    end do

    aamax = 0.

    do i = j, n
      temp = a (i, j)
      do k = 1, j - 1
       temp = temp - a (i, k) * a (k, j)
      end do

      a (i, j) = temp

      dum = vv (i) * abs (temp)
      if (dum>=aamax) then
        imax = i
        aamax = dum
      end if

    end do

    if (j/=imax) then
      do k = 1, n
        dum = a (imax, k)
        a (imax, k) = a (j, k)
        a (j, k) = dum
      end do
      vv (imax) = vv (j)
    end if
    indx (j) = imax

    if (a (j, j)==0.) a (j, j) = zero
    if (j/=n) then
      dum = 1. / a (j, j)
      do i = j + 1, n
        a (i, j) = a (i, j) * dum
      end do
    end if
  end do
END subroutine ludecomp
!======================================
SUBROUTINE lusolver (a,n,b,indx)
  !
  !Subroutine lusolver -  forward substitution and backsubstitution
  !                   A · x = (L · U) · x = L · (U · x) = b
  IMPLICIT NONE
  INTEGER,INTENT(IN)              :: n           !Size of quadratic system
  INTEGER,DIMENSION(n),INTENT(IN) :: indx        !Index for pivot obtained from LU decomposition
  REAL,DIMENSION(n,n),INTENT(IN)  :: a           !Array obtained from LU decomposition
  REAL,DIMENSION(n),INTENT(INOUT) :: b           !b - Independent terms of equations (in input)
  INTEGER :: i, ii, j, ll
  DOUBLE PRECISION :: temp
  ii = 0
  do i = 1, n

    ll = indx (i)
    temp= b (ll)
    b (ll) = b (i)
    if (ii/=0) then
       do j = ii, i - 1
         temp= temp- a (i, j) * b (j)
       end do
    else if (temp/=0.) then
       ii = i
    end if
    b (i) = temp
  end do
  do i = n, 1, - 1
    temp= b (i)
    do j = i + 1, n
      temp= temp- a (i, j) * b (j)
    end do
    b (i) = temp/ a (i, i)
  end do
END subroutine lusolver
!----------------------------------------------------------------------
   SUBROUTINE cup_up_cape(aa0,z,zu,dby,GAMMA_CUP,t_cup,  &
              kbcon,ktop,ierr,tempco,qco,qrco, qo_cup,   &
              itf,jtf,ktf,its,ite, jts,jte, kts,kte      )

   IMPLICIT NONE
     integer ,intent (in   )                   ::        &
        itf,jtf,ktf, its,ite, jts,jte, kts,kte

  ! aa0 = dummy array for CAPE
  ! gamma_cup = gamma on model cloud levels
  ! t_cup = temperature (Kelvin) on model cloud levels
  ! dby = buoancy term
  ! zu= normalized updraft mass flux
  ! z = heights of model levels
  ! ierr = error value, maybe modified in this routine
  ! tempco = in-cloud temperature (Kelvin) on model cloud levels
  ! qco    = in-cloud water vapor mixing ratio on model cloud levels
  ! qo_cup = environ water vapor mixing ratio on model cloud levels
  ! qrco   = in-cloud liquid water mixing ratio on model cloud levels

     real,    dimension (its:ite,kts:kte) ,intent (in   )    ::        &
        z,zu,gamma_cup,t_cup,dby,tempco,qco,qrco, qo_cup
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        kbcon,ktop
!
! input and output
     integer, dimension (its:ite),intent (inout)   ::                  &
        ierr
     real,    dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
        aa0
!
!  local variables in this routine
   integer                              ::   i,k
   real                                 ::   dz,daa0
!
   integer, parameter :: cape_formulation = 2 ! = 1  like cloud work function
                                              ! = 2  traditional formulation
   DO i=its,itf
         AA0(i)=0.
   ENDDO
   IF(cape_formulation == 1) then

       DO k=kts+1,ktf
        DO i=its,itf
         IF(ierr(i).ne.0) cycle!GO TO 100
         IF(K.LT.KBCON(I)) cycle!GO TO 100
         IF(K.Gt.KTOP(I)) cycle!GO TO 100
         DZ=Z(I,K)-Z(I,K-1)
         daa0=  DZ*(9.81/(1004.*( &
                (T_cup(I,K)))))*DBY(I,K-1)/ &
             (1.+GAMMA_CUP(I,K))
         IF(K.eq.KTOP(I).and.daa0.le.0.) cycle !go to 100
         AA0(I)=AA0(I)+daa0
         if(aa0(i).lt.0.)aa0(i)=0.
       end do
       end do

   ELSEIF(cape_formulation == 2) then

        DO i=its,itf
          IF(ierr(i) == 0) then

             DO k=max(kbcon(i),2),ktop(i)
              DZ=Z(I,K)-Z(I,K-1)
              daa0=9.81*DZ*( (tempco(i,k)- t_cup(i,k) )/ t_cup(i,k) + &
                     0.608*  (qco   (i,k)-qo_cup(i,k) )/qo_cup(i,k) - &
                              qrco  (i,k)                             )
              AA0(I)=AA0(I)+daa0
              !print*,"cape",k,AA0(I),tempco(i,k),t_cup(i,k), qrco  (i,k),dz
            ENDDO
          ENDIF
        ENDDO
   ELSE
         stop "wrong option for cape calculation - must be 1 or 2"
   ENDIF
 END SUBROUTINE cup_up_cape
 !====================================================================
    SUBROUTINE cup_kbcon(ierrc,cap_inc,iloop,k22,kbcon,he_cup,hes_cup, &
              hkb,ierr,kbmax,p_cup,cap_max,                                &
              xl,cp,ztexec,zqexec,use_excess,                              &
              jprnt,itf,jtf,ktf,its,ite, jts,jte, kts,kte,                  &
              z_cup,entr_rate,heo)

   IMPLICIT NONE
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
        use_excess,jprnt,itf,jtf,ktf,           &
        its,ite, jts,jte, kts,kte
  !
  !
  !
  ! ierr error value, maybe modified in this routine
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        he_cup,hes_cup,p_cup
     real,    dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        ztexec,zqexec,cap_max,cap_inc
     real,intent (in   )                  ::                           &
        xl,cp
     real,    dimension (its:ite)                                      &
        ,intent (inout   )                   ::                           &
        hkb
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        kbmax
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        kbcon,k22,ierr
     integer                                                           &
        ,intent (in   )                   ::                           &
        iloop
     character*50 :: ierrc(its:ite)

     real, dimension (its:ite,kts:kte),intent (in) :: z_cup,heo
     real,intent (in) :: entr_rate
!
!  local variables in this routine
!

     integer                              ::                           &
        i,k,k1,k2
     real                                 ::                           &
        pbcdif,plus,hetest,dz
     real, dimension (its:ite,kts:kte) :: hcot
!
!--- DETERMINE THE LEVEL OF CONVECTIVE CLOUD BASE  - KBCON
!
      DO 27 i=its,itf
        kbcon(i)=1
        IF(ierr(I).ne.0)GO TO 27
        KBCON(I)=K22(I)+1
        if(iloop.eq.5)KBCON(I)=K22(I)

       !== including entrainment for hetest
        hcot(i,1:k22(i)) = HKB(I)
        do k=k22(i)+1,KBMAX(i)+3
           dz=z_cup(i,k)-z_cup(i,k-1)

           hcot(i,k)= ( (1.-0.5*entr_rate*dz)*hcot(i,k-1)   &
                         + entr_rate*dz*heo(i,k-1)       )/ &
                      (1.+0.5*entr_rate*dz)
        enddo
       !==

        GO TO 32
 31     CONTINUE
        KBCON(I)=KBCON(I)+1

        IF(KBCON(I).GT.KBMAX(i)+2)THEN
           if(iloop.ne.4)then
                ierr(i)=3
                ierrc(i)="could not find reasonable kbcon in cup_kbcon"
           endif
           GO TO 27
        ENDIF
 32     CONTINUE

       !==
       !hetest= hkb(i)
        hetest= hcot(i,KBCON(I))
       !==

       if(iloop.eq.5)then
           hetest=HKB(I)
!          do k=1,k22(i)
!          hetest=max(hetest,he_cup(i,k))
!          enddo
        endif

        IF(HETEST.LT.HES_cup(I,KBCON(I)))then
           if(jprnt.eq.1)write(0,*)'htest',k22(i),kbcon(i),HETEST,-P_cup(I,KBCON(I))+P_cup(I,K22(I))
           GO TO 31
        ENDIF

!       cloud base pressure and max moist static energy pressure
!       i.e., the depth (in mb) of the layer of negative buoyancy
        if(KBCON(I)-K22(I).eq.1)go to 27

        if(iloop.eq.5 .and. (KBCON(I)-K22(I)).le.2)go to 27

        PBCDIF=-P_cup(I,KBCON(I))+P_cup(I,K22(I))
        plus=max(25.,cap_max(i)-float(iloop-1)*cap_inc(i))

!----------
        if(iloop.eq.4)plus=cap_max(i)
        !
        ! for shallow convection, if cap_max is greater than 25, it is the pressure at pbltop
        if(iloop.eq.5)plus=150.
        if(iloop.eq.5.and.cap_max(i).gt.25)pbcdif=-P_cup(I,KBCON(I))+cap_max(i)
!----------

        IF(PBCDIF.GT.plus)THEN

            if(jprnt.eq.1)write(0,*)'htest',k22(i),kbcon(i),plus,-P_cup(I,KBCON(I))+P_cup(I,K22(I))

            K22  (I)=K22(I)+1
            KBCON(I)=K22(I)+1

            !==     including entrainment for hetest
            hcot(i,1:k22(i)) = HKB(I)
            do k=k22(i)+1,KBMAX(i)+3
               dz=z_cup(i,k)-z_cup(i,k-1)

               hcot(i,k)= ( (1.-0.5*entr_rate*dz)*hcot(i,k-1)        &
                                  + entr_rate*dz* heo (i,k-1)        )/ &
                            (1.+0.5*entr_rate*dz)
            enddo
            !==

            if(iloop.eq.5)KBCON(I)=K22(I)

            IF(KBCON(I).GT.KBMAX(i)+2)THEN
               if(iloop.ne.4)then
                 ierr(i)=3
                 ierrc(i)="could not find reasonable kbcon in cup_kbcon"
               endif
               GO TO 27
            ENDIF
            GO TO 32
        ENDIF
 27   CONTINUE

   END SUBROUTINE cup_kbcon

END MODULE module_cu_gf2
