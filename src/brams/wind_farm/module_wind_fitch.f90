!>
!!@package BRAMS-WindFram Main module
!!@author Luiz Flavio
!!@date 05/22/2017
!!@brief Set of routines to calc power produced by turbines and tke's disturbs
!!@copyright Under CC-GPL license
!!@link https://creativecommons.org/licenses/GPL/2.0/legalcode.pt
!!@param see especif routines
!!@par History and Info
!!
!! NOTICE
!! The following paper should be cited whenever presenting results using this scheme
!! (using either the original version or any modified versions of the scheme):
!!@citeFitch, A. C. et al. 2012: Local and Mesoscale Impacts of Wind Farms as Parameterized in a
!! Mesoscale NWP Model. Monthly Weather Review, doi:http://dx.doi.org/10.1175/MWR-D-11-00352.1
!!
!! Anna C. Fitch, National Center for Atmospheric Research (formerly University of Bergen)
!!
!!
!! History of changes:
!!
!! WRFV3.5.1:
!! WRFV3.6:   Modified by Pedro A. Jimenez to include:
!!             - Initialize the wind turbines in this module.
!!             - Introduce z_at_walls to avoid instabilities due to neglecting
!!                the perturbation of the geopotential height.
!!             - User friendly interface to introduce the technical characteritics of
!!                the wind turbines.
!!             - Only uses one set of turbine coefficients using the wind speed at hub height
!!             - Two standing coefficients.
!!             - Calculates the power produced by the wind turbines.
!! BRAMS 5.2: - modified by Luiz Flavio to adapt in BRAMS Model
!!
!!@par References:
!!
!!@cite Fitch, A. C. et al. 2012: Local and Mesoscale Impacts of Wind Farms as Parameterized in a
!!    Mesoscale NWP Model. Monthly Weather Review, doi:http://dx.doi.org/10.1175/MWR-D-11-00352.1
!!@cite Fitch, A. C. et al. 2013: Mesoscale Influences of Wind Farms Throughout a Diurnal Cycle.
!!    Monthly Weather Review, doi:http://dx.doi.org/10.1175/MWR-D-12-00185.1
!!@cite Fitch, A. C. et al. 2013: Parameterization of Wind Farms in Climate Models.
!!    Journal of Climate, doi:http://dx.doi.org/10.1175/JCLI-D-12-00376.1
!!@cite Jimenez, P.A., J. Navarro, A.M. Palomares and J. Dudhia:  Mesoscale modeling of offshore wind turbines
!!    wakes at the wind farm resolving scale: a composite-based analysis with the WRF model over Horns Rev.
!!    Wind Energy, (In Press.).
!!
!! @version 1.0
!!
!!@warning This routines needs a windturbines.txt with turbines data and turbine infos.
MODULE module_wind_fitch

  IMPLICIT NONE

  integer :: max_domains
  real :: piconst
  INTEGER, PARAMETER :: MAXVALS  = 100
  INTEGER, PARAMETER :: MAXVALS2 = 100
!
  INTEGER           :: nt
  INTEGER, DIMENSION(:), ALLOCATABLE :: NKIND
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: ival,jval
  real, DIMENSION(:), ALLOCATABLE :: lat,lon
  character(len=32), ALLOCATABLE, dimension(:) :: wfname
  character(len=3), ALLOCATABLE, dimension(:) :: turbine_id


  REAL, DIMENSION(:), ALLOCATABLE :: hubheight,diameter,stc,stc2,cutin,cutout,npower
!
  REAL :: turbws(maxvals,maxvals2),turbtc(maxvals,maxvals2),turbpw(maxvals,maxvals2)

  integer, parameter :: windfarm_ij=0 !lfr substituindo grid_config_rec_type
!
CONTAINS

  SUBROUTINE  dragforce(&         ! 01
       id,                      & ! 02
       z_at_w,                  & ! 03
       u,                       & ! 04
       v,                       & ! 05
       dx,                      & ! 06
       dy,                      & !
       dz,                      & !
       dt,                      & !
       qke,                     & !
       du,                      & !
       dv,                      & !
       windfarm_opt,            & !
       power,                   & !
       ids,ide,jds,jde,kds,kde, &
       ims,ime,jms,jme,kms,kme, &
       its,ite,jts,jte,kts,kte, &
       mynum,istp,time, &
       wind_speed,u_speed,v_speed &
       )
!
!
!
  INTEGER, INTENT(IN) :: id,windfarm_opt,mynum,istp
  INTEGER, INTENT(IN) :: its,ite,jts,jte,kts,kte
  INTEGER, INTENT(IN) :: ims,ime,jms,jme,kms,kme
  INTEGER, INTENT(IN) :: ids,ide,jds,jde,kds,kde
  REAL, DIMENSION(ide,jde), INTENT(IN) :: dx,dy
  REAL, INTENT(IN) :: dt,time
  REAL, DIMENSION(kde,ide,jde), INTENT(IN) :: dz,u,v,z_at_w
  REAL, DIMENSION(kde,ide,jde), INTENT(INOUT) :: du,dv,qke

  REAL, DIMENSION(nt), INTENT(OUT) :: power
  REAL, DIMENSION(nt), INTENT(OUT) :: wind_speed
  REAL, DIMENSION(nt), INTENT(OUT) :: u_speed
  REAL, DIMENSION(nt), INTENT(OUT) :: v_speed
!
! Local
!
  REAL     blade_l_point,blade_u_point,zheightl,zheightu,z1,z2,tarea
  REAL     speed,tkecof,powcof,thrcof,wfdensity
  INTEGER  itf,jtf,ktf
  INTEGER  i,j,k,n
  INTEGER  k_turbine_bot, k_turbine_top,fnum

  real, dimension(nt) :: wind_dir

  LOGICAL :: kfound
!
! ... PAJ: more variables ...
!
  REAL :: speedhub,speed1,speed2
  real :: power1,power2,area,ec
  INTEGER :: kbot,ktop,kt

  itf=MIN0(ite,ide-1)
  jtf=MIN0(jte,jde-1)
  ktf=MIN0(kte,kde-1)


    power=0.

    DO kt = 1,nt
      IF ( windfarm_opt .eq. 1 ) THEN
!
! vertical layers cut by turbine blades
!
        k_turbine_bot=0      !bottom level
        k_turbine_top=-1     !top level
        i = ival(kt,id)
        j = jval(kt,id)
!	print *,'ktloop: ', mynum,kt,i,j,(( its .LE. i .AND. i .LE. itf ) .AND. ( jts .LE. j .AND. j .LE. jtf ))
!
         if (i>0 .and. j>0) then
            IF (( its .LE. i .AND. i .LE. itf ) .AND. &
            ( jts .LE. j .AND. j .LE. jtf )  ) THEN
!
               blade_l_point=hubheight(kt)-diameter(kt)/2. ! height of lower blade tip above ground (m)
               blade_u_point=hubheight(kt)+diameter(kt)/2. ! height of upper blade tip above ground (m)
               wfdensity = dx(i,j)*dy(i,j)   !  per turbine, so numerator is 1
!
!write (*,fmt='("Bld: ",2(I2.2,1X),2(F10.5,1X))') mynum,kt,blade_l_point,blade_u_point

               kfound = .false.
               zheightl=0.0
               ! find vertical levels cut by turbine blades
               DO k=kts,ktf
                  IF(.NOT. kfound) THEN
                     zheightu = zheightl + dz(k,i,j) ! increment height
!write (*,fmt='("Lyr: ",3(I2.2,1X),2(F10.5,1X))') mynum,kt,k,zheightl,zheightu,dz(k,i,j)
                     IF(blade_l_point .GE. zheightl .AND. blade_l_point .LE. zheightu) THEN
                        k_turbine_bot=k ! lower blade tip cuts this level
                     ENDIF
                     IF(blade_u_point .GE. zheightl .AND. blade_u_point .LE. zheightu) THEN
                        k_turbine_top=k ! upper blade tip cuts this level
                        kfound = .TRUE.
                     ENDIF
                     zheightl = zheightu
                  ENDIF
               ENDDO
!print *,'LFR->kfound: ',mynum,kt,kfound
               IF ( kfound ) THEN
                  !
                  ! ... PAJ: Changes introduced to compute only one set of turbine coefficients ...
                  !          First computes the wind speed at the hub height.
                  !
                  kfound = .false.
                  zheightl=0.
                  ! find vertical levels (half levels) within the hub height
                  DO k=kts,ktf
                     IF(.NOT. kfound) THEN
                        z2 = zheightl + 0.5*dz(k,i,j)
!
                        IF(hubheight(kt) .GE. z2 ) THEN
                           kbot=k
                        ELSE
                           ktop=k
                           kfound = .TRUE.
                        ENDIF
!
                        if (.NOT. kfound) z1=z2
                        zheightl = z2 + 0.5*dz(k,i,j)
                     ENDIF
                   ENDDO
!
                   speed1=0.
                   speed2=0.
                   if (ktop.eq.1) then
                      speedhub=sqrt(u(1,i,j)**2.+v(1,i,j)**2.)*hubheight(kt)/z1
                      u_speed(kt)=u(1,i,j)
                      v_speed(kt)=v(1,i,j)
                   else
                      speed1=sqrt(u(kbot,i,j)**2.+v(kbot,i,j)**2.)
                      speed2=sqrt(u(ktop,i,j)**2.+v(ktop,i,j)**2.)
                      speedhub=speed1+((speed2-speed1)/(z2-z1))*(hubheight(kt)-z1)
                      u_speed(kt)=(u(kbot,i,j)+u(ktop,i,j))*.5
                      v_speed(kt)=(v(kbot,i,j)+v(ktop,i,j))*.5
		      !if(i==2 .and. j==9) print *,'LFR->2',kbot,ktop,z1,z2,u(kbot,i,j),u(ktop,i,j),v(kbot,i,j),v(ktop,i,j); call flush(6)
                   endif
                   wind_speed(kt)=speedhub

!write (*,fmt='("Hub: ",3(I2.2,1X),8(F10.5,1X))') mynum,ktop,kbot,u(kbot,i,j),v(kbot,i,j),u(ktop,i,j),v(ktop,i,j),hubheight(kt),u(1,i,j),v(1,i,j),z1,z2

!write (*,fmt='("In : ",3(I2.2,1X),7(F10.5,1X))') mynum,kt,nkind(kt),npower(kt),diameter(kt),stc(kt),stc2(kt),speedhub,cutin(kt),cutout(kt)
                   !
                   ! ... calculate TKE, power and thrust coeffs
                   CALL dragcof(tkecof,powcof,thrcof,               &
                           speedhub,cutin(kt),cutout(kt),   &
                           npower(kt),diameter(kt),stc(kt),stc2(kt),nkind(kt))
                   !
                   ! ... PAJ: Computation of power generated by the wind turbine ...
                   area=piconst/4.*diameter(kt)**2.          ! area swept by turbine blades
                   power1=0.5*1.23*speedhub**3.*area*powcof
                   power(kt)=power1+power(kt)
                   power2=0.
!write (*,fmt='("Out: ",2(I2.2,1X),3(F20.5,1X))') mynum,kt,area,powcof,power(i,j); call flush(6)
!
                   DO k=k_turbine_bot,k_turbine_top ! loop over turbine blade levels
                      z1=z_at_w(k,i,j)-blade_l_point-z_at_w(1,i,j)  ! distance between k level and lower blade tip
                      z2=z_at_w(k+1,i,j)-blade_l_point-z_at_w(1,i,j) ! distance between k+1 level and lower blade tip
                      IF(z1 .LT. 0.) z1=0.0 ! k level lower than lower blade tip
                      IF(z2 .GT. diameter(kt)) z2=diameter(kt) ! k+1 level higher than turbine upper blade tip
                      CALL turbine_area(z1,z2,diameter(kt),wfdensity,tarea)
!
                      speed=sqrt(u(k,i,j)**2.+v(k,i,j)**2.)
                      power2=power2+0.5*powcof*1.23*(speed**3.)*tarea/wfdensity
                   ENDDO
                   !
                   ! ... PAJ: Computes the tendencies of TKE and momentum ...
                   !
                   DO k=k_turbine_bot,k_turbine_top ! loop over turbine blade levels
                      z1=z_at_w(k,i,j)-blade_l_point-z_at_w(1,i,j)  ! distance between k lev and lower blade tip
                      z2=z_at_w(k+1,i,j)-blade_l_point-z_at_w(1,i,j) !distance between k+1 lev and lower blade tip
                      IF(z1 .LT. 0.) z1=0.0 ! k level lower than lower blade tip
                      IF(z2 .GT. diameter(kt)) z2=diameter(kt) ! k+1 level higher than turbine upper blade tip
!
                      CALL turbine_area(z1,z2,diameter(kt),wfdensity,tarea)
!
                      speed=sqrt(u(k,i,j)**2.+v(k,i,j)**2.)
                      !`
                      ! ... PAJ: normalization introduced to conserve energy ...
                      !
                      if (power1.eq.0.or.power2.eq.0) then
                         ec=1.
                      else
                         ec=power1/power2
                      endif
!
                      ! output TKE
                      qke(k,i,j) = qke(k,i,j)+speed**3.*tarea*tkecof*dt/dz(k,i,j)*ec
                      ! output u tendency
                      du(k,i,j) = du(k,i,j)-.5*u(k,i,j)*thrcof*speed*tarea/dz(k,i,j)*ec
                      ! output v tendency
                      dv(k,i,j) = dv(k,i,j)-.5*v(k,i,j)*thrcof*speed*tarea/dz(k,i,j)*ec
                   ENDDO
                 ENDIF
               ENDIF
             endif
        ENDIF
    ENDDO

  END SUBROUTINE dragforce

! This subroutine calculates area of turbine between two vertical levels
! Input variables :
!            z1 = distance between k level and lower blade tip
!            z2 = distance between k+1 level and lower blade tip
!            wfdensity = wind farm density in m^-2
!     tdiameter = turbine diameter
! Output variable :
!         tarea = area of turbine between two levels * wfdensity
  SUBROUTINE turbine_area(z1,z2,tdiameter,wfdensity,tarea)

  REAL, INTENT(IN) ::tdiameter,wfdensity
  REAL, INTENT(INOUT) ::z1,z2
  REAL, INTENT(OUT):: tarea
  REAL r,zc1,zc2

  r=tdiameter/2.              !r = turbine radius
  z1=r-z1                   !distance of kth level from turbine center
  z2=r-z2                   !distance of k+1 th level from turbine center
  zc1=abs(z1)
  zc2=abs(z2)
  !turbine area between z1 and z2
  IF(z1 .GT. 0. .AND. z2 .GT. 0.) THEN
     tarea=zc1*sqrt(r*r-zc1*zc1)+r*r*asin(zc1/r)- &
     (zc2*sqrt(r*r-zc2*zc2)+r*r*asin(zc2/r))
  ELSE IF(z1 .LT. 0. .AND. z2 .LT. 0.) THEN
     tarea=zc2*sqrt(r*r-zc2*zc2)+r*r*asin(zc2/r)- &
     (zc1*sqrt(r*r-zc1*zc1)+r*r*asin(zc1/r))
  ELSE
     tarea=zc2*sqrt(r*r-zc2*zc2)+r*r*asin(zc2/r)+ &
     zc1*sqrt(r*r-zc1*zc1)+r*r*asin(zc1/r)
  ENDIF
  tarea=tarea*wfdensity      !turbine area * wind farm density

  END SUBROUTINE turbine_area


  SUBROUTINE dragcof(tkecof,powcof,thrcof,speed,cispeed,cospeed, &
                     tpower,tdiameter,stdthrcoef,stdthrcoef2,nkind)


  REAL, INTENT(IN):: speed, cispeed, cospeed, tpower,tdiameter,stdthrcoef,stdthrcoef2
  REAL, INTENT(OUT):: tkecof,powcof,thrcof
  REAL :: power,area,mspeed,hspeed
!
! ... PAJ ...
!
   INTEGER :: nkind,k,nu,nb
   LOGICAL :: vfound
   REAL :: fac1,fac2

  area=piconst/4.*tdiameter**2.          ! area swept by turbine blades

      vfound=.false.
      DO k=1,maxvals2
            IF(.NOT. vfound) THEN
              IF(turbws(nkind,k).GT.speed) THEN
                nu=k
                nb=k-1
                vfound=.true.
              ENDIF
            ENDIF
      ENDDO
!
  IF (speed .LE. cispeed) THEN
     thrcof = stdthrcoef
  ELSE
    IF (speed .GE. cospeed) THEN
     thrcof = stdthrcoef2
     ELSE
     thrcof = turbtc(nkind,nb)+(turbtc(nkind,nu)-turbtc(nkind,nb))/(turbws(nkind,nu)-turbws(nkind,nb))*(speed-turbws(nkind,nb))
    ENDIF
  ENDIF
!
! ... power coeficient ...
!
  IF(speed .LE. cispeed .OR. speed .GE. cospeed) THEN
     power=0.
     powcof=0.
  ELSE
      fac1=1000./(0.5*1.23*turbws(nkind,nb)**3.*area)
      fac2=1000./(0.5*1.23*turbws(nkind,nu)**3.*area)
!print *, 'LFR-> inside dragcof: ',fac1,turbws(nkind,nb),fac2,turbws(nkind,nu),area,piconst,tdiameter
      power = turbpw(nkind,nb)+(turbpw(nkind,nu)-turbpw(nkind,nb))/(turbws(nkind,nu)-turbws(nkind,nb)) &
                               *(speed-turbws(nkind,nb))
!print *, 'LFR-> inside dragcof: ',turbpw(nkind,nb),fac1,turbpw(nkind,nu),fac2,turbws(nkind,nu),turbws(nkind,nb),speed
      powcof = turbpw(nkind,nb)*fac1+(turbpw(nkind,nu)*fac2-turbpw(nkind,nb)*fac1)/(turbws(nkind,nu)-turbws(nkind,nb)) &
                                     *(speed-turbws(nkind,nb))
  ENDIF
!
  ! tke coefficient calculation

  tkecof=thrcof-powcof
  IF(tkecof .LT. 0.) tkecof=0.
!
  END SUBROUTINE dragcof
!
  SUBROUTINE init_module_wind_fitch(xlong,xlat,windfarm_initialized,&
                                            m2,m3,ia,iz,ja,jz, &
                                            mynum,master_num,ngrids,nmachs,ngrid)
!
  INCLUDE 'mpif.h'
!
   integer,intent(in) :: m2,m3,ia,iz,ja,jz
   integer,intent(in) :: mynum,master_num,ngrids
   integer,intent(in) :: nmachs,ngrid
   REAL,DIMENSION(m2,m3) , INTENT(IN) :: xlong,xlat

   logical :: windfarm_initialized
!
   CHARACTER*256 :: num,input,message_wind,chead
   real :: ts_rx,ts_ry
   REAL :: known_lat, known_lon
   INTEGER :: i,j,nval,k,ierr,proc,ihead
   integer :: posi,posj

   IF (mynum==master_num) THEN
!
! ... PAJ: Opens the file with the location of the wind turbines ...
!
        if ( windfarm_ij .eq. 1 ) then
          !open(70,file='windturbines-ij.txt',form='formatted',status='old')
	        stop 'No prepared for i,j turbines position'
        else
          open(70,file='./tables/windturbines/windturbines.txt',form='formatted',status='old')
        end if
!
! ... PAJ: Counts the turbines ...
!
       nt=0
 10    read(70,*,end=100)
       nt=nt+1
       goto 10
!
 100   continue
       rewind (70)
     endif
!
    !Sending number of turbines for all processors
     call MPI_BCAST(nt,1, MPI_INTEGER, master_num, MPI_COMM_WORLD, ierr)
!
! ... PAJ: Initializes the configuration of the wind farm(s) ...
!
     if (.not. windfarm_initialized) then
       max_domains=ngrids
       allocate (nkind(nt),ival(nt,max_domains),jval(nt,max_domains))
       allocate (hubheight(nt),stc(nt),stc2(nt),cutin(nt),cutout(nt),diameter(nt),npower(nt))
       allocate (lat(nt),lon(nt))
       allocate (wfname(nt))
       allocate (turbine_id(nt))
       windfarm_initialized=.true.
     endif

    !Readind position of all turbines
    if (mynum==master_num) then
     write (*,fmt='(A)') ' --- Windturbine information --- '
     write (*,fmt='(A3,1X,2(A6,1X),A4,1X,A,1X,A)') &
	'#tn','lat','lon','Type','File Name','Turbine ID'

     do k=1,nt
       if ( windfarm_ij .eq. 1 ) then
         !read(70,*) ival(k,id), jval(k,id), nkind(k)
         stop 'No prepared for i,j turbines position'
       else
         read(70,*)lat(k),lon(k),nkind(k),wfname(k),turbine_id(k)
	       write (*,fmt='(I3.3,1X,2(F8.2,1X),I4.4,1X,A,1X,A,1X,A)') k,lat(k),lon(k),nkind(k),wfname(k),turbine_id(k)
       endif
      enddo
      close(70)
    endif

    ! sending the position and type of all turbines to all processors
    CALL MPI_BCAST(lat,nt, MPI_REAL, master_num, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST(lon,nt, MPI_REAL, master_num, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST(nkind,nt, MPI_INTEGER, master_num, MPI_COMM_WORLD, ierr)
    DO k=1,nt
      call MPI_BCAST(wfname(k),32, MPI_CHARACTER, master_num, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(turbine_id(k),3, MPI_CHARACTER, master_num, MPI_COMM_WORLD, ierr)
    enddo

    !Searching for position in x,y in local processor
    do k=1,nt
        call LatLon2xy(lat(k),lon(k),xlat,xlong,ival(k,ngrid),jval(k,ngrid),m2,m3,ia,iz,ja,jz,mynum)

!         if(ival(k,ngrid)>0 .and. jval(k,ngrid)>0) &
         !write (*,fmt='(A,I3.3,A,F6.2,A,F6.2,A,I4,A,I4,A,I1)') &
 	!"WINDFARM Turbine #",k," | Lat,Lon =",lat(k),",",lon(k)," | (i,j) = (",ival(k,ngrid),",",jval(k,ngrid),"), | Type = ",nkind(k)
 	call flush(6)
    end do

    if (mynum==master_num) then
!
! ... PAJ: Read the tables with the turbine's characteristics ...
!
         turbws=0.
         turbtc=0.
         turbpw=0.
         DO i=1,nt
          write(num,*) nkind(i)
          num=adjustl(num)
          input="./tables/windturbines/wind-turbine-"//trim(num)//".tbl"
          OPEN(file=TRIM(input),unit=19,FORM='FORMATTED',STATUS='OLD')
          do ihead=1,9
             READ (19,*,ERR=132) chead
          end do
          READ (19,*,ERR=132)nval
          READ(19,*,ERR=132)hubheight(i),diameter(i),stc(i),npower(i)
            DO k=1,nval
              READ(19,*,ERR=132)turbws(nkind(i),k),turbtc(nkind(i),k),turbpw(nkind(i),k)
            ENDDO
          cutin(i)  = turbws(nkind(i),1)
          cutout(i) = turbws(nkind(i),nval)
          stc2(i) = turbtc(nkind(i),nval)
	  !print *, 'LFR->stc: ',i,nval,nkind(i),turbtc(nkind(i),nval)
          close (19)
         ENDDO

 132   continue
!
    endif

    !Sending turbine's characteristics to all processors
    call MPI_BCAST(hubheight,nt, MPI_REAL, master_num, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(diameter,nt, MPI_REAL, master_num, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(stc,nt, MPI_REAL, master_num, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(stc2,nt, MPI_REAL, master_num, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(npower,nt, MPI_REAL, master_num, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(cutin,nt, MPI_REAL, master_num, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(cutout,nt, MPI_REAL, master_num, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(turbws,maxvals*maxvals2, MPI_REAL, master_num, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(turbtc,maxvals*maxvals2, MPI_REAL, master_num, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(turbpw,maxvals*maxvals2, MPI_REAL, master_num, MPI_COMM_WORLD, ierr)


  end subroutine init_module_wind_fitch


  subroutine LatLon2xy(lat,lon,latin,lonin,iPos,jPos,m2,m3,ia,iz,ja,jz,mchnum)
    !Find x and y pos from a especific lat and lon inside the model tile

    integer, intent(in) :: m2,m3,mchnum !max size #points x and y of tile
    integer, intent(in) :: ia,iz,ja,jz
    !#points without border
    real,intent(in) :: lat ! in degrees
    real,intent(in) :: lon ! lon in degrees
    real,intent(in) :: latin(m2,m3) ! array of lat's
    real,intent(in) :: lonin(m2,m3) ! array of lon's
    integer,intent(out) :: iPos !i Position of point in lat,lon
    integer,intent(out) :: jPos !j Position of point in lat,lon

    integer :: i,j
    iPos=-1
    jPos=-1

    do i=1,iz-1
      do j=1,jz-1
	!Compare the position with lat and lon of each point. If found fill the Ipos, Jpos
        if((lonin(i,j)<=lon .and. lonin(i+1,j)>lon) .and. &
	      (latin(i,j)<=lat .and. latin(i,j+1)>lat)) then
	      !write (80+mchnum,fmt='(6(F6.2,1X),2(I2.2,1X))') lonin(i,j),lon,lonin(i+1,j),latin(i,j),lat,latin(i,j+1),i+1,j+1; call flush(80+mchnum)
           iPos=i+1
           jPos=j+1
        endif
      enddo
    enddo

   end subroutine LatLon2xy

END MODULE module_wind_fitch
