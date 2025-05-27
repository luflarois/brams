!>
!!@package BRAMS-WindFarm
!!@author Luiz Flavio
!!@date 05/22/2017
!!@brief Interface to call windfarms routines and plots the results at BRAMS' end.
!!@copyright Under CC-GPL license
!!@link https://creativecommons.org/licenses/GPL/2.0/legalcode.pt
!!@param ngrid,m1,m2,m3,ia,iz,ja,jz
!!@version 1.0
!!
!!@warning This routines needs a windturbines.txt with turbines data.
MODULE wind_Farm
  
  use ModNamelistFile, only: namelistFile
  use io_params, only :           & !INTENT(IN)
       frqanl, &
       afilout

  implicit none
  include "files.h"

  integer :: windfarm
  logical :: windfarm_initialized=.false.
  real, allocatable, dimension(:,:,:) :: power
  real, allocatable, dimension(:,:,:) :: wind_speed
  real, allocatable, dimension(:,:,:) :: u_speed
  real, allocatable, dimension(:,:,:) :: v_speed
  integer :: nTicsTime
  logical, allocatable :: firstwf(:)
  character(len=f_name_length) :: turbFilename
  character(len=f_name_length) :: wfFile
  integer :: timeInterval_wf

CONTAINS


subroutine wind_farm_driver(ngrid,m1,m2,m3,ia,iz,ja,jz)

  use mem_grid, only     : ngrids,grid_g,dzt,zm,istp,time,deltaxn,dtlongn,timmax, &
                           idate1,imonth1,iyear1,ihour1,itime1
  use mem_basic, only: basic_g
  use mem_turb, only: turb_g
  use mem_tend, only: tend
  use rconstants, only: pi
  use node_mod, only: master_num,mchnum,nmachs
  use module_wind_fitch, only: dragforce,init_module_wind_fitch,piconst,nt,ival,jval,lat,lon, &
                               hubheight,diameter,stc,stc2,cutin,cutout,npower,wfname,turbine_id


  integer, intent(in) :: ngrid,m1,m2,m3,ia,iz,ja,jz

  real, dimension(m1,m2,m3) :: dz       !dz between full levels (m)
  real, dimension(m1,m2,m3) :: z_at_w   !height above sea level at interface (m)
  REAL, DIMENSION(m1,m2,m3) :: du,dv,qke

  integer, parameter :: lid=1
  integer, parameter :: windfarm_opt=1
  logical, parameter :: windfarm_dump=.true. ! T - Write a file with turbine output
  character(len=3) :: cturb
  real :: latlon(2)

  integer :: i,j,k,nk,fnum,ntt
  
  if(windfarm==0) return
  
  piconst = pi

  if(.not. windfarm_initialized) then
    nTicsTime=int(timmax/dtlongn(ngrid))
    !if (mchnum==0) print *,'--- Initializing Windfarm ---'; call flush(6)
    call init_module_wind_fitch(       &
        grid_g(ngrid)%glon(:,:),       &
        grid_g(ngrid)%glat(:,:),       &
        windfarm_initialized,          &
        m2,m3,ia,iz,ja,jz,             &
        mchnum,master_num,ngrids,nmachs, &
        ngrid)

        allocate(power(nt,nTicsTime,ngrids),wind_speed(nt,nTicsTime,ngrids))
        allocate(u_speed(nt,nTicsTime,ngrids),v_speed(nt,nTicsTime,ngrids))
				power=0.0; wind_speed=0.0; u_speed=0.0; v_speed=0.0


	      if (windfarm_dump) then
	        allocate(firstwf(nt))
	        firstwf=.true.
	      endif

        print *,"LFR windfarm_initialized=",windfarm_initialized

  endif


  do j=ja,jz
     do i=ia,iz
       do k=1,m1
          dz(k,i,j) = 1./dzt(k)
          z_at_w(k,i,j) = grid_g(ngrid)%rtgt(i,j)*zm(k) + grid_g(ngrid)%topt(i,j)
       enddo
    enddo
  enddo

  du=0.0;dv=0.0;qke=0.0

  call dragforce(                 & ! 01
       lid,                       & ! 02
       z_at_w,  & ! 03
       basic_g(ngrid)%up, & ! 04
       basic_g(ngrid)%vp, & ! 05
       grid_g(ngrid)%dxt, &
       grid_g(ngrid)%dyt, &
       dz,      & ! 07
       dtlongn(ngrid),            & ! 08
       turb_g(ngrid)%tkep, & ! 09
       tend%ut               , & ! 10
       tend%vt,                & ! 11
       windfarm_opt,              & ! 12
       power(:,istp,ngrid),       & ! 13
       1,m2,1,m3,1,m1,            & !
       ia,iz,ja,jz,1,m1,          & !
       1,m2,1,m3,1,m1,            & !
       mchnum,istp,time,          & !
       wind_speed (:,istp,ngrid), & !
       u_speed(:,istp,ngrid),     & !
       v_speed(:,istp,ngrid)      &
       )

  end subroutine wind_farm_driver

  subroutine StoreNamelistFileAtWindFarm(oneNamelistFile)
    type(namelistFile), pointer :: oneNamelistFile

    windfarm = oneNamelistFile%windfarm
    wfFile=oneNamelistFile%wfFile
    timeInterval_wf = oneNamelistFile%timeInterval_wf

  end subroutine StoreNamelistFileAtWindfarm


  subroutine output_windFarms()
    use mem_grid, only     : ngrids,grid_g,dzt,zm,istp,time,deltaxn,dtlongn,timmax, &
                             idate1,imonth1,iyear1,ihour1,itime1
    use mem_basic, only: basic_g
    use mem_turb, only: turb_g
    use mem_tend, only: tend
    use rconstants, only: pi
    use node_mod, only: master_num,mchnum,nmachs
    use module_wind_fitch, only: dragforce,init_module_wind_fitch,piconst,nt,ival,jval,lat,lon, &
                                 hubheight,diameter,stc,stc2,cutin,cutout,npower,wfname,turbine_id

    integer :: k
    character(len=3) :: cturb
    real :: latlon(2)
    integer,parameter :: ngrid=1

    !--- check if it is time to write the hdf turbine file
    if((mod(time,real(timeInterval_wf)) < 1.0 .or. time==timmax) .and. time>1.0) then
      !print *, 'LFR Dentro da escrita',timeInterval_wf
      do k=1,nt
        if(ival(k,ngrid)<0 .or. jval(k,ngrid)<0) cycle
        write(cturb,fmt='(I3.3)') k
        latlon(1)=grid_g(ngrid)%glat(ival(k,ngrid),jval(k,ngrid))
        latlon(2)=grid_g(ngrid)%glon(ival(k,ngrid),jval(k,ngrid))
        !write (80+mchnum,*) trim(wfname(k))
        !print *, "LFR END TURBINES wfname(k)=",trim(wfname(k)),turbine_id(k)
        call writeTurbines(power(k,:,ngrid),wind_speed(k,:,ngrid), &
                           u_speed(k,:,ngrid),v_speed(k,:,ngrid),wfname(k),&
                           cturb,nTicsTime,latlon,lat(k),lon(k),idate1,imonth1, &
                             iyear1,ihour1,itime1,time,dtlongn(ngrid),piconst, &
                             ngrid,wffile,frqanl,istp,turbine_id(k))

      end do
    end if

  end subroutine output_windFarms

SUBROUTINE writeTurbines(power,wind_speed,u_speed,v_speed,windfarm, &
                           cTurbine,ixDim,LatLon,lat,lon,day,month, &
                           year,hour,itime,time,deltaT,piconst,ngrid,afilout,&
													 frqanl,istp,tid)
      !Subroutine for write power & others
      !  The size of filename+fileExtension must be less than 253 characters
      ! Author: Luiz Flavio
			use ModDateUtils

      IMPLICIT NONE

      character(len=*),parameter :: wfmt='(I4.4,"-",I2.2,"-",I2.2,"T",I2.2,":",I2.2,":",I2.2,",",F10.1,",",3(F5.1,","),F5.1)'

      INTEGER, INTENT(IN) :: ixdim,ngrid !Dims of variable
      REAL, INTENT(IN) :: power(ixdim) !
      REAL, INTENT(IN) :: wind_speed(ixdim) !
      REAL, INTENT(IN) :: u_speed(ixdim) !
      REAL, INTENT(IN) :: v_speed(ixdim) !

      real, intent(in) :: latlon(2)
      real, intent(in) :: lat,lon,piconst

      integer,intent(in) :: day,month,year,hour,itime,istp
      real,intent(in) :: time,frqanl
			CHARACTER(LEN=*), INTENT(IN) :: windfarm !Filename to be composed
      CHARACTER(LEN=*), INTENT(IN) :: cTurbine !Filenamen extension to be added
      character(len=*), intent(in) :: afilout
      character(len=*), intent(in) :: tid

      CHARACTER(LEN=256) :: filename ! File name

      INTEGER     ::   error ! Error flag
      INTEGER     ::   i,j,k,i1,i2
			integer     :: ticsPerFrq,outyear,outmonth,outdate,outhour,inyear,inmonth,indate,inhour
			integer     :: ipos,begdata,enddata,ih,im,is
			real :: deltaT,dv
			character :: cgrid

			!ticsPerFrq=frqanl/deltaT
      ticsPerFrq=timeInterval_wf/deltaT
			begdata=istp-ticsPerFrq+1
			enddata=istp

			CALL date_add_to (year,month,day,hour,  & !Determine the current model date
       time,'s',inyear,inmonth,indate,inhour)

			CALL date_add_to (inyear,inmonth,indate,inhour,  & !Determine the initial date for this write
       (real(timeInterval_wf))*(-1),'s',outyear,outmonth,outdate,outhour)

      !coxilha_negra-T-2022-08-01-040000-072.dat
      !
       write(filename,fmt='(A,"/",A,"-",A,"-",2(i4.4,I2.2,i2.2,"T",i6.6,"."),"csv")') trim(afilout),trim(windfarm),tid &
            ,outyear,outmonth,outdate,outhour,inyear,inmonth,indate,inhour

      !write(cgrid,fmt='(I1.1)') ngrid
      !CALL makefnam(fileName,afilout, time, &
      !      year, month, day, itime*100,'T', cTurbine, 'dat')
      !print *,"filename=",trim(filename)
			call rams_f_open(25,trim(filename),'FORMATTED','REPLACE','WRITE',1)

			!write (25,fmt='(A,1X,A)') windfarm,cTurbine
			!write (25,fmt='(i4.4,1X,I2.2,1X,i2.2,1X,i6.6)') inyear,inmonth,indate,inhour
			!write (25,fmt='(i4.4,1X,I2.2,1X,i2.2,1X,i6.6)') outyear,outmonth,outdate,outhour
			!write (25,*) frqanl,deltaT,lat,lon
			!write (25,*) enddata-begdata+1
      write(25,*) "time,power[W],wind_speed[m/s],u_speed[m/s],v_speed[m/s],direction[gr]"
			do i=begdata,enddata
        CALL date_add_to (outyear,outmonth,outdate,outhour,  & !Determine the current model date
        real(i-1)*deltaT,'s',inyear,inmonth,indate,inhour)
				dv=(180.0*atan2(u_speed(i),v_speed(i)))/piconst
        if(dv<0.0) dv=360.0+dv
        !2023-06-05T00:00:00
        
        ih = int(inhour/10000.)
        im = int((inhour-ih*10000)/100)
        is = int(mod(inhour-ih*10000,100))
        !wfmt='(I4.4,"-",I2.2,"-",I2.2,"T",I2.2,":",I2.2,":",I2.2,",",F10.1,",",3(5.1,","),5.1)'
				write(25,fmt=wfmt) inyear,inmonth,indate,ih,im,is,power(i),wind_speed(i),u_speed(i),v_speed(i),dv
			end do

			close(25)

end subroutine writeTurbines


END MODULE wind_Farm
