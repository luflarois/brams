!=============================================================================================
module filesMod
   !# Module to manipulate namelist file
   !#
   !# @note
   !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
   !#
   !# **Brief**: manipulate namelist file
   !#
   !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
   !#
   !# **Author(s)**: Luiz Flavio Rodrigues **&#9993;**<luiz.rodrigues@inpe.br>
   !#
   !# **Date**: 26 August 2020 (Wednesday)
   !# @endnote
   !#
   !# @changes
   !# &#9744; <br/>
   !# @endchanges
   !# @bug
   !#
   !#@endbug
   !#
   !#@todo
   !#  &#9744; <br/>
   !# @endtodo
   !#
   !# @warning
   !# Now is under CC-GPL License, please see
   !# &copy; <https://creativecommons.org/licenses/GPL/2.0/legalcode.pt>
   !# @endwarning
   !#

   !Use area
   use dump

   implicit none

   include "constants.f90"
   character(len=*),parameter :: sourceName='filesMod.f90' !Name of source code
   character(len=*),parameter :: procedureName='**namelist**' !Name of this procedure
   !
   !Local Parameters

   !Local variables
   integer :: chemjIni,chemjFin,chemnpj
   integer :: chemiIni,chemiFin,chemnpi

contains

   !=============================================================================================
   subroutine readNamelist(namelistFile)
      !# read namelist file from input
      !#
      !# @note
      !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
      !#
      !# **Brief**: Read the namelist file, pre.nml and return data
      !#
      !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
      !#
      !# **Author(s)**: Luiz Flavio Rodrigues **&#9993;**<luiz.rodrigues@inpe.br>
      !#
      !# **Date**: 26 August 2020 (Wednesday)
      !# @endnote
      !#
      !# @changes
      !# &#9744; <br/>
      !# @endchanges
      !# @bug
      !#
      !#@endbug
      !#
      !#@todo
      !#  &#9744; <br/>
      !# @endtodo
      !#
      !# @warning
      !# Now is under CC-GPL License, please see
      !# &copy; <https://creativecommons.org/licenses/GPL/2.0/legalcode.pt>
      !# @endwarning
      !#

      !Use area
      use dump
      use memoryMod, only: &
         init_year ,  init_month,  init_day,  init_hour,  &
         final_year, final_month, final_day, final_hour,  &
         step,                                            &
         atmos_type, atmos_prefix,atmos_sufix,atmos_idir, &
         levels, initial_latitude,final_latitude,         &
         initial_longitude,final_longitude,               &
         chem_type, chem_idir, &
         chem1_prefix,chem1_sufix,                        &
         out_type, out_prefix, out_sufix, out_dir

      use utilsMod, only: &
         getUnit, &
         releaseUnit

      implicit none

      include "constants.f90"
      character(len=*),parameter :: procedureName='**readNamelist**' !Name of this procedure
      !
      !Local Parameters

      !Input/Output variables
      character(len=*), intent(in) :: namelistFile

      !Local variables
      integer :: lunit

      !Code
      NAMELIST/args_input/init_year ,  init_month,  init_day,  init_hour,  &
         final_year, final_month, final_day, final_hour,  &
         step,                                            &
         atmos_type, atmos_prefix,atmos_sufix,atmos_idir, &
         levels, initial_latitude,final_latitude,         &
         initial_longitude,final_longitude,               &
         chem_type, chem_idir, &
         chem1_prefix,chem1_sufix,&
         out_type, out_prefix, out_sufix, out_dir

      lunit=getUnit()
      open(unit=lunit, file=trim(namelistFile), status='old')
      read(lunit, args_input)
      ierrNumber=releaseUnit(lunit)


   end subroutine readNamelist

   !=============================================================================================
   subroutine createGrib2FilesNames(numberOfSteps)
      !# Create a set of files Names for Grib2
      !#
      !# @note
      !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
      !#
      !# **Brief**: Creates a set of filenames for GribFiles accordingly the initial date and
      !# others informations from Namelist
      !#
      !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
      !#
      !# **Author(s)**: Luiz Flavio Rodrigues **&#9993;**<luiz.rodrigues@inpe.br>
      !#
      !# **Date**: 26 August 2020 (Wednesday)
      !# @endnote
      !#
      !# @changes
      !# &#9744; <br/>
      !# @endchanges
      !# @bug
      !#
      !#@endbug
      !#
      !#@todo
      !#  &#9744; <br/>
      !# @endtodo
      !#
      !# @warning
      !# Now is under CC-GPL License, please see
      !# &copy; <https://creativecommons.org/licenses/GPL/2.0/legalcode.pt>
      !# @endwarning
      !#

      !Use area
      use dump
      use memoryMod, only: &
         atmos_prefix, &
         atmos_sufix, &
         atmos_idir, &
         step, &
         grib2FilesNames, &
         grib2InvFilesNames, &
         atmosDate, &
         atmosHoursCount

      implicit none

      include "constants.f90"
      character(len=*),parameter :: procedureName='**createGrib2FilesNames**' !Name of this procedure
      !
      !Local Parameters
      character(len=69) :: line='+--+----------------------------------------------------------------+'

      !Input/Output variables
      integer, intent(in) :: numberOfSteps

      !Local variables
      integer :: i

      if(allocated(grib2FilesNames)) iErrNumber=dumpMessage(c_tty,c_yes,sourceName,procedureName &
         ,c_fatal,' grib2FilesNames array already allocated with size=',size(grib2FilesNames),"I2.2")

      allocate(grib2FilesNames(numberOfSteps))

      if(allocated(grib2InvFilesNames)) iErrNumber=dumpMessage(c_tty,c_yes,sourceName,procedureName &
         ,c_fatal,' grib2InvFilesNames array already allocated with size=',size(grib2InvFilesNames),"I2.2")

      allocate(grib2InvFilesNames(numberOfSteps))
      allocate(atmosHoursCount(numberOfSteps))

      !Code

      do i=1,numberOfSteps
         write(grib2FilesNames(i),fmt='(A,I3.3,A)') trim(atmos_prefix),(i-1)*step,trim(atmos_sufix)
         write(grib2InvFilesNames(i),fmt='(A,I3.3,A)') trim(atmos_prefix),(i-1)*step,trim(atmos_sufix)//'.inv'
         atmosHoursCount(i)=(i-1)*step
         !print *,trim(grib2FilesNames(i))
         !print *,trim(grib2InvFilesNames(i))
      enddo

      iErrNumber=dumpMessage(c_tty,c_yes,'','',c_notice,'Grib2 file names inventory:')
      write(*,fmt='(A)') line
      do i=1,numberOfSteps
         write(*,fmt='("|",I2.2,"|",A64,"|")') i,trim(grib2FilesNames(i))
         write(*,fmt='(A)') line
      enddo

      if(allocated(atmosDate)) iErrNumber=dumpMessage(c_tty,c_yes,sourceName,procedureName &
         ,c_fatal,' atmosDate array already allocated with size=',size(atmosDate),"I2.2")
      allocate(atmosDate(numberOfSteps))


   end subroutine createGrib2FilesNames

   !=============================================================================================
   subroutine createEra5FilesNames(numberOfSteps)
      !# Create a set of files Names for era5
      !#
      !# @note
      !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
      !#
      !# **Brief**: Creates a set of filenames for era5 reanalisys files accordingly the
      !# initial date and others informations from Namelist
      !#
      !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
      !#
      !# **Author(s)**: Luiz Flavio Rodrigues **&#9993;**<luiz.rodrigues@inpe.br>
      !#
      !# **Date**: 17 February 2021 (Wednesday)
      !# @endnote
      !#
      !# @changes
      !# &#9744; <br/>
      !# @endchanges
      !# @bug
      !#
      !#@endbug
      !#
      !#@todo
      !#  &#9744; <br/>
      !# @endtodo
      !#
      !# @warning
      !# Now is under CC-GPL License, please see
      !# &copy; <https://creativecommons.org/licenses/GPL/2.0/legalcode.pt>
      !# @endwarning
      !#

      !Use area
      use dump
      use memoryMod, only: &
         atmos_prefix, &
         atmos_sufix, &
         atmos_idir, &
         step, &
         era5FilesNames, &
         atmosDate, &
         atmosHoursCount, &
         init_year, &
         init_month, &
         init_day, &
         init_hour
      use utilsMod, only: &
         date_add_to_dble

      implicit none

      include "constants.f90"
      character(len=*),parameter :: procedureName='**createEra5FilesNames**' !Name of this procedure
      !
      !Local Parameters
      character(len=69) :: line='+--+----------------------------------------------------------------+'
      character, parameter :: cTime='h'

      !Input/Output variables
      integer, intent(in) :: numberOfSteps

      !Local variables
      integer :: i,iyy,imm,idd,ihh
      integer :: thisYear
      integer :: thisHour
      integer :: thisMonth
      integer :: thisDay


      if(allocated(era5FilesNames)) iErrNumber=dumpMessage(c_tty,c_yes,sourceName,procedureName &
         ,c_fatal,' era5FilesNames array already allocated with size=',size(era5FilesNames),"I2.2")

      allocate(era5FilesNames(numberOfSteps))

      allocate(atmosHoursCount(numberOfSteps))

      !Code

      thisYear=init_year
      thisMonth=init_month
      thisDay=init_day
      thisHour=init_Hour

      do i=1,numberOfSteps
         write(era5FilesNames(i),fmt='(A,I4.4,I2.2,I2.2,I2.2,A)') trim(atmos_prefix) &
            ,thisYear,thisMonth,thisDay,thisHour,trim(atmos_sufix)
!            print *,'LFR-DBG:',thisYear,thisMonth,thisDay,thisHour*10000
         call date_add_to_dble(thisYear,thisMonth,thisDay,thisHour*10000,dble(step),cTime &
            ,iyy,imm,idd,ihh)
         thisYear=iyy
         thisMonth=imm
         thisDay=idd
         thisHour=ihh/10000
!            print *,'LFR-DBG:',iyy,imm,idd,ihh
      enddo

      iErrNumber=dumpMessage(c_tty,c_yes,'','',c_notice,'Era5 file names inventory:')
      write(*,fmt='(A)') line
      do i=1,numberOfSteps
         write(*,fmt='("|",I2.2,"|",A64,"|")') i,trim(era5FilesNames(i))
         write(*,fmt='(A)') line
      enddo

      if(allocated(atmosDate)) iErrNumber=dumpMessage(c_tty,c_yes,sourceName,procedureName &
         ,c_fatal,' atmosDate array already allocated with size=',size(atmosDate),"I2.2")
      allocate(atmosDate(numberOfSteps))

   end subroutine createEra5FilesNames



   !=============================================================================================
   subroutine readAtmosGrib2(iStep,stepsBetDates)
      !# Open a grib2 file and read the contents
      !#
      !# @note
      !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
      !#
      !# **Brief**: Open GFS Grib2 file, read the contents
      !#
      !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
      !#
      !# **Author(s)**: Luiz Flavio Rodrigues **&#9993;**<luiz.rodrigues@inpe.br>
      !#
      !# **Date**: 31 July 2020 (Friday)
      !# @endnote
      !#
      !# @changes
      !# &#9744; <br/>
      !# @endchanges
      !# @bug
      !#
      !#@endbug
      !#
      !#@todo
      !#  &#9744; <br/>
      !# @endtodo
      !#
      !# @warning
      !# Now is under CC-GPL License, please see
      !# &copy; <https://creativecommons.org/licenses/GPL/2.0/legalcode.pt>
      !# @endwarning
      !#

      !Use area
      use dump

      use wgrib2api

      use memoryMod, only: &
         atmos_idir, &
         grib2FilesNames, &
         grib2InvFilesNames, &
         initial_latitude, &
         final_latitude, &
         initial_longitude, &
         final_longitude, &
         atmosNz, &
         atmosNx, &
         atmosNy, &
         atmosLon, &
         atmosLat, &
         atmosXs, &
         atmosYs, &
         atmosNv, &
         atmosLevels, &
         atmosVarNames,&
         atmosValues  , &
         atmosDate, &
         atmosHoursCount, &
         levels

      use utilsMod, only: &
         getUnit, &
         releaseUnit, &
         fileExist, &
         date_add_to_dble


      implicit none

      include "constants.f90"
      character(len=*),parameter :: procedureName='**openatmosGrib2**'
      !# Name of this procedure
      !
      !Local Parameters
      character(len=*),parameter :: wind_u_varname='UGRD'
      !# Name of wind u var in grib2 (1)
      character(len=*),parameter :: wind_v_varname='VGRD'
      !# Name of wind v var in grib2 (2)
      character(len=*),parameter :: temperature_varname='TMP'
      !# Name of temperature var in grib2 (3)
      character(len=*),parameter :: geo_varname='HGT'
      !# Name of geopotential var in grib2 (4)
      character(len=*),parameter :: ur_varname='RH'
      !# Name of relative humidity var in grib2 (5)
      real,dimension(5),parameter :: scale_factor = (/1.,1.,1.,1.,0.01/)
      !# Factor to multiplie each var (1,2,3,4,5)
      real,dimension(5), parameter :: limi=(/-100.,-100.,180., -500.,  0./)
      !# inferior limits for each variable
      real,dimension(5), parameter :: lims=(/ 100., 100.,333.,80000.,100./)
      !# inferior limits for each variable

      integer, parameter :: maxpr=100

      !Input/Output variables
      integer, intent(in)         :: iStep
      !# timestep
      integer, intent(in)         :: stepsBetDates
      !# total of timesteps

      !Local variables
      integer :: nprz
      !# Total of vertical levels
      character(len=32) :: varName(5)
      !# Name os GFS variables
      character(len=200) :: metadata
      !# Metadata from grib2 information
      character(len=30), allocatable :: slevs(:)
      !# Levels description from grib2 file
      character(len=30), allocatable :: slevsInv(:)
      !# Levels description INverted
      real,allocatable :: levpr_grib2(:)
      !# levels values from grib2 file
      character (len=300) :: grid_info
      !# grid information of grib2 file
      character (len=99) :: invline
      !# Var information from grib2 file
      real, allocatable :: var(:,:)
      !# Variable value for each (x,y) point
      real, allocatable :: lat(:,:)
      !# lattitude value for each (x,y) point
      real, allocatable :: lon(:,:)
      !# longitude value for each (x,y) point
      integer :: nxGrib
      !# Size of grib data in lon direction
      integer :: nyGrib
      !# Size of grib data in lat direction
      real :: dx
      !# deltaX longitude increment
      real :: dy
      !# deltaY longitude increment
      integer :: iyear,imonth,iday,ihour
      !# year, month, day and hour read from grib2 file
      real :: latIni,latFin
      !# Initial and final valid latittudes
      real :: lonIni,lonFin
      !# Initial and final valid longitudes
      integer :: jIni,jFin
      !# first and last latittude valid position
      integer :: iIni,Ifin
      !# first and last longitude valid position
      integer :: npi,npj
      !# Total of valid points i,j
      character(len=3) :: clvi
      integer :: iyy,imm,idd,ihh

      integer :: i,j,k,lv,nvar,iCount,jCount,ii,jj,lvi
      integer :: lunit

      character(len=256) :: fName,fNameI,fNameB
      character(len=3) :: ncl

      fName=trim(atmos_idir)//trim(grib2FilesNames(iStep))
      fNameI=trim(grib2InvFilesNames(iStep))
      fNameB=trim(grib2FilesNames(iStep))//'.blow'

      ! Read the header of input pressure file.
      varName(1)=':'//temperature_varname//':'

      !Code

      if(.not. fileExist(trim(fName))) iErrNumber=dumpMessage(c_tty,c_yes,sourceName,procedureName &
         ,c_fatal,'File '//trim(fName) &
         //' not found. Please, verify and solve it!')

      iErrNumber=dumpMessage(c_tty,c_yes,'','',c_notice,'Opening/reading Grib2 file: '//trim(fName))

      !Making inventory using a wgrib2 function grb2_mk_inv
      iErrNumber = grb2_mk_inv(fName,fNameI)
      !print *,iErrNumber
      !if (iErrNumber.ne.1) iErrNumber=dumpMessage(c_tty,c_yes,sourceName,procedureName &
      !           ,c_fatal,' Error making grib2 inventory, check it!')

      ! get number of first var levels using grb2_inq from wgrib2 lib
      nprz = grb2_inq(fName,fNameI,trim(varname(1)),' mb:')
      allocate(slevs(nprz),levpr_grib2(nprz),slevsInv(nprz))
      !
      !print *,'nprz=',nprz
      !
      do i = 1,nprz
         ! get pressure leves and levels labels using grb2_inq from wgrib2 lib
         iErrNumber=grb2_inq(fName,fNameI,trim(varname(1)),' mb:' &
            ,sequential=i-1,desc=metadata)
         if (iErrNumber.ne.1) iErrNumber=dumpMessage(c_tty,c_yes,sourceName,procedureName &
            ,c_fatal,' Error getting pressure levels from grib2 file, check it!')

         j = index(metadata,trim(varname(1))) + len(trim(varname(1)))
         k = index(metadata," mb:") + len(" mb:")-1
         !Transform text values to real for each level
         read(metadata(j:),*) levpr_grib2(i)
         !Fill array with level name+suffix
         slevs(i) = metadata(j-1:k)
      enddo
      !
      !print *,index(metadata," mb:"),index(metadata,'hour')
      !print *,'hour: ',metadata(index(metadata," mb:")+4:index(metadata,'hour')-1)
      !print *,'slevs=',slevs
      !print *,'levpr=',levpr_grib2
      !
      !Getting size and geo information in GRIB2 file. Using the first variable (U)
      iErrNumber = grb2_inq(fName,fNameI,trim(varName(1)) &
         , data2=var &
         , lat=lat, lon=lon, grid_desc=grid_info, desc=invline)

      !Closing grib2 file to save memory
      iErrNumber = grb2_free_file(fName)

      !print *,grid_info
      !print *,invline

      !print *,var

      nxGrib=size(var,1)
      nyGrib=size(var,2)

      !print *,nxGrib,nyGrib

      ! Increments lat and lon
      !dy=(lat(1,nyGrib)-lat(1,1))/(nyGrib-1)
      !dx=(lon(nxGrib,1)-lon(1,1))/(nxGrib-1)
      
      dy =abs(lat(1,2)-lat(1,1))
      dx = lon(2,1)-lon(1,1)
      !print *,dx,dy

      do j=1,nyGrib
         if(lat(1,j)>initial_latitude) then
            latIni=lat(1,j)
            jIni=j
            exit
         endif
      enddo
      do j=jINi,nyGrib
         if(lat(1,j)<final_latitude) then
            latFin=lat(1,j)
            jFin=j
         endif
      enddo
      do i=1,nxGrib
         if(lon(i,1)>initial_longitude) then
            lonIni=lon(i,1)
            iIni=i
            exit
         endif
      enddo
      do i=iINi,nxGrib
         if(lon(i,1)<final_longitude) then
            lonFin=lon(i,1)
            iFin=i
         endif
      enddo
      npi=iFin-iIni+1
      npj=jFin-jIni+1
      print *,jIni,jFin,latIni,latFin,npj
      print *,iIni,iFin,lonIni,lonFin,npi

      !convert date string got from grib2 file to integer
      read(invline( 3: 6),*,iostat=iErrNumber)  iyear
      read(invline( 7: 8),*,iostat=iErrNumber)  imonth
      read(invline( 9:10),*,iostat=iErrNumber)  iday
      read(invline(11:12),*,iostat=iErrNumber)  ihour

      call date_add_to_dble(iyear,imonth,iday,ihour,dble(atmosHoursCount(iStep)),'h' &
         ,iyy,imm,idd,ihh)
      atmosDate(iStep)%year=iyy
      atmosDate(iStep)%month=imm
      atmosDate(iStep)%day=idd
      atmosDate(iStep)%hour=ihh/10000

      !print *,iyear,imonth,iday,ihour

      !Adjustin dpatmos type vars
      atmosNz=levels
      atmosNx=npi
      atmosNy=npj
      atmosLon(1)=lonIni-360.0
      atmosLat(1)=latIni
      atmosXs=dx
      atmosYs=dy
      if(.not. allocated(atmosLevels)) allocate(atmosLevels(atmosNz))
      do lv=1,nprz
         lvi=nprz-lv+1
         if(lvi>atmosNz) cycle
         atmosLevels(lvi)=levpr_grib2(lv)
         slevsInv(lvi)=slevs(lv)
      enddo

      if(iStep==1) then
         atmosLat(2) = atmosLat(1) + (atmosYs*atmosNy)
         atmosLon(2) = atmosLon(1) + (atmosXs*atmosNx)
         iErrNumber=dumpMessage(c_tty,c_yes,'','',c_notice,'Grib2 sizes and levels inventory: ')
         write(*,fmt='("Ntimes: ",I3.3," Nvars: ",I3.3)') 1,5
         write(*,fmt='("NLons : ",I3.3," NLats: ",I3.3," NLevs : ",I3.3)') atmosNx,atmosNy,atmosNz
         write(*,fmt='("LatI  : ",F8.2," LatF : ",F8.2," DeltaY: ",F6.3)') atmosLat(1),atmosLat(2),atmosYs
         write(*,fmt='("LonI  : ",F8.2," LonF : ",F8.2," DeltaX: ",F6.3)') atmosLon(1),atmosLon(2),atmosXs
         write(ncl,fmt='(I3.3)') atmosNz
         write(*,fmt='("Levels: ",'//ncl//'(F6.1,1X))') atmosLevels
      endif

      atmosVarNames = (/wind_u_varname,wind_v_varname,temperature_varname &
         ,geo_varname,ur_varname/)

      if(.not. allocated(atmosValues)) allocate(atmosValues(stepsBetDates,atmosNv,atmosNx, atmosNy, atmosNz))

      lunit=getUnit()
      open(unit=lunit,file=fNameB,action='write',status='replace')

      do lv=1,levels!nprz
         !lvi=lv !nprz-lv+1
         !if(lvi>atmosnZ) cycle
         write(clvi,fmt='(I3.3)') lv
         !Compose name of var to search inside grib2
         varName(1)=':'//trim(wind_u_varname)//trim(slevsInv(lv))
         varName(2)=':'//trim(wind_v_varname)//trim(slevsInv(lv))
         varName(3)=':'//trim(temperature_varname)//trim(slevsInv(lv))
         varName(4)=':'//trim(geo_varname)//trim(slevsInv(lv))
         varName(5)=':'//trim(ur_varname)//trim(slevsInv(lv))
         write(*,fmt='(I2.2,5(A15,1X))') lv,trim(varName(1)),trim(varname(2)),trim(varName(3)),trim(varName(4)),trim(varName(5))
         do nvar=1,5
            !iErrNumber=dumpMessage(c_tty,c_yes,'','',c_notice,'Getting: '//varName(nvar))
            !Getting values for variable nvar
            iErrNumber = grb2_inq(fName,fNameI,trim(varName(nVar)) &
               ,data2=var, lat=lat, lon=lon)
            ii=0
            jj=0
            !Fill values only inside valid area
            do i=iIni,iFin
               ii=ii+1
               jj=0
               do j=jIni,jFin
                  jj=jj+1
                  atmosValues(iStep,nVar,ii,jj,lv)=var(i,j)
                  !Verify var limits
                  if(var(i,j)<limi(nvar)) &
                     iErrNumber=dumpMessage(lunit,c_yes,sourceName,procedureName &
                     ,c_warning,''//' - lev='//clvi//' - Estouro de limite inf. da var. ' &
                     //trim(atmosVarNames(nvar)),(/var(i,j),limi(nvar)/),'E13.4')
                  if(var(i,j)>lims(nvar)) &
                     iErrNumber=dumpMessage(lunit,c_yes,sourceName,procedureName &
                     ,c_warning,''//' - lev='//clvi//' - Estouro de limite sup. da var. ' &
                     //trim(atmosVarNames(nvar)),(/var(i,j),lims(nvar)/),'E13.4')
               enddo
            enddo
         enddo
      end do

      !Closing grib2 file to save memory
      !print *, 'Closing ',trim(file)
      iErrNumber = grb2_free_file(fName)
      iErrNumber = releaseUnit(lunit)

   end subroutine readAtmosGrib2


   !=============================================================================================
   subroutine readAtmosEra5(iStep,stepsBetDates)
      !# Open a era5 file and read the contents
      !#
      !# @note
      !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
      !#
      !# **Brief**: Open era5 ECMWF file, read the contents
      !#
      !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
      !#
      !# **Author(s)**: Luiz Flavio Rodrigues **&#9993;**<luiz.rodrigues@inpe.br>
      !#
      !# **Date**: 17 February 2021 (wednesday)
      !# @endnote
      !#
      !# @changes
      !# &#9744; <br/>
      !# @endchanges
      !# @bug
      !#
      !#@endbug
      !#
      !#@todo
      !#  &#9744; <br/>
      !# @endtodo
      !#
      !# @warning
      !# Now is under CC-GPL License, please see
      !# &copy; <https://creativecommons.org/licenses/GPL/2.0/legalcode.pt>
      !# @endwarning
      !#

      !Use area
      use dump

      use netcdf

      use memoryMod, only: &
         atmos_idir, &
         era5FilesNames, &
         initial_latitude, &
         final_latitude, &
         initial_longitude, &
         final_longitude, &
         atmosNz, &
         atmosNx, &
         atmosNy, &
         atmosLon, &
         atmosLat, &
         atmosXs, &
         atmosYs, &
         atmosNv, &
         atmosLevels, &
         atmosVarNames,&
         atmosValues  , &
         atmosDate, &
         atmosHoursCount, &
         levels

      use utilsMod, only: &
         getUnit, &
         releaseUnit, &
         fileExist, &
         date_add_to_dble


      implicit none

      include "constants.f90"
      include "netcdf.inc"

      character(len=*),parameter :: procedureName='**openatmosEra5**'
      !# Name of this procedure
      !
      !Local Parameters
      character(len=*),parameter :: wind_u_varname='u'
      !# Name of wind u var in grib2 (1)
      character(len=*),parameter :: wind_v_varname='v'
      !# Name of wind v var in grib2 (2)
      character(len=*),parameter :: temperature_varname='t'
      !# Name of temperature var in grib2 (3)
      character(len=*),parameter :: geo_varname='z'
      !# Name of geopotential var in grib2 (4)
      character(len=*),parameter :: ur_varname='r'
      !# Name of relative humidity var in grib2 (5)
      real,dimension(5),parameter :: scale_factor = (/1.,1.,1.,1.,0.01/)
      !# Factor to multiplie each var (1,2,3,4,5)
      real,dimension(5), parameter :: limi=(/-100.,-100.,180., -500.,  0./)
      !# inferior limits for each variable
      real,dimension(5), parameter :: lims=(/ 100., 100.,333.,80000.,100./)
      !# inferior limits for each variable


      integer, parameter :: ncTotal=5
      !# Total number of NetCDF vars

      !Input/Output variables
      integer, intent(in)         :: iStep
      !# timestep
      integer, intent(in)         :: stepsBetDates
      !# total of timesteps

      !Local variables
      integer :: i,j,an,ii,jj,nVar
      integer :: ncid
      integer :: ndims
      integer :: nvars
      integer :: nxEra5
      integer :: nyEra5
      integer :: nzEra5
      integer :: nTimes
      integer :: levelVarn
      integer :: latVarN
      integer :: lonVarN
      integer :: timeN
      integer :: xtype
      integer :: iLen
      integer :: uN
      integer :: vN
      integer :: tN
      integer :: zN
      integer :: qN
      integer :: npi
      integer :: npj
      integer :: jIni
      integer :: iIni
      integer :: jFin
      integer :: iFin
      integer :: iyy,imm,idd,ihh
      integer :: lvi,lv,lunit,ji


      real :: xValue
      real :: missingValue(32)
      real :: mValue
      real :: deltaLat
      real :: deltaLon
      real :: latIni
      real :: latFin
      real :: lonIni
      real :: lonFin

      character(len=256) :: fName,fNameB
      character(len=32),allocatable :: varName(:)
      character(len=32) :: name,atName,vunit
      integer, allocatable :: nat(:)
      integer, allocatable :: varDim(:)
      integer, allocatable :: lenDim(:)
      real, allocatable :: hoursSince1900(:)
      real, allocatable :: PressLevs(:)
      real, allocatable :: llat(:)
      real, allocatable :: llon(:)
      real :: scaleFactor(32)
      real :: addOfset(32)
      character(len=16) :: units(32)
      character(len=3) :: ncl,clvi

      type nct
         real, allocatable  :: values(:,:,:) !i,j,k
      end type nct
      type (nct) :: ncVar(ncTotal)

      !Code
      initial_longitude=initial_longitude-360.0
      final_longitude=final_longitude-360.0


      fName=trim(atmos_idir)//trim(era5FilesNames(iStep))
      fNameB=trim(era5FilesNames(iStep))//'.blow'

      if(.not. fileExist(trim(fName))) iErrNumber=dumpMessage(c_tty,c_yes,sourceName,procedureName &
         ,c_fatal,'File '//trim(fName) &
         //' not found. Please, verify and solve it!')

      iErrNumber=dumpMessage(c_tty,c_yes,'','',c_notice,'Opening/reading era5 file: '//trim(fName))

      !Open NetCDF file
      !write(*,fmt='(A,A,A)') '{',trim(innpr),'}'
      iErrNumber = nf90_open(path = trim(fname), mode = nf90_nowrite, ncid = ncid)
      ! get info about netCDF file: number of vars and number of dims
      iErrNumber=nf90_inquire(ncid, ndims, nvars)

      allocate(varName(nvars),varDim(nvars),nat(nvars),lenDim(nvars))

      ! Get position in data for levels, lons, lats
      nzEra5=0
      nxEra5=0
      nyEra5=0
      do i=1,nvars
         iErrNumber = nf90_Inquire_Dimension(ncid, i, name, lenDim(i))
         !print *,name
         select case (trim(name))
          case ('level')
            nzEra5=lenDim(i)
            levelVarN=i
          case ('longitude')
            nxEra5=lenDIm(i)
            lonVarN=i
          case ('latitude')
            nyEra5=lenDIm(i)
            latVarN=i
          case ('time')
            nTimes=lendim(i)
            timeN=i
         end select
      enddo

      !Allocate all vars
      allocate(PressLevs(nzEra5))
      allocate(llat(nyEra5),llon(nxEra5))
      allocate(hoursSince1900(nTimes))

      !Get pressure levels
      iErrNumber=nf90_get_var(ncid, levelVarN, PressLevs)

      !Get names of vars in netCDF file and total of atributes
      do i=1,nvars
         iErrNumber=nf90_Inquire_Variable(ncid, i, name=varName(i) &
            ,ndims=varDim(i), nAtts=nat(i))
         select case (trim(varName(i)))
          case (wind_u_varname)
            !iErrNumber=dumpMessage(c_tty,c_yes,'','',c_notice,'# wind U, u, present!')
            uN=i
          case (wind_v_varname)
            !iErrNumber=dumpMessage(c_tty,c_yes,'','',c_notice,'# wind V, v, present!')
            vN=i
          case (temperature_varname)
            !iErrNumber=dumpMessage(c_tty,c_yes,'','',c_notice,'# Temperature, t, present!')
            tN=i
          case (geo_varname)
            !iErrNumber=dumpMessage(c_tty,c_yes,'','',c_notice,'# Geopotencial, z, present!')
            zN=i
          case (ur_varname)
            !iErrNumber=dumpMessage(c_tty,c_yes,'','',c_notice,'# Moisture, q, present!')
            qN=i
         end select
      enddo

      do i=1,nvars
         do an=1,nat(i)
            iErrNumber=nf90_inq_attname(ncid, i, an, atName)
            iErrNumber=nf90_Inquire_Attribute(ncid, i, trim(atName), xtype=xtype, len=iLen)
            if(trim(atName)=='scale_factor') then
               iErrNumber = nf90_get_att(ncid, i, atName, xValue)
               scaleFactor(i)=xvalue
            elseif(trim(atName)=='add_offset') then
               iErrNumber = nf90_get_att(ncid, i, atName, xValue)
               addOfset(i)=xvalue
            elseif(trim(atName)=='units') then
               iErrNumber = nf90_get_att(ncid, i, atName, vunit)
               units(i)=vunit
            elseif(trim(atName)=='missing_value') then
               iErrNumber = nf90_get_att(ncid, i, atName, mValue)
               missingValue(i)=mvalue
            endif
         enddo
         !print *,i,scaleFactor(i),addOfset(i),units(i),missingValue(i)
      end do

      !get lon and lat info
      iErrNumber = nf90_get_var(ncid, latVarN, llat)
      iErrNumber = nf90_get_var(ncid, lonVarN, llon)
      iErrNumber = nf90_get_var(ncid, timeN  , hoursSince1900)

      !write (*,fmt='(A,F16.0)') "Hours past since 1900-01-01 00h: ",HoursSince1900(1)
      call date_add_to_dble(1900,1,1,0,dble(HoursSince1900(1)),'h' &
         ,iyy,imm,idd,ihh)
      write (*,fmt='("Year: ",I4.4,"- Month: ",I2.2,"- Day: ",I2.2,"- Hour: ",I2.2)') iyy,imm &
         ,idd,ihh/10000


      do i=1,ncTotal
         allocate(ncVar(i)%values(nxEra5,nyEra5,nzEra5))
      enddo

      atmosDate(iStep)%year =iyy
      atmosDate(iStep)%month=imm
      atmosDate(iStep)%day  =idd
      atmosDate(iStep)%hour =ihh/10000



      if(iStep==1) then
         deltaLat=llat(1)-llat(2) !ECMWF data is from +90 to -90
         deltaLon=llon(2)-llon(1)

         do j=1,nyEra5
            if(llat(j)>=final_latitude) then
               latFin=llat(j)
               jFin=j
            endif
         enddo
         do j=jFin,nyEra5
            if(llat(j)<=initial_latitude) then
               latIni=llat(j)
               jIni=j
               exit
            endif
         enddo
         do i=1,nxEra5
            !print *,i,llon(i),initial_longitude
            if(llon(i)>=initial_longitude) then
               lonIni=llon(i)
               iIni=i
               exit
            endif
         enddo
         do i=iINi,nxEra5
            if(llon(i)<=final_longitude) then
               lonFin=llon(i)
               iFin=i
            endif
         enddo
         npi=iFin-iIni+1
         npj=jIni-jFin+1

         atmosVarNames = (/wind_u_varname,wind_v_varname,temperature_varname &
            ,geo_varname,ur_varname/)

         !Adjustin dpatmos type vars
         atmosNz=levels
         atmosNx=npi
         atmosNy=npj
         atmosLon(1)=lonIni
         atmosLat(1)=latIni
         atmosXs=deltaLon
         atmosYs=deltaLat
         if(.not. allocated(atmosLevels)) allocate(atmosLevels(atmosNz))
         do lv=1,nzEra5
            lvi=nzEra5-lv+1
            if(lvi>atmosNz) cycle
            atmosLevels(lvi)=PressLevs(lv)
            !slevsInv(lvi)=slevs(lv)
         enddo

         atmosLat(2) = atmosLat(1) + (atmosYs*atmosNy)
         atmosLon(2) = atmosLon(1) + (atmosXs*atmosNx)
         iErrNumber=dumpMessage(c_tty,c_yes,'','',c_notice,'Grib2 sizes and levels inventory: ')
         write(*,fmt='("Ntimes: ",I3.3," Nvars: ",I3.3)') 1,5
         write(*,fmt='("NLons : ",I3.3," NLats: ",I3.3," NLevs : ",I3.3)') atmosNx,atmosNy,atmosNz
         write(*,fmt='("LatI  : ",F8.2," LatF : ",F8.2," DeltaY: ",F6.3)') atmosLat(1),atmosLat(2),atmosYs
         write(*,fmt='("LonI  : ",F8.2," LonF : ",F8.2," DeltaX: ",F6.3)') atmosLon(1),atmosLon(2),atmosXs
         write(ncl,fmt='(I3.3)') atmosNz
         write(*,fmt='("Levels: ",'//ncl//'(F6.1,1X))') atmosLevels

      endif

      if(.not. allocated(atmosValues)) allocate(atmosValues(stepsBetDates,atmosNv,atmosNx, atmosNy, atmosNz))

      !Read all upper air values for u,v,t,z and q
      do i=1,ncTotal
         select case (i)
          case(1)
            iErrNumber = nf90_get_var(ncid, uN, ncVar(i)%values)
            ncVar(i)%values=ncVar(i)%values*scaleFactor(uN) &
               +addOfset(uN)
          case(2)
            iErrNumber = nf90_get_var(ncid, vN, ncVar(i)%values)
            ncVar(i)%values=ncVar(i)%values*scaleFactor(vN) &
               +addOfset(vN)
          case(3)
            iErrNumber = nf90_get_var(ncid, tN, ncVar(i)%values)
            ncVar(i)%values=ncVar(i)%values*scaleFactor(tN) &
               +addOfset(tN)
          case(4)
            iErrNumber = nf90_get_var(ncid, zN, ncVar(i)%values)
            ncVar(i)%values=(ncVar(i)%values*scaleFactor(zN) &
               +addOfset(zN))/c_grav
          case(5)
            iErrNumber = nf90_get_var(ncid, qN, ncVar(i)%values)
            ncVar(i)%values=ncVar(i)%values*scaleFactor(qN) &
               +addOfset(qN)
         end select
      enddo

      iErrNumber = nf90_close(ncid)

      lunit=getUnit()
      open(unit=lunit,file=fNameB,action='write',status='replace')

      do lv=1,nzEra5
         lvi=nzEra5-lv+1
         if(lvi>atmosNz) cycle
         write(clvi,fmt='(I3.3)') lvi
         write(*,fmt='(F6.0,":",5(A1,":"))') atmosLevels(lvi),trim(atmosVarNames(1)),trim(atmosVarNames(2)) &
            ,trim(atmosVarNames(3)) &
            ,trim(atmosVarNames(4)),trim(atmosVarNames(5))
         do nvar=1,5
            !iErrNumber=dumpMessage(c_tty,c_yes,'','',c_notice,'Getting: '//varName(nvar))
            !Getting values for variable nvar
            ii=0
            !Fill values only inside valid area
            do i=iIni,iFin
               ii=ii+1
               jj=0
               do j=jFin,jIni
                  jj=jj+1
                  ji=atmosNy-jj+1
                  if (ji>atmosNy .or. ii>atmosNx) then
                     iErrNumber=dumpMessage(c_tty,c_yes,sourceName,procedureName &
                        ,c_warning,''//' - lev='//clvi//' - Estouro de limites (i,j,atmosNx,atmosNy): ' &
                        ,(/ii,ji,atmosNx,atmosNx/),'I5')
                  endif
                  atmosValues(iStep,nVar,ii,ji,lvi)=ncVar(nVar)%values(i,j,lv)
                  !Verify var limits
                  if(ncVar(nVar)%values(i,j,lv)<limi(nvar)) &
                     iErrNumber=dumpMessage(lunit,c_yes,sourceName,procedureName &
                     ,c_warning,''//' - lev='//clvi//' - Estouro de limite inf. da var. ' &
                     //trim(atmosVarNames(nvar)),(/ncVar(nVar)%values(i,j,lv),limi(nvar)/),'E13.4')
                  if(ncVar(nVar)%values(i,j,lv)>lims(nvar)) &
                     iErrNumber=dumpMessage(lunit,c_yes,sourceName,procedureName &
                     ,c_warning,''//' - lev='//clvi//' - Estouro de limite sup. da var. ' &
                     //trim(atmosVarNames(nvar)),(/ncVar(nVar)%values(i,j,lv),lims(nvar)/),'E13.4')
               enddo
            enddo
         enddo
      enddo

      iErrNumber = releaseUnit(lunit)

   end subroutine readAtmosEra5



   !=============================================================================================
   subroutine createCamsFilesNames(monthBetDates)
      !# Creates the set of Cams chem climatology names
      !#
      !# @note
      !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
      !#
      !# **Brief**: Creates the set of Cams chem climatology names and store in memoryMod
      !#
      !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
      !#
      !# **Author(s)**: Luiz Flavio Rodrigues **&#9993;**<luiz.rodrigues@inpe.br>
      !#
      !# **Date**: 27 August 2020 (Thursday)
      !# @endnote
      !#
      !# @changes
      !# &#9744; <br/>
      !# @endchanges
      !# @bug
      !#
      !#@endbug
      !#
      !#@todo
      !#  &#9744; <br/>
      !# @endtodo
      !#
      !# @warning
      !# Now is under CC-GPL License, please see
      !# &copy; <https://creativecommons.org/licenses/GPL/2.0/legalcode.pt>
      !# @endwarning
      !#

      !Use area
      use dump

      use memoryMod, only: &
         monthCount, &
         cams1FilesNames, &
         chem1_prefix, &
         chem1_sufix, &
         chemDate

      use utilsMod, only: &
         to_upper

      implicit none

      include "constants.f90"
      character(len=*),parameter :: procedureName='**createCamsFilesNames**' !Name of this procedure
      character(len=71) :: line='+--+-+----------------------------------------------------------------+'
      !
      !Local Parameters

      !Input/Output variables
      logical, intent(in) :: monthBetDates(12)
      !# Array with true for each month present

      !Local variables
      integer :: i,iCount
      character(len=3) :: monthCamsName(12)

      !Code
      monthCount=0
      do i=1,12
         if(monthBetDates(i)) then
            monthCount=monthCount+1
            monthCamsName(monthCount)=to_upper(month_Name(i))
         endif
      enddo
      if(allocated(cams1FilesNames)) iErrNumber=dumpMessage(c_tty,c_yes,sourceName,procedureName &
         ,c_fatal,' cams1FilesNames array already allocated with size=',size(cams1FilesNames),"I2.2")
      allocate(cams1FilesNames(monthCount))
      do i=1,monthCount
         write(cams1FilesNames(i),fmt='(A,A,A)') trim(chem1_prefix),monthCamsName(i),trim(chem1_sufix)
      enddo

      iErrNumber=dumpMessage(c_tty,c_yes,'','',c_notice,'Cams file names inventory:')
      write(*,fmt='(A)') line
      do i=1,monthCount
         write(*,fmt='("|",I2.2,"|",I1,"|",A64,"|")') i,1,trim(cams1FilesNames(i))
         !write(*,fmt='("|",A2,"|",I1,"|",A64,"|")') "",2,trim(cams2FilesNames(i))
         write(*,fmt='(A)') line
      enddo

   end subroutine createCamsFilesNames

   !=============================================================================================
   subroutine readCamsFile()
      !# Read and store data from cams files
      !#
      !# @note
      !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
      !#
      !# **Brief**: REad and store data from camsFiles
      !#
      !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
      !#
      !# **Author(s)**: Luiz Flavio Rodrigues **&#9993;**<luiz.rodrigues@inpe.br>
      !#
      !# **Date**: 27 August 2020 (Thursday)
      !# @endnote
      !#
      !# @changes
      !# &#9744; <br/>
      !# @endchanges
      !# @bug
      !#
      !#@endbug
      !#
      !#@todo
      !#  &#9744; <br/>
      !# @endtodo
      !#
      !# @warning
      !# Now is under CC-GPL License, please see
      !# &copy; <https://creativecommons.org/licenses/GPL/2.0/legalcode.pt>
      !# @endwarning
      !#

      !Use area
      use dump

      use netcdf

      use memoryMod, only: &
         monthCount, &
         cams1FilesNames, &
         chem_idir, &
         chemDate, &
         spcCamsName, &
         whichCams, &
         cams

      use utilsMod, only: &
         fileExist, &
         getUnit, &
         releaseUnit, &
         date_add_to_dble, &
         outRealSize

      use chemMod, only: &
         nSpc

      implicit none

      include "constants.f90"
      include "netcdf.inc"
      character(len=*),parameter :: procedureName='**readCamsFile**' !Name of this procedure
      !
      !Local Parameters

      !Input/Output variables

      !Local variables
      integer :: nFile
      integer :: ncid
      character(len=256) :: fName
      character(len=32), allocatable :: varName(:)
      character(len=32) :: name,atName,cValueAt
      integer, allocatable :: varDim(:),nat(:),lenDim(:)
      integer :: levelVarn,lonVarN,latVarN,timeVarN
      integer :: nVars,nDims,an,nnt,ii,jj,nvar
      real, allocatable :: llat(:),llon(:)
      integer, allocatable :: HoursFrom1900(:)
      integer :: nx,ny,nz,nt,ki,k
      integer :: j
      integer :: i
      real :: latIni,latFin
      real :: lonIni,LonFin
      real, allocatable :: Values(:,:,:,:,:)
      double precision :: xValue
      integer :: iLen,xtype,lunit
      integer :: iyy,imm,idd,ihh,nlev,levI
      real, allocatable :: localLevels(:),offset(:),factor(:)
      character(len=3) :: cnz
      integer :: recordLen, iRec
      character(len=256) :: localFileName,cValue
      character(len=15) :: tDef
      integer, allocatable :: varNum(:)
      integer :: varCount,ji


      !Code
      do nfile=1,monthCount
         fName=trim(chem_idir)//trim(cams1FilesNames(nFile))//'.nc'
         if(.not. fileExist(trim(fName))) iErrNumber=dumpMessage(c_tty,c_yes,sourceName,procedureName &
            ,c_fatal,'File '//trim(fName) &
            //' not found. Please, verify and solve it!')
         iErrNumber=dumpMessage(c_tty,c_yes,'','',c_notice,'Opening/reading cams file: '//trim(fName))

         !Open the NetCDF file
         iErrNumber = nf90_open(path = trim(fName), mode = nf90_nowrite, ncid = ncid)
         ! get info about netCDF file
         iErrNumber=nf90_inquire(ncid, ndims, nvars)
         !print *, ndims,nvars
         if(nFile==1) then
            allocate(varName(nvars),varDim(nvars),nat(nvars),lenDim(nvars))
         endif

         do i=1,nvars
            iErrNumber = nf90_Inquire_Dimension(ncid, i, name, lenDim(i))
            select case (trim(name))
             case ('level')
               levelVarn=i
               nz=lenDim(i)
             case ('longitude')
               lonVarN=i
               nx=lenDIm(i)
             case ('latitude')
               latVarN=i
               ny=lenDIm(i)
             case ('time')
               timeVarN=i
               nt=lenDIm(i)
            end select
            !print *,i,name,lenDIm(i)
         enddo
         if(nFile==1) then
            allocate(cams%lat(ny),cams%lon(nx))
            allocate(llat(ny))
            allocate(values(nt,nVars,nx,ny,nz))
            if(allocated(chemDate)) iErrNumber=dumpMessage(c_tty,c_yes,sourceName,procedureName &
               ,c_fatal,' chemDate array already allocated with size=',size(chemDate),"I2.2")
            allocate(chemDate(monthCount,nt))
         endif
         iErrNumber=nf90_get_var(ncid, latVarN, llat)
         iErrNumber=nf90_get_var(ncid, lonVarN, cams%lon)

         cams%nx=nx
         cams%ny=ny
         cams%nz=nz
         cams%nt=nt

         !Invertendo latitudes (cams eh de +90 -> -90)
         do j=1,ny
            ji=ny-j+1
            cams%lat(ji)=llat(j)
         enddo
         do i=1,nx
            cams%lon(i)=cams%lon(i)-360.
         enddo

         if(nFile==1) then
            allocate(localLevels(nz))
            allocate(cams%HoursFrom1900(nt))
         endif

         !Read levels, dates, latitudes and longitudes
         iErrNumber=nf90_get_var(ncid, levelVarN, localLevels)
         iErrNumber=nf90_get_var(ncid, timeVarN, cams%HoursFrom1900)


         if(nfile==1) then
            allocate(varNum(nvars))
            allocate(cams%Levels(Nz))
         endif

         !Invertendo os niveis (cams eh de top->bottom)
         do i=1,nz
            ii=nz-i+1
            cams%Levels(ii)=localLevels(i)
         enddo
         !Deltas lon e lat
         cams%dlon=cams%lon(2)-cams%lon(1)
         cams%dlat=cams%lat(2)-cams%lat(1)


         do i=1,nt
            !write (*,fmt='(A,I2,A,I16)') "Hours for ",i," is ",HoursFrom1900(i)
            call date_add_to_dble(1900,1,1,0,dble(cams%HoursFrom1900(i)),'h' &
               ,iyy,imm,idd,ihh)
            chemDate(nFile,i)%month=imm
            !print *,'LFR: ',i,chemDate(nFile,i)%month,HoursFrom1900(i)
            chemDate(nFile,i)%hour=ihh/10000
            !print *,'LFR: ',i,chemDate(nFile,i)%hour,HoursFrom1900(i)
         enddo

         !Get var Names
         varCount=0
         do i=1,nvars
            iErrNumber=nf90_Inquire_Variable(ncid, i, name=varName(i) &
               ,ndims=varDim(i), nAtts=nat(i))
            !Percorre as variaveis do chem1_list lidas na tabela em gas_equivalence
            do j=1,nSpc
               !Percorre até 5 membros que compoe a especie
               do k=1,5
                  if(varName(i)==spcCamsName(j,k)) then !se eh a especie
                     whichCams(j,k)=1 ! Marca que foi encontrada no arquivo
                     varCount=varCount+1 !ACrescenta uma variavel
                     varNum(varCount)=i !Guarda o numero dessa variavel
                     !print *,i,varName(i),varCount
                  endif
               enddo
            enddo
         enddo
         cams%Nv=varCount


         if(nFile==1) then
            !Allocate vars to read
            allocate(Offset(cams%nv))
            allocate(Factor(cams%nv))
            Factor=1.0
            Offset=0.0
            allocate(cams%Values(monthCount,nt,cams%Nv,Nx,Ny,Nz))
            allocate(cams%LongName(cams%Nv),cams%Units(cams%Nv))
            allocate(cams%varName(cams%Nv))
         endif

         do i=1,varCount
            cams%varName(i)=varName(varNum(i))
            do an=1,nat(varNum(i))
               iErrNumber=nf90_inq_attname(ncid, varNum(i), an, atName)
               iErrNumber=nf90_Inquire_Attribute(ncid, varNum(i), trim(atName), xtype=xtype, len=iLen)
               if(trim(atName)=='scale_factor') then
                  iErrNumber = nf90_get_att(ncid, varNum(i), atName, xValue)
                  Factor(i)=xvalue
               elseif(trim(atName)=='add_offset') then
                  iErrNumber = nf90_get_att(ncid, varNum(i), atName, xValue)
                  Offset(i)=xvalue
               elseif(trim(atName)=='long_name') then
                  iErrNumber = nf90_get_att(ncid, varNum(i), atName, cValue)
                  cams%LongName(i)=cValue
               elseif(trim(atName)=='units') then
                  iErrNumber = nf90_get_att(ncid,varNum(i), atName, cValue)
                  cams%Units(i)=cValue
               endif
            enddo
         end do
         iErrNumber=dumpMessage(c_tty,c_yes,'','',c_notice,'Cams sizes and levels inventory for file #1: ')
         write(*,fmt='("Ntimes: ",I3.3," Nvars: ",I3.3)') cams%Nt,cams%Nv
         write(*,fmt='("NLons : ",I3.3," NLats: ",I3.3," NLevs : ",I3.3)') Nx,Ny,Nz
         write(*,fmt='("LatI  : ",F8.2," LatF : ",F8.2," DeltaY: ",F6.3)') cams%Lat(1),cams%Lat(ny),cams%dlat
         write(*,fmt='("LonI  : ",F8.2," LonF : ",F8.2," DeltaX: ",F6.3)') cams%Lon(1),cams%Lon(nx),cams%dlon
         write(cnz,fmt='(I3.3)') cams%Nz
         write(*,fmt='("Levels: ",'//cnz//'(F6.1,1X))') cams%Levels
         write(*,fmt='("Dates : ",8(F7.0,1X))') cams%HoursFrom1900

         !Reading the data.
         do nnt=1,cams%Nt
            iErrNumber=dumpMessage(c_tty,c_yes,'','',c_notice,'Reading Cams #1 variables for t: ',nnt,'I2.2')
            write(*,fmt='(A2,1X,A16,1X,A40,1X,A16,1X,4(A18,1X))') '#N','Name','Long Name' &
               ,'Units','x factor','+ Offset','Max Value','Min value'
            do i=1,cams%nv
               iErrNumber = nf90_get_var(ncid, varNum(i), Values(nnt,varNum(i),:,:,:) &
                  ,start = (/ 1, 1, 1,nnt /))
               write(*,fmt='(I2.2,1X,A16,1X,A40,1X,A16,1X,4(E18.4,1X))') i,cams%VarName(i),cams%LongName(i) &
                  ,cams%Units(i),factor(i),Offset(i),maxval(Values(nnt,varNum(i),:,:,:)),minval(Values(nnt,varNum(i),:,:,:))
            enddo
         enddo

         do nnt=1,cams%Nt
            do nvar=1,cams%Nv
               do j=1,cams%Ny
                  ji=cams%Ny-j+1 !Invertendo a latitude
                  do k=1,cams%Nz
                     ki=cams%Nz-k+1 !Invertendo os niveis
                     cams%Values(nfile,nnt,nVar,:,ji,ki)=Values(nnt,varNum(nVar),:,j,k) &
                        *Factor(nVar)+Offset(nVar)
                  enddo
               enddo
            enddo
         enddo

         iErrNumber = nf90_close(ncid)

         write(tDef,fmt='(I2.2,":00z",I2.2,A3,I4.4)')  0,1 &
            ,month_Name(chemDate(nFile,1)%month),2019
         write(localFileName,fmt='("Leitura1-",I2.2,"-",I2.2)') nFile,1
         lunit=getUnit()
         open(unit=lunit,file=trim(localFileName)//'.ctl',action='write',status='replace')
         !writing the name of grads file
         write(lunit,*) 'dset ^'//trim(localFileName)//'.gra'
         !writing others infos to ctl
         write(lunit,*) 'undef -0.9990000E+34'
         write(lunit,*) 'title '//procedureName
         write(lunit,*) 'xdef ',cams%Nx,' linear ',cams%Lon(1),cams%dlon
         write(lunit,*) 'ydef ',cams%Ny,' linear ',cams%Lat(1),cams%dlat
         write(lunit,*) 'zdef ',cams%Nz,'levels',cams%Levels
         write(lunit,*) 'tdef ',nnt,' linear '//tDef//' 3hr'
         write(lunit,*) 'vars ',cams%Nv
         do i=1,cams%Nv
            write(lunit,*) trim(cams%VarName(i)),cams%Nz,'99 ',trim(cams%LongName(i))//' '//trim(cams%Units(i))
         enddo
         write(lunit,*) 'endvars'
         ierrNumber=releaseUnit(lunit)

         recordLen=outRealSize()*cams%Nx*cams%Ny

         lunit=getUnit()
         open(unit=lunit,file=trim(localFileName)//'.gra',action='WRITE',status='REPLACE' &
            ,form='UNFORMATTED',access='DIRECT',recl=recordLen)

         irec=1
         do nnt=1,cams%Nt
            do nvar=1,cams%Nv
               do k=1,cams%Nz
                  write(lunit,rec=irec) cams%Values(nfile,nnt,nvar,:,:,k)
                  irec=irec+1
               enddo
            enddo
         enddo

      enddo

   end subroutine readCamsFile


   !=============================================================================================
   subroutine writeGradsCtlFile(iTime,nSpc,values,nx,ny,nz,lon,lat,dlon,dlat,levs)
      !# Write the output in grads and ctl files
      !#
      !# @note
      !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
      !#
      !# **Brief**: write output in grads and ctl file with meteorological and chem data
      !#
      !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
      !#
      !# **Author(s)**: Luiz Flavio Rodrigues **&#9993;**<luiz.rodrigues@inpe.br>
      !#
      !# **Date**: 27 August 2020 (Thursday)
      !# @endnote
      !#
      !# @changes
      !# &#9744; <br/>
      !# @endchanges
      !# @bug
      !#
      !#@endbug
      !#
      !#@todo
      !#  &#9744; <br/>
      !# @endtodo
      !#
      !# @warning
      !# Now is under CC-GPL License, please see
      !# &copy; <https://creativecommons.org/licenses/GPL/2.0/legalcode.pt>
      !# @endwarning
      !#

      !Use area
      use dump

      use memoryMod, only: &
         out_type, &
         out_prefix, &
         out_sufix, &
         out_dir, &
         atmosNx, &
         atmosNy, &
         atmosNz, &
         atmosLat, &
         atmosLon, &
         atmosXs, &
         atmosYs, &
         atmosValues, &
         atmosVarNames, &
         atmosDate, &
         atmosNv, &
         atmosLevels , &
         spcName, &
         chem_type

      use utilsMod, only: &
         getUnit , &
         releaseUnit, &
         outRealSize

      implicit none

      include "constants.f90"
      character(len=*),parameter :: procedureName='**writeGradsCtlFile**' !Name of this procedure
      !
      !Local Parameters

      !Input/Output variables
      integer :: iTime
      !#
      integer, intent(in) :: nSpc
      !#
      real, intent(in) :: values(nSpc,atmosNx,atmosNy,atmosNz)
      !#
      integer, intent(in) :: nx,ny,nz
      !#
      real,intent(in) :: lon
      !#
      real,intent(in) :: lat
      !#
      real,intent(in) :: dlon
      !#
      real,intent(in) :: dlat
      !#
      real,intent(in) :: levs(nz)
      !#

      !Local variables
      character(len=256) :: ctlFileName
      character(len=256) :: binFileName
      character(len=15) :: tDef
      integer :: TotalVars,nVar,lunit,recordLen
      integer :: iRec,nnz

      !Code
      !SUm the amount of vars
      if(chem_type>=1) then
         TotalVars=atmosNv+nSpc
      else
         TotalVars=atmosNv
      endif
      !Creating output names
      write(ctlFileName,fmt='(A,I4.4,I2.2,I2.2,I2.2,A)') trim(out_prefix),atmosDate(ITime)%year,atmosDate(Itime)%month &
         ,atmosDate(iTime)%day,atmosDate(iTime)%hour,trim(out_sufix)//'.ctl'
      write(binFileName,fmt='(A,I4.4,I2.2,I2.2,I2.2,A)') trim(out_prefix),atmosDate(ITime)%year,atmosDate(Itime)%month &
         ,atmosDate(iTime)%day,atmosDate(iTime)%hour,trim(out_sufix)//'.gra'

      write(tDef,fmt='(I2.2,":00z",I2.2,A3,I4.4)')  atmosDate(iTime)%hour,atmosDate(iTime)%day &
         ,month_Name(atmosDate(Itime)%month),atmosDate(ITime)%year
      !Writing ctl file
      lunit=getUnit()
      open(unit=lunit, file=trim(out_dir)//trim(ctlFileName), action='write', status='replace')

      !writing the name of grads file
      write(lunit,*) 'dset ^'//trim(binFileName)
      !writing others infos to ctl
      write(lunit,*) 'undef -0.9990000E+34'
      write(lunit,*) 'title PRE BRAMS data'
      write(lunit,*) 'xdef ',Nx-1,' linear ',lon+dlon,dlon
      write(lunit,*) 'ydef ',Ny,' linear ',lat,dlat
      write(lunit,*) 'zdef ',Nz,'levels',levs
      write(lunit,*) 'tdef 1 linear '//tDef//' 1mo'
      write(lunit,*) 'vars ',TotalVars


      do nvar=1,atmosNv
         write(lunit,*) atmosVarNames(nvar),atmosNz,'99 ',atmosVarNames(nvar)
      enddo
      if(chem_type>=1) then
         do nvar=1,nSpc
            write(lunit,*) spcName(nvar),nz,'99 ',spcName(nvar)
         enddo
      endif
      write(lunit,*) 'endvars'
      ierrNumber=releaseUnit(lunit)

      iErrNumber=dumpMessage(c_tty,c_yes,'','',c_notice,'writing grads/ctl file for '//tDef)
      recordLen=outRealSize()*(nx-1)*ny
      iErrNumber=dumpMessage(c_tty,c_yes,'','',c_notice,'dimensions nx,ny,recordLen ',(/nx,ny,recordLen/),"I10")
      lunit=getUnit()
      open(unit=lunit,file=trim(out_dir)//trim(binFileName),action='WRITE',status='REPLACE' &
         ,form='UNFORMATTED',access='DIRECT',recl=recordLen)

      irec=1

      do nvar=1,atmosNv
         do nnz=1,atmosNz
            write(lunit,rec=irec) atmosValues(iTime,nvar,2:nx,:,nnz)
            irec=irec+1
         enddo
      enddo

      if(chem_type>=1) then
         do nvar=1,nSpc
            do nnz=1,atmosNz
               write(lunit,rec=irec) values(nvar,2:nx,:,nnz)
               irec=irec+1
            enddo
         enddo
      endif

   end subroutine writeGradsCtlFile

end module filesMod
