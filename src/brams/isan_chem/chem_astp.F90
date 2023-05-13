!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine chem_pressure_stage(n1, n2, nhem, glat, glon, glat2, glon2)
   use isan_coms

   !--(DMK-CCATT-INI)----------------------------------------------------------------
   !srf-chem
   use chem_isan_coms
   use chem1_list, only: chemical_mechanism ! intent(in)
   use mem_chem1, only: CHEM_ASSIM, CHEMISTRY ! intent(in)
   !srf-chem-end
   !--(DMK-CCATT-FIM)----------------------------------------------------------------

   use rconstants

   implicit none

   integer :: n1, n2, nhem
   real :: glat(n1, n2), glon(n1, n2), glat2(n1, n2), glon2(n1, n2)

   real :: fnprx, fnpry, grx, gry, gglat, gglon, thmax, thmin &
           , xswlon_east, cntlon_east, rr, each_levp_real(maxpr)
   integer :: i, j, k, lv, ifm, loop, n, iunit
   real, external :: rs

   !--(DMK-CCATT-INI)----------------------------------------------------------------
   !srf-chem
   character(len=32) :: chemical_mechanism_test
   !srf-chem-end
   !--(DMK-CCATT-FIM)----------------------------------------------------------------

   ! Read the header of input pressure file.
   xnelat = 0.0; xnelon = 0.0; cntlat = 0.0; cntlon = 0.0; secondlat = 0.0

   print 91, innpr(1:len_trim(innpr))
91 format(//, ' Reading pressure gridded data', /, a, //)

   !write(*,fmt='(2(A3,1X),4(A7,1X))') 'n1','n2','glat','glon','glat2','glon2'
   !write(*,fmt='(2(I3.3,1X),4(F7.2,1X))') n1,n2,glat(n1,n2),glon(n1,n2),glat2(n1,n2),glon2(n1,n2)

   open (11, file=innpr(1:len_trim(innpr)))
   read (11, *) marker, isversion
   if (marker .ne. 999999) isversion = 1
   !  print *,'File to be read:',innpr(1:len_trim(innpr)),isversion
   select case (isversion)
   case (1)
      rewind 11
      read (11, *) iyy, imm, idd, ihh, nprz, nprx, npry, xswlon, xswlat, gdatdx, gdatdy
      read (11, *) (each_levp_real(n), n=1, nprz)
      !dmae - geraDp generating integer, dprep generating real
      do n = 1, nprz
         levpr(n) = int(each_levp_real(n))
      end do

      write (*, *) "Info 1: ", iyy, imm, idd, ihh, nprz, nprx, npry, xswlon, xswlat, gdatdx, gdatdy
      write (*, *) "Levels: ", (levpr(n), n=1, nprz)

      inproj = 1
      if (iyy .lt. 100) iyy = iyy + 1900
      ihh = ihh*100
   case (2)
      print *, 'doing RALPH 2 format'
      read (11, *) iyy, imm, idd, ihh, itinc, nprz, nprx, npry
      read (11, *) inproj, gdatdx, gdatdy, xswlat, xswlon &
         , xnelat, xnelon, cntlat, cntlon, secondlat
      read (11, *) ivertcoord, (levpr(lv), lv=1, nprz)
   end select

   !--(DMK-CCATT-INI)----------------------------------------------------------------
   !srf-chem
   if (CHEM_ASSIM == 1 .and. CHEMISTRY >= 0) then
      !   read chemical mechanism that dp file is prepared for
      read (11, *) chemical_mechanism_test
      !   check consistency
      if (trim(chemical_mechanism_test) /= trim(chemical_mechanism)) then
         print *, 'chemical mechanism not the same as expected'
         print *, 'expected=', trim(chemical_mechanism(1:len_trim(chemical_mechanism)))
         print *, 'read    =', trim(chemical_mechanism_test(1:len_trim(chemical_mechanism_test)))
         stop 'wrong chem mechanism at chem_astp'
      else
         print *, 'ISAN-CHEM data assimilation for chemical mechanism: ', &
            trim(chemical_mechanism(1:len_trim(chemical_mechanism)))
      end if
   end if
   !srf-chem-end
   !--(DMK-CCATT-FIM)----------------------------------------------------------------

   ! Check for consistency between file parameters and namelist parameters

   if (iyy .ne. iyear .or. imm .ne. imonth &
       .or. idd .ne. idate .or. ihh .ne. ihour) then
      print *, 'Pressure file dates not the same as namelist!'
      print *, 'Year :', iyy, iyear
      print *, 'Month:', imm, imonth
      print *, 'Day  :', idd, idate
      print *, 'Hour :', ihh, ihour
      stop 'pr_dates'
   end if

   ! Check pressure data domain size and location

   if (inproj /= 1 .and. nhem > 0) then
      print *, 'You must input a lat-lon pressure grid ' &
         , 'to run a global simulation !!!!'
      stop 'glob-no-press'
   end if

   if (inproj .eq. 1) then

      ! If necessary, convert longitude specification to values in the range
      ! [-180.,180.].

      xswlon = mod(xswlon + 900., 360.) - 180.
      xnelon = mod(xnelon - 900., -360.) + 180.

      if (xswlon .lt. -180. .or. xswlon .ge. 180. .or. &
          xnelon .lt. -180. .or. xnelon .gt. 180.01 .or. &
          xswlat .lt. -90. .or. xnelat .gt. 90.01) then
         print *, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         print *, '!!! BAD DOMAIN SPECIFICATION   !!'
         print *, '!!! xswlat,xswlon - ', xswlat, xswlon
         print *, '!!! xnelat,xnelon-xswlon - ', xnelat, xnelon - xswlon
         print *, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         stop 'bad_domain'
      end if

      fnprx = float(nprx)
      fnpry = float(npry)

      ! Set "global domain" flags to determine whether 2 extra rows will be
      ! added to borders of the input gridded data.  iglobew = 1 if doing
      ! full 360 degrees of longitude.  iglobs = 1 if including the south
      ! pole and if iglobew = 1, and iglobn = 1 if including
      ! the north pole and if iglobew = 1.

      iglobew = 0
      iglobn = 0
      iglobs = 0
      if (fnprx*gdatdx .gt. 359.9) then
         if (abs(xswlon + 180.) .gt. 1.e-3) then
            print *, 'When extracting all longitudes of pressure data,'
            print *, 'xswlon must be set to -180.'
            !stop 'xswlon'
         end if
         iglobew = 1
      end if
      if (iglobew .eq. 1 .and. xswlat .lt. -89.9) iglobs = 1
      if (iglobew .eq. 1 .and. xnelat .gt. 89.9) iglobn = 1
      idatelin = 0
      if (xnelon .lt. xswlon + .01) idatelin = 1

      ! Search coarse grid points for any that are outside the bounds of
      ! the pressure data.  If any are found, stop.

      loop = 1
      if (nhem > 0) loop = 2
      do ifm = 1, loop
         do j = 1, n2
            do i = 1, n1
               if (ifm == 1) then
                  gglat = glat(i, j)
                  gglon = glon(i, j)
               elseif (ifm == 2) then
                  gglat = glat2(i, j)
                  gglon = glon2(i, j)
               end if
               gry = (gglat - xswlat)/gdatdy + 1.
               if (gry .lt. 2 - iglobs .or. gry .gt. fnpry + iglobn - 1) then
                  print *, 'Model grid point latitude must be'
                  print *, 'at least 1 gridpoint inside'
                  print *, 'pressure data area'
                  print *, 'ifm,i,j,glat,xswlat,xnelat' &
                     , ifm, i, j, gglat, xswlat, xnelat
                  stop 'isan 1: outside p-data 1'
               end if
               if (idatelin .eq. 1 .and. gglon .lt. xswlon) gglon = gglon + 360.
               grx = (gglon - xswlon)/gdatdx + 1.
               if (grx .lt. 2 - iglobew .or. grx .gt. fnprx + iglobew - 1) then
                  print *, 'Model grid point longitude must be'
                  print *, 'at least 1 gridpoint inside'
                  print *, 'pressure data area'
                  print *, 'ifm,i,j,glon,xswlon,xnelon' &
                     , ifm, i, j, gglon, xswlon, xnelon, iglobew, iglobs, iglobn
                  stop 'isan 1: outside p-data 2'
               end if
            end do
         end do
      end do
   end if

   ! Deallocate memory for the pressure data

   if (allocated(p_u)) deallocate (p_u)
   if (allocated(p_v)) deallocate (p_v)
   if (allocated(p_t)) deallocate (p_t)
   if (allocated(p_z)) deallocate (p_z)
   if (allocated(p_r)) deallocate (p_r)
   if (allocated(p_ur)) deallocate (p_ur)
   if (allocated(p_vr)) deallocate (p_vr)
   if (allocated(p_lat)) deallocate (p_lat)
   if (allocated(p_lon)) deallocate (p_lon)
   if (allocated(p_slp)) deallocate (p_slp)
   if (allocated(p_sfp)) deallocate (p_sfp)
   if (allocated(p_sft)) deallocate (p_sft)
   if (allocated(p_snow)) deallocate (p_snow)
   if (allocated(p_sst)) deallocate (p_sst)

   !--(DMK-CCATT-INI)----------------------------------------------------------------
   !srf-chem
   if (allocated(p_sc)) deallocate (p_sc)
   !srf-chem-end
   !--(DMK-CCATT-FIM)----------------------------------------------------------------

   ! Allocate memory for the pressure data

   allocate (p_u(nprx, npry, nprz))
   allocate (p_v(nprx, npry, nprz))
   allocate (p_t(nprx, npry, nprz))
   allocate (p_z(nprx, npry, nprz))
   allocate (p_r(nprx, npry, nprz))

   !--(DMK-CCATT-INI)----------------------------------------------------------------
   !srf-chem
   if (CHEM_ASSIM == 1 .and. nspecies > 0) allocate (p_sc(nprx, npry, nprz, nspecies))
   !srf-chem-end
   !--(DMK-CCATT-FIM)----------------------------------------------------------------

   ! p_ur,p_vr arrays for rotated winds used later

   allocate (p_ur(nprx, npry, nprz))
   allocate (p_vr(nprx, npry, nprz))

   allocate (p_lat(nprx, npry))
   allocate (p_lon(nprx, npry))

   allocate (p_slp(nprx, npry))
   allocate (p_sfp(nprx, npry))
   allocate (p_sft(nprx, npry))
   allocate (p_snow(nprx, npry))
   allocate (p_sst(nprx, npry))

   ! Fill with missing in case they are not present
   p_slp(1:nprx, 1:npry) = 1.e30
   p_sfp(1:nprx, 1:npry) = 1.e30
   p_sft(1:nprx, 1:npry) = 1.e30
   p_snow(1:nprx, 1:npry) = 1.e30
   p_sst(1:nprx, 1:npry) = 1.e30

   iunit = 11

   print 92
92 format(////, 1x, 70('*')/ &
           , '  Access coarse resolution pressure data', /, 1x, 70('*')///)

   do k = 1, nprz
      pnpr(k) = levpr(k)*100.
   end do

   ! Call routine to fill pressure arrays from the chosen dataset.

   !--(DMK-CCATT-INI)----------------------------------------------------------------
   call chem_get_press(iunit)

   !--(DMK-CCATT-OLD)----------------------------------------------------------------
   !  call get_press (iunit)
   !--(DMK-CCATT-FIM)----------------------------------------------------------------

    !!!!!!!! Be careful !!!!!!!!!
   !  Check input humidity variable. Assume that if the max of the field is greater
   !  than 1.1 (allow for some machine roundoff),
   !  it is specific humidity (in g/kg) which needs to be converted to rh

    !!$if(maxval(p_r(1:nprx,1:npry,1:nprz)) > 1.1)then
    !!$   print*,'------------------------------------------------'
    !!$   print*,' Converting specific humidity to rh',nprx,npry,nprz
    !!$   print*,'------------------------------------------------'
    !!$   do k=1,nprz
    !!$      do j=1,npry
    !!$         do i=1,nprx
    !!$            rr=p_r(i,j,k)*.001/(1.-p_r(i,j,k)*.001)
    !!$            p_r(i,j,k)=rr/rs(pnpr(k),p_t(i,j,k))
    !!$            if(p_r(i,j,k) <= 0.1) p_r(i,j,k)=.1
    !!$!   if(i.eq.1.and.j.eq.1)print *,k,tn(i,j,k),pnpr(k),rr,rn(i,j,k),zn(i,j,k)
    !!$         enddo
    !!$      enddo
    !!$   enddo
    !!$endif

   ! Find max-min theta at bottom and top levels
   thmax = 1.
   thmin = 1000.
   do j = 1, npry
      do i = 1, nprx
         if (p_t(i, j, nprz) .lt. 1.e20) &
            thmax = max(thmax, p_t(i, j, nprz)*(p00/pnpr(nprz))**rocp)
         if (p_t(i, j, 1) .lt. 1.e20) &
            thmin = min(thmin, p_t(i, j, 1)*(p00/pnpr(1))**rocp)
      end do
   end do

   print 300, levpr(1), thmin, levpr(nprz), thmax
300 format(//, ' Minimum THETA at ', I4, ' mb - ', F8.2/ &
           , ' Maximum THETA at ', I4, ' mb - ', F8.2)

   !  Compute lat-lon at input pressure data points

   select case (inproj)
   case (1)
      do j = 1, npry
         do i = 1, nprx
            p_lat(i, j) = xswlat + (j - 1)*gdatdy
            p_lon(i, j) = xswlon + (i - 1)*gdatdx
         end do
      end do
   case (2)
      print *, 'lamb-con:', xswlat, xswlon, cntlat, cntlon, gdatdx, gdatdy
      do j = 1, npry
         do i = 1, nprx
            call lc2ll(cntlat, cntlon, p_lat(i, j), p_lon(i, j) &
                       , float(i), float(j), xswlat, xswlon, gdatdx)
            !         print*,i,j,p_lat(i,j),p_lon(i,j)
         end do
      end do
   case (3)
      print *, 'polar:', cntlat, cntlon, gdatdx, gdatdy, xswlat, xswlon
      xswlon_east = xswlon
      if (xswlon < 0.) xswlon_east = xswlon + 360.
      cntlon_east = cntlon
      if (cntlon < 0.) cntlon_east = cntlon + 360.
      do j = 1, npry
         do i = 1, nprx
            call w3fb07(float(i), float(j), xswlat, xswlon_east, gdatdx, cntlon_east &
                        , p_lat(i, j), p_lon(i, j))
            !         print*,i,j,p_lat(i,j),p_lon(i,j)
         end do
      end do
   end select

   close (iunit)

   return

!--(DMK-CCATT-INI)----------------------------------------------------------------
end subroutine chem_pressure_stage
!--(DMK-CCATT-OLD)----------------------------------------------------------------
!end subroutine pressure_stage
!--(DMK-CCATT-FIM)----------------------------------------------------------------

subroutine chem_pressure_stage_grib2(n1, n2, nhem, glat, glon, glat2, glon2)
   !--(DMK-CCATT-OLD)----------------------------------------------------------------
   !subroutine pressure_stage(n1,n2,nhem,glat,glon,glat2,glon2)
   !--(DMK-CCATT-FIM)----------------------------------------------------------------

   use isan_coms

   use ModDateUtils, only: &
      date_add_to

   !--(DMK-CCATT-INI)----------------------------------------------------------------
   !srf-chem
   use chem_isan_coms
   use chem1_list, only: chemical_mechanism ! intent(in)
   use mem_chem1, only: CHEM_ASSIM, CHEMISTRY ! intent(in)
   !srf-chem-end
   !--(DMK-CCATT-FIM)----------------------------------------------------------------

   use rconstants

   !use dump, only: &
   !  dumpMessage
#ifdef GRIB2
   use wgrib2api
#endif

   implicit none

   include "constants.f90"
   integer :: n1, n2, nhem
   real :: glat(n1, n2), glon(n1, n2), glat2(n1, n2), glon2(n1, n2)

   real :: fnprx, fnpry, grx, gry, gglat, gglon, thmax, thmin &
           , xswlon_east, cntlon_east, rr, each_levp_real(maxpr)
   integer :: i, j, k, lv, ifm, loop, n, iunit, stat
   real, external :: rs
   character(len=200) :: metadata
   character(len=10) :: varName(5)
   !# composed name of levels (:press:)
   character(len=300) :: grid_info
   character(len=99) :: invline
   real, allocatable :: var(:, :), lat(:, :), lon(:, :)
   integer :: nx, ny, thisHour
   integer :: iyyl, imml, iddl, ihhl
   character(len=3) :: cnz

   !--(DMK-CCATT-INI)----------------------------------------------------------------
   !srf-chem
   character(len=32) :: chemical_mechanism_test
   !srf-chem-end
   !--(DMK-CCATT-FIM)----------------------------------------------------------------

   call str2int(innpr(len_trim(innpr) - 2:len_trim(innpr)), thisHour, stat)
#ifdef GRIB2
   ! iErrNumber = nf90_open(path = trim(innpr), mode = nf90_nowrite, ncid = ncid)
   ! if (status /= nf90_noerr) iErrNumber=dumpMessage(c_tty,c_yes,header,modelVersion &
   !               ,c_fatal,trim(innpr)//' not found. Please, check it!')

   ! Read the header of input pressure file.
   varName(1) = ':'//trim(temperature_varname)//':'

   print 91, innpr(1:len_trim(innpr))
91 format(//, ' Reading pressure gridded netCDF data', /, a, //)

   !print 92,innpr(1:len_trim(innpr))//'.inv'
   !92 format(//,' using the inventory GRIB2 data',/,a,//)

   !write(*,fmt='(2(A3,1X),4(A7,1X))') 'n1','n2','glat','glon','glat2','glon2'
   !write(*,fmt='(2(I3.3,1X),4(F7.2,1X))') n1,n2,glat(n1,n2),glon(n1,n2),glat2(n1,n2),glon2(n1,n2)

   !Making inventory using a wgrib2 function grb2_mk_inv
   iErrNumber = grb2_mk_inv(innpr(1:len_trim(innpr)), innpr(1:len_trim(innpr))//'.inv')

   ! get number of first var levels using grb2_inq from wgrib2 lib
   nprz_grib2 = grb2_inq(innpr(1:len_trim(innpr)), innpr(1:len_trim(innpr))//'.inv' &
                         , trim(varname(1)), ' mb:')
   allocate (slevs(nprz_grib2))

   !If the total of levels ir greater than defined in RAMSIN use the defined
   if (nprz_grib2 > z_max_level) then
      nprz = z_max_level
   else
      nprz = nprz_grib2
   end if

   do i = 1, nprz_grib2
      ! get pressure leves and levels labels using grb2_inq from wgrib2 lib
      iErrNumber = grb2_inq(trim(innpr), trim(innpr)//'.inv', trim(varname(1)), ' mb:' &
                            , sequential=i - 1, desc=metadata)
      if (iErrNumber .ne. 1) stop 2

      j = index(metadata, trim(varname(1))) + len(trim(varname(1)))
      k = index(metadata, " mb:") + len(" mb:") - 1
      read (metadata(j:), *) levpr_grib2(i)
      slevs(i) = metadata(j - 1:k)
      !print *,'metadata(j:): ',metadata(j:),slevs(i),levpr_grib2(i)
   end do

   !Getting information in GRIB2 file. Using the first variable (U)
   iErrNumber = grb2_inq(trim(innpr), trim(innpr)//'.inv', trim(varName(1)) &
                         , data2=var &
                         , lat=lat, lon=lon, grid_desc=grid_info, desc=invline)

   !Closing grib2 file to save memory
   iErrNumber = grb2_free_file(trim(innpr))

   nxGrib = size(var, 1)
   nyGrib = size(var, 2)

   allocate (mask(nxGrib, nyGrib))

   ! Increments lat and lon
   gdatdy = (lat(1, nyGrib) - lat(1, 1))/(nyGrib - 1)
   gdatdx = (lon(nxGrib, 1) - lon(1, 1))/(nxGrib - 1)
   !
   !Make number of points lon and lat accordingly setup in RAMSIN
   nprx = ((final_longitude - initial_longitude)/gdatdx) + 1
   npry = ((final_latitude - initial_latitude)/gdatdy) + 1

   write (cnz, fmt='(I3.3)') nprz_grib2
   write (*, fmt='("Inventory(1): ",3(I5,1X))') nxGrib, nyGrib, nprz_grib2
   write (*, fmt='("Inventory(2): ",4(E18.4,1X))') lon(1, 1), gdatdx, lat(1, 1), gdatdy
   write (*, fmt='("Inventory(3): "'//cnz//'(F6.1,1X))') (levpr_grib2(i), i=1, nprz_grib2)
   !
   !convert date string got from grib2 file to integer
   call str2int(invline(3:6), iyyl, stat)
   call str2int(invline(7:8), imml, stat)
   call str2int(invline(9:10), iddl, stat)
   call str2int(invline(11:12), ihhl, stat)

   !Setting initial lon and lat
   xswlon = initial_longitude - 360.0
   xswlat = initial_latitude
   !
   inproj = 1
   !if(iyy.lt.100) iyy=iyy+1900
   !ihh=(ihh+thisHour)*100
   !print *,iyyl,imml,iddl,ihhl,real(thisHour*3600)
   call date_add_to(iyyl, imml, iddl, ihhl, real(thisHour*3600), 's' &
                    , iyy, imm, idd, ihh)
   ihh = ihh/100
   !
   !  write(*,*) "Info 1: ",iyy,imm,idd,ihh,nprz,nprx,npry,xswlon,xswlat,gdatdx,gdatdy
   !  write(*,*) "Levels: ",(levpr(n),n=1,nprz)

   ! Check for consistency between file parameters and namelist parameters
   if (iyy .ne. iyear .or. imm .ne. imonth &
       .or. idd .ne. idate .or. ihh .ne. ihour) then
      print *, 'Pressure file dates not the same as namelist!'
      print *, 'Year :', iyy, iyear
      print *, 'Month:', imm, imonth
      print *, 'Day  :', idd, idate
      print *, 'Hour :', ihh, ihour
      stop 'pr_dates'
   end if

   ! Check pressure data domain size and location
   if (inproj /= 1 .and. nhem > 0) then
      print *, 'You must input a lat-lon pressure grid ' &
         , 'to run a global simulation !!!!'
      stop 'glob-no-press'
   end if

   if (inproj .eq. 1) then
      ! If necessary, convert longitude specification to values in the range
      ! [-180.,180.].
      xswlon = mod(xswlon + 900., 360.) - 180.
      xnelon = mod(xnelon - 900., -360.) + 180.
      if (xswlon .lt. -180. .or. xswlon .ge. 180. .or. &
          xnelon .lt. -180. .or. xnelon .gt. 180.01 .or. &
          xswlat .lt. -90. .or. xnelat .gt. 90.01) then
         print *, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         print *, '!!! BAD DOMAIN SPECIFICATION   !!'
         print *, '!!! xswlat,xswlon - ', xswlat, xswlon
         print *, '!!! xnelat,xnelon-xswlon - ', xnelat, xnelon - xswlon
         print *, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         stop 'bad_domain'
      end if
      fnprx = float(nprx)
      fnpry = float(npry)
      ! Set "global domain" flags to determine whether 2 extra rows will be
      ! added to borders of the input gridded data.  iglobew = 1 if doing
      ! full 360 degrees of longitude.  iglobs = 1 if including the south
      ! pole and if iglobew = 1, and iglobn = 1 if including
      ! the north pole and if iglobew = 1.
      iglobew = 0
      iglobn = 0
      iglobs = 0
      if (fnprx*gdatdx .gt. 359.9) then
         if (abs(xswlon + 180.) .gt. 1.e-3) then
            print *, 'When extracting all longitudes of pressure data,'
            print *, 'xswlon must be set to -180.'
            !stop 'xswlon'
         end if
         iglobew = 1
      end if
      if (iglobew .eq. 1 .and. xswlat .lt. -89.9) iglobs = 1
      if (iglobew .eq. 1 .and. xnelat .gt. 89.9) iglobn = 1
      idatelin = 0
      if (xnelon .lt. xswlon + .01) idatelin = 1
      ! Search coarse grid points for any that are outside the bounds of
      ! the pressure data.  If any are found, stop.
      loop = 1
      if (nhem > 0) loop = 2
      do ifm = 1, loop
         do j = 1, n2
            do i = 1, n1
               if (ifm == 1) then
                  gglat = glat(i, j)
                  gglon = glon(i, j)
               elseif (ifm == 2) then
                  gglat = glat2(i, j)
                  gglon = glon2(i, j)
               end if
               gry = (gglat - xswlat)/gdatdy + 1.
               if (gry .lt. 2 - iglobs .or. gry .gt. fnpry + iglobn - 1) then
                  print *, 'Model grid point latitude must be'
                  print *, 'at least 1 gridpoint inside'
                  print *, 'pressure data area'
                  print *, 'ifm,i,j,glat,xswlat,xnelat' &
                     , ifm, i, j, gglat, xswlat, xnelat
                  stop 'isan 2: outside p-data 1'
               end if
               if (idatelin .eq. 1 .and. gglon .lt. xswlon) gglon = gglon + 360.
               grx = (gglon - xswlon)/gdatdx + 1.
               if (grx .lt. 2 - iglobew .or. grx .gt. fnprx + iglobew - 1) then
                  print *, 'Model grid point longitude must be'
                  print *, 'at least 1 gridpoint inside'
                  print *, 'pressure data area'
                  print *, 'ifm,i,j,glon,xswlon,xnelon' &
                     , ifm, i, j, gglon, xswlon, xnelon, iglobew, iglobs, iglobn
                  stop 'isan 2: outside p-data 2'
               end if
            end do
         end do
      end do
   end if

   ! Deallocate memory for the pressure data
   if (allocated(p_u)) deallocate (p_u)
   if (allocated(p_v)) deallocate (p_v)
   if (allocated(p_t)) deallocate (p_t)
   if (allocated(p_z)) deallocate (p_z)
   if (allocated(p_r)) deallocate (p_r)
   if (allocated(p_ur)) deallocate (p_ur)
   if (allocated(p_vr)) deallocate (p_vr)
   if (allocated(p_lat)) deallocate (p_lat)
   if (allocated(p_lon)) deallocate (p_lon)
   if (allocated(p_slp)) deallocate (p_slp)
   if (allocated(p_sfp)) deallocate (p_sfp)
   if (allocated(p_sft)) deallocate (p_sft)
   if (allocated(p_snow)) deallocate (p_snow)
   if (allocated(p_sst)) deallocate (p_sst)
   if (allocated(p_sc)) deallocate (p_sc)

   ! Allocate memory for the pressure data
   allocate (p_u(nprx, npry, nprz))
   allocate (p_v(nprx, npry, nprz))
   allocate (p_t(nprx, npry, nprz))
   allocate (p_z(nprx, npry, nprz))
   allocate (p_r(nprx, npry, nprz))

   if (CHEM_ASSIM == 1 .and. nspecies > 0) allocate (p_sc(nprx, npry, nprz, nspecies))

   ! p_ur,p_vr arrays for rotated winds used later
   allocate (p_ur(nprx, npry, nprz))
   allocate (p_vr(nprx, npry, nprz))
   allocate (p_lat(nprx, npry))
   allocate (p_lon(nprx, npry))
   allocate (p_slp(nprx, npry))
   allocate (p_sfp(nprx, npry))
   allocate (p_sft(nprx, npry))
   allocate (p_snow(nprx, npry))
   allocate (p_sst(nprx, npry))

   ! Fill with missing in case they are not present
   p_slp(1:nprx, 1:npry) = 1.e30
   p_sfp(1:nprx, 1:npry) = 1.e30
   p_sft(1:nprx, 1:npry) = 1.e30
   p_snow(1:nprx, 1:npry) = 1.e30
   p_sst(1:nprx, 1:npry) = 1.e30

   print 192
192 format(////, 1x, 70('*')/ &
           , '  Access coarse resolution pressure data', /, 1x, 70('*')///)
   !do k=1,nprz
   !   pnpr(k)=levpr(k)*100.
   !enddo

   ! Call routine to fill pressure arrays from the chosen dataset.
   call chem_get_press_grib2()

   ! Find max-min theta at bottom and top levels
   thmax = 1.
   thmin = 1000.
   do j = 1, npry
      do i = 1, nprx
         if (p_t(i, j, nprz) .lt. 1.e20) &
            thmax = max(thmax, p_t(i, j, nprz)*(p00/pnpr(nprz))**rocp)
         if (p_t(i, j, 1) .lt. 1.e20) &
            thmin = min(thmin, p_t(i, j, 1)*(p00/pnpr(1))**rocp)
      end do
   end do

   print 300, levpr(1), thmin, levpr(nprz), thmax
300 format(//, ' Minimum THETA at ', I4, ' mb - ', F8.2/ &
           , ' Maximum THETA at ', I4, ' mb - ', F8.2)

   !  Compute lat-lon at input pressure data points
   select case (inproj)
   case (1)
      do j = 1, npry
         do i = 1, nprx
            p_lat(i, j) = xswlat + (j - 1)*gdatdy
            p_lon(i, j) = xswlon + (i - 1)*gdatdx
         end do
      end do
   case (2)
      print *, 'lamb-con:', xswlat, xswlon, cntlat, cntlon, gdatdx, gdatdy
      do j = 1, npry
         do i = 1, nprx
            call lc2ll(cntlat, cntlon, p_lat(i, j), p_lon(i, j) &
                       , float(i), float(j), xswlat, xswlon, gdatdx)
            !         print*,i,j,p_lat(i,j),p_lon(i,j)
         end do
      end do
   case (3)
      print *, 'polar:', cntlat, cntlon, gdatdx, gdatdy, xswlat, xswlon
      xswlon_east = xswlon
      if (xswlon < 0.) xswlon_east = xswlon + 360.
      cntlon_east = cntlon
      if (cntlon < 0.) cntlon_east = cntlon + 360.
      do j = 1, npry
         do i = 1, nprx
            call w3fb07(float(i), float(j), xswlat, xswlon_east, gdatdx, cntlon_east &
                        , p_lat(i, j), p_lon(i, j))
            !         print*,i,j,p_lat(i,j),p_lon(i,j)
         end do
      end do
   end select
#else
   write (*, fmt='(A)') 'To use grib 2 the model must be compiled with grib2!'
   stop
#endif

end subroutine chem_pressure_stage_grib2

subroutine chem_get_press_grib2()
   use isan_coms
   use chem_isan_coms
#ifdef GRIB2
   use wgrib2api
#endif
   use dump, only: &
      dumpMessage

   implicit none

   include "constants.f90"
   !
   logical, external :: checkInside

   integer :: ithere(maxpr, 5), isfthere(5)
   character*4 field, idat(5)
   data idat/'T', 'R', 'U', 'V', 'H'/
   character(len=60) :: varName(5)
   real, allocatable::as(:, :)
   character(len=200) :: metadata
   character(len=300) :: grid_info
   character(len=99) :: invline
   real, allocatable :: lat(:, :), lon(:, :)
   real, allocatable :: aux_as(:, :)

   integer :: i, j, k, nv, nvar, misstot, lv, n, lv_aux
   integer :: irec, recordLen, iCount, jCount
   integer :: iAbove(5), iBellow(5)
   real :: vMax(5), vMin(5)
   integer, external :: outRealSize

   !Allocate array for trimmer area
   allocate (aux_as(nprx, npry))
   iAbove = 0
   iBellow = 0
   vMax = uLimit
   vMin = dlimit
#ifdef GRIB2
   allocate (lat(nxGrib, nyGrib), lon(nxGrib, nyGrib), as(nxGrib, nyGrib))
   !  Read upper air fields
   do lv = 1, nprz_grib2

      !Compose name of var to search inside grib2
      varName(1) = ':'//trim(wind_u_varname)//trim(slevs(lv))
      varName(2) = ':'//trim(wind_v_varname)//trim(slevs(lv))
      varName(3) = ':'//trim(temperature_varname)//trim(slevs(lv))
      varName(4) = ':'//trim(geo_varname)//trim(slevs(lv))
      varName(5) = ':'//trim(ur_varname)//trim(slevs(lv))

      do nvar = 1, 5 !+nspecies
         print *, trim(innpr), trim(varName(nVar)); call flush (6)
         iErrNumber = grb2_inq(trim(innpr), trim(innpr)//'.inv', trim(varName(nVar)) &
                               , data2=as, lat=lat, lon=lon)
         iCount = 0
         jCount = 1

         do j = 1, nyGrib
            do i = 1, nxGrib
               !Teste if i,j is inside trimmer area for regional model
               if (checkInside(lat(i, j), lon(i, j))) then

                  ! Mapping the position inside trimmer area
                  iCount = iCount + 1
                  if (iCount > nprx) then
                     iCount = 1
                     jCount = jCount + 1
                  end if

                  !Copying inside trimmer area
                  aux_as(iCount, jCount) = as(i, j)*scale_factor(nvar)

                  call limitVariableWithCount(varName(nVar), aux_as(iCount, jCount) &
                                              , dlimit(nVar), uLimit(nVar), iAbove(nVar), iBellow(nVar) &
                                              , vMax(nVar), vMin(nVar))
                  !
               end if
            end do
         end do
         !Invert pressure leves (Grib2 -> BRAMS) and fill
         lv_aux = nprz_grib2 - lv + 1
         levpr(lv_aux) = levpr_grib2(lv)
         pnpr(lv_aux) = levpr(lv_aux)*100

         !Do not process levels above z_max_level
         if (lv_aux > z_max_level) cycle

         select case (nvar)
         case (1)
            call prfill(nprx, npry, aux_as, p_u(1, 1, lv_aux))
         case (2)
            call prfill(nprx, npry, aux_as, p_v(1, 1, lv_aux))
         case (3)
            call prfill(nprx, npry, aux_as, p_t(1, 1, lv_aux))
         case (4)
            call prfill(nprx, npry, aux_as, p_z(1, 1, lv_aux))
         case (5)
            call prfill(nprx, npry, aux_as, p_r(1, 1, lv_aux))
         end select

         print 555, int(levpr_grib2(lv)), nvar, imonth, idate, iyear, ihour
555      format(' ==  Read pressure field  ', 2i4, ' at ', I2, '/', I2 &
                , '/', I4, I6.4, ' UTC')
      end do
   end do

   do i = 1, 5
      !Print information about values bellow and above limits.
      call dumpLimits(iAbove(i), iBellow(i), varName(i), uLimit(i), dLimit(i), vMin(i), vMax(i))
   end do

   if (ccGradsWrite == 1) call writeGradsSub()

   !Fill surface fields
   do nvar = 1, 5
      as = -999.0

      select case (nvar)
      case (1)
         call prfill(nprx, npry, as, p_slp(1, 1))
      case (2)
         call prfill(nprx, npry, as, p_sfp(1, 1))
      case (3)
         call prfill(nprx, npry, as, p_sft(1, 1))
      case (4)
         call prfill(nprx, npry, as, p_snow(1, 1))
      case (5)
         call prfill(nprx, npry, as, p_sst(1, 1))
      end select
   end do

   ! Check for levels that may be all missing
   ithere(1:nprz, 1:5) = 0
   isfthere(1:5) = 0

   do k = 1, nprz
      do j = 1, npry
         do i = 1, nprx
            if (p_t(i, j, k) > 1e20) ithere(k, 1) = ithere(k, 1) + 1
            if (p_r(i, j, k) > 1e20) ithere(k, 2) = ithere(k, 2) + 1
            if (p_u(i, j, k) > 1e20) ithere(k, 3) = ithere(k, 3) + 1
            if (p_v(i, j, k) > 1e20) ithere(k, 4) = ithere(k, 4) + 1
            if (p_z(i, j, k) > 1e20) ithere(k, 5) = ithere(k, 5) + 1
         end do
      end do
   end do
   do nv = 1, 5
      do k = 1, nprz
         if (ithere(k, nv) < nprx*npry) then
            ithere(k, nv) = 1
         else
            ithere(k, nv) = 0
         end if
      end do
   end do

   misstot = 0
   do nv = 1, 5
      do k = 1, nprz
         if (ithere(k, nv) == 0) misstot = misstot + 1
      end do
   end do

   do j = 1, npry
      do i = 1, nprx
         if (p_slp(i, j) > 1e20) isfthere(1) = isfthere(1) + 1
         if (p_sfp(i, j) > 1e20) isfthere(2) = isfthere(2) + 1
         if (p_sft(i, j) > 1e20) isfthere(3) = isfthere(3) + 1
         if (p_snow(i, j) > 1e20) isfthere(4) = isfthere(4) + 1
         if (p_sst(i, j) > 1e20) isfthere(5) = isfthere(5) + 1
      end do
   end do
   do nv = 1, 5
      if (isfthere(nv) < nprx*npry) then
         isfthere(nv) = 1
      else
         isfthere(nv) = 0
      end if
   end do

   print *, '------------------------------------------------'
   print *, ' Missing parameters: 0 = all missing'
   print *, '------------------------------------------------'
   print '(t20,5(a1,6x))', (idat(n), n=1, 5)

   do k = 1, nprz
      print '(f10.1,t14,5(i7))', pnpr(k), (ithere(k, n), n=1, 5)
   end do

   print *, '------------------------------------------------'
   print *, ' Missing surface parameters: 0 = all missing'
   print *, '------------------------------------------------'
   print *, '              SLP    SFP    SFT    SNOW   SST'

   print '(t14,5(i7))', (isfthere(n), n=1, 5)

   if (misstot .gt. 0) then
      ! Let's see if we can get creative and make up data for the missing fields

      call press_miss(nprx, npry, nprz, p_u, p_v, p_t, p_z, p_r &
                      , ithere, maxpr, levpr)

      ! Check again for missing fields

      misstot = 0
      do nv = 1, 5
         do k = 1, nprz
            if (ithere(k, nv) .eq. 0) misstot = misstot + 1
         end do
      end do

      print *, '------------------------------------------------'
      print *, ' After missing parameters check: 0 = all missing'
      print *, '------------------------------------------------'
      print '(t20,5(a1,6x))', (idat(n), n=1, 5)

      do k = 1, nprz
         print '(f10.1,t14,5(i7))', pnpr(k), (ithere(k, n), n=1, 5)
      end do

   end if

   deallocate (slevs)
   deallocate (mask)
   !write(*,fmt='(A,I4.4,2(I2.2),I4.4)') 'All done for ',iyear,imonth,idate,ihour
   return
#endif
!--(DMK-CCATT-INI)----------------------------------------------------------------
end subroutine chem_get_press_grib2

!--(DMK-CCATT-INI)----------------------------------------------------------------
subroutine chem_get_press(iunit)
   !--(DMK-CCATT-OLD)----------------------------------------------------------------
   !subroutine get_press (iunit)
   !--(DMK-CCATT-FIM)----------------------------------------------------------------

   use isan_coms

   !--(DMK-CCATT-INI)----------------------------------------------------------------
   !srf-chem
   use chem_isan_coms
    !!!!  use chem1_list, only : spc_name
   !srf-chem-end
   !--(DMK-CCATT-FIM)----------------------------------------------------------------

   implicit none

   integer :: iunit

   real, allocatable::as(:, :)
   integer :: ithere(maxpr, 5), isfthere(5)
   character*4 field, idat(5)
   data idat/'T', 'R', 'U', 'V', 'H'/

   integer :: i, j, k, nv, nvar, misstot, lv, n
   integer :: recordLen, irec

   ! Allocate array to hold one level data for one variable

   allocate (as(nprx, npry))

   !  Read upper air fields
   do lv = 1, nprz
      !print*," ===" ,lv,nprz,nspecies;call flush(6)
      !--(DMK-CCATT-INI)----------------------------------------------------------------
      !srf-chem
      !do nvar=1,5
      do nvar = 1, 5 + nspecies
         !srf-chem-end

         !srf-chem
         !- define the type of format that DP was writed and and reads the data
         if (trim(innpr(len_trim(innpr) - 2:len_trim(innpr))) == 'vfm') then
            !print*,'reading vfm file'
            call vfirec(iunit, as, nprx*npry, 'LIN')
         else
            !read(iunit,*,end=70,err=70) ((as(i,j),i=1,nprx),j=1,npry)
            print *, 'reading text file *********'
            read (iunit, *, end=70, err=70) ((as(i, j), i=1, nprx), j=1, npry)
            !if(lv==1 .and. nvar==3) then
            !    do i=1,nprx
            !        do j=1,npry
            !            write (51,fmt='(3(F18.2,1X))') real(i),real(j) &
            !                ,as(i,j)
            !        enddo
            !    enddo
            !endif
         end if
         !
         if (nvar .eq. 1) then
            call prfill(nprx, npry, as, p_u(1, 1, lv))
         elseif (nvar .eq. 2) then
            call prfill(nprx, npry, as, p_v(1, 1, lv))
         elseif (nvar .eq. 3) then
            call prfill(nprx, npry, as, p_t(1, 1, lv))
         elseif (nvar .eq. 4) then
            call prfill(nprx, npry, as, p_z(1, 1, lv))
         elseif (nvar .eq. 5) then
            call prfill(nprx, npry, as, p_r(1, 1, lv))
            !--(DMK-CCATT-INI)----------------------------------------------------------------
            !srf-chem
         elseif (nvar > 5) then
            !do j=1,npry;do i=1,nprx; if(as(i,j)>0.01) print*,'as=',as(i,j);enddo;enddo
            call prfill(nprx, npry, as, p_sc(1, 1, lv, nvar - 5))
            !srf-chem
            !--(DMK-CCATT-FIM)----------------------------------------------------------------
         end if

         print 555, levpr(lv), nvar, imonth, idate, iyear, ihour
555      format(' ==  Read pressure field  ', 2i4, ' at ', I2, '/', I2 &
                , '/', I4, I6.4, ' UTC')
      end do
   end do

   if (ccGradsWrite == 1) call writeGradsSub()

   !  Read surface fields

   do nvar = 1, 5
      read (iunit, *, end=70, err=70) ((as(i, j), i=1, nprx), j=1, npry)
      !     print *,'MaxMin(as): ',maxval(as),minval(as)

      select case (nvar)
      case (1)
         call prfill(nprx, npry, as, p_slp(1, 1))
      case (2)
         call prfill(nprx, npry, as, p_sfp(1, 1))
      case (3)
         call prfill(nprx, npry, as, p_sft(1, 1))
      case (4)
         call prfill(nprx, npry, as, p_snow(1, 1))
      case (5)
         call prfill(nprx, npry, as, p_sst(1, 1))
      end select
   end do

   goto 71

70 continue
   print *, 'Premature end of file or error in pressure input file!'
   print *, 'We''ll close our eyes and pretend it didn''t happen!'
71 continue

   deallocate (as)

   ! Check for levels that may be all missing

   ithere(1:nprz, 1:5) = 0
   isfthere(1:5) = 0

   do k = 1, nprz
      do j = 1, npry
         do i = 1, nprx
            if (p_t(i, j, k) > 1e20) ithere(k, 1) = ithere(k, 1) + 1
            if (p_r(i, j, k) > 1e20) ithere(k, 2) = ithere(k, 2) + 1
            if (p_u(i, j, k) > 1e20) ithere(k, 3) = ithere(k, 3) + 1
            if (p_v(i, j, k) > 1e20) ithere(k, 4) = ithere(k, 4) + 1
            if (p_z(i, j, k) > 1e20) ithere(k, 5) = ithere(k, 5) + 1
         end do
      end do
   end do
   do nv = 1, 5
      do k = 1, nprz
         if (ithere(k, nv) < nprx*npry) then
            ithere(k, nv) = 1
         else
            ithere(k, nv) = 0
         end if
      end do
   end do

   misstot = 0
   do nv = 1, 5
      do k = 1, nprz
         if (ithere(k, nv) == 0) misstot = misstot + 1
      end do
   end do

   do j = 1, npry
      do i = 1, nprx
         if (p_slp(i, j) > 1e20) isfthere(1) = isfthere(1) + 1
         if (p_sfp(i, j) > 1e20) isfthere(2) = isfthere(2) + 1
         if (p_sft(i, j) > 1e20) isfthere(3) = isfthere(3) + 1
         if (p_snow(i, j) > 1e20) isfthere(4) = isfthere(4) + 1
         if (p_sst(i, j) > 1e20) isfthere(5) = isfthere(5) + 1
      end do
   end do
   do nv = 1, 5
      if (isfthere(nv) < nprx*npry) then
         isfthere(nv) = 1
      else
         isfthere(nv) = 0
      end if
   end do

   print *, '------------------------------------------------'
   print *, ' Missing parameters: 0 = all missing'
   print *, '------------------------------------------------'
   print '(t20,5(a1,6x))', (idat(n), n=1, 5)

   do k = 1, nprz
      print '(f10.1,t14,5(i7))', pnpr(k), (ithere(k, n), n=1, 5)
   end do

   print *, '------------------------------------------------'
   print *, ' Missing surface parameters: 0 = all missing'
   print *, '------------------------------------------------'
   print *, '              SLP    SFP    SFT    SNOW   SST'

   print '(t14,5(i7))', (isfthere(n), n=1, 5)

   if (misstot .gt. 0) then
      ! Let's see if we can get creative and make up data for the missing fields

      call press_miss(nprx, npry, nprz, p_u, p_v, p_t, p_z, p_r &
                      , ithere, maxpr, levpr)

      ! Check again for missing fields

      misstot = 0
      do nv = 1, 5
         do k = 1, nprz
            if (ithere(k, nv) .eq. 0) misstot = misstot + 1
         end do
      end do

      print *, '------------------------------------------------'
      print *, ' After missing parameters check: 0 = all missing'
      print *, '------------------------------------------------'
      print '(t20,5(a1,6x))', (idat(n), n=1, 5)

      do k = 1, nprz
         print '(f10.1,t14,5(i7))', pnpr(k), (ithere(k, n), n=1, 5)
      end do

   end if

   return

   !--(DMK-CCATT-INI)----------------------------------------------------------------
end subroutine chem_get_press

#ifdef cdf
subroutine chem_pressure_stage_netCDF(n1, n2, nhem, glat, glon, glat2, glon2)
   !--(DMK-CCATT-OLD)----------------------------------------------------------------
   !subroutine pressure_stage(n1,n2,nhem,glat,glon,glat2,glon2)
   !--(DMK-CCATT-FIM)----------------------------------------------------------------

   use isan_coms

   use ModDateUtils, only: &
      date_add_to, date_add_to_dble

   !--(DMK-CCATT-INI)----------------------------------------------------------------
   !srf-chem
   use chem_isan_coms
   use chem1_list, only: chemical_mechanism ! intent(in)
   use mem_chem1, only: CHEM_ASSIM, CHEMISTRY ! intent(in)
   !srf-chem-end
   !--(DMK-CCATT-FIM)----------------------------------------------------------------

   use rconstants

   use netcdf

   use dump, only: &
      dumpMessage

   implicit none

   include "constants.f90"
   include "netcdf.inc"

   character(len=*), parameter :: header = '**(chem_pressure_stage_netCDF)**'

   logical, external :: checkInside

   integer :: n1, n2, nhem
   real :: glat(n1, n2), glon(n1, n2), glat2(n1, n2), glon2(n1, n2)

   real :: fnprx, fnpry, grx, gry, gglat, gglon, thmax, thmin &
           , xswlon_east, cntlon_east, rr, each_levp_real(maxpr)
   integer :: i, j, k, lv, ifm, loop, n, iunit, stat
   real, external :: rs
   character(len=200) :: metadata
   character(len=32), allocatable :: varName(:)
   character(len=32) :: name, atName
   !# composed name of levels (:press:)
   character(len=300) :: grid_info
   character(len=99) :: invline
   real, allocatable :: var(:, :), lat(:, :), lon(:, :)
   real, allocatable :: llat(:), llon(:), pLev(:)
   integer :: nx, ny, thisHour, iFound, iFactor
   integer :: iyyl, imml, iddl, ihhl, an, iLen, xtype
   integer :: ndims, nvars, nglobalatts, unlimdimid, ncid
   integer, allocatable :: nat(:), varDim(:), lenDim(:)
   character(len=32), allocatable :: atValue(:)
   integer :: levelVarN, lonVarN, latVarN, uVarN, vVarN, tVarN, zVarN, qVarN
   integer :: nTimes, timeVarN, validTime
   integer, allocatable :: HoursFrom1900(:)
   real, allocatable :: latTest(:)
   real :: xvalue
   real :: scaleFactor(32), addOfset(32)
   integer :: ii, jj, kk
   integer :: iAbove(5), iBellow(5)
   real :: vMax(5), vMin(5)
   integer :: beginDate, beginTime
   character(len=33) :: timeLong
   real :: missingValue

   !--(DMK-CCATT-INI)----------------------------------------------------------------
   !srf-chem
   character(len=32) :: chemical_mechanism_test
   !srf-chem-end
   !--(DMK-CCATT-FIM)----------------------------------------------------------------

   call str2int(innpr(len_trim(innpr) - 2:len_trim(innpr)), thisHour, stat)

   !Open NetCDF file
   !write(*,fmt='(A,A,A)') '{',trim(innpr),'}'
   iErrNumber = nf90_open(path=trim(innpr), mode=nf90_nowrite, ncid=ncid)
   if (iErrNumber /= nf90_noerr) iErrNumber = dumpMessage(c_tty, c_yes, header, modelVersion &
                                                          , c_fatal, trim(innpr)//' not found. Please, check it!')

   print 91, innpr(1:len_trim(innpr))
91 format(//, ' Reading pressure gridded netCDF data', /, a, //)

   ! get info about netCDF file
   iErrNumber = nf90_inquire(ncid, ndims, nvars)
   !print *,'File open. Total of variables are ',nvars
   allocate (varName(nvars), varDim(nvars), nat(nvars), lenDim(nvars))

   ! Get dimensions in levels, lons, lats and times
   nprz_grib2 = 0
   nxGrib = 0
   nyGrib = 0
   do i = 1, nvars
      iErrNumber = nf90_Inquire_Dimension(ncid, i, name, lenDim(i))
      !print *,'Variable ',i,'name ',trim(name)
      select case (trim(name))
      case ('level')
         nprz_grib2 = lenDim(i)
         levelVarN = i
      case ('lev')
         nprz_grib2 = lenDim(i)
         levelVarN = i
      case ('longitude')
         nxGrib = lenDIm(i)
         lonVarN = i
      case ('lon')
         nxGrib = lenDIm(i)
         lonVarN = i
      case ('latitude')
         nyGrib = lenDIm(i)
         latVarN = i
      case ('lat')
         nyGrib = lenDIm(i)
         latVarN = i
      case ('time')
         nTimes = lenDIm(i)
         timeVarN = i
      end select
   end do
   !If one of variables not found is error
   if (nprz_grib2 == 0 .or. nxGrib == 0 .or. nyGrib == 0) then
      print *, 'Vars expected: level, longitude, latitude and time.'
      iErrNumber = dumpMessage(c_tty, c_yes, header, modelVersion &
                               , c_fatal, 'Var expected for dimension not found in netCDF file. Please check it!')
   end if

   write (*, fmt='("number of lats points: ",I4)') nygrib
   write (*, fmt='("number of lons points: ",I4)') nxgrib
   write (*, fmt='("number of levs       : ",I4)') nprz_grib2

   !Allocate all vars
   allocate (PressLevs(nprz_grib2))
   allocate (mask(nxGrib, nyGrib))
   allocate (llat(nyGrib), llon(nxGrib))
   allocate (HoursFrom1900(nTimes))
   do i = 1, 5
      allocate (ncVar(i)%valuesRead(nxGrib, nyGrib, nprz_grib2))
   end do

   !Get pressure levels
   iErrNumber = nf90_get_var(ncid, levelVarN, PressLevs)
   do i = 1, nprz_grib2
      levpr_grib2(i) = PressLevs(i)
   end do

   !If the total of levels ir greater than defined in RAMSIN use the defined
   if (nprz_grib2 > z_max_level) then
      nprz = z_max_level
   else
      nprz = nprz_grib2
   end if

   !Get names of vars in netCDF file and total of atributes
   do i = 1, nvars
      iErrNumber = nf90_Inquire_Variable(ncid, i, name=varName(i) &
                                         , ndims=varDim(i), nAtts=nat(i))
      !print *,i,trim(varName(i)),varDim(i),nat(i)
   end do

   !checking if required varName is present
   iFound = 0
   do i = 1, nvars
      if (trim(varName(i)) == trim(WIND_U_VARNAME)) then
         iFound = iFound + 1
         uVarN = i
      elseif (trim(varName(i)) == trim(WIND_V_VARNAME)) then
         iFound = iFound + 1
         vVarN = i
      elseif (trim(varName(i)) == trim(TEMPERATURE_VARNAME)) then
         iFound = iFound + 1
         tVarN = i
      elseif (trim(varName(i)) == trim(GEO_VARNAME)) then
         iFound = iFound + 1
         zVarN = i
      elseif (trim(varName(i)) == trim(UR_VARNAME)) then
         iFound = iFound + 1
         qVarN = i
      end if
   end do
   if (iFound /= 5) then !Not found all needs
      print *, 'Available names of var in netCDF file:'
      do i = 1, nvars
         print *, trim(varName(i))
      end do
      print *, 'Expected name os vars from RAMSIN:'
      print *, trim(WIND_U_VARNAME)
      print *, trim(WIND_V_VARNAME)
      print *, trim(TEMPERATURE_VARNAME)
      print *, trim(GEO_VARNAME)
      print *, trim(UR_VARNAME)
      iErrNumber = dumpMessage(c_tty, c_yes, header, modelVersion &
                               , c_fatal, 'Incompatible var names. See the list above and RAMSIN vars!')
   end if

   missingValue = 1.0e+20
   !Getting variables atributes
   do i = 1, nvars
      do an = 1, nat(i)
         iErrNumber = nf90_inq_attname(ncid, i, an, atName)
         iErrNumber = nf90_Inquire_Attribute(ncid, i, trim(atName), xtype=xtype, len=iLen)
         if (trim(atName) == 'scale_factor') then
            iErrNumber = nf90_get_att(ncid, i, atName, xValue)
            scaleFactor(i) = xvalue
         elseif (trim(atName) == 'add_offset') then
            iErrNumber = nf90_get_att(ncid, i, atName, xValue)
            addOfset(i) = xvalue
         elseif (trim(atName) == 'begin_date') then
            iErrNumber = nf90_get_att(ncid, i, atName, xValue)
            beginDate = xValue
         elseif (trim(atName) == 'begin_time') then
            iErrNumber = nf90_get_att(ncid, i, atName, xValue)
            beginTime = xValue
         elseif (trim(atName) == 'units' .and. i == timeVarN) then
            iErrNumber = nf90_get_att(ncid, i, atName, timeLong)
         elseif (trim(atName) == 'missing_value') then
            iErrNumber = nf90_get_att(ncid, i, atName, missingValue)
         end if
      end do
   end do

   !get lon and lat info
   iErrNumber = nf90_get_var(ncid, latVarN, llat)
   iErrNumber = nf90_get_var(ncid, lonVarN, llon)
   if (ICFILETYPE == 3) then !NASA data is from -90 to +90
      gdatdy = llat(2) - llat(1)
   else
      gdatdy = llat(1) - llat(2) !ECMWF data is from +90 to -90
   end if
   gdatdx = llon(2) - llon(1)

   write (*, fmt='("First lat point      : ",F8.2)') llat(1)
   write (*, fmt='("First lon point      : ",F8.2)') llon(1)
   write (*, fmt='("Delta lat            : ",F8.4)') gdatdy
   write (*, fmt='("Delta lon            : ",F8.4)') gdatdx
   write (*, fmt='("Missing Value        : ",E18.4)') missingValue

   nprx = ((final_longitude - initial_longitude)/gdatdx) + 1
   npry = ((final_latitude - initial_latitude)/gdatdy) + 1

   !Setting initial lon and lat
   xswlon = initial_longitude - 360.0
   xswlat = initial_latitude

   ! Check for date consistency between file parameters and namelist parameters
   iErrNumber = nf90_get_var(ncid, timeVarN, HoursFrom1900)
   if (iErrNumber /= nf90_noerr) print *, 'Err: ', iErrNumber

   !iErrNumber = nf90_get_var(ncid, latVarN, latTest)
   !if(iErrNumber /= nf90_noerr) print *, 'Err: ',iErrNumber
   !print *, "LatTest: ",size(latTest),latTest

   write (*, fmt='(A,I2,A,I2,A)') "Var time from netCDF is #", timeVarN, " with ", size(HoursFrom1900), ' times!'
   print *, HoursFrom1900, nTimes

   validTime = 0
   if (ICFILETYPE == 3) then
      call splitDate(timeLong, iyy, imm, idd, ihh)
      !print *,timelong,iyy,imm,idd,ihh
      print *, 'Time: ', iyy, imm, idd, ihh, iyear, imonth, idate, ihour
      if (iyy == iyear .and. imm == imonth .and. idd == idate &
          .and. ihh == ihour) then
         validTime = 1
      end if
   else
      do i = 1, nTimes
         write (*, fmt='(A,I2,A,I16)') "Hours for ", i, " is ", HoursFrom1900(i)
         call date_add_to_dble(1900, 1, 1, 0, dble(HoursFrom1900(i)), 'h' &
                               , iyy, imm, idd, ihh)
         if (iyy == iyear .and. imm == imonth .and. idd == idate &
             .and. ihh/100 == ihour) then
            validTime = i
            exit
         end if
      end do
   end if

   if (validTime == 0) then !Not found the date inside netCDF
      print *, 'Pressure file dates not the same as namelist!'
      print *, 'Namelist require: ', iyear, imonth, idate, ihour
      print *, 'The netCDF file has:'
      do i = 1, nTimes
         call date_add_to_dble(1900, 1, 1, 0, dble(HoursFrom1900(i)), 'h' &
                               , iyy, imm, idd, ihh)
         write (*, fmt='(I2.2,1X,I4.4,1X,I2.2,1X,I2.2,1X,I8)') i, iyy, imm, idd, ihh, HoursFrom1900(i)
      end do
      iErrNumber = dumpMessage(c_tty, c_yes, header, modelVersion &
                               , c_fatal, 'Please, check the dates!')
   end if

   !Read all upper air values for u,v,t,z and r
   do i = 1, 5
      select case (i)
      case (1)
         iErrNumber = nf90_get_var(ncid, uVarN, ncVar(i)%valuesRead &
                                   , start=(/1, 1, 1, validTime/))
         call adjustMissingValues(icFileType, i, nxGrib, nyGrib, nprz_grib2, missingValue &
                                  , levpr_grib2, ncVar(i)%valuesRead)
         ncVar(i)%valuesRead = ncVar(i)%valuesRead*scaleFactor(uVarN) &
                               + addOfset(uVarN)
      case (2)
         iErrNumber = nf90_get_var(ncid, vVarN, ncVar(i)%valuesRead &
                                   , start=(/1, 1, 1, validTime/))
         call adjustMissingValues(icFileType, i, nxGrib, nyGrib, nprz_grib2, missingValue &
                                  , levpr_grib2, ncVar(i)%valuesRead)
         ncVar(i)%valuesRead = ncVar(i)%valuesRead*scaleFactor(vVarN) &
                               + addOfset(vVarN)
      case (3)
         iErrNumber = nf90_get_var(ncid, tVarN, ncVar(i)%valuesRead &
                                   , start=(/1, 1, 1, validTime/))
         call adjustMissingValues(icFileType, i, nxGrib, nyGrib, nprz_grib2, missingValue &
                                  , levpr_grib2, ncVar(i)%valuesRead)
         ncVar(i)%valuesRead = ncVar(i)%valuesRead*scaleFactor(tVarN) &
                               + addOfset(tVarN)
      case (4)
         iErrNumber = nf90_get_var(ncid, zVarN, ncVar(i)%valuesRead &
                                   , start=(/1, 1, 1, validTime/))
         call adjustMissingValues(icFileType, i, nxGrib, nyGrib, nprz_grib2, missingValue &
                                  , levpr_grib2, ncVar(i)%valuesRead)
         ncVar(i)%valuesRead = ncVar(i)%valuesRead*scaleFactor(zVarN) &
                               + addOfset(zVarN)
      case (5)
         iErrNumber = nf90_get_var(ncid, qVarN, ncVar(i)%valuesRead &
                                   , start=(/1, 1, 1, validTime/))
         call adjustMissingValues(icFileType, i, nxGrib, nyGrib, nprz_grib2, missingValue &
                                  , levpr_grib2, ncVar(i)%valuesRead)
         ncVar(i)%valuesRead = ncVar(i)%valuesRead*scaleFactor(qVarN) &
                               + addOfset(qVarN)

      end select

      !Adjust the value accordingly RAMSIN scale_factor
      ncVar(i)%valuesRead = ncVar(i)%valuesRead*scale_factor(i)
      !If values out of valid values setup to limits
      iAbove(i) = 0
      iBellow(i) = 0
      vMax(i) = uLimit(i)
      vMin(i) = dlimit(i)
      do ii = 1, nxGrib
         do jj = 1, nyGrib
            if (.not. checkInside(llat(jj), llon(ii))) cycle
            do kk = 1, nprz_grib2
               call limitVariableWithCount(varName(i), ncVar(i)%valuesRead(ii, jj, kk) &
                                           , dlimit(i), uLimit(i), iAbove(i), iBellow(i), vMax(i), vMin(i))
            end do
         end do
      end do
      !Print information about values bellow and above limits.
      call dumpLimits(iAbove(i), iBellow(i), varName(i), uLimit(i), dLimit(i), vMin(i), vMax(i))
   end do

   if (ICFILETYPE == 3) then
      inproj = 1 !NASA lons is from -180 to +180
   else
      inproj = 1 !ECMWF lons is from 0 to 360
   end if

   ! Check pressure data domain size and location
   if (inproj /= 1 .and. nhem > 0) then
      print *, 'You must input a lat-lon pressure grid ' &
         , 'to run a global simulation !!!!'
      stop 'glob-no-press'
   end if

   if (inproj .eq. 1) then
      ! If necessary, convert longitude specification to values in the range
      ! [-180.,180.].
      xswlon = mod(xswlon + 900., 360.) - 180.
      xnelon = mod(xnelon - 900., -360.) + 180.
      if (xswlon .lt. -180. .or. xswlon .ge. 180. .or. &
          xnelon .lt. -180. .or. xnelon .gt. 180.01 .or. &
          xswlat .lt. -90. .or. xnelat .gt. 90.01) then
         print *, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         print *, '!!! BAD DOMAIN SPECIFICATION   !!'
         print *, '!!! xswlat,xswlon - ', xswlat, xswlon
         print *, '!!! xnelat,xnelon-xswlon - ', xnelat, xnelon - xswlon
         print *, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         stop 'bad_domain'
      end if
      fnprx = float(nprx)
      fnpry = float(npry)
      ! Set "global domain" flags to determine whether 2 extra rows will be
      ! added to borders of the input gridded data.  iglobew = 1 if doing
      ! full 360 degrees of longitude.  iglobs = 1 if including the south
      ! pole and if iglobew = 1, and iglobn = 1 if including
      ! the north pole and if iglobew = 1.
      iglobew = 0
      iglobn = 0
      iglobs = 0
      if (fnprx*gdatdx .gt. 359.9) then
         if (abs(xswlon + 180.) .gt. 1.e-3) then
            print *, 'When extracting all longitudes of pressure data,'
            print *, 'xswlon must be set to -180.'
            !stop 'xswlon'
         end if
         iglobew = 1
      end if
      if (iglobew .eq. 1 .and. xswlat .lt. -89.9) iglobs = 1
      if (iglobew .eq. 1 .and. xnelat .gt. 89.9) iglobn = 1
      idatelin = 0
      if (xnelon .lt. xswlon + .01) idatelin = 1
      ! Search coarse grid points for any that are outside the bounds of
      ! the pressure data.  If any are found, stop.
      loop = 1
      if (nhem > 0) loop = 2
      do ifm = 1, loop
         do j = 1, n2
            do i = 1, n1
               if (ifm == 1) then
                  gglat = glat(i, j)
                  gglon = glon(i, j)
               elseif (ifm == 2) then
                  gglat = glat2(i, j)
                  gglon = glon2(i, j)
               end if

               gry = (gglat - xswlat)/gdatdy + 1.

               if (gry .lt. 2 - iglobs .or. gry .gt. fnpry + iglobn - 1) then
                  write (*, fmt='("i      : ",I4.4)') i
                  write (*, fmt='("j      : ",I4.4)') j
                  write (*, fmt='("xswlat : ",E18.4)') xswlat
                  write (*, fmt='("xnelat : ",E18.4)') xnelat
                  write (*, fmt='("gLat   : ",E18.4)') gglat
                  write (*, fmt='("xnelat : ",E18.4)') xnelat
                  write (*, fmt='("iglobs : ",  I1)') iglobs
                  write (*, fmt='("iglobn : ",  I1)') iglobn
                  write (*, fmt='("fnpry  : ",E18.4)') fnpry
                  write (*, fmt='("gdatdy : ",E18.4)') gdatdy
                  write (*, fmt='("gry    : ",E18.4)') gry
                  if (gry .lt. 2 - iglobs) write (*, fmt='(A)') "gry is < 2-iglobs!"
                  if (gry .gt. fnpry + iglobn - 1) write (*, fmt='(A)') "gry is > fnpry+iglobn-1!"
                  iErrNumber = dumpMessage(c_tty, c_yes, header, 'Pos 3' &
                                           , c_fatal, 'Model grid point latitude must be at least 1' &
                                           //'gridpoint inside pressure data area, see info above!')
               end if

               if (idatelin .eq. 1 .and. gglon .lt. xswlon) gglon = gglon + 360.
               grx = (gglon - xswlon)/gdatdx + 1.
               if (grx .lt. 2 - iglobew .or. grx .gt. fnprx + iglobew - 1) then
                  print *, 'Model grid point longitude must be'
                  print *, 'at least 1 gridpoint inside'
                  print *, 'pressure data area'
                  print *, 'ifm,i,j,glon,xswlon,xnelon' &
                     , ifm, i, j, gglon, xswlon, xnelon, iglobew, iglobs, iglobn
                  stop 'isan 3: outside p-data 2'
               end if
            end do
         end do
      end do
   end if

   ! Deallocate memory for the pressure data
   if (allocated(p_u)) deallocate (p_u)
   if (allocated(p_v)) deallocate (p_v)
   if (allocated(p_t)) deallocate (p_t)
   if (allocated(p_z)) deallocate (p_z)
   if (allocated(p_r)) deallocate (p_r)
   if (allocated(p_ur)) deallocate (p_ur)
   if (allocated(p_vr)) deallocate (p_vr)
   if (allocated(p_lat)) deallocate (p_lat)
   if (allocated(p_lon)) deallocate (p_lon)
   if (allocated(p_slp)) deallocate (p_slp)
   if (allocated(p_sfp)) deallocate (p_sfp)
   if (allocated(p_sft)) deallocate (p_sft)
   if (allocated(p_snow)) deallocate (p_snow)
   if (allocated(p_sst)) deallocate (p_sst)
   if (allocated(p_sc)) deallocate (p_sc)

   ! Allocate memory for the pressure data
   allocate (p_u(nprx, npry, nprz))
   allocate (p_v(nprx, npry, nprz))
   allocate (p_t(nprx, npry, nprz))
   allocate (p_z(nprx, npry, nprz))
   allocate (p_r(nprx, npry, nprz))

   if (CHEM_ASSIM == 1 .and. nspecies > 0) allocate (p_sc(nprx, npry, nprz, nspecies))

   ! p_ur,p_vr arrays for rotated winds used later
   allocate (p_ur(nprx, npry, nprz))
   allocate (p_vr(nprx, npry, nprz))
   allocate (p_lat(nprx, npry))
   allocate (p_lon(nprx, npry))
   allocate (p_slp(nprx, npry))
   allocate (p_sfp(nprx, npry))
   allocate (p_sft(nprx, npry))
   allocate (p_snow(nprx, npry))
   allocate (p_sst(nprx, npry))

   ! Fill with missing in case they are not present
   p_slp(1:nprx, 1:npry) = 1.e30
   p_sfp(1:nprx, 1:npry) = 1.e30
   p_sft(1:nprx, 1:npry) = 1.e30
   p_snow(1:nprx, 1:npry) = 1.e30
   p_sst(1:nprx, 1:npry) = 1.e30

   print 192
192 format(////, 1x, 70('*')/ &
           , '  Access coarse resolution pressure data', /, 1x, 70('*')///)

   !Call routine to fill pressure arrays from the chosen dataset.
   call chem_get_press_netCDF(uVarN, vVarN, tVarN, zVarN, qVarN &
                              , llat, llon, nTimes)

   ! Find max-min theta at bottom and top levels
   thmax = 1.
   thmin = 1000.
   do j = 1, npry
      do i = 1, nprx
         if (p_t(i, j, nprz) .lt. 1.e20) &
            thmax = max(thmax, p_t(i, j, nprz)*(p00/pnpr(nprz))**rocp)
         if (p_t(i, j, 1) .lt. 1.e20) &
            thmin = min(thmin, p_t(i, j, 1)*(p00/pnpr(1))**rocp)
      end do
   end do

   print 300, levpr(1), thmin, levpr(nprz), thmax
300 format(//, ' Minimum THETA at ', I4, ' mb - ', F8.2/ &
           , ' Maximum THETA at ', I4, ' mb - ', F8.2)

   !  Compute lat-lon at input pressure data points
   select case (inproj)
   case (1)
      do j = 1, npry
         do i = 1, nprx
            p_lat(i, j) = xswlat + (j - 1)*gdatdy
            p_lon(i, j) = xswlon + (i - 1)*gdatdx
         end do
      end do
   case (2)
      print *, 'lamb-con:', xswlat, xswlon, cntlat, cntlon, gdatdx, gdatdy
      do j = 1, npry
         do i = 1, nprx
            call lc2ll(cntlat, cntlon, p_lat(i, j), p_lon(i, j) &
                       , float(i), float(j), xswlat, xswlon, gdatdx)
            !         print*,i,j,p_lat(i,j),p_lon(i,j)
         end do
      end do
   case (3)
      print *, 'polar:', cntlat, cntlon, gdatdx, gdatdy, xswlat, xswlon
      xswlon_east = xswlon
      if (xswlon < 0.) xswlon_east = xswlon + 360.
      cntlon_east = cntlon
      if (cntlon < 0.) cntlon_east = cntlon + 360.
      do j = 1, npry
         do i = 1, nprx
            call w3fb07(float(i), float(j), xswlat, xswlon_east, gdatdx, cntlon_east &
                        , p_lat(i, j), p_lon(i, j))
            !         print*,i,j,p_lat(i,j),p_lon(i,j)
         end do
      end do
   end select

   iErrNumber = nf90_close(ncid)

end subroutine chem_pressure_stage_netCDF

subroutine chem_get_press_netCDF(uVarN, vVarN, tVarN, zVarN, qVarN &
                                 , llat, llon, nTimes)
   use isan_coms
   use chem_isan_coms
   use dump, only: &
      dumpMessage

   use netcdf

   implicit none

   include "constants.f90"
   !

       character(len=*), parameter :: cFMT='(" ==  Read pressure field  ",I4," mBar for var=",A1," at ",I2,"/",I2,"/",I4,I6.4," UTC, maxval= ",E18.6,", minval= ",E18.6)'
   integer, intent(in) :: uVarN, vVarN, tVarN, zVarN, qVarN, nTimes
   real, intent(in) :: llat(nyGrib), llon(nxGrib)
   integer :: ncid

   logical, external :: checkInside

   integer :: ithere(maxpr, 5), isfthere(5)
   character*4 field, idat(5)
   data idat/'T', 'R', 'U', 'V', 'H'/
   character, parameter, dimension(5) :: varName = (/'U', 'V', 'T', 'Z', 'R'/)
   real, allocatable::as(:, :)
   character(len=200) :: metadata
   character(len=300) :: grid_info
   character(len=99) :: invline
   real, allocatable :: lat(:, :), lon(:, :)
   real, allocatable :: aux_as(:, :), aux_as_(:, :)

   integer :: i, j, k, nv, nvar, misstot, lv, n, lv_aux, auxJ
   integer :: irec, recordLen, iCount, jCount, iFactor

   !Allocate array for trimmer area
   allocate (aux_as(nprx, npry))
   allocate (aux_as_(nprx, npry))
   allocate (as(nxGrib, nyGrib))

   ! Put the upper air info using correct lats and lons and invert levels
   do lv = nprz_grib2, 1, -1
      do nvar = 1, 5
         as = ncVar(nvar)%valuesRead(:, :, lv)
         !if(nvar==3) write(55,fmt='(I2.2,1X,6(F8.2,1X),10(E18.3,1X))') lv &
         !                                                     ,llat(151),llon(216)    &
         !                                                     ,initial_latitude,final_latitude,initial_longitude,final_longitude &
         !                                                     ,as(251,161),as(216,151) &
         !                                                     ,as(251,162),as(216,152) &
         !                                                     ,as(251,163),as(216,153) &
         !                                                     ,as(251,164),as(216,154) &
         !                                                     ,as(251,165),as(216,155)
         iCount = 0
         jCount = 1
         do j = nyGrib, 1, -1
            do i = 1, nxGrib
               !Teste if i,j is inside trimmer area for regional model
               if (checkInside(llat(j), llon(i))) then
                  ! Mapping the position inside trimmer area
                  iCount = iCount + 1
                  if (iCount > nprx) then
                     iCount = 1
                     jCount = jCount + 1
                  end if
                  !Copying inside trimmer area
                  aux_as(iCount, jCount) = as(i, j)
                  !     if(nvar==3 .and. i==216 .and. j==151) write(55,fmt='(2(I2.2,1X),E18.3)') &
                  !                       icount,jcount,aux_as(iCount,jCount)
               end if
            end do
         end do

         if (ICFILETYPE == 3) then !Inverting Latitudes (NASA is from -90 to +90)
            lv_aux = lv
            do j = jCount, 1, -1
               auxJ = jCount - j + 1
               aux_as_(:, j) = aux_as(:, auxJ)
            end do
            aux_as = aux_as_
         else !Inverting pressure levels (netCDF -> BRAMS) and fill
            lv_aux = nprz_grib2 - lv + 1 !Others is from up to down
         end if
         levpr(lv_aux) = levpr_grib2(lv)
         pnpr(lv_aux) = levpr(lv_aux)*100
         !
         !Do not process levels above z_max_level
         if (lv_aux > z_max_level) cycle

         select case (nvar)
         case (1)
            call prfill(nprx, npry, aux_as, p_u(1, 1, lv_aux))
         case (2)
            call prfill(nprx, npry, aux_as, p_v(1, 1, lv_aux))
         case (3)
            call prfill(nprx, npry, aux_as, p_t(1, 1, lv_aux))
         case (4)
            call prfill(nprx, npry, aux_as, p_z(1, 1, lv_aux))
         case (5)
            call prfill(nprx, npry, aux_as, p_r(1, 1, lv_aux))
         end select

         write (*, fmt=trim(cFMT)) int(levpr_grib2(lv)), varName(nvar), imonth, idate, iyear, ihour &
            , maxval(aux_as), minval(aux_as)
      end do
   end do

   if (ccGradsWrite == 1) call writeGradsSub()

   !Fill surface fields
   do nvar = 1, 5
      as = -999.0
      select case (nvar)
      case (1)
         call prfill(nprx, npry, as, p_slp(1, 1))
      case (2)
         call prfill(nprx, npry, as, p_sfp(1, 1))
      case (3)
         call prfill(nprx, npry, as, p_sft(1, 1))
      case (4)
         call prfill(nprx, npry, as, p_snow(1, 1))
      case (5)
         call prfill(nprx, npry, as, p_sst(1, 1))
      end select
   end do

   ! Check for levels that may be all missing
   ithere(1:nprz, 1:5) = 0
   isfthere(1:5) = 0

   do k = 1, nprz
      do j = 1, npry
         do i = 1, nprx
            if (p_t(i, j, k) > 1e20) ithere(k, 1) = ithere(k, 1) + 1
            if (p_r(i, j, k) > 1e20) ithere(k, 2) = ithere(k, 2) + 1
            if (p_u(i, j, k) > 1e20) ithere(k, 3) = ithere(k, 3) + 1
            if (p_v(i, j, k) > 1e20) ithere(k, 4) = ithere(k, 4) + 1
            if (p_z(i, j, k) > 1e20) ithere(k, 5) = ithere(k, 5) + 1
         end do
      end do
   end do
   do nv = 1, 5
      do k = 1, nprz
         if (ithere(k, nv) < nprx*npry) then
            ithere(k, nv) = 1
         else
            ithere(k, nv) = 0
         end if
      end do
   end do

   misstot = 0
   do nv = 1, 5
      do k = 1, nprz
         if (ithere(k, nv) == 0) misstot = misstot + 1
      end do
   end do

   do j = 1, npry
      do i = 1, nprx
         if (p_slp(i, j) > 1e20) isfthere(1) = isfthere(1) + 1
         if (p_sfp(i, j) > 1e20) isfthere(2) = isfthere(2) + 1
         if (p_sft(i, j) > 1e20) isfthere(3) = isfthere(3) + 1
         if (p_snow(i, j) > 1e20) isfthere(4) = isfthere(4) + 1
         if (p_sst(i, j) > 1e20) isfthere(5) = isfthere(5) + 1
      end do
   end do
   do nv = 1, 5
      if (isfthere(nv) < nprx*npry) then
         isfthere(nv) = 1
      else
         isfthere(nv) = 0
      end if
   end do

   print *, '------------------------------------------------'
   print *, ' Missing parameters: 0 = all missing'
   print *, '------------------------------------------------'
   print '(t20,5(a1,6x))', (idat(n), n=1, 5)

   do k = 1, nprz
      print '(f10.1,t14,5(i7))', pnpr(k), (ithere(k, n), n=1, 5)
   end do

   print *, '------------------------------------------------'
   print *, ' Missing surface parameters: 0 = all missing'
   print *, '------------------------------------------------'
   print *, '              SLP    SFP    SFT    SNOW   SST'

   print '(t14,5(i7))', (isfthere(n), n=1, 5)

   if (misstot .gt. 0) then
      ! Let's see if we can get creative and make up data for the missing fields

      call press_miss(nprx, npry, nprz, p_u, p_v, p_t, p_z, p_r &
                      , ithere, maxpr, levpr)

      ! Check again for missing fields

      misstot = 0
      do nv = 1, 5
         do k = 1, nprz
            if (ithere(k, nv) .eq. 0) misstot = misstot + 1
         end do
      end do

      print *, '------------------------------------------------'
      print *, ' After missing parameters check: 0 = all missing'
      print *, '------------------------------------------------'
      print '(t20,5(a1,6x))', (idat(n), n=1, 5)

      do k = 1, nprz
         print '(f10.1,t14,5(i7))', pnpr(k), (ithere(k, n), n=1, 5)
      end do

   end if

   deallocate (mask)
   deallocate (pressLevs)
   do nVar = 1, 5
      deallocate (ncVar(nvar)%valuesRead)
   end do
   !write(*,fmt='(A,I4.4,2(I2.2),I4.4)') 'All done for ',iyear,imonth,idate,ihour
   return

end subroutine chem_get_press_netCDF
#endif

   !!===============================

!Utils functions and subroutines

logical function checkInside(xlat, xlon)
   use isan_coms

   real, intent(in) :: xlat, xlon
   checkInside = .false.

   if (xlat >= initial_latitude .and. xlat <= final_latitude &
       .and. xlon >= initial_longitude .and. xlon &
       <= final_longitude) checkInside = .true.

end function checkInside

subroutine limitVariableWithCount(variableName, variable, inferiorLim &
                                  , superiorLim, above, bellow, vMax, vMin)

   use dump, only: &
      dumpMessage

   include "constants.f90"
   character(len=*), parameter :: header = '**(limitVariable)**'
   real, intent(inout) :: variable
   integer, intent(inout) :: above
   integer, intent(inout) :: bellow
   real, intent(inout) :: vMax
   real, intent(inout) :: vMin
   real, intent(in) :: inferiorLim, superiorLim
   character(len=*), intent(in) :: variableName

   if (variable .lt. inferiorLim) then
      if (variable < vmin) vmin = variable
      variable = inferiorLim
      bellow = bellow + 1
   else if (variable .gt. superiorLim) then
      if (variable > vMax) vmax = variable
      variable = superiorLim
      above = above + 1
   end if

end subroutine limitVariableWithCount

subroutine dumpLimits(iAbove, iBellow, varName, uLimit, dLimit, vMin, vMax)
   !Header_vai_aqui!
   use dump, only: &
      dumpMessage
   implicit none

   include "constants.f90"

   !Parameters (constants)

   ! Input/Output variables
   integer, intent(in)    :: iAbove, iBellow
   !#
   real, intent(in)    :: uLimit, dLimit, vMin, vMax
   !#
   character(len=*), intent(in) :: varName
   !#
   !Code
   if (iAbove > 0) then
      iErrNumber = dumpMessage(c_tty, c_yes, '', '', c_warning &
                               , 'The maximum value allowed for '//trim(varName) &
                               //' is:', uLimit, 'E18.4')
      iErrNumber = dumpMessage(c_tty, c_yes, '', '', c_warning &
                               , 'There are values above the limit ' &
                               //' - total count:', iabove, 'I8')
      iErrNumber = dumpMessage(c_tty, c_yes, '', '', c_warning &
                               , 'The greater value found is:', vMax, 'E18.4')
   end if
   if (iBellow > 0) then
      iErrNumber = dumpMessage(c_tty, c_yes, '', '', c_warning &
                               , 'The minimum value allowed for '//trim(varName) &
                               //' is:', dLimit, 'E18.4')
      iErrNumber = dumpMessage(c_tty, c_yes, '', '', c_warning &
                               , 'There are values bellow the limit ' &
                               //' - total count:', iBellow, 'I8')
      iErrNumber = dumpMessage(c_tty, c_yes, '', '', c_warning &
                               , 'The minor value is:', vMin, 'E18.4')
   end if

end subroutine dumpLimits

subroutine writeGradsSub()
   !Header_vai_aqui!

   use isan_coms

   implicit none

   include "constants.f90"

   !Parameters (constants)

   ! Input/Output variables

   !Local variables
   integer :: recordLen, irec, nvar, k
   character(len=4) :: cYear, cHour
   character(len=2) :: cDay, cMonth
   character(len=256) :: gradsFileName
   character(len=19) :: fileName
   integer, external :: outRealSize

   !Code
   write (cYear, fmt='(I4.4)') iyear
   write (cMonth, fmt='(I2.2)') iMonth
   write (cDay, fmt='(I2.2)') idate
   write (cHour, fmt='(I4.4)') ihour

   fileName = 'ciREad.'//cYear//cMonth//cDay//cHour

   gradsFileName = trim(ICGRADSPREFIX)//fileName

   recordLen = outRealSize()*nprx*npry
!  recordLen=nprx*npry
   open (unit=33, file=trim(gradsFileName)//'.gra', &
         action='WRITE', status='REPLACE', form='UNFORMATTED', access='DIRECT', &
         recl=recordLen)

   !# writing grads binary and fill variables
   irec = 1
   do nvar = 1, 5
      do k = 1, nprz
         select case (nvar)
         case (1)
            write (33, rec=irec) p_u(1:nprx, 1:npry, k)
         case (2)
            write (33, rec=irec) p_v(1:nprx, 1:npry, k)
         case (3)
            write (33, rec=irec) p_t(1:nprx, 1:npry, k)
         case (4)
            write (33, rec=irec) p_z(1:nprx, 1:npry, k)
         case (5)
            write (33, rec=irec) p_r(1:nprx, 1:npry, k)
         end select
         irec = irec + 1
      end do
   end do

   close (33)

   open (unit=33, file=trim(gradsFileName)//'.ctl' &
         , action='WRITE', status='replace', form='FORMATTED')

   !writing the name of grads file
   write (33, *) 'dset ^'//fileName//'.gra'
   !writing others infos to ctl
   write (33, *) 'undef -0.9990000E+34'
   write (33, *) 'title GRIB NCEP GFS'
   write (33, *) 'xdef ', nprx, ' linear ', xswlon, gdatdx
   write (33, *) 'ydef ', npry, ' linear ', xswlat, gdatdy
   write (33, *) 'zdef ', nprz, 'levels', levpr(1:nprz)
   write (33, *) 'tdef 1 linear 00:00z01jan2018     1mo'
   write (33, *) 'vars 5'
   write (33, *) 'U ', nprz, '99 U wind'
   write (33, *) 'V ', nprz, '99 V Wind'
   write (33, *) 'T ', nprz, '99 Temp'
   write (33, *) 'Z ', nprz, '99 Geopot'
   write (33, *) 'R ', nprz, '99 Rel Umid'
   write (33, *) 'endvars'

   close (33)

end subroutine writeGradsSub

subroutine str2int(str, int, stat)
   implicit none
   ! Arguments
   character(len=*), intent(in) :: str
   integer, intent(out)         :: int
   integer, intent(out)         :: stat
   read (str, *, iostat=stat) int
end subroutine str2int

subroutine limitVariable(variableName, variable, inferiorLim, superiorLim)

   use dump, only: &
      dumpMessage

   include "constants.f90"
   character(len=*), parameter :: header = '**(limitVariable)**'
   real, intent(inout) :: variable
   real, intent(in) :: inferiorLim, superiorLim
   character(len=*), intent(in) :: variableName

   if (variable .lt. inferiorLim) then
      write (12, fmt='("Variable, Limit: ",2(F18.2,1X))') variable, inferiorLim
      iErrNumber = dumpMessage(12, c_yes, header, modelVersion &
                               , c_warning, 'Value of ' &
                               //trim(variableName)//' is inferior than limit.' &
                               //' Using limit.')
      variable = inferiorLim
   else if (variable .gt. superiorLim) then
      write (12, fmt='("Variable, Limit: ",2(F18.2,1X))') variable, superiorLim
      iErrNumber = dumpMessage(12, c_yes, header, modelVersion &
                               , c_warning, 'Value of ' &
                               //trim(variableName)//' is superior than limit.' &
                               //' Using limit.')
      variable = superiorLim
   end if

end subroutine limitVariable

subroutine splitDate(di, yy, mm, dd, hh)
   !# Subroutine to convert GEOS NetCDF date long unit to
   !# year, month, day and hour
   !#
   !# @note
   !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
   !#
   !# **Brief**: This subroutine trim the lond unit date description
   !#  minutes since 2020-02-13 00:00:00 into year, month, day and hour
   !#
   !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
   !#
   !# **Author**: Luiz Flavio Rodrigues **&#9993;**<mailto:luiz.rodrigues@inpe.br>
   !#
   !# **Date**: 2020-02-17
   !# @endnote
   !#
   !# @changes
   !#
   !# +
   !# @endchanges
   !# @bug
   !# No active bugs reported now
   !# @endbug
   !#
   !# @todo
   !#  &#9744; <br/>
   !# @endtodo
   !#
   !# @warning
   !# Now is under CC-GPL License, please see
   !# &copy; <https://creativecommons.org/licenses/GPL/2.0/legalcode.pt>
   !# @endwarning
   !#
   !#--- ----------------------------------------------------------------------------------------
   !
   character(len=*), intent(in) :: di
   !# date in long format
   integer, intent(out) :: yy
   !# year of date
   integer, intent(out) :: mm
   !# month of date
   integer, intent(out) :: dd
   !# day of date
   integer, intent(out) :: hh
   !# hour of date

   character(len=4) :: cYear
   character(len=2) :: cMonth, cDay, cHour
   integer :: stat

   !0        1         2         3
   !123456789012345678901234567890123
   !minutes since 2020-02-13 00:00:00
   !
   cYear = di(15:18)
   cMonth = di(20:21)
   cDay = di(23:24)
   cHour = di(26:27)

   !write (*,fmt='(A4,"/",A2,"/",A2," : ",A2)') cYear,cMonth,cDay,cHour
   read (cYear, *, iostat=stat) yy
   read (cMonth, *, iostat=stat) mm
   read (cDay, *, iostat=stat) dd
   read (cHour, *, iostat=stat) hh
   !write (*,fmt='(I4.4,"/",I2.2,"/",I2.2," : ",I2.2)') yy,mm,dd,hh
   hh = hh*100

end subroutine splitDate

subroutine adjustMissingValues(icFileType, nVar, nxGrib, nyGrib, nprz_grib2, missingValue &
                               , levpr_grib2, valuesRead)
   !# Subroutine to check and fix variable from input
   !#
   !# @note
   !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
   !#
   !# **Brief**: This subroutine looks for missing values inside each input variable
   !# u,v,t,z and r and fix the problem.
   !#
   !# If variable is T (Temperature) the subroutine does the replacement using an equation:
   !#
   !# \[$T_z=\frac{{P_{z+1}-P_z}}{\rho.C_p}+T_{z+1}$\]
   !#
   !# where:
   !#
   !# \[$\rho=\frac{P_z}{R.T_{z+1}}$\]
   !#
   !# otherwise it just copy the value of variable of level above for the current
   !# Z position
   !#
   !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
   !#
   !# **Author**: Luiz Flavio Rodrigues **&#9993;**<mailto:luiz.rodrigues@inpe.br>
   !#
   !# **Date**: 2020-02-19
   !# @endnote
   !#
   !# @changes
   !#
   !# +
   !# @endchanges
   !# @bug
   !# No active bugs reported now
   !# @endbug
   !#
   !# @todo
   !#  &#9744; <br/>
   !# @endtodo
   !#
   !# @warning
   !# Now is under CC-GPL License, please see
   !# &copy; <https://creativecommons.org/licenses/GPL/2.0/legalcode.pt>
   !# @endwarning
   !#
   !#--- ----------------------------------------------------------------------------------------
   !
   include "constants.f90"
   integer, parameter :: p_tVar = 3

   integer, intent(in) :: icFileType
   !# Type of input data
   integer, intent(in) :: nVar
   !# Number of incoming variable
   integer, intent(in) :: nxGrib
   !# Size of longitudes
   integer, intent(in) :: nyGrib
   !# Size of latitudes
   integer, intent(in) :: nprz_grib2
   !# Size of levels
   real, intent(in) :: levpr_grib2(nprz_grib2)
   !# Value of level pressure
   real, intent(in) :: missingValue
   !# Value for data missing
   real, intent(inout) :: valuesRead(nxGrib, nyGrib, nprz_grib2)
   !# Array with variable

   integer :: i, j, k
   !# Just counters
   real :: rho
   !# $\rho=\frac{P_z}{R.T_{z+1}}$
   real :: rhoCpInv
   !# $rhoCpInv=\frac{1.}{\rho.C_p}$
   real :: diffPress
   !# $diffPress=P_{k}-P_{k+1}$

   select case (icFileType)
   case (3)
      do k = nprz_grib2 - 1, 1, -1 !Cycle through all levels from highest to surface
         do j = 1, nyGrib !Cycle through all lats and lons
            do i = 1, nxGrib
               if (valuesRead(i, j, k) == missingValue) then !If var value is missing
                  if (nvar == p_tVar) then !If var is Temperature
                     !Equation
                     rho = levpr_grib2(k)/(c_rgas*valuesRead(i, j, k + 1))
                     rhoCpInv = 1./(rho*c_cp)
                     diffPress = levpr_grib2(k) - levpr_grib2(k + 1)
                     valuesRead(i, j, k) = valuesRead(i, j, k + 1) + diffPress*rhoCpInv
                  else !If var is U,V,Z or R
                     valuesRead(i, j, k) = valuesRead(i, j, k + 1) !Just copy the value from level above
                  end if
               end if
            end do
         end do
      end do
   case (2)
      do k = 2, nprz_grib2
         do j = 1, nyGrib
            do i = 1, nxGrib
               if (valuesRead(i, j, k) == missingValue) then
                  if (nvar == p_tVar) then
                     rho = levpr_grib2(k)/(c_rgas*valuesRead(i, j, k - 1))
                     rhoCpInv = 1./(rho*c_cp)
                     diffPress = levpr_grib2(k - 1) - levpr_grib2(k)
                     valuesRead(i, j, k) = valuesRead(i, j, k - 1) + diffPress*rhoCpInv
                  else
                     valuesRead(i, j, k) = valuesRead(i, j, k - 1)
                  end if
               end if
            end do
         end do
      end do
   end select

end subroutine adjustMissingValues

subroutine chem_pressure_stage_grads(n1, n2, n3, nhem, glat, glon, glat2, glon2)

   use isan_coms

   use ModDateUtils, only: &
      date_add_to

   use chem_isan_coms
   use chem1_list, only: chemical_mechanism ! intent(in)
   use mem_chem1, only: CHEM_ASSIM, CHEMISTRY ! intent(in)

   use rconstants

   use dump

   use mem_varinit, only: &
      nudlat

   implicit none

   include "constants.f90"

   character(len=*), parameter :: sourceName = 'chem_astp.F90' !Name of source code
   character(len=*), parameter :: procedureName = '**chem_pressure_stage_grads**' !Name of this procedure

   integer :: n1, n2, n3, nhem, nt
   real :: glat(n1, n2), glon(n1, n2), glat2(n1, n2), glon2(n1, n2)

   real :: fnprx, fnpry, grx, gry, gglat, gglon, thmax, thmin &
           , xswlon_east, cntlon_east, rr, each_levp_real(maxpr)
   integer :: i, j, k, lv, ifm, loop, n, iunit, stat
   real, allocatable :: var(:, :), lat(:, :), lon(:, :)
   integer :: nx, ny, nz, thisHour
   integer :: iyyl, imml, iddl, ihhl
   character(len=256) :: lixo, ctlFileName
   real :: lonini, dx, latini, dy
   integer, parameter :: nzmax = 256
   real :: levs(nzmax)
   character(len=3) :: cnz
   integer :: nGradsVars
   character(len=8), allocatable :: gradsVarsNames(:)
   character(len=15) :: dataCtl

   xnelat = 0.0; xnelon = 0.0; cntlat = 0.0; cntlon = 0.0; secondlat = 0.0

   ctlFileName = innpr(1:len_trim(innpr) - 3)//'ctl'

   iErrNumber = dumpMessage(c_tty, c_yes, '', '' &
                            , c_notice, 'Reading pressure gridded CTL Grads data '//innpr(1:len_trim(ctlFileName)))

   open (unit=33, file=trim(ctlFileName) &
         , action='read', status='old', form='FORMATTED')
   read (33, *) lixo
   read (33, *) lixo
   read (33, *) lixo
   read (33, *) lixo, nx, lixo, lonini, dx
   read (33, *) lixo, ny, lixo, latini, dy
   read (33, *) lixo, nz, lixo, levs(1:nz)
   read (33, *) lixo, nt, lixo, dataCtl
   read (33, *) lixo, nGradsVars
   allocate (gradsVarsNames(nGradsVars))
   do i = 1, nGradsVars
      read (33, *) gradsVarsNames(i)
      print *, 'LFR->', i, gradsVarsNames(i)
   end do

   close (33)

   allocate (slevs(nz))
   do i = 1, nz
      write (slevs(i), fmt='(F5.1)') levs(i)
      levpr_grib2(i) = levs(i)
   end do

   write (cnz, fmt='(I3.3)') nz
   write (*, fmt='(A)') ''
   write (*, fmt='(A)') '¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨'
   write (*, fmt='("Grads input data inventory for ",A)') DataCtl
   write (*, fmt='("Nx ......: ",I5.5)') nx
   write (*, fmt='("Ny ......: ",I5.5)') ny
   write (*, fmt='("Nz ......: ",I5.5)') nz
   write (*, fmt='("Init. Lon: ",F8.2,", delta Lon: ",F8.2)') lonini, dx
   write (*, fmt='("Init. Lat: ",F8.2,", delta Lat: ",F8.2)') latini, dy
   write (*, fmt='("Final Lon: ",F8.2)') lonini + (nx - 1)*dx
   write (*, fmt='("Final Lat: ",F8.2)') latini + (ny - 1)*dy
   write (*, fmt='("Levels ..: "'//cnz//'(F6.1,1X))') (levpr_grib2(i), i=1, nz)
   write (*, fmt='(A)') ''
   write (*, fmt='(A)') '¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨'
   write (*, fmt='("Model limits inventory",A)') ''
   write (*, fmt='("Nx ......: ",I5.5)') n1
   write (*, fmt='("Ny ......: ",I5.5)') n2
   write (*, fmt='("Nz ......: ",I5.5)') n3
   write (*, fmt='("Init. Lon: ",F8.2,", delta Lon: ",F8.2)') glon(1, 1), glon(2, 1) - glon(1, 1)
   write (*, fmt='("Init. Lat: ",F8.2,", delta Lat: ",F8.2)') glat(1, 1), glat(1, 2) - glat(1, 1)
   write (*, fmt='("Final Lon: ",F8.2)') glon(n1, 1)
   write (*, fmt='("Final Lat: ",F8.2)') glat(1, n2)
   write (*, fmt='(A)') '¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨'
   write (*, fmt='(A)') c_pearlWhite
   write (*, fmt='("Model MINIMUM limits with lateral nuding",A)') ''
   write (*, fmt='("NudLat   : ",I5.5)') nudlat
   write (*, fmt='("Nx+NudLat: ",I5.5)') n1 + nudlat
   write (*, fmt='("Ny+NudLat: ",I5.5)') n2 + nudlat
   write (*, fmt='("Nz ......: ",I5.5)') n3
   write (*, fmt='("Init. Lon: ",F8.2,", delta Lon: ",F8.2)') glon(1, 1) - nudlat*(glon(2, 1) - glon(1, 1)), glon(2, 1) - glon(1, 1)
   write (*, fmt='("Init. Lat: ",F8.2,", delta Lat: ",F8.2)') glat(1, 1) - nudlat*(glat(1, 2) - glat(1, 1)), glat(1, 2) - glat(1, 1)
   write (*, fmt='("Final Lon: ",F8.2)') glon(n1, 1) + nudlat*(glon(2, 1) - glon(1, 1))
   write (*, fmt='("Final Lat: ",F8.2)') glat(1, n2) + nudlat*(glat(1, 2) - glat(1, 1))
   write (*, fmt='(A)') 'c_noColor'
   if (glon(1, 1) < lonini .or. glon(n1, 1) > lonini + (nx - 1)*dx &
       .or. glat(1, 1) < latini .or. glat(1, n2) > latini + (ny - 1)*dy) then
      iErrNumber = dumpMessage(c_tty, c_yes, sourceName, procedureName &
                               , c_fatal, 'Bad domain specification, see limits above!')
   end if

   !If the total of levels ir greater than defined in RAMSIN use the defined
   if (nz > z_max_level) then
      nprz = z_max_level
   else
      nprz = nz
   end if
   nprz_grib2 = nz
   nxGrib = nx
   nyGrib = ny
   !
   ! Increments lat and lon
   gdatdy = dy
   gdatdx = dx

   !Make number of points lon and lat accordingly setup in RAMSIN
   nprx = nx
   npry = ny

   !Setting initial lon and lat
   xswlon = lonini
   xswlat = latini

   inproj = 1

   if (inproj .eq. 1) then
      ! If necessary, convert longitude specification to values in the range
      ! [-180.,180.].
      xswlon = mod(xswlon + 900., 360.) - 180.
      xnelon = mod(xnelon - 900., -360.) + 180.
      if (xswlon < -180.00 .or. xswlon > 180.00 &
          .or. xnelon < -180.00 .or. xnelon > 180.01 &
          .or. xswlat < -90.00 .or. xnelat > 90.01) then
         write (*, fmt='("xswlon: ",F8.2)') xswlon
         write (*, fmt='("xnelon: ",F8.2)') xnelon
         write (*, fmt='("xswlat: ",F8.2)') xswlat
         iErrNumber = dumpMessage(c_tty, c_yes, sourceName, procedureName &
                                  , c_fatal, 'Bad domain specification, see info above!')
      end if

      fnprx = float(nprx)
      fnpry = float(npry)
      ! Set "global domain" flags to determine whether 2 extra rows will be
      ! added to borders of the input gridded data.  iglobew = 1 if doing
      ! full 360 degrees of longitude.  iglobs = 1 if including the south
      ! pole and if iglobew = 1, and iglobn = 1 if including
      ! the north pole and if iglobew = 1.
      iglobew = 0
      iglobn = 0
      iglobs = 0
      if (fnprx*gdatdx > 359.9) then
         if (abs(xswlon + 180.0) < 1.e-3) then
            print *, 'When extracting all longitudes of pressure data,'
            print *, 'xswlon must be set to -180.'
         end if
         iglobew = 1
      end if

      if (iglobew == 1 .and. xswlat < -89.9) iglobs = 1
      if (iglobew == 1 .and. xnelat > 89.9) iglobn = 1

      idatelin = 0
      if (xnelon .lt. xswlon + .01) idatelin = 1
      ! Search coarse grid points for any that are outside the bounds of
      ! the pressure data.  If any are found, stop.
      loop = 1
      if (nhem > 0) loop = 2
      do ifm = 1, loop
         do j = 1, n2
            do i = 1, n1
               if (ifm == 1) then
                  gglat = glat(i, j)
                  gglon = glon(i, j)
               elseif (ifm == 2) then
                  gglat = glat2(i, j)
                  gglon = glon2(i, j)
               end if
               gry = (gglat - xswlat)/gdatdy + 1.
               if (gry < 2 - iglobs .or. gry > fnpry + iglobn - 1) then
                  write (*, fmt='("i      : ",I4.4)') i
                  write (*, fmt='("j      : ",I4.4)') j
                  write (*, fmt='("xswlat : ",F8.2)') xswlat
                  write (*, fmt='("xnelat : ",F8.2)') xnelat
                  write (*, fmt='("gLat   : ",F8.2)') gglat
                  write (*, fmt='("iglobs : ",  I1)') iglobs
                  write (*, fmt='("iglobn : ",  I1)') iglobn
                  write (*, fmt='("fnpry  : ",F8.2)') fnpry
                  write (*, fmt='("gdatdy : ",F8.2)') gdatdy
                  write (*, fmt='("gry    : ",F8.2)') gry
                  if (gry .lt. 2 - iglobs) write (*, fmt='(A)') "gry is < 2-iglobs!"
                  if (gry .gt. fnpry + iglobn - 1) write (*, fmt='(A)') "gry is > fnpry+iglobn-1!"
                  iErrNumber = dumpMessage(c_tty, c_yes, sourceName, procedureName &
                                           , c_fatal, 'Model grid point latitude must be at least 1 ' &
                                           //'gridpoint inside pressure data area, see info above!')
               end if

               if (idatelin .eq. 1 .and. gglon .lt. xswlon) gglon = gglon + 360.
               grx = (gglon - xswlon)/gdatdx + 1.
               if (grx < 2 - iglobew .or. grx > fnprx + iglobew - 1) then
                  write (*, fmt='("i      : ",I4.4)') i
                  write (*, fmt='("j      : ",I4.4)') j
                  write (*, fmt='("xswlon : ",F8.2)') xswlon
                  write (*, fmt='("xnelon : ",F8.2)') xnelon
                  write (*, fmt='("gLon   : ",F8.2)') gglon
                  write (*, fmt='("iglobew: ",  I1)') iglobew
                  write (*, fmt='("iglobn : ",  I1)') iglobn
                  write (*, fmt='("fnprx  : ",F8.2)') fnprx
                  write (*, fmt='("gdatdx : ",F8.2)') gdatdx
                  write (*, fmt='("grx    : ",F8.2)') grx
                  if (grx < 2 - iglobew) write (*, fmt='(A)') "grx is < 2-iglobew!"
                  if (grx > fnprx + iglobew - 1) write (*, fmt='(A)') "grx is > fnprx+iglobew-1!"
                  iErrNumber = dumpMessage(c_tty, c_yes, sourceName, procedureName &
                                           , c_fatal, 'Model grid point longitude must be at least 1 ' &
                                           //'gridpoint inside pressure data area, see info above!')
               end if
            end do
         end do
      end do
   end if

   ! Allocate memory for the pressure data
   if (.not. allocated(p_u)) allocate (p_u(nprx, npry, nprz))
   if (.not. allocated(p_v)) allocate (p_v(nprx, npry, nprz))
   if (.not. allocated(p_t)) allocate (p_t(nprx, npry, nprz))
   if (.not. allocated(p_z)) allocate (p_z(nprx, npry, nprz))
   if (.not. allocated(p_r)) allocate (p_r(nprx, npry, nprz))

   !print *,'LFR: ',nprx,npry,nprz,nspecies_aer_in
   if (CHEM_ASSIM == 1 .and. nspecies > 0 .and. .not. allocated(p_sc)) allocate (p_sc(nprx, npry, nprz, nspecies))
   if(CHEM_ASSIM == 1 .and. nspecies_aer_in>0 .and. .not. allocated(p_aer_sc)) allocate(p_aer_sc(nprx,npry,nprz,nspecies_aer_in))

   ! p_ur,p_vr arrays for rotated winds used later
   if (.not. allocated(p_ur)) allocate (p_ur(nprx, npry, nprz))
   if (.not. allocated(p_vr)) allocate (p_vr(nprx, npry, nprz))
   if (.not. allocated(p_lat)) allocate (p_lat(nprx, npry))
   if (.not. allocated(p_lon)) allocate (p_lon(nprx, npry))
   if (.not. allocated(p_slp)) allocate (p_slp(nprx, npry))
   if (.not. allocated(p_sfp)) allocate (p_sfp(nprx, npry))
   if (.not. allocated(p_sft)) allocate (p_sft(nprx, npry))
   if (.not. allocated(p_snow)) allocate (p_snow(nprx, npry))
   if (.not. allocated(p_sst)) allocate (p_sst(nprx, npry))

   ! Fill with missing in case they are not present
   p_slp(1:nprx, 1:npry) = 1.e30
   p_sfp(1:nprx, 1:npry) = 1.e30
   p_sft(1:nprx, 1:npry) = 1.e30
   p_snow(1:nprx, 1:npry) = 1.e30
   p_sst(1:nprx, 1:npry) = 1.e30

   iErrNumber = dumpMessage(c_tty, c_yes, '', '' &
                            , c_notice, 'Access coarse resolution pressure data')
   write (*, *) ''

   ! Call routine to fill pressure arrays from the chosen dataset.
   call chem_get_press_grads(lonini, dx, latini, dy, nGradsVars, gradsVarsNames)

   ! Find max-min theta at bottom and top levels
   thmax = 1.
   thmin = 1000.
   do j = 1, npry
      do i = 1, nprx
         if (p_t(i, j, nprz) .lt. 1.e20) &
            thmax = max(thmax, p_t(i, j, nprz)*(p00/pnpr(nprz))**rocp)
         if (p_t(i, j, 1) .lt. 1.e20) &
            thmin = min(thmin, p_t(i, j, 1)*(p00/pnpr(1))**rocp)
      end do
   end do

   print 300, levpr(1), thmin, levpr(nprz), thmax
300 format(//, ' Minimum THETA at ', I4, ' mb - ', F8.2/ &
           , ' Maximum THETA at ', I4, ' mb - ', F8.2)

   !  Compute lat-lon at input pressure data points
   select case (inproj)
   case (1)
      do j = 1, npry
         do i = 1, nprx
            p_lat(i, j) = xswlat + (j - 1)*gdatdy
            p_lon(i, j) = xswlon + (i - 1)*gdatdx
         end do
      end do
   case (2)
      print *, 'lamb-con:', xswlat, xswlon, cntlat, cntlon, gdatdx, gdatdy
      do j = 1, npry
         do i = 1, nprx
            call lc2ll(cntlat, cntlon, p_lat(i, j), p_lon(i, j) &
                       , float(i), float(j), xswlat, xswlon, gdatdx)
            !         print*,i,j,p_lat(i,j),p_lon(i,j)
         end do
      end do
   case (3)
      print *, 'polar:', cntlat, cntlon, gdatdx, gdatdy, xswlat, xswlon
      xswlon_east = xswlon
      if (xswlon < 0.) xswlon_east = xswlon + 360.
      cntlon_east = cntlon
      if (cntlon < 0.) cntlon_east = cntlon + 360.
      do j = 1, npry
         do i = 1, nprx
            call w3fb07(float(i), float(j), xswlat, xswlon_east, gdatdx, cntlon_east &
                        , p_lat(i, j), p_lon(i, j))
            !         print*,i,j,p_lat(i,j),p_lon(i,j)
         end do
      end do
   end select

end subroutine chem_pressure_stage_grads

subroutine chem_get_press_grads(lonini, dx, latini, dy, nGradsVars, gradsVarsNames)
   use isan_coms
   use chem_isan_coms
#ifdef GRIB2
   use wgrib2api
#endif
   use dump, only: &
      dumpMessage

   implicit none

   include "constants.f90"

   real, intent(in) :: lonini, dx, latini, dy
   !
   logical, external :: checkInside
   integer, intent(in) :: nGradsVars
   character(len=8), intent(in) :: gradsVarsNames(nGradsVars)

   integer :: ithere(maxpr, 5), isfthere(5)
   character*4 field, idat(5)
   data idat/'T', 'R', 'U', 'V', 'H'/
   character(len=60) :: varName(5)
   real, allocatable::as(:, :)
   real, allocatable :: lat(:, :), lon(:, :)
   real, allocatable :: aux_as(:, :)

   integer :: i, j, k, nv, nvar, misstot, lv, n, lv_aux
   integer :: irec, recordLen, iCount, jCount
   integer :: iAbove(5), iBellow(5)
   real :: vMax(5), vMin(5)
   real, allocatable :: dados(:, :, :, :)
   character(len=256) :: message
   integer, external :: outRealSize

   !Allocate array for trimmer area
   allocate (aux_as(nprx, npry))
   iAbove = 0
   iBellow = 0
   vMax = uLimit
   vMin = dlimit

   allocate (lat(nxGrib, nyGrib), lon(nxGrib, nyGrib), as(nxGrib, nyGrib))
   allocate (dados(5 + nspecies + nspecies_aer_in, nxGrib, nyGrib, nprz_grib2))

   !Making lat and lon arrays
   do i = 1, nxGrib
      lon(i, :) = lonini + (i - 1)*dx
      !write (88,*) 'Lon:',i,lon(i,1)
   end do
   do j = 1, nyGrib
      lat(:, j) = latini + (j - 1)*dy
      !write (89,*) "Lat:",j,lat(1,j)
   end do

   !print *,'INITIAL_LATITUDE: ',INITIAL_LATITUDE,minval(lat)
   !print *,'FINAL_LATITUDE: ',FINAL_LATITUDE,maxval(lat)
   !print *,'INITIAL_LONGITUDE:',INITIAL_LONGITUDE,minval(lon)
   !print *,'FINAL_LONGITUDE: ',FINAL_LONGITUDE,maxval(lon)

   !  Read upper air fields
   recordLen = outRealSize()*nxGrib*nyGrib
!srf  recordLen=outRealSize()*nxGrib*nyGrib
   open (unit=33, file=trim(innpr), &
         action='read', status='old', form='UNFORMATTED', access='DIRECT', &
         recl=recordLen)

   !# reading grads binary and fill variables
   irec = 1
   do nvar = 1, 5 + nspecies + nspecies_aer_in
      do k = 1, nprz_grib2
         read (33, rec=irec) dados(nvar, 1:nxGrib, 1:nyGrib, k)
         irec = irec + 1
      end do
   end do

   close (33)

   do lv = 1, nprz_grib2

      levpr(lv) = levpr_grib2(lv)
      lv_aux = lv
      pnpr(lv_aux) = levpr_grib2(lv)*100
      !Do not process levels above z_max_level
      !    if(lv_aux>z_max_level) cycle
      do nvar = 1, 5 + nspecies + nspecies_aer_in

         if (nvar <= 5) then
            aux_as = dados(nvar, :, :, lv)*scale_factor(nvar)
         else
            aux_as = dados(nvar, :, :, lv)
         end if

         write (message, fmt='("Reading pressure field ",I4," mBar for ",A8)') &
            int(levpr_grib2(lv)), gradsVarsNames(nvar)
         iErrNumber = dumpMessage(c_tty, c_yes, '', '' &
                                  , c_notice, trim(message))

         if (nvar .eq. 1) then
            call prfill(nprx, npry, aux_as, p_u(1, 1, lv_aux))
         elseif (nvar .eq. 2) then
            call prfill(nprx, npry, aux_as, p_v(1, 1, lv_aux))
         elseif (nvar .eq. 3) then
            call prfill(nprx, npry, aux_as, p_t(1, 1, lv_aux))
         elseif (nvar .eq. 4) then
            call prfill(nprx, npry, aux_as, p_z(1, 1, lv_aux))
         elseif (nvar .eq. 5) then
            call prfill(nprx, npry, aux_as, p_r(1, 1, lv_aux))
         elseif (nvar > 5 .and. nvar <= 5 + nspecies) then
            print *, 'LFR-Chem>nvar,ini', nvar, nvar - 5
            call prfill(nprx, npry, aux_as, p_sc(1, 1, lv, nvar - 5))
         elseif (nvar > 5 + nspecies) then
            print *, 'LFR-Aer>nvar,ini', nvar, nvar - 5 - nspecies
            call prfill(nprx, npry, aux_as, p_aer_sc(1, 1, lv, nvar - 5 - nspecies))
         end if
      end do
   end do

   if (ccGradsWrite == 1) call writeGradsSub4Grads(nGradsVars, gradsVarsNames)

   !Fill surface fields
   do nvar = 1, 5
      as = -999.0

      select case (nvar)
      case (1)
         call prfill(nprx, npry, as, p_slp(1, 1))
      case (2)
         call prfill(nprx, npry, as, p_sfp(1, 1))
      case (3)
         call prfill(nprx, npry, as, p_sft(1, 1))
      case (4)
         call prfill(nprx, npry, as, p_snow(1, 1))
      case (5)
         call prfill(nprx, npry, as, p_sst(1, 1))
      end select
   end do

   ! Check for levels that may be all missing
   ithere(1:nprz, 1:5) = 0
   isfthere(1:5) = 0

   do k = 1, nprz
      do j = 1, npry
         do i = 1, nprx
            if (p_t(i, j, k) > 1e20) ithere(k, 1) = ithere(k, 1) + 1
            if (p_r(i, j, k) > 1e20) ithere(k, 2) = ithere(k, 2) + 1
            if (p_u(i, j, k) > 1e20) ithere(k, 3) = ithere(k, 3) + 1
            if (p_v(i, j, k) > 1e20) ithere(k, 4) = ithere(k, 4) + 1
            if (p_z(i, j, k) > 1e20) ithere(k, 5) = ithere(k, 5) + 1
         end do
      end do
   end do
   do nv = 1, 5
      do k = 1, nprz
         if (ithere(k, nv) < nprx*npry) then
            ithere(k, nv) = 1
         else
            ithere(k, nv) = 0
         end if
      end do
   end do

   misstot = 0
   do nv = 1, 5
      do k = 1, nprz
         if (ithere(k, nv) == 0) misstot = misstot + 1
      end do
   end do

   do j = 1, npry
      do i = 1, nprx
         if (p_slp(i, j) > 1e20) isfthere(1) = isfthere(1) + 1
         if (p_sfp(i, j) > 1e20) isfthere(2) = isfthere(2) + 1
         if (p_sft(i, j) > 1e20) isfthere(3) = isfthere(3) + 1
         if (p_snow(i, j) > 1e20) isfthere(4) = isfthere(4) + 1
         if (p_sst(i, j) > 1e20) isfthere(5) = isfthere(5) + 1
      end do
   end do
   do nv = 1, 5
      if (isfthere(nv) < nprx*npry) then
         isfthere(nv) = 1
      else
         isfthere(nv) = 0
      end if
   end do

   print *, '------------------------------------------------'
   print *, ' Missing parameters: 0 = all missing'
   print *, '------------------------------------------------'
   print '(t20,5(a1,6x))', (idat(n), n=1, 5)

   do k = 1, nprz
      print '(f10.1,t14,5(i7))', pnpr(k), (ithere(k, n), n=1, 5)
   end do

   print *, '------------------------------------------------'
   print *, ' Missing surface parameters: 0 = all missing'
   print *, '------------------------------------------------'
   print *, '              SLP    SFP    SFT    SNOW   SST'

   print '(t14,5(i7))', (isfthere(n), n=1, 5)

   if (misstot .gt. 0) then
      ! Let's see if we can get creative and make up data for the missing fields

      call press_miss(nprx, npry, nprz, p_u, p_v, p_t, p_z, p_r &
                      , ithere, maxpr, levpr)

      ! Check again for missing fields

      misstot = 0
      do nv = 1, 5
         do k = 1, nprz
            if (ithere(k, nv) .eq. 0) misstot = misstot + 1
         end do
      end do

      print *, '------------------------------------------------'
      print *, ' After missing parameters check: 0 = all missing'
      print *, '------------------------------------------------'
      print '(t20,5(a1,6x))', (idat(n), n=1, 5)

      do k = 1, nprz
         print '(f10.1,t14,5(i7))', pnpr(k), (ithere(k, n), n=1, 5)
      end do

   end if

   deallocate (slevs)
   !deallocate(mask)
   !write(*,fmt='(A,I4.4,2(I2.2),I4.4)') 'All done for ',iyear,imonth,idate,ihour
   return

   !--(DMK-CCATT-INI)----------------------------------------------------------------
end subroutine chem_get_press_grads

subroutine writeGradsSub4Grads(nGradsVars, gradsVarsNames)
   !Header_vai_aqui!

   use isan_coms
   use chem_isan_coms

   implicit none

   include "constants.f90"

   !Parameters (constants)

   ! Input/Output variables
   integer, intent(in) :: nGradsVars !,nspecies
   character(len=8) :: gradsVarsNames(nGradsVars)

   !Local variables
   integer :: recordLen, irec, nvar, k
   character(len=4) :: cYear, cHour
   character(len=2) :: cDay, cMonth
   character(len=256) :: gradsFileName
   character(len=19) :: fileName
   integer, external :: outRealSize

   real :: cs(nprx, npry)

   !Code
   write (cYear, fmt='(I4.4)') iyear
   write (cMonth, fmt='(I2.2)') iMonth
   write (cDay, fmt='(I2.2)') idate
   write (cHour, fmt='(I4.4)') ihour

   fileName = 'ci.'//cYear//cMonth//cDay//cHour

   gradsFileName = trim(ICGRADSPREFIX)//fileName

   recordLen = outRealSize()*nprx*npry
   open (unit=33, file=trim(gradsFileName)//'.gra', &
         action='WRITE', status='REPLACE', form='UNFORMATTED', access='DIRECT', &
         recl=recordLen)

   !# writing grads binary and fill variables
   irec = 1
   do nvar = 1, 10
      do k = 1, nprz
         select case (nvar)
         case (1)
            write (33, rec=irec) p_u(1:nprx, 1:npry, k)
         case (2)
            write (33, rec=irec) p_v(1:nprx, 1:npry, k)
         case (3)
            write (33, rec=irec) p_t(1:nprx, 1:npry, k)
         case (4)
            write (33, rec=irec) p_z(1:nprx, 1:npry, k)
         case (5)
            write (33, rec=irec) p_r(1:nprx, 1:npry, k)
         case (6)
            write (33, rec=irec) p_slp(1:nprx, 1:npry)
         case (7)
            write (33, rec=irec) p_sfp(1:nprx, 1:npry)
         case (8)
            write (33, rec=irec) p_sft(1:nprx, 1:npry)
         case (9)
            write (33, rec=irec) p_snow(1:nprx, 1:npry)
         case (10)
            write (33, rec=irec) p_sst(1:nprx, 1:npry)
         end select
         irec = irec + 1
      end do
   end do
   !do nvar=1,nspecies
   !    do k=1,nprz
   !        cs=p_sc(:,:,k,nvar)
   !        write (33,rec=irec)  cs
   !        irec=irec+1
   !    enddo
   !enddo
   !do nvar=1,nspecies_aer_in
   !    do k=1,nprz
   !        cs=p_aer_sc(:,:,k,nvar)
   !        write (33,rec=irec)  cs
   !        irec=irec+1
   !    enddo
   !enddo

   close (33)

   open (unit=33, file=trim(gradsFileName)//'.ctl' &
         , action='WRITE', status='replace', form='FORMATTED')

   !writing the name of grads file
   write (33, *) 'dset ^./'//trim(fileName)//'.gra'
   !writing others infos to ctl
   write (33, *) 'undef -0.9990000E+34'
   write (33, *) 'title GRIB NCEP GFS'
   write (33, *) 'xdef ', nprx, ' linear ', xswlon, gdatdx
   write (33, *) 'ydef ', npry, ' linear ', xswlat, gdatdy
   write (33, *) 'zdef ', nprz, 'levels', levpr(1:nprz)
   write (33, *) 'tdef 1 linear 00:00z01jan2018     1mo'
   write (33, *) 'vars ', 10
   do nvar = 1, 5
      write (33, *) gradsVarsNames(nvar), nprz, '99 ', gradsVarsNames(nvar)
   end do
   write (33, *) 'slp', nprz, '99 ', 'slp'
   write (33, *) 'sfp', nprz, '99 ', 'sfp'
   write (33, *) 'sft', nprz, '99 ', 'sft'
   write (33, *) 'snow', nprz, '99 ', 'snow'
   write (33, *) 'sst', nprz, '99 ', 'sst'

   write (33, *) 'endvars'

   close (33)

end subroutine writeGradsSub4Grads
