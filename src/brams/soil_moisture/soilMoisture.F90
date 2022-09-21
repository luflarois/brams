module soilMoisture
   !! Soil Moisture Initialization for NWP Models
   !!
   !! @note
   !!
   !! **Project**: BRAMS
   !! **Author(s)**: Freitas, S. R. [SRF], Gevaerd, R. [GVR]
   !! **e-mail**: <mailto:saulo.freitas@inpe.br>
   !! **Date**:  2007
   !!
   !! **Full description**:
   !! Soil Moisture Initialization for NWP Models
   !! Coded and implemented by Rodrigo Gevaerd and Saulo Freitas
   !! Ref.: Gevaerd, R. e S. R. Freitas, Estimativa operacional da umidade
   !! do solo para iniciacao de modelos de previsao numerica da atmosfera.
   !! Parte I: Descricao da metodologia e validacao. Rev. Bras. Meteo.,
   !! volume especial do LBA, 2007.
   !! @endnote
   !!
   !! @warning
   !!
   !!  [](https://www.gnu.org/graphics/gplv3-127x51.png'')
   !!
   !!     This program is free software: you can redistribute it and/or modify
   !!     it under the terms of the GNU General Public License as published by
   !!     the  Free  Software  Foundation, either version 3 of the License, or
   !!     (at your option) any later version.
   !!
   !!     This program is distributed in the hope that it  will be useful, but
   !!     WITHOUT  ANY  WARRANTY;  without  even  the   implied   warranty  of
   !!     MERCHANTABILITY or FITNESS FOR A  PARTICULAR PURPOSE.  See  the, GNU
   !!     GNU General Public License for more details.
   !!
   !!     You should have received a copy  of the GNU General  Public  License
   !!     along with this program.  If not, see <https://www.gnu.org/licenses/>.
   !!
   !! @endwarning

   use dump

   use ModNamelistFile, only: namelistFile

   use memSoilMoisture, only: &
      soil_moist, soil_moist_fail, usdata_in, usmodel_in ! INTENT(IN)

   implicit none
   include 'constants.f90'

   character(len=*), parameter :: sourceName = 'soilMoisture.F90' ! Nome do arquivo fonte

   public  :: soilMoistureInit, StoreNamelistFileAtSoilMoisture
   real ::  factorsm = 1.0

!!$  PRIVATE :: gatherData
!!$  PRIVATE :: gatherData2D
!!$  PRIVATE :: gatherData4D
!  PUBLIC :: apiPrlatlon
!  PUBLIC :: interpolacao
!  PRIVATE :: changeDay

!!$  interface gatherData
!!$     module procedure gatherData2D, gatherData4D
!!$  end interface

contains

   subroutine soilmoistureinit(n1, n2, n3, mzg, mzs, npat, ifm,   &
      theta, pi0, pp,                                           &
      soil_water, soil_energy, soil_text,                       &
      glat, glon,                                               &
      lpw_r,seatp,seatf                                         &
      )
      !# _
      !#
      !# @note
      !#
      !#
      !# **Brief**: _
      !#
      !# **Documentation**: <LF>
      !#
      !# **Author**: Luiz Flavio **&#9993;**<mailto:lufla@enterprise>
      !#
      !# **Date**: 2020-04-30
      !# @endnote
      !#
      !# @changes
      !#**Changelogs:**
      !# +
      !# @endchanges
      !# @bug
      !# **Open bugs:**
      !#
      !# @endbug
      !#
      !# @todo
      !# **Todo list:**n    !#
      !# @endtodo
      !#
      !# @warning
      !# Now is under CC-GPL License, please see
      !# &copy; <https://creativecommons.org/licenses/GPL/2.0/legalcode.pt>
      !# @endwarning
      !#
      !#--- ----------------------------------------------------------------------------------------
      use dump
      use mem_grid, only : runtype,  &          ! intent(in)
         iyear1, imonth1, idate1, itime1, &   ! intent(in)
         nnxp, nnyp, nnzp,                &   ! intent(in)
         globalsizes                          ! subroutine

      use io_params, only : timstr     ! intent(in)

      use rconstants, only : cpi       ! intent(in)

      use leaf_coms, only : soilcp,  & ! intent(in)
         slmsts,                   & ! intent(in)
         slcpd                       ! intent(in)

      use mem_leaf, only : stgoff,   & ! intent(in)
         slmstr,                   & ! intent(in)
         slz                         ! intent(in)

      use node_mod, only: &
         nodei0, nodej0, & ! intent(in)
         nodemxp, nodemyp, & ! intent(in)
         nmachs,  & ! intent(in)
         mynum,  &  ! intent(in)
         mchnum, &  ! intent(in)
         master_num ! intent(in)

      !use parlib, only: &
      use parlib, only: &
         parf_bcast ! subroutine

      use mem_aerad, only: &
         nwave ! intent(in)

      use readbcst, only: &
         gatherdata

      !!!!!!DSM {
#ifdef cdf

      use netcdf, only: nf90_nowrite, nf90_open, nf90_put_att

#endif
      !!!!!!DSM}

      implicit none

      character(len=*),parameter :: procedureName="soilmoistureinit"
      character(len=*),parameter :: srcName="soilMoisture.f90"
        include "constants.f90"
        include "i8.h"
        include "files.h"

      ! arguments:
      integer, intent(in) :: n1, n2, n3, mzg, mzs, npat, ifm
      real, intent(in)    :: theta(:,:,:)         !(n1,n2,n3)
      real, intent(in)    :: pi0(:,:,:)           !(n1,n2,n3)
      real, intent(in)    :: pp(:,:,:)            !(n1,n2,n3)
      real, intent(in)    :: glat(:,:)            !(n2,n3)
      real, intent(in)    :: glon(:,:)            !(n2,n3)
      real, intent(in)    :: lpw_r(:,:)           !(n2,n3)
      real, intent(inout) :: soil_water(:,:,:,:)  !(mzg,n2,n3,npat)
      real, intent(inout) :: soil_energy(:,:,:,:) !(mzg,n2,n3,npat)
      real, intent(in)    :: soil_text(:,:,:,:)   !(mzg,n2,n3,npat)

      real, intent(in)    :: seatp(:,:) ,seatf(:,:)

      integer :: lpw(n2,n3)             !(n2,n3)
      ! local variables:

      integer            :: i, j, k, ipat, nveg, nsoil,it
      real               :: c1, airtemp, pis, soiltemp
      real               :: globalsoilwater (mzg, nnxp(ifm), nnyp(ifm), npat)
      real               :: globalsoilenergy(mzg, nnxp(ifm), nnyp(ifm), npat)
      real               :: globalsoiltext(mzg, nnxp(ifm), nnyp(ifm), npat)
      real               :: globalglon(nnxp(ifm), nnyp(ifm))
      real               :: globalglat(nnxp(ifm), nnyp(ifm))
      integer            :: qi1, qi2, qj1, qj2, ncount,               &
         ii, jj,iii,jjj,i_ini,j_ini,i_fim,j_fim,skip,range,is,js, jc, ic, i1, j1, i2, j2, kk, ifname, k2, &
         ipref, ipref_start, icihourmin
      integer            :: n4us
      integer            :: nlat, nlon, STATUSi,ncid,var_id
      integer            :: nf_INQ_DIMLEN,nf_inq_varid,nf_get_var,nf_close &
         ,nf_def_var,nf_enddef,nf_put_var, nf_INQ_DIMID
      real, allocatable  :: slz_us(:)
      real               :: latni, latnf, lonni, lonnf, ilatn, ilonn, &
         ilats, ilons, latn, lonn, lats, lons, dlatr, dlonr
      logical            :: there, theref
      character(len=f_name_length) :: usdata, usmodel
      character(len=50)  :: pref
      character(len=2)   :: cidate, cimon
      character(len=1)   :: cgrid
      character(len=4)   :: ciyear
      character(len=4)   :: cihourmin
      real, allocatable  :: api_us(:,:,:), prlat(:,:), prlon(:,:), usdum(:)&
         ,api_temp(:,:,:), tempdum(:),smcl(:,:,:,:),t_soil(:,:,:,:),t_soilJ(:,:,:)
      real               :: dif_time, dummy
      integer            :: int_dif_time, idate2, imonth2, iyear2, hourmin, &
         da, sair
      integer            :: ierr
      !!$    integer, parameter :: idim_type_min = 2
      !!$    integer, parameter :: idim_type_max = 7
      integer, parameter :: idim_type     = 4
      !!$    integer            :: il1(nmachs)
      !!$    integer            :: ir2(nmachs)
      !!$    integer            :: jb1(nmachs)
      !!$    integer            :: jt2(nmachs)
      !!$    integer            :: localsize(nmachs,idim_type_min:idim_type_max)
      !!$    integer            :: disp(nmachs,idim_type_min:idim_type_max)
      !!$    integer            :: maxlocalsize
      !!$    integer            :: sizegathered(idim_type_min:idim_type_max)
      !!$    integer            :: maxsizegathered
      !!$    integer            :: sizefullfield(idim_type_min:idim_type_max)
      !!$    integer            :: maxsizefullfield
      !!$    integer            :: globalsize(idim_type_min:idim_type_max)
      !!$    real, allocatable  :: localchunk(:)
      !!$    real, allocatable  :: gathered(:)
      !!$    real, allocatable  :: fullfield(:)
      real :: usdum_before(200,200,50)
      character(len=16)  :: varn
      integer            :: i0, j0, ia, iz, ja, jz,irec,nrec
      logical            :: general
      character(len=256) :: dset,title,tdef,lixo
      character(len=32),allocatable :: soilVarName(:)
      integer, allocatable :: slevels(:)
      real,allocatable :: us_geos(:,:,:)
      real,allocatable :: soilM_GFS(:,:,:),soilT_GFS(:,:,:)
      real,allocatable :: lats_(:),lons_(:)
      real :: undef, shift_lon
      integer :: vars,ki
      integer :: USE_SMOOTH_PROF     = 0
      real,    dimension (1:20,2)     ::  tend2d
      real,    dimension (2)          ::  tend1d
      integer, external :: outRealSize
      character(LEN(usdata)+100)  :: outs  !DSM

      namelist /gradeumso/ latni, latnf, lonni, lonnf, ilatn, ilonn, nlat, nlon

      lpw=int(lpw_r)
      ! initial grid definitions for grid ifm
      i0 = nodei0(mynum,ifm)
      j0 = nodej0(mynum,ifm)
      ia = nodei0(mynum,ifm) + 1
      iz = nodei0(mynum,ifm) + nodemxp(mynum,ifm)
      ja = nodej0(mynum,ifm) + 1
      jz = nodej0(mynum,ifm) + nodemyp(mynum,ifm)

      iyear2  = iyear1
      imonth2 = imonth1
      idate2  = idate1

      ! determinacao do tipo de produto de umidade
      !DSM     do i=256,1,-1
      !DSM        if (usdata_in(i:i)=='/') then
      !DSM           ipref_start = i + 1
      !DSM           exit
      !DSM        endif
      !DSM     enddo
      ipref_start=index(usdata_in,'/',BACK = .true.)+1

      shift_lon = 0.

      ! definicao da espessura das camadas
      ipref = len_trim(usdata_in)
      pref = usdata_in(ipref_start:ipref)

      if (trim(pref)=='SM_V2.') then
         n4us = 6                    ! modelo v2 com 6 camadas
         allocate(slz_us(0:n4us))
         slz_us = (/-3.0, -2.0, -1.0, -0.5, -0.25, -0.1, 0. /)
      elseif (trim(pref)=='GL_SM.GPCP.' .or.  trim(pref)=='GL_SM.GPNR.' ) then
         n4us = 8                    ! modelo glsm v2 com 8 camadas
         allocate(slz_us(0:n4us))
         slz_us = (/-4.5, -2.5, -1.75, -1.0, -0.5, -0.25, -0.13, -0.05, 0./)
      elseif (trim(pref)=='GFS025GR.SOILM.') then
         n4us = 4                    ! modelo do gfs 4 camadas
         allocate(slz_us(0:n4us))
         slz_us = (/-2.0, -1.0,  -0.40, -0.10, 0./)
      elseif (trim(pref)=='SM.GEOS.') then
         !do nothing

      elseif (usdata_in(ipref-10:ipref)=='YYYYMMDD.nc') then  !DSM - Lendo a umidade do solo proveniente do JULES
#ifdef cdf
         write (cidate, '(i2.2)') idate2
         write (cimon,  '(i2.2)') imonth2
         write (ciyear, '(i4)')   iyear2
         write (cgrid,  '(i1)')   ifm

         !usdata = usdata_in(1:ipref-11)//ciyear//cimon//cidate//'.nc'
         usdata = trim(usdata_in)
         call replace_text(usdata,'YYYY',ciyear,outs);usdata=trim(outs)
         call replace_text(usdata,'MM',cimon,outs);usdata=trim(outs)
         call replace_text(usdata,'DD',cidate,outs);usdata=trim(outs)

         inquire(file=usdata,exist=theref)
         if (.not. theref) then
            print*;print*;print*, 'ERROR! Not found the file:'//trim(usdata);print*
            stop 'In:  src/brams/soil_moisture/soilMoisture.F90'
         endif

         ifname = len_trim(usdata)
         STATUSi = nf90_open(usdata(1:ifname), NF90_NOWRITE, ncid); if (STATUSi/=0) stop 'ERROR in: nf90_open'
         STATUSi = nf_INQ_DIMLEN(ncid, 4, n4us); if (STATUSi/=0) stop 'ERROR in: nf_INQ_DIMLEN - n4us'

         if (n4us /= 7) then
            print*, 'nsoil of '//usdata(1:ifname)//' is /= 7. (nsoil=',n4us,')'
            stop
         endif
         allocate(slz_us(0:n4us))
         slz_us = (/-7.0, -5.0, -3.4, -2.0, -1.0, -0.4, -0.1, 0./)
#else

         print*
         print*, '=> To use JULES soil moisture ('//usdata_in(ipref-10:ipref)//') is necessary compiler the model with NetCDF'
         stop 'STOP in: src/brams/soil_moisture/soilMoisture.F90'

#endif
      elseif (trim(pref)=='GFS.SOIL:UMID_TEMP.') then
         n4us = 4                    ! modelo GFS com 4 camadas
         allocate(slz_us(0:n4us))
         slz_us = (/ -2.0, -1.0, -0.4, -0.1, 0. /)
         shift_lon = 360.
         !print*,"UMID1=",trim(pref),slz_us;call flush(6)

      elseif (trim(pref)=='ERA5.SOIL:UMID_TEMP.') then
         n4us = 4                    ! modelo ERA5 com 4 camadas
         allocate(slz_us(0:n4us))
         slz_us = (/ -2.89, -1.0, -0.28, -0.07, 0. /)
         shift_lon = 0.
      else
         n4us = 4                    ! modelo original com 4 camadas
         allocate(slz_us(n4us))
         slz_us = (/-2.4, -0.4, -0.1, 0./)

      endif

      ! composicao do nome dos arquivos de entrada e saida

      if ( (runtype(1:7)=='history') .and.                 &
         ((soil_moist=='h') .or. (soil_moist=='H') .or.  &
         (soil_moist=='a') .or. (soil_moist=='A'))) then

         dif_time = timstr

         !if (timeunit == 'h') then
         !int_dif_time = floor(dif_time/24.)
         !elseif (timeunit == 'm') then
         !int_dif_time = floor(dif_time/1440.)
         !elseif (timeunit == 's') then
         int_dif_time = floor(dif_time/86400.)
         !endif

         call changeday(idate1, imonth1, iyear1, int_dif_time,   &
            idate2, imonth2, iyear2)

      else

         int_dif_time = 0

      endif

      if ((soil_moist_fail=='l') .or. (soil_moist_fail=='L')) then
         ! looking for until 5 days old file
         da = 5
      else
         da = 1
      endif

      sair = 0

      do i=1,da

         write (cidate, '(i2.2)') idate2
         write (cimon,  '(i2.2)') imonth2
         write (ciyear, '(i4)')   iyear2
         write (cgrid,  '(i1)')   ifm

         ! calculating the hour of simulation
         if ((itime1>=0000) .and. (itime1<1200)) then
            hourmin = 0000
            if (trim(pref)=='GL_SM.GPCP.' .or. trim(pref)=='GL_SM.GPNR.' .or. trim(pref)=='GFS.SOIL:UMID_TEMP.' &
               .or. trim(pref)=='ERA5.SOIL:UMID_TEMP.' ) hourmin = 00
         else
            hourmin = 1200
            if (trim(pref)=='GL_SM.GPCP.'         .or. trim(pref)=='GL_SM.GPNR.'          ) hourmin = 12
            if (trim(pref)=='GFS.SOIL:UMID_TEMP.' .or. trim(pref)=='ERA5.SOIL:UMID_TEMP.' ) hourmin = 21
         endif

         if (trim(pref)=='GL_SM.GPCP.' .or. trim(pref)=='GL_SM.GPNR.' .or.  trim(pref)=='GFS.SOIL:UMID_TEMP.' &
            .or. trim(pref)=='ERA5.SOIL:UMID_TEMP.' ) then
            write (cihourmin, '(i2.2)') hourmin
            icihourmin = 2
         else
            write (cihourmin, '(i4.4)') hourmin
            icihourmin = 4
         endif

         if (trim(pref)=='us') then
            cihourmin  = ''
            icihourmin = 0
         endif

         if (trim(pref)=='GFS.SOIL:UMID_TEMP.' .or. trim(pref)=='ERA5.SOIL:UMID_TEMP.' ) then
            usdata = trim(usdata_in)//ciyear//cimon//cidate//cihourmin(1:icihourmin)//'.bin'
         elseif (usdata_in(ipref-10:ipref)/='YYYYMMDD.nc') then
            usdata = trim(usdata_in)//ciyear//cimon//cidate//cihourmin(1:icihourmin)//'.vfm'
         endif
         ifname = len_trim(usdata)

         if (mchnum==master_num) inquire(file=usdata(1:ifname),exist=theref)
         call parf_bcast(theref, master_num)

         if (.not.theref) &
            usdata = trim(usdata_in)//ciyear//cimon//cidate//cihourmin(1:icihourmin)//'.gra'

         usmodel = trim(usmodel_in)//ciyear//cimon//cidate//cihourmin(1:icihourmin)//'_g'//cgrid//'.mod'

         if (mchnum==master_num) then
            print *, 'looking up soil moisture date files: '
            print *, '  usdata : ', usdata(1:len_trim(usdata))
            print *, '  usmodel: ', usmodel(1:len_trim(usmodel))

            ifname = len_trim(usmodel)
            inquire(file=usmodel(1:ifname), exist=there)

            if(.not.there) then

               ifname = len_trim(usdata)
               inquire(file=usdata(1:ifname), exist=there)

               if (there) then
                  sair = 1
               else
                  print *, 'files:', usdata(1:ifname), &
                     ' and: ', usmodel(1:ifname), ' not found!'
               endif
            else
               sair = 1
            endif
         endif !(mchnum==master_num)

         call parf_bcast(sair, master_num)

         if (sair==1) exit

         call changeday(idate1, imonth1, iyear1, (int_dif_time - i),   &
            idate2, imonth2, iyear2)

      enddo

      ifname = len_trim(usmodel)

      if (mchnum==master_num) inquire(file=usmodel(1:ifname), exist=there)
      call parf_bcast(there, master_num)

      !    because the the soil moisture dataset interpolated to the model grid
      !    is never checked again, we will allways require the interpolation to
      !    be done every model initialization. this will prevents soil moisture
      !    interpolated for given grid specification being erroneously be used
      !    for a different model configuration

      ifname = len_trim(usdata)
      if (mchnum==master_num) inquire(file=usdata(1:ifname), exist=there)
      call parf_bcast(there, master_num)

      if (.not.there) then
         if (mchnum==master_num) then
            print *, '  usdata : ', usdata(1:ifname)
            print *, '  usmodel: ', usmodel(1:len_trim(usmodel))
         end if
         if ((soil_moist_fail=='s') .or. (soil_moist_fail=='')) then
            call fatal_error('* heterogenous soil moisture init. error! '//&
               '  program will stop!')
            !!stop
         else
            if (mchnum==master_num) then
               print *, '  homogeneous soil moisture initialization.'
            end if
            return
         endif
      endif

      c1 = 0.5*cpi

      do j=1,n3
         do i=1,n2

            k2  = lpw(i,j)
            pis = c1*(pi0(k2-1,i,j) + pi0(k2,i,j) + pp(k2-1,i,j) + pp(k2,i,j))
            airtemp = theta(k2,i,j)*pis

            do ipat=2,npat
               do k=1,mzg
                  nsoil = nint(soil_text(k,i,j,ipat))
                  soil_water(k,i,j,ipat) = &
                     max(soilcp(nsoil), slmstr(k)*slmsts(nsoil))
                  soil_energy(k,i,j,ipat) = &
                     (airtemp - 273.15 + stgoff(k))*  &
                     (slcpd(nsoil) + soil_water(k,i,j,ipat)*4.186e6) + &
                     soil_water(k,i,j,ipat)*3.34e8
               enddo
            enddo
         enddo
      enddo

      !--------------------------------------------------------- heterogeneous
      if (mchnum==master_num) then
         print*,'----------------------------------------------'
         print*,'       soil moisture initialization'
         print*,'----------------------------------------------'
      endif

      ! dados da grade de precipitacao

      if (mchnum==master_num) &
         inquire(file=trim(usdata_in)//'_ent', exist=general)

      call parf_bcast(general, master_num)

      if (general) then
         if (mchnum==master_num) then
            open (unit=93,file=usdata_in(1:len_trim(usdata_in))//'_ent',status='old')
            read (unit=93,nml=gradeumso) !latni, latnf, lonni, lonnf, ilatn, ilonn, nlat, nlon
            close (unit=93, status='keep')
         endif
         call parf_bcast(latni, master_num)
         call parf_bcast(latnf, master_num)
         call parf_bcast(lonni, master_num)
         call parf_bcast(lonnf, master_num)
         call parf_bcast(ilatn, master_num)
         call parf_bcast(ilonn, master_num)
         call parf_bcast(nlat,  master_num)
         call parf_bcast(nlon,  master_num)

      elseif (trim(pref)=='GL_SM.GPCP.') then
         latni =  -89.5            !  gpcp global
         latnf =   89.5
         lonni = -179.5
         lonnf =  179.5
         ilatn =    1.
         ilonn =    1.
         nlat  =  180
         nlon  =  360
      elseif (trim(pref)=='GL_SM.GPNR.') then
         latni =  -89.875             !  trmm/navy + gpcp global
         latnf =   89.875
         lonni = -179.875
         lonnf =  179.875
         ilatn =    0.250
         ilonn =    0.250
         nlat  =  720
         nlon  = 1440


         !====== GFS soil moisture (volumetric: [0-1]) and soil temperature (Kelvin)
      elseif (trim(pref)=='GFS.SOIL:UMID_TEMP.') then
         latni =  -89.875
         latnf =   89.875
         lonni =   0. !-179.875
         lonnf =   360. !179.875
         ilatn =   0.250
         ilonn =   0.250
         nlat  =  721
         nlon  = 1440
         allocate( soilM_GFS(nlon,nlat,n4us),  soilT_GFS(nlon,nlat,n4us))

         if(mchnum==master_num) then

            write(*,*) '===================================='
            write(*,*) 'Input Soil Moisture File Inventory: '
            write(*,*) '===================================='
            write(*,*) 'Filename    : ',usdata(1:len_trim(usdata)-3)//'bin'


            open(2, status='old', form='unformatted', access='direct', &
               recl=    outRealSize()*nlat*nlon, file=usdata(1:len_trim(usdata)))
            !--- check this later         recl=4 * nlat*nlon, file=usdata(1:len_trim(usdata)))
            nrec=1
            do k=1,n4us
               read(2,rec=nrec) soilM_GFS(1:nlon,1:nlat,k)
               print*,"reading SoilMoist GFS :",k,maxval(soilM_GFS(1:nlon,1:nlat,k)),minval(abs(soilM_GFS(1:nlon,1:nlat,k)))
               nrec=nrec+1
            enddo
            do k=1,n4us
               read(2,rec=nrec) soilT_GFS(1:nlon,1:nlat,k)
               print*,"reading SoilTemp GFS :",k,"of",n4us,maxval(soilT_GFS(1:nlon,1:nlat,k)),minval(abs(soilT_GFS(1:nlon,1:nlat,k)))
               nrec=nrec+1
            enddo
            close(2)
            print*,"end of reading SoilTemp GFS :",n4us
         endif

         call parf_bcast(soilM_GFS, int(nlon,i8), int(nlat,i8), int(n4us,i8), master_num)
         call parf_bcast(soilT_GFS, int(nlon,i8), int(nlat,i8), int(n4us,i8), master_num)


         !====== ERA5 soil moisture (wetness: [0-1]) and soil temperature (Kelvin)
      elseif (trim(pref)=='ERA5.SOIL:UMID_TEMP.') then
         latni =  -59.875
         latnf =   29.875
         lonni =  -110. !-179.875
         lonnf =   0. !179.875
         ilatn =   0.250
         ilonn =   0.250
         nlat  =  361
         nlon  =  441
         allocate( soilM_GFS(nlon,nlat,n4us),  soilT_GFS(nlon,nlat,n4us))

         if(mchnum==master_num) then

            write(*,*) '===================================='
            write(*,*) 'Input Soil Moisture File Inventory: '
            write(*,*) '===================================='
            write(*,*) 'Filename    : ',usdata(1:len_trim(usdata)-3)//'bin'


            open(2, status='old', form='unformatted', access='direct', &
               recl=    outRealSize()*nlat*nlon, file=usdata(1:len_trim(usdata)))
            !--- check this later         recl=4 * nlat*nlon, file=usdata(1:len_trim(usdata)))
            nrec=1
            do k=1,n4us
               read(2,rec=nrec) soilM_GFS(1:nlon,1:nlat,k)
               print*,"reading SoilMoist ERA5:",k,maxval(soilM_GFS(1:nlon,1:nlat,k)),minval(abs(soilM_GFS(1:nlon,1:nlat,k)))
               nrec=nrec+1
            enddo
            do k=1,n4us
               read(2,rec=nrec) soilT_GFS(1:nlon,1:nlat,k)
               print*,"reading SoilTemp  ERA5:",k,maxval(soilT_GFS(1:nlon,1:nlat,k)),minval(abs(soilT_GFS(1:nlon,1:nlat,k)))
               nrec=nrec+1
            enddo
            close(2)
         endif

         call parf_bcast(soilM_GFS, int(nlon,i8), int(nlat,i8), int(n4us,i8), master_num)
         call parf_bcast(soilT_GFS, int(nlon,i8), int(nlat,i8), int(n4us,i8), master_num)






         !====== GEOS soil moisture
      elseif(trim(pref)=='SM.GEOS.') then

         !Opening the descriptor file to read info
         print *, usdata(1:len_trim(usdata)-3)//'ctl'
         open(2, status='old', form='formatted', action='read', &
            file=usdata(1:len_trim(usdata)-3)//'ctl')

         !First time: read variables and allocate all
         read(2,*) lixo,dset
         read(2,*) lixo,undef
         read(2,*) lixo,title
         read(2,*) lixo,nlon,lixo,lonni,ilonn
         read(2,*) lixo,nlat,lixo,latni,ilatn
         read(2,*) lixo,n4us
         n4us=n4us-1
         allocate(slz_us(0:n4us))
         read(2,*) lixo,lixo,lixo,tdef
         read(2,*) lixo,vars
         allocate(SoilVarName(vars),Slevels(vars))


         !Goto to begin of file and read all variables
         rewind(unit=2)
         read(2,*) lixo,dset
         read(2,*) lixo,undef
         read(2,*) lixo,title
         read(2,*) lixo,nlon,lixo,lonni,ilonn
         read(2,*) lixo,nlat,lixo,latni,ilatn
         read(2,*) lixo,lixo,lixo,slz_us(0:n4us)
         read(2,*) lixo,lixo,lixo,tdef
         read(2,*) lixo,vars
         do k=1,vars
            read(2,*)  SoilVarName(k),Slevels(k)
         enddo

         close(unit=2)

         !Determine the end of lats and lons
         latnf=latni+(nlat-1)*ilatn
         lonnf=lonni+(nlon-1)*ilonn

         if(mchnum==master_num) then
            write(*,*) c_lightYellow
            write(*,*) '===================================='
            write(*,*) 'Input Soil Moisture File Inventory: '
            write(*,*) '===================================='
            write(*,*) 'Filename    : ',usdata(1:len_trim(usdata)-3)//'ctl'
            write(*,*) 'Bin filename: ',trim(dset)
            write(*,*) 'Undef Value : ',undef
            write(*,*) 'File Title  : ' ,trim(title)
            write(*,*) 'Lons        : ',nlon,lonni,lonnf,ilonn
            write(*,*) 'Lats        : ',nlat,latni,latnf,ilatn
            write(*,fmt='(A,I2.2,A)') ' # of levels : (00:',n4us,')'
            write(*,*) 'Levels      : ',slz_us(0:n4us)
            write(*,*) 'Time of data: ',trim(tdef)
            write(*,*) '# of vars   : ',vars
            do k=1,vars
               write(*,fmt='(A,I1,A,A,1X,I2.2)') ' Var(',k,')   : ',SoilVarName(k),Slevels(k)
            enddo
            write(*,*) '===================================='
            write(*,*) c_noColor

            !allocate(lats_(nlat),lons_(nlon))
            !do j=1,nlat
            !       lats_(i)=latni+(i-1)*ilatn
            !enddo
            !do i=1,nlon
            !       lons_(j)=lonni+(j-1)*ilonn
            !enddo

         endif

         allocate(us_geos(nlon,nlat,0:n4us))

         if(mchnum==master_num) then
            !Reading the file
            open(2, status='old', form='unformatted', access='direct', &
               recl=4*nlat*nlon, file=usdata(1:len_trim(usdata)))
            nrec=1
            do k=0,n4us,1
               !Geos SM data is upside down
               ki=n4us-k
               read(2,rec=nrec) us_geos(1:nlon,1:nlat,ki)
               nrec=nrec+1
            enddo
            close(2)

            !open(unit=2,file='./us.gra',action='WRITE',status='REPLACE' &
            !   ,form='UNFORMATTED',access='DIRECT', &
            !recl=4*nlat*nlon)
            !nrec=1
            !do k=0,n4us,1
            !   write (2,rec=nrec) us_geos(1:nlon,1:nlat,k)
            !   nrec=nrec+1
            !enddo
            !close(2)
            !
            !open(unit=2,file='./us.ctl',action='WRITE',status='replace' &
            !   ,form='FORMATTED')
            !write(2,*) 'dset ^./us.gra'
            !write(2,*) 'undef',undef
            !write(2,*) 'title us'
            !write(2,*) 'xdef',nlon,'linear',lonni,ilonn
            !write(2,*) 'ydef',nlat,'linear',latni,ilatn
            !write(2,*) 'zdef',n4us+1,'levels',slz_us(0:n4us)
            !write(2,*) 'tdef   1 linear ',tdef,'  1hr'
            !write(2,*) 'vars 1'
            !write(2,*) 'us',n4us+1,'99 Us'
            !write(2,*) 'endvars'
            !close(2)
         endif

         call parf_bcast(us_geos, int(nlon,i8), int(nlat,i8), int(n4us+1,i8), master_num)

      elseif (usdata_in(ipref-10:ipref)=='YYYYMMDD.nc') then  !DSM - Lendo a umidade do solo proveniente do JULES
#ifdef cdf
         ! --- Lendo as dimensoes do arquivo .nc ---
         STATUSi = nf_INQ_DIMLEN(ncid, 1, nlon); if (STATUSi/=0) stop 'ERROR in: nf_INQ_DIMLEN - nlon'
         STATUSi = nf_INQ_DIMLEN(ncid, 2, nlat); if (STATUSi/=0) stop 'ERROR in: nf_INQ_DIMLEN - nlat'
         allocate (prlat(nlon,nlat),prlon(nlon,nlat))
         STATUSi = nf_inq_varid(ncid, 'latitude', var_id) ; if (STATUSi/=0) stop 'ERROR in: nf_inq_varid - latitude'
         STATUSi = nf_get_var(ncid, var_id, prlat) ; if (STATUSi/=0) stop 'ERROR in: nf_get_var - prlat'
         STATUSi = nf_inq_varid(ncid, 'longitude', var_id) ; if (STATUSi/=0) stop 'ERROR in: nf_inq_varid - longitude'
         STATUSi = nf_get_var(ncid, var_id, prlon) ; if (STATUSi/=0) stop 'ERROR in: nf_get_var - prlon'

         undef=-1e+20

         lonni=prlon(1,1)
         lonnf=prlon(nlon,1)
         ilonn=prlon(2,1)-prlon(1,1)
         latni=prlat(1,1)
         latnf=prlat(1,nlat)
         ilatn=prlat(1,2)-prlat(1,1)
#endif
      else
         call fatal_error('unexpected prefix for soil moisture')
      endif

      if (usdata_in(ipref-10:ipref)/='YYYYMMDD.nc') then  !DSM - Lendo a umidade do solo proveniente do JULES
         allocate(prlat(nlon,nlat), stat=ierr)
         if (ierr/=0) call fatal_error("error allocating prlats (soilmoistureinit)")
         allocate(prlon(nlon,nlat), stat=ierr)
         if (ierr/=0) call fatal_error("error allocating prlon (soilmoistureinit)")

         call apiprlatlon(nlon, nlat, prlat, prlon, ilatn, ilonn, latni, lonni)
      endif

      allocate(api_us(nlon,nlat,n4us), stat=ierr)
      if (ierr/=0) call fatal_error("error allocating api_us (soilmoistureinit)")


      allocate(api_temp(nlon,nlat,n4us), stat=ierr)
      if (ierr/=0) call fatal_error("error allocating api_temp (soilmoistureinit)")

      allocate(usdum(n4us), stat=ierr)
      if (ierr/=0) call fatal_error("error allocating usdum (soilmoistureinit)")

      allocate(tempdum(n4us), stat=ierr)
      if (ierr/=0) call fatal_error("error allocating tempdum (soilmoistureinit)")


      if (trim(pref) .ne. 'GFS.SOIL:UMID_TEMP.' .and. trim(pref) .ne. 'ERA5.SOIL:UMID_TEMP.') then


         if (mchnum==master_num) then
            print *, '-------------------------------------- grid=', ifm
            print *, 'opening soil moisture data= ', trim(usdata)
         endif

         if (.not.theref) then ! arquivo .gra
            if ( trim(pref)=='SM.GEOS.') then   ! arquivo .gra acesso direto
               if (mchnum==master_num) then
                  do k=1,n4us,1
                     api_us(1:nlon,1:nlat,k)=us_geos(1:nlon,1:nlat,k) ! wetness
                  enddo
               endif

               call parf_bcast(api_us, int(nlon,i8), int(nlat,i8), int(n4us,i8), &
                  master_num)
            endif

         elseif (usdata_in(ipref-10:ipref)=='YYYYMMDD.nc') then  !DSM - Lendo a temperatura e umidade do solo proveniente do JULES
#ifdef cdf

            allocate (smcl(nlon,nlat,n4us,4),t_soil(nlon,nlat,n4us,4),t_soilJ(nlon,nlat,n4us))
            STATUSi = nf_inq_varid(ncid, 'smcl', var_id) ; if (STATUSi/=0) stop 'ERROR in: nf_inq_varid - smcl'
            STATUSi = nf_get_var(ncid, var_id, smcl) ; if (STATUSi/=0) stop 'ERROR in: nf_get_var - smcl'
            STATUSi = nf_inq_varid(ncid, 't_soil', var_id) ; if (STATUSi/=0) stop 'ERROR in: nf_inq_varid - smcl'
            STATUSi = nf_get_var(ncid, var_id, t_soil) ; if (STATUSi/=0) stop 'ERROR in: nf_get_var - t_soil'

            !( x,  y, z,t)
            !print*, 4,n4us,nlon,nlat, smcl(312,403,4,2)
            !slz_us = (/-7.0, -5.0, -3.4, -2.0, -1.0, -0.4, -0.1, 0./)

            select case (cihourmin(1:icihourmin-2))
               case ('00')
                  it=1
               case ('06')
                  it=2
               case ('12')
                  it=3
               case ('18')
                  it=4
               case default
                  stop 'ERROR! Not found this hour in JULES soil moisture file'
            end select

            do k=1,n4us
               api_us(:,:,k)=smcl(:,:,n4us-k+1,it)/abs(slz_us(k-1)-slz_us(k))/1000.0
               t_soilJ(:,:,k)=t_soil(:,:,n4us-k+1,it)   !utilizando a variavel soil_energy para armazenar t_soil
            enddo
            !print*,api_us(312,403,:)
#endif
         else  ! arquivo .vfm
            if (mchnum==master_num) then
               open(unit=2, file=usdata(1:len_trim(usdata)), form='formatted', status='old')
               call vfirec(2, api_us, (nlat*nlon*n4us), 'lin')
               close(2)
            endif

            call parf_bcast(api_us, int(nlon,i8), int(nlat,i8), int(n4us,i8), master_num)

         endif

         ! gathering data
         varn = 'soil_water'
         call gatherdata(idim_type, varn, ifm, mzg, nnxp(ifm), nnyp(ifm), &
            npat, nmachs, mchnum, mynum, master_num,                    &
            soil_water, globalsoilwater)
         varn = 'soil_text'
         call gatherdata(idim_type, varn, ifm, mzg, nnxp(ifm), nnyp(ifm), &
            npat, nmachs, mchnum, mynum, master_num,                    &
            soil_text, globalsoiltext)
         varn = 'glon'
         call gatherdata(2, varn, ifm, nnxp(ifm), nnyp(ifm), &
            nmachs, mchnum, mynum, master_num,             &
            glon, globalglon)
         varn = 'glat'
         call gatherdata(2, varn, ifm, nnxp(ifm), nnyp(ifm), &
            nmachs, mchnum, mynum, master_num,             &
            glat, globalglat)


         ! loop no dominio global do modelo
         globalsoilenergy=0.0
         do j=1,nnyp(ifm)
            do i=1,nnxp(ifm)
               ! evitando pontos fora do dominio da grade de umidade de solo
               if ( globalglat(i,j)<latni .or. globalglat(i,j)>latnf .or. &
                  globalglon(i,j)<lonni .or. globalglon(i,j)>lonnf) cycle

               call interpolacao(globalglon(i,j), globalglat(i,j), nlon, nlat, &
                  prlat, prlon, i1, i2, ic, j1, j2, jc)
               if (ic>=0 .and. jc>=0) then
                  dlonr = 0.5*(globalglon(nnxp(ifm)-1,j) - globalglon(1,j))/&
                     float(nnxp(ifm)-1)
                  dlatr = 0.5*(globalglat(i,nnyp(ifm)-1) - globalglat(i,1))/&
                     float(nnyp(ifm)-1)
                  qi1   = int(dlonr/ilonn+0.5)
                  qi2   = int(dlonr/ilonn+0.5)
                  qj1   = int(dlatr/ilatn+0.5)
                  qj2   = int(dlatr/ilatn+0.5)

                  do k=1,n4us
                     ncount = 0
                     usdum(k)=0.
                     i_ini=max(1,ic-qi1)
                     i_fim=min(nlon,ic+qi2)
                     j_ini=max(1,jc-qj1)
                     j_fim=min(nlat,jc+qj2)
                     do jj=j_ini,j_fim
                        do ii=i_ini,i_fim
                           if (usdata_in(ipref-10:ipref)=='YYYYMMDD.nc') then  !DSM {
                              !--- Encontrando o valor mais proximo que nao seja indefinido - como a grade dos dados eh
                              !--- diferente ocorre de pegar pontos que estao sobre a agua (nao calculado pelo jules offline).
                              iii=ii
                              jjj=jj
                              !print*,'oi0',1,1,'|',iii,jjj,t_soilJ(iii,jjj,k)
                              if (t_soilJ(iii,jjj,k) < 0.) then
                                 do range=1,2 ! max(iii,nlon-iii+1) --- Considerando que eh oceano caso for acima de 4 pontos
                                    do js=-1*range,range
                                       if (jj+js>=1 .and. jj+js<=nlat) jjj=jj+js
                                       if (abs(js)==range) then
                                          skip=1
                                       else
                                          skip=2*range
                                       endif
                                       do is=-1*range,range,skip
                                          if (ii+is>=1 .and. ii+is<=nlon) iii=ii+is
                                          !print*,'oi1',is,js,'|',iii,jjj,t_soilJ(iii,jjj,k)
                                          if (t_soilJ(iii,jjj,k) > 0.) exit
                                       enddo
                                       if (t_soilJ(iii,jjj,k) > 0.) exit
                                    enddo
                                    if (t_soilJ(iii,jjj,k) > 0.) exit
                                 enddo
                              endif

                              do ipat=2,npat
                                 globalsoilenergy(k,i,j,ipat)=t_soilJ(iii,jjj,k)
                              enddo
                              api_us(ii,jj,k)=api_us(iii,jjj,k)
                           endif  !DSM }


                           if (api_us(ii,jj,k)>1.e-5) then
                              !If the value is undef probably is ocean
                              if(trim(pref)=='SM.GEOS.') then
                                 if(api_us(ii,jj,k)>=undef) then
                                    api_us(ii,jj,k)=0.0
                                 endif
                              endif
                              do ipat=2,npat
                                 ncount = ncount + 1
                                 usdum(k) = usdum(k) + api_us(ii,jj,k)
                              enddo
                           endif

                        enddo
                     enddo

                     usdum(k) = usdum(k)/(float(ncount) + 1.e-10)
                  enddo

                  !usdum_before(i,j,k)=usdum(k)

                  do k=mzg,1,-1
                     do kk=n4us,1,-1
                        if (slz(k)>=slz_us(kk)) then
                           do ipat=2,npat
                              nsoil = nint(globalsoiltext(k,i,j,ipat))
                              !- umidade em mm3/mm3 (gfs case)
                              globalsoilwater(k,i,j,ipat) = usdum(kk+1)
                              !write(90,fmt='(6(I4.4,1X),3(F18.6,1X))') k,kk,i,j,ipat,nsoil,globalsoilwater(k,i,j,ipat),usdum(kk+1),slmsts(nsoil)
                              ! umidade lida em % (armazenamento)
                              if ( trim(pref)=='sm_v2.'      .or. &
                                 trim(pref)=='gl_sm.gpcp.' .or. &
                                 trim(pref)=='gl_sm.gpnr.' .or. &
                                 trim(pref)=='SM.GEOS.') then
                                 globalsoilwater(k,i,j,ipat) = &
                                    globalsoilwater(k,i,j,ipat)*slmsts(nsoil)
                              endif
                              if (usdum(kk+1)<1.e-5) &
                                 globalsoilwater(k,i,j,ipat) = slmstr(k)*slmsts(nsoil) !oceano
                              globalsoilwater(k,i,j,ipat) = max(soilcp(nsoil), &
                                 min(globalsoilwater(k,i,j,ipat), &
                                 slmsts(nsoil)))
                           !write(91,fmt='(6(I4.4,1X),3(F18.6,1X))') k,kk,i,j,ipat,nsoil,globalsoilwater(k,i,j,ipat),usdum(kk+1),slmsts(nsoil)
                           enddo
                           exit
                        elseif (slz(k)<slz_us(1)) then
                           do ipat=2,npat
                              nsoil = nint(globalsoiltext(k,i,j,ipat))
                              !- umidade em mm3/mm3 (gfs case)
                              globalsoilwater(k,i,j,ipat) = usdum(1)
                              ! umidade lida em % (armazenamento) - v2 (sm_v2.)
                              if ( trim(pref)=='sm_v2.'      .or. &
                                 trim(pref)=='gl_sm.gpcp.' .or. &
                                 trim(pref)=='gl_sm.gpnr.' .or. &
                                 trim(pref)=='SM.GEOS.') then
                                 globalsoilwater(k,i,j,ipat) = &
                                    globalsoilwater(k,i,j,ipat)*slmsts(nsoil)
                              endif
                              if (usdum(1)<1.e-5)   &
                                 globalsoilwater(k,i,j,ipat) = slmstr(k)*slmsts(nsoil)!oceano
                              globalsoilwater(k,i,j,ipat) = max(soilcp(nsoil), &
                                 min(globalsoilwater(k,i,j,ipat), &
                                 slmsts(nsoil)))
                           !write(92,fmt='(6(I4.4,1X),3(F18.6,1X))') k,kk,i,j,ipat,nsoil,globalsoilwater(k,i,j,ipat),usdum(kk+1),slmsts(nsoil)

                           enddo
                           exit
                        endif
                     enddo
                  enddo
               endif
            enddo
         enddo

         !- GFS/ERA5 soil moisture and soil temperature
      elseif (trim(pref) == 'GFS.SOIL:UMID_TEMP.' .or. trim(pref) == 'ERA5.SOIL:UMID_TEMP.') then
         do k=1,n4us
            api_us  (:,:,k)= soilM_GFS(:,:,n4us-k+1)
            api_temp(:,:,k)= soilT_GFS(:,:,n4us-k+1)
         enddo
         where(api_us   < 0.) api_us  = 0.2  !ocean
         where(api_temp < 0.) api_temp= 300. !ocean

         ! gathering data
         varn = 'soil_water'
         call gatherdata(idim_type, varn, ifm, mzg, nnxp(ifm), nnyp(ifm), &
            npat, nmachs, mchnum, mynum, master_num,             &
            soil_water, globalsoilwater)
         varn = 'soil_energy'
         call gatherdata(idim_type, varn, ifm, mzg, nnxp(ifm), nnyp(ifm), &
            npat, nmachs, mchnum, mynum, master_num,             &
            soil_energy, globalsoilenergy)
         varn = 'soil_text'
         call gatherdata(idim_type, varn, ifm, mzg, nnxp(ifm), nnyp(ifm), &
            npat, nmachs, mchnum, mynum, master_num,             &
            soil_text, globalsoiltext)
         varn = 'glon'
         call gatherdata(2, varn, ifm, nnxp(ifm), nnyp(ifm), &
            nmachs, mchnum, mynum, master_num,      &
            glon, globalglon)
         varn = 'glat'
         call gatherdata(2, varn, ifm, nnxp(ifm), nnyp(ifm), &
            nmachs, mchnum, mynum, master_num,      &
            glat, globalglat)

         ! loop no dominio global do modelo
         do j=1,nnyp(ifm)
            do i=1,nnxp(ifm)

               call interpolacao(globalglon(i,j)+shift_lon, globalglat(i,j), nlon, nlat, prlat, prlon, i1, i2, ic, j1, j2, jc)

               !print*,"lonlat=",globalglon(i,j), globalglat(i,j),ic,jc , prlat(ic,jc), prlon(ic,jc),ic,jc

               if (ic>=0 .and. jc>=0) then

                  dlonr = 0.5*(globalglon(nnxp(ifm)-1,j) - globalglon(1,j))/float(nnxp(ifm)-1)
                  dlatr = 0.5*(globalglat(i,nnyp(ifm)-1) - globalglat(i,1))/float(nnyp(ifm)-1)

                  qi1   = int(dlonr/ilonn+0.5); qi2   = int(dlonr/ilonn+0.5)
                  qj1   = int(dlatr/ilatn+0.5); qj2   = int(dlatr/ilatn+0.5)

                  do k=1,n4us
                     !---reset to zero
                     ncount = 0; usdum(k)=0.; tempdum(k)=0.

                     do jj=max(1,jc-qj1),min(nlat,jc+qj2)
                        do ii=max(1,ic-qi1),min(nlon,ic+qi2)

                           !            if (api_us(ii,jj,k)>0.) then

                           do ipat=2,npat
                              ncount = ncount + 1
                              usdum  (k) = usdum  (k) + api_us  (ii,jj,k)
                              tempdum(k) = tempdum(k) + api_temp(ii,jj,k)
                           enddo
                        !                         endif
                        enddo
                     enddo

                     usdum  (k) = usdum  (k)/(float(ncount) + 1.e-10)
                     tempdum(k) = tempdum(k)/(float(ncount) + 1.e-10)

                  enddo

                  !--- interpolation to BRAMS grid
                  do k=mzg,1,-1
                     kk = minloc( abs(slz_us(0:n4us)-slz(k)), 1) - 1

                     globalsoilwater (k,i,j,2:npat) = usdum  (kk+1) !*slmsts(nsoil), vol soil moisture
                     globalsoilenergy(k,i,j,2:npat) = tempdum(kk+1) !- this is actually soil temp

                  enddo

                  do k=1,mzg
                     do ipat=2,npat
                        nsoil = nint(globalsoiltext(k,i,j,ipat))
                        globalsoilwater(k,i,j,ipat) = max(soilcp(nsoil)              , &
                           min(globalsoilwater(k,i,j,ipat), &
                           slmsts(nsoil)))
                     enddo
                  enddo
                  !---smoothness
                  if(USE_SMOOTH_PROF > 0) then
                     ipat=2
                     do k=1,mzg
                        ncount = 0 ; tend1d=0.
                        do kk= max(1,k-USE_SMOOTH_PROF),min(mzg,k+USE_SMOOTH_PROF)
                           ncount =ncount +1
                           tend1d(1)  = tend1d(1)  +  globalsoilwater (kk,i,j,ipat)
                           tend1d(2)  = tend1d(2)  +  globalsoilenergy(kk,i,j,ipat)
                        enddo
                        tend2d(k,1:2)  = tend1d(1:2) /ncount
                     enddo
                     !--- get the final/smoother fields for all land patches
                     do k=mzg,1,-1
                        globalsoilwater (k,i,j,:)=tend2d(k,1)
                        globalsoilenergy(k,i,j,:)=tend2d(k,2)
                     enddo
                  endif
                  !if(maxval(usdum)>0.) print*,"soilM=",maxval(usdum),minval(globalsoilwater(:,i,j,:))


               endif ! ic>=0 .and. jc>=0
            enddo ! loop i
         enddo  ! loop j
         ! --- adjust soil moisture, if desired
         if(factorsm .ne. 1.) globalsoilenergy = globalsoilenergy * factorsm
      else

         stop 'unknown soil moisture data type'

      endif

      !-- scattering local data
      !-- 1) soil water
      call mk_4_buff(globalsoilwater(:,:,:,:), soil_water(:,:,:,:), &
         mzg, nnxp(ifm), nnyp(ifm), npat, mzg, n2, n3, npat, ia, iz, ja, jz)

      if (trim(pref) .ne. 'GFS.SOIL:UMID_TEMP.' .and. &
         trim(pref) .ne. 'ERA5.SOIL:UMID_TEMP.'.and. &
         usdata_in(ipref-10:ipref) .ne. 'YYYYMMDD.nc' ) then
         !----- recalculate soil_energy
         c1 = 0.5*cpi

         do j=1,n3
            do i=1,n2

               k2      = lpw(i,j)
               pis     = c1*(pi0(k2-1,i,j) + pi0(k2,i,j) + pp(k2-1,i,j) + pp(k2,i,j))
               airtemp = theta(k2,i,j)*pis

               do ipat=2,npat
                  do k=1,mzg
                     nsoil = nint(soil_text(k,i,j,ipat))
                     soil_energy(k,i,j,ipat) = (airtemp - 273.15 + stgoff(k))*   &
                        (slcpd(nsoil) + soil_water(k,i,j,ipat)*4.186e6)      + &
                        soil_water(k,i,j,ipat)*3.34e8
                  enddo
               enddo
            enddo
         enddo

         if (mchnum==master_num) then
            open(2, status='unknown', form='unformatted', access='direct', &
               recl=4*nnxp(ifm)*nnyp(ifm)*mzg*npat, file=usmodel(1:len_trim(usmodel)))
            write(unit=2,rec=1) globalsoilwater  !soil moisture
            close(2)
         endif


         !- GFS/ERA5 soil moisture and soil temperature
      elseif (trim(pref) == 'GFS.SOIL:UMID_TEMP.' .or. trim(pref) == 'ERA5.SOIL:UMID_TEMP.') then


         !-- scattering local data
         !-- 2) soil temperature (using soil_energy array)
         call mk_4_buff(globalsoilenergy(:,:,:,:), soil_energy(:,:,:,:), &
            mzg, nnxp(ifm), nnyp(ifm), npat, mzg, n2, n3, npat, ia, iz, ja, jz)

         do j=1,n3
            do i=1,n2
               !k2      = lpw(i,j)
               !pis     = c1*(pi0(k2-1,i,j) + pi0(k2,i,j) + pp(k2-1,i,j) + pp(k2,i,j))
               !airtemp = theta(k2,i,j)*pis

               !--- land
               do ipat=2,npat
                  do k=1,mzg
                     nsoil = nint(soil_text(k,i,j,ipat))

                     !soiltemp = airtemp -273.15 ! using air temperature for initial soil temperature

                     soiltemp = soil_energy(k,i,j,ipat)-273.15 !here soil_energy is actually soil temperature (k)

                     !---temporary restriction
                     soiltemp = max(10., min(soiltemp, 50.)) !Celsius

                     soil_energy(k,i,j,ipat) = (soiltemp)*(slcpd(nsoil) + soil_water(k,i,j,ipat)*4.186e6)+ &
                        soil_water(k,i,j,ipat)*3.34e8


                     !if(soiltemp<20. .or. soiltemp> 40.) print*,"soil=",k,soiltemp,soil_energy(k,i,j,ipat)

                  enddo
               enddo

               !--- ocean (ipat=1)
               soil_energy(1:mzg,i,j,1) = 334000. + 4186. * (0.5*(seatf(i,j) + seatp(i,j))-273.15)


            enddo
         enddo
         !- this is no longer necessary since the model will always reinitiate the soil moist/temp fields
         !- when it starts.
         !if (mchnum==master_num) then
         !    open(2, status='unknown', form='unformatted', access='direct', &
         !         recl=4*nnxp(ifm)*nnyp(ifm)*mzg*npat, file=usmodel(1:len_trim(usmodel)))
         !    write(unit=2,rec=1) globalsoilwater  !soil moisture
         !    write(unit=2,rec=1) globalsoilenergy !soil temp
         !    close(2)
         !endif


      elseif (usdata_in(ipref-10:ipref)=='YYYYMMDD.nc') then  !DSM - Lendo a temperatura e umidade do solo proveniente do JULES

         !-- 2) soil temperature (using soil_energy array)
         call mk_4_buff(globalsoilenergy(:,:,:,:), soil_energy(:,:,:,:), &
            mzg, nnxp(ifm), nnyp(ifm), npat, mzg, n2, n3, npat, ia, iz, ja, jz)

         do j=1,n3
            do i=1,n2
               !--- land
               do ipat=2,npat
                  do k=1,mzg
                     nsoil = nint(soil_text(k,i,j,ipat))
                     if (soil_energy(k,i,j,ipat)>223. .and. soil_energy(k,i,j,ipat)<323.) then
                        soiltemp = soil_energy(k,i,j,ipat)-273.15        !- this is actually soil temp
                     else
                        soiltemp = 0.5*(seatf(i,j) + seatp(i,j))-273.15
                     endif

                     soil_energy(k,i,j,ipat) = (soiltemp)*(slcpd(nsoil) + soil_water(k,i,j,ipat)*4.186e6)+ &
                        soil_water(k,i,j,ipat)*3.34e8
                  enddo
               enddo
               !--- ocean (ipat=1)
               soil_energy(1:mzg,i,j,1) = 334000. + 4186. * (0.5*(seatf(i,j) + seatp(i,j))-273.15)
            enddo
         enddo
         deallocate(t_soil,t_soilJ)

      endif

      deallocate(api_us, usdum,  tempdum, prlat, prlon,api_temp)



   end subroutine soilMoistureInit


   !!$  ! Recreating Global Information (Gathering data)
   !!$  SUBROUTINE gatherData2D(idim_type, varn, ifm, mzg, nnxp, nnyp, npat, nwave, &
   !!$       nmachs, mchnum, mynum, master_num,                                   &
   !!$       localData2D, globalData2D)
   !!$
   !!$    USE ReadBcst, ONLY: &
   !!$         PreProcAndGather, & ! Subroutine
   !!$         LocalSizesAndDisp   ! Subroutine
   !!$    USE mem_grid, ONLY : &
   !!$         GlobalSizes         ! Subroutine
   !!$    USE ParLib, ONLY: &
   !!$         parf_bcast ! Subroutine
   !!$
   !!$    IMPLICIT NONE
   !!$    INCLUDE "i8.h"
   !!$    ! Arguments:
   !!$    INTEGER, INTENT(IN)           :: idim_type, ifm, mzg, nnxp, nnyp, npat, &
   !!$         nwave, nmachs, mchnum, mynum, master_num
   !!$    CHARACTER(LEN=16), INTENT(IN) :: varn
   !!$    REAL, INTENT(IN)              :: localData2D(:,:)
   !!$    REAL, INTENT(OUT)             :: globalData2D(:,:)
   !!$    ! Local Variables:
   !!$    CHARACTER(LEN=16)  :: localVarn
   !!$    INTEGER            :: ierr
   !!$    INTEGER, PARAMETER :: idim_type_min = 2
   !!$    INTEGER, PARAMETER :: idim_type_max = 7
   !!$    INTEGER            :: il1(nmachs)
   !!$    INTEGER            :: ir2(nmachs)
   !!$    INTEGER            :: jb1(nmachs)
   !!$    INTEGER            :: jt2(nmachs)
   !!$    INTEGER            :: localSize(nmachs,idim_type_min:idim_type_max)
   !!$    INTEGER            :: disp(nmachs,idim_type_min:idim_type_max)
   !!$    INTEGER            :: maxLocalSize
   !!$    INTEGER            :: sizeGathered(idim_type_min:idim_type_max)
   !!$    INTEGER            :: maxSizeGathered
   !!$    INTEGER            :: sizeFullField(idim_type_min:idim_type_max)
   !!$    INTEGER            :: maxsizeFullField
   !!$    INTEGER            :: globalSize(idim_type_min:idim_type_max)
   !!$    REAL, ALLOCATABLE  :: localChunk(:)
   !!$    REAL, ALLOCATABLE  :: gathered(:)
   !!$    REAL, ALLOCATABLE  :: fullField(:)
   !!$
   !!$    ! Recreating Global information about Soil Water
   !!$    ! grid dependent, field independent constants for gather and unpacking
   !!$    ! as a function of idim_type
   !!$    CALL LocalSizesAndDisp(ifm, il1, ir2, jb1, jt2, localSize, disp)
   !!$    maxLocalSize = MAXVAL(localSize(mynum,:))
   !!$    ALLOCATE(localChunk(maxLocalSize), stat=ierr)
   !!$    IF (ierr/=0) THEN
   !!$       CALL fatal_error("Error allocating localChunk (gatherData)")
   !!$    ENDIF
   !!$    CALL CopyLocalChunk(localData2D(1,1), localChunk, &
   !!$         LocalSize(mynum,idim_type))
   !!$    sizeGathered(:) = disp(nmachs,:) + localSize(nmachs,:)
   !!$    maxSizeGathered = MAXVAL(sizeGathered)
   !!$    ALLOCATE(gathered(maxSizeGathered), stat=ierr)
   !!$    IF (ierr/=0) THEN
   !!$       CALL fatal_error("Error allocating gathered (gatherData)")
   !!$    ENDIF
   !!$    ! grid dependent field sizes as a function of idim_type
   !!$    CALL GlobalSizes(ifm, nmachs, nwave, globalSize)
   !!$    IF (mchnum==master_num) THEN
   !!$       sizeFullField(:) = globalSize(:)
   !!$    ELSE
   !!$       sizeFullField(:) = 1
   !!$    END IF
   !!$    maxSizeFullField = MAXVAL(sizeFullField)
   !!$    ALLOCATE(fullField(sizeFullField(idim_type)), stat=ierr)
   !!$    IF (ierr/=0) THEN
   !!$       CALL fatal_error("Error allocating fullField (gatherData)")
   !!$    ENDIF
   !!$    localVarn = trim(varn)
   !!$    CALL PreProcAndGather(.FALSE., ifm, idim_type, localVarn, &
   !!$         il1, ir2, jb1, jt2, localSize, disp,                 &
   !!$         localSize(mynum,idim_type), LocalChunk,              &
   !!$         sizeGathered(idim_type), gathered,                   &
   !!$         sizeFullField(idim_type), fullField                  )
   !!$
   !!$    IF (mchnum==master_num) THEN
   !!$       globalData2D = RESHAPE(fullField, (/nnxp, nnyp/))
   !!$    ENDIF
   !!$
   !!$    call parf_bcast(globalData2D, int(nnxp,i8), int(nnyp,i8), master_num)
   !!$
   !!$    DEALLOCATE(fullField)
   !!$    DEALLOCATE(gathered)
   !!$    DEALLOCATE(localChunk)
   !!$
   !!$  END SUBROUTINE gatherData2D
   !!$
   !!$  ! Recreating Global Information (Gathering data)
   !!$  SUBROUTINE gatherData4D(idim_type, varn, ifm, mzg, nnxp, nnyp, npat, nwave, &
   !!$       nmachs, mchnum, mynum, master_num,                                   &
   !!$       localData4D, globalData4D)
   !!$
   !!$    USE ReadBcst, ONLY: &
   !!$         PreProcAndGather, & ! Subroutine
   !!$         LocalSizesAndDisp   ! Subroutine
   !!$    USE mem_grid, ONLY : &
   !!$         GlobalSizes         ! Subroutine
   !!$    USE ParLib, ONLY: &
   !!$         parf_bcast ! Subroutine
   !!$
   !!$    IMPLICIT NONE
   !!$    INCLUDE "i8.h"
   !!$    ! Arguments:
   !!$    INTEGER, INTENT(IN)           :: idim_type, ifm, mzg, nnxp, nnyp, npat, &
   !!$         nwave, nmachs, mchnum, mynum, master_num
   !!$    CHARACTER(LEN=16), INTENT(IN) :: varn
   !!$    REAL, INTENT(IN)              :: localData4D(:,:,:,:)
   !!$    REAL, INTENT(OUT)             :: globalData4D(:,:,:,:)
   !!$    ! Local Variables:
   !!$    CHARACTER(LEN=16)  :: localVarn
   !!$    INTEGER            :: ierr
   !!$    INTEGER, PARAMETER :: idim_type_min = 2
   !!$    INTEGER, PARAMETER :: idim_type_max = 7
   !!$    INTEGER            :: il1(nmachs)
   !!$    INTEGER            :: ir2(nmachs)
   !!$    INTEGER            :: jb1(nmachs)
   !!$    INTEGER            :: jt2(nmachs)
   !!$    INTEGER            :: localSize(nmachs,idim_type_min:idim_type_max)
   !!$    INTEGER            :: disp(nmachs,idim_type_min:idim_type_max)
   !!$    INTEGER            :: maxLocalSize
   !!$    INTEGER            :: sizeGathered(idim_type_min:idim_type_max)
   !!$    INTEGER            :: maxSizeGathered
   !!$    INTEGER            :: sizeFullField(idim_type_min:idim_type_max)
   !!$    INTEGER            :: maxsizeFullField
   !!$    INTEGER            :: globalSize(idim_type_min:idim_type_max)
   !!$    REAL, ALLOCATABLE  :: localChunk(:)
   !!$    REAL, ALLOCATABLE  :: gathered(:)
   !!$    REAL, ALLOCATABLE  :: fullField(:)
   !!$
   !!$    ! Recreating Global information about Soil Water
   !!$    ! grid dependent, field independent constants for gather and unpacking
   !!$    ! as a function of idim_type
   !!$    CALL LocalSizesAndDisp(ifm, il1, ir2, jb1, jt2, localSize, disp)
   !!$    maxLocalSize = MAXVAL(localSize(mynum,:))
   !!$    ALLOCATE(localChunk(maxLocalSize), stat=ierr)
   !!$    IF (ierr/=0) THEN
   !!$       CALL fatal_error("Error allocating localChunk (gatherData)")
   !!$    ENDIF
   !!$    CALL CopyLocalChunk(localData4D(1,1,1,1), localChunk, &
   !!$         LocalSize(mynum,idim_type))
   !!$    sizeGathered(:) = disp(nmachs,:) + localSize(nmachs,:)
   !!$    maxSizeGathered = MAXVAL(sizeGathered)
   !!$    ALLOCATE(gathered(maxSizeGathered), stat=ierr)
   !!$    IF (ierr/=0) THEN
   !!$       CALL fatal_error("Error allocating gathered (gatherData)")
   !!$    ENDIF
   !!$    ! grid dependent field sizes as a function of idim_type
   !!$    CALL GlobalSizes(ifm, nmachs, nwave, globalSize)
   !!$    IF (mchnum==master_num) THEN
   !!$       sizeFullField(:) = globalSize(:)
   !!$    ELSE
   !!$       sizeFullField(:) = 1
   !!$    END IF
   !!$    maxSizeFullField = MAXVAL(sizeFullField)
   !!$    ALLOCATE(fullField(sizeFullField(idim_type)), stat=ierr)
   !!$    IF (ierr/=0) THEN
   !!$       CALL fatal_error("Error allocating fullField (gatherData)")
   !!$    ENDIF
   !!$    localVarn = trim(varn)
   !!$    CALL PreProcAndGather(.FALSE., ifm, idim_type, localVarn, &
   !!$         il1, ir2, jb1, jt2, localSize, disp,                 &
   !!$         localSize(mynum,idim_type), LocalChunk,              &
   !!$         sizeGathered(idim_type), gathered,                   &
   !!$         sizeFullField(idim_type), fullField                  )
   !!$
   !!$    IF (mchnum==master_num) THEN
   !!$       globalData4D = RESHAPE(fullField, (/mzg, nnxp, nnyp, npat/))
   !!$    ENDIF
   !!$
   !!$    call parf_bcast(globalData4D, &
   !!$         int(mzg,i8), int(nnxp,i8), int(nnyp,i8), int(npat,i8), master_num)
   !!$
   !!$    DEALLOCATE(fullField)
   !!$    DEALLOCATE(gathered)
   !!$    DEALLOCATE(localChunk)
   !!$
   !!$  END SUBROUTINE gatherData4D
   subroutine StoreNamelistFileAtSoilMoisture(oneNamelistFile)
      implicit none
      type(namelistFile), pointer :: oneNamelistFile
      soil_moist = oneNamelistFile%soil_moist
      soil_moist_fail = oneNamelistFile%soil_moist_fail
      usdata_in = oneNamelistFile%usdata_in
      usmodel_in = oneNamelistFile%usmodel_in
   end subroutine StoreNamelistFileAtSoilMoisture
end module soilMoisture

 !
 ! prlatlon
 !----------------------------------------------------------------
 ! SUB-ROTINA QUE ESTABELECE LATITUDES E LONGITUDES DOS PONTOS DE
 ! GRADE DO CAMPO DE PRECIPITACAO
subroutine apiPrlatlon(nlon, nlat, prlat, prlon, ilatn, ilonn, latni, lonni)
   implicit none
   ! Arguments:
   integer, intent(in) :: nlon, nlat
   real, intent(out)   :: prlat(nlon,nlat) !(nlon,nlat)
   real, intent(out)   :: prlon(nlon,nlat) !(nlon,nlat)
   real, intent(in)    ::  ilatn, ilonn, latni, lonni
   ! Local Variables:
   integer :: i, j

   do j=1,nlat
      do i=1,nlon
         prlon(i,j) = lonni + (i-1)*ilonn
         prlat(i,j) = latni + (j-1)*ilatn
      enddo
   enddo

end subroutine apiPrlatlon

 !----------------------------------------------------------------
 ! interpolacao
 !----------------------------------------------------------------
 ! SUB-ROTINA QUE REALIZA INTERPOLACAO ENTRE GRADES (RAMS E UMIDADE DO SOLO)
subroutine interpolacao(glon, glat, nlon, nlat, prlat, prlon, &
   i1, i2, ic, j1, j2, jc)

   implicit none
   ! Arguments:
   integer, intent(out) :: i1, i2, ic, j1, j2, jc
   real, intent(in)     :: glat, glon
   integer, intent(in)  :: nlon, nlat
   real, intent(in)     :: prlat(nlon,nlat)
   real, intent(in)     :: prlon(nlon,nlat)
   ! Local Variables
   real    :: diffx1, diffx2, diffy1, diffy2
   integer :: i, j

   do i=1,nlon
      if (glon<=prlon(i,1)) exit
   enddo
   i2 = i
   i1 = i-1

   do j=1,nlat
      if (glat<=prlat(1,j)) exit
   enddo
   j2 = j
   j1 = j-1

   diffx1 =   glon - prlon(i1,j1)
   diffx2 = -(glon - prlon(i1,j2))
   diffy1 =   glat - prlat(i1,j1)
   diffy2=  -(glat - prlat(i2,j1))

   jc = j1
   ic = i1
   if (diffx1>diffx2) ic = i2
   if (diffy1>diffy2) jc = j2

   if (i1<1 .or. i1>nlon .or. j1<1 .or. j1>nlat) then
      ic = -9999
      jc = -9999
   endif

end subroutine interpolacao

 !----------------------------------------------------------------

subroutine changeDay(idate1, imonth1, iyear1, INT_DIF_TIME, &
   idate2, imonth2, iyear2)
   implicit none
   ! Arguments:
   integer, intent(in)  :: idate1, imonth1, iyear1, INT_DIF_TIME
   integer, intent(out) :: idate2, imonth2, iyear2
   ! Local Variables:
   integer :: i, increm, DMES(12)

   ! Initiate DMES
   DMES = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)

   iyear2  = iyear1
   imonth2 = imonth1
   idate2  = idate1

   increm  = 1

   if (INT_DIF_TIME<1) increm = increm*(-1)

   do i=1,ABS(INT_DIF_TIME)
      idate2 = idate2 + increm
      if (idate2<1) then
         imonth2 = imonth2 + increm
         if (imonth2<1) then
            imonth2 = 12
            iyear2  = iyear2 - 1
         endif
         idate2 = DMES(imonth2)
      elseif (idate2>DMES(imonth2)) then
         imonth2 = imonth2 + increm
         if (imonth2>12) then
            imonth2 = 1
            iyear2  = iyear2 + 1
         endif
         idate2 = 1
      endif
   enddo

end subroutine changeDay


!-srf END MODULE soilMoisture

!============================================================================

subroutine Swap32(A, N)
   !
   !      REVERSE ORDER OF BYTES IN INTEGER*4 WORD, or REAL*4
   !
   implicit none
   ! Arguments:
   integer, intent(in)            :: n
   integer(kind=4), intent(inout) :: a(n)
   ! Local Varaibles:
   character(len=1) :: jtemp(4)
   character(len=1) :: ktemp
   !
   ! Local variables
   integer :: i, itemp

   equivalence (jtemp(1), itemp)
   !
   save
   !
   do i=1,n
      itemp    = a(i)
      ktemp    = jtemp(4)
      jtemp(4) = jtemp(1)
      jtemp(1) = ktemp
      ktemp    = jtemp(3)
      jtemp(3) = jtemp(2)
      jtemp(2) = ktemp
      a(i)     = itemp
   enddo

end subroutine Swap32

!DSM - Para trocar por exempo AAAAMMDD por 20050425 {
subroutine Replace_Text (s,text,rep,outs)
   implicit none
   character(*)        :: s,text,rep
   character(LEN(s)+100) :: outs     ! provide outs with extra 100 char len
   integer             :: i, nt, nr

   outs = s ; nt = LEN_TRIM(text) ; nr = LEN_TRIM(rep)
   do
      i = INDEX(outs,text(:nt)) ; if (i == 0) exit
      outs = outs(:i-1) // rep(:nr) // outs(i+nt:)
   end do
end subroutine Replace_Text
!DSM }

