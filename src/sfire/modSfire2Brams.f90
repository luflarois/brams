module modSfire2Brams

    use mem_grid, only: &
        nxtnest, nnxp, nnyp, npatch, nzp, nzg, ipm, jpm, platn, plonn, &
        iyear1, imonth1, idate1, ihour1, xtn, ytn, deltaxn,  deltayn, time, &
        deltax, deltay, zmn

    use node_mod, only: &
        mynum,         &
        mchnum,        &
        master_num

    use mem_mksfc, only: &
        dato, datp, cdato, nump, numpind, numpind1, numpind2, &
        ptable, datq_patch, npq

    use aer1_list, only: &
        nspecies

    use mem_aer1, only: aer1_g ! aer1_g(imode,ispc,1)%sc_p)
    !use mem_chem, only: chem1_g !chem1_g(ispc)%sc_p

    use node_mod, only: &
        nodei0, & ! intent(in)
        nodej0    ! intent(in)

    use ParLib, only: parf_barrier, parf_bcast
    
    implicit none
    private


    include "i8.h"

    logical, parameter :: vegDump =.false.
    character(len=256), parameter :: h5name = " "
    integer, parameter :: maxmiss=1000, ifm=1
    integer, parameter :: nveg=17, nveg_olson=97, nveg_agreg = 4
    real, parameter :: fx=1.0
    integer, parameter :: nspecies_sfire = 5
    integer, parameter :: max_sfire_points_by_brams = 400
    
    real,    dimension(7,nveg) :: queiparms
    integer, allocatable :: vegType(:,:,:,:)
    real,    dimension(nveg)   :: fc,bas_queiparms
    real :: burntArea
    real, dimension(:,:), allocatable :: gridVolume
    real, dimension(:,:,:,:), allocatable :: glatp_l, glonp_l
    integer :: npx, npy
    real, dimension(:), allocatable :: lat, lon
    character(len=256) :: vegFile
    logical :: test_only
    real :: time_2_test

    integer, parameter :: catb(17) = (/&
          1, 1, 1, 1                   & !floresta: 1 a 4
        , 1, 3, 3, 3, 3                & !cerrado/woody savanna :6 a 9
        , 4, 4, 4, 4, 4, 0, 4, 0     /)  !pastagem/lavouras: 10 ...

    real, parameter, dimension(0:nveg_agreg) :: flaming = (/ &
       0.00 &  ! 
      ,0.30 &  ! % biomass burned at flaming phase : floresta              igbp 1 a 4
      ,0.30 &  ! % biomass burned at flaming phase : floresta              igbp 1 a 4
      ,0.50 &  ! % biomass burned at flaming phase : cerrado/woody savanna igbp 5 a 9
      ,0.00 /) ! % biomass burned at flaming phase : pastagem/lavoura:     igbp 10 a 17

    !#Fire info data
    real, allocatable :: sfire_info_area(:)
    !# burn area
    real, allocatable :: sfire_info_time(:)
    !# burn time
    real, allocatable :: sfire_info_frp (:)
    !Fire radiate power
    real, allocatable :: sfire_info_lat (:)
    !Latitude of point
    real, allocatable :: sfire_info_lon (:)
    !Longitude of point
    real, allocatable :: flam_frac(:,:)
    real, allocatable :: mean_frp (:,:)
    real, allocatable :: std_frp  (:,:)
    real, allocatable :: mean_size(:,:)
    real, allocatable :: std_size (:,:)
    integer, allocatable :: valid_ij(:,:)

    integer, allocatable :: count_in_grid(:,:)
    real, allocatable :: qsc(:,:,:,:)
    real, allocatable :: b_fires_area(:,:,:)
    real, allocatable :: b_fires_time(:,:,:)
    real, allocatable :: b_fires_frp(:,:,:)
    real, allocatable :: b_fires_lat(:,:,:)
    real, allocatable :: b_fires_lon(:,:,:)

    !-------------------------------------------------------------------------------
    data queiparms/  &
    ! dados de Andreae & Merlet(A&M), Ward et al., Kauffman (unpubl.), Sinha et al, 2003, CAPOS 2004
    ! Palacios-Orueta et al, Barbosa e Fearnside
    !--------------------------------------------------------------------------------------------------------
    !fator  !Biomassa   !EF    !EF     !EF    !EF    !EF   !IGBP Land Cover              ! observation
    !de com-!acima do   !g/kg  !g/kg   !g/kg  !g/kg  !g/kg !Legend and                   ! reference
    !bustao !solo(kg/m2)!CO2   !CO     !BC    !OC   !SO2  !description                  ! 
    !-------------------------------------------------------------------------------------------------------
    0.50,	 29.24,     1568., 106.0,  .56,  9.1,  1.00, &! 1 Evergreen Needleleaf Fores !A&M
    0.50,	 29.24,     1527., 128.8,  .63,  5.2,  0.57, &! 2 Evergreen Broadleaf Forest !CAPOS 2004
    0.43,	 12.14,     1639.,  99.0,  .46,  3.2,  0.35, &! 3 Deciduous Needleleaf Forest!  
    0.43,	 12.14,     1639.,  99.0,  .46,  3.2,  0.35, &! 4 Deciduous Broadleaf Forest !  
    0.43,	 12.14,     1639.,  99.0,  .46,  3.2,  0.35, &! 5 Mixed Forest		   !		     
    0.87,	  7.40,     1700.,  68.0,  .46,  3.2,  0.35, &! 6 Closed Shrublands =Caatinga!Kauffman, Sinha et al, 2003
    0.72,	  0.86,     1700.,  68.0,  .46,  3.2,  0.35, &! 7 Open Shrublands	     !Sinha et al, 2003
    0.45,	  10.0,     1700.,  68.0,  .46,  3.2,  0.35, &! 8 Woody Savannas	     !Sinha et al, 2003
    0.52,	  6.20,     1700.,  68.0,  .46,  3.2,  0.35, &! 9 Savannas		     !Sinha et al, 2003 	 
    1.00,	  0.71,     1700.,  68.0,  .46,  3.2,  0.35, &!10 Grasslands		     !Sinha et al, 2003   
    0.40,	  3.80,     1700.,  68.0,  .46,  3.2,  0.35, &!11 Permanent Wetlands	
    0.40,	  3.80,     1700.,  68.0,  .46,  3.2,  0.35, &!12 Croplands	 Palacios-orueta.
    0.40,	  3.80,     1700.,  68.0,  .46,  3.2,  0.35, &!13 Urban and Built-Up	
    0.40,	  3.80,     1700.,  68.0,  .46,  3.2,  0.35, &!14 Cropland/Natural Veg.Mosaic Palacios-orueta
    0.00,	  0.00,     0000.,  00.0,  .00,  0.0,  0.00, &!15 Snow and Ice
    0.84,	  1.00,     1723.,  57.5,  .46,  3.2,  0.35, &!16 Barren or Sparsely Vegetated 
    0.00,	  0.00,     0000.,  00.0,  .00,  0.0,  0.00  /!17 Water Bodies
    !-------------------------------------------------------------------------------------------------------

    real,    dimension(10,nspecies_sfire) :: emission_factor

    data emission_factor/  &
    !---------------------------------------------------------------------------------------------------
    ! Emission Factor
    !---------------------------------------------------------------------------------------------------
    ! TropFor! Extra !Savanna! pasture | Biofuel! CharcMak ! CharcBurn ! AgResid ! Laboratory ! Molec. !
    !	 !tropF  !	 ! cropland|	    !  	       !           !         !            ! Weight !
    !---------------------------------------------------------------------------------------------------
      8.3,    15.7,  7.5,	   7.5,    6.0,       2.1,        1.6,        7.4,     5.1,    -9999, &! PM2.5
      11.8,   17.6,  8.2,	   8.2,    6.9,       4.3,        2.4,        10.8,    7.2,    -9999, &! PM10 (old TPM)
      6.0,    8.3,   4.0,	   4.0,    3.4,       -9999,      3.3,        4.4,     3.6,    12.00, &! TC
      4.3,    7.7,   3.3,	   3.3,    3.1,       -9999,      4.8,        4.0,     2.5,    12.00, &! OC
      0.56,   0.58,  0.61,     0.61,   0.54,      -9999,      1.50,       0.41,    1.08,   12.00 / ! BC

    character(len=4), dimension(nspecies_sfire) :: sfire_species_name = (/"PM25","PM10","TC  ","OC  ","BC  "/)

      public :: readModisVeg, get_emission_in_global_brams_grid, nspecies_sfire &
                , sfire_species_name, sfire_info_area, sfire_info_time, sfire_info_frp, sfire_info_lat &
                , sfire_info_lon, flam_frac, mean_frp, std_frp, mean_size, std_size, qsc,count_in_grid &
                , alloc_dealloc_sfire_info, valid_ij, comm_aer_data, testModSfire2Brams, test_only &
                , fill_sfire_info_and_get_emissions

contains

    function alloc_dealloc_sfire_info(size_info) result (ok)
        implicit none

        integer, intent(in), optional :: size_info
        logical :: ok
        if (present(size_info)) then
            if (.not. allocated(sfire_info_area)) allocate(sfire_info_area(size_info))
            if (.not. allocated(sfire_info_time)) allocate(sfire_info_time (size_info))
            if (.not. allocated(sfire_info_frp )) allocate(sfire_info_frp (size_info))
            if (.not. allocated(sfire_info_lat )) allocate(sfire_info_lat (size_info))
            if (.not. allocated(sfire_info_lon )) allocate(sfire_info_lon (size_info))  

        else
            if (allocated(sfire_info_area)) deallocate(sfire_info_area)
            if (allocated(sfire_info_time)) deallocate(sfire_info_time )
            if (allocated(sfire_info_frp )) deallocate(sfire_info_frp )
            if (allocated(sfire_info_lat )) deallocate(sfire_info_lat )
            if (allocated(sfire_info_lon )) deallocate(sfire_info_lon )       
        end if     
        ok = .true.

    end function alloc_dealloc_sfire_info

    subroutine readModisVeg(rtgt,razaoBSF)
        implicit none

        character(len=*), parameter :: iaction = "veg"
        real, intent(in) :: rtgt(nnxp(1), nnyp(1))
        real, intent(in) :: razaoBSF

        integer :: datq,lsp  &
        ,datsoil,soil_count,nlev_soil,nmiss      &
        ,jr,jp,ir,ip,iqv,ing,maxdq,jng,jng1,ngstor1    &
        ,ngstor2,npatpixs,nwat,ipat,idatq,isoil,jsoil,k      &
        ,no_veg ,iblksizo_veg ,isbego_veg ,iwbego_veg

        real :: deltallo_veg ,deltallo_soil ,deltallo_ndvi  &
        ,offlat_veg   ,offlat_soil   ,offlat_ndvi    &
        ,offlon_veg   ,offlon_soil   ,offlon_ndvi    &
        ,xp,yp,fracwat,plpp

        character(len=256) :: fnmiss(maxmiss)
        character(len=256) :: ivegtfn 

        integer :: i,j
        real, allocatable, dimension(:) :: xtn_l, ytn_l



        if (mchnum==master_num .and. time==0.0) then

            call setup_fireVeg(rtgt,razaoBSF)

            ivegtfn = trim(vegFile)

            npx = deltaxn(ifm)/1000*nnxp(ifm)
            npy = deltayn(ifm)/1000*nnyp(ifm)

            npq = 1

            !Allocate arrays that need npq dimension',npq
            allocate (glatp_l(npq,npq,npx,npy),  &
            glonp_l(npq,npq,npx,npy),       &
            datq_patch(npq,npq,npx,npy),  &
            datp(npq,npq,npx,npy),xtn_l(npx),ytn_l(npy) &
            , lat(npy), lon(npx))

            xtn_l(1) = -(npx/2*1000.)
            do i = 2,npx
                xtn_l(i) = xtn_l(i-1)+1000.0
            end do 
            ytn_l(1) = -(npy/2*1000.)
            do i = 2,npy
                ytn_l(i) = ytn_l(i-1)+1000.0
            end do             

            !Allocate the array that will have the vegetation type
            allocate (vegType(npq,npq,npx,npy))

            call patch_latlon_fix1km(npq, npx, npy, xtn_l, ytn_l, 1000., 1000., platn(ifm), plonn(ifm))


            lat = glatp_l(1,1,1,:)
            lon = glonp_l(1,1,:,1)

            !print *, 'Reading header of MODIS vegetation types'
            call read_header(ivegtfn,iblksizo_veg,no_veg,isbego_veg  &
            ,iwbego_veg,offlat_veg,offlon_veg,deltallo_veg,'veg', h5name)

            call fill_datp_fix1km(npx,npy,no_veg,iblksizo_veg,isbego_veg,iwbego_veg  &
            ,platn(ifm),plonn(ifm),offlat_veg,offlon_veg,deltallo_veg,ivegtfn,iaction  &
            ,nmiss,fnmiss, h5name)

            vegType = datp

        end if

    end subroutine readModisVeg

    subroutine fill_datp_fix1km(n2,n3,no,iblksizo,isbego,iwbego  &
        ,platn,plonn,offlat,offlon,deltallo,ofn,iaction,nmiss,fnmiss, h5name)
   
     use mem_mksfc, only: &
          dato, datp, cdato, nump, numpind, numpind1, numpind2, &
          ptable
   
     implicit none
   
     include "files.h"
   
     integer, parameter :: maxmiss=1000
   
     character(len=f_name_length), intent(IN)    :: ofn
     character(len=f_name_length), intent(INOUT) :: fnmiss(maxmiss)
     character(len=*), intent(IN)                :: iaction
     character(len=*), intent(IN)                :: h5name
   
     integer :: nmiss
   
     integer :: n2,n3,no,iblksizo,isbego,iwbego,isoc,iwoc,iofr  &
          ,isocpt,isocpo,iwocph,iwocpt,iwocpo,lb,io,jo  &
          ,ir,jr,ip,jp,ind,nc3,nc2,j3d,j2d,j1d,ind1,ind2,io1,jo1  &
          ,ifile_max,jfile_max,ifile,jfile,missing,ptab,ptab0,idatp,nn
   
     real :: rio,rjo,rno,platn,plonn,offlat,offlon  &
          ,glatp1,glonp1,deltallo,wio2,wjo2,wio1,wjo1
   
   !!$  character :: title1*3,title2*4,title3*80
     character(len=3)   :: title1
     character(len=4)   :: title2
     character(len=f_name_length) :: title3
   
     logical l1,l2
   
   !!$  include 'interface.h'
   
     ! Compute number of files in input dataset that span all latitudes and
     ! longitudes on earth.  Allocate nump and numpind arrays to this size and
     ! initialize to 0.
     ! Allocate ptable array.
   
     rno = float(no)
     ifile_max = 360 / iblksizo
     jfile_max = 180 / iblksizo
   

    allocate (cdato(no,no))
   
     allocate (nump    (ifile_max,jfile_max)  &
          ,numpind (ifile_max,jfile_max)  &
          ,numpind1(ifile_max,jfile_max)  &
          ,numpind2(ifile_max,jfile_max)  &
          ,ptable  (npq*npq*n2*n3))
   
     do jfile = 1,jfile_max
        do ifile = 1,ifile_max
           nump(ifile,jfile) = 0
           numpind(ifile,jfile) = 0
        enddo
     enddo
   
     ! Get file index (ifile,jfile) within full dataset and count number of p
     ! points (nump) that occur in each file
     !print *,'LFR-DBG: n2,n3,npq: ',n2,n3,npq
     do jr = 1,n3
        do jp = 1,npq
           do ir = 1,n2
              do ip = 1,npq
   
                 glatp1 = max(-89.9999,min(89.9999,glatp_l(ip,jp,ir,jr) - offlat))
                 glonp1 = glonp_l(ip,jp,ir,jr) - offlon
                 if (glonp1 .ge.  180.) glonp1 = glonp1 - 360.
                 if (glonp1 .le. -180.) glonp1 = glonp1 + 360.
   
                 ifile = int((glonp1 - float(iwbego)) / float(iblksizo)) + 1
                 jfile = int((glatp1 - float(isbego)) / float(iblksizo)) + 1
   
                 nump(ifile,jfile) = nump(ifile,jfile) + 1
   
              enddo
           enddo
        enddo
     enddo
   
     ! Set up array index values for ptable array
   
     ind = 1
     do jfile = 1,jfile_max
        do ifile = 1,ifile_max
           numpind1(ifile,jfile) = ind
           numpind2(ifile,jfile) = ind
           ind = ind + nump(ifile,jfile)
        enddo
     enddo
   
     ! Fill ptable array
   
     nc3 = n2 * npq * npq
     nc2 = npq * npq
   
     do jr = 1,n3
        j3d = (jr - 1) * nc3
        do ir = 1,n2
           j2d = (ir - 1) * nc2
           do jp = 1,npq
              j1d = (jp - 1) * npq
              do ip = 1,npq
   
                 glatp1 = max(-89.9999,min(89.9999,glatp_l(ip,jp,ir,jr) - offlat))
                 glonp1 = glonp_l(ip,jp,ir,jr) - offlon
                 if (glonp1 .ge.  180.) glonp1 = glonp1 - 360.
                 if (glonp1 .le. -180.) glonp1 = glonp1 + 360.
   
                 ifile = int((glonp1 - float(iwbego)) / float(iblksizo)) + 1
                 jfile = int((glatp1 - float(isbego)) / float(iblksizo)) + 1
   
                 ind = numpind2(ifile,jfile)
                 ptable(ind) = j3d + j2d + j1d + ip
                 numpind2(ifile,jfile) = numpind2(ifile,jfile) + 1
   
              enddo
           enddo
        enddo
     enddo
   
     !print *,'LFR-DBG: Read files and extract data'
   
     do jfile = 1,jfile_max
        do ifile = 1,ifile_max
   
           ind1 = numpind1(ifile,jfile)
           ind2 = numpind2(ifile,jfile)
           !if(ind2 .gt. ind1) print *,'LFR-DBG: ind1,ind2: ',ind2 .gt. ind1
   
   
           if (ind2 .gt. ind1) then
              isoc = (jfile - 1) * iblksizo + isbego
              iwoc = (ifile - 1) * iblksizo + iwbego
   
              ! Construct filename
   
              isocpt = abs(isoc) / 10
              isocpo = abs(isoc) - isocpt*10
              iwocph = abs(iwoc) / 100
              iwocpt = (abs(iwoc) - iwocph * 100) / 10
              iwocpo = abs(iwoc) - iwocph * 100 - iwocpt * 10
   
   
              if (isoc .ge. 0) then
                 write(title1,'(2i1,a1)') isocpt,isocpo,'N'
              else
                 write(title1,'(2i1,a1)') isocpt,isocpo,'S'
              endif
   
   
              if (iwoc .ge. 0) then
                 write(title2,'(3i1,a1)') iwocph,iwocpt,iwocpo,'E'
              else
                 write(title2,'(3i1,a1)') iwocph,iwocpt,iwocpo,'W'
              endif
   
              lb = len_trim(ofn)
              title3 = ofn(1:lb)//title1//title2
   
              lb = len_trim(title3)
   
              inquire(file=title3(1:lb),exist=l1,opened=l2)
              print *,'LFR-DBG: lendo ',title3(1:lb),l1,l2
              if(l1) then
                missing = 0
                print*, 'getting file ', title3(1:lb)
                  call rams_c_open(title3(1:lb)//char(0), 'rb'//char(0))
                  call rams_c_read_char(4, no*no, cdato(1,1))
                  call rams_c_close()
              else
                do nn=1,nmiss
                  if(trim(title3(1:lb)) == trim(fnmiss(nn))) goto 302
                enddo
                nmiss         = nmiss + 1
                fnmiss(nmiss) = title3(1:lb)
   302          continue
                missing       = 1
              endif
   ! ------------------------------------------------------------------------------------------------
   
              do ind=ind1,ind2-1
   
   
                 ptab  = ptable(ind)
                 ptab0 = ptab - 1
                 jr    = ptab0 / nc3 + 1
                 j3d   = (jr - 1) * nc3
                 ir    = (ptab0 - j3d) / nc2 + 1
                 j2d   = (ir - 1) * nc2
                 jp    = (ptab0 - j3d - j2d) / npq + 1
                 j1d   = (jp - 1) * npq
                 ip    = ptab - j3d - j2d - j1d
   
                 glatp1 = max(-89.9999,min(89.9999,glatp_l(ip,jp,ir,jr) - offlat))
                 glonp1 = glonp_l(ip,jp,ir,jr) - offlon
                 if (glonp1 .ge.  180.) glonp1 = glonp1 - 360.
                 if (glonp1 .le. -180.) glonp1 = glonp1 + 360.
   
                 rio = (glonp1 - float(iwoc)) / deltallo + 1.
                 rjo = (glatp1 - float(isoc)) / deltallo + 1.
   
                 if (rio .lt. .9 .or. rio .gt. rno+.1 .or.  &
                      rjo .lt. .9 .or. rjo .gt. rno+.1) then
                    print*, 'rio,rjo out of range ',rio,rjo,ip,jp,ir,jr
                    stop 45
                 endif
   
                 if (missing .eq. 0) then
   
                       io = nint(rio)
                       jo = nint(rjo)
                       idatp = ichar(cdato(io,jo))
                       datp(ip,jp,ir,jr) = float(mod(idatp+256,256))
                    
                 endif
   
              enddo
           endif
        enddo
     enddo
   
     if (iaction == 'ndvi') then
        deallocate (dato)
     else
        deallocate (cdato)
     endif
   
     deallocate(nump,numpind,numpind1,numpind2,ptable)
   
     return
   end subroutine fill_datp_fix1km

    subroutine patch_latlon_fix1km(npq,n2,n3,xt,yt,deltax,deltay,platn,plonn)

        implicit none
      
        integer :: n2,n3,npq
        real :: xt(n2),yt(n3),platn,plonn
      
        integer :: jr,jp,ip,ir
        real :: yp,xp,deltax,deltay,deltaxp,deltayp
      
        ! Fill arrays with offset latitudes and longitudes of all p points
      
        deltaxp = deltax / float(npq)
        deltayp = deltay / float(npq)
        !print *,'LFR-DBG deltaxp,deltayp,npq: ',deltaxp,deltayp,npq
        !print *,'LFR-DBG npq,npq,n2,n3: ',npq,npq,n2,n3,platn,plonn
        do jr = 1,n3
           do jp = 1,npq
              yp = yt(jr) + (float(jp) - .5 * float(npq+1)) * deltayp
              do ir = 1,n2
                 do ip = 1,npq
                    xp = xt(ir) + (float(ip) - .5 * float(npq+1)) * deltaxp

                    call xy_ll(glatp_l(ip,jp,ir,jr),glonp_l(ip,jp,ir,jr)  &
                         ,platn,plonn,xp,yp)

                 enddo
              enddo
           enddo
        enddo
      end subroutine patch_latlon_fix1km



    subroutine setup_fireVeg(rtgt,razaoBSF)
        implicit none

        real, intent(in) :: rtgt(nnxp(1), nnyp(1))
        real, intent(in) :: razaoBSF
        !# razão entre as resolução das grades BRAMS e Sfire
        integer :: i,j,iesp,iUnit, err
        real :: factArea, qarea
        logical :: ex

        namelist /VEG/ vegFile, test_only, time_2_test

        do i=1,nveg
            fc(i)             = queiparms(1,i)        ! fracao de biomassa queimada
            bas_queiparms(i)  = queiparms(2,i)        ! kg/m2
        end do  

        allocate(flam_frac(nnxp(ifm),nnyp(ifm)))
        allocate(mean_frp (nnxp(ifm),nnyp(ifm)))
        allocate(std_frp  (nnxp(ifm),nnyp(ifm)))
        allocate(mean_size(nnxp(ifm),nnyp(ifm)))
        allocate(std_size (nnxp(ifm),nnyp(ifm)))
        allocate(valid_ij (nnxp(ifm),nnyp(ifm)))
        allocate(count_in_grid(nnxp(ifm),nnyp(ifm)))
        allocate(qsc(nspecies_sfire,2,nnxp(ifm),nnyp(ifm)))

        allocate(b_fires_area(max_sfire_points_by_brams,nnxp(ifm),nnyp(ifm)))
        allocate(b_fires_time(max_sfire_points_by_brams,nnxp(ifm),nnyp(ifm)))
        allocate(b_fires_frp(max_sfire_points_by_brams,nnxp(ifm),nnyp(ifm)))
        allocate(b_fires_lat(max_sfire_points_by_brams,nnxp(ifm),nnyp(ifm)))
        allocate(b_fires_lon(max_sfire_points_by_brams,nnxp(ifm),nnyp(ifm)))

        flam_frac = 0.0
        mean_frp  = 0.0
        std_frp   = 0.0
        mean_size = 0.0
        std_size  = 0.0
        count_in_grid = 0
        valid_ij = 0
        do iesp = 1,nspecies_sfire
            qsc(iesp,:,:,:) = 0.0
        end do        
        allocate(gridVolume(nnxp(1), nnyp(1)))
        do j = 1, nnyp(ifm)
            do i = 1, nnxp(ifm)
                gridVolume(i,j) = (deltax * deltay * zmn(2,1) * rtgt(i,j))/razaoBSF
            end do
        end do

        inquire(file='fireVeg.nml', exist=ex)
        if (.not. ex) then
            print *, "Error: file fireVeg.nml don't exist!"
            print *, "Please, verify!"
            stop
        end if
        open (newunit=iUnit, FILE='fireVeg.nml', STATUS='OLD')
        read (iunit, iostat=err, NML=VEG)
        if (err /= 0) then
            print *, "Error: reading fireVeg.nml!"
            print *, "Please, check sintax!"
            stop   
        else
            print *,"Modis Veg File is ", vegFile, " ", test_only, " ", time_2_test
        end if
        

    end subroutine setup_fireVeg

    subroutine get_emission_in_global_brams_grid(glat,glon,dn0)
        implicit none

        real, intent(in) :: dn0(:,:,:)
        real, intent(in) :: glat(:,:)
        real, intent(in) :: glon(:,:)

        integer :: nv, ip, jp, iesp
        integer :: iveg(max_sfire_points_by_brams,nnxp(ifm),nnyp(ifm)), iveg_agreg
        integer, parameter :: itpm = 2 !PM10 old TPM

        real :: bas_by_fire, biomass_burned, qflam, emission_rate
        real :: desvioQ_frp, desvioQ_siz
        integer :: x ,y
        integer :: i,j,pcx,pcy,pt
        integer :: nfocos
        integer, dimension(2) :: pos


        nfocos = size(sfire_info_lat)

        flam_frac = 0.0
        mean_frp  = 0.0
        std_frp   = 0.0
        mean_size = 0.0
        std_size  = 0.0
        valid_ij = 0
        count_in_grid = 0
        do iesp = 1,nspecies_sfire
            qsc(iesp,:,:,:) = 0.0
        end do


        do i = 1, nfocos
            !Localiza o ponto na grade do BRAMS
            pos = localizar_ponto(sfire_info_lat(i), sfire_info_lon(i), glat(1,:),glon(:,1),nnyp(ifm),nnxp(ifm))
            x = pos(2)
            y = pos(1)
            ! Incrementa a contagem de pontos no ponto de grade
            count_in_grid(x,y) = count_in_grid(x,y) + 1
            valid_ij(x,y) = count_in_grid(x,y)
            !Copia os valores de fire obtidos para a grade do BRAMS
            b_fires_area(count_in_grid(x,y),x,y) = sfire_info_area(i)
            b_fires_frp (count_in_grid(x,y),x,y) = sfire_info_frp(i)
            b_fires_time(count_in_grid(x,y),x,y) = sfire_info_time(i)
            b_fires_lat (count_in_grid(x,y),x,y) = sfire_info_lat(i)
            b_fires_lon (count_in_grid(x,y),x,y) = sfire_info_lon(i)
        end do

        do j = 1,nnyp(ifm)
            do i = 1,nnxp(ifm)
                if(count_in_grid(i,j) == 0) cycle !Sem pontos sfire na grade
                do pt = 1, count_in_grid(i,j)
                    !Localiza o ponto na grade de 1km para pegar a vegetação
                    pos = localizar_ponto(b_fires_lat(pt,i,j), b_fires_lon(pt,i,j), lat, lon, npy, npx)
                    !if (i==32 .and. j==134) print *,'pos_1: ',b_fires_frp(pt,i,j)
                    if (b_fires_frp(pt,i,j)<=0) then !Sem queima
                        valid_ij(i,j) = valid_ij(i,j)-1
                        cycle 
                    end if
                    iveg(pt,i,j) = vegType(1,1,pos(2),pos(1))
                    if (iveg(pt,i,j)==0 .or. iveg(pt,i,j)==15 .or. iveg(pt,i,j)>16) then ! Não é vegetação
                        valid_ij(i,j) = valid_ij(i,j)-1
                        cycle
                    end if
                    iveg_agreg= catb(iveg(pt,i,j))
                    biomass_burned =  b_fires_frp(pt,i,j) ! kg[dry biomass burned]
                    qflam =  flaming(iveg_agreg)
                    biomass_burned =  fx*biomass_burned
                    do iesp = 1,nspecies_sfire
                        emission_rate = biomass_burned * (emission_factor(iveg_agreg,iesp) &
                                        /emission_factor(iveg_agreg,itpm))
                        emission_rate = emission_rate/b_fires_time(pt,i,j)/dn0(2,i,j)/gridVolume(i,j)
                        qsc(iesp,2,i,j) =  qsc(iesp,2,i,j) + emission_rate * qflam !ACrescentado a fração Level 2
                        qsc(iesp,1,i,j) =  qsc(iesp,1,i,j) + emission_rate * (1-qflam) !ACrescentado a fração  Level 1
                    end do
                    !Acumulando os valores para cálculo da média
                    !flam_frac(i,j) = 1.0 !flam_frac(i,j) + qflam !Fixado em 1.0 para evitar o cálculo em plumerise
                    mean_frp(i,j)  = mean_frp(i,j)  + b_fires_frp(pt,i,j)
                    !if (i==32 .and. j==134) print *,'pos_2: ',mean_frp(i,j)
                    mean_size(i,j) = mean_size(i,j) + b_fires_area(pt,i,j)
                end do
                ! Fazendo a média
                mean_frp(i,j)  = mean_frp(i,j) /valid_ij(i,j)
                !if (i==32 .and. j==134) print *,'pos_3: ',mean_frp(i,j)
                mean_size(i,j) = mean_size(i,j)/valid_ij(i,j)
                ! Precisa verificar se flam_frac é média ou não!
                flam_frac(i,j) = 1.0 !flam_frac(i,j)/valid_ij(i,j)
            end do
        end do

        !Calculando o desvio padrão
        do j = 1,nnyp(ifm)
            do i = 1,nnxp(ifm)
                desvioQ_frp = 0.0 
                desvioQ_siz = 0.0
                if(valid_ij(i,j) <= 0) cycle !Sem pontos válidos sfire na grade
                !Soma os desvios quadráticos
                do pt = 1, count_in_grid(i,j)
                    if (b_fires_frp(pt,i,j)<=0.0) cycle
                    if (iveg(pt,i,j)==0 .or. iveg(pt,i,j)==15 .or. iveg(pt,i,j)>16) cycle
                    desvioQ_frp = desvioQ_frp + (b_fires_frp(pt,i,j)  - mean_frp(i,j))**2
                    !if (i==32 .and. j==134) print *,'pos_4: ',mean_frp(i,j),b_fires_frp(pt,i,j),desvioQ_frp
                    desvioQ_siz = desvioQ_siz + (b_fires_area(pt,i,j) - mean_size(i,j))**2
                end do
                !Calcula o desvio padrão
                std_frp(i,j)  = sqrt(desvioQ_frp/valid_ij(i,j))
                !if (i==32 .and. j==134) print *,'pos_5: ',std_frp(i,j)
                std_size(i,j) = sqrt(desvioQ_siz/valid_ij(i,j))
            end do
        end do
            
    end subroutine get_emission_in_global_brams_grid


    function localizar_ponto(lat, lon, lat_array, lon_array, nx, ny) result(pos)

        implicit none
        real, intent(in) :: lat, lon                  ! Latitude e longitude do ponto a ser localizado
        real, intent(in) :: lat_array(nx), lon_array(ny) ! Arrays de latitudes e longitudes
        integer, intent(in) :: nx, ny                     ! Dimensões dos arrays
        integer :: i, j
        real(8) :: dist_min, dist
        integer, dimension(2) :: pos                      ! Posição relativa do ponto (x, y)
        
        ! Inicializando valores
        dist_min = 1.0e30   ! Um valor muito alto de distância para iniciar
        pos = (-1, -1)      ! Valores default caso não sejam encontrados
    
        ! Laço para localizar a posição mais próxima do ponto (lat, lon)
        do i = 1, nx
            do j = 1, ny
                ! Calcula a distância quadrada entre o ponto e o array (métrica simples de distância euclidiana 2D)
                dist = (lat - lat_array(i))**2 + (lon - lon_array(j))**2
    
                ! Se a distância atual for menor que a menor distância encontrada até agora
                if (dist < dist_min) then
                    dist_min = dist
                    pos = (/i, j/)  ! Atualiza a posição x, y
                end if
            end do
        end do
    
    end function localizar_ponto

    subroutine comm_aer_data(mzp, mxp, myp, ia, iz, ja, jz)
        use node_mod, only: &
            mynum, &
            mchnum, &
            master_num, &
            nodei0, nodej0, &
            nmachs, & 
            nodemxp, &
            nodemyp

        use mem_aer1, only: aer1_vars, aer1_g

        use mem_plume_chem1, only:  &
               plume_fre_g, &
               iflam_frac, &
               istd_frp, &
               imean_frp, &
               imean_size, &
               istd_size

        use aer1_list, only : spc_alloc,spc_name, src, ddp, wdp, fdda, offline, on ,off &
               ,mode_alloc, mode_name, aer_name, nspecies_aer=> nspecies,nmodes, bburn

        implicit none

        integer, intent(in) :: mzp, mxp, myp, ia, iz, ja, jz

        real :: aer_L1(nspecies_sfire,nnxp(1),nnyp(1))
        real :: aer_L2(nspecies_sfire,nnxp(1),nnyp(1))
        real :: aerLocal(nspecies_sfire,mxp,myp)
        real :: pvar(mxp,myp)
        integer :: i1,i2,j1,j2

        real :: flam_frac_com(nnxp(1),nnyp(1))
        real :: mean_frp_com(nnxp(1),nnyp(1))
        real :: std_frp_com(nnxp(1),nnyp(1))
        real :: mean_size_com(nnxp(1),nnyp(1))
        real :: std_size_com(nnxp(1),nnyp(1))
        character(len=16) :: register_name

        integer :: k, ispc, imode

        if (mchnum == master_num) then
            aer_L1 = qsc(:,1,:,:)
            aer_L2 = qsc(:,2,:,:)
            flam_frac_com  = flam_frac
            mean_frp_com   = mean_frp
            std_frp_com    = std_frp
            mean_size_com  = mean_size
            std_size_com   = std_size
        end if

        i1 = nodei0(mynum,ifm)+1
        i2 = nodei0(mynum,ifm)+nodemxp(mynum,ifm)
        j1 = nodej0(mynum,ifm)+1
        j2 = nodej0(mynum,ifm)+nodemyp(mynum,ifm)

        !Send frp data to all processors
        call parf_bcast(flam_frac_com, int(nnxp(ifm),i8), int(nnyp(ifm),i8),master_num)
        call parf_bcast(mean_frp_com , int(nnxp(ifm),i8), int(nnyp(ifm),i8),master_num)
        call parf_bcast(std_frp_com  , int(nnxp(ifm),i8), int(nnyp(ifm),i8),master_num)
        call parf_bcast(mean_size_com, int(nnxp(ifm),i8), int(nnyp(ifm),i8),master_num)
        call parf_bcast(std_size_com , int(nnxp(ifm),i8), int(nnyp(ifm),i8),master_num)

        !Send the aer from master to all processors 
        do k=1,nspecies_sfire
            !aer = qsc(k,:,:)
            call parf_bcast(aer_L1(k,:,:), int(nnxp(ifm),i8), int(nnyp(ifm),i8),master_num)
        end do
        !Copy the geographic part of aer from global aer to local array, aerLocal
        aerLocal = 0.0
        call copy_buff(a=aer_L1(1,1,1), b=aerLocal(1,1,1), &
        n1=nspecies_sfire, n2=nnxp(ifm), n3=nnyp(ifm), m1=nspecies_sfire, m2=nodemxp(mynum,ifm), m3=nodemyp(mynum,ifm), &
        i1=i1, i2=i2, j1=j1, j2=j2)

        do ispc=1,nspecies_aer
            do imode=1,nmodes
               register_name=trim(spc_name(ispc))//trim(mode_name(imode))
               if (trim(register_name)=='bburn2') then
                  aer1_g(imode, bburn, 1)%sc_src(1, :,:) = &
                  aer1_g(imode, bburn, 1)%sc_src(1, :,:) + aerLocal(1,:,:)
               end if
               if (trim(register_name)=='bburn3') then
                  aer1_g(imode, bburn, 1)%sc_src(1, :,:) = &
                  aer1_g(imode, bburn, 1)%sc_src(1, :,:) + aerLocal(2,:,:)
               end if
            end do
         end do

        !Send the aer from master to all processors 
        do k=1,nspecies_sfire
            !aer = qsc(k,:,:)
            call parf_bcast(aer_L2(k,:,:), int(nnxp(ifm),i8), int(nnyp(ifm),i8),master_num)
        end do
        !Copy the geographic part of aer from global aer to local array, aerLocal
        aerLocal = 0.0
        call copy_buff(a=aer_L2(1,1,1), b=aerLocal(1,1,1), &
        n1=nspecies_sfire, n2=nnxp(ifm), n3=nnyp(ifm), m1=nspecies_sfire, m2=nodemxp(mynum,ifm), m3=nodemyp(mynum,ifm), &
        i1=i1, i2=i2, j1=j1, j2=j2)

        do ispc=1,nspecies_aer
            do imode=1,nmodes
               register_name=trim(spc_name(ispc))//trim(mode_name(imode))
               if (trim(register_name)=='bburn2') then
                  aer1_g(imode, bburn, 1)%sc_src(2, :,:) = &
                  aer1_g(imode, bburn, 1)%sc_src(2, :,:) + aerLocal(1,:,:)
               end if
               if (trim(register_name)=='bburn3') then
                  aer1_g(imode, bburn, 1)%sc_src(2, :,:) = &
                  aer1_g(imode, bburn, 1)%sc_src(2, :,:) + aerLocal(2,:,:)
               end if
            end do
         end do


        call copy_buff2D(a=flam_frac_com, b=pvar, &
        n2=nnxp(ifm), n3=nnyp(ifm), m2=nodemxp(mynum,ifm), m3=nodemyp(mynum,ifm), &
        i1=i1, i2=i2, j1=j1, j2=j2)
        plume_fre_g(iflam_frac,ifm)%pvar = pvar

        call copy_buff2D(a=mean_frp_com, b=pvar, &
        n2=nnxp(ifm), n3=nnyp(ifm), m2=nodemxp(mynum,ifm), m3=nodemyp(mynum,ifm), &
        i1=i1, i2=i2, j1=j1, j2=j2)
        plume_fre_g(imean_frp,ifm)%pvar = pvar 

        call copy_buff2D(a=std_frp_com, b=pvar, &
        n2=nnxp(ifm), n3=nnyp(ifm), m2=nodemxp(mynum,ifm), m3=nodemyp(mynum,ifm), &
        i1=i1, i2=i2, j1=j1, j2=j2)
        plume_fre_g(istd_frp,ifm)%pvar = pvar   
     
        call copy_buff2D(a=mean_size_com, b=pvar, &
        n2=nnxp(ifm), n3=nnyp(ifm), m2=nodemxp(mynum,ifm), m3=nodemyp(mynum,ifm), &
        i1=i1, i2=i2, j1=j1, j2=j2)
        plume_fre_g(imean_size,ifm)%pvar = pvar  

        call copy_buff2D(a=std_size_com, b=pvar, &
        n2=nnxp(ifm), n3=nnyp(ifm), m2=nodemxp(mynum,ifm), m3=nodemyp(mynum,ifm), &
        i1=i1, i2=i2, j1=j1, j2=j2)
        plume_fre_g(istd_size,ifm)%pvar = pvar     



    end subroutine comm_aer_data

    subroutine copy_buff(a,b,n1,n2,n3,m1,m2,m3,i1,i2,j1,j2)
        implicit none
        integer, intent(in) :: n1,n2,n3
        integer, intent(in) :: m1,m2,m3
        integer, intent(in) :: i1,i2,j1,j2
        real, intent(in)    :: a(n1,n2,n3)
        real, intent(out)   :: b(m1,m2,m3)
      
        b(1:m1,1:m2,1:m3) = a(1:n1,i1:i2,j1:j2)
      
      end subroutine copy_buff
  
      subroutine copy_buff2D(a,b,n2,n3,m2,m3,i1,i2,j1,j2)
        implicit none
        integer, intent(in) :: n2,n3
        integer, intent(in) :: m2,m3
        integer, intent(in) :: i1,i2,j1,j2
        real, intent(in)    :: a(n2,n3)
        real, intent(out)   :: b(m2,m3)
      
        b(1:m2,1:m3) = a(i1:i2,j1:j2)
      
      end subroutine copy_buff2D
    
      subroutine testModSfire2Brams(glat,glon,dn0,time)
        !Simula um caso de focos de fogo ativos
        !Lê uma tabela contendo as informações abaixo:
        !
        !-------------------------------------------------
        !lat    lon      frp[W]      area [m2]   time[s]
        !38
        !40.26  -7.73    0.05e+08    0.340e+06   30
        !40.26  -7.73    0.08e+10    0.440e+06   30
        !...
        !-------------------------------------------------
        !
        ! Na segunda linha está definido o número de focos e depois,
        ! nas linhas seguintes, segue a sequência de informações sobre cada foco.
        !
        !Preenche os arrays e chama a rotina que gera as emissões e
        !estatísticas

        implicit none

        real, intent(in) :: dn0(:,:,:)
        real, intent(in) :: glat(:,:)
        real, intent(in) :: glon(:,:)
        real, intent(in) :: time

        character(len=256) :: lixo
        integer :: nfocos
        integer :: k
        logical :: ok

        if (time /= time_2_test) return 

        !Le uma tabela de focos de fogo simulando uma condição real
        open(unit=22,file="tabela_focos.txt",status="old",action="read")
        read(22,*) lixo
        read(22,*) nfocos

        !Aloca os arrays a serem usados
        ok = alloc_dealloc_sfire_info(nfocos)
        !Preenche os arrays com os focos lidos do arquivo
        do k=1,nfocos
           read(22,*) sfire_info_lat (k), sfire_info_lon (k), sfire_info_frp (k), &
              sfire_info_area(k), sfire_info_time(k) 
        end do
        close(22)

        !Chama a rotina que calcula as emissões segundo os dados vindos de sfire
        !A  rotina calcula as emissões e estatísticas - Deve estar dentro do nó mestre
        call get_emission_in_global_brams_grid(glat,glon,dn0)

        !Dealoca os arrays sem uso
        ok = alloc_dealloc_sfire_info()

    end subroutine testModSfire2Brams

    subroutine fill_sfire_info_and_get_emissions(ifds,ifde,jfds,jfde,frp,fire_area,latf,lonf,tburn,glat,glon,dn0)
        implicit none

        integer, intent(in) :: ifde,ifds,jfds,jfde
        real, intent(in) :: frp(:,:)
        real, intent(in) :: fire_area(:,:)
        real, intent(in) :: latf(:,:)
        real, intent(in) :: lonf(:,:)
        real, intent(in) :: tburn(:,:)
        real, intent(in) :: dn0(:,:,:)
        real, intent(in) :: glat(:,:)
        real, intent(in) :: glon(:,:)
        
        integer :: i, j, ifoc, isize, jsize, k
        integer, allocatable :: ival(:), jval(:)
        logical :: ok

        isize = size(latf,1)
        jsize = size(latf,2)
        allocate(ival(isize),jval(jsize))

        !Determinando focos válidos
        ifoc = 0
        do j = jfds, ifde
            do i = ifde ,ifds
                if(frp(i,j)>0.0 .and. fire_area(i,j)>0.0) then
                    ifoc = ifoc+1
                    ival(ifoc) = i
                    jval(ifoc) = j
                end if
            end do
        end do

        ok = alloc_dealloc_sfire_info(ifoc)
        !Percorrendo focos válidos e preenchendo sfire_info
        do i = k ,ifoc
            sfire_info_lat (k) = latf(ival(k),jval(k))
            sfire_info_lon (k) = lonf(ival(k),jval(k))
            sfire_info_frp (k) = frp(ival(k),jval(k))
            sfire_info_area(k) = fire_area(ival(k),jval(k))
            sfire_info_time(k) = tburn(ival(k),jval(k))
        end do        
        
        !Chama a rotina que calcula as emissões segundo os dados vindos de sfire
        !A  rotina calcula as emissões e estatísticas - Deve estar dentro do nó mestre
        call get_emission_in_global_brams_grid(glat,glon,dn0)

        !Dealoca os arrays sem uso
        ok = alloc_dealloc_sfire_info()

    end subroutine fill_sfire_info_and_get_emissions

end module modSfire2Brams