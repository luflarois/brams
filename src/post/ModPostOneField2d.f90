module ModPostOneField2d

   USE mem_grid, ONLY : time   ! INTENT(IN)  !DSM
   use io_params, only: frqanl
   use micphys, only : &
         mcphys_type, &! INTENT(IN)
         level          ! INTENT(IN)
   !
   !use mem_grid  , only: ngrid ! INTENT(IN)

   use ModOutputUtils, only : GetVarFromMemToOutput

   use ModPostGrid, only : OutputGradsField
   use ModBramsGrid, only : BramsGrid
   use ModPostTypes, only : PostGrid
   use ModPostOneFieldUtils, only : PrepareAndOutputGradsField
   use ModPostOneFieldUtils, only : PrepareGradsField
   !use ModPostOneFieldUtils, only : PostVarType
   use ModPostTypes, only: PostVarType

   use ModPostUtils, only : undef
   use ModPostUtils, only : rams_comp_tempc
   use ModPostUtils, only : rams_comp_tempk
   use ModPostUtils, only : rams_comp_dewk
   use ModPostUtils, only : rams_comp_thetv
   use ModPostUtils, only : rams_get_surface
   use ModPostUtils, only : rams_comp_1minus
   use ModPostUtils, only : rams_comp_slpmm5
   use ModPostUtils, only : rams_comp_rh
   use ModPostUtils, only : rams_comp_press
   use ModPostUtils, only : rams_comp_pbl
   use ModPostUtils, only : cape_cine
   use ModPostUtils, only : rams_comp_dn0
   use ModPostUtils, only : get_ZItheta
   use ModPostUtils, only : rams_comp_speed
   use ModPostUtils, only : rams_comp_sfc_press
   use ModPostUtils, only : rams_reduced_temp
   use ModPostUtils, only : rams_fill_sst
   use ModPostUtils, only : get_leaf_soil
   use ModPostUtils, only : comp_vertint
   use ModPostUtils, only : comp_vertint_press
   use ModPostUtils, only : comp_slp_metar
   use ModPostUtils, only : rams_reduced_rv
   use ModPostUtils, only : rams_comp_dewk_2m
   use ModPostUtils, only : rams_reduced_wind
   use ModPostUtils, only : rams_comp_dir
   use ModPostUtils, only : calc_u10m
   use ModPostUtils, only : calc_v10m
   use ModPostUtils, only : relative_humidity_2m
   use ModPostUtils, only : calc_poda_index
   use ModPostUtils, only : copy_x_to_y
   use ModPostUtils, only : checkUsingJules

   !LFR
   use modTimeLineFRN, only: writeTimeLineFRN

   use io_params, only : & ! 
      IPOS

   implicit none
   real, allocatable, dimension(:,:) :: convprec,totprec
   logical :: alloc2d = .true. 
   private

   public :: Brams2Post_2d

contains


   subroutine Brams2Post_2d (one_post_variable, oneBramsGrid, onePostGrid)
      use dump, only: &
        dumpMessage
      use mem_radiate, only: &
        ISWRTYP, &
        ILWRTYP
      use node_mod, only:  &
             mynum, &
             master_num
      use mem_aerad, only: nwave  !INTENT(IN)

      include "constants.f90"
      type(PostVarType) :: one_post_variable
      type(BramsGrid), pointer :: oneBramsGrid
      type(PostGrid), pointer :: onePostGrid

      real :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp)
      real :: OutputField3d(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp) ! for use in a special case
      real :: ScrT1N01
      real :: ScrT1N02

      real :: ScrT2N01(oneBramsGrid%mxp, oneBramsGrid%myp)
      real :: ScrT2N02(oneBramsGrid%mxp, oneBramsGrid%myp)
      real :: ScrT2N03(oneBramsGrid%mxp, oneBramsGrid%myp)
      real :: ScrT2N04(oneBramsGrid%mxp, oneBramsGrid%myp)
      real :: ScrT2N05(oneBramsGrid%mxp, oneBramsGrid%myp)

      real :: ScrT3N01(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)
      real :: ScrT3N02(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)
      real :: ScrT3N03(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)
      real :: ScrT3N04(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)
      real :: ScrT3N05(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)

      real :: ScrT4N01(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%nzg, oneBramsGrid%npatch)

      real :: ScrT6N01(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%npatch)
      real :: ScrT6N02(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%npatch)
      real :: ScrT6N03(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%npatch)
      real :: ScrT6N04(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%npatch)
      real :: ScrT6N05(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%npatch)
      real :: ScrT6N06(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%npatch)
      
      real :: ScrT7N01(oneBramsGrid%mxp, oneBramsGrid%myp, nwave)
      
      character(len=*),parameter :: header='***(Brams2Post_2d)***'
      character(len=*),parameter :: version='5.4'

      integer :: i, j, k,ierr

      if(alloc2d) then 
	 !if(one_post_variable%fieldName == 'CONVPREC' .or.  &
         !   one_post_variable%fieldName == 'TOTPREC'   ) then 
	    call alloc2d_routine(oneBramsGrid%mxp, oneBramsGrid%myp)
            alloc2d = .false. 
	 !endif
      endif      
      

      select case (one_post_variable%fieldName)
      case ('TOTPCP')
         call getAccComponents(oneBramsGrid, OutputField)
         OutputField = max(OutputField, 0.0)
      case ('ACCCON')
         OutputField = getAconpr(oneBramsGrid)
         OutputField = max(OutputField, 0.0)
      case ('CONVPREC')
         OutputField = getAconpr(oneBramsGrid)-convprec
	 convprec    = getAconpr(oneBramsGrid) ! save for the next step
         !--- output mm/hour
         OutputField = max(OutputField, 0.0)/(frqanl/3600.)

      case ('RSHORT')
         call GetVarFromMemToOutput ('RSHORT', oneBramsGrid%currGrid, OutputField)
      case ('RLONG')
         call GetVarFromMemToOutput ('RLONG', oneBramsGrid%currGrid, OutputField)
      case ('SEA_PRESS')
         call GetVarFromMemToOutput ('TOPT', oneBramsGrid%currGrid, ScrT2N01)
         call GetVarFromMemToOutput ('PI', oneBramsGrid%currGrid, ScrT3N01)
         call GetVarFromMemToOutput ('THETA', oneBramsGrid%currGrid, ScrT3N02)
         call GetVarFromMemToOutput ('RV', oneBramsGrid%currGrid, ScrT3N03)
         call rams_comp_thetv (ScrT3N02, ScrT3N03)
         OutputField = ScrT3N03(:, :, 1)
         call rams_comp_slpmm5 (ScrT3N02, ScrT3N01, ScrT2N01, OutputField)
      case ('CAPE')
         call GetVarFromMemToOutput ('RV', oneBramsGrid%currGrid, ScrT3N01)
         call GetVarFromMemToOutput ('PI', oneBramsGrid%currGrid, ScrT3N02)
         call GetVarFromMemToOutput ('THETA', oneBramsGrid%currGrid, ScrT3N03)
         call rams_comp_rh (ScrT3N01, ScrT3N02, ScrT3N03)
         ScrT3N01 = max(ScrT3N01, 0.0)
         call GetVarFromMemToOutput ('THETA', oneBramsGrid%currGrid, ScrT3N03)
         call GetVarFromMemToOutput ('PI', oneBramsGrid%currGrid, ScrT3N02)
         call rams_comp_tempk (ScrT3N03, ScrT3N02)
         call GetVarFromMemToOutput ('PI', oneBramsGrid%currGrid, ScrT3N02)
         call rams_comp_press(ScrT3N02)
         call cape_cine (ScrT3N02, ScrT3N03, ScrT3N01, OutputField, 'cape')
      case ('CINE')
         call GetVarFromMemToOutput ('RV', oneBramsGrid%currGrid, ScrT3N01)
         call GetVarFromMemToOutput ('PI', oneBramsGrid%currGrid, ScrT3N02)
         call GetVarFromMemToOutput ('THETA', oneBramsGrid%currGrid, ScrT3N03)
         call rams_comp_rh (ScrT3N01, ScrT3N02, ScrT3N03)
         ScrT3N01 = max(ScrT3N01, 0.0)
         call GetVarFromMemToOutput ('THETA', oneBramsGrid%currGrid, ScrT3N03)
         call GetVarFromMemToOutput ('PI', oneBramsGrid%currGrid, ScrT3N02)
         call rams_comp_tempk (ScrT3N03, ScrT3N02)
         call GetVarFromMemToOutput ('PI', oneBramsGrid%currGrid, ScrT3N02)
         call rams_comp_press(ScrT3N02)
         call cape_cine (ScrT3N02, ScrT3N03, ScrT3N01, OutputField, 'cine')
      case ('TOPO')
         call GetVarFromMemToOutput ('TOPT', oneBramsGrid%currGrid, OutputField)
      case ('PRECIP')
         call getAccComponents(oneBramsGrid, OutputField)
         ScrT2N01 = getAconpr(oneBramsGrid)
         OutputField = OutputField + ScrT2N01
         OutputField = max(OutputField, 0.0)
      case ('TOTPREC')
         call getAccComponents(oneBramsGrid, ScrT2N02)
         ScrT2N01 = getAconpr(oneBramsGrid)
         OutputField = ScrT2N02 + ScrT2N01 - TotPrec
	 TotPrec= ScrT2N02 + ScrT2N01  ! save for the next step
	 
         !--- output mm/hour
         OutputField = max(OutputField, 0.0)/(frqanl/3600.)
      case ('LE')
         call GetVarFromMemToOutput ('SFLUX_R', oneBramsGrid%currGrid, OutputField)
         call GetVarFromMemToOutput ('TOPT', oneBramsGrid%currGrid, ScrT2N01)
         call rams_comp_dn0 (ScrT3N01, ScrT3N02, ScrT3N03, ScrT2N01, oneBramsGrid%pi01dn, &
               oneBramsGrid%th01dn, oneBramsGrid%ztn, oneBramsGrid%ztop, oneBramsGrid%dzmn)
         call rams_get_surface(ScrT2N01, ScrT3N03)
         OutputField = OutputField * ScrT2N01
         OutputField = OutputField * 2.5e6
      case ('H')
         call GetVarFromMemToOutput ('SFLUX_T', oneBramsGrid%currGrid, OutputField)
         call GetVarFromMemToOutput ('TOPT', oneBramsGrid%currGrid, ScrT2N01)
         call rams_comp_dn0 (ScrT3N01, ScrT3N02, ScrT3N03, ScrT2N01, oneBramsGrid%pi01dn, &
               oneBramsGrid%th01dn, oneBramsGrid%ztn, oneBramsGrid%ztop, oneBramsGrid%dzmn)
         call rams_get_surface(ScrT2N01, ScrT3N03)
         OutputField = OutputField * ScrT2N01
         OutputField = OutputField * 1004.
      case ('RLONGUP')
         call GetVarFromMemToOutput ('RLONGUP', oneBramsGrid%currGrid, OutputField)
      case ('ALBEDT')
         call GetVarFromMemToOutput ('ALBEDT', oneBramsGrid%currGrid, OutputField)
      case ('TEMPC2M')
         call GetVarFromMemToOutput ('UP', oneBramsGrid%currGrid, ScrT3N01)
         call GetVarFromMemToOutput ('VP', oneBramsGrid%currGrid, ScrT3N02)
         call rams_comp_speed (ScrT3N01, ScrT3N02)
         call GetVarFromMemToOutput ('THETA', oneBramsGrid%currGrid, ScrT3N02)
         call GetVarFromMemToOutput ('PI', oneBramsGrid%currGrid, ScrT3N03)
         call GetVarFromMemToOutput ('TOPT', oneBramsGrid%currGrid, ScrT2N01)
         call GetVarFromMemToOutput ('USTAR', oneBramsGrid%currGrid, ScrT6N01)
         call GetVarFromMemToOutput ('PATCH_ROUGH', oneBramsGrid%currGrid, ScrT6N02)
         call GetVarFromMemToOutput ('CAN_TEMP', oneBramsGrid%currGrid, ScrT6N03)
         call GetVarFromMemToOutput ('PATCH_AREA', oneBramsGrid%currGrid, ScrT6N04)
         call GetVarFromMemToOutput ('TSTAR', oneBramsGrid%currGrid, ScrT6N05)
         call rams_reduced_temp (OutputField, ScrT3N01, ScrT6N01, ScrT6N05, 2., oneBramsGrid%ztn(2), &
               ScrT6N02, ScrT6N04, ScrT6N03, ScrT3N02, ScrT3N03, ScrT2N01, oneBramsGrid%ztop)
         call rams_get_surface(ScrT2N01, ScrT3N03)
         call rams_comp_tempk (OutputField, ScrT2N01)
         call rams_comp_tempc (OutputField)
      case ('RH2M')
         ! ToDo - refactoring for this block and Brams2Post_td2mj to use the same code
         !from Brams2Post_td2mj
         call checkUsingJules(one_post_variable%fieldName)
         !real :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp)
         call GetVarFromMemToOutput ('RV2MJ', oneBramsGrid%currGrid, ScrT2N02)
         ScrT2N02 = ScrT2N02 * 1.e3
         call GetVarFromMemToOutput ('T2MJ', oneBramsGrid%currGrid, ScrT2N03)
         call GetVarFromMemToOutput ('PI', oneBramsGrid%currGrid, ScrT3N01)
         ScrT2N03 = (ScrT3N01(:, :, 1) + ScrT3N01(:, :, 2)) * 0.5
         call RAMS_comp_dewK_2m(ScrT2N02, ScrT2N03, ScrT2N03)
         do j = 1, size(ScrT2N02, 2)
            do i = 1, size(ScrT2N02, 1)
               ScrT2N02(i, j) = ScrT2N02(i, j) - 273.16
               ScrT2N03(i, j) = ScrT2N03(i, j) - 273.16
            end do
         end do
         call relative_humidity_2m(OutputField, ScrT2N02, ScrT2N03)
         OutputField = max(OutputField, 0.0)
      case ('SST')
         call GetVarFromMemToOutput ('SOIL_ENERGY', oneBramsGrid%currGrid, ScrT4N01)
         call rams_fill_sst (oneBramsGrid%nzg, OutputField, ScrT4N01)
      case ('LAND')
         call GetVarFromMemToOutput ('PATCH_AREA', oneBramsGrid%currGrid, ScrT6N01)
         call rams_comp_1minus(OutputField, ScrT6N01)
      case ('PWT')
         ScrT3N05 = 0.0
         call GetVarFromMemToOutput ('RV', oneBramsGrid%currGrid, ScrT3N04)
         ScrT3N05 = ScrT3N05 + ScrT3N04
         call GetVarFromMemToOutput ('RCP', oneBramsGrid%currGrid, ScrT3N04)
         ScrT3N05 = ScrT3N05 + ScrT3N04
         call GetVarFromMemToOutput ('RRP', oneBramsGrid%currGrid, ScrT3N04)
         ScrT3N05 = ScrT3N05 + ScrT3N04
         call GetVarFromMemToOutput ('RPP', oneBramsGrid%currGrid, ScrT3N04)
         ScrT3N05 = ScrT3N05 + ScrT3N04
         call GetVarFromMemToOutput ('RSP', oneBramsGrid%currGrid, ScrT3N04)
         ScrT3N05 = ScrT3N05 + ScrT3N04
         !-For GThompson microphysics micphys_type>1 : RAP does not exist
         if(mcphys_type .le. 1) then
            call GetVarFromMemToOutput ('RAP', oneBramsGrid%currGrid, ScrT3N04)
            ScrT3N05 = ScrT3N05 + ScrT3N04
            call GetVarFromMemToOutput ('RHP', oneBramsGrid%currGrid, ScrT3N04)
            ScrT3N05 = ScrT3N05 + ScrT3N04
         endif
         call GetVarFromMemToOutput ('RGP', oneBramsGrid%currGrid, ScrT3N04)
         ScrT3N05 = ScrT3N05 + ScrT3N04
         ScrT3N05 = max(ScrT3N05, 0.0)
         call GetVarFromMemToOutput ('TOPT', oneBramsGrid%currGrid, ScrT2N01)
         call rams_comp_dn0 (ScrT3N01, ScrT3N02, ScrT3N03, ScrT2N01, oneBramsGrid%pi01dn, &
               oneBramsGrid%th01dn, oneBramsGrid%ztn, oneBramsGrid%ztop, oneBramsGrid%dzmn)
         ScrT3N05 = ScrT3N05 * ScrT3N03
         call comp_vertint (ScrT3N05, ScrT3N05, ScrT2N01, oneBramsGrid%ztop, oneBramsGrid%zmn)
         call rams_get_surface(OutputField, ScrT3N05)
         OutputField = OutputField * 0.1
      case ('SLP_METAR')
         call GetVarFromMemToOutput ('TOPT', oneBramsGrid%currGrid, ScrT2N01)
         call GetVarFromMemToOutput ('PI', oneBramsGrid%currGrid, ScrT3N01)
         call RAMS_comp_press(ScrT3N01)
         call GetVarFromMemToOutput ('THETA', oneBramsGrid%currGrid, ScrT3N02)
         call rams_comp_tempk(ScrT3N02, ScrT3N01)
         call comp_slp_metar(OutputField, ScrT3N01, ScrT2N01, ScrT3N02, oneBramsGrid%ztn)
      case ('TD2M')
         call GetVarFromMemToOutput ('UP', oneBramsGrid%currGrid, ScrT3N01)
         call GetVarFromMemToOutput ('VP', oneBramsGrid%currGrid, ScrT3N02)
         call rams_comp_speed (ScrT3N01, ScrT3N02)
         call GetVarFromMemToOutput ('RV', oneBramsGrid%currGrid, ScrT3N02)
         call GetVarFromMemToOutput ('PI', oneBramsGrid%currGrid, ScrT3N03)
         call GetVarFromMemToOutput ('THETA', oneBramsGrid%currGrid, ScrT3N04)
         call GetVarFromMemToOutput ('TOPT', oneBramsGrid%currGrid, ScrT2N01)
         call GetVarFromMemToOutput ('USTAR', oneBramsGrid%currGrid, ScrT6N01)
         call GetVarFromMemToOutput ('PATCH_ROUGH', oneBramsGrid%currGrid, ScrT6N02)
         call GetVarFromMemToOutput ('CAN_RVAP', oneBramsGrid%currGrid, ScrT6N03)
         call GetVarFromMemToOutput ('PATCH_AREA', oneBramsGrid%currGrid, ScrT6N04)
         call GetVarFromMemToOutput ('RSTAR', oneBramsGrid%currGrid, ScrT6N05)
         call GetVarFromMemToOutput ('CAN_TEMP', oneBramsGrid%currGrid, ScrT6N06)
         call rams_reduced_rv (OutputField, ScrT3N01, ScrT6N01, ScrT6N05, 2., oneBramsGrid%ztn(2), &
               ScrT6N02, ScrT6N04, ScrT6N03, ScrT3N02, ScrT3N03, ScrT2N01, &
               oneBramsGrid%ztop, ScrT6N06, ScrT3N04)
         OutputField = OutputField * 1.e3
         OutputField = max(OutputField, 0.0)
         call GetVarFromMemToOutput ('THETA', oneBramsGrid%currGrid, ScrT3N02)
         call GetVarFromMemToOutput ('CAN_TEMP', oneBramsGrid%currGrid, ScrT6N03)
         call GetVarFromMemToOutput ('TSTAR', oneBramsGrid%currGrid, ScrT6N05)
         call rams_reduced_temp (ScrT2N02, ScrT3N01, ScrT6N01, ScrT6N05, 2., oneBramsGrid%ztn(2), &
               ScrT6N02, ScrT6N04, ScrT6N03, ScrT3N02, ScrT3N03, ScrT2N01, &
               oneBramsGrid%ztop)
         ScrT2N03 = (ScrT3N03(:, :, 1) + ScrT3N03(:, :, 2)) * 0.5
         call rams_comp_tempk (ScrT2N02, ScrT2N03)
         call RAMS_comp_dewK_2m (OutputField, ScrT2N03, ScrT2N02)
         call RAMS_comp_tempC (OutputField)
      case ('U10M')
         call getWinds10m(oneBramsGrid, 'U', OutputField)
      case ('V10M')
         call getWinds10m(oneBramsGrid, 'V', OutputField)
      case ('U10MJ')
         call checkUsingJules(one_post_variable%fieldName)
         if (time==0.) then
            call getWinds10m(oneBramsGrid, 'U', OutputField)
         else
            call GetVarFromMemToOutput ('U10MJ', oneBramsGrid%currGrid, OutputField)
         endif
      case ('V10MJ')
         call checkUsingJules(one_post_variable%fieldName)
         if (time==0.) then
            call getWinds10m(oneBramsGrid, 'V', OutputField)
         else
            call GetVarFromMemToOutput ('V10MJ', oneBramsGrid%currGrid, OutputField)
         endif
      case ('U10MJ1HR')
         call checkUsingJules(one_post_variable%fieldName)
         onePostGrid%ivar_type = one_post_variable%ivar_type
         onePostGrid%fieldUnits = one_post_variable%fieldUnits

         call GetVarFromMemToOutput ('U10MJ1hr', oneBramsGrid%currGrid, ScrT3N01)
         ! output binary field
         onePostGrid%fieldName = 'U10M_6HR_A'
         onePostGrid%fieldDescription = 'Zonal Wind at 10m - 6 hours ago'
         OutputField(:, :) = ScrT3N01(:, :, 1)
         call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
         onePostGrid%fieldName = 'U10M_5HR_A'
         onePostGrid%fieldDescription = 'Zonal Wind at 10m - 5 hours ago'
         OutputField(:, :) = ScrT3N01(:, :, 2)
         call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
         onePostGrid%fieldName = 'U10M_4HR_A'
         onePostGrid%fieldDescription = 'Zonal Wind at 10m - 4 hours ago'
         OutputField(:, :) = ScrT3N01(:, :, 3)
         call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
         onePostGrid%fieldName = 'U10M_3HR_A'
         onePostGrid%fieldDescription = 'Zonal Wind at 10m - 3 hours ago'
         OutputField(:, :) = ScrT3N01(:, :, 4)
         call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
         onePostGrid%fieldName = 'U10M_2HR_A'
         onePostGrid%fieldDescription = 'Zonal Wind at 10m - 2 hours ago'
         OutputField(:, :) = ScrT3N01(:, :, 5)
         call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
         onePostGrid%fieldName = 'U10M_1HR_A'
         onePostGrid%fieldDescription = 'Zonal Wind at 10m - 1 hours ago'
         OutputField(:, :) = ScrT3N01(:, :, 6)
         call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
         !return
      case ('V10MJ1HR')
         call checkUsingJules(one_post_variable%fieldName)
         onePostGrid%ivar_type = one_post_variable%ivar_type
         onePostGrid%fieldUnits = one_post_variable%fieldUnits

         onePostGrid%fieldName = 'V10M_6HR_A'
         onePostGrid%fieldDescription = 'Meridional Wind at 10m - 6 hours ago'
         OutputField(:, :) = ScrT3N01(:, :, 1)
         call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)

         onePostGrid%fieldName = 'V10M_5HR_A'
         onePostGrid%fieldDescription = 'Meridional Wind at 10m - 5 hours ago'
         OutputField(:, :) = ScrT3N01(:, :, 2)
         call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)

         onePostGrid%fieldName = 'V10M_4HR_A'
         onePostGrid%fieldDescription = 'Meridional Wind at 10m - 4 hours ago'
         OutputField(:, :) = ScrT3N01(:, :, 3)
         call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)

         onePostGrid%fieldName = 'V10M_3HR_A'
         onePostGrid%fieldDescription = 'Meridional Wind at 10m - 3 hours ago'
         OutputField(:, :) = ScrT3N01(:, :, 4)
         call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)

         onePostGrid%fieldName = 'V10M_2HR_A'
         onePostGrid%fieldDescription = 'Meridional Wind at 10m - 2 hours ago'
         OutputField(:, :) = ScrT3N01(:, :, 5)
         call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)

         onePostGrid%fieldName = 'V10M_1HR_A'
         onePostGrid%fieldDescription = 'Meridional Wind at 10m - 1 hours ago'
         OutputField(:, :) = ScrT3N01(:, :, 6)
         call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
         !return
      case ('T2MJ')
         call checkUsingJules(one_post_variable%fieldName)
         onePostGrid%fieldName = one_post_variable%fieldName
         onePostGrid%ivar_type = one_post_variable%ivar_type
         onePostGrid%fieldUnits = one_post_variable%fieldUnits
         if (time==0.) then
            call GetVarFromMemToOutput ('UP', oneBramsGrid%currGrid, ScrT3N01)
            call GetVarFromMemToOutput ('VP', oneBramsGrid%currGrid, ScrT3N02)
            call rams_comp_speed (ScrT3N01, ScrT3N02)
            call GetVarFromMemToOutput ('THETA', oneBramsGrid%currGrid, ScrT3N02)
            call GetVarFromMemToOutput ('PI', oneBramsGrid%currGrid, ScrT3N03)
            call GetVarFromMemToOutput ('TOPT', oneBramsGrid%currGrid, ScrT2N01)
            call GetVarFromMemToOutput ('USTAR', oneBramsGrid%currGrid, ScrT6N01)
            call GetVarFromMemToOutput ('PATCH_ROUGH', oneBramsGrid%currGrid, ScrT6N02)
            call GetVarFromMemToOutput ('CAN_TEMP', oneBramsGrid%currGrid, ScrT6N03)
            call GetVarFromMemToOutput ('PATCH_AREA', oneBramsGrid%currGrid, ScrT6N04)
            call GetVarFromMemToOutput ('TSTAR', oneBramsGrid%currGrid, ScrT6N05)
            call rams_reduced_temp (OutputField, ScrT3N01, ScrT6N01, ScrT6N05, 2., oneBramsGrid%ztn(2), &
                  ScrT6N02, ScrT6N04, ScrT6N03, ScrT3N02, ScrT3N03, ScrT2N01, oneBramsGrid%ztop)
            call rams_get_surface(ScrT2N01, ScrT3N03)
            call rams_comp_tempk (OutputField, ScrT2N01)
            call rams_comp_tempc (OutputField)
            onePostGrid%fieldDescription = 'temp - 2m AGL;'
         else
            call GetVarFromMemToOutput ('T2MJ', oneBramsGrid%currGrid, OutputField)
            OutputField = OutputField - 273.16
            where (OutputField<-70.) OutputField = undef
            onePostGrid%fieldDescription = 'Temperature at 2m - from JULES'
         endif
         call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
         !return
      case ('T2MJ_MAX')
         call checkUsingJules(one_post_variable%fieldName)
         onePostGrid%fieldName = one_post_variable%fieldName
         onePostGrid%ivar_type = one_post_variable%ivar_type
         onePostGrid%fieldUnits = one_post_variable%fieldUnits
         if (time==0.) then
            call GetVarFromMemToOutput ('UP', oneBramsGrid%currGrid, ScrT3N01)
            call GetVarFromMemToOutput ('VP', oneBramsGrid%currGrid, ScrT3N02)
            call rams_comp_speed (ScrT3N01, ScrT3N02)
            call GetVarFromMemToOutput ('THETA', oneBramsGrid%currGrid, ScrT3N02)
            call GetVarFromMemToOutput ('PI', oneBramsGrid%currGrid, ScrT3N03)
            call GetVarFromMemToOutput ('TOPT', oneBramsGrid%currGrid, ScrT2N01)
            call GetVarFromMemToOutput ('USTAR', oneBramsGrid%currGrid, ScrT6N01)
            call GetVarFromMemToOutput ('PATCH_ROUGH', oneBramsGrid%currGrid, ScrT6N02)
            call GetVarFromMemToOutput ('CAN_TEMP', oneBramsGrid%currGrid, ScrT6N03)
            call GetVarFromMemToOutput ('PATCH_AREA', oneBramsGrid%currGrid, ScrT6N04)
            call GetVarFromMemToOutput ('TSTAR', oneBramsGrid%currGrid, ScrT6N05)
            call rams_reduced_temp (OutputField, ScrT3N01, ScrT6N01, ScrT6N05, 2., oneBramsGrid%ztn(2), &
                  ScrT6N02, ScrT6N04, ScrT6N03, ScrT3N02, ScrT3N03, ScrT2N01, oneBramsGrid%ztop)
            call rams_get_surface(ScrT2N01, ScrT3N03)
            call rams_comp_tempk (OutputField, ScrT2N01)
            call rams_comp_tempc (OutputField)
            onePostGrid%fieldDescription = 'temp - 2m AGL;'
         else
            call GetVarFromMemToOutput ('T2MJ_MAX', oneBramsGrid%currGrid, OutputField)
            OutputField = OutputField - 273.16
            where (OutputField<-70.) OutputField = undef
            onePostGrid%fieldDescription = 'Max Temp at 2m - from JULES'
         endif
         call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
         !return
      case ('T2MJ_MIN')
         call checkUsingJules(one_post_variable%fieldName)
         onePostGrid%fieldName = one_post_variable%fieldName
         onePostGrid%ivar_type = one_post_variable%ivar_type
         onePostGrid%fieldUnits = one_post_variable%fieldUnits
         if (time==0.) then
            call GetVarFromMemToOutput ('UP', oneBramsGrid%currGrid, ScrT3N01)
            call GetVarFromMemToOutput ('VP', oneBramsGrid%currGrid, ScrT3N02)
            call rams_comp_speed (ScrT3N01, ScrT3N02)
            call GetVarFromMemToOutput ('THETA', oneBramsGrid%currGrid, ScrT3N02)
            call GetVarFromMemToOutput ('PI', oneBramsGrid%currGrid, ScrT3N03)
            call GetVarFromMemToOutput ('TOPT', oneBramsGrid%currGrid, ScrT2N01)
            call GetVarFromMemToOutput ('USTAR', oneBramsGrid%currGrid, ScrT6N01)
            call GetVarFromMemToOutput ('PATCH_ROUGH', oneBramsGrid%currGrid, ScrT6N02)
            call GetVarFromMemToOutput ('CAN_TEMP', oneBramsGrid%currGrid, ScrT6N03)
            call GetVarFromMemToOutput ('PATCH_AREA', oneBramsGrid%currGrid, ScrT6N04)
            call GetVarFromMemToOutput ('TSTAR', oneBramsGrid%currGrid, ScrT6N05)
            call rams_reduced_temp (OutputField, ScrT3N01, ScrT6N01, ScrT6N05, 2., oneBramsGrid%ztn(2), &
                  ScrT6N02, ScrT6N04, ScrT6N03, ScrT3N02, ScrT3N03, ScrT2N01, oneBramsGrid%ztop)
            call rams_get_surface(ScrT2N01, ScrT3N03)
            call rams_comp_tempk (OutputField, ScrT2N01)
            call rams_comp_tempc (OutputField)
            onePostGrid%fieldDescription = 'temp - 2m AGL;'
         else
            call GetVarFromMemToOutput ('T2MJ_MIN', oneBramsGrid%currGrid, OutputField)
            OutputField = OutputField - 273.16
            where (OutputField<-70.) OutputField = undef
            onePostGrid%fieldDescription = 'Min Temp at 2m - from JULES'
         endif
         call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
         !return
      case ('ZI')
         ! because is 3d in a 2d variable definition, need to set onePostGrid and call OutputGradsField
         call PrepareGradsField (one_post_variable, onePostGrid)
         call GetVarFromMemToOutput ('TKEP', oneBramsGrid%currGrid, OutputField3d)
         OutputField3d = max(OutputField3d, 0.0)
         call GetVarFromMemToOutput ('TOPT', oneBramsGrid%currGrid, ScrT2N01)
         call rams_comp_pbl (OutputField3d, ScrT2N01, oneBramsGrid%ztn, oneBramsGrid%ztop)
         call OutputGradsField (oneBramsGrid, onePostGrid, OutputField3d)
         !return
      case ('HT_FLUXJ')
         call GetVarFromMemToOutput ('ht_fluxj', oneBramsGrid%currGrid, OutputField)
      case ('CSJ')
         call checkUsingJules(one_post_variable%fieldName)
         onePostGrid%fieldName = one_post_variable%fieldName
         onePostGrid%ivar_type = one_post_variable%ivar_type
         if (time==0.) then
            call GetVarFromMemToOutput ('UP', oneBramsGrid%currGrid, ScrT3N01)
            call GetVarFromMemToOutput ('VP', oneBramsGrid%currGrid, ScrT3N02)
            call rams_comp_speed (ScrT3N01, ScrT3N02)
            call GetVarFromMemToOutput ('THETA', oneBramsGrid%currGrid, ScrT3N02)
            call GetVarFromMemToOutput ('PI', oneBramsGrid%currGrid, ScrT3N03)
            call GetVarFromMemToOutput ('TOPT', oneBramsGrid%currGrid, ScrT2N01)
            call GetVarFromMemToOutput ('USTAR', oneBramsGrid%currGrid, ScrT6N01)
            call GetVarFromMemToOutput ('PATCH_ROUGH', oneBramsGrid%currGrid, ScrT6N02)
            call GetVarFromMemToOutput ('CAN_TEMP', oneBramsGrid%currGrid, ScrT6N03)
            call GetVarFromMemToOutput ('PATCH_AREA', oneBramsGrid%currGrid, ScrT6N04)
            call GetVarFromMemToOutput ('TSTAR', oneBramsGrid%currGrid, ScrT6N05)
            call rams_reduced_temp (OutputField, ScrT3N01, ScrT6N01, ScrT6N05, 2., oneBramsGrid%ztn(2), &
                  ScrT6N02, ScrT6N04, ScrT6N03, ScrT3N02, ScrT3N03, ScrT2N01, oneBramsGrid%ztop)
            call rams_get_surface(ScrT2N01, ScrT3N03)
            call rams_comp_tempk (OutputField, ScrT2N01)
            call rams_comp_tempc (OutputField)
            onePostGrid%fieldDescription = 'temp - 2m AGL;'
            onePostGrid%fieldUnits = 'C'
         else
            call GetVarFromMemToOutput ('CSJ', oneBramsGrid%currGrid, OutputField)
            onePostGrid%fieldDescription = 'Soil Carbon - from JULES'
            onePostGrid%fieldUnits = 'kg/m2'
         endif
         call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
         !return
      case ('RV2MJ')
         call checkUsingJules(one_post_variable%fieldName)
         onePostGrid%fieldName = one_post_variable%fieldName
         onePostGrid%fieldUnits = one_post_variable%fieldUnits
         if (time==0.) then
            ! because is 3d in a 2d variable definition, need to set onePostGrid and call OutputGradsField

            call GetVarFromMemToOutput ('RV', oneBramsGrid%currGrid, OutputField3d)
            Outputfield(:, :) = OutputField3d(:, :, 2) * 1.e3
            OutputField = max(OutputField, 0.0)
            onePostGrid%ivar_type = 2
            onePostGrid%fieldDescription = 'Mixing rate at 2m'
            call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
         else
            call GetVarFromMemToOutput ('RV2MJ', oneBramsGrid%currGrid, OutputField)
            OutputField = OutputField * 1000.
            onePostGrid%ivar_type = 2
            onePostGrid%fieldDescription = 'Mixing rate at 2m - from JULES'
            call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
         endif
         !return
      case ('ZITHETA')
         call GetVarFromMemToOutput ('THETA', oneBramsGrid%currGrid, ScrT3N01)
         call GetVarFromMemToOutput ('RCP', oneBramsGrid%currGrid, ScrT3N02)
         call get_ZItheta(OutputField, ScrT3N01, ScrT3N02, oneBramsGrid%ztn)
      case ('GPP')
         call GetVarFromMemToOutput ('GPP', oneBramsGrid%currGrid, OutputField)
      case ('RESP_S')
         call GetVarFromMemToOutput ('RESP_S', oneBramsGrid%currGrid, OutputField)
      case ('RESP_P')
         call GetVarFromMemToOutput ('RESP_P', oneBramsGrid%currGrid, OutputField)
      case ('NPP')
         call GetVarFromMemToOutput ('NPP', oneBramsGrid%currGrid, OutputField)
      case ('SFC_PRESS')
         call GetVarFromMemToOutput ('PI', oneBramsGrid%currGrid, ScrT3N01)
         call rams_comp_sfc_press (OutputField, ScrT3N01)
      case ('AOT500')
         call GetVarFromMemToOutput ('AOT', oneBramsGrid%currGrid, ScrT7N01)
         Outputfield(:, :) = ScrT7N01(:, :, 11)
      case ('AOT550')
         call GetVarFromMemToOutput ('AOT', oneBramsGrid%currGrid, ScrT7N01)
         Outputfield(:, :) = ScrT7N01(:, :, 10)
      case ('LMO')
         call GetVarFromMemToOutput ('LMO', oneBramsGrid%currGrid, OutputField)
      case ('PBLHGT')
         call GetVarFromMemToOutput ('PBLHGT', oneBramsGrid%currGrid, OutputField)
      case ('TD2MJ')
         call checkUsingJules(one_post_variable%fieldName)
         onePostGrid%fieldName = one_post_variable%fieldName
         onePostGrid%fieldUnits = one_post_variable%fieldUnits
         onePostGrid%ivar_type = 2
         if (time==0.) then
            call GetVarFromMemToOutput ('UP', oneBramsGrid%currGrid, ScrT3N01)
            call GetVarFromMemToOutput ('VP', oneBramsGrid%currGrid, ScrT3N02)
            call rams_comp_speed (ScrT3N01, ScrT3N02)
            call GetVarFromMemToOutput ('RV', oneBramsGrid%currGrid, ScrT3N02)
            call GetVarFromMemToOutput ('PI', oneBramsGrid%currGrid, ScrT3N03)
            call GetVarFromMemToOutput ('THETA', oneBramsGrid%currGrid, ScrT3N04)
            call GetVarFromMemToOutput ('TOPT', oneBramsGrid%currGrid, ScrT2N01)
            call GetVarFromMemToOutput ('USTAR', oneBramsGrid%currGrid, ScrT6N01)
            call GetVarFromMemToOutput ('PATCH_ROUGH', oneBramsGrid%currGrid, ScrT6N02)
            call GetVarFromMemToOutput ('CAN_RVAP', oneBramsGrid%currGrid, ScrT6N03)
            call GetVarFromMemToOutput ('PATCH_AREA', oneBramsGrid%currGrid, ScrT6N04)
            call GetVarFromMemToOutput ('RSTAR', oneBramsGrid%currGrid, ScrT6N05)
            call GetVarFromMemToOutput ('CAN_TEMP', oneBramsGrid%currGrid, ScrT6N06)
            call rams_reduced_rv (OutputField, ScrT3N01, ScrT6N01, ScrT6N05, 2., oneBramsGrid%ztn(2), &
                  ScrT6N02, ScrT6N04, ScrT6N03, ScrT3N02, ScrT3N03, ScrT2N01, &
                  oneBramsGrid%ztop, ScrT6N06, ScrT3N04)
            OutputField = OutputField * 1.e3
            OutputField = max(OutputField, 0.0)
            call GetVarFromMemToOutput ('THETA', oneBramsGrid%currGrid, ScrT3N02)
            call GetVarFromMemToOutput ('CAN_TEMP', oneBramsGrid%currGrid, ScrT6N03)
            call GetVarFromMemToOutput ('TSTAR', oneBramsGrid%currGrid, ScrT6N05)
            call rams_reduced_temp (ScrT2N02, ScrT3N01, ScrT6N01, ScrT6N05, 2., oneBramsGrid%ztn(2), &
                  ScrT6N02, ScrT6N04, ScrT6N03, ScrT3N02, ScrT3N03, ScrT2N01, &
                  oneBramsGrid%ztop)
            ScrT2N03 = (ScrT3N03(:, :, 1) + ScrT3N03(:, :, 2)) * 0.5
            call rams_comp_tempk (ScrT2N02, ScrT2N03)
            call RAMS_comp_dewK_2m (OutputField, ScrT2N03, ScrT2N02)
            call RAMS_comp_tempC (OutputField)
            onePostGrid%fieldDescription = 'Dewpoint temp in 2m'
         else
            call GetVarFromMemToOutput ('RV2MJ', oneBramsGrid%currGrid, OutputField)
            OutputField = OutputField * 1.e3
            call GetVarFromMemToOutput ('T2MJ', oneBramsGrid%currGrid, ScrT2N01)
            call GetVarFromMemToOutput ('PI', oneBramsGrid%currGrid, ScrT3N01)
            ScrT2N03 = (ScrT3N01(:, :, 1) + ScrT3N01(:, :, 2)) * 0.5
            call RAMS_comp_dewK_2m(OutputField, ScrT2N03, ScrT2N01)
            OutputField = OutputField - 273.16
            onePostGrid%fieldDescription = 'Dewpoint temp at 2m - from JULES'
         end if
         call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
         !return
!LFR New chemical outputs on surface
      case ('CO_SFC')
         call PrepareGradsField (one_post_variable, onePostGrid)
         call GetVarFromMemToOutput ('COP', oneBramsGrid%currGrid, OutputField3d)
         Outputfield(:, :) = OutputField3d(:, :, 2)
         OutputField = max(OutputField, 0.0)
         OutputField = OutputField * (28.96 / 28.)
      case ('O3_SFC')
         call PrepareGradsField (one_post_variable, onePostGrid)
         call GetVarFromMemToOutput ('O3P', oneBramsGrid%currGrid, OutputField3d)
         Outputfield(:, :) = OutputField3d(:, :, 2)
         OutputField = max(OutputField, 0.0)
         OutputField = OutputField * (28.96 / 48.)
      case ('NO_SFC')
         call PrepareGradsField (one_post_variable, onePostGrid)
         call GetVarFromMemToOutput ('NOP', oneBramsGrid%currGrid, OutputField3d)
         Outputfield(:, :) = OutputField3d(:, :, 2)
         OutputField = max(OutputField, 0.0)
         OutputField = OutputField * (28.96 / 30.)
      case ('HNO3_SFC')
         call PrepareGradsField (one_post_variable, onePostGrid)
         call GetVarFromMemToOutput ('HNO3P', oneBramsGrid%currGrid, OutputField3d)
         Outputfield(:, :) = OutputField3d(:, :, 2)
         OutputField = max(OutputField, 0.0)
         OutputField = OutputField * (28.96 / 63.)
      case ('NO2_SFC')
         call PrepareGradsField (one_post_variable, onePostGrid)
         call GetVarFromMemToOutput ('NO2P', oneBramsGrid%currGrid, OutputField3d)
         Outputfield(:, :) = OutputField3d(:, :, 2)
         OutputField = max(OutputField, 0.0)
         OutputField = OutputField * (28.96 / 46.)
      case ('PM25_SFC')
         call PrepareGradsField (one_post_variable, onePostGrid)
         call GetVarFromMemToOutput ('bburn2P', oneBramsGrid%currGrid, OutputField3d)
         call GetVarFromMemToOutput ('urban2P', oneBramsGrid%currGrid, ScrT3N04)
         call GetVarFromMemToOutput ('urban3P', oneBramsGrid%currGrid, ScrT3N05)
         ! call RAMS_comp_noneg(n1,n2,n3,a)
         OutputField = max(OutputField3d(:,:,2) + ScrT3N04(:,:,2) + ScrT3N05(:,:,2), 0.0)
         ! call RAMS_comp_mults(n1,n2,n3,a,1.e-9)  ! converte de 1e-9 kg/kg para 1 kg/kg
         OutputField = OutputField * 1.e-9
         ! ierr= RAMS_getvar('TOPT',idim_type,ngrd,e,b,flnm)
         call GetVarFromMemToOutput ('TOPT', oneBramsGrid%currGrid, ScrT2N01)
         ! call RAMS_comp_dn0(n1,n2,n3,b,c,d,e,ngrd)
         call rams_comp_dn0 (ScrT3N01, ScrT3N02, ScrT3N03, ScrT2N01, oneBramsGrid%pi01dn, &
               oneBramsGrid%th01dn, oneBramsGrid%ztn, oneBramsGrid%ztop, oneBramsGrid%dzmn)
         ! call RAMS_transf_ugm3(n1,n2,n3,a,d)
         OutputField = OutputField * scrT3N03(:,:,2) * 1.e+9
         ! call RAMS_comp_noneg(n1,n2,n3,a)
         !OutputField = max(OutputField, 0.0)
!LFR
      case ('CO2_BURN2D')
         ! because is 3d in a 2d variable definition, need to set onePostGrid and call OutputGradsField
         call PrepareGradsField (one_post_variable, onePostGrid)
         call GetVarFromMemToOutput ('TOPT', oneBramsGrid%currGrid, ScrT2N01)
         call GetVarFromMemToOutput ('CO2_bburn_SRC', oneBramsGrid%currGrid, OutputField3d)
         OutputField3d(:, :, 2) = OutputField3d(:, :, 1) * (1. - ScrT2N01(:, :) / oneBramsGrid%ztop) &
               / oneBramsGrid%dztn(2) / 86400. * 1.e-3  ! convertendo de kg/m3/dia para mg/m2/s
         call OutputGradsField (oneBramsGrid, onePostGrid, OutputField3d)
         !return
      case ('CO2_bioge')
         ! because is 3d in a 2d variable definition, need to set onePostGrid and call OutputGradsField
         call PrepareGradsField (one_post_variable, onePostGrid)
         call GetVarFromMemToOutput ('TOPT', oneBramsGrid%currGrid, ScrT2N01)
         call GetVarFromMemToOutput ('CO2_bioge_SRC', oneBramsGrid%currGrid, OutputField3d)
         OutputField3d(:, :, 2) = OutputField3d(:, :, 1) * (1. - ScrT2N01(:, :) / oneBramsGrid%ztop) &
               / oneBramsGrid%dztn(2) * 1.e-3  ! convertendo de kg/m3/dia para mg/m2/s
         call OutputGradsField (oneBramsGrid, onePostGrid, OutputField3d)
         !return
      case ('CO2_TOTAL')
         call GetVarFromMemToOutput ('CO2P', oneBramsGrid%currGrid, ScrT3N04)
         ! call RAMS_comp_noneg(n1,n2,n3,a)
         ScrT3N04 = max(ScrT3N04, 0.0)
         call GetVarFromMemToOutput ('TOPT', oneBramsGrid%currGrid, ScrT2N01)
         call rams_comp_dn0 (ScrT3N01, ScrT3N02, ScrT3N03, ScrT2N01, oneBramsGrid%pi01dn, &
               oneBramsGrid%th01dn, oneBramsGrid%ztn, oneBramsGrid%ztop, oneBramsGrid%dzmn)
         ScrT3N04 = ScrT3N04 * ScrT3N03
         call GetVarFromMemToOutput ('PI', oneBramsGrid%currGrid, ScrT3N01)
         call rams_comp_press(ScrT3N01)
         !Ate aqui acho que esta OK
         call comp_vertint_press (ScrT2N02, ScrT3N04, ScrT3N01, ScrT2N01, oneBramsGrid%ztop, oneBramsGrid%zmn)
         OutputField = ScrT2N02 * 1e-9 * 1000. * 1. / 44. * 1e-4 * 6.022e23
      case ('PM25WD')
         ! ierr= RAMS_getvar('bburn2WD',idim_type,ngrd,a,b,flnm)
         call GetVarFromMemToOutput ('bburn2WD', oneBramsGrid%currGrid, OutputField)
         ! call RAMS_comp_mults(n1,n2,n3,a,1.e-9)  ! converte de 1e-9 kg/kg para 1 kg/kg
         OutputField = OutputField * 1.e-9
      case ('PMINT')
         ! because is 3d in a 2d variable definition, need to set onePostGrid and call OutputGradsField
         call PrepareGradsField (one_post_variable, onePostGrid)
         ! ierr= RAMS_getvar('bburn2P',idim_type,ngrd,a,b,flnm)
         call GetVarFromMemToOutput ('bburn2P', oneBramsGrid%currGrid, OutputField3d)
         call GetVarFromMemToOutput ('urban2P', oneBramsGrid%currGrid, ScrT3N04)
         call GetVarFromMemToOutput ('urban3P', oneBramsGrid%currGrid, ScrT3N05)
         ! call RAMS_comp_noneg(n1,n2,n3,a)
         !OutputField3d = max(OutputField3d, 0.0)
         OutputField3d = max(OutputField3d + ScrT3N04 + ScrT3N05, 0.0)
         ! ierr= RAMS_getvar('TOPT',idim_type,ngrd,e,b,flnm)
         call GetVarFromMemToOutput ('TOPT', oneBramsGrid%currGrid, ScrT2N01)
         ! call RAMS_comp_dn0(n1,n2,n3,b,c,d,e,ngrd)
         call rams_comp_dn0 (ScrT3N01, ScrT3N02, ScrT3N03, ScrT2N01, oneBramsGrid%pi01dn, &
               oneBramsGrid%th01dn, oneBramsGrid%ztn, oneBramsGrid%ztop, oneBramsGrid%dzmn)
         ! call RAMS_comp_mult(n1,n2,n3,a,d) !Unit: kg[pm25]/m3
         OutputField3d = OutputField3d * ScrT3N03
         !call RAMS_comp_vertint(n1,n2,n3,a,e,ngrd) ! Unit: kg[pm25]/m2
         call comp_vertint (OutputField3d, OutputField3d, ScrT2N01, oneBramsGrid%ztop, oneBramsGrid%zmn)
         !call RAMS_comp_mults(n1,n2,n3,a,1.e+6*1.e-9 )  ! converte de kg/m2 para mg/m2
         OutputField3d = OutputField3d * 1.e+6 * 1.e-9        ! converte de kg/m2 para mg/m2
         call OutputGradsField (oneBramsGrid, onePostGrid, OutputField3d)
         !return
      case ('APRGR')
         call GetVarFromMemToOutput ('accpr_gr', oneBramsGrid%currGrid, OutputField)
         OutputField = max(OutputField, 0.0)
      case ('APRST')
         call GetVarFromMemToOutput ('accpr_st', oneBramsGrid%currGrid, OutputField)
         OutputField = max(OutputField, 0.0)
      case ('APRMC')
         call GetVarFromMemToOutput ('accpr_mc', oneBramsGrid%currGrid, OutputField)
         OutputField = max(OutputField, 0.0)
      case ('APRW')
         call GetVarFromMemToOutput ('accpr_w', oneBramsGrid%currGrid, OutputField)
         OutputField = max(OutputField, 0.0)
      case ('APRAS')
         call GetVarFromMemToOutput ('accpr_as', oneBramsGrid%currGrid, OutputField)
         OutputField = max(OutputField, 0.0)
      case ('PODA')
         call checkUsingJules(one_post_variable%fieldName)
         !--- rv_2.0m ---
         call GetVarFromMemToOutput ('RV2MJ', oneBramsGrid%currGrid, ScrT2N01)
         ScrT2N01 = ScrT2N01 * 1.e3
         !--- tempk2m ---
         call GetVarFromMemToOutput ('T2MJ', oneBramsGrid%currGrid, ScrT2N02)
         !--- PI ---
         call GetVarFromMemToOutput ('PI', oneBramsGrid%currGrid, ScrT3N03)
         ScrT2N03 = (ScrT3N03(:, :, 1) + ScrT3N03(:, :, 2)) * 0.5
         call RAMS_comp_dewK_2m(ScrT2N01, ScrT2N03, ScrT2N02)
         !call undef (n1,n2,n3,a,1,-80.0+273,80.0+273.)
         !cdname='dewpoint temp at 2.0m height'
         !cdunits='K'

         !-compute slp (hPa)
         call GetVarFromMemToOutput ('TOPT', oneBramsGrid%currGrid, ScrT2N04)
         call GetVarFromMemToOutput ('THETA', oneBramsGrid%currGrid, ScrT3N01)
         call GetVarFromMemToOutput ('RV', oneBramsGrid%currGrid, ScrT3N02)

         call RAMS_comp_thetv(ScrT3N01, ScrT3N02)
         ScrT2N05 = ScrT3N02(:, :, 1)
         call RAMS_comp_slpmm5(ScrT3N01, ScrT3N03, ScrT2N04, ScrT2N05)
         !cdname='sea level pressure;'
         !cdunits='hPa'
         call copy_x_to_y(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp, ScrT2N05, ScrT3N02)

         !- compute  air temperature at 1.5m height
         !cdunits='K'
         !call undef (n1,n2,n3,c,1,193.0,353.0)
         call calc_poda_index(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp, ScrT2N02, ScrT3N02, ScrT2N01, OutputField)
         !--(tests for wbiol of UV)---

      case ('P47UM')
         call GetVarFromMemToOutput ('P47UM', oneBramsGrid%currGrid, OutputField)
      case ('EXD14')
         call GetVarFromMemToOutput ('EXD14', oneBramsGrid%currGrid, OutputField)
      case ('DNADV')
         call GetVarFromMemToOutput ('DNADV', oneBramsGrid%currGrid, OutputField)
      case ('SCUPM')
         call GetVarFromMemToOutput ('SCUPM', oneBramsGrid%currGrid, OutputField)
      case ('SCUPH')
         call GetVarFromMemToOutput ('SCUPH', oneBramsGrid%currGrid, OutputField)
      case ('CIEHE')
         call GetVarFromMemToOutput ('CIEHE', oneBramsGrid%currGrid, OutputField)
      case ('UVIND')
         call GetVarFromMemToOutput ('UVIND', oneBramsGrid%currGrid, OutputField)
      case ('ERYTH')
         call GetVarFromMemToOutput ('ERYTH', oneBramsGrid%currGrid, OutputField)
      case ('OCTLV')
         call GetVarFromMemToOutput ('OCTLV', oneBramsGrid%currGrid, OutputField)
      case ('PHYTO')
         call GetVarFromMemToOutput ('PHYTO', oneBramsGrid%currGrid, OutputField)
      case ('PHYPH')
         call GetVarFromMemToOutput ('PHYPH', oneBramsGrid%currGrid, OutputField)
      case ('PHYPR')
         call GetVarFromMemToOutput ('PHYPR', oneBramsGrid%currGrid, OutputField)
      case ('CAPIG')
         call GetVarFromMemToOutput ('CAPIG', oneBramsGrid%currGrid, OutputField)
      case ('PLDCW')
         call GetVarFromMemToOutput ('PLDCW', oneBramsGrid%currGrid, OutputField)
      case ('PLDFL')
         call GetVarFromMemToOutput ('PLDFL', oneBramsGrid%currGrid, OutputField)
      case ('MFUP')
         call GetVarFromMemToOutput ('MFUP', oneBramsGrid%currGrid, OutputField)
      case ('MFDD')
         call GetVarFromMemToOutput ('MFDD', oneBramsGrid%currGrid, OutputField)
      case ('MFSH')
         call GetVarFromMemToOutput ('MFSH', oneBramsGrid%currGrid, OutputField)

      case ('AA0')
         call GetVarFromMemToOutput ('accpr_gr', oneBramsGrid%currGrid, OutputField)
      case ('AA1')
         call GetVarFromMemToOutput ('accpr_w', oneBramsGrid%currGrid, OutputField)
      case ('AA2')
         call GetVarFromMemToOutput ('accpr_mc', oneBramsGrid%currGrid, OutputField)
      case ('AA3')
         call GetVarFromMemToOutput ('accpr_st', oneBramsGrid%currGrid, OutputField)
      case ('AA1_BL')
         call GetVarFromMemToOutput ('accpr_as', oneBramsGrid%currGrid, OutputField)
      case ('L_CLD_FRAC')
         call GetVarFromMemToOutput ('CLOUD_FRACTION', oneBramsGrid%currGrid, ScrT3N01)
         call GetVarFromMemToOutput ('PI', oneBramsGrid%currGrid, ScrT3N02)
         call rams_comp_press(ScrT3N02) 
         call cloudFract(ScrT3N02,ScrT3N01,680.0,1000.0,Outputfield)
      case ('M_CLD_FRAC')
         call GetVarFromMemToOutput ('CLOUD_FRACTION', oneBramsGrid%currGrid, ScrT3N01)
         call GetVarFromMemToOutput ('PI', oneBramsGrid%currGrid, ScrT3N02)
         call rams_comp_press(ScrT3N02) 
         call cloudFract(ScrT3N02,ScrT3N01,440.0,680.0,Outputfield)
      case ('H_CLD_FRAC')
         call GetVarFromMemToOutput ('CLOUD_FRACTION', oneBramsGrid%currGrid, ScrT3N01)
         call GetVarFromMemToOutput ('PI', oneBramsGrid%currGrid, ScrT3N02)
         call rams_comp_press(ScrT3N02) 
         call cloudFract(ScrT3N02,ScrT3N01,0.0,440.0,Outputfield)
      case default
         write(*, "(a)") "**(OnePostField)** Post field 2d " // one_post_variable%fieldName // " not implemented!"
         stop
      end select


      if(IPOS == 11 .or. IPOS == 10) call writeTimeLineFRN(one_post_variable%fieldName,OutputField,oneBramsGrid%mxp &
         , oneBramsGrid%myp,time)

      if(IPOS==11) return
      ! most of cases run this. Some just returns before
      call PrepareAndOutputGradsField(one_post_variable, oneBramsGrid, onePostGrid, OutPutField)

   end subroutine Brams2Post_2d

   !=============================================================================================
   subroutine CloudFract(press,totalFract,minLevel,maxLevel,cFract)
       !# compute the cloud fract
       !#
       !# @note
       !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
       !#
       !# **Brief**: compute cloud fract beteween x1 minLevel and maxLevel [mBar]
       !#
       !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
       !#
       !# **Author(s)**: Luiz Flavio Rodrigues **&#9993;**<luiz.rodrigues@inpe.br>
       !#
       !# **Date**: 2021Apr08
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
       character(len=*),parameter :: sourceName='ModPostOneField2d.f90' !Name of source code
       character(len=*),parameter :: procedureName='**lowCloudFract**' !Name of this procedure
       !
       !Local Parameters
   
       !Input/Output variables
       real, intent(in) :: press(:,:,:)
       !# Input pressure values
       real, intent(in) :: totalFract(:,:,:)
       !# Input cloud_fract for all levels
       real, intent(in) :: minLevel
       !# min pressure value
       real, intent(in) :: maxLevel
       !# Max pressure value
       real, intent(out) :: cFract(:,:)
       
       !Local variables
       integer :: i,j,k
       integer :: nLevels
       real :: p1, p2, deltaPress
       logical :: firstTime
       
       !Code  
       cFract=0.0
       do j = 1, size(press, 2)
         do i = 1, size(press, 1)
            firstTime=.true.
            do k = 1, size(press, 3)
               !write(88,*) k,press(i,j,k)
               if(press(i,j,k)>minLevel .and. press(i,j,k)<maxLevel) then
                  deltaPress=press(i,j,k-1)-press(i,j,k)
                  if (firstTime) then
                     p1=press(i,j,k)
                     firstTime=.false.
                  endif
                  cFract(i,j)=cFract(i,j)+totalFract(i,j,k)*deltaPress
                  p2=press(i,j,k)
               endif
            enddo
            cFract(i,j)=cFract(i,j)/(p1-p2)
         enddo
      enddo
   
   end subroutine CloudFract 


   function getAconpr(oneBramsGrid) result(aconpr)
      use mem_cuparm, only : hasAconpr

      type(BramsGrid), intent(in), pointer :: oneBramsGrid
      real :: aconpr(oneBramsGrid%mxp, oneBramsGrid%myp)

      if(hasAconpr(oneBramsGrid%currGrid)) then
         call GetVarFromMemToOutput ('ACONPR', oneBramsGrid%currGrid, aconpr)
      else
         aconpr = 0
      end if

   end function getAconpr

   subroutine getAccComponents(oneBramsGrid, OutputField)

      type(BramsGrid), pointer :: oneBramsGrid
      real, intent(inout) :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp)

      real :: ScrT2N01(oneBramsGrid%mxp, oneBramsGrid%myp)
      ! OutputField <- 0.0
      ! ScrT2N01 <- field 'ACCPR' from var_table
      ! OutputField <- OutputField + ScrT2N01
      ! ScrT2N01 <- field 'ACCPP' from var_table
      ! OutputField <- OutputField + ScrT2N01
      ! ScrT2N01 <- field 'ACCPS' from var_table
      ! OutputField <- OutputField + ScrT2N01
      ! ScrT2N01 <- field 'ACCPA' from var_table
      ! OutputField <- OutputField + ScrT2N01
      ! ScrT2N01 <- field 'ACCPG' from var_table
      ! OutputField <- OutputField + ScrT2N01
      ! ScrT2N01 <- field 'ACCPH' from var_table
      ! OutputField <- OutputField + ScrT2N01
      ! OutputField <- max(OutputField, 0.0)

      OutputField = 0.0
      call GetVarFromMemToOutput ('ACCPR', oneBramsGrid%currGrid, ScrT2N01)
      OutputField = OutputField + ScrT2N01
      !-For GThompson microphysics micphys_type>1 : ACCPR is already the total precip
      if(mcphys_type .le. 1) then
         call GetVarFromMemToOutput ('ACCPP', oneBramsGrid%currGrid, ScrT2N01)
         OutputField = OutputField + ScrT2N01
         call GetVarFromMemToOutput ('ACCPS', oneBramsGrid%currGrid, ScrT2N01)
         OutputField = OutputField + ScrT2N01
         call GetVarFromMemToOutput ('ACCPA', oneBramsGrid%currGrid, ScrT2N01)
         OutputField = OutputField + ScrT2N01
         call GetVarFromMemToOutput ('ACCPG', oneBramsGrid%currGrid, ScrT2N01)
         OutputField = OutputField + ScrT2N01
         call GetVarFromMemToOutput ('ACCPH', oneBramsGrid%currGrid, ScrT2N01)
         OutputField = OutputField + ScrT2N01
      endif

   end subroutine getAccComponents


   subroutine getWinds10m(oneBramsGrid, u_or_v, OutputField)
      type(BramsGrid), pointer :: oneBramsGrid
      character(len = *), intent(in) :: u_or_v
      real, intent(out) :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp)
      integer :: firstX, lastX, firstY, lastY

      real :: ScrT2N01(oneBramsGrid%mxp, oneBramsGrid%myp)
      real :: ScrT3N01(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)
      real :: ScrT3N02(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)
      real :: ScrT3N03(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)

      real :: ScrT6N01(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%npatch)
      real :: ScrT6N02(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%npatch)
      real :: ScrT6N03(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%npatch)
      real :: ScrT6N04(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%npatch)

      call GetVarFromMemToOutput ('UP', oneBramsGrid%currGrid, ScrT3N01)
      call GetVarFromMemToOutput ('VP', oneBramsGrid%currGrid, ScrT3N02)
      call rams_comp_speed (ScrT3N01, ScrT3N02)
      call GetVarFromMemToOutput ('THETA', oneBramsGrid%currGrid, ScrT3N02)
      call GetVarFromMemToOutput ('PI', oneBramsGrid%currGrid, ScrT3N03)
      call GetVarFromMemToOutput ('TOPT', oneBramsGrid%currGrid, ScrT2N01)
      call GetVarFromMemToOutput ('USTAR', oneBramsGrid%currGrid, ScrT6N01)
      call GetVarFromMemToOutput ('PATCH_ROUGH', oneBramsGrid%currGrid, ScrT6N02)
      call GetVarFromMemToOutput ('CAN_TEMP', oneBramsGrid%currGrid, ScrT6N03)
      call GetVarFromMemToOutput ('PATCH_AREA', oneBramsGrid%currGrid, ScrT6N04)
      call rams_reduced_wind (OutputField, ScrT3N01, ScrT6N01, 10., oneBramsGrid%ztn(2), ScrT6N02, ScrT6N04, &
            ScrT6N03, ScrT3N02, ScrT3N03, ScrT2N01, oneBramsGrid%ztop)
      call GetVarFromMemToOutput ('UP', oneBramsGrid%currGrid, ScrT3N01)
      call GetVarFromMemToOutput ('VP', oneBramsGrid%currGrid, ScrT3N02)
      firstX = oneBramsGrid%nodei0(oneBramsGrid%mynum) + 1
      lastX = oneBramsGrid%nodei0(oneBramsGrid%mynum) + oneBramsGrid%mxp
      firstY = oneBramsGrid%nodej0(oneBramsGrid%mynum) + 1
      lastY = oneBramsGrid%nodej0(oneBramsGrid%mynum) + oneBramsGrid%myp
      call RAMS_comp_dir (ScrT3N01, ScrT3N02, oneBramsGrid%xtn(firstX : lastX), oneBramsGrid%ytn(firstY : lastY), &
            oneBramsGrid%polelat, oneBramsGrid%polelon)
      if(u_or_v .eq. 'U') then
         call calc_u10m (OutputField, ScrT3N01)
      else
         call calc_v10m (OutputField, ScrT3N02)
      end if
      !--- Indefinindo campos "zerados" - ex: analise ---
      if (maxval(OutputField)==0.0 .and. minval(OutputField)==0.0) OutputField = -9.99e+33
   end subroutine getWinds10m


   subroutine alloc2d_routine(mxp, myp)
      implicit none
      integer, intent(in) :: mxp,myp	        
  
  
      allocate(totprec (mxp,myp)); totprec= 0.0
      allocate(convprec(mxp,myp)); convprec=0.0
   end subroutine alloc2d_routine
end module ModPostOneField2d
