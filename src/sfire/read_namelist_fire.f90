  module ModNamelistsfireFile

  !use grid_dims, only: ngrid

   
   public :: namelistsfireFile
   public :: CreateNamelistsfireFile
   public :: DestroyNamelistsfireFile
   public :: GetNamelistsfireFileName
   public :: ReadNamelistsfireFile
  !public :: BroadcastNamelistsfireFile
   public :: DumpNamelistsfireFile
  !public :: TimeUnitsToSeconds


  !include "files.h"
  integer, parameter :: ngrid=1
  integer, parameter :: f_name_length=11
  type grid_config_rec_type
 
 character(len=f_name_length) :: fileName 

     ! namelist RAMSIN_FIRE

     include "namelist_defines2.inc"

     integer ims,ime, kms,kme, jms,jme &
			,ifds,ifde, jfds,jfde, kfds, kfde &
			,ifms,ifme, jfms,jfme, kfms, kfme &
			,ifps,ifpe, jfps,jfpe, kfps, kfpe &
			,ifts,ifte, jfts,jfte, kfts, kfte &
			,ids,ide, jds,jde, kds,kde                        & 
        		,ips,ipe, jps,jpe, kps,kpe			   & 
        		,its,ite, jts,jte, kts,kte
    
    	
     LOGICAL  :: start=.TRUE.
 
 end type grid_config_rec_type
  
contains


subroutine CreateNamelistsfireFile(oneNamelistFile)
    implicit none
    type(grid_config_rec_type), pointer :: oneNamelistFile
    if (associated(oneNamelistFile)) then
       deallocate(oneNamelistFile)
    end if
    allocate(oneNamelistFile)
    print*,'createnamelist'
  end subroutine CreateNamelistsfireFile




  subroutine DestroyNamelistsfireFile(oneNamelistFile)
    implicit none
    type(grid_config_rec_type), pointer :: oneNamelistFile
    if (associated(oneNamelistFile)) then
       deallocate(oneNamelistFile)
    end if
    nullify(oneNamelistFile)
    print*,'destroinamelist'
  end subroutine DestroyNamelistsfireFile


   subroutine GetNamelistsfireFileName(oneNamelistFile)
    implicit none
    type(grid_config_rec_type), pointer :: oneNamelistFile

    integer :: nargs
    integer :: iarg
    integer :: lenArg
    integer :: status
    logical :: flagName ! true iff arg="-f"; next arg is file name
    integer, parameter :: f_name_length=11
    character(len=f_name_length) :: arg
    character(len=f_name_length) ::fileName
    oneNamelistFile%fileName="sfire.in" ! default namelist

    ! search command line for "-f " <namelist file name> 
    ! return default if not found

  !  nargs = command_argument_count()
  !  if (nargs >= 0) then
   !    flagName = .false.
    !   do iarg = 0, nargs
     !     call get_command_argument(iarg, arg, lenArg, status)
     !     if (status == 0) then
     !        if (flagName) then
     !           oneNamelistFile%fileName = arg(1:lenArg)
      !          exit
      !       else
      !          flagName = arg(1:lenArg) == "-f"
      !       end if
      !    end if
      ! end do
   ! end if
  end subroutine GetNamelistsfireFileName


  ! ReadNamelistFile:
  !    open, reads and close namelist file
  !    implements defaults for namelist variables
  !    check input options consistency


  subroutine ReadNamelistsfireFile(oneNamelistFile)

    implicit none
    type(grid_config_rec_type), pointer :: oneNamelistFile

  !  include "files.h"


    integer :: iunit                    ! io unit number
    integer, parameter :: firstUnit=14  ! lowest io unit number available
    integer, parameter :: lastUnit=15
   ! integer, parameter :: firstUnit=20  ! lowest io unit number available
   ! integer, parameter :: lastUnit=99   ! highest io unit number available
    logical :: op                       ! io unit number opened or not
    logical :: ex                       ! namelist file exists?
    integer :: err                      ! return code on iostat
    character(len=10) :: c0             ! scratch
    character(len=*), parameter :: h="**(ReadNamelistsfireFile)**"  ! program unit name

    ! namelist /FIRE/

    real:: dx
    real:: dy
    real:: dt
    integer ::tracer_opt
    real ::cen_lat
    real ::cen_lon
    logical ::restart
    integer ::sr_x
    integer ::sr_y
    logical :: fmoist_run
    logical :: fmoist_interp
    logical :: fmoist_only
    integer :: fmoist_freq
    real :: fmoist_dt
    integer :: ifire
    logical :: crown
    integer :: fire_boundary_guard
    integer :: fire_num_ignitions
    real :: fire_ignition_ros1
    real :: fire_ignition_start_lon1
    real :: fire_ignition_start_lat1
    real :: fire_ignition_end_lon1
    real :: fire_ignition_end_lat1
    real :: fire_ignition_radius1
    real :: fire_ignition_start_time1
    real :: fire_ignition_end_time1
    real :: fire_ignition_ros2
    real :: fire_ignition_start_lon2
    real :: fire_ignition_start_lat2
    real :: fire_ignition_end_lon2
    real :: fire_ignition_end_lat2
    real :: fire_ignition_radius2
    real :: fire_ignition_start_time2
    real :: fire_ignition_end_time2
    real :: fire_ignition_ros3
    real :: fire_ignition_start_lon3
    real :: fire_ignition_start_lat3
    real :: fire_ignition_end_lon3
    real :: fire_ignition_end_lat3
    real :: fire_ignition_radius3
    real :: fire_ignition_start_time3
    real :: fire_ignition_end_time3
    real :: fire_ignition_ros4
    real :: fire_ignition_start_lon4
    real :: fire_ignition_start_lat4
    real :: fire_ignition_end_lon4 
    real :: fire_ignition_end_lat4
    real :: fire_ignition_radius4
    real :: fire_ignition_start_time4
    real :: fire_ignition_end_time4
    real :: fire_ignition_ros5 
    real :: fire_ignition_start_lon5
    real :: fire_ignition_start_lat5
    real :: fire_ignition_end_lon5
    real :: fire_ignition_end_lat5
    real :: fire_ignition_radius5
    real :: fire_ignition_start_time5
    real :: fire_ignition_end_time5
    real :: fire_ignition_start_x1
    real :: fire_ignition_start_y1
    real :: fire_ignition_end_x1
    real :: fire_ignition_end_y1
    real :: fire_ignition_start_x2
    real :: fire_ignition_start_y2
    real :: fire_ignition_end_x2
    real :: fire_ignition_end_y2
    real :: fire_ignition_start_x3
    real :: fire_ignition_start_y3
    real :: fire_ignition_end_x3
    real :: fire_ignition_end_y3
    real :: fire_ignition_start_x4
    real :: fire_ignition_start_y4
    real :: fire_ignition_end_x4
    real :: fire_ignition_end_y4
    real :: fire_ignition_start_x5
    real :: fire_ignition_start_y5
    real :: fire_ignition_end_x5
    real :: fire_ignition_end_y5
    real :: fire_perimeter_time
    real :: fire_lat_init
    real :: fire_lon_init
    real :: fire_ign_time
    integer ::fire_shape
    integer ::fire_sprd_mdl
    real ::fire_crwn_hgt
    real ::fire_ext_grnd
    real ::fire_ext_crwn
    integer ::fire_wind_log_interp
    integer ::fire_use_windrf
    integer ::fire_fuel_read
    integer ::fire_fmc_read
    integer ::fire_fuel_cat
    integer ::fire_print_msg
    integer ::fire_print_file
    logical ::fire_restart
    integer ::fire_time_step_ratio
    integer ::fire_debug_hook_sec
    integer ::fire_fuel_left_method
    integer ::fire_fuel_left_irl
    integer ::fire_fuel_left_jrl
    real :: fire_back_weight
    integer ::fire_grows_only
    integer ::fire_upwinding
    real :: fire_viscosity
    real :: fire_lfn_ext_up
    integer :: fire_topo_from_atm
    integer :: fire_advection
    integer :: fire_test_steps
    real :: fire_const_time
    real :: fire_const_grnhfx
    real :: fire_const_grnqfx
    integer :: fire_hfx_given
    integer :: fire_hfx_num_lines
    real :: fire_hfx_latent_part
    real :: fire_hfx_value1
    real :: fire_hfx_start_time1
    real :: fire_hfx_end_time1
    real :: fire_hfx_trans_time1
    real :: fire_hfx_radius1
    real :: fire_hfx_start_x1
    real :: fire_hfx_end_x1
    real :: fire_hfx_start_lat1
    real :: fire_hfx_end_lat1
    real :: fire_hfx_start_y1
    real :: fire_hfx_end_y1
    real :: fire_hfx_start_lon1
    real :: fire_hfx_end_lon1
    real :: fire_atm_feedback
    integer :: chem_opt
    integer :: nfmc
    integer :: fire_ignition_clamp
    integer :: fire_update_fuel_frac
    integer :: fndwi_from_ndwi
    integer :: kfmc_ndwi
    integer :: fire_can_top_read
    integer :: sfire_upwinding

   namelist /FIRE_DEFAULT/                                        &
         dx,dy,dt,tracer_opt,cen_lat,cen_lon,restart,sr_x,sr_y, &
         fmoist_run,fmoist_interp,fmoist_only,fmoist_freq,fmoist_dt, &
         fire_ignition_ros4,fire_ignition_start_lon4, &
         fire_ignition_start_lat4,fire_ignition_end_lon4, &
         fire_ignition_end_lat4,fire_ignition_radius4, &
         fire_ignition_start_time4,fire_ignition_end_time4, &
         fire_ignition_ros5,fire_ignition_start_lon5, &
         fire_ignition_start_lat5,fire_ignition_end_lon5, &
         fire_ignition_end_lat5,fire_ignition_radius5, &
         fire_ignition_start_time5,fire_ignition_end_time5, &
         fire_ignition_start_x1,fire_ignition_start_y1, &
         fire_ignition_end_x1,fire_ignition_end_y1, &
         fire_ignition_start_x2,fire_ignition_start_y2, &
         fire_ignition_end_x2,fire_ignition_end_y2, &
         fire_ignition_start_x3,fire_ignition_start_y3, &
         fire_ignition_end_x3,fire_ignition_end_y3, &
         fire_ignition_start_x4,fire_ignition_start_y4, &
         fire_ignition_end_x4,fire_ignition_end_y4, &
         fire_ignition_start_x5,fire_ignition_start_y5, &
         fire_ignition_end_x5,fire_ignition_end_y5, &
         fire_perimeter_time,fire_lat_init,fire_lon_init, &
         fire_ign_time,fire_shape,fire_sprd_mdl,fire_crwn_hgt, &
         fire_ext_grnd,fire_ext_crwn,fire_wind_log_interp, &
         fire_use_windrf,fire_fmc_read, &
         fire_restart,fire_time_step_ratio,fire_debug_hook_sec, &
         fire_back_weight,fire_advection, &
         fire_const_time,fire_const_grnhfx,fire_const_grnqfx, &
         fire_hfx_given,fire_hfx_num_lines,fire_hfx_latent_part, &
         fire_hfx_value1,fire_hfx_start_time1, &
         fire_hfx_end_time1,fire_hfx_trans_time1, &
         fire_hfx_radius1,fire_hfx_start_x1,fire_hfx_end_x1, &
         fire_hfx_start_lat1,fire_hfx_end_lat1,fire_hfx_start_y1, &
         fire_hfx_end_y1,fire_hfx_start_lon1,fire_hfx_end_lon1, &
         chem_opt,nfmc,fire_ignition_clamp,fire_update_fuel_frac,&
         fndwi_from_ndwi,kfmc_ndwi,fire_can_top_read,sfire_upwinding


     namelist /FIRE_SIMUL_INFO/ &
         ifire,crown,fire_fuel_read,fire_fuel_cat,fire_num_ignitions, &       
         fire_ignition_start_lon1,fire_ignition_start_lat1, & 
         fire_ignition_end_lon1,fire_ignition_end_lat1, &
         fire_ignition_radius1,fire_ignition_ros1, & 
         fire_ignition_start_time1,fire_ignition_end_time1, & 
         fire_ignition_start_lon2,fire_ignition_start_lat2, &
         fire_ignition_end_lon2,fire_ignition_end_lat2, & 
         fire_ignition_radius2,fire_ignition_ros2, & 
         fire_ignition_start_time2,fire_ignition_end_time2, & 
         fire_ignition_start_lon3,fire_ignition_start_lat3, & 
         fire_ignition_end_lon3,fire_ignition_end_lat3, &  
         fire_ignition_radius3,fire_ignition_ros3, &
         fire_ignition_start_time3,fire_ignition_end_time3, & 
         fire_print_msg,fire_print_file,fire_boundary_guard, &   
         fire_fuel_left_method,fire_fuel_left_irl, &      
         fire_fuel_left_jrl,fire_atm_feedback,fire_grows_only, &        
         fire_viscosity,fire_upwinding,fire_lfn_ext_up, &        
         fire_test_steps,fire_topo_from_atm


    !FIRE_DEFAULT

    dx                                  = 200.
    dy                                  = 200.
    dt                                  = 2.
    tracer_opt                          = 0
    cen_lat                             = 39.7053
    cen_lon                             = -107.2907
    restart                             = .false.
    sr_x                                = 0
    sr_y                                = 0
    fmoist_run                          = .false.
    fmoist_interp                       = .false.
    fmoist_only                         = .false.
    fmoist_freq                         = 0
    fmoist_dt                           = 600
    fire_ignition_ros4                  = 0.01
    fire_ignition_start_lon4            = 0.
    fire_ignition_start_lat4            = 0.
    fire_ignition_end_lon4              = 0.
    fire_ignition_end_lat4              = 0.
    fire_ignition_radius4               = 0.
    fire_ignition_start_time4           = 0.
    fire_ignition_end_time4             = 0.
    fire_ignition_ros5                  = 0.01
    fire_ignition_start_lon5            = 0.
    fire_ignition_start_lat5            = 0.
    fire_ignition_end_lon5              = 0.
    fire_ignition_end_lat5              = 0.
    fire_ignition_radius5               = 0.
    fire_ignition_start_time5           = 0.
    fire_ignition_end_time5             = 0.
    fire_ignition_start_x1              = 0.
    fire_ignition_start_y1              = 0.
    fire_ignition_end_x1                = 0.
    fire_ignition_end_y1                = 0.
    fire_ignition_start_x2              = 0.
    fire_ignition_start_y2              = 0.
    fire_ignition_end_x2                = 0.
    fire_ignition_end_y2                = 0.
    fire_ignition_start_x3              = 0.
    fire_ignition_start_y3              = 0.
    fire_ignition_end_x3                = 0.
    fire_ignition_end_y3                = 0.
    fire_ignition_start_x4              = 0.
    fire_ignition_start_y4              = 0.
    fire_ignition_end_x4                = 0.
    fire_ignition_end_y4                = 0.
    fire_ignition_start_x5              = 0.
    fire_ignition_start_y5              = 0.
    fire_ignition_end_x5                = 0.
    fire_ignition_end_y5                = 0.
    fire_perimeter_time                 = 0.
    fire_lat_init                       = 0.
    fire_lon_init                       = 0.
    fire_ign_time                       = 0.
    fire_shape                          = 0
    fire_sprd_mdl                       = 1
    fire_crwn_hgt                       = 15.
    fire_ext_grnd                       = 50.
    fire_ext_crwn                       = 50.
    fire_wind_log_interp                = 4
    fire_use_windrf                     = 0
    fire_fmc_read                       = 1
    fire_restart                        = .false.
    fire_time_step_ratio                = 1
    fire_debug_hook_sec                 = 0
    fire_back_weight                    = 0.5
    fire_advection                      = 1
    fire_const_time                     = -1.
    fire_const_grnhfx                   = 0.
    fire_const_grnqfx                   = 0.
    fire_hfx_given                      = 0
    fire_hfx_num_lines                  = 0
    fire_hfx_latent_part                = 0.084
    fire_hfx_value1                     = 0.
    fire_hfx_start_time1                = 0.
    fire_hfx_end_time1                  = 0.
    fire_hfx_trans_time1                = 0.
    fire_hfx_radius1                    = 0.
    fire_hfx_start_x1                   = 0.
    fire_hfx_end_x1                     = 0.
    fire_hfx_start_lat1                 = 0.
    fire_hfx_end_lat1                   = 0.
    fire_hfx_start_y1                   = 0.
    fire_hfx_end_y1                     = 0.
    fire_hfx_start_lon1                 = 0.
    fire_hfx_end_lon1                   = 0.
    chem_opt                            = 0
    nfmc                                = 5
    fire_ignition_clamp                 = 0
    fire_update_fuel_frac               = 1
    fndwi_from_ndwi                     = 1
    kfmc_ndwi                           = 1
    fire_can_top_read                   = 1
    sfire_upwinding                     = 3

    !FIRE_SIMUL_INFO
 
    ifire                               = 2 
    crown                               = .true.
    fire_fuel_read                      = -1
    fire_fuel_cat                       = 3
    fire_num_ignitions                  = 3       
    fire_ignition_start_lon1            = -107.293664
    fire_ignition_start_lat1            = 39.698696
    fire_ignition_end_lon1              = -107.293664 
    fire_ignition_end_lat1              = 39.710990 
    fire_ignition_radius1               = 370.
    fire_ignition_ros1 = 10. 
    fire_ignition_start_time1           = 2
    fire_ignition_end_time1             = 2 
    fire_ignition_start_lon2            = -107.287954  
    fire_ignition_start_lat2            =  39.698696 
    fire_ignition_end_lon2              = -107.287954 
    fire_ignition_end_lat2              = 39.71099 
    fire_ignition_radius2               = 370.
    fire_ignition_ros2 = 10. 
    fire_ignition_start_time2           = 3 
    fire_ignition_end_time2             = 3 
    fire_ignition_start_lon3            = -107.289096 
    fire_ignition_start_lat3            = 39.706599 
    fire_ignition_end_lon3              = -107.289096  
    fire_ignition_end_lat3              = 39.706599  
    fire_ignition_radius3               = 400. 
    fire_ignition_ros3                  = 10.
    fire_ignition_start_time3           = 4
    fire_ignition_end_time3             = 4 
    fire_print_msg                      = 1  
    fire_print_file                     = 1   
    fire_boundary_guard                 = -1   
    fire_fuel_left_method               = 1  
    fire_fuel_left_irl                  = 2      
    fire_fuel_left_jrl                  = 2      
    fire_atm_feedback                   = 1.        
    fire_grows_only                     = 1       
    fire_viscosity                      = 0.4          
    fire_upwinding                      = 3        
    fire_lfn_ext_up                     = 1.0        
    fire_test_steps                     = 0         
    fire_topo_from_atm                  = 0                  


    print*,'readnamelist 1'
! select unused i/o unit

    do iunit = firstUnit, lastUnit
       inquire(iunit,opened=op)
       if (.not. op) exit
    end do

    if (iunit > lastUnit) then
    !   call fatal_error(h//" all i/o units in use")
      print*, "all i/o units in use"
    end if

    ! if namelist file exists, open, read each section and close

    inquire(file=trim(oneNamelistFile%fileName), exist=ex)
    if (.not. ex) then
     !  call fatal_error(h//" namelist file "//trim(oneNamelistFile%fileName)//&
     !       " does not exist")
     PRINT*, "ERRO :: ARQUIVO NAO EXISTE"
    end if

    open(iunit, file=trim(oneNamelistFile%fileName), status="old", action="read",&
         iostat=err)
    if (err /= 0) then
  !     write(c0,"(i10)") err
  !     call fatal_error(h//" open namelist file "//trim(oneNamelistFile%fileName)//&
   !         " returned iostat="//trim(adjustl(c0)))

      PRINT*, "ERRO :: NAO LIDO "//trim(oneNamelistFile%fileName)
    end if
    print*,'fire_default'
    read (iunit, iostat=err, NML=FIRE_DEFAULT)
    print*,'ERR=',err
    if (err /= 0) then
       write(*,"(a)") h//"**(ERROR)** reading section FIRE_DEFAULT "//&
            &"of namelist file "//trim(oneNamelistFile%fileName)
       write(*,"(a)") h//" compare values read with file contents:"

       write(*,*) "dx=",dx
       write(*,*) "dy=",dy
       write(*,*) "dt=",dt
       write(*,*) "tracer_opt=",tracer_opt
       write(*,*) "cen_lat=",cen_lat
       write(*,*) "cen_lon=",cen_lon
       write(*,*) "restart=",restart
       write(*,*) "sr_x=",sr_x
       write(*,*) "sr_y=",sr_y
       write(*,*) "fmoist_run=",fmoist_run
       write(*,*) "fmoist_interp=",fmoist_interp
       write(*,*) "fmoist_only=",fmoist_only
       write(*,*) "fmoist_freq=",fmoist_freq
       write(*,*) "fmoist_dt=", fmoist_dt
       write(*,*) "fire_ignition_ros4=",fire_ignition_ros4
       write(*,*) "fire_ignition_start_lon4=",fire_ignition_start_lon4
       write(*,*) "fire_ignition_start_lat4=",fire_ignition_start_lat4
       write(*,*) "fire_ignition_end_lon4=",fire_ignition_end_lon4
       write(*,*) "fire_ignition_end_lat4=",fire_ignition_end_lat4
       write(*,*) "fire_ignition_radius4=",fire_ignition_radius4
       write(*,*) "fire_ignition_start_time4=",fire_ignition_start_time4
       write(*,*) "fire_ignition_end_time4=",fire_ignition_end_time4
       write(*,*) "fire_ignition_ros5=",fire_ignition_ros5
       write(*,*) "fire_ignition_start_lon5=",fire_ignition_start_lon5
       write(*,*) "fire_ignition_start_lat5=",fire_ignition_start_lat5
       write(*,*) "fire_ignition_end_lon5=",fire_ignition_end_lon5
       write(*,*) "fire_ignition_end_lat5=",fire_ignition_end_lat5
       write(*,*) "fire_ignition_radius5=",fire_ignition_radius5
       write(*,*) "fire_ignition_start_time5=",fire_ignition_start_time5
       write(*,*) "fire_ignition_end_time5=",fire_ignition_end_time5
       write(*,*) "fire_ignition_start_x1=",fire_ignition_start_x1
       write(*,*) "fire_ignition_start_y1=",fire_ignition_start_y1
       write(*,*) "fire_ignition_end_x1=",fire_ignition_end_x1
       write(*,*) "fire_ignition_end_y1=",fire_ignition_end_y1
       write(*,*) "fire_ignition_start_x2=",fire_ignition_start_x2
       write(*,*) "fire_ignition_start_y2=",fire_ignition_start_y2
       write(*,*) "fire_ignition_end_x2=",fire_ignition_end_x2
       write(*,*) "fire_ignition_end_y2=",fire_ignition_end_y2
       write(*,*) "fire_ignition_start_x3=",fire_ignition_start_x3
       write(*,*) "fire_ignition_start_y3=",fire_ignition_start_y3
       write(*,*) "fire_ignition_end_x3=",fire_ignition_end_x3
       write(*,*) "fire_ignition_end_y3=",fire_ignition_end_y3
       write(*,*) "fire_ignition_start_x4=",fire_ignition_start_x4  
       write(*,*) "fire_ignition_start_y4=", fire_ignition_start_y4
       write(*,*) "fire_ignition_end_x4=",fire_ignition_end_x4
       write(*,*) "fire_ignition_end_y4=",fire_ignition_end_y4
       write(*,*) "fire_ignition_start_x5=",fire_ignition_start_x5
       write(*,*) "fire_ignition_start_y5=",fire_ignition_start_y5
       write(*,*) "fire_ignition_end_x5=",fire_ignition_end_x5
       write(*,*) "fire_ignition_end_y5=",fire_ignition_end_y5
       write(*,*) "fire_perimeter_time=",fire_perimeter_time
       write(*,*) "fire_lat_init=",fire_lat_init
       write(*,*) "fire_lon_init=",fire_lon_init
       write(*,*) "fire_ign_time=",fire_ign_time
       write(*,*) "fire_shape=",fire_shape
       write(*,*) "fire_sprd_mdl=",fire_sprd_mdl
       write(*,*) "fire_crwn_hgt=",fire_crwn_hgt
       write(*,*) "fire_ext_grnd=",fire_ext_grnd
       write(*,*) "fire_ext_crwn=",fire_ext_crwn
       write(*,*) "fire_wind_log_interp=",fire_wind_log_interp
       write(*,*) "fire_use_windrf=",fire_use_windrf
       write(*,*) "fire_fmc_read=",fire_fmc_read
       write(*,*) "fire_restart=",fire_restart
       write(*,*) "fire_time_step_ratio=",fire_time_step_ratio
       write(*,*) "fire_debug_hook_sec=",fire_debug_hook_sec
       write(*,*) "fire_back_weight=",fire_back_weight
       write(*,*) "fire_advection=",fire_advection
       write(*,*) "fire_const_time=",fire_const_time
       write(*,*) "fire_const_grnhfx=",fire_const_grnhfx
       write(*,*) "fire_const_grnqfx=",fire_const_grnqfx
       write(*,*) "fire_hfx_given=",fire_hfx_given 
       write(*,*) "fire_hfx_num_lines=",fire_hfx_num_lines
       write(*,*) "fire_hfx_latent_part=",fire_hfx_latent_part
       write(*,*) "fire_hfx_value1=",fire_hfx_value1
       write(*,*) "fire_hfx_start_time1=", fire_hfx_start_time1
       write(*,*) "fire_hfx_end_time1=",fire_hfx_end_time1
       write(*,*) "fire_hfx_trans_time1=",fire_hfx_trans_time1
       write(*,*) "fire_hfx_radius1=",fire_hfx_radius1
       write(*,*) "fire_hfx_start_x1=",fire_hfx_start_x1
       write(*,*) "fire_hfx_end_x1=",fire_hfx_end_x1
       write(*,*) "fire_hfx_start_lat1=",fire_hfx_start_lat1
       write(*,*) "fire_hfx_end_lat1=",fire_hfx_end_lat1
       write(*,*) "fire_hfx_start_y1=",fire_hfx_start_y1
       write(*,*) "fire_hfx_end_y1=",fire_hfx_end_y1
       write(*,*) "fire_hfx_start_lon1=",fire_hfx_start_lon1
       write(*,*) "fire_hfx_end_lon1=",fire_hfx_end_lon1
       write(*,*) "chem_opt=",chem_opt
       write(*,*) "nfmc=",nfmc
       write(*,*) "fire_ignition_clamp=",fire_ignition_clamp
       write(*,*) "fire_update_fuel_frac=",fire_update_fuel_frac
       write(*,*) "fndwi_from_ndwi=",fndwi_from_ndwi
       write(*,*) "kfmc_ndwi=",kfmc_ndwi
       write(*,*) "fire_can_top_read=",fire_can_top_read
       write(*,*) "sfire_upwinding=",sfire_upwinding

      ! call fatal_error(h//" reading namelist")
    else
       ! namelist /FIRE_DEFAULT/
     print*, 'isilda1'
     oneNamelistFile%dx = dx
     oneNamelistFile%dy = dy
     oneNamelistFile%dt = dt
     oneNamelistFile%tracer_opt = tracer_opt
     oneNamelistFile%cen_lat = cen_lat
     oneNamelistFile%cen_lon = cen_lon
     oneNamelistFile%restart = restart
     oneNamelistFile%sr_x = sr_x
     oneNamelistFile%sr_y = sr_y
     oneNamelistFile%fmoist_run = fmoist_run
     oneNamelistFile%fmoist_interp = fmoist_interp
     oneNamelistFile%fmoist_only =fmoist_only
     oneNamelistFile%fmoist_freq =fmoist_freq
     oneNamelistFile%fmoist_dt =fmoist_dt
     oneNamelistFile%fire_ignition_ros4 = fire_ignition_ros4
     oneNamelistFile%fire_ignition_start_lon4 = fire_ignition_start_lon4
     oneNamelistFile%fire_ignition_start_lat4 = fire_ignition_start_lat4
     oneNamelistFile%fire_ignition_end_lon4 = fire_ignition_end_lon4
     oneNamelistFile%fire_ignition_end_lat4 = fire_ignition_end_lat4
     oneNamelistFile%fire_ignition_radius4 = fire_ignition_radius4
     oneNamelistFile%fire_ignition_start_time4 = fire_ignition_start_time4
     oneNamelistFile%fire_ignition_end_time4 = fire_ignition_end_time4
     oneNamelistFile%fire_ignition_ros5 = fire_ignition_ros5
     oneNamelistFile%fire_ignition_start_lon5 = fire_ignition_start_lon5
     oneNamelistFile%fire_ignition_start_lat5 = fire_ignition_start_lat5
     oneNamelistFile%fire_ignition_end_lon5 = fire_ignition_end_lon5
     oneNamelistFile%fire_ignition_end_lat5 = fire_ignition_end_lat5
     oneNamelistFile%fire_ignition_radius5 = fire_ignition_radius5
     oneNamelistFile%fire_ignition_start_time5 = fire_ignition_start_time5
     oneNamelistFile%fire_ignition_end_time5 = fire_ignition_end_time5
     oneNamelistFile%fire_ignition_start_x1 = fire_ignition_start_x1
     oneNamelistFile%fire_ignition_start_y1 = fire_ignition_start_y1
     oneNamelistFile%fire_ignition_end_x1 = fire_ignition_end_x1
     oneNamelistFile%fire_ignition_end_y1 = fire_ignition_end_y1
     oneNamelistFile%fire_ignition_start_x2 = fire_ignition_start_x2
     oneNamelistFile%fire_ignition_start_y2 = fire_ignition_start_y2
     oneNamelistFile%fire_ignition_end_x2 = fire_ignition_end_x2
     oneNamelistFile%fire_ignition_end_y2 = fire_ignition_end_y2
     oneNamelistFile%fire_ignition_start_x3 = fire_ignition_start_x3
     oneNamelistFile%fire_ignition_start_y3 = fire_ignition_start_y3
     oneNamelistFile%fire_ignition_end_x3 = fire_ignition_end_x3
     oneNamelistFile%fire_ignition_end_y3 = fire_ignition_end_y3
     oneNamelistFile%fire_ignition_start_x4 = fire_ignition_start_x4
     oneNamelistFile%fire_ignition_start_y4 = fire_ignition_start_y4
     oneNamelistFile%fire_ignition_end_x4 = fire_ignition_end_x4
     oneNamelistFile%fire_ignition_end_y4 = fire_ignition_end_y4
     oneNamelistFile%fire_ignition_start_x5 = fire_ignition_start_x5
     oneNamelistFile%fire_ignition_start_y5 = fire_ignition_start_y5
     oneNamelistFile%fire_ignition_end_x5 = fire_ignition_end_x5
     oneNamelistFile%fire_ignition_end_y5 = fire_ignition_end_y5
     oneNamelistFile%fire_perimeter_time = fire_perimeter_time
     oneNamelistFile%fire_lat_init = fire_lat_init
     oneNamelistFile%fire_lon_init = fire_lon_init
     oneNamelistFile%fire_ign_time = fire_ign_time
     oneNamelistFile%fire_shape = fire_shape
     oneNamelistFile%fire_sprd_mdl = fire_sprd_mdl
     oneNamelistFile%fire_crwn_hgt = fire_crwn_hgt
     oneNamelistFile%fire_ext_grnd = fire_ext_grnd
     oneNamelistFile%fire_ext_crwn = fire_ext_crwn
     oneNamelistFile%fire_wind_log_interp = fire_wind_log_interp
     oneNamelistFile%fire_use_windrf = fire_use_windrf
     oneNamelistFile%fire_fmc_read = fire_fmc_read
     oneNamelistFile%fire_restart = fire_restart
     oneNamelistFile%fire_time_step_ratio = fire_time_step_ratio
     oneNamelistFile%fire_debug_hook_sec = fire_debug_hook_sec
     oneNamelistFile%fire_back_weight = fire_back_weight
     oneNamelistFile%fire_advection = fire_advection
     oneNamelistFile%fire_const_time = fire_const_time
     oneNamelistFile%fire_const_grnhfx = fire_const_grnhfx
     oneNamelistFile%fire_const_grnqfx = fire_const_grnqfx
     oneNamelistFile%fire_hfx_given = fire_hfx_given
     oneNamelistFile%fire_hfx_num_lines = fire_hfx_num_lines
     oneNamelistFile%fire_hfx_latent_part = fire_hfx_latent_part
     oneNamelistFile%fire_hfx_value1 = fire_hfx_value1
     oneNamelistFile%fire_hfx_start_time1 = fire_hfx_start_time1
     oneNamelistFile%fire_hfx_end_time1 = fire_hfx_end_time1
     oneNamelistFile%fire_hfx_trans_time1 = fire_hfx_trans_time1
     oneNamelistFile%fire_hfx_radius1 = fire_hfx_radius1
     oneNamelistFile%fire_hfx_start_x1 = fire_hfx_start_x1
     oneNamelistFile%fire_hfx_end_x1 = fire_hfx_end_x1
     oneNamelistFile%fire_hfx_start_lat1 = fire_hfx_start_lat1
     oneNamelistFile%fire_hfx_end_lat1 = fire_hfx_end_lat1
     oneNamelistFile%fire_hfx_start_y1 = fire_hfx_start_y1
     oneNamelistFile%fire_hfx_end_y1 = fire_hfx_end_y1
     oneNamelistFile%fire_hfx_start_lon1 = fire_hfx_start_lon1
     oneNamelistFile%fire_hfx_end_lon1 = fire_hfx_end_lon1
     oneNamelistFile%chem_opt = chem_opt
     oneNamelistFile%nfmc = nfmc
     oneNamelistFile%fire_ignition_clamp = fire_ignition_clamp
     oneNamelistFile%fire_update_fuel_frac = fire_update_fuel_frac
     oneNamelistFile%fndwi_from_ndwi = fndwi_from_ndwi
     oneNamelistFile%kfmc_ndwi = kfmc_ndwi
     oneNamelistFile%fire_can_top_read = fire_can_top_read
     oneNamelistFile%sfire_upwinding = sfire_upwinding

    end if

     print*,'fire_simul_info'
     read (iunit, iostat=err, NML=FIRE_SIMUL_INFO)
    if (err /= 0) then
       write(*,"(a)") h//"**(ERROR)** reading section FIRE_SIMUL_INFO "//&
            &"of namelist file "//trim(oneNamelistFile%fileName)
       write(*,"(a)") h//" compare values read with file contents:"
       
       
       write(*,*) "ifire=", ifire
       write(*,*) "crown=", crown
       write(*,*) "fire_fuel_read=", fire_fuel_read
       write(*,*) "fire_fuel_cat=", fire_fuel_cat
       write(*,*) "fire_num_ignitions=", fire_num_ignitions
       write(*,*) "fire_ignition_start_lon1=", fire_ignition_start_lon1
       write(*,*) "fire_ignition_start_lat1=", fire_ignition_start_lat1
       write(*,*) "fire_ignition_end_lon1=", fire_ignition_end_lon1
       write(*,*) "fire_ignition_end_lat1=", fire_ignition_end_lat1
       write(*,*) "fire_ignition_radius1=", fire_ignition_radius1
       write(*,*) "fire_ignition_ros1=", fire_ignition_ros1
       write(*,*) "fire_ignition_start_time1=", fire_ignition_start_time1
       write(*,*) "fire_ignition_end_time1=", fire_ignition_end_time1
       write(*,*) "fire_ignition_start_lon2=", fire_ignition_start_lon2
       write(*,*) "fire_ignition_start_lat2=", fire_ignition_start_lat2
       write(*,*) "fire_ignition_end_lon2=", fire_ignition_end_lon2
       write(*,*) "fire_ignition_end_lat2=", fire_ignition_end_lat2
       write(*,*) "fire_ignition_radius2=", fire_ignition_radius2
       write(*,*) "fire_ignition_ros2=", fire_ignition_ros2 
       write(*,*) "fire_ignition_start_time2=", fire_ignition_start_time2
       write(*,*) "fire_ignition_end_time2=", fire_ignition_end_time2
       write(*,*) "fire_ignition_start_lon3=", fire_ignition_start_lon3
       write(*,*) "fire_ignition_start_lat3=", fire_ignition_start_lat3
       write(*,*) "fire_ignition_end_lon3=", fire_ignition_end_lon3
       write(*,*) "fire_ignition_end_lat3=", fire_ignition_end_lat3
       write(*,*) "fire_ignition_radius3=", fire_ignition_radius3
       write(*,*) "fire_ignition_ros3=", fire_ignition_ros3
       write(*,*) "fire_ignition_start_time3=", fire_ignition_start_time3
       write(*,*) "fire_ignition_end_time3=", fire_ignition_end_time3
       write(*,*) "fire_print_msg=", fire_print_msg
       write(*,*) "fire_print_file=", fire_print_file
       write(*,*) "fire_boundary_guard=", fire_boundary_guard
       write(*,*) "fire_fuel_left_method=", fire_fuel_left_method
       write(*,*) "fire_fuel_left_irl=", fire_fuel_left_irl
       write(*,*) "fire_fuel_left_jrl=", fire_fuel_left_jrl
       write(*,*) "fire_atm_feedback=", fire_atm_feedback
       write(*,*) "fire_grows_only=", fire_grows_only       
       write(*,*) "fire_viscosity=", fire_viscosity
       write(*,*) "fire_upwinding=", fire_upwinding
       write(*,*) "fire_lfn_ext_up=", fire_lfn_ext_up
       write(*,*) "fire_test_steps=", fire_test_steps
       write(*,*) "fire_test_steps=", fire_test_steps
       write(*,*) "fire_topo_from_atm=", fire_topo_from_atm

      
       
     !  call fatal_error(h//" reading namelist")
    else
      print*,'isilda2'
       
       oneNamelistFile%ifire = ifire
       oneNamelistFile%crown = crown 
       oneNamelistFile%fire_fuel_read = fire_fuel_read
       oneNamelistFile%fire_fuel_cat = fire_fuel_cat
       oneNamelistFile%fire_num_ignitions = fire_num_ignitions     
       oneNamelistFile%fire_ignition_start_lon1 = fire_ignition_start_lon1
       oneNamelistFile%fire_ignition_start_lat1 = fire_ignition_start_lat1
       oneNamelistFile%fire_ignition_end_lon1 = fire_ignition_end_lon1
       oneNamelistFile%fire_ignition_end_lat1 = fire_ignition_end_lat1
       oneNamelistFile%fire_ignition_radius1 = fire_ignition_radius1
       oneNamelistFile%fire_ignition_ros1 = fire_ignition_ros1
       oneNamelistFile%fire_ignition_start_time1 = fire_ignition_start_time1
       oneNamelistFile%fire_ignition_end_time1 = fire_ignition_end_time1
       oneNamelistFile%fire_ignition_start_lon2 = fire_ignition_start_lon2
       oneNamelistFile%fire_ignition_start_lat2 = fire_ignition_start_lat2
       oneNamelistFile%fire_ignition_end_lon2 = fire_ignition_end_lon2 
       oneNamelistFile%fire_ignition_end_lat2 = fire_ignition_end_lat2
       oneNamelistFile%fire_ignition_radius2 = fire_ignition_radius2
       oneNamelistFile%fire_ignition_ros2 = fire_ignition_ros2
       oneNamelistFile%fire_ignition_start_time2 = fire_ignition_start_time2 
       oneNamelistFile%fire_ignition_end_time2 = fire_ignition_end_time2
       oneNamelistFile%fire_ignition_start_lon3 = fire_ignition_start_lon3
       oneNamelistFile%fire_ignition_start_lat3 = fire_ignition_start_lat3
       oneNamelistFile%fire_ignition_end_lon3 = fire_ignition_end_lon3 
       oneNamelistFile%fire_ignition_end_lat3 = fire_ignition_end_lat3 
       oneNamelistFile%fire_ignition_radius3 = fire_ignition_radius3
       oneNamelistFile%fire_ignition_ros3 = fire_ignition_ros3
       oneNamelistFile%fire_ignition_start_time3 = fire_ignition_start_time3
       oneNamelistFile%fire_ignition_end_time3 = fire_ignition_end_time3
       oneNamelistFile%fire_print_msg = fire_print_msg 
       oneNamelistFile%fire_print_file = fire_print_file   
       oneNamelistFile%fire_boundary_guard = fire_boundary_guard  
       oneNamelistFile%fire_fuel_left_method = fire_fuel_left_method 
       oneNamelistFile%fire_fuel_left_irl = fire_fuel_left_irl     
       oneNamelistFile%fire_fuel_left_jrl = fire_fuel_left_jrl    
       oneNamelistFile%fire_atm_feedback = fire_atm_feedback       
       oneNamelistFile%fire_grows_only = fire_grows_only      
       oneNamelistFile%fire_viscosity = fire_viscosity         
       oneNamelistFile%fire_upwinding = fire_upwinding      
       oneNamelistFile%fire_lfn_ext_up = fire_lfn_ext_up       
       oneNamelistFile%fire_test_steps = fire_test_steps      
       oneNamelistFile%fire_topo_from_atm = fire_topo_from_atm                

     end if

     close(iunit, iostat=err)
   ! if (err /= 0) then
    !   write(c0,"(i10)") err
     !  call fatal_error(h//" closing file "//&
      !      trim(oneNamelistFile%fileName)//" returned iostat="//&
       !     trim(adjustl(c0)))
    !end if

     return
     end subroutine ReadNamelistsfireFile

    subroutine DumpNamelistsfireFile(oneNamelistFile)
    implicit none
    integer, parameter :: ngrid= 1
   type(grid_config_rec_type), pointer :: oneNamelistFile
    real:: dx
    real:: dy
    real:: dt
    integer ::tracer_opt
    real ::cen_lat
    real ::cen_lon
    logical ::restart
    integer ::sr_x
    integer ::sr_y
    logical :: fmoist_run
    logical :: fmoist_interp
    logical :: fmoist_only
    integer :: fmoist_freq
    real :: fmoist_dt
    integer :: ifire
    logical :: crown
    integer :: fire_boundary_guard
    integer :: fire_num_ignitions
    real :: fire_ignition_ros1
    real :: fire_ignition_start_lon1
    real :: fire_ignition_start_lat1
    real :: fire_ignition_end_lon1
    real :: fire_ignition_end_lat1
    real :: fire_ignition_radius1
    real :: fire_ignition_start_time1
    real :: fire_ignition_end_time1
    real :: fire_ignition_ros2
    real :: fire_ignition_start_lon2
    real :: fire_ignition_start_lat2
    real :: fire_ignition_end_lon2
    real :: fire_ignition_end_lat2
    real :: fire_ignition_radius2
    real :: fire_ignition_start_time2
    real :: fire_ignition_end_time2
    real :: fire_ignition_ros3
    real :: fire_ignition_start_lon3
    real :: fire_ignition_start_lat3
    real :: fire_ignition_end_lon3
    real :: fire_ignition_end_lat3
    real :: fire_ignition_radius3
    real :: fire_ignition_start_time3
    real :: fire_ignition_end_time3
    real :: fire_ignition_ros4
    real :: fire_ignition_start_lon4
    real :: fire_ignition_start_lat4
    real :: fire_ignition_end_lon4 
    real :: fire_ignition_end_lat4
    real :: fire_ignition_radius4
    real :: fire_ignition_start_time4
    real :: fire_ignition_end_time4
    real :: fire_ignition_ros5 
    real :: fire_ignition_start_lon5
    real :: fire_ignition_start_lat5
    real :: fire_ignition_end_lon5
    real :: fire_ignition_end_lat5
    real :: fire_ignition_radius5
    real :: fire_ignition_start_time5
    real :: fire_ignition_end_time5
    real :: fire_ignition_start_x1
    real :: fire_ignition_start_y1
    real :: fire_ignition_end_x1
    real :: fire_ignition_end_y1
    real :: fire_ignition_start_x2
    real :: fire_ignition_start_y2
    real :: fire_ignition_end_x2
    real :: fire_ignition_end_y2
    real :: fire_ignition_start_x3
    real :: fire_ignition_start_y3
    real :: fire_ignition_end_x3
    real :: fire_ignition_end_y3
    real :: fire_ignition_start_x4
    real :: fire_ignition_start_y4
    real :: fire_ignition_end_x4
    real :: fire_ignition_end_y4
    real :: fire_ignition_start_x5
    real :: fire_ignition_start_y5
    real :: fire_ignition_end_x5
    real :: fire_ignition_end_y5
    real :: fire_perimeter_time
    real :: fire_lat_init
    real :: fire_lon_init
    real :: fire_ign_time
    integer ::fire_shape
    integer ::fire_sprd_mdl
    real ::fire_crwn_hgt
    real ::fire_ext_grnd
    real ::fire_ext_crwn
    integer ::fire_wind_log_interp
    integer ::fire_use_windrf
    integer ::fire_fuel_read
    integer ::fire_fmc_read
    integer ::fire_fuel_cat
    integer ::fire_print_msg
    integer ::fire_print_file
    logical ::fire_restart
    integer ::fire_time_step_ratio
    integer ::fire_debug_hook_sec
    integer ::fire_fuel_left_method
    integer ::fire_fuel_left_irl
    integer ::fire_fuel_left_jrl
    real :: fire_back_weight
    integer ::fire_grows_only
    integer ::fire_upwinding
    real :: fire_viscosity
    real :: fire_lfn_ext_up
    integer :: fire_topo_from_atm
    integer :: fire_advection
    integer :: fire_test_steps
    real :: fire_const_time
    real :: fire_const_grnhfx
    real :: fire_const_grnqfx
    integer :: fire_hfx_given
    integer :: fire_hfx_num_lines
    real :: fire_hfx_latent_part
    real :: fire_hfx_value1
    real :: fire_hfx_start_time1
    real :: fire_hfx_end_time1
    real :: fire_hfx_trans_time1
    real :: fire_hfx_radius1
    real :: fire_hfx_start_x1
    real :: fire_hfx_end_x1
    real :: fire_hfx_start_lat1
    real :: fire_hfx_end_lat1
    real :: fire_hfx_start_y1
    real :: fire_hfx_end_y1
    real :: fire_hfx_start_lon1
    real :: fire_hfx_end_lon1
    real :: fire_atm_feedback
    integer :: chem_opt
    integer :: nfmc
    integer :: fire_ignition_clamp
    integer :: fire_update_fuel_frac
    integer :: fndwi_from_ndwi
    integer :: kfmc_ndwi
    integer :: fire_can_top_read
    integer :: sfire_upwinding
    integer :: ng
    open(18, file='copy_RAMSIN_FIRE.dat', action='write', form='formatted')

    do ng = 1, ngrid

     write(18,*)'ng=',ng,'FIRE_DEFAULT'

       write(18,*) "dx=",oneNamelistFile%dx
       write(18,*) "dy=",oneNamelistFile%dy
       write(18,*) "dt=",oneNamelistFile%dt
       write(18,*) "tracer_opt=",oneNamelistFile%tracer_opt
       write(18,*) "cen_lat=",oneNamelistFile%cen_lat
       write(18,*) "cen_lon=",oneNamelistFile%cen_lon
       write(18,*) "restart=",oneNamelistFile%restart
       write(18,*) "sr_x=",oneNamelistFile%sr_x
       write(18,*) "sr_y=",oneNamelistFile%sr_y
       write(18,*) "fmoist_run=",oneNamelistFile%fmoist_run
       write(18,*) "fmoist_interp=",oneNamelistFile%fmoist_interp
       write(18,*) "fmoist_only=",oneNamelistFile%fmoist_only
       write(18,*) "fmoist_freq=",oneNamelistFile%fmoist_freq
       write(18,*) "fmoist_dt=", oneNamelistFile%fmoist_dt
       write(18,*) "fire_ignition_ros4=",oneNamelistFile%fire_ignition_ros4
       write(18,*) "fire_ignition_start_lon4=",&
      oneNamelistFile%fire_ignition_start_lon4
       write(18,*) "fire_ignition_start_lat4=",&
       oneNamelistFile%fire_ignition_start_lat4
       write(18,*) "fire_ignition_end_lon4=",&
       oneNamelistFile%fire_ignition_end_lon4
       write(18,*) "fire_ignition_end_lat4=",&
       oneNamelistFile%fire_ignition_end_lat4
       write(18,*) "fire_ignition_radius4=",&
       oneNamelistFile%fire_ignition_radius4
       write(18,*) "fire_ignition_start_time4=",&
       oneNamelistFile%fire_ignition_start_time4
       write(18,*) "fire_ignition_end_time4=",&
       oneNamelistFile%fire_ignition_end_time4
       write(18,*) "fire_ignition_ros5=",oneNamelistFile%fire_ignition_ros5
       write(18,*) "fire_ignition_start_lon5=",&
       oneNamelistFile%fire_ignition_start_lon5
       write(18,*) "fire_ignition_start_lat5=",&
       oneNamelistFile%fire_ignition_start_lat5
       write(18,*) "fire_ignition_end_lon5=",&
       oneNamelistFile%fire_ignition_end_lon5
       write(18,*) "fire_ignition_end_lat5=",&
       oneNamelistFile%fire_ignition_end_lat5
       write(18,*) "fire_ignition_radius5=",&
       oneNamelistFile%fire_ignition_radius5
       write(18,*) "fire_ignition_start_time5=",&
       oneNamelistFile%fire_ignition_start_time5
       write(18,*) "fire_ignition_end_time5=",&
       oneNamelistFile%fire_ignition_end_time5
       write(18,*) "fire_ignition_start_x1=", &
		oneNamelistFile%fire_ignition_start_x1
       write(18,*) "fire_ignition_start_y1=",&
      oneNamelistFile%fire_ignition_start_y1
       write(18,*) "fire_ignition_end_x1=",&
      oneNamelistFile%fire_ignition_end_x1
       write(18,*) "fire_ignition_end_y1=",&
      oneNamelistFile%fire_ignition_end_y1
       write(18,*) "fire_ignition_start_x2=",&
      oneNamelistFile%fire_ignition_start_x2
       write(18,*) "fire_ignition_start_y2=",&
      oneNamelistFile%fire_ignition_start_y2
       write(18,*) "fire_ignition_end_x2=",&
      oneNamelistFile%fire_ignition_end_x2
       write(18,*) "fire_ignition_end_y2=",&
      oneNamelistFile%fire_ignition_end_y2
       write(18,*) "fire_ignition_start_x3=",&
      oneNamelistFile%fire_ignition_start_x3
       write(18,*) "fire_ignition_start_y3=",&
      oneNamelistFile%fire_ignition_start_y3
       write(18,*) "fire_ignition_end_x3=",&
      oneNamelistFile%fire_ignition_end_x3
       write(18,*) "fire_ignition_end_y3=",&
      oneNamelistFile%fire_ignition_end_y3
       write(18,*) "fire_ignition_start_x4=",&
      oneNamelistFile%fire_ignition_start_x4  
       write(18,*) "fire_ignition_start_y4=", &
      oneNamelistFile%fire_ignition_start_y4
       write(18,*) "fire_ignition_end_x4=",&
      oneNamelistFile%fire_ignition_end_x4
       write(18,*) "fire_ignition_end_y4=",&
       oneNamelistFile%fire_ignition_end_y4
       write(18,*) "fire_ignition_start_x5="&
       ,oneNamelistFile%fire_ignition_start_x5
       write(18,*) "fire_ignition_start_y5=",&
      oneNamelistFile%fire_ignition_start_y5
       write(18,*) "fire_ignition_end_x5=",&
       oneNamelistFile%fire_ignition_end_x5
       write(18,*) "fire_ignition_end_y5=",&
       oneNamelistFile%fire_ignition_end_y5
       write(18,*) "fire_perimeter_time=",&
       oneNamelistFile%fire_perimeter_time
       write(18,*) "fire_lat_init=",oneNamelistFile%fire_lat_init
       write(18,*) "fire_lon_init=",oneNamelistFile%fire_lon_init
       write(18,*) "fire_ign_time=",oneNamelistFile%fire_ign_time
       write(18,*) "fire_shape=",oneNamelistFile%fire_shape
       write(18,*) "fire_sprd_mdl=",oneNamelistFile%fire_sprd_mdl
       write(18,*) "fire_crwn_hgt=",oneNamelistFile%fire_crwn_hgt
       write(18,*) "fire_ext_grnd=",oneNamelistFile%fire_ext_grnd
       write(18,*) "fire_ext_crwn=",oneNamelistFile%fire_ext_crwn
       write(18,*) "fire_wind_log_interp=",oneNamelistFile%fire_wind_log_interp
       write(18,*) "fire_use_windrf=",oneNamelistFile%fire_use_windrf
       write(18,*) "fire_fmc_read=",oneNamelistFile%fire_fmc_read
       write(18,*) "fire_restart=",oneNamelistFile%fire_restart
       write(18,*) "fire_time_step_ratio=",oneNamelistFile%fire_time_step_ratio
       write(18,*) "fire_debug_hook_sec=",oneNamelistFile%fire_debug_hook_sec
       write(18,*) "fire_back_weight=",oneNamelistFile%fire_back_weight
       write(18,*) "fire_advection=",oneNamelistFile%fire_advection
       write(18,*) "fire_const_time=",oneNamelistFile%fire_const_time
       write(18,*) "fire_const_grnhfx=",oneNamelistFile%fire_const_grnhfx
       write(18,*) "fire_const_grnqfx=",oneNamelistFile%fire_const_grnqfx
       write(18,*) "fire_hfx_given=",oneNamelistFile%fire_hfx_given 
       write(18,*) "fire_hfx_num_lines=",oneNamelistFile%fire_hfx_num_lines
       write(18,*) "fire_hfx_latent_part=",oneNamelistFile%fire_hfx_latent_part
       write(18,*) "fire_hfx_value1=",oneNamelistFile%fire_hfx_value1
       write(18,*) "fire_hfx_start_time1=", oneNamelistFile%fire_hfx_start_time1
       write(18,*) "fire_hfx_end_time1=",oneNamelistFile%fire_hfx_end_time1
       write(18,*) "fire_hfx_trans_time1=",oneNamelistFile%fire_hfx_trans_time1
       write(18,*) "fire_hfx_radius1=",oneNamelistFile%fire_hfx_radius1
       write(18,*) "fire_hfx_start_x1=",oneNamelistFile%fire_hfx_start_x1
       write(18,*) "fire_hfx_end_x1=",oneNamelistFile%fire_hfx_end_x1
       write(18,*) "fire_hfx_start_lat1=",oneNamelistFile%fire_hfx_start_lat1
       write(18,*) "fire_hfx_end_lat1=",oneNamelistFile%fire_hfx_end_lat1
       write(18,*) "fire_hfx_start_y1=",oneNamelistFile%fire_hfx_start_y1
       write(18,*) "fire_hfx_end_y1=",oneNamelistFile%fire_hfx_end_y1
       write(18,*) "fire_hfx_start_lon1=",oneNamelistFile%fire_hfx_start_lon1
       write(18,*) "fire_hfx_end_lon1=",oneNamelistFile%fire_hfx_end_lon1
       write(18,*) "chem_opt=",oneNamelistFile%chem_opt
       write(18,*) "nfmc=",oneNamelistFile%nfmc
       write(18,*) "fire_ignition_clamp=",oneNamelistFile%fire_ignition_clamp
       write(18,*)"fire_update_fuel_frac=",&
       oneNamelistFile%fire_update_fuel_frac
       write(18,*) "fndwi_from_ndwi=",oneNamelistFile%fndwi_from_ndwi
       write(18,*) "kfmc_ndwi=",oneNamelistFile%kfmc_ndwi
       write(18,*) "fire_can_top_read=",oneNamelistFile%fire_can_top_read
       write(18,*) "sfire_upwinding=",oneNamelistFile%sfire_upwinding




       write(18,*)'ng=',ng,'FIRE_SIMUL_INFO'

       write(18,*) "ifire=", oneNamelistFile%ifire
       write(18,*) "crown=", oneNamelistFile%crown
       write(18,*) "fire_fuel_read=", oneNamelistFile%fire_fuel_read
       write(18,*) "fire_fuel_cat=", oneNamelistFile%fire_fuel_cat
       write(18,*) "fire_num_ignitions=", oneNamelistFile%fire_num_ignitions
       write(18,*) "fire_ignition_start_lon1=", &
       oneNamelistFile%fire_ignition_start_lon1
       write(18,*) "fire_ignition_start_lat1=",&
       oneNamelistFile%fire_ignition_start_lat1
       write(18,*) "fire_ignition_end_lon1=",&
       oneNamelistFile%fire_ignition_end_lon1
       write(18,*) "fire_ignition_end_lat1=",&
       oneNamelistFile%fire_ignition_end_lat1
       write(18,*) "fire_ignition_radius1=",&
       oneNamelistFile%fire_ignition_radius1
       write(18,*) "fire_ignition_ros1=", oneNamelistFile%fire_ignition_ros1
       write(18,*) "fire_ignition_start_time1=",&
       oneNamelistFile%fire_ignition_start_time1
       write(18,*) "fire_ignition_end_time1=", oneNamelistFile%fire_ignition_end_time1
       write(18,*) "fire_ignition_start_lon2=",&
       oneNamelistFile%fire_ignition_start_lon2
       write(18,*) "fire_ignition_start_lat2=", &
       oneNamelistFile%fire_ignition_start_lat2
       write(18,*) "fire_ignition_end_lon2=", oneNamelistFile%fire_ignition_end_lon2
       write(18,*) "fire_ignition_end_lat2=", oneNamelistFile%fire_ignition_end_lat2
       write(18,*) "fire_ignition_radius2=", oneNamelistFile%fire_ignition_radius2
       write(18,*) "fire_ignition_ros2=", oneNamelistFile%fire_ignition_ros2 
       write(18,*) "fire_ignition_start_time2=", oneNamelistFile%fire_ignition_start_time2
       write(18,*) "fire_ignition_end_time2=", oneNamelistFile%fire_ignition_end_time2
       write(18,*) "fire_ignition_start_lon3=", oneNamelistFile%fire_ignition_start_lon3
       write(18,*) "fire_ignition_start_lat3=", oneNamelistFile%fire_ignition_start_lat3
       write(18,*) "fire_ignition_end_lon3=", oneNamelistFile%fire_ignition_end_lon3
       write(18,*) "fire_ignition_end_lat3=", oneNamelistFile%fire_ignition_end_lat3
       write(18,*) "fire_ignition_radius3=", oneNamelistFile%fire_ignition_radius3
       write(18,*) "fire_ignition_ros3=", oneNamelistFile%fire_ignition_ros3
       write(18,*) "fire_ignition_start_time3=", oneNamelistFile%fire_ignition_start_time3
       write(18,*) "fire_ignition_end_time3=", oneNamelistFile%fire_ignition_end_time3
       write(18,*) "fire_print_msg=", oneNamelistFile%fire_print_msg
       write(18,*) "fire_print_file=", oneNamelistFile%fire_print_file
       write(18,*) "fire_boundary_guard=", oneNamelistFile%fire_boundary_guard
       write(18,*) "fire_fuel_left_method=", oneNamelistFile%fire_fuel_left_method
       write(18,*) "fire_fuel_left_irl=", oneNamelistFile%fire_fuel_left_irl
       write(18,*) "fire_fuel_left_jrl=", oneNamelistFile%fire_fuel_left_jrl
       write(18,*) "fire_atm_feedback=", oneNamelistFile%fire_atm_feedback
       write(18,*) "fire_grows_only=", oneNamelistFile%fire_grows_only      
       write(18,*) "fire_viscosity=", oneNamelistFile%fire_viscosity
       write(18,*) "fire_upwinding=", oneNamelistFile%fire_upwinding
       write(18,*) "fire_lfn_ext_up=", oneNamelistFile%fire_lfn_ext_up
       write(18,*) "fire_test_steps=", oneNamelistFile%fire_test_steps
       write(18,*) "fire_test_steps=", oneNamelistFile%fire_test_steps
       write(18,*) "fire_topo_from_atm=",oneNamelistFile%fire_topo_from_atm
     enddo
     close(18)
     return
     END SUBROUTINE DumpNamelistsfireFile

    end module ModNamelistsfireFile
