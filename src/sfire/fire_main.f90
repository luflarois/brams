program fire_main
	use byteswap, only: SWAP_F4
	use module_domain_type
	use config_rec_type
	use module_fr_sfire_driver_brams
	USE module_tiles, ONLY : set_tiles
	USE module_machine
	use netcdf
	implicit none


	type(domain) :: grid
	type(grid_config_rec_type) :: config_flags
	
	INTEGER    ::  its, ite, jts, jte, kts, kte
	
	integer, parameter :: ids = 1
	integer, parameter :: ide = 43 
	integer, parameter :: kds = 1
	integer, parameter :: kde = 43
	integer, parameter :: jds = 1
	integer, parameter :: jde = 43
	integer, parameter :: ims = -4
	integer, parameter :: ime = 48
	integer, parameter :: kms = 1
	integer, parameter :: kme = 43
	integer, parameter :: jms = -4
	integer, parameter :: jme = 48           
        integer, parameter :: ips = 1
	integer, parameter :: ipe = 43
	integer, parameter :: kps = 1
	integer, parameter :: kpe = 43
	integer, parameter :: jps = 1
	integer, parameter :: jpe = 43
	
	integer, parameter :: ifds =1
    	integer, parameter :: ifde = 430
    	integer, parameter :: kfds = 1
    	integer, parameter :: kfde = 43
    	integer, parameter :: jfds = 1
    	integer, parameter :: jfde = 430
    	integer, parameter :: ifms = -49
    	integer, parameter :: ifme = 480
    	integer, parameter :: kfms = 1 
    	integer, parameter :: kfme = 43
    	integer, parameter :: jfms = -49
    	integer, parameter :: jfme = 480
    	integer, parameter :: ifps = 1
    	integer, parameter :: ifpe = 430
    	integer, parameter :: kfps = 1
    	integer, parameter :: kfpe = 43
    	integer, parameter :: jfps = 1
    	integer, parameter :: jfpe = 430
	
	integer, parameter :: borders = 5
	
	!real, parameter :: g = 9.8 !m/s
	
	integer :: nz
	integer :: ntime
	integer :: nlat
	integer :: nlon 
	integer :: nlatSubgrid
	integer :: nlonSubgrid
	!netcdf
	
	integer :: ncid
	integer :: varid
	
	integer :: i
	integer :: j
	integer :: k
	integer :: t
	integer :: nrec
	
	integer, parameter :: num_tiles = 1
	
	
	
	character (len=255) :: filename, tc 
	
	real, allocatable, dimension(:,:,:)  :: xlatFromFile
	real, allocatable, dimension(:,:,:)  :: xlonFromFile
	
	real, allocatable, dimension(:,:,:,:)  :: phbFromFile
	real, allocatable, dimension(:,:,:,:)  :: ph_2FromFile
	
	real, allocatable, dimension(:,:,:)  :: mutFromFile
	
	
	real,dimension(ims:ime, kms:kme, jms:jme)::rho      !densidade do ar
	real,dimension(ims:ime, kms:kme, jms:jme)::dz8w     !altura da camada (dz across w-lvl)

	!****grid%
	real,dimension(ims:ime, jms:jme) :: xlong
	real,dimension(ims:ime, jms:jme) :: xlat            ! inout because of extension at bdry 
	
	!****grid%
	real,allocatable, dimension(:,:,:) :: fxlong
	real,allocatable, dimension(:,:,:) :: fxlat            ! inout because of extension at bdry 
	
	integer :: i_start
	integer ::i_end
	integer :: k_start
	integer ::k_end
	integer ::j_start
	integer ::j_end  ! atm grid tiling 
	
	real:: dx
	
	
	integer :: tbyte              
	!real,dimension(ims:ime,kms:kme,jms:jme)::u_2
	!real,dimension(ims:ime,kms:kme,jms:jme)::v_2   ! wind velocity (m/s) (staggered atm grid) 
	!real,dimension(ims:ime,kms:kme,jms:jme)::ph
	!real,dimension(ims:ime,kms:kme,jms:jme)::ph_2
	!real,dimension(ims:ime,kms:kme,jms:jme)::phb                      ! geopotential (w-points atm grid) 
	!real,dimension(ims:ime, jms:jme)::z0                ! roughness height 
	!real,dimension(ims:ime, jms:jme)::zs                ! terrain height  
	!real,dimension(ims:ime,jms:jme)::rain_old      ! rain_old(i,j) = rainc(i,j) + rainnc(i,j )
	!real,dimension(ims:ime,jms:jme):: t2
	!real,dimension(ims:ime,jms:jme):: q2
	!real,dimension(ims:ime,jms:jme):: psfc
	!real,dimension(ims:ime,jms:jme):: rainc
	!real,dimension(ims:ime,jms:jme):: rainnc       !temperature (K), vapor contents (kg/kg), pr
	!real,dimension(ims:ime,kms:kme,jms:jme):: fmc_gc       !temperature (K), vapor contents (kg/kg), pr
	
	
	
	INQUIRE (IOLENGTH=tbyte) dx
	
	kts = kps ; kte = kpe     
   	its = ips ; ite = ipe     
   	jts = jps ; jte = jpe 
   
	i_start = its
      	i_end   = min(ite,ide-1)
      	j_start = jts
      	j_end   = min(jte,jde-1)
      	k_start = kts
      	k_end   = kte-1
	
	config_flags%dx = 200
config_flags%dy = 200
config_flags%dt = 2.
config_flags%tracer_opt = 0
config_flags%cen_lat = 39.7053
config_flags%cen_lon = -107.2907
config_flags%restart = .false.
config_flags%sr_x = 0
config_flags%sr_y = 0
config_flags%fmoist_run = .false.
config_flags%fmoist_interp = .false.
config_flags%fmoist_only = .false.
config_flags%fmoist_freq = 0
config_flags%fmoist_dt = 600
config_flags%ifire = 0
config_flags%fire_boundary_guard = 2
config_flags%fire_num_ignitions = 0.
config_flags%fire_ignition_ros1 = 0.01
config_flags%fire_ignition_start_lon1 = 0.
config_flags%fire_ignition_start_lat1 = 0.
config_flags%fire_ignition_end_lon1 = 0.
config_flags%fire_ignition_end_lat1 = 0.
config_flags%fire_ignition_radius1 = 400.
config_flags%fire_ignition_start_time1 = 0.
config_flags%fire_ignition_end_time1 = 0.
config_flags%fire_ignition_ros2 = 0.01
config_flags%fire_ignition_start_lon2 = 0.
config_flags%fire_ignition_start_lat2 = 0.
config_flags%fire_ignition_end_lon2 = 0.
config_flags%fire_ignition_end_lat2 = 0.
config_flags%fire_ignition_radius2 = 400.
config_flags%fire_ignition_start_time2 = 0.
config_flags%fire_ignition_end_time2 = 0.
config_flags%fire_ignition_ros3 = 0.01
config_flags%fire_ignition_start_lon3 = 0.
config_flags%fire_ignition_start_lat3 = 0.
config_flags%fire_ignition_end_lon3 = 0.
config_flags%fire_ignition_end_lat3 = 0.
config_flags%fire_ignition_radius3 = 600.
config_flags%fire_ignition_start_time3 = 0.
config_flags%fire_ignition_end_time3 = 0.
config_flags%fire_ignition_ros4 = 0.01
config_flags%fire_ignition_start_lon4 = 0.
config_flags%fire_ignition_start_lat4 = 0.
config_flags%fire_ignition_end_lon4 = 0.
config_flags%fire_ignition_end_lat4 = 0.
config_flags%fire_ignition_radius4 = 0.
config_flags%fire_ignition_start_time4 = 0.
config_flags%fire_ignition_end_time4 = 0.
config_flags%fire_ignition_ros5 = 0.01
config_flags%fire_ignition_start_lon5 = 0.
config_flags%fire_ignition_start_lat5 = 0.
config_flags%fire_ignition_end_lon5 = 0.
config_flags%fire_ignition_end_lat5 = 0.
config_flags%fire_ignition_radius5 = 0.
config_flags%fire_ignition_start_time5 = 0.
config_flags%fire_ignition_end_time5 = 0.
config_flags%fire_ignition_start_x1 = 0.
config_flags%fire_ignition_start_y1 = 0.
config_flags%fire_ignition_end_x1 = 0.
config_flags%fire_ignition_end_y1 = 0.
config_flags%fire_ignition_start_x2 = 0.
config_flags%fire_ignition_start_y2 = 0.
config_flags%fire_ignition_end_x2 = 0.
config_flags%fire_ignition_end_y2 = 0.
config_flags%fire_ignition_start_x3 = 0.
config_flags%fire_ignition_start_y3 = 0.
config_flags%fire_ignition_end_x3 = 0.
config_flags%fire_ignition_end_y3 = 0.
config_flags%fire_ignition_start_x4 = 0.
config_flags%fire_ignition_start_y4 = 0.
config_flags%fire_ignition_end_x4 = 0.
config_flags%fire_ignition_end_y4 = 0.
config_flags%fire_ignition_start_x5 = 0.
config_flags%fire_ignition_start_y5 = 0.
config_flags%fire_ignition_end_x5 = 0.
config_flags%fire_ignition_end_y5 = 0.
config_flags%fire_perimeter_time = 0.
config_flags%fire_lat_init = 0.
config_flags%fire_lon_init = 0.
config_flags%fire_ign_time = 0.
config_flags%fire_shape = 0
config_flags%fire_sprd_mdl = 1
config_flags%fire_crwn_hgt = 15.
config_flags%fire_ext_grnd = 50.
config_flags%fire_ext_crwn = 50.
config_flags%fire_wind_log_interp = 4
config_flags%fire_use_windrf = 0
config_flags%fire_fuel_read = -1
config_flags%fire_fmc_read = 1
config_flags%fire_fuel_cat = 1
config_flags%fire_print_msg = 0
config_flags%fire_print_file = 0
config_flags%fire_restart = .false.
config_flags%fire_time_step_ratio = 1
config_flags%fire_debug_hook_sec = 0
config_flags%fire_fuel_left_method = 1
config_flags%fire_fuel_left_irl = 2
config_flags%fire_fuel_left_jrl = 2
config_flags%fire_back_weight = 0.5
config_flags%fire_grows_only = 1
config_flags%fire_upwinding = 3
config_flags%fire_viscosity = 0.4
config_flags%fire_lfn_ext_up = 1.0
config_flags%fire_topo_from_atm = 1
config_flags%fire_advection = 1
config_flags%fire_test_steps = 0
config_flags%fire_const_time = -1.
config_flags%fire_const_grnhfx = 0.
config_flags%fire_const_grnqfx = 0.
config_flags%fire_hfx_given = 0
config_flags%fire_hfx_num_lines = 0
config_flags%fire_hfx_latent_part = 0.084
config_flags%fire_hfx_value1 = 0.
config_flags%fire_hfx_start_time1 = 0.
config_flags%fire_hfx_end_time1 = 0.
config_flags%fire_hfx_trans_time1 = 0.
config_flags%fire_hfx_radius1 = 0.
config_flags%fire_hfx_start_x1 = 0.
config_flags%fire_hfx_end_x1 = 0.
config_flags%fire_hfx_start_lat1 = 0.
config_flags%fire_hfx_end_lat1 = 0.
config_flags%fire_hfx_start_y1 = 0.
config_flags%fire_hfx_end_y1 = 0.
config_flags%fire_hfx_start_lon1 = 0.
config_flags%fire_hfx_end_lon1 = 0.
config_flags%fire_atm_feedback = 1.
config_flags%chem_opt = 0

	
	
	
	! START - config flags from namelist
	
	config_flags%run_days                            = 0
 	config_flags%run_hours                           = 2
 config_flags%run_minutes                         = 0
 config_flags%run_seconds                         = 0
 config_flags%start_year                          = 2013
 config_flags%start_month                         = 05
 config_flags%start_day                           = 10
 config_flags%start_hour                          = 00
 config_flags%start_minute                        = 00
 config_flags%start_second                        = 00
 config_flags%end_year                            = 2013
 config_flags%end_month                           = 05
 config_flags%end_day                             = 12
 config_flags%end_hour                            = 00
 config_flags%end_minute                          = 00
 config_flags%end_second                          = 00
 config_flags%interval_seconds                    = 21600
 config_flags%input_from_file                     = .true.
 config_flags%history_interval_s                  = 30
 config_flags%frames_per_outfile                  = 1000
 config_flags%restart                             = .false.
 config_flags%restart_interval                    = 1
 config_flags%io_form_history                     = 2
 config_flags%io_form_restart                     = 2
 config_flags%io_form_input                       = 2
 config_flags%io_form_boundary                    = 2
 config_flags%debug_level                         = 0
 
 config_flags%time_step                           = 0
 config_flags%time_step_fract_num                 = 5
 config_flags%time_step_fract_den                 = 10
 config_flags%max_dom                             = 1
 config_flags%s_we                                = 1
 config_flags%e_we                                = 43
 config_flags%s_sn                                = 1
 config_flags%e_sn                                = 43
 config_flags%s_vert                              = 1
 config_flags%e_vert                              = 43
 config_flags%num_metgrid_levels                  = 27
 config_flags%num_metgrid_soil_levels             = 4
 config_flags%dx                                  = 60
 config_flags%dy                                  = 60
 config_flags%grid_id                             = 1
 config_flags%parent_id                           = 0
 config_flags%i_parent_start                      = 0
 config_flags%j_parent_start                      = 0
 config_flags%parent_grid_ratio                   = 1
 config_flags%parent_time_step_ratio              = 1
 config_flags%feedback                            = 1
 config_flags%smooth_option                       = 0
 config_flags%sr_x                                = 10
 config_flags%sr_y                                = 10
 config_flags%sfcp_to_sfcp                        = .true.
 config_flags%p_top_requested                     = 10000
 
 config_flags%mp_physics                          = 2
 config_flags%ra_lw_physics                       = 1
 config_flags%ra_sw_physics                       = 1
 config_flags%radt                                = 10
 config_flags%sf_sfclay_physics                   = 1
 config_flags%sf_surface_physics                  = 1
 config_flags%bl_pbl_physics                      = 1
 config_flags%bldt                                = 0
 config_flags%cu_physics                          = 1
 config_flags%cudt                                = 5
 config_flags%isfflx                              = 1
 config_flags%ifsnow                              = 0
 config_flags%icloud                              = 0
 config_flags%surface_input_source                = 1
 config_flags%num_soil_layers                     = 5
 config_flags%sf_urban_physics                    = 0
 config_flags%maxiens                             = 1
 config_flags%maxens                              = 3
 config_flags%maxens2                             = 3
 config_flags%maxens3                             = 16
 config_flags%ensdim                              = 144
 
 config_flags%scm_th_adv                          = .true.
 config_flags%scm_wind_adv                        = .true.
 config_flags%scm_qv_adv                          = .true.
 config_flags%scm_ql_adv                          = .true.
 config_flags%scm_vert_adv                        = .true.
 
 config_flags%rk_ord                              = 3
 config_flags%w_damping                           = 0
 config_flags%diff_opt                            = 2
 config_flags%km_opt                              = 2
 config_flags%damp_opt                            = 0
 config_flags%base_temp                           = 290.
 config_flags%zdamp                               = 5000.
 config_flags%dampcoef                            = 0.2
 config_flags%khdif                               = 0.05
 config_flags%kvdif                               = 0.05
 config_flags%smdiv                               = 0.1
 config_flags%emdiv                               = 0.01
 config_flags%epssm                               = 0.1      
 config_flags%pert_coriolis                       = .false.
 config_flags%tke_drag_coefficient                = 0.
 config_flags%tke_heat_flux                       = 0.
 config_flags%time_step_sound                     = 20
 config_flags%h_mom_adv_order                     = 5
 config_flags%v_mom_adv_order                     = 3
 config_flags%h_sca_adv_order                     = 5
 config_flags%v_sca_adv_order                     = 3
 config_flags%non_hydrostatic                     = .true.
 config_flags%momentum_adv_opt                    = 1
 config_flags%moist_adv_opt                       = 1
 config_flags%scalar_adv_opt                      = 1
 config_flags%chem_adv_opt                        = 1     
 config_flags%tke_adv_opt                         = 1  
 

 
 config_flags%ifire              = 2    

 config_flags%fire_fuel_read     = -1   
 config_flags%fire_fuel_cat      = 3   
 
 config_flags%nfmc = 5
 
 config_flags%dt = 2.


 config_flags%fire_num_ignitions = 3       
 config_flags%fire_ignition_start_lon1=-107.293664
 config_flags%fire_ignition_start_lat1 =  39.698696
 config_flags%fire_ignition_end_lon1 = -107.293664 
 config_flags%fire_ignition_end_lat1 =    39.710990 
 config_flags%fire_ignition_radius1 =    370. 
 config_flags%fire_ignition_start_time1  =      2
 config_flags%fire_ignition_end_time1  =      2 
 config_flags%fire_ignition_start_lon2=-107.287954  
 config_flags%fire_ignition_start_lat2 =  39.698696 
 config_flags%fire_ignition_end_lon2 = -107.287954 
 config_flags%fire_ignition_end_lat2 =    39.71099 
 config_flags%fire_ignition_radius2 =    370. 
 config_flags%fire_ignition_start_time2  =      3 
 config_flags%fire_ignition_end_time2  =      3 
 config_flags%fire_ignition_start_lon3=-107.289096 
 config_flags%fire_ignition_start_lat3 =  39.706599 
 config_flags%fire_ignition_end_lon3 =   -107.289096  
 config_flags%fire_ignition_end_lat3 =    39.706599  
 config_flags%fire_ignition_radius3 =    400. 
 config_flags%fire_ignition_start_time3  =      4
 config_flags%fire_ignition_end_time3  =      4 

 config_flags%fire_print_msg     = 1    
 config_flags%fire_print_file    = 1   
 
 config_flags%fire_wind_log_interp    = 4
 config_flags%fire_use_windrf    = 0      

 config_flags%fire_boundary_guard = -1      
 config_flags%fire_fuel_left_method=1      
 config_flags%fire_fuel_left_irl=2          
 config_flags%fire_fuel_left_jrl=2          
 config_flags%fire_atm_feedback=1.          
 config_flags%fire_grows_only=1             
 config_flags%fire_viscosity=0.4            
 config_flags%fire_upwinding=3             
 config_flags%fire_lfn_ext_up=1.0          
 config_flags%fire_test_steps=0             
 config_flags%fire_topo_from_atm=0         
 
 config_flags%spec_bdy_width                      = 5
 config_flags%spec_zone                           = 1
 config_flags%relax_zone                          = 4
 config_flags%specified                           = .true.
 config_flags%periodic_x                          = .false.
 config_flags%symmetric_xs                        = .false.
 config_flags%symmetric_xe                        = .false.
 config_flags%open_xs                             = .false.
 config_flags%open_xe                             = .false.
 config_flags%periodic_y                          = .false.
 config_flags%symmetric_ys                        = .false.
 config_flags%symmetric_ye                        = .false.
 config_flags%open_ys                             = .false.
 config_flags%open_ye                             = .false.
 config_flags%nested                              = .false.
 
 
 config_flags%fmoist_run = .false.
 config_flags%fmoist_interp = .false.
	

 !Ramsin default
 
  config_flags%tracer_opt                         = 0
  config_flags%cen_lat                         = 39.7053
  config_flags%cen_lon                         = -107.2907	
  config_flags%fmoist_only = .false.
  config_flags%fmoist_freq = 0
  config_flags%fmoist_dt = 600	
  
  config_flags%fire_ignition_ros1 = 10.
  config_flags%fire_ignition_ros2 = 10.
  config_flags%fire_ignition_ros3 = 10.
  
  
  
  
	
	! END - config flags from namelist
	
	
grid%time_step                           = 0
 grid%time_step_fract_num                 = 5
 grid%time_step_fract_den                 = 10
 grid%max_dom                             = 1
 grid%s_we                                = 1
 grid%e_we                                = 43
 grid%s_sn                                = 1
 grid%e_sn                                = 43
 grid%s_vert                              = 1
 grid%e_vert                              = 43

 grid%dx                                  = 60
 grid%dy                                  = 60
 grid%grid_id                             = 1
 grid%parent_id                           = 0
 grid%i_parent_start                      = 0
 grid%j_parent_start                      = 0
 grid%parent_grid_ratio                   = 1
 grid%parent_time_step_ratio              = 1


 grid%sr_x                                = 10
 grid%sr_y                                = 10

	
 grid%fire_time_step_ratio                                = 2	
 
  grid%u_frame = 0.
 grid%v_frame = 0.
 !grid%uf
 !grid%vf
 !grid%zsf
 !grid%dzdxf
 !grid%dzdyf
 !grid%bbb
 !grid%phisc
 !grid%phiwc
 !grid%r_0
 !grid%fgip
 !grid%ischap
 !grid%fuel_time
 !grid%fmc_g
 !grid%nfuel_cat
 !
 config_flags%fire_time_step_ratio = 2	
	
	allocate(grid%xlong(ims:ime,jms:jme))
	allocate(grid%xlat(ims:ime,jms:jme))
	
	allocate(grid%u_2(ims:ime,kms:kme,jms:jme))
	allocate(grid%v_2(ims:ime,kms:kme,jms:jme))
	allocate(grid%ph_2(ims:ime,kms:kme,jms:jme))
	allocate(grid%phb(ims:ime,kms:kme,jms:jme))
	allocate(grid%z0(ims:ime,jms:jme))
	
	allocate(grid%ht(ims:ime,jms:jme))
	allocate(grid%rain_old(ims:ime,jms:jme))
	allocate(grid%t2(ims:ime,jms:jme))
	allocate(grid%q2(ims:ime,jms:jme))
	allocate(grid%psfc(ims:ime,jms:jme))
	allocate(grid%rainc(ims:ime,jms:jme))
	allocate(grid%rainnc(ims:ime,jms:jme))
	
	allocate(grid%t2_old(ims:ime,jms:jme))
	allocate(grid%q2_old(ims:ime,jms:jme))
	allocate(grid%psfc_old(ims:ime,jms:jme))
	
	allocate(grid%rh_fire(ims:ime,jms:jme))
	
	allocate(grid%fmc_gc(ims:ime,1:config_flags%nfmc,jms:jme))
	allocate(grid%fmc_equi(ims:ime,1:config_flags%nfmc,jms:jme))
	allocate(grid%fmc_lag(ims:ime,1:config_flags%nfmc,jms:jme))
	
	allocate(grid%lfn(ifms:ifme,jfms:jfme))
	allocate(grid%tign_g(ifms:ifme,jfms:jfme))
	allocate(grid%fire_area(ifms:ifme,jfms:jfme))
	allocate(grid%fuel_frac(ifms:ifme,jfms:jfme))
	
	allocate(grid%zsf(ifms:ifme,jfms:jfme))
	allocate(grid%dzdxf(ifms:ifme,jfms:jfme))
	allocate(grid%dzdyf(ifms:ifme,jfms:jfme))
	allocate(grid%fmc_g(ifms:ifme,jfms:jfme))
	allocate(grid%nfuel_cat(ifms:ifme,jfms:jfme))
	allocate(grid%fuel_time(ifms:ifme,jfms:jfme))
	allocate(grid%fuel_frac_burnt(ifms:ifme,jfms:jfme))
	
	allocate(grid%uf(ifms:ifme,jfms:jfme))
	allocate(grid%vf(ifms:ifme,jfms:jfme))
	allocate(grid%bbb(ifms:ifme,jfms:jfme))
	allocate(grid%phisc(ifms:ifme,jfms:jfme))
	allocate(grid%phiwc(ifms:ifme,jfms:jfme))
	allocate(grid%r_0(ifms:ifme,jfms:jfme))
	allocate(grid%fgip(ifms:ifme,jfms:jfme))
	allocate(grid%ischap(ifms:ifme,jfms:jfme))

        allocate(grid%phisc_FM10(ifms:ifme,jfms:jfme))!INTRODUZIDO POR ISILDA MENEZES
        allocate(grid%phiwc_FM10(ifms:ifme,jfms:jfme))!INTRODUZIDO POR ISILDA MENEZES
        allocate(grid%r_0_FM10(ifms:ifme,jfms:jfme))!INTRODUZIDO POR ISILDA MENEZES
        allocate(grid%bbb_FM10(ifms:ifme,jfms:jfme))!INTRODUZIDO POR ISILDA MENEZES
        allocate(grid%CFB(ifms:ifme,jfms:jfme))!INTRODUZIDO POR ISILDA MENEZES
        allocate(grid%HPA(ifms:ifme,jfms:jfme))!INTRODUZIDO POR ISILDA MENEZES

	allocate(grid%fz0(ifms:ifme,jfms:jfme))
	allocate(grid%fwh(ifms:ifme,jfms:jfme))
	
	allocate(grid%avg_fuel_frac(ims:ime, jms:jme))
	allocate(grid%grnhfx(ims:ime, jms:jme))
	allocate(grid%grnqfx(ims:ime, jms:jme))
	allocate(grid%canhfx(ims:ime, jms:jme))
	allocate(grid%canqfx(ims:ime, jms:jme))
	allocate(grid%uah(ims:ime,jms:jme))
	allocate(grid%vah(ims:ime,jms:jme))
	
	allocate(grid%fgrnhfx(ifms:ifme, jfms:jfme))
	allocate(grid%fgrnqfx(ifms:ifme, jfms:jfme))
	allocate(grid%fcanhfx(ifms:ifme, jfms:jfme))
	allocate(grid%fcanqfx(ifms:ifme, jfms:jfme))
	allocate(grid%ros(ifms:ifme,jfms:jfme))
	allocate(grid%flineint(ifms:ifme,jfms:jfme))
	
	allocate(grid%flineint2(ifms:ifme,jfms:jfme))
	allocate(grid%f_ros0(ifms:ifme,jfms:jfme))
	allocate(grid%f_rosx(ifms:ifme,jfms:jfme))
	allocate(grid%f_rosy(ifms:ifme,jfms:jfme))
	allocate(grid%f_ros(ifms:ifme,jfms:jfme))
	allocate(grid%f_int(ifms:ifme,jfms:jfme))
	
	allocate(grid%f_lineint(ifms:ifme,jfms:jfme))
	allocate(grid%f_lineint2(ifms:ifme,jfms:jfme))
	
	
	allocate(grid%fxlong(ifms:ifme, jfms:jfme))
	allocate(grid%fxlat(ifms:ifme, jfms:jfme))
	allocate(grid%fire_hfx(ifms:ifme, jfms:jfme))
	
	allocate(grid%z_at_w(ims:ime, kms:kme, jms:jme))
	
	allocate(grid%mut(ims:ime, jms:jme))
	
	allocate(grid%rthfrten(ims:ime, kms:kme, jms:jme))
	allocate(grid%rqvfrten(ims:ime, kms:kme, jms:jme))

	
	
	grid%num_tiles = num_tiles
	grid%num_tiles_spec = num_tiles
	
	
	
	CALL init_module_machine()
	
	
	CALL set_tiles ( grid , ids , ide , jds , jde , ips , ipe , jps , jpe )


	! Read namelist fire with grid parameters 
	! (Time, south_north, west_east) ;
	! 
	
	filename = '/scratchin/grupos/catt-brams/home/isilda.menezes/wrf-fire/wrfv2_fire/run/wrfout_d01_2013-05-10_00:00:00'

	
! Read netcdf file for xlat, xlon ...
! 

	call check( nf90_open(trim(filename), NF90_NOWRITE, ncid) )
	
	call check( nf90_inq_dimid(ncid, 'south_north', varid) )
	call check( nf90_inquire_dimension(ncid, varid, len = nlat) )
	call check( nf90_inq_dimid(ncid, 'south_north_subgrid', varid) )
	call check( nf90_inquire_dimension(ncid, varid, len = nlatSubGrid) )
	call check( nf90_inq_dimid(ncid, 'west_east', varid) )
	call check( nf90_inquire_dimension(ncid, varid, len = nlon) )
	call check( nf90_inq_dimid(ncid, 'west_east_subgrid', varid) )
	call check( nf90_inquire_dimension(ncid, varid, len = nlonSubgrid) )
	call check( nf90_inq_dimid(ncid, 'Time', varid) )
	call check( nf90_inquire_dimension(ncid, varid, len = ntime) )	
	call check( nf90_inq_dimid(ncid, 'bottom_top_stag', varid) )
	call check( nf90_inquire_dimension(ncid, varid, len = nz) )	
	
	allocate(xlatFromFile(nlon, nlat, 1), &
		xlonFromFile(nlon, nlat, 1), &
		fxlat(nlonSubgrid, nlatSubgrid, 1), &
		fxlong(nlonSubgrid, nlatSubgrid, 1), &
		phbFromFile(nlon, nlat, nz, 1), &
		mutFromFile(nlon, nlat, 1))
	
	
	call check( nf90_inq_varid(ncid, 'XLAT', varid) )
	call check( nf90_get_var(ncid, varid, xlatFromFile, start = (/ 1, 1, 241 /)) )
	
	call check( nf90_inq_varid(ncid, 'XLONG', varid) )
	call check( nf90_get_var(ncid, varid, xlonFromFile, start = (/ 1, 1, 241 /)) )
	
	call check( nf90_inq_varid(ncid, 'FXLAT', varid) )
	call check( nf90_get_var(ncid, varid, fxlat, start = (/ 1, 1, 241 /)) )
	
	call check( nf90_inq_varid(ncid, 'FXLONG', varid) )
	call check( nf90_get_var(ncid, varid, fxlong, start = (/ 1, 1, 241 /)) )
	
	
	
	!call check( nf90_inq_varid(ncid, 'PHB', varid) )
	!call check( nf90_get_var(ncid, varid, phb, start = (/ 1, 1, 241 /)) )
	
	!call check( nf90_close(ncid) )
	
	do j=jps,max(jpe,jde-1)
            do i=ips,max(ipe,ide-1)
                 grid%xlong(i,j)=xlonFromFile(i, j, 1)
		 grid%xlat(i,j)=xlatFromFile(i, j, 1)
	    enddo
        enddo
	
	do j=jfps,max(jfpe,jfde-1)
            do i=ifps,max(ifpe,ifde-1)
                 grid%fxlong(i,j)=fxlong(i, j, 1)
		 grid%fxlat(i,j)=fxlat(i, j, 1)
	    enddo
        enddo

!if(config_flags%ifire.eq.2)then

 
   
!endif	
	
	!do j=jts,min(jte,jde-1)
        !    do i=its,min(ite,ide-1)
        !       do k=kts,kte
        !         phb(i,k,j)=phbFromFile(i, j, k, 1)
        !       enddo
	!    enddo
        !enddo
	
	!do j=jts,min(jte,jde-1)
        !    do i=its,min(ite,ide-1)
        !       do k=kts,kte
        !         ph_2(i,k,j)=ph_2FromFile(i, j, k, 1)
        !       enddo
	!    enddo
        !enddo
	
	!do j=jts,jte
        !    do i=its,ite
        !        grid%ph_2(i,1,j) = 0.
        !    enddo
        !enddo
	
	!do j=jts,min(jte,jde-1)
        !    do i=its,min(ite,ide-1)
        !       do k=kts,kte
        !         z_at_w(i,k,j)=(grid%ph_2(i,k,j)+grid%phb(i,k,j))/g
        !       enddo
	!    enddo
        !enddo
	
	!grid%ph_2(i,1,j) = grid%phb(i,1,j)
        !do k = 2,kte
        !	pfu = grid%mu0(i,j)*grid%znw(k)   + grid%p_top
        !        pfd = grid%mu0(i,j)*grid%znw(k-1) + grid%p_top
        !        phm = grid%mu0(i,j)*grid%znu(k-1) + grid%p_top
        !        grid%ph_2(i,k,j) = grid%ph_2(i,k-1,j) + grid%alt(i,k-1,j)*phm*LOG(pfd/pfu)
        !END DO

        !DO k = 1,kte
        !	grid%ph_2(i,k,j) = grid%ph_2(i,k,j) - grid%phb(i,k,j)
        !END DO
	
        !tbyte = 1
	
	print*, tbyte, (ime+borders)*kme*(jme+borders)*tbyte, ime, borders, kme, jme
!*isilda*****
         open(unit=13,file="/stornext/home/rafael.stockler/projetos/BRAMS/data/u_2.dat",&
	 	action="read",form="unformatted",recl=(ime+borders)*kme*(jme+borders)*tbyte,access="direct")
         read(13, rec=1) (((grid%u_2(i,k,j), i=ims,ime),k=kms,kme),j=jms,jme)
         close(unit=13)
	 
	 do j=1, jme
	 	do k=1, kme
	 		do i=1, ime
	 			call SWAP_F4(grid%u_2(i,k,j))
			enddo
		enddo
	 enddo
	 
	 !where( isnan(grid%u_2)) grid%u_2 = -999.
	 
	  
	 
!*isilda*****
         open(unit=13,file="/stornext/home/rafael.stockler/projetos/BRAMS/data/v_2.dat",&
	 	action="read",form="unformatted",recl=(ime+borders)*kme*(jme+borders)*tbyte,access="direct")
         read(13, rec=1) (((grid%v_2(i,k,j), i=ims,ime),k=kms,kme),j=jms,jme)
         close(unit=13)
	 
	 do j=1, jme
	 	do k=1, kme
	 		do i=1, ime
	 			call SWAP_F4(grid%v_2(i,k,j))
			enddo
		enddo
	 enddo
	 
	 !where( isnan(grid%v_2)) grid%v_2 = -999.
	 
	 
!*isilda******
!*isilda*****
         open(unit=13,file="/stornext/home/rafael.stockler/projetos/BRAMS/data/phb.dat",&
	 	action="read",form="unformatted",recl=(ime+borders)*kme*(jme+borders)*tbyte,access="direct")
         read(13, rec=1) (((grid%phb(i,k,j), i=ims,ime),k=kms,kme),j=jms,jme)
         close(unit=13)
	 
	 do j=1, jme
	 	do k=1, kme
	 		do i=1, ime
	 			call SWAP_F4(grid%phb(i,k,j))
			enddo
		enddo
	 enddo
	 
	 !where( isnan(grid%phb)) grid%phb = -999.
	 
	 
!*isilda******
!*isilda*****
!         open(13,file='z_at_w.dat',action="write",form='unformatted')
 !        write(13) z_at_w
  !       close(13)
!*isilda******

!*isilda*****
         open(unit=13,file="/stornext/home/rafael.stockler/projetos/BRAMS/data/t2.dat",&
	 	action="read",form="unformatted",recl=(ime+borders)*(jme+borders)*tbyte,access="direct")
         read(13, rec=1) ((grid%t2(i,j), i=ims,ime),j=jms,jme)
         close(unit=13)
	 
	 do j=1, jme
	 	do i=1, ime
	 		call SWAP_F4(grid%t2(i,j))
		enddo
	 enddo
	 
	 
	 
	 !where( isnan(grid%t2)) grid%t2 = -999.
	 where( grid%t2 <= 200.) grid%t2 = 200.
!*isilda******
         open(unit=13,file="/stornext/home/rafael.stockler/projetos/BRAMS/data/q2.dat",&
	 	action="read",form="unformatted",recl=(ime+borders)*(jme+borders)*tbyte,access="direct")
         read(13, rec=1) ((grid%q2(i,j), i=ims,ime),j=jms,jme)
         close(unit=13)
	 
	 do j=1, jme
	 	do i=1, ime
	 		call SWAP_F4(grid%q2(i,j))
		enddo
	 enddo
	 
	 !where( isnan(grid%q2)) grid%q2 = -999.
	 where( grid%q2 <= 0.0000001) grid%q2 = 0.0000001
	 
	 
!*isilda******
         open(unit=13,file="/stornext/home/rafael.stockler/projetos/BRAMS/data/psfc.dat",&
	 	action="read",form="unformatted",recl=(ime+borders)*(jme+borders)*tbyte,access="direct")
         read(13, rec=1) ((grid%psfc(i,j), i=ims,ime),j=jms,jme)
         close(unit=13)
	 
	 do j=1, jme
	 	do i=1, ime
	 		call SWAP_F4(grid%psfc(i,j))
		enddo
	 enddo
	 
	 !where( isnan(grid%psfc)) grid%psfc = -999.
	 where( grid%psfc <= 1000.) grid%psfc = 1000.
	 
	 
!*isilda******
         open(unit=13,file="/stornext/home/rafael.stockler/projetos/BRAMS/data/rainc.dat",&
	 	action="read",form="unformatted",recl=(ime+borders)*(jme+borders)*tbyte,access="direct")
         read(13, rec=1) ((grid%rainc(i,j), i=ims,ime),j=jms,jme)
         close(unit=13)
	 
	do j=1, jme
	 	do i=1, ime
	 		call SWAP_F4(grid%rainc(i,j))
		enddo
	 enddo
	 
	 !where( isnan(grid%rainc)) grid%rainc = -999.
	 
	  grid%rainc = 1.
	 
	  
!*isilda******
         open(unit=13,file="/stornext/home/rafael.stockler/projetos/BRAMS/data/rainnc.dat",&
	 	action="read",form="unformatted",recl=(ime+borders)*(jme+borders)*tbyte,access="direct")
         read(13, rec=1) ((grid%rainnc(i,j), i=ims,ime),j=jms,jme)
         close(unit=13)
	 
	 do j=1, jme
	 	do i=1, ime
	 		call SWAP_F4(grid%rainnc(i,j))
		enddo
	 enddo
	 !where( isnan(grid%rainnc)) grid%rainnc = -999.
	 
	 grid%rainnc = 1.2
	 
	 
!*isilda******
!*isilda******
         open(unit=13,file="/stornext/home/rafael.stockler/projetos/BRAMS/data/zsf.dat",&
	 	action="read",form="unformatted",recl=(ifme+borders*grid%sr_x)*(jfme+borders*grid%sr_y)*tbyte,access="direct")
         read(13, rec=1) ((grid%zsf(i,j), i=ifms,ifme),j=jfms,jfme)
         close(unit=13)
	 
	 do j=1, jfme
	 	do i=1, ifme
	 		call SWAP_F4(grid%zsf(i,j))
		enddo
	enddo
	!where( isnan(grid%zsf)) grid%zsf = -999.
	

	 
!*isilda******
!*isilda******
         open(unit=13,file="/stornext/home/rafael.stockler/projetos/BRAMS/data/nfuel_cat.dat",&
	 	action="read",form="unformatted",recl=(ifme+borders*grid%sr_x)*(jfme+borders*grid%sr_y)*tbyte,access="direct")
         read(13, rec=1) ((grid%nfuel_cat(i,j), i=ifms,ifme),j=jfms,jfme)
         close(unit=13)
	 
	 do j=1, jfme
	 	do i=1, ifme
	 		call SWAP_F4(grid%nfuel_cat(i,j))
		enddo
	enddo
	!where( isnan(grid%nfuel_cat)) grid%nfuel_cat = -999.

   	
!*isilda******
!*isilda******
         open(unit=13,file="/stornext/home/rafael.stockler/projetos/BRAMS/data/ph_2.dat",&
	 	action="read",form="unformatted",recl=(ime+borders)*kme*(jme+borders)*tbyte,access="direct")
         read(13, rec=1) (((grid%ph_2(i,k,j), i=ims,ime),k=kms,kme),j=jms,jme)
         close(unit=13)
	 
	 do j=1, jme
	 	do k=1, kme
	 		do i=1, ime
	 			call SWAP_F4(grid%ph_2(i,k,j))
			enddo
		enddo
	 enddo
	 !where( isnan(grid%ph_2)) grid%ph_2 = -999.
	 
	
!*isilda******
!*isilda******
         open(unit=13,file="/stornext/home/rafael.stockler/projetos/BRAMS/data/z0.dat",&
	 	action="read",form="unformatted",recl=(ime+borders)*(jme+borders)*tbyte,access="direct")
         read(13, rec=1) ((grid%z0(i,j), i=ims,ime),j=jms,jme)
         close(unit=13)
	 
	 do j=1, jme
	 	do i=1, ime
	 		call SWAP_F4(grid%z0(i,j))
		enddo
	 enddo
	 !where( isnan(grid%z0)) grid%z0 = -999.
	 
	
!*isilda******
!*isilda******
        
open(unit=13,file="/stornext/home/rafael.stockler/projetos/BRAMS/data/ht.dat",&
	action="read",form="unformatted",recl=(ime+borders)*(jme+borders)*tbyte,access="direct")
         read(13, rec=1) ((grid%ht(i,j), i=ims,ime),j=jms,jme)
         close(unit=13)
	 
	 do j=1, jme
	 	do i=1, ime
	 		call SWAP_F4(grid%ht(i,j))
		enddo
	 enddo
	 !where( isnan(grid%ht)) grid%ht = -999.
	 
	 

	
print*, 'INIT SFIRE 1'


   grid%itimestep=0
   

   call sfire_driver_em_init ( grid , config_flags    &
            ,ids,ide, kds,kde, jds,jde                &
            ,ims,ime, kms,kme, jms,jme                &
            ,ips,ipe, kps,kpe, jps,jpe ) 

   !!CALL wrf_debug ( 100 , 'start_domain_em: After call to sfire_driver_em_init' )


nrec=0




!nrec=nrec+1
!   open(unit=13,file="/stornext/home/rafael.stockler/projetos/BRAMS/data/u_2_braco2.dat",action="read",form="unformatted",recl=ime*kme*jme*tbyte,access="direct", &
!	 convert="big_endian")
!         read(13, rec=nrec) grid%u_2
!         close(unit=13)
!print*, grid%u_2
!nrec=nrec+1
!   open(unit=13,file="/stornext/home/rafael.stockler/projetos/BRAMS/data/u_2_braco2.dat",action="read",form="unformatted",recl=ime*kme*jme*tbyte,access="direct", &
!	 convert="big_endian")
!         read(13, rec=nrec) grid%v_2
!         close(unit=13)
!
!print*, '-------'
!print*, grid%v_2

!do j=1, jme
!	 	do k=1, kme
!	 		do i=1, ime
!	 			if (grid%u_2(i,k,j) /= grid%v_2(i,k,j)) then
!					print*, 'DIFF :: ', grid%u_2(i,k,j), grid%v_2(i,k,j)
!				endif
!			enddo
!		enddo
!	 enddo


!stop

nrec=nrec+1
do t=1, 5  
   
   grid%itimestep = grid%itimestep + 1.
   write(tc, '(I4)') t
   
   print*, tc
   
   open(unit=13,file="/stornext/home/rafael.stockler/projetos/BRAMS/data/u_2_braco2_"// trim(adjustl(tc)) //".dat",&
	 	action="read",form="unformatted",recl=(ime+borders)*kme*(jme+borders)*tbyte,access="direct")
         read(13, rec=nrec) (((grid%u_2(i,k,j), i=ims,ime),k=kms,kme),j=jms,jme)
         close(unit=13)
	 
	 do j=1, jme
	 	do k=1, kme
	 		do i=1, ime
	 			call SWAP_F4(grid%u_2(i,k,j))
			enddo
		enddo
	 enddo
	 
	 !where( isnan(grid%u_2)) grid%u_2 = -999.
	 
	  
	 
	open(unit=13,file="/stornext/home/rafael.stockler/projetos/BRAMS/data/mut_braco2_"//trim(adjustl(tc)) //".dat",&
		action="read",form="unformatted",recl=(ime+borders)*(jme+borders)*tbyte,access="direct")

  	read(13, rec=nrec) ((grid%mut(i,j), i=ims,ime),j=jms,jme)

  	close(unit=13)
	
	do j=1, jme	
		do i=1, ime
			call SWAP_F4(grid%mut(i,j))

                enddo

        enddo

        !where( isnan(grid%mut)) grid%mut = -999.
	
	
		 
!*isilda*****
         open(unit=13,file="/stornext/home/rafael.stockler/projetos/BRAMS/data/v_2_braco2_"// trim(adjustl(tc)) //".dat",&
	 	action="read",form="unformatted",recl=(ime+borders)*kme*(jme+borders)*tbyte,access="direct")
         read(13, rec=nrec) (((grid%v_2(i,k,j), i=ims,ime),k=kms,kme),j=jms,jme)
         close(unit=13)
	 
	 do j=1, jme
	 	do k=1, kme
	 		do i=1, ime
	 			call SWAP_F4(grid%v_2(i,k,j))
			enddo
		enddo
	 enddo
	 
	 !where( isnan(grid%v_2)) grid%v_2 = -999.
	 
	
!*isilda******
!*isilda*****
         open(unit=13,file="/stornext/home/rafael.stockler/projetos/BRAMS/data/phb_braco2_"// trim(adjustl(tc)) //".dat",&
	 	action="read",form="unformatted",recl=(ime+borders)*kme*(jme+borders)*tbyte,access="direct")
         read(13, rec=nrec) (((grid%phb(i,k,j), i=ims,ime),k=kms,kme),j=jms,jme)
         close(unit=13)
	 
	 do j=1, jme
	 	do k=1, kme
	 		do i=1, ime
	 			call SWAP_F4(grid%phb(i,k,j))
			enddo
		enddo
	 enddo
	 
	 !where( isnan(grid%phb)) grid%phb = -999.
	 
	
!*isilda******
!*isilda*****
!         open(13,file='z_at_w.dat',action="write",form='unformatted')
 !        write(13) z_at_w
  !       close(13)
!*isilda******

!*isilda*****
         open(unit=13,file="/stornext/home/rafael.stockler/projetos/BRAMS/data/t2_braco2_"// trim(adjustl(tc)) //".dat",&
	 	action="read",form="unformatted",recl=(ime+borders)*(jme+borders)*tbyte,access="direct")
         read(13, rec=nrec) ((grid%t2(i,j), i=ims,ime),j=jms,jme)
         close(unit=13)
	 
	 do j=1, jme
	 	do i=1, ime
	 		call SWAP_F4(grid%t2(i,j))
		enddo
	 enddo
	 
	 !where( isnan(grid%t2)) grid%t2 = -999.
	 where( grid%t2 <= 200.) grid%t2 = 200.
	 
	 
!*isilda******
         open(unit=13,file="/stornext/home/rafael.stockler/projetos/BRAMS/data/q2_braco2_"// trim(adjustl(tc)) //".dat",&
	 	action="read",form="unformatted",recl=(ime+borders)*(jme+borders)*tbyte,access="direct")
         read(13, rec=nrec) ((grid%q2(i,j), i=ims,ime),j=jms,jme)
         close(unit=13)
	 
	 do j=1, jme
	 	do i=1, ime
	 		call SWAP_F4(grid%q2(i,j))
		enddo
	 enddo
	 
	 !where( isnan(grid%q2)) grid%q2 = -999.
	 where( grid%q2 <= 0.0000001) grid%q2 = 0.0000001
	 
	 
!*isilda******
         open(unit=13,file="/stornext/home/rafael.stockler/projetos/BRAMS/data/psfc_braco2_"// trim(adjustl(tc)) //".dat",&
	 	action="read",form="unformatted",recl=(ime+borders)*(jme+borders)*tbyte,access="direct")
         read(13, rec=nrec) ((grid%psfc(i,j), i=ims,ime),j=jms,jme)
         close(unit=13)
	 
	 do j=1, jme
	 	do i=1, ime
	 		call SWAP_F4(grid%psfc(i,j))
		enddo
	 enddo
	 
	 !where( isnan(grid%psfc)) grid%psfc = -999.
	 where( grid%psfc <= 1000.) grid%psfc = 1000.
	 
	
!*isilda******
         open(unit=13,file="/stornext/home/rafael.stockler/projetos/BRAMS/data/rainc_braco2_"// trim(adjustl(tc)) //".dat",&
	 	action="read",form="unformatted",recl=(ime+borders)*(jme+borders)*tbyte,access="direct")
         read(13, rec=nrec) ((grid%rainc(i,j), i=ims,ime),j=jms,jme)
         close(unit=13)
	 
	do j=1, jme
	 	do i=1, ime
	 		call SWAP_F4(grid%rainc(i,j))
		enddo
	 enddo
	 
	 !where( isnan(grid%rainc)) grid%rainc = -999.
	 
	  
	 
	 print*, 'rainc :: '
	 
	 do j=1, jme
	 	
	 		do i=1, ime
	 			if (grid%rainc(i,j) > 0.) then
					print*, grid%rainc(i,j)
				endif
			enddo
		
	 enddo
	 
	 print*, 'rainc :: '
	 grid%rainc = 1.
!*isilda******
         open(unit=13,file="/stornext/home/rafael.stockler/projetos/BRAMS/data/rainnc_braco2_"// trim(adjustl(tc)) //".dat",&
	 	action="read",form="unformatted",recl=(ime+borders)*(jme+borders)*tbyte,access="direct")
         read(13, rec=nrec) ((grid%rainnc(i,j), i=ims,ime),j=jms,jme)
         close(unit=13)
	 
	 do j=1, jme
	 	do i=1, ime
	 		call SWAP_F4(grid%rainnc(i,j))
		enddo
	 enddo
	 !where( isnan(grid%rainnc)) grid%rainnc = -999.
	 
	 
	 print*, 'rainnc :: '
	 
	 do j=1, jme
	 	
	 		do i=1, ime
	 			if (grid%rainnc(i,j) > 0.) then
					print*, grid%rainnc(i,j)
				endif
			enddo
		
	 enddo
	 
	 print*, 'rainnc :: '
	 grid%rainnc = 1.2
!*isilda******
!*isilda******
         open(unit=13,file="/stornext/home/rafael.stockler/projetos/BRAMS/data/zsf_braco2_"// trim(adjustl(tc)) //".dat",&
	 	action="read",form="unformatted",recl=(ifme+borders*grid%sr_x)*(jfme+borders*grid%sr_y)*tbyte,access="direct")
         read(13, rec=nrec) ((grid%zsf(i,j), i=ifms,ifme),j=jfms,jfme)
         close(unit=13)
	 
	 do j=1, jfme
	 	do i=1, ifme
	 		call SWAP_F4(grid%zsf(i,j))
		enddo
	enddo
	!where( isnan(grid%zsf)) grid%zsf = -999.
	
	
	 !open(unit=13,file="/stornext/home/rafael.stockler/projetos/BRAMS/data/fmc_g.dat",action="read",form="unformatted",recl=ifpe*jfpe*tbyte,access="direct")
         !read(13, rec=1) grid%fmc_g
        ! close(unit=13)
	 
!*isilda******
!*isilda******
         open(unit=13,file="/stornext/home/rafael.stockler/projetos/BRAMS/data/nfuel_cat_braco2_"// trim(adjustl(tc)) //".dat",&
	 	action="read",form="unformatted",recl=(ifme+borders*grid%sr_x)*(jfme+borders*grid%sr_y)*tbyte,access="direct")
         read(13, rec=nrec) ((grid%nfuel_cat(i,j), i=ifms,ifme),j=jfms,jfme)
         close(unit=13)
	 
	 do j=1, jfme
	 	do i=1, ifme
	 		call SWAP_F4(grid%nfuel_cat(i,j))
		enddo
	enddo
	!where( isnan(grid%nfuel_cat)) grid%nfuel_cat = -999.
	
	
!*isilda******
!*isilda******
 !        open(13,file='dz8w.dat',action="write",form='unformatted')
  !       write(13) grid%dz8w
   !      close(13)
!*isilda******
!*isilda******
         open(unit=13,file="/stornext/home/rafael.stockler/projetos/BRAMS/data/ph_2_braco2_"// trim(adjustl(tc)) //".dat",&
	 	action="read",form="unformatted",recl=(ime+borders)*kpe*(jme+borders)*tbyte,access="direct")
         read(13, rec=nrec) (((grid%ph_2(i,k,j), i=ims,ime),k=kms,kme),j=jms,jme)
         close(unit=13)
	 
	 do j=1, jme
	 	do k=1, kme
	 		do i=1, ime
	 			call SWAP_F4(grid%ph_2(i,k,j))
			enddo
		enddo
	 enddo
	 !where( isnan(grid%ph_2)) grid%ph_2 = -999.
	 
	
!*isilda******
!*isilda******
         open(unit=13,file="/stornext/home/rafael.stockler/projetos/BRAMS/data/z0_braco2_"// trim(adjustl(tc)) //".dat",&
	 	action="read",form="unformatted",recl=(ime+borders)*(jme+borders)*tbyte,access="direct")
         read(13, rec=nrec) ((grid%z0(i,j), i=ims,ime),j=jms,jme)
         close(unit=13)
	 
	 do j=1, jme
	 	do i=1, ime
	 		call SWAP_F4(grid%z0(i,j))
		enddo
	 enddo
	 !where( isnan(grid%z0)) grid%z0 = -999.
	 
	 
!*isilda******
!*isilda******
        
open(unit=13,file="/stornext/home/rafael.stockler/projetos/BRAMS/data/ht_braco2_"// trim(adjustl(tc)) //".dat",&
	action="read",form="unformatted",recl=(ime+borders)*(jme+borders)*tbyte,access="direct")
         read(13, rec=nrec) ((grid%ht(i,j), i=ims,ime),j=jms,jme)
         close(unit=13)
	 
	 do j=1, jme
	 	do i=1, ime
	 		call SWAP_F4(grid%ht(i,j))
		enddo
	 enddo
	 !where( isnan(grid%ht)) grid%ht = -999.
	 
	 
!*isilda******
!*isilda******
!         open(13,file='rho.dat',action="write",form='unformatted')
!         write(13) grid%rho
!         close(13)
!*isilda*****
!*isilda******
         open(unit=13,file="/stornext/home/rafael.stockler/projetos/BRAMS/data/rho_braco2_"// trim(adjustl(tc)) //".dat",&
	 	action="read",form="unformatted",recl=(ime+borders)*kme*(jme+borders)*tbyte,access="direct")
         read(13, rec=nrec) (((rho(i,k,j), i=ims,ime),k=kms,kme),j=jms,jme)
         close(unit=13)
	 
	 do j=1, jme
	 	do k=1, kme
	 		do i=1, ime
	 			call SWAP_F4(rho(i,k,j))
			enddo
		enddo
	 enddo
	 !where( isnan(rho)) rho = -999.
	 where( rho <= 0.00000000000000000000001) rho = 0.00001
	 
	
!*isilda******
!*isilda******
         open(unit=13,file="/stornext/home/rafael.stockler/projetos/BRAMS/data/dz8w_braco2_"// trim(adjustl(tc)) //".dat",&
	 	action="read",form="unformatted",recl=(ime+borders)*kme*(jme+borders)*tbyte,access="direct")
         read(13, rec=nrec) (((dz8w(i,k,j), i=ims,ime),k=kms,kme),j=jms,jme)
         close(unit=13)
	 
	 do j=1, jme
	 	do k=1, kme
	 		do i=1, ime
	 			call SWAP_F4(dz8w(i,k,j))
				
			enddo
		enddo
	enddo
	 
	 !where( isnan(dz8w)) dz8w = -999.
	 where( dz8w < 1.) dz8w = 1.
	 
	
	 
!*isilda******	
!*isilda*****
         open(unit=13,file="/stornext/home/rafael.stockler/projetos/BRAMS/data/z_at_w_braco2_"// trim(adjustl(tc)) //".dat",&
	 	action="read",form="unformatted",recl=(ime+borders)*kme*(jme+borders)*tbyte,access="direct")
         read(13, rec=nrec) (((grid%z_at_w(i,k,j), i=ims,ime),k=kms,kme),j=jms,jme)
         close(unit=13)
	 
	 do j=1, jme
	 	do k=1, kme
	 		do i=1, ime
	 			call SWAP_F4(grid%z_at_w(i,k,j))
			enddo
		enddo
	 enddo
	 !where( isnan(grid%z_at_w)) grid%z_at_w = -999.
	 
	 


print*, 'INIT SFIRE 2'


   	

!if(config_flags%ifire.eq.2)then 
!	! initialization moved to start_em:start_domain_em
!
!        ! one timestep of the fire model
        call sfire_driver_em_step (  grid , config_flags    &
            ,ids,ide, kds,kde, jds,jde                              &
            ,ims,ime, kms,kme, jms,jme                              &
            ,ips,ipe, kps,kpe, jps,jpe                              &
            ,rho,grid%z_at_w,dz8w)
	    
	
	
	
enddo
contains
  subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine check

end program fire_main
