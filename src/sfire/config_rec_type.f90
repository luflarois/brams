module config_rec_type

	use module_driver_constants
	implicit none

	!TYPE model_config_rec_type
      	!SEQUENCE
	!	include "namelist_defines.inc"

   	!END TYPE model_config_rec_type

   	TYPE grid_config_rec_type
		include "namelist_defines2.inc"
		
   	END TYPE grid_config_rec_type

end module config_rec_type
