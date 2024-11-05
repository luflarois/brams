module_fr_sfire_model.o: $(SFIRE)/module_fr_sfire_model.f90 module_fr_sfire_core.o \
	module_fr_sfire_phys.o module_fr_sfire_util.o 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

module_fr_sfire_util.o: $(SFIRE)/module_fr_sfire_util.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

module_fr_sfire_phys.o: $(SFIRE)/module_fr_sfire_phys.f90 module_fr_sfire_util.o \
	module_model_constants.o read_namelist_fire.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

module_fr_sfire_core.o: $(SFIRE)/module_fr_sfire_core.f90 module_fr_sfire_phys.o \
	module_fr_sfire_util.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

module_fr_sfire_atm.o: $(SFIRE)/module_fr_sfire_atm.f90 module_fr_sfire_phys.o \
	module_fr_sfire_util.o module_model_constants.o module_state_description.o \
	module_configure.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

module_fr_sfire_driver.o: $(SFIRE)/module_fr_sfire_driver.f90 	module_configure.o \
	module_domain_type.o module_fr_sfire_atm.o module_fr_sfire_model.o module_fr_sfire_phys.o \
	module_fr_sfire_util.o module_model_constants.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src	

byteswap.o: $(SFIRE)/byteswap.f90
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src	

module_state_description.o: $(SFIRE)/module_state_description.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

module_domain_type.o: $(SFIRE)/module_domain_type.f90 
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

read_namelist_fire.o: $(SFIRE)/read_namelist_fire.f90
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

module_configure.o: $(SFIRE)/module_configure.f90 module_domain_type.o \
	module_driver_constants.o module_state_description.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

module_machine.o: $(SFIRE)/module_machine.f90 module_driver_constants.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

module_tiles.o: $(SFIRE)/module_tiles.f90 module_configure.o module_domain_type.o \
	module_driver_constants.o module_machine.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

module_model_constants.o: $(SFIRE)/module_model_constants.f90
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

module_driver_constants.o: $(SFIRE)/module_driver_constants.f90
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

mem_sfire.o: $(SFIRE)/mem_sfire.f90 module_domain_type.o read_namelist_fire.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

read_namelist_fire.o: $(SFIRE)/read_namelist_fire.f90
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

module_fr_sfire_driver_brams.o: $(SFIRE)/module_fr_sfire_driver_brams.f90 module_domain_type.o \
	module_fr_sfire_atm.o module_fr_sfire_driver.o module_fr_sfire_util.o module_state_description.o \
	read_namelist_fire.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

modSfire.o: $(SFIRE)/modSfire.f90 read_namelist_fire.o mem_sfire.o module_fr_sfire_driver_brams.o \
	modSfire2Brams.o
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND_LIGHT) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src

modSfire2Brams.o : $(SFIRE)/modSfire2Brams.f90
	@cp -f $< $(<F:.f90=.f90)
	$(F_COMMAND) $(<F:.f90=.f90) $(EXTRAFLAGSF)
	@mv -f $(<F:.f90=.f90) ../doc/src
