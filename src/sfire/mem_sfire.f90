MODULE mem_sfire
  USE module_domain_type
  USE ModNamelistsfireFile 
  use ModNamelistFile, only: namelistFile
    
  implicit none
  TYPE(domain) :: sfire_g
  TYPE(grid_config_rec_type), pointer :: config_flags => null()

  integer :: sfire
  

CONTAINS

  SUBROUTINE alloc_sfire_brams(sfire &
  	,ims,ime, kms,kme, jms,jme &
	,ifms,ifme, jfms,jfme, nfmc)

    IMPLICIT NONE
    TYPE (domain) :: sfire
    INTEGER, INTENT(in) :: ims,ime, kms,kme, jms,jme &
	,ifms,ifme, jfms,jfme, nfmc

    ! Allocate arrays based on options (if necessary)
    
    ALLOCATE(sfire%xlong(ims:ime,jms:jme))
    ALLOCATE(sfire%xlat(ims:ime,jms:jme))   
    ALLOCATE(sfire%rho(ims:ime,kms:kme,jms:jme))
    ALLOCATE(sfire%dz8w(ims:ime,kms:kme,jms:jme)) 
    ALLOCATE(sfire%u_2(ims:ime,kms:kme,jms:jme))
    ALLOCATE(sfire%v_2(ims:ime,kms:kme,jms:jme))
    ALLOCATE(sfire%ph_2(ims:ime,kms:kme,jms:jme))
    ALLOCATE(sfire%phb(ims:ime,kms:kme,jms:jme))
    ALLOCATE(sfire%z0(ims:ime,jms:jme))
    ALLOCATE(sfire%ht(ims:ime,jms:jme))
    ALLOCATE(sfire%rain_old(ims:ime,jms:jme))
    ALLOCATE(sfire%t2(ims:ime,jms:jme))
    ALLOCATE(sfire%q2(ims:ime,jms:jme))
    ALLOCATE(sfire%psfc(ims:ime,jms:jme))
    ALLOCATE(sfire%rainc(ims:ime,jms:jme))
    ALLOCATE(sfire%rainnc(ims:ime,jms:jme))
    ALLOCATE(sfire%t2_old(ims:ime,jms:jme))
    ALLOCATE(sfire%q2_old(ims:ime,jms:jme))
    ALLOCATE(sfire%psfc_old(ims:ime,jms:jme))
    ALLOCATE(sfire%rh_fire(ims:ime,jms:jme))
    !LFR
    ALLOCATE(sfire%fndwi(ifms:ifme,jfms:jfme))
    ALLOCATE(sfire%ndwi(ifms:ifme,jfms:jfme))
    ALLOCATE(sfire%tign_in(ifms:ifme, jfms:jfme))
    ALLOCATE(sfire%can_top(ifms:ifme, jfms:jfme))
    ALLOCATE(sfire%cuf(ifms:ifme, jfms:jfme))
    ALLOCATE(sfire%cvf(ifms:ifme, jfms:jfme))
    !LFR
    ALLOCATE(sfire%fmc_gc(ims:ime,1:nfmc,jms:jme))
    ALLOCATE(sfire%fmc_equi(ims:ime,1:nfmc,jms:jme))
    ALLOCATE(sfire%fmc_lag(ims:ime,1:nfmc,jms:jme))
    ALLOCATE(sfire%lfn(ifms:ifme,jfms:jfme))
    ALLOCATE(sfire%lfn_out(ifms:ifme,jfms:jfme)) !INTRODUZIDO POR ISILDA CM
    ALLOCATE(sfire%tign_g(ifms:ifme,jfms:jfme))
    ALLOCATE(sfire%fire_area(ifms:ifme,jfms:jfme))
    ALLOCATE(sfire%fuel_frac(ifms:ifme,jfms:jfme))
    ALLOCATE(sfire%zsf(ifms:ifme,jfms:jfme))
    ALLOCATE(sfire%dzdxf(ifms:ifme,jfms:jfme))
    ALLOCATE(sfire%dzdyf(ifms:ifme,jfms:jfme))
    ALLOCATE(sfire%fmc_g(ifms:ifme,jfms:jfme))
    ALLOCATE(sfire%nfuel_cat(ifms:ifme,jfms:jfme))
    ALLOCATE(sfire%fuel_time(ifms:ifme,jfms:jfme))
    ALLOCATE(sfire%fuel_frac_burnt(ifms:ifme,jfms:jfme))
    ALLOCATE(sfire%uf(ifms:ifme,jfms:jfme))
    ALLOCATE(sfire%vf(ifms:ifme,jfms:jfme))
    ALLOCATE(sfire%bbb(ifms:ifme,jfms:jfme))
    ALLOCATE(sfire%phisc(ifms:ifme,jfms:jfme))
    ALLOCATE(sfire%phiwc(ifms:ifme,jfms:jfme))
    ALLOCATE(sfire%r_0(ifms:ifme,jfms:jfme))
    ALLOCATE(sfire%bbb_FM10(ifms:ifme,jfms:jfme)) !INTRODUZIDO POR ISILDA CM
    ALLOCATE(sfire%phisc_FM10(ifms:ifme,jfms:jfme))!INTRODUZIDO POR ISILDA CM
    ALLOCATE(sfire%phiwc_FM10(ifms:ifme,jfms:jfme))!INTRODUZIDO POR ISILDA CM
    ALLOCATE(sfire%r_0_FM10(ifms:ifme,jfms:jfme))!INTRODUZIDO POR ISILDA CM
    ALLOCATE(sfire%CFB(ifms:ifme,jfms:jfme))!INTRODUZIDO POR ISILDA CM
    ALLOCATE(sfire%HPA(ifms:ifme,jfms:jfme))!INTRODUZIDO POR ISILDA CM
    ALLOCATE(sfire%fgip(ifms:ifme,jfms:jfme))
    ALLOCATE(sfire%ischap(ifms:ifme,jfms:jfme))
    ALLOCATE(sfire%fz0(ifms:ifme,jfms:jfme))
    ALLOCATE(sfire%fwh(ifms:ifme,jfms:jfme))
    ALLOCATE(sfire%wz0(ifms:ifme,jfms:jfme)) !INTRODUZIDO POR ISILDA CM
    ALLOCATE(sfire%avg_fuel_frac(ims:ime, jms:jme))
    ALLOCATE(sfire%grnhfx(ims:ime, jms:jme))
    ALLOCATE(sfire%grnqfx(ims:ime, jms:jme))
    ALLOCATE(sfire%canhfx(ims:ime, jms:jme))
    ALLOCATE(sfire%canqfx(ims:ime, jms:jme))
    ALLOCATE(sfire%uah(ims:ime,jms:jme))
    ALLOCATE(sfire%vah(ims:ime,jms:jme))
    ALLOCATE(sfire%fgrnhfx(ifms:ifme, jfms:jfme))
    ALLOCATE(sfire%fgrnqfx(ifms:ifme, jfms:jfme))
    ALLOCATE(sfire%fcanhfx(ifms:ifme, jfms:jfme))
    ALLOCATE(sfire%fcanqfx(ifms:ifme, jfms:jfme))
    ALLOCATE(sfire%ros(ifms:ifme,jfms:jfme))
    ALLOCATE(sfire%flineint(ifms:ifme,jfms:jfme))
    ALLOCATE(sfire%flineint2(ifms:ifme,jfms:jfme))
    ALLOCATE(sfire%flineint_total(ifms:ifme,jfms:jfme))  !!!!INTRODUZIDO POR ISILDA CUNHA MENEZES
    ALLOCATE(sfire%f_ros0(ifms:ifme,jfms:jfme))
    ALLOCATE(sfire%f_rosx(ifms:ifme,jfms:jfme))
    ALLOCATE(sfire%f_rosy(ifms:ifme,jfms:jfme))
    ALLOCATE(sfire%f_ros(ifms:ifme,jfms:jfme))
    ALLOCATE(sfire%f_ros11(ifms:ifme,jfms:jfme))!!!!INTRODUZIDO POR ISILDA CUNHA MENEZES
    ALLOCATE(sfire%f_ros12(ifms:ifme,jfms:jfme))!!!!INTRODUZIDO POR ISILDA CUNHA MENEZES
    ALLOCATE(sfire%f_ros13(ifms:ifme,jfms:jfme))!!!!INTRODUZIDO POR ISILDA CUNHA MENEZES
    ALLOCATE(sfire%f_ros21(ifms:ifme,jfms:jfme))!!!!INTRODUZIDO POR ISILDA CUNHA MENEZES
    ALLOCATE(sfire%f_ros23(ifms:ifme,jfms:jfme))!!!!INTRODUZIDO POR ISILDA CUNHA MENEZES
    ALLOCATE(sfire%f_ros31(ifms:ifme,jfms:jfme))!!!!INTRODUZIDO POR ISILDA CUNHA MENEZES
    ALLOCATE(sfire%f_ros32(ifms:ifme,jfms:jfme))!!!!INTRODUZIDO POR ISILDA CUNHA MENEZES
    ALLOCATE(sfire%f_ros33(ifms:ifme,jfms:jfme))!!!!INTRODUZIDO POR ISILDA CUNHA MENEZES
    ALLOCATE(sfire%f_int(ifms:ifme,jfms:jfme))
    ALLOCATE(sfire%f_lineint(ifms:ifme,jfms:jfme))
    ALLOCATE(sfire%f_lineint2(ifms:ifme,jfms:jfme))
    ALLOCATE(sfire%fxlong(ifms:ifme, jfms:jfme))
    ALLOCATE(sfire%fxlat(ifms:ifme, jfms:jfme))
    ALLOCATE(sfire%fire_hfx(ifms:ifme, jfms:jfme))
    ALLOCATE(sfire%z_at_w(ims:ime, kms:kme, jms:jme))	
    ALLOCATE(sfire%mut(ims:ime, jms:jme))
    ALLOCATE(sfire%rthfrten(ims:ime, kms:kme, jms:jme))
    ALLOCATE(sfire%rqvfrten(ims:ime, kms:kme, jms:jme))
    ALLOCATE(sfire%fmep(ims:ime,1:2,jms:jme)) !!!!INTRODUZIDO POR ISILDA CUNHA MENEZES
    RETURN
  END SUBROUTINE alloc_sfire_brams

  SUBROUTINE nullify_sfire_brams(sfire)

    IMPLICIT NONE
    
    TYPE (domain) :: sfire
    


    IF (ASSOCIATED(sfire%xlong))  NULLIFY (sfire%xlong)
    IF (ASSOCIATED(sfire%xlat))  NULLIFY (sfire%xlat)   
    IF (ASSOCIATED(sfire%u_2))  NULLIFY (sfire%u_2)
    IF (ASSOCIATED(sfire%v_2))  NULLIFY (sfire%v_2)
    IF (ASSOCIATED(sfire%ph_2))  NULLIFY (sfire%ph_2)
    IF (ASSOCIATED(sfire%phb))  NULLIFY (sfire%phb)
    IF (ASSOCIATED(sfire%z0))  NULLIFY (sfire%z0)
    IF (ASSOCIATED(sfire%ht))  NULLIFY (sfire%ht)
    IF (ASSOCIATED(sfire%rain_old))  NULLIFY (sfire%rain_old)
    IF (ASSOCIATED(sfire%t2))  NULLIFY (sfire%t2)
    IF (ASSOCIATED(sfire%q2))  NULLIFY (sfire%q2)
    IF (ASSOCIATED(sfire%psfc))  NULLIFY (sfire%psfc)
    IF (ASSOCIATED(sfire%rainc))  NULLIFY (sfire%rainc)
    IF (ASSOCIATED(sfire%rainnc))  NULLIFY (sfire%rainnc)
    IF (ASSOCIATED(sfire%t2_old))  NULLIFY (sfire%t2_old)
    IF (ASSOCIATED(sfire%q2_old))  NULLIFY (sfire%q2_old)
    IF (ASSOCIATED(sfire%psfc_old))  NULLIFY (sfire%psfc_old)
    IF (ASSOCIATED(sfire%rh_fire))  NULLIFY (sfire%rh_fire)
    IF (ASSOCIATED(sfire%fmc_gc))  NULLIFY (sfire%fmc_gc)
    IF (ASSOCIATED(sfire%fmc_equi))  NULLIFY (sfire%fmc_equi)
    IF (ASSOCIATED(sfire%fmc_lag))  NULLIFY (sfire%fmc_lag)
    IF (ASSOCIATED(sfire%lfn))  NULLIFY (sfire%lfn)
    IF (ASSOCIATED(sfire%lfn_out))  NULLIFY (sfire%lfn_out)!INTRODUZIDO POR ISILDA CM
    IF (ASSOCIATED(sfire%tign_g))  NULLIFY (sfire%tign_g)
    IF (ASSOCIATED(sfire%fire_area))  NULLIFY (sfire%fire_area)
    IF (ASSOCIATED(sfire%fuel_frac))  NULLIFY (sfire%fuel_frac)
    IF (ASSOCIATED(sfire%zsf))  NULLIFY (sfire%zsf)
    IF (ASSOCIATED(sfire%dzdxf))  NULLIFY (sfire%dzdxf)
    IF (ASSOCIATED(sfire%dzdyf))  NULLIFY (sfire%dzdyf)
    IF (ASSOCIATED(sfire%fmc_g))  NULLIFY (sfire%fmc_g)
    IF (ASSOCIATED(sfire%nfuel_cat))  NULLIFY (sfire%nfuel_cat)
    IF (ASSOCIATED(sfire%fuel_time))  NULLIFY (sfire%fuel_time)
    IF (ASSOCIATED(sfire%fuel_frac_burnt))  NULLIFY (sfire%fuel_frac_burnt)
    IF (ASSOCIATED(sfire%uf))  NULLIFY (sfire%uf)
    IF (ASSOCIATED(sfire%vf))  NULLIFY (sfire%vf)
    IF (ASSOCIATED(sfire%bbb))  NULLIFY (sfire%bbb)
    IF (ASSOCIATED(sfire%phisc))  NULLIFY (sfire%phisc)
    IF (ASSOCIATED(sfire%phiwc))  NULLIFY (sfire%phiwc)
    IF (ASSOCIATED(sfire%r_0))  NULLIFY (sfire%r_0)
    IF (ASSOCIATED(sfire%bbb_FM10))  NULLIFY (sfire%bbb_FM10)!INTRODUZIDO POR ISILDA CM
    IF (ASSOCIATED(sfire%phisc_FM10))  NULLIFY (sfire%phisc_FM10)!INTRODUZIDO POR ISILDA CM
    IF (ASSOCIATED(sfire%phiwc_FM10))  NULLIFY (sfire%phiwc_FM10)!INTRODUZIDO POR ISILDA CM
    IF (ASSOCIATED(sfire%r_0_FM10))  NULLIFY (sfire%r_0_FM10)!INTRODUZIDO POR ISILDA CM
    IF (ASSOCIATED(sfire%CFB))  NULLIFY (sfire%CFB)!INTRODUZIDO POR ISILDA CM
    IF (ASSOCIATED(sfire%HPA))  NULLIFY (sfire%HPA)!INTRODUZIDO POR ISILDA CM
    IF (ASSOCIATED(sfire%fgip))  NULLIFY (sfire%fgip)
    IF (ASSOCIATED(sfire%ischap))  NULLIFY (sfire%ischap)
    IF (ASSOCIATED(sfire%fz0))  NULLIFY (sfire%fz0)
    IF (ASSOCIATED(sfire%fwh))  NULLIFY (sfire%fwh)
    IF (ASSOCIATED(sfire%wz0))  NULLIFY (sfire%wz0) !INTRODUZIDO POR ISILDA CM
    IF (ASSOCIATED(sfire%avg_fuel_frac))  NULLIFY (sfire%avg_fuel_frac)
    IF (ASSOCIATED(sfire%grnhfx))  NULLIFY (sfire%grnhfx)
    IF (ASSOCIATED(sfire%grnqfx))  NULLIFY (sfire%grnqfx)
    IF (ASSOCIATED(sfire%canhfx))  NULLIFY (sfire%canhfx)
    IF (ASSOCIATED(sfire%canqfx))  NULLIFY (sfire%canqfx)
    IF (ASSOCIATED(sfire%uah))  NULLIFY (sfire%uah)
    IF (ASSOCIATED(sfire%vah))  NULLIFY (sfire%vah)
    IF (ASSOCIATED(sfire%fgrnhfx))  NULLIFY (sfire%fgrnhfx)
    IF (ASSOCIATED(sfire%fgrnqfx))  NULLIFY (sfire%fgrnqfx)
    IF (ASSOCIATED(sfire%fcanhfx))  NULLIFY (sfire%fcanhfx)
    IF (ASSOCIATED(sfire%fcanqfx))  NULLIFY (sfire%fcanqfx)
    IF (ASSOCIATED(sfire%ros))  NULLIFY (sfire%ros)
    IF (ASSOCIATED(sfire%flineint))  NULLIFY (sfire%flineint)
    IF (ASSOCIATED(sfire%flineint2))  NULLIFY (sfire%flineint2)
	IF (ASSOCIATED(sfire%flineint_total))  NULLIFY (sfire%flineint_total)  !!!!INTRODUZIDO POR ISILDA CUNHA MENEZES
    IF (ASSOCIATED(sfire%f_ros0))  NULLIFY (sfire%f_ros0)
    IF (ASSOCIATED(sfire%f_rosx))  NULLIFY (sfire%f_rosx)
    IF (ASSOCIATED(sfire%f_rosy))  NULLIFY (sfire%f_rosy)
    IF (ASSOCIATED(sfire%f_ros))  NULLIFY (sfire%f_ros)
    IF (ASSOCIATED(sfire%f_ros11))  NULLIFY (sfire%f_ros11)!!!!INTRODUZIDO POR ISILDA CUNHA MENEZES
    IF (ASSOCIATED(sfire%f_ros12))  NULLIFY (sfire%f_ros12)!!!!INTRODUZIDO POR ISILDA CUNHA MENEZES
    IF (ASSOCIATED(sfire%f_ros13))  NULLIFY (sfire%f_ros13)!!!!INTRODUZIDO POR ISILDA CUNHA MENEZES
    IF (ASSOCIATED(sfire%f_ros21))  NULLIFY (sfire%f_ros21)!!!!INTRODUZIDO POR ISILDA CUNHA MENEZES
    IF (ASSOCIATED(sfire%f_ros23))  NULLIFY (sfire%f_ros23)!!!!INTRODUZIDO POR ISILDA CUNHA MENEZES
    IF (ASSOCIATED(sfire%f_ros31))  NULLIFY (sfire%f_ros31)!!!!INTRODUZIDO POR ISILDA CUNHA MENEZES
    IF (ASSOCIATED(sfire%f_ros32))  NULLIFY (sfire%f_ros32)!!!!INTRODUZIDO POR ISILDA CUNHA MENEZES
    IF (ASSOCIATED(sfire%f_ros33))  NULLIFY (sfire%f_ros33)!!!!INTRODUZIDO POR ISILDA CUNHA MENEZES
    IF (ASSOCIATED(sfire%f_int))  NULLIFY (sfire%f_int)
    IF (ASSOCIATED(sfire%f_lineint))  NULLIFY (sfire%f_lineint)
    IF (ASSOCIATED(sfire%f_lineint2))  NULLIFY (sfire%f_lineint2)
    IF (ASSOCIATED(sfire%fxlong))  NULLIFY (sfire%fxlong)
    IF (ASSOCIATED(sfire%fxlat))  NULLIFY (sfire%fxlat)
    IF (ASSOCIATED(sfire%fire_hfx))  NULLIFY (sfire%fire_hfx)
    IF (ASSOCIATED(sfire%z_at_w))  NULLIFY (sfire%z_at_w)
    IF (ASSOCIATED(sfire%mut))  NULLIFY (sfire%mut)
    IF (ASSOCIATED(sfire%rthfrten))  NULLIFY (sfire%rthfrten)
    IF (ASSOCIATED(sfire%rqvfrten))  NULLIFY (sfire%rqvfrten)
    IF (ASSOCIATED(sfire%fmep))  NULLIFY (sfire%fmep) !!!!INTRODUZIDO POR ISILDA CUNHA MENEZES
    RETURN
  END SUBROUTINE nullify_sfire_brams

  SUBROUTINE dealloc_sfire_brams(sfire)

    IMPLICIT NONE
    
    TYPE (domain) :: sfire
    
    IF (ASSOCIATED(sfire%xlong))  DEALLOCATE (sfire%xlong)
    IF (ASSOCIATED(sfire%xlat))  DEALLOCATE (sfire%xlat)   
    IF (ASSOCIATED(sfire%u_2))  DEALLOCATE (sfire%u_2)
    IF (ASSOCIATED(sfire%v_2))  DEALLOCATE (sfire%v_2)
    IF (ASSOCIATED(sfire%ph_2))  DEALLOCATE (sfire%ph_2)
    IF (ASSOCIATED(sfire%phb))  DEALLOCATE (sfire%phb)
    IF (ASSOCIATED(sfire%z0))  DEALLOCATE (sfire%z0)
    IF (ASSOCIATED(sfire%ht))  DEALLOCATE (sfire%ht)
    IF (ASSOCIATED(sfire%rain_old))  DEALLOCATE (sfire%rain_old)
    IF (ASSOCIATED(sfire%t2))  DEALLOCATE (sfire%t2)
    IF (ASSOCIATED(sfire%q2))  DEALLOCATE (sfire%q2)
    IF (ASSOCIATED(sfire%psfc))  DEALLOCATE (sfire%psfc)
    IF (ASSOCIATED(sfire%rainc))  DEALLOCATE (sfire%rainc)
    IF (ASSOCIATED(sfire%rainnc))  DEALLOCATE (sfire%rainnc)
    IF (ASSOCIATED(sfire%t2_old))  DEALLOCATE (sfire%t2_old)
    IF (ASSOCIATED(sfire%q2_old))  DEALLOCATE (sfire%q2_old)
    IF (ASSOCIATED(sfire%psfc_old))  DEALLOCATE (sfire%psfc_old)
    IF (ASSOCIATED(sfire%rh_fire))  DEALLOCATE (sfire%rh_fire)
    IF (ASSOCIATED(sfire%fmc_gc))  DEALLOCATE (sfire%fmc_gc)
    IF (ASSOCIATED(sfire%fmc_equi))  DEALLOCATE (sfire%fmc_equi)
    IF (ASSOCIATED(sfire%fmc_lag))  DEALLOCATE (sfire%fmc_lag)
    IF (ASSOCIATED(sfire%lfn))  DEALLOCATE (sfire%lfn)
    IF (ASSOCIATED(sfire%lfn_out))  DEALLOCATE (sfire%lfn_out) !INTRODUZIDO POR ISILDA CUNHA MENEZES
    IF (ASSOCIATED(sfire%tign_g))  DEALLOCATE (sfire%tign_g)
    IF (ASSOCIATED(sfire%fire_area))  DEALLOCATE (sfire%fire_area)
    IF (ASSOCIATED(sfire%fuel_frac))  DEALLOCATE (sfire%fuel_frac)
    IF (ASSOCIATED(sfire%zsf))  DEALLOCATE (sfire%zsf)
    IF (ASSOCIATED(sfire%dzdxf))  DEALLOCATE (sfire%dzdxf)
    IF (ASSOCIATED(sfire%dzdyf))  DEALLOCATE (sfire%dzdyf)
    IF (ASSOCIATED(sfire%fmc_g))  DEALLOCATE (sfire%fmc_g)
    IF (ASSOCIATED(sfire%nfuel_cat))  DEALLOCATE (sfire%nfuel_cat)
    IF (ASSOCIATED(sfire%fuel_time))  DEALLOCATE (sfire%fuel_time)
    IF (ASSOCIATED(sfire%fuel_frac_burnt))  DEALLOCATE (sfire%fuel_frac_burnt)
    IF (ASSOCIATED(sfire%uf))  DEALLOCATE (sfire%uf)
    IF (ASSOCIATED(sfire%vf))  DEALLOCATE (sfire%vf)
    IF (ASSOCIATED(sfire%bbb))  DEALLOCATE (sfire%bbb)
    IF (ASSOCIATED(sfire%phisc))  DEALLOCATE (sfire%phisc)
    IF (ASSOCIATED(sfire%phiwc))  DEALLOCATE (sfire%phiwc)
    IF (ASSOCIATED(sfire%r_0))  DEALLOCATE (sfire%r_0)
    IF (ASSOCIATED(sfire%bbb_FM10))  DEALLOCATE (sfire%bbb)!INTRODUZIDO POR ISILDA CM 
    IF (ASSOCIATED(sfire%phisc_FM10))  DEALLOCATE (sfire%phisc)!INTRODUZIDO POR ISILDA CM
    IF (ASSOCIATED(sfire%phiwc_FM10))  DEALLOCATE (sfire%phiwc)!INTRODUZIDO POR ISILDA CM
    IF (ASSOCIATED(sfire%r_0_FM10))  DEALLOCATE (sfire%r_0)!INTRODUZIDO POR ISILDA CM
    IF (ASSOCIATED(sfire%CFB))  DEALLOCATE (sfire%CFB)!INTRODUZIDO POR ISILDA CM
    IF (ASSOCIATED(sfire%HPA))  DEALLOCATE (sfire%HPA)!INTRODUZIDO POR ISILDA CM
    IF (ASSOCIATED(sfire%fgip))  DEALLOCATE (sfire%fgip)
    IF (ASSOCIATED(sfire%ischap))  DEALLOCATE (sfire%ischap)
    IF (ASSOCIATED(sfire%fz0))  DEALLOCATE (sfire%fz0)
    IF (ASSOCIATED(sfire%fwh))  DEALLOCATE (sfire%fwh)
    IF (ASSOCIATED(sfire%wz0))  DEALLOCATE (sfire%wz0) !INTRODUZIDO POR ISILDA CUNHA MENEZES
    IF (ASSOCIATED(sfire%avg_fuel_frac))  DEALLOCATE (sfire%avg_fuel_frac)
    IF (ASSOCIATED(sfire%grnhfx))  DEALLOCATE (sfire%grnhfx)
    IF (ASSOCIATED(sfire%grnqfx))  DEALLOCATE (sfire%grnqfx)
    IF (ASSOCIATED(sfire%canhfx))  DEALLOCATE (sfire%canhfx)
    IF (ASSOCIATED(sfire%canqfx))  DEALLOCATE (sfire%canqfx)
    IF (ASSOCIATED(sfire%uah))  DEALLOCATE (sfire%uah)
    IF (ASSOCIATED(sfire%vah))  DEALLOCATE (sfire%vah)
    IF (ASSOCIATED(sfire%fgrnhfx))  DEALLOCATE (sfire%fgrnhfx)
    IF (ASSOCIATED(sfire%fgrnqfx))  DEALLOCATE (sfire%fgrnqfx)
    IF (ASSOCIATED(sfire%fcanhfx))  DEALLOCATE (sfire%fcanhfx)
    IF (ASSOCIATED(sfire%fcanqfx))  DEALLOCATE (sfire%fcanqfx)
    IF (ASSOCIATED(sfire%ros))  DEALLOCATE (sfire%ros)
    IF (ASSOCIATED(sfire%flineint))  DEALLOCATE (sfire%flineint)
    IF (ASSOCIATED(sfire%flineint2)) DEALLOCATE  (sfire%flineint2)
	IF (ASSOCIATED(sfire%flineint_total)) DEALLOCATE (sfire%flineint_total)!!!!INTRODUZIDO POR ISILDA CUNHA MENEZES
    IF (ASSOCIATED(sfire%f_ros0))  DEALLOCATE (sfire%f_ros0)
    IF (ASSOCIATED(sfire%f_rosx))  DEALLOCATE (sfire%f_rosx)
    IF (ASSOCIATED(sfire%f_rosy))  DEALLOCATE (sfire%f_rosy)
    IF (ASSOCIATED(sfire%f_ros))  DEALLOCATE (sfire%f_ros)
    IF (ASSOCIATED(sfire%f_ros11))  DEALLOCATE (sfire%f_ros11)!!!!INTRODUZIDO POR ISILDA CUNHA MENEZES
    IF (ASSOCIATED(sfire%f_ros12))  DEALLOCATE (sfire%f_ros12)!!!!INTRODUZIDO POR ISILDA CUNHA MENEZES
    IF (ASSOCIATED(sfire%f_ros13))  DEALLOCATE (sfire%f_ros13)!!!!INTRODUZIDO POR ISILDA CUNHA MENEZES
    IF (ASSOCIATED(sfire%f_ros21))  DEALLOCATE (sfire%f_ros21)!!!!INTRODUZIDO POR ISILDA CUNHA MENEZES
    IF (ASSOCIATED(sfire%f_ros23))  DEALLOCATE (sfire%f_ros23)!!!!INTRODUZIDO POR ISILDA CUNHA MENEZES
    IF (ASSOCIATED(sfire%f_ros31))  DEALLOCATE (sfire%f_ros31)!!!!INTRODUZIDO POR ISILDA CUNHA MENEZES
    IF (ASSOCIATED(sfire%f_ros32))  DEALLOCATE (sfire%f_ros32)!!!!INTRODUZIDO POR ISILDA CUNHA MENEZES
    IF (ASSOCIATED(sfire%f_ros33))  DEALLOCATE (sfire%f_ros33)!!!!INTRODUZIDO POR ISILDA CUNHA MENEZES
    IF (ASSOCIATED(sfire%f_int))  DEALLOCATE (sfire%f_int)
    IF (ASSOCIATED(sfire%f_lineint))  DEALLOCATE (sfire%f_lineint)
    IF (ASSOCIATED(sfire%f_lineint2))  DEALLOCATE (sfire%f_lineint2)
    IF (ASSOCIATED(sfire%fxlong))  DEALLOCATE (sfire%fxlong)
    IF (ASSOCIATED(sfire%fxlat))  DEALLOCATE (sfire%fxlat)
    IF (ASSOCIATED(sfire%fire_hfx))  DEALLOCATE (sfire%fire_hfx)
    IF (ASSOCIATED(sfire%z_at_w))  DEALLOCATE (sfire%z_at_w)
    IF (ASSOCIATED(sfire%mut))  DEALLOCATE (sfire%mut)
    IF (ASSOCIATED(sfire%rthfrten))  DEALLOCATE (sfire%rthfrten)
    IF (ASSOCIATED(sfire%rqvfrten))  DEALLOCATE (sfire%rqvfrten)
    IF (ASSOCIATED(sfire%fmep))  DEALLOCATE (sfire%fmep) !!!!INTRODUZIDO POR ISILDA CUNHA
    RETURN
  END SUBROUTINE dealloc_sfire_brams

  

  SUBROUTINE zero_sfire_brams(sfire)

    IMPLICIT NONE
   
    TYPE (domain) :: sfire

    !LFR
    sfire%fndwi(:,:) = 0.0
    sfire%ndwi(:,:) = 0.0
    sfire%tign_in(:,:) = 0.0
    sfire%can_top(:,:) = 0.0
    sfire%cuf(:,:) = 0.0
    sfire%cvf(:,:) = 0.0
    !LFR
    sfire%lfn(:,:)= 0.
    sfire%lfn_out(:,:)=0.  !INTRODUZIDO POR ISILDA CUNHA MENEZES
    sfire%tign_g(:,:)= 0.
    sfire%fire_area(:,:)= 0.
    sfire%fuel_frac(:,:)= 0.
    sfire%zsf(:,:)= 0.
    sfire%dzdxf(:,:)= 0.
    sfire%dzdyf(:,:)= 0.
    sfire%fmc_g(:,:)= 0.
    sfire%nfuel_cat(:,:)= 0.
    sfire%fuel_time(:,:) = 0.
    sfire%fuel_frac_burnt(:,:)= 0.
    sfire%uf(:,:)= 0.
    sfire%vf(:,:)= 0.
    sfire%bbb(:,:)= 0.
    sfire%phisc(:,:)= 0.
    sfire%phiwc(:,:)= 0.
    sfire%r_0(:,:)= 0.
    sfire%bbb_FM10(:,:)= 0. !INTRODUZIDO POR ISILDA CM
    sfire%phisc_FM10(:,:)= 0.!INTRODUZIDO POR ISILDA CM
    sfire%phiwc_FM10(:,:)= 0.!INTRODUZIDO POR ISILDA CM
    sfire%r_0_FM10(:,:)= 0.!INTRODUZIDO POR ISILDA CM
    sfire%CFB(:,:)= 0.!INTRODUZIDO POR ISILDA CM
    sfire%HPA(:,:)= 0.!INTRODUZIDO POR ISILDA CM
    sfire%fgip(:,:)= 0.
    sfire%ischap(:,:)= 0.
    sfire%fz0(:,:)= 0.
    sfire%fwh(:,:)= 0.
    sfire%wz0(:,:)= 0. !INTRODUZIDO POR ISILDA CUNHA MENEZES 
    sfire%fgrnhfx(:, :)= 0.
    sfire%fgrnqfx(:, :)= 0.
    sfire%fcanhfx(:, :)= 0.
    sfire%fcanqfx(:, :)= 0.
    sfire%ros(:,:)= 0.
    sfire%flineint(:,:)= 0.
    sfire%flineint2(:,:)= 0.
    sfire%flineint_total(:,:)=0. !!!!INTRODUZIDO POR ISILDA CUNHA MENEZES
    sfire%f_ros0(:,:)= 0.
    sfire%f_rosx(:,:)= 0.
    sfire%f_rosy(:,:)= 0.
    sfire%f_ros(:,:)= 0.
    sfire%f_ros11(:,:)= 0.
    sfire%f_ros12(:,:)= 0.
    sfire%f_ros13(:,:)= 0.
    sfire%f_ros21(:,:)= 0.
    sfire%f_ros23(:,:)= 0.
    sfire%f_ros31(:,:)= 0.
    sfire%f_ros32(:,:)= 0.
    sfire%f_ros33(:,:)= 0.
    sfire%f_int(:,:)= 0.
    sfire%f_lineint(:,:)= 0.
    sfire%f_lineint2(:,:)= 0.
    sfire%fxlong(:, :)= 0.
    sfire%fxlat(:, :)= 0.
    sfire%fire_hfx(:, :)= 0.



    sfire%fmc_gc(:,:,:)= 0.
    sfire%fmc_equi(:,:,:)= 0.
    sfire%fmc_lag(:,:,:)= 0.



    sfire%xlong(:,:)= 0.
    sfire%xlat(:,:)= 0.    
    
    sfire%z0(:,:)= 0.
    sfire%ht(:,:)= 0.
    sfire%rain_old(:,:)= 0.
    sfire%t2(:,:)= 0.
    sfire%q2(:,:)= 0.
    sfire%psfc(:,:)= 0.
    sfire%rainc(:,:)= 0.
    sfire%rainnc(:,:)= 0.
    sfire%t2_old(:,:)= 0.
    sfire%q2_old(:,:)= 0.
    sfire%psfc_old(:,:)= 0.
    sfire%rh_fire(:,:)= 0.
    sfire%avg_fuel_frac(:, :)= 0.
    sfire%grnhfx(:, :)= 0.
    sfire%grnqfx(:, :)= 0.
    sfire%canhfx(:, :)= 0.
    sfire%canqfx(:, :)= 0.
    sfire%uah(:,:)= 0.
    sfire%vah(:,:)= 0.
    sfire%mut(:, :)= 0.

 

    sfire%u_2(:,:,:)= 0.
    sfire%v_2(:,:,:)= 0.
    sfire%ph_2(:,:,:)= 0.
    sfire%phb(:,:,:)= 0
    sfire%z_at_w(:, :, :)= 0.
    sfire%rthfrten(:, :, :)= 0.
    sfire%rqvfrten(:, :, :)= 0.
       
    sfire%fmep(:, :, :)= 0.  !!!!INTRODUZIDO POR ISILDA CUNHA MENEZES 

  END SUBROUTINE zero_sfire_brams

  subroutine StoreNamelistFileAtMem_sfire(oneNamelistFile)
    implicit none
    type(namelistFile), pointer :: oneNamelistFile
    sfire = oneNamelistFile%sfire

  end subroutine StoreNamelistFileAtMem_sfire


END MODULE mem_sfire
