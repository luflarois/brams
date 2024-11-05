MODULE ChemDryDepDriver


  USE rconstants, ONLY: &
       cpi,             &
       cpor,            &
       p00,             &
       g,               &
       vonk
       
  USE mem_grid, ONLY: &
       grid_g,        &
       jdim,          &
       dzt,           &
       zt,            &
       nzpmax,        &
       npatch,        &
       dtlt,          &
       imonth1,       &
       idate1,        &
       iyear1,        &
       ngrid

  USE micphys, ONLY: &
       level

  USE mem_cuparm, ONLY: &
       cuparm_g,        &
       nnqparm

  USE mem_basic, ONLY: &
       basic_g

  USE mem_turb, ONLY: &
       turb_g

  USE mem_leaf, ONLY: &
       leaf_g

  USE mem_micro, ONLY: &
       micro_g

  USE mem_radiate, ONLY: &
       radiate_g

  USE mem_chem1, ONLY: &
       chem1_g,        &
       chemistry
  
  USE mem_aer1, ONLY: &
       aerosol, &
       aer1_g

  USE module_dry_dep, ONLY: &
       dd_sedim,            &
       dry_dep                 ! Subroutine



  IMPLICIT NONE

  PRIVATE


  PUBLIC :: drydep_driver



CONTAINS


  !========================================================================
  SUBROUTINE drydep_driver(m1,m2,m3,ia,iz,ja,jz)

    INTEGER,              INTENT(IN)    :: m1
    INTEGER,              INTENT(IN)    :: m2
    INTEGER,              INTENT(IN)    :: m3
    INTEGER,              INTENT(IN)    :: ia
    INTEGER,              INTENT(IN)    :: iz
    INTEGER,              INTENT(IN)    :: ja
    INTEGER,              INTENT(IN)    :: jz

    REAL, POINTER :: conprr_dummy(:,:)    

       
    IF (associated(cuparm_g(ngrid)%conprr)) THEN
       conprr_dummy=>cuparm_g(ngrid)%conprr
    ELSE
       ALLOCATE(conprr_dummy(m2,m3))
       !LFR para testar NULLIFY(conprr_dummy)
    ENDIF

    CALL dry_dep(m1,m2,m3,ia,iz,ja,jz           & 
                ,cpi,cpor,p00,g,vonk            &
                ,jdim,dzt,zt,nzpmax,npatch,dtlt &
                ,level,nnqparm(ngrid)           &
                ,imonth1,idate1,iyear1          &
                ,chemistry,aerosol              &
                ,basic_g(ngrid)%theta	        &
                ,basic_g(ngrid)%rv	        &
                ,basic_g(ngrid)%pp	        &
                ,basic_g(ngrid)%dn0	        &
                ,basic_g(ngrid)%pi0	        &
                ,basic_g(ngrid)%up              &
                ,basic_g(ngrid)%vp   	        &
                ,turb_g(ngrid)%tkep	        &
                ,turb_g(ngrid)%sflux_t          &
                ,turb_g(ngrid)%sflux_r          &
                ,turb_g(ngrid)%sflux_u          &
                ,turb_g(ngrid)%sflux_v          &
                ,leaf_g(ngrid)%r_aer            &
                ,leaf_g(ngrid)%ustar	        &
                ,leaf_g(ngrid)%tstar	        &
                ,leaf_g(ngrid)%patch_area       &
                ,leaf_g(ngrid)%leaf_class       &
                ,leaf_g(ngrid)%patch_rough      & 
                ,micro_g(ngrid)%rcp	        & 
                ,micro_g(ngrid)%pcpg            & 
                ,grid_g(ngrid)%rtgt	        & 
                ,radiate_g(ngrid)%rshort        & 
!                ,cuparm_g(ngrid)%conprr         & 
                ,conprr_dummy                   & 
!-srf-27jan2015
!-    changed (:,:,ngrid) to (:,:,:) to avoid
!-    seg violation when AEROSOL=0
!               ,aer1_g(:,:,ngrid)             & 
!               ,chem1_g(:,ngrid)              & 
!               ,dd_sedim(:,ngrid))    
                ,aer1_g  (:,:,:)              & 
                ,chem1_g (:,:)                & 
                ,dd_sedim(:,:)                )    

    RETURN
  END SUBROUTINE drydep_driver
  !========================================================================


END MODULE ChemDryDepDriver
