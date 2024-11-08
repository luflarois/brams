module module_fr_sfire_driver_brams

   use module_fr_sfire_driver
   use module_fr_sfire_atm
   implicit none

contains

   subroutine sfire_driver_em_init(grid, config_flags, time_step_start, dt)
   
      use module_domain_type
      use ModNamelistsfireFile
      
      implicit none

      type(domain), target          :: grid
      type(grid_config_rec_type), intent(IN)          :: config_flags

      integer :: ids, ide, kds, kde, jds, jde
      integer :: ims, ime, kms, kme, jms, jme
      integer :: ifds,ifde,kfds, kfde,jfds,jfde
      integer :: ifms,ifme,kfms,kfme,jfms,jfme   
      
      real, intent(in) ::time_step_start, dt

      print *, 'sfire_driver_em_init: SFIRE initialization start'
      call flush(6)

      ! Ajustar os limites da grelha principal e subgrelhas
      ids = config_flags%ids
      ide = config_flags%ide
      jds = config_flags%jds
      jde = config_flags%jde
      kds = config_flags%kds
      kde = config_flags%kde
      ims = config_flags%ims
      ime = config_flags%ime
      kms = config_flags%kms
      kme = config_flags%kme
      jms = config_flags%jms
      jme = config_flags%jme

      ! Ajustar os índices das subgrelhas (subgrids)
      ifds = config_flags%ifds
      ifde = config_flags%ifde
      jfds = config_flags%jfds
      jfde = config_flags%jfde
      kfds = config_flags%kfds
      kfde = config_flags%kfde
      ifms = config_flags%ifms
      ifme = config_flags%ifme
      kfms = config_flags%kfms
      kfme = config_flags%kfme
      jfms = config_flags%jfms
      jfme = config_flags%jfme
      
      ! Chamada à rotina principal, mantendo os índices principais e de subgrids
      call sfire_driver_em(grid, config_flags &
                          , time_step_start, dt &
                           , ifun_beg, ifun_step - 1, 0 &
                           , ids, ide, kds, kde, jds, jde &
                           , ims, ime, kms, kme, jms, jme &
                           , ifds, ifde, jfds, jfde &
                           , ifms, ifme, jfms, jfme &
                           )
                                                  
                           
      call message('sfire_driver_em_init: SFIRE initialization complete')
 
   end subroutine sfire_driver_em_init


   subroutine sfire_driver_em_step(grid, config_flags, time_step_start, dt)
   
      use module_domain_type
      use ModNamelistsfireFile
      use module_fr_sfire_util, only: fire_test_steps
      use module_state_description, only: num_tracer
 
      implicit none
      type(domain), target :: grid
      type(grid_config_rec_type), intent(in) :: config_flags

      integer :: &
         ids, ide, kds, kde, jds, jde &
         , ims, ime, kms, kme, jms, jme

      real, intent(in) ::time_step_start, dt
      integer::fire_time_step_ratio, itime_step, i, j
      real, dimension(config_flags%ids:config_flags%ide, &
                      config_flags%jds:config_flags%jde) :: grnhfx_save, grnqfx_save, &
                                            canhfx_save, canqfx_save
      integer:: ij, ipe1, jpe1, kpe1

      character(len=128)::msg
      integer :: ifds
      integer :: ifde
      integer :: kfds
      integer :: kfde
      integer :: jfds
      integer :: jfde
      integer :: ifms
      integer :: ifme
      integer :: kfms
      integer :: kfme
      integer :: jfms
      integer :: jfme
 
      call message('sfire_driver_em_step: SFIRE step start')

      fire_time_step_ratio = config_flags%fire_time_step_ratio

      if (fire_time_step_ratio .lt. 1) then
         call crash('fire_time_step_ratio must be >= 1')
      end if

      !time_step_start=TimeInterval2Sec(domain_get_time_since_sim_start(grid))
      !dt=TimeInterval2Sec(domain_get_time_step(grid))/fire_time_step_ratio
      !time_step_start=(grid%itimestep+1)/2.
      
      ! Ajustar os limites da grelha principal e subgrelhas
      ids = config_flags%ids
      ide = config_flags%ide
      jds = config_flags%jds
      jde = config_flags%jde
      kds = config_flags%kds
      kde = config_flags%kde
      ims = config_flags%ims
      ime = config_flags%ime
      kms = config_flags%kms
      kme = config_flags%kme
      jms = config_flags%jms
      jme = config_flags%jme

      ! Ajustar os índices das subgrelhas (subgrids)
      ifds = config_flags%ifds
      ifde = config_flags%ifde
      jfds = config_flags%jfds
      jfde = config_flags%jfde
      kfds = config_flags%kfds
      kfde = config_flags%kfde
      ifms = config_flags%ifms
      ifme = config_flags%ifme
      kfms = config_flags%kfms
      kfme = config_flags%kfme
      jfms = config_flags%jfms
      jfme = config_flags%jfme
      

      
      ! Inicializar os arrays de saída para armazenar os valores de calor e umidade
      grnhfx_save(:, :) = 0.
      grnqfx_save(:, :) = 0.
      canhfx_save(:, :) = 0.
      canqfx_save(:, :) = 0.
      
      
      do itime_step = 1, fire_time_step_ratio

      call sfire_driver_em(grid, config_flags &
                              , time_step_start, dt &
                              , ifun_step, ifun_end, fire_test_steps &
                              , ids, ide, kds, kde, jds, jde &
                              , ims, ime, kms, kme, jms, jme &
                              , ifds, ifde, jfds, jfde &
                              , ifms, ifme, jfms, jfme &
                              , grid%rho, grid%z_at_w, grid%dz8w &
                              )

  

      ! Loop principal sobre a grelha principal
      do j = jds, jde
         do i = ids, ide
            grnhfx_save(i, j) = grnhfx_save(i, j) + grid%grnhfx(i, j)
            grnqfx_save(i, j) = grnqfx_save(i, j) + grid%grnqfx(i, j)
            canhfx_save(i, j) = canhfx_save(i, j) + grid%canhfx(i, j)
            canqfx_save(i, j) = canqfx_save(i, j) + grid%canqfx(i, j)
         end do
      end do
     enddo
      ! Atualizar os valores no grid para os limites das subgrelhas
      do j = jfds, jfde
         do i = ifds, ifde
            grid%grnhfx(i, j) = grnhfx_save(i, j)/fire_time_step_ratio
            grid%grnqfx(i, j) = grnqfx_save(i, j)/fire_time_step_ratio
            grid%canhfx(i, j) = canhfx_save(i, j)/fire_time_step_ratio
            grid%canqfx(i, j) = canqfx_save(i, j)/fire_time_step_ratio
         end do
      end do
      
      call print_chsum(0,ims,ime,kms,kme,jms,jme,ids,ide,kds,kde,jds,jde,grid%z_at_w,'z_at_w')
      call print_chsum(0,ims,ime,kms,kme,jms,jme,ids,ide,kds,kde,jds,jde,grid%dz8w,'dz8w')
      call print_chsum(0,ims,ime,kms,kme,jms,jme,ids,ide,kds,kde,jds,jde,grid%rho,'rho')
      call print_chsum(0, ims, ime, 1, 1, jms, jme, ids, ide, 1, 1, jds, jde,grid%mut, 'mu')
      
      call fire_tendency(&
            ids, ide - 1, kds, kde, jds, jde - 1, &
            ims, ime, kms, kme, jms, jme, &
            grid%grnhfx, grid%grnqfx, grid%canhfx, grid%canqfx, &
            config_flags%fire_ext_grnd, config_flags%fire_ext_crwn, config_flags%fire_crwn_hgt, &
            grid%ht, grid%z_at_w, grid%dz8w, grid%mut, grid%rho, &
            grid%rthfrten, grid%rqvfrten)
            
      call print_3d_stats(ids, ide, kds, kde, jds, jde,ims, ime, kms, kme, jms, jme, grid%rthfrten, 'fire_driver_phys:rthfrten')
      call print_3d_stats(ids, ide, kds, kde, jds, jde,ims, ime, kms, kme, jms, jme, grid%rqvfrten, 'fire_driver_phys:rqvfrten')

      call message('sfire_driver_em_step: SFIRE step complete')
      
   end subroutine sfire_driver_em_step

end module module_fr_sfire_driver_brams