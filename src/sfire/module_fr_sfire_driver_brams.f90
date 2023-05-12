module module_fr_sfire_driver_brams

   use module_fr_sfire_driver
   use module_fr_sfire_atm
!USE module_utility, only: WRFU_TimeInterval,WRFU_TimeIntervalGet, WRFU_SUCCESS
   implicit none

contains

   subroutine sfire_driver_em_init(grid, config_flags, time_step_start, dt)
      use module_domain_type
      use ModNamelistsfireFile
!ver isto
      !  USE module_configure , only : grid_config_rec_type
      implicit none

      type(domain), target          :: grid
      type(grid_config_rec_type), intent(IN)          :: config_flags

      integer :: &
         ids, ide, kds, kde, jds, jde &
         , ims, ime, kms, kme, jms, jme &
         , ips, ipe, kps, kpe, jps, jpe

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
      integer :: ifps
      integer :: ifpe
      integer :: kfps
      integer :: kfpe
      integer :: jfps
      integer :: jfpe

      real, intent(in) ::time_step_start, dt

      print *, 'sfire_driver_em_init: SFIRE initialization start'
      call flush (6)
      call message('sfire_driver_em_init: SFIRE initialization start')

      ids = config_flags%ids
      ide = config_flags%ide
      kds = config_flags%kds
      kde = config_flags%kde
      jds = config_flags%jds
      jde = config_flags%jde

      ims = config_flags%ims
      ime = config_flags%ime
      kms = config_flags%kms
      kme = config_flags%kme
      jms = config_flags%jms
      jme = config_flags%jme
      ips = config_flags%ips
      ipe = config_flags%ipe
      kps = config_flags%kps
      kpe = config_flags%kpe
      jps = config_flags%jps
      jpe = config_flags%jpe

      ifds = config_flags%ifds
      ifde = config_flags%ifde
      kfds = config_flags%kfds
      kfde = config_flags%kfde
      jfds = config_flags%jfds
      jfde = config_flags%jfde
      ifms = config_flags%ifms
      ifme = config_flags%ifme
      kfms = config_flags%kfms
      kfme = config_flags%kfme
      jfms = config_flags%jfms
      jfme = config_flags%jfme
      ifps = config_flags%ifps
      ifpe = config_flags%ifpe
      kfps = config_flags%kfps
      kfpe = config_flags%kfpe
      jfps = config_flags%jfps
      jfpe = config_flags%jfpe

      !  time_step_start=TimeInterval2Sec(domain_get_time_since_sim_start(grid))
      !  dt=TimeInterval2Sec(domain_get_time_step(grid))
      !print *, "estou no module_fr_sfire_driver_brams na rotina sfire_drive_em_init vou chamar sfire_driver_em"
      call sfire_driver_em(grid, config_flags &
                           , time_step_start, dt &
                           , ifun_beg, ifun_step - 1, 0 &
                           , ids, ide, kds, kde, jds, jde &
                           , ims, ime, kms, kme, jms, jme &
                           , ips, ipe, kps, kpe, jps, jpe &
                           , ifds, ifde, jfds, jfde &
                           , ifms, ifme, jfms, jfme &
                           , ifps, ifpe, jfps, jfpe &
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
         , ims, ime, kms, kme, jms, jme &
         , ips, ipe, kps, kpe, jps, jpe

      integer :: its, ite, jts, jte, kts, kte
      real, intent(in) ::time_step_start, dt
      integer::fire_time_step_ratio, itime_step, i, j
      real, dimension(config_flags%ips:config_flags%ipe, &
                      config_flags%jps:config_flags%jpe) :: grnhfx_save, grnqfx_save, &
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
      integer :: ifps
      integer :: ifpe
      integer :: kfps
      integer :: kfpe
      integer :: jfps
      integer :: jfpe
      !print *, "estou no module_fr_sfire_driver_brams na rotina sfire_drive_em_step vou chamar sfire_driver_em"
      flush (6)

      !print *, config_flags%ips, config_flags%ipe, config_flags%jps, &
      !   config_flags%jpe
      !flush (6)

      call message('sfire_driver_em_step: SFIRE step start')

      fire_time_step_ratio = config_flags%fire_time_step_ratio

      if (fire_time_step_ratio .lt. 1) then
         call crash('fire_time_step_ratio must be >= 1')
      end if

      !time_step_start=TimeInterval2Sec(domain_get_time_since_sim_start(grid))
      !dt=TimeInterval2Sec(domain_get_time_step(grid))/fire_time_step_ratio
      !time_step_start=(grid%itimestep+1)/2.

      ids = config_flags%ids
      ide = config_flags%ide
      kds = config_flags%kds
      kde = config_flags%kde
      jds = config_flags%jds
      jde = config_flags%jde

      ims = config_flags%ims
      ime = config_flags%ime
      kms = config_flags%kms
      kme = config_flags%kme
      jms = config_flags%jms
      jme = config_flags%jme
      ips = config_flags%ips
      ipe = config_flags%ipe
      kps = config_flags%kps
      kpe = config_flags%kpe
      jps = config_flags%jps
      jpe = config_flags%jpe

      ifds = config_flags%ifds
      ifde = config_flags%ifde
      kfds = config_flags%kfds
      kfde = config_flags%kfde
      jfds = config_flags%jfds
      jfde = config_flags%jfde
      ifms = config_flags%ifms
      ifme = config_flags%ifme
      kfms = config_flags%kfms
      kfme = config_flags%kfme
      jfms = config_flags%jfms
      jfme = config_flags%jfme
      ifps = config_flags%ifps
      ifpe = config_flags%ifpe
      kfps = config_flags%kfps
      kfpe = config_flags%kfpe
      jfps = config_flags%jfps
      jfpe = config_flags%jfpe

      grnhfx_save(:, :) = 0.
      grnqfx_save(:, :) = 0.
      canhfx_save(:, :) = 0.
      canqfx_save(:, :) = 0.

      !dt=0.5

      print *,'LFR-DBG: time_step_start, dt:', time_step_start, dt

      ipe1 = min(ipe, ide - 1)
      kpe1 = kpe - 1

      do itime_step = 1, fire_time_step_ratio

         call sfire_driver_em(grid, config_flags &
                              , time_step_start, dt &
                              , ifun_step, ifun_end, fire_test_steps &
                              , ids, ide, kds, kde, jds, jde &
                              , ims, ime, kms, kme, jms, jme &
                              , ips, ipe, kps, kpe, jps, jpe &
                              , ifds, ifde, jfds, jfde &
                              , ifms, ifme, jfms, jfme &
                              , ifps, ifpe, jfps, jfpe &
                              , grid%rho, grid%z_at_w, grid%dz8w &
                              )

        !print *, 'LFR-DBG: ', "config_flags%ips,config_flags%ipe:",config_flags%ips,config_flags%ipe
        !print *, 'LFR-DBG: ', "config_flags%jps,config_flags%jpe:",config_flags%jps,config_flags%jpe

         !print *, 'LFR-DBG: ', "jps,jpe1:", jps, jpe1
         !print *, 'LFR-DBG: ', "ips,ipe1:", ips, ipe1
         !print *, 'LFR-DBG: ', size(grnhfx_save, 1), size(grnhfx_save, 2)
         !print *, 'LFR-DBG: ', size(grid%grnhfx, 1), size(grid%grnhfx, 2)

         do j = jps, jpe !LFR 1
            do i = ips, ipe !LFR 1
               grnhfx_save(i, j) = grnhfx_save(i, j) + grid%grnhfx(i, j)
               grnqfx_save(i, j) = grnqfx_save(i, j) + grid%grnqfx(i, j)
               canhfx_save(i, j) = canhfx_save(i, j) + grid%canhfx(i, j)
               canqfx_save(i, j) = canqfx_save(i, j) + grid%canqfx(i, j)
            end do
         end do

         ! time_step_start=time_step_start+dt ESTE COMENTARIO E PARA TIRAR DEPOIS
      end do

      do j = jps, jpe !LFR 1
         do i = ips, ipe !LFR 1
            grid%grnhfx(i, j) = grnhfx_save(i, j)/fire_time_step_ratio
            grid%grnqfx(i, j) = grnqfx_save(i, j)/fire_time_step_ratio
            grid%canhfx(i, j) = canhfx_save(i, j)/fire_time_step_ratio
            grid%canqfx(i, j) = canqfx_save(i, j)/fire_time_step_ratio
         end do
      end do
      !print *, "estou no module_fr_sfire_driver_brams na rotina sfire_drive_em_step some o calor"

      call print_chsum(0,ims,ime,kms,kme,jms,jme,ids,ide,kds,kde,jds,jde,ips,ipe1,kps,kpe1,jps,jpe1,0,0,0,grid%z_at_w,'z_at_w')
      call print_chsum(0,ims,ime,kms,kme,jms,jme,ids,ide,kds,kde,jds,jde,ips,ipe1,kps,kpe1,jps,jpe1,0,0,0,grid%dz8w,'dz8w')
      call print_chsum(0,ims,ime,kms,kme,jms,jme,ids,ide,kds,kde,jds,jde,ips,ipe1,kps,kpe1,jps,jpe1,0,0,0,grid%rho,'rho')
      call print_chsum(0, ims, ime, 1, 1, jms, jme, ids, ide, 1, 1, jds, jde, ips, ipe1, 1, 1, jps, jpe1, 0, 0, 0, grid%mut, 'mu')
      call print_3d_stats(ips, ipe1, kps, kpe1, jps, jpe1, ims, ime, kms, kme, jms, jme, grid%rho, 'rho')
      call print_3d_stats(ips, ipe1, kps, kpe1, jps, jpe1, ims, ime, kms, kme, jms, jme, grid%z_at_w, 'z_at_w')
      call print_3d_stats(ips, ipe1, kps, kpe1, jps, jpe1, ims, ime, kms, kme, jms, jme, grid%dz8w, 'dz8w')

      do ij = 1, grid%num_tiles

         its = grid%i_start(ij)
         ite = min(grid%i_end(ij), ide - 1)
         jts = grid%j_start(ij)
         jte = min(grid%j_end(ij), jde - 1)
         kts = kds
         kte = kde
         !print *, "estou no module_fr_sfire_driver_brams na rotina sfire_drive_em_step vou chamar fire_tendency"
         call fire_tendency( &
            ids, ide - 1, kds, kde, jds, jde - 1, &
            ims, ime, kms, kme, jms, jme, &
            its, ite, kts, kte, jts, jte, &
            grid%grnhfx, grid%grnqfx, grid%canhfx, grid%canqfx, &
            config_flags%fire_ext_grnd, config_flags%fire_ext_crwn, config_flags%fire_crwn_hgt, &
            grid%ht, grid%z_at_w, grid%dz8w, grid%mut, grid%rho, &
            grid%rthfrten, grid%rqvfrten)

      end do

      write (msg, 993) lbound(grid%tracer, 4), ubound(grid%tracer, 4)
993   format('tracer array dimensions ', i3, ':', i3)
      call message(msg)
      write (msg, 994) num_tracer, config_flags%tracer_opt
994   format('number of tracers:', i3, ' tracer_opt=', i3)
      call message(msg)

      call print_3d_stats(its, ite, kts, kte, jts, jte, ims, ime, kms, kme, jms, jme, grid%rthfrten, 'fire_driver_phys:rthfrten')
      call print_3d_stats(its, ite, kts, kte, jts, jte, ims, ime, kms, kme, jms, jme, grid%rqvfrten, 'fire_driver_phys:rqvfrten')

      call message('sfire_driver_em_step: SFIRE step complete')
      !print *, "estou no module_fr_sfire_driver_brams e vou sair da rotina sfire_drive_em_step"

   end subroutine sfire_driver_em_step

!double precision function TimeInterval2Sec(time)

!    TYPE(WRFU_TimeInterval), intent(in) :: time

!    integer::rc,S,Sn,Sd

!    call WRFU_TimeIntervalGet(time,S=S,Sd=Sd,Sn=Sn,rc=rc)
!    if(rc.ne.WRFU_SUCCESS)call crash('TimeInterval2Sec: WRFU_TimeIntervalGet failed')

!    if(Sd.ne.0)then
!         TimeInterval2Sec=dble(S)+dble(Sn)/dble(Sd)
!    else
!         TimeInterval2Sec=dble(S)
!    endif
!end function TimeInterval2Sec

end module module_fr_sfire_driver_brams

