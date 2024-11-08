module sfireMod
!############################# Change Log ##################################
! Interface entre BRAMS e SFIRE
   use ReadBcst, only: &
      gatherData, Broadcast

   use node_mod, only: &
      mynum, &
      mchnum, &
      master_num, &
      nmachs, & !INTRODUZIDO ISILDA
      nodemxp, & !INTRODUZIDO ISILDA
      nodemyp, & !INTRODUZIDO ISILDA
      nxbeg, & !INTRODUZIDO ISILDA
      nxend, & !INTRODUZIDO ISILDA
      nybeg, & !INTRODUZIDO ISILDA
      nyend !INTRODUZIDO ISILDA

   use ParLib, only: parf_barrier
   private

   public sfclyr_sfire

contains

   subroutine sfclyr_sfire(mzp, mxp, myp, ia, iz, ja, jz)

      use ModNamelistsfireFile, only: &
         namelistsfireFile, &
         CreateNamelistsfireFile, &
         DestroyNamelistsfireFile, &
         GetNamelistsfireFileName, &
         ReadNamelistsfireFile, &
         DumpNamelistsfireFile

      use mem_grid, only: ngrid, deltaxn, deltax, deltayn, dtlongn, time, istp, &
                          nnzp, nnxp, nnyp, grid_g, npatch, zt, nzpmax ! dtlong, ESTES SO INTROD POR ISILDA

      use mem_sfire, only: &
         sfire_g, &
         alloc_sfire_brams, &
         nullify_sfire_brams, &
         dealloc_sfire_brams, &
         zero_sfire_brams, &
         config_flags

      use module_domain_type
      use module_fr_sfire_driver_brams
      use module_model_constants

      use, intrinsic :: ieee_arithmetic

      ! USE mem_cuparm, only: cuparm_g
      use mem_stilt, only: stilt_g
      use mem_micro, only: micro_g
      use mem_basic, only: basic_g
      use mem_leaf, only: leaf_g
      use mem_jules, only: jules_g
      use mem_varinit, only: varinit_g
      use mem_turb, only: turb_g
      use rconstants, only: cpi, p00, cpor
      use mem_basic, only: basic_g
      use micphys, only: mcphys_type

      use modSfire2Brams, only: readModisVeg, get_emission_in_global_brams_grid, nspecies_sfire &
                              , sfire_species_name, flam_frac, mean_frp, std_frp &
                              , mean_size, std_size, qsc, count_in_grid, sfire_info_area &
                              , sfire_info_time, sfire_info_frp, sfire_info_lat &
                              , sfire_info_lon,alloc_dealloc_sfire_info, valid_ij &
                              , comm_aer_data, testModSfire2Brams, test_only &
                              , fill_sfire_info_and_get_emissions


      implicit none

      integer, intent(IN) :: mzp, mxp, myp, ia, iz, ja, jz
      ! REAL, INTENT(IN) :: time
      real :: fdx, fdy

      !type(grid_config_rec_type), pointer :: config_flags => null()

      !Local Variables
      integer            :: ng, a_step, i, j, k, ierr, np
      real ::hf(nnxp(1), nnyp(1))
      real ::ch(nnxp(1), nnyp(1))
      real ::qf(nnxp(1), nnyp(1))
      real ::cq(nnxp(1), nnyp(1))
      real ::temp_sflux_t(nnxp(1), nnyp(1))
      real ::temp_sflux_r(nnxp(1), nnyp(1))
      ! real ::temp_aconpr(nnxp(1),nnyp(1))
      real ::temp_accpr(nnxp(1), nnyp(1))
      real ::temp_accpp(nnxp(1), nnyp(1))
      real ::temp_accps(nnxp(1), nnyp(1))
      real ::temp_accpa(nnxp(1), nnyp(1))
      real ::temp_accpg(nnxp(1), nnyp(1))
      real ::temp_accph(nnxp(1), nnyp(1))
      real ::temp_t2mj(nnxp(1), nnyp(1))
      real ::temp_topt(nnxp(1), nnyp(1))
      real ::temp_rv2mj(nnzp(1), nnxp(1), nnyp(1))
      real ::temp_glat(nnxp(1), nnyp(1))
      real ::temp_glon(nnxp(1), nnyp(1))
      real ::temp_veg_rough(nnxp(1), nnyp(1), npatch)
      real ::temp_veg_class(nnxp(1), nnyp(1), npatch)
      ! real ::zt(nzpmax)
      real ::temp_rtgt(nnxp(1), nnyp(1))
      real ::temp_h(nnzp(1), nnxp(1), nnyp(1))
!  real ::temp_dnp(nnzp(1),nnxp(1),nnyp(1))
      real ::temp_up(nnzp(1), nnxp(1), nnyp(1))
      real ::temp_vp(nnzp(1), nnxp(1), nnyp(1))
      real ::temp_dn0(nnzp(1), nnxp(1), nnyp(1))
      character(len=16)  :: varn
      real :: reduz
      real ::temp_pi0(nnzp(1), nnxp(1), nnyp(1))
      real ::temp_pp(nnzp(1), nnxp(1), nnyp(1))
      real ::temp_press(nnxp(1), nnyp(1))
      real ::temp_theta(nnzp(1), nnxp(1), nnyp(1))
      real ::picpi
      real ::a, b, step_isil, dxf,dyf
      integer, parameter :: ifm = 1 ! MUDAR ESTE PARAMETRO PARA UM LOOP SE EXISTIREM MAIS GRIDS

      if (mcphys_type /= 0) then

         print *, " SFIRE PRECISA DE MCPHYS_TYPE=0"
         stop
      end if

      varn = 'SFLUX_T'
      call gatherData(2, varn, 1, nnxp(1), nnyp(1), &
                      nmachs, mchnum, mynum, master_num, &
                      turb_g(ngrid)%sflux_t, temp_sflux_t)


      varn = 'SFLUX_R'
      call gatherData(2, varn, 1, nnxp(1), nnyp(1), &
                      nmachs, mchnum, mynum, master_num, &
                      turb_g(ngrid)%sflux_r, temp_sflux_r)

      varn = 'ACCPR'
      call gatherData(2, varn, 1, nnxp(1), nnyp(1), &
                      nmachs, mchnum, mynum, master_num, &
                      micro_g(ngrid)%accpr, temp_accpr)
      varn = 'ACCPP'
      call gatherData(2, varn, 1, nnxp(1), nnyp(1), &
                      nmachs, mchnum, mynum, master_num, &
                      micro_g(ngrid)%accpp, temp_accpp)
      varn = 'ACCPS'
      call gatherData(2, varn, 1, nnxp(1), nnyp(1), &
                      nmachs, mchnum, mynum, master_num, &
                      micro_g(ngrid)%accps, temp_accps)
      varn = 'ACCPA'
      call gatherData(2, varn, 1, nnxp(1), nnyp(1), &
                      nmachs, mchnum, mynum, master_num, &
                      micro_g(ngrid)%accpa, temp_accpa)
      varn = 'ACCPG'
      call gatherData(2, varn, 1, nnxp(1), nnyp(1), &
                      nmachs, mchnum, mynum, master_num, &
                      micro_g(ngrid)%accpg, temp_accpg)
      varn = 'ACCPH'
      call gatherData(2, varn, 1, nnxp(1), nnyp(1), &
                      nmachs, mchnum, mynum, master_num, &
                      micro_g(ngrid)%accph, temp_accph)
      varn = 'TOPT'
      call gatherData(2, varn, 1, nnxp(1), nnyp(1), &
                      nmachs, mchnum, mynum, master_num, &
                      grid_g(ngrid)%topt, temp_topt)
      varn = 'RV' !***ATENCAO acoplamos o leaf
      call gatherData(3, varn, 1, nnzp(1), nnxp(1), nnyp(1), &
                      nmachs, mchnum, mynum, master_num, &
                      basic_g(ngrid)%rv, temp_rv2mj)
      varn = 'PP'
      call gatherData(3, varn, 1, nnzp(1), nnxp(1), nnyp(1), &
                      nmachs, mchnum, mynum, master_num, &
                      basic_g(ngrid)%pp, temp_pp)
      varn = 'GLAT'
      call gatherData(2, varn, 1, nnxp(1), nnyp(1), &
                      nmachs, mchnum, mynum, master_num, &
                      grid_g(ngrid)%glat, temp_glat)
      varn = 'GLON'
      call gatherData(2, varn, 1, nnxp(1), nnyp(1), &
                      nmachs, mchnum, mynum, master_num, &
                      grid_g(ngrid)%glon, temp_glon)

      !======================================================================
      do np = 1, npatch
         varn = 'VEG_ROUGH'
         call gatherData(2, varn, 1, nnxp(1), nnyp(1), &
                         nmachs, mchnum, mynum, master_num, &
                         leaf_g(ngrid)%veg_rough(:, :, np), &
                         temp_veg_rough(:, :, np))
      end do
      ! print*,"passei VEG_ROUGH"
      do np = 1, npatch
         varn = 'LEAF_CLASS'
         call gatherData(2, varn, 1, nnxp(1), nnyp(1), &
                         nmachs, mchnum, mynum, master_num, &
                         leaf_g(ngrid)%leaf_class(:, :, np), &
                         temp_veg_class(:, :, np))
      end do

      varn = 'RTGT'
      call gatherData(2, varn, 1, nnxp(1), nnyp(1), &
                      nmachs, mchnum, mynum, master_num, &
                      grid_g(ngrid)%rtgt, temp_rtgt)
      varn = 'UP'
      call gatherData(3, varn, 1, nnzp(1), nnxp(1), nnyp(1), &
                      nmachs, mchnum, mynum, master_num, &
                      basic_g(ngrid)%up, temp_up)
      ! print*,"passei UP"
      varn = 'VP'
      call gatherData(3, varn, 1, nnzp(1), nnxp(1), nnyp(1), &
                      nmachs, mchnum, mynum, master_num, &
                      basic_g(ngrid)%vp, temp_vp)
      ! print*,"passei VP"
      varn = 'DN0'
      call gatherData(3, varn, 1, nnzp(1), nnxp(1), nnyp(1), &
                      nmachs, mchnum, mynum, master_num, &
                      basic_g(ngrid)%dn0, temp_dn0)
      !  print*, "passei DN0"
      varn = 'THETA'
      call gatherData(3, varn, 1, nnzp(1), nnxp(1), nnyp(1), &
                      nmachs, mchnum, mynum, master_num, &
                      basic_g(ngrid)%theta, temp_theta)
      ! print*,"passei THETA"
      varn = 'PI0'
      call gatherData(3, varn, 1, nnzp(1), nnxp(1), nnyp(1), &
                      nmachs, mchnum, mynum, master_num, &
                      basic_g(ngrid)%pi0, temp_pi0)

      if (mchnum == master_num) then
         
         print *, 'Inicializando a vegetacao MODIS e comuncação com o BRAMS'
         if (time ==0.0) call readModisVeg(temp_rtgt,deltax/1000.0) !Setup inicial  do BRAMS<->SFIRE
         if (.not. test_only) then

            do k = 1, nnzp(1)
               do j = 1, nnyp(1)
                  do i = 1, nnxp(1)
                     temp_h(k, i, j) = zt(k)*temp_rtgt(i, j)
                  end do
               end do
            end do

            k = 2
            !MFF - Troquei nnxp e nnyp de lugar
            do j = 1, nnyp(1)
               do i = 1, nnxp(1)
                  picpi = (temp_pi0(k, i, j) + temp_pp(k, i, j))*cpi
                  temp_press(i, j) = p00*picpi**cpor
                  temp_t2mj(i, j) = temp_theta(k, i, j)*picpi !ACOPLAMENTO DO LEAF
               end do
            end do

            !print *,'LFR-DBG 02: se for 1 timestep, inicialização: ', istp; call flush(6)
            if ((istp - 1) == 0) then

               print *, "LFR-DBG - FIST ISTP ", istp, 'Inicializando...' 
               call GlobalLatLonGrid(temp_glon, temp_glat)

               !print *, " Criando namelist do SFIRE"
               call flush (6)
               call CreateNamelistsfireFile(config_flags)

               !print *, " Pegando namelist do SFIRE"
               call flush (6)
               call GetNamelistsfireFileName(config_flags)

               !print *, " Lendo namelist do SFIRE"
               call flush (6)
               call ReadNamelistsfireFile(config_flags)

               write (*, "(50('-'),/,a)") 'Reading model control file (./sfire.in)...'

               call get_ijk_from_subgrid(config_flags)

               sfire_g%sr_x = config_flags%sr_x
               sfire_g%sr_y = config_flags%sr_y


               sfire_g%dx = deltaxn(ngrid)
               sfire_g%dy = deltayn(ngrid)

               sfire_g%itimestep = istp - 1
               sfire_g%num_tiles = 1
               sfire_g%num_tiles_spec = 1

               call nullify_sfire_brams(sfire_g)
               !print *,'LFR-DBG 04: Alocando sfrire_brams (sfire_g)- passando config_flags (size)'; call flush(6)
               call alloc_sfire_brams(sfire_g, config_flags%ims, config_flags%ime, config_flags%kms, config_flags%kme &
                                      , config_flags%jms, config_flags%jme, config_flags%ifms, config_flags%ifme &
                                      , config_flags%jfms, config_flags%jfme, config_flags%nfmc)

               ! Putting zero on all values
               !print *,'LFR-DBG 05: Zerando sfire_g que receberá o BRAMS'; call flush(6)
               call zero_sfire_brams(sfire_g)

               if (associated(sfire_g%i_start)) then; deallocate (sfire_g%i_start); nullify (sfire_g%i_start); end if
               if (associated(sfire_g%i_end)) then; deallocate (sfire_g%i_end); nullify (sfire_g%i_end); end if
               if (associated(sfire_g%j_start)) then; deallocate (sfire_g%j_start); nullify (sfire_g%j_start); end if
               if (associated(sfire_g%j_end)) then; deallocate (sfire_g%j_end); nullify (sfire_g%j_end); end if

               allocate (sfire_g%i_start(sfire_g%num_tiles))
               allocate (sfire_g%i_end(sfire_g%num_tiles))
               allocate (sfire_g%j_start(sfire_g%num_tiles))
               allocate (sfire_g%j_end(sfire_g%num_tiles))
            
                sfire_g%i_start(1) = config_flags%ids
                sfire_g%i_end(1) = config_flags%ide
                sfire_g%j_start(1) = config_flags%jds
                sfire_g%j_end(1) = config_flags%jde
            
               !print *,'LFR-DBG 07: Fazendo swap BRAMS - Leva do BRAMS para o sfire_g'; call flush(6)
               call swap_brams_sfire(sfire_g, & !temp_aconpr,&
                                     temp_accpr, temp_accpp, temp_accps, temp_accpa, &
                                     temp_accpg, temp_accph, temp_t2mj,temp_topt, &
                                     temp_rv2mj, temp_pp, temp_glat, temp_glon, &
                                     temp_veg_rough, temp_h, temp_dn0, temp_up, &
                                     temp_vp, temp_press)


               !print *, " Iniciando o SFIRE"
               !call flush (6)

               step_isil = (dtlongn(ngrid)*((config_flags%dx/1000.)*6.))/dtlongn(ngrid)

               !print *,'LFR-DBG 08: Chamando o sfire_driver_init para...?'; call flush(6)
               call sfire_driver_em_init(sfire_g, config_flags, time, step_isil)

               where (ieee_is_nan(sfire_g%rho)) sfire_g%rho = 0.00001
               where (sfire_g%rho <= 0.00000000000000000000001) sfire_g%rho = 0.00001

               where (ieee_is_nan(sfire_g%dz8w)) sfire_g%dz8w = 1.
               where (sfire_g%dz8w < 1.) sfire_g%dz8w = 1.

            end if !(time==0)

            !print *,'LFR-DBG 09: Marco!!! Passou pela inicialização do tempo ZERO !!! '; call flush(6)

            sfire_g%itimestep = istp
            print *, "sfire_g%itimestep =", sfire_g%itimestep
            !print *, " Convertendo indices das variaveis para o SFIRE - simu run"
            !call flush (6)

            ! get fire mesh dimensions
            !print *,'LFR-DBG 10: Pega ijk da subgrade fire'; call flush(6)

            call get_ijk_from_subgrid(config_flags)

            !print *, " sai do get_ijk"
            !call flush (6)

            sfire_g%sr_x = config_flags%sr_x
            sfire_g%sr_y = config_flags%sr_y
            print *, 'LFR-DBG: sr_x,sr_y:',sfire_g%sr_x ,sfire_g%sr_y

            sfire_g%dx = deltaxn(ngrid)
            sfire_g%dy = deltayn(ngrid)

            sfire_g%i_start(1) = config_flags%ids
            sfire_g%i_end(1) = config_flags%ide
            sfire_g%j_start(1) = config_flags%jds
            sfire_g%j_end(1) = config_flags%jde



            !print *,'LFR-DBG 11: Chamando o swap_brams -> fire'; call flush(6)
            call swap_B_brams_sfire(sfire_g, & !temp_aconpr,&
                                  temp_accpr, temp_accpp, temp_accps, temp_accpa, &
                                  temp_accpg, temp_accph, temp_t2mj,temp_topt, &
                                  temp_rv2mj, temp_pp, temp_glat, temp_glon, &
                                  temp_veg_rough, temp_h, temp_dn0, temp_up, &
                                  temp_vp, temp_press)

            !print *, " sai do swap_brams_sfire"
            !call flush (6)


            sfire_g%rainc = 0.

            where (ieee_is_nan(sfire_g%rho)) sfire_g%rho = 0.00001
            where (sfire_g%rho <= 0.00000000000000000000001) sfire_g%rho = 0.00001

            where (ieee_is_nan(sfire_g%dz8w)) sfire_g%dz8w = 1.
            where (sfire_g%dz8w < 1.) sfire_g%dz8w = 1.

            !print *, " Integrando o SFIRE "
            !call flush (6)

            ! Comentando apenas para testar a inicializacao do modulo e compilacao
            step_isil = (dtlongn(ngrid)*((config_flags%dx/1000.)*6.))/dtlongn(ngrid)
            print *,'LFR-DBG 13: iniciando a integração do sfire, step',step_isil, 'time: ',time; call flush(6)

            call sfire_driver_em_step(sfire_g, config_flags, time, step_isil)

        !   ###########################################################################
        !   passagem da area em fogo e do FRP

            do j=config_flags%jfds,config_flags%jfde
              do i=config_flags%ifds,config_flags%ifde
                  if(sfire_g%FRP(i,j) > 0.) then
                  write(51,*)'FIRE_AREA',sfire_g%fire_area(i,j)
                  write(52,*)'FRP',sfire_g%FRP(i,j)
                  endif
               enddo
            enddo 
        !   ############################################################################

            !passagem das unidades W/m^2 para K m/s

            hf(:, :) = 0.0
            ch(:, :) = 0.0
            do j = config_flags%jds, config_flags%jde
               do i = config_flags%ids, config_flags%ide
                  ! reduz = dn0(2,i,j) * 1004.
                  reduz = sfire_g%rho(i, 2, j)*1004
                  hf(i, j) = sfire_g%grnhfx(i, j)/reduz
                  ch(i, j) = sfire_g%canhfx(i, j)/reduz
                  !if(hf(i,j) > 0.)then
                  !print*,"ISILDA4 - h"
                  !print*,hf(i,j)
                  !endif
               end do
            end do

            !passagem das unidades de W/m^2 para Kg/kg m/s

            qf(:, :) = 0.0
            cq(:, :) = 0.0
            do j = config_flags%jds, config_flags%jde
               do i = config_flags%ids, config_flags%ide
                  !reduz = dn0(2,i,j) * 2.5e6
                  reduz = sfire_g%rho(i, 2, j)*2.5e6
                  qf(i, j) = sfire_g%grnqfx(i, j)/reduz
                  cq(i, j) = sfire_g%canqfx(i, j)/reduz
                !   if(sfire_g%grnqfx(i, j) > 0.)then
                !   print*,"ISILDA4 - h"
                !   print*,'IM-DBG:FORA:i,j,qf(i,j),sfire_g%grnqfx(i, j)= ',i,j,qf(i,j),sfire_g%grnqfx(i, j)
                !   endif
               end do
            end do

            !print *, 'LFR-DBG 14: somando ao calor do  Jules'; call flush (6)
            ! soma ao calor do  Jules

            do j = 1, config_flags%jds - 1
               do i = 1, nnxp(1)

                  temp_sflux_t(i, j) = temp_sflux_t(i, j)
                  temp_sflux_r(i, j) = temp_sflux_r(i, j)
               end do
            end do
            do j = config_flags%jde + 1, nnyp(1)
               do i = 1, nnxp(1)

                  temp_sflux_t(i, j) = temp_sflux_t(i, j)
                  temp_sflux_r(i, j) = temp_sflux_r(i, j)
               end do
            end do
            do j = config_flags%jds, config_flags%jde
               do i = 1, config_flags%ids - 1

                  temp_sflux_t(i, j) = temp_sflux_t(i, j)
                  temp_sflux_r(i, j) = temp_sflux_r(i, j)
               end do
            end do
            do j = config_flags%jds, config_flags%jde
               do i = config_flags%ide + 1, nnxp(1)

                  temp_sflux_t(i, j) = temp_sflux_t(i, j)
                  temp_sflux_r(i, j) = temp_sflux_r(i, j)
                  !LFR-DEB beg
                  if(istp>5) then
                     write (78,*) i,j,temp_sflux_t(i, j),hf(i, j),ch(i, j)
                  endif
                  !LFR-DEB end
               end do
            end do
            do j = config_flags%jds, config_flags%jde
               do i = config_flags%ids, config_flags%ide
                  temp_sflux_t(i, j) = temp_sflux_t(i, j) + hf(i, j) + ch(i, j)
                  !  print*,"temp_sflux_t",temp_sflux_t(i,j)
                  !  call flush(6)
                  temp_sflux_r(i, j) = temp_sflux_r(i, j) + qf(i, j) + cq(i, j)
                  !  print*,"temp_sflux_r",temp_sflux_r(i,j)
                  !  call flush(6)

               end do
            end do

         end if
      
         
         if (test_only) then 
            !Se for um teste preencha os dados de teste para enviar ao BRAMS
            call testModSfire2Brams(temp_glat,temp_glon,temp_dn0,time) 
         else
            !Caso contrário, preencher os valores de sfire_info e gerar emissões
            call fill_sfire_info_and_get_emissions(config_flags%ifds,config_flags%ifde &
                                ,config_flags%jfds,config_flags%jfde &
                                ,sfire_g%FRP,sfire_g%fire_area &
                                ,sfire_g%xlat,sfire_g%xlong &
                                ,sfire_g%tburn,temp_glat,temp_glon,temp_dn0)
         end if
         
         !print *, 'LFR-DBG 15: fim do setor serial'; call flush (6)
      end if
      !call Broadcast(mynum, master_num, "sync")
      !print *, 'LFR-DBG 16: ', mynum, 'Saindo do setor serial'; call flush (6)

      ! Comunica os pontos calculados para todos os processadores
      ! E determina aer1_g vars, valores médios e desvio padrão
      call comm_aer_data(mzp, mxp, myp, ia, iz, ja, jz)

      !call swap_sfire_brams(temp_sflux_t, temp_sflux_r, &
      !                      turb_g(ngrid)%sflux_t, turb_g(ngrid)%sflux_r)

      !  print*,"VOU SAIR DO SFIRE"
      return

   end subroutine sfclyr_sfire
   
   

   subroutine swap_brams_sfire(sfire, & !temp_aconpr,&
                               temp_accpr, temp_accpp, temp_accps, temp_accpa, &
                               temp_accpg, temp_accph, temp_t2mj,temp_topt, &
                               temp_rv2mj, temp_pp, temp_glat, temp_glon, &
                               temp_veg_rough, temp_h, temp_dn0, temp_up, &
                               temp_vp, temp_press)

      use module_domain_type
      use mem_grid, only: nnzp, nnxp, nnyp, &
                          polelat, polelon, npatch

      implicit none

      !real, INTENT(IN) ::temp_aconpr(nnxp(1),nnyp(1))
      real, intent(IN) ::temp_accpr(nnxp(1), nnyp(1))
      real, intent(IN) ::temp_accpp(nnxp(1), nnyp(1))
      real, intent(IN) ::temp_accps(nnxp(1), nnyp(1))
      real, intent(IN) ::temp_accpa(nnxp(1), nnyp(1))
      real, intent(IN) ::temp_accpg(nnxp(1), nnyp(1))
      real, intent(IN) ::temp_accph(nnxp(1), nnyp(1))
      real, intent(IN) ::temp_t2mj(nnxp(1), nnyp(1))
      real, intent(IN) ::temp_topt(nnxp(1), nnyp(1))
      real, intent(IN) ::temp_rv2mj(nnzp(1), nnxp(1), nnyp(1))
      real, intent(IN) ::temp_pp(nnzp(1), nnxp(1), nnyp(1))
      real, intent(IN) ::temp_glat(nnxp(1), nnyp(1))
      real, intent(IN) ::temp_glon(nnxp(1), nnyp(1))
      real, intent(IN) ::temp_veg_rough(nnxp(1), nnyp(1), npatch)
      real, intent(IN) ::temp_h(nnzp(1), nnxp(1), nnyp(1))
      real, intent(IN) ::temp_dn0(nnzp(1), nnxp(1), nnyp(1))
      real, intent(IN) ::temp_up(nnzp(1), nnxp(1), nnyp(1))
      real, intent(IN) ::temp_vp(nnzp(1), nnxp(1), nnyp(1))
      real, intent(IN) ::temp_press(nnxp(1), nnyp(1))
      real :: x, y

      type(domain), target :: sfire
      integer :: k
      real, dimension(nnzp(1), nnxp(1), nnyp(1)) :: temp

      sfire%rainc(:, :) = 0.0
      sfire%rainnc(:, :) = sfire%rainc(:, :) + temp_accpr(:, :) &
                           + temp_accpp(:, :) + temp_accps(:, :) &
                           + temp_accpa(:, :) + temp_accpg(:, :) &
                           + temp_accph(:, :)

      sfire%t2(:, :) = temp_t2mj(:, :)

      sfire%ht(:, :) = temp_topt(:, :)
      sfire%zsf(:, :) = 0.
      sfire%dzdxf(:,:) = 0.
      sfire%dzdyf(:,:) = 0.
      sfire%nfuel_cat(:,:) = 14
      sfire%q2(:, :) = temp_rv2mj(2, :, :)
      sfire%mut(:, :) = temp_pp(1, :, :)
      sfire%xlat(:, :) = temp_glat(:, :)
      sfire%xlong(:, :) = temp_glon(:, :)
      sfire%z0(:, :) = temp_veg_rough(:, :, 2)
      sfire%psfc(:, :) = temp_press(:, :)

      do k = 1, nnzp(1) - 1
         temp(k, :, :) = temp_h(k + 1, :, :) - temp_h(k, :, :)
      end do

      call swap_kij_to_ikj(sfire%rho, temp_dn0)
      call swap_kij_to_ikj(sfire%dz8w, temp)
      call swap_kij_to_ikj(sfire%z_at_w, temp_h)

      sfire%ph_2(:, :, :) = sfire%z_at_w(:, :, :)
      sfire%phb(:, :, :) = sfire%z_at_w(:, :, :)

      call swap_kij_to_ikj(sfire%u_2, temp_up)
      call swap_kij_to_ikj(sfire%v_2, temp_vp)

      return
   end subroutine swap_brams_sfire

subroutine swap_B_brams_sfire(sfire, & !temp_aconpr,&
                               temp_accpr, temp_accpp, temp_accps, temp_accpa, &
                               temp_accpg, temp_accph, temp_t2mj, temp_topt, &
                               temp_rv2mj, temp_pp, temp_glat, temp_glon, &
                               temp_veg_rough, temp_h, temp_dn0, temp_up, &
                               temp_vp, temp_press)

      use module_domain_type
      use mem_grid, only: nnzp, nnxp, nnyp, &
                          polelat, polelon, npatch

      implicit none

      !real, INTENT(IN) ::temp_aconpr(nnxp(1),nnyp(1))
      real, intent(IN) ::temp_accpr(nnxp(1), nnyp(1))
      real, intent(IN) ::temp_accpp(nnxp(1), nnyp(1))
      real, intent(IN) ::temp_accps(nnxp(1), nnyp(1))
      real, intent(IN) ::temp_accpa(nnxp(1), nnyp(1))
      real, intent(IN) ::temp_accpg(nnxp(1), nnyp(1))
      real, intent(IN) ::temp_accph(nnxp(1), nnyp(1))
      real, intent(IN) ::temp_t2mj(nnxp(1), nnyp(1))
      real, intent(IN) ::temp_topt(nnxp(1), nnyp(1))
      real, intent(IN) ::temp_rv2mj(nnzp(1), nnxp(1), nnyp(1))
      real, intent(IN) ::temp_pp(nnzp(1), nnxp(1), nnyp(1))
      real, intent(IN) ::temp_glat(nnxp(1), nnyp(1))
      real, intent(IN) ::temp_glon(nnxp(1), nnyp(1))
      real, intent(IN) ::temp_veg_rough(nnxp(1), nnyp(1), npatch)
      real, intent(IN) ::temp_h(nnzp(1), nnxp(1), nnyp(1))
      real, intent(IN) ::temp_dn0(nnzp(1), nnxp(1), nnyp(1))
      real, intent(IN) ::temp_up(nnzp(1), nnxp(1), nnyp(1))
      real, intent(IN) ::temp_vp(nnzp(1), nnxp(1), nnyp(1))
      real, intent(IN) ::temp_press(nnxp(1), nnyp(1))
      real :: x, y

      type(domain), target :: sfire
      integer :: k
      real, dimension(nnzp(1), nnxp(1), nnyp(1)) :: temp

      sfire%rainc(:, :) = 0.0
      sfire%rainnc(:, :) = sfire%rainc(:, :) + temp_accpr(:, :) &
                           + temp_accpp(:, :) + temp_accps(:, :) &
                           + temp_accpa(:, :) + temp_accpg(:, :) &
                           + temp_accph(:, :)

      sfire%t2(:, :) = temp_t2mj(:, :)
      sfire%ht(:, :) = temp_topt(:, :)
      sfire%q2(:, :) = temp_rv2mj(2, :, :)
      sfire%mut(:, :) = temp_pp(1, :, :)
      sfire%xlat(:, :) = temp_glat(:, :)
      sfire%xlong(:, :) = temp_glon(:, :)
      sfire%z0(:, :) = temp_veg_rough(:, :, 2)
      sfire%psfc(:, :) = temp_press(:, :)

      do k = 1, nnzp(1) - 1
         temp(k, :, :) = temp_h(k + 1, :, :) - temp_h(k, :, :)
      end do

      call swap_kij_to_ikj(sfire%rho, temp_dn0)
      call swap_kij_to_ikj(sfire%dz8w, temp)
      call swap_kij_to_ikj(sfire%z_at_w, temp_h)

      sfire%ph_2(:, :, :) = sfire%z_at_w(:, :, :)
      sfire%phb(:, :, :) = sfire%z_at_w(:, :, :)
      call swap_kij_to_ikj(sfire%u_2, temp_up)
      call swap_kij_to_ikj(sfire%v_2, temp_vp)

      return
   end subroutine swap_B_brams_sfire

   subroutine swap_sfire_brams(temp_sflux_t, temp_sflux_r, &
                               said_sflux_t, said_sflux_r)

      use mem_grid, only: nnxp, nnyp

      use node_mod, only: &
         mynum, &
         master_num, &
         nodemxp, &
         nodemyp, &
         ia, iz, &
         ja, jz
      use ParLib, only: &
         parf_bcast

      implicit none
      real, intent(IN) ::temp_sflux_t(nnxp(1), nnyp(1))
      real, intent(IN) ::temp_sflux_r(nnxp(1), nnyp(1))
      ! real, INtENT(OUT) ::said_sflux_t(nnxp(1), nnyp(1))
      real, intent(OUT) ::said_sflux_t(nodemxp(mynum, 1), nodemyp(mynum, 1))
      ! real, INTENT(OUT) ::said_sflux_r(nnxp(1), nnyp(1))
      real, intent(OUT) ::said_sflux_r(nodemxp(mynum, 1), nodemyp(mynum, 1))
      integer:: n1, n2

      n1 = nnxp(1)
      n2 = nnyp(1)
      said_sflux_t(:, :) = temp_sflux_t(:, :)

      said_sflux_r(:, :) = temp_sflux_r(:, :)
      return

   end subroutine swap_sfire_brams

   subroutine swap_kij_to_ikj(a, b)
      ! USE node_mod, only: mxp, myp, mzp
      use mem_grid, only: nnzp, nnxp, nnyp, ngrid
      implicit none

      integer :: i, j, k
      real :: b(nnzp(ngrid), nnxp(ngrid), nnyp(ngrid)), a(nnxp(ngrid), nnzp(ngrid), nnyp(ngrid))

      do j = 1, nnyp(ngrid)
         do k = 1, nnzp(ngrid)
            do i = 1, nnxp(ngrid)
               a(i, k, j) = b(k, i, j)
            end do
         end do
      end do

   end subroutine swap_kij_to_ikj



   subroutine get_ijk_from_subgrid(config_flags)
      use mem_grid, only: nnzp, nnxp, nnyp
      use ModNamelistsfireFile

      type(grid_config_rec_type), target :: config_flags

      !- converting WRF setting to BRAMS
      !alterei o indice inicial ids e jds
      config_flags%ids = 3; config_flags%ide = nnxp(1) - 2; config_flags%jds = 3; config_flags%jde = nnyp(1) - 2
      config_flags%kds = 2; config_flags%kde = nnzp(1)
      config_flags%ims = 1; config_flags%ime = nnxp(1); config_flags%jms = 1; config_flags%jme = nnyp(1); config_flags%kms = 1
      config_flags%kme = nnzp(1)

      config_flags%ifds = config_flags%ids
      config_flags%ifde = config_flags%ide*config_flags%sr_x
      config_flags%jfds = config_flags%jds
      config_flags%jfde = config_flags%jde*config_flags%sr_y
      config_flags%kfds = config_flags%kds
      config_flags%kfde = config_flags%kde

      config_flags%ifms = (config_flags%ims - 1)*config_flags%sr_x + 1
      config_flags%ifme = config_flags%ime*config_flags%sr_x
      config_flags%jfms = (config_flags%jms - 1)*config_flags%sr_y + 1
      config_flags%jfme = config_flags%jme*config_flags%sr_y
      config_flags%kfms = config_flags%kms
      config_flags%kfme = config_flags%kme
      
         
      return
   end subroutine get_ijk_from_subgrid

   subroutine GlobalLatLonGrid(temp_glon, temp_glat)

      use mem_grid, only: nnxp, nnyp, ngrid

      implicit none

      integer :: i
      integer :: ierr
      real, intent(IN) ::temp_glat(nnxp(ngrid), nnyp(ngrid))
      real, intent(IN) ::temp_glon(nnxp(ngrid), nnyp(ngrid))
      real :: delLon, delLat
      real :: deltaSum
      integer :: nLon, nLat
      real :: firstLon, firstLat
      real :: latMin, latMax, latMinBrams, latMaxBrams
      real :: lonMin, lonMax, lonMinBrams, lonMaxBrams
      real :: loni(ngrid)
      real :: lonf(ngrid)
      real :: lati(ngrid)
      real :: latf(ngrid)
      character(len=8) :: c0
      character(len=*), parameter :: h = "**(GlobalLatLonGrid)**"
      logical, parameter :: dumpLocal = .true.
      logical, parameter :: project = .true.
      real, allocatable :: lon(:) ! longitudes
      real, allocatable :: lat(:) ! longitudes

      print *, "DEBUG GlobalLatLonGrid :: ", h

      deltaSum = sum &
                 ( &
                 (temp_glon(nnxp(ngrid), :) - temp_glon(1, :)) &
                 /(nnxp(ngrid) - 1) &
                 )
      delLon = deltaSum/nnyp(ngrid)
      print *, "DEBUG GlobalLatLonGrid - delLon ::", delLon
      ! same procedure for latitudes:
      ! sum of average delta latitude over all longitudes
      ! divided by number of longitudes

      deltaSum = sum &
                 ( &
                 (temp_glat(:, nnyp(ngrid)) - temp_glat(:, 1)) &
                 /(nnyp(ngrid) - 1) &
                 )
      delLat = deltaSum/nnxp(ngrid)
      print *, "DEBUG GlobalLatLonGrid - delLat ::", delLat
      ! Step 2: Compute origin and number of points, restricted to
      ! user selected area and projection  specified at namelist file

      if (project) then

         ! if projection:
         !   find BRAMS grid extremes
         !   intersect with user selected area
         !   define origin and number of points keeping original delta

         ! envelope BRAMS grid area at earth's surface

         lonMinBrams = minval(temp_glon)
         lonMaxBrams = maxval(temp_glon)
         latMinBrams = minval(temp_glat)
         latMaxBrams = maxval(temp_glat)

         print *, "DEBUG:: lonMinBrams = minval(temp_glon)", lonMinBrams
         print *, "DEBUG:: lonMaxBrams = maxval(temp_glon)", lonMaxBrams
         print *, "DEBUG:: latMinBrams = minval(temp_glat)", latMinBrams
         print *, "DEBUG:: latMaxBrams = maxval(temp_glat)", latMaxBrams

         ! intersection of BRAMS envelope and user defined area
         loni(ngrid) = -180.
         lonf(ngrid) = 180.
         lati(ngrid) = -90.
         latf(ngrid) = 90.

         !ATENCAO::: VER SE MUDAMOS MAIS TARDE PARA ISTO::oneNamelistFile%loni(ngrid)), oneNamelistFile%lonf(ngrid)),oneNamelistFile%lati(ngrid)), oneNamelistFile%latf(ngrid))
         lonMin = max(lonMinBrams, loni(ngrid))
         lonMax = min(lonMaxBrams, lonf(ngrid))
         latMin = max(latMinBrams, lati(ngrid))
         latMax = min(latMaxBrams, latf(ngrid))
         print *, "DEBUG:: lonMin =", lonMin
         print *, "DEBUG:: lonMax =", lonMax
         print *, "DEBUG:: latMin =", latMin
         print *, "DEBUG:: latMax =", latMax
         ! post grid origin

         firstLon = lonMin
         firstLat = latMin

         ! post grid number of points (#intervals + 1)

         nLon = 1 + &
                ceiling((lonMax - lonMin)/delLon)
         nLat = 1 + &
                ceiling((latMax - latMin)/delLat)

         print *, "nlon,nlat - project", nLon, nLat, project
      else

         ! if no projection:
         !   find BRAMS grid origin as average first latitude and average first longitude
         !   compute last latitude and last longitude keeping delta and # points
         !   intersect with user selected area
         !   define origin and number of points keeping original delta

         ! intersection of BRAMS envelope and user defined area
         loni(ngrid) = -180.
         lonf(ngrid) = 180.
         lati(ngrid) = -90.
         latf(ngrid) = 90.

         lonMinBrams = sum(temp_glon(1, :))/nnyp(ngrid)
         latMinBrams = sum(temp_glat(:, 1))/nnxp(ngrid)

         lonMaxBrams = lonMinBrams + delLon*real(nnxp(ngrid) - 1)
         latMaxBrams = latMinBrams + delLat*real(nnyp(ngrid) - 1)

         !ATENCAO::: VER SE MUDAMOS MAIS TARDE PARA ISTO::oneNamelistFile%loni(ngrid)), oneNamelistFile%lonf(ngrid)),oneNamelistFile%lati(ngrid)), oneNamelistFile%latf(ngrid))
         lonMin = max(lonMinBrams, loni(ngrid))
         lonMax = min(lonMaxBrams, lonf(ngrid))
         latMin = max(latMinBrams, lati(ngrid))
         latMax = min(latMaxBrams, latf(ngrid))
         firstLon = lonMin
         firstLat = latMin
         print *, "lonMin, lonMax,latMin,latMax =", lonMin, lonMax, latMin, latMax
         print *, "firstLon, firstLat =", firstLon, firstLat
         if ((lonMin == lonMinBrams) .and. &
             (lonMax == lonMaxBrams) .and. &
             (latMin == latMinBrams) .and. &
             (latMax == latMaxBrams)) then

            nLon = nnxp(ngrid)
            nLat = nnyp(ngrid)
            print *, "vim pelo NNXP E NNYP coorde"
         else
            print *, "VOU PELO CEILING"
            nLon = 1 + &
                   ceiling((lonMax - lonMin)/delLon)
            nLat = 1 + &
                   ceiling((latMax - latMin)/delLat)
            print *, "nLon ,nLat =", nLon, nLat
         end if
      end if

      ! Step 3: Compute longitude and latitude points
      ! given first point, delta and number of points

      allocate (lon(nLon), stat=ierr)
      if (ierr /= 0) then
         write (c0, "(i8)") nLon
         call fatal_error(h//" allocate lon("// &
                          trim(adjustl(c0))//") fails")
      end if

      call Axis(nLon, firstLon, &
                delLon, lon)

      allocate (lat(nLat), stat=ierr)
      if (ierr /= 0) then
         write (c0, "(i8)") nLat
         call fatal_error(h//" allocate lat("// &
                          trim(adjustl(c0))//") fails")
      end if

      call Axis(nLat, firstLat, &
                delLat, lat)

      if (dumpLocal) then
         open (15, file='coordenadas.dat', form='formatted', status='replace')
         write (15, "(a6,i8)") "nLon =", nLon
         write (15, "(a6,i8)") "nLat =", nLat
         write (15, "(a10,f12.7)") "firstLon =", firstLon
         write (15, "(f12.7)") lon
         write (15, "(a8,f12.7)") "delLon =", delLon
         write (15, "(a10,f12.7)") "firstLat =", firstLat
         write (15, "(f12.7)") lat
         write (15, "(a8,f12.7)") "delLat =", delLat
         close (15)
         !call MsgDump (h//" lat-lon has size"//&
         !     " ("//trim(adjustl(c0))//","//trim(adjustl(c1))//")"//&
         !     " with lon,lat = "//&
         !     " ("//trim(adjustl(c2))//":"//trim(adjustl(c3))//":"//trim(adjustl(c4))//", "//&
         !     trim(adjustl(c5))//":"//trim(adjustl(c6))//":"//trim(adjustl(c7))//")")

      end if
   end subroutine GlobalLatLonGrid

   subroutine Axis(nVal, firstVal, delVal, val)
      integer, intent(in) :: nVal
      real, intent(in) :: firstVal
      real, intent(in) :: delVal
      real, intent(out) :: val(nVal)

      integer :: i

      print *, 'DEBUG :: func axis'
      do i = 1, nVal
         val(i) = firstVal + real(i - 1)*delVal
      end do

      return
   end subroutine Axis

end module sfireMod

