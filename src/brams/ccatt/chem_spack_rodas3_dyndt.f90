module mod_chem_spack_rodas3_dyndt

   use mem_chem1, only: &
      chem1_vars        ! Type

   use mem_spack, only: &
      spack_type, & ! Type
      spack_type_2d     ! Type

   use Extras, only: &
      ext3d             ! Type

!  USE solve_sparse, ONLY: &
!       Solve_linear      & ! Subroutine, OUT: spack_2d(ijk,inob)%DLk1

   use mod_chem_spack_ros, only: &
      get_always_zeros, &  ! Subroutine
      get_number_nonzeros  ! Subroutine

   use mod_chem_spack_jacdchemdc, only: &
      jacdchemdc        ! Subroutine

   use mod_chem_spack_kinetic, only: &
      kinetic           ! Subroutine

   use mod_chem_spack_fexchem, only: &
      fexchem           ! Subroutine

   use node_mod, only: &
          mynum

   implicit none

   private

   public chem_rodas3_dyndt

contains

   !========================================================================================
   subroutine chem_rodas3_dyndt(m1, m2, m3, dtlt, pp, pi0, theta, rv, dn0, cosz, rcp, &
                                nspecies, nr, nr_photo, weight, PhotojMethod, maxnspecies, &
                                nspecies_chem_transported, nspecies_chem_no_transported, &
                                transp_chem_index, no_transp_chem_index, chem1_g, nob, maxblock_size, &
                                block_end, indexk, indexi, indexj, kij_index, last_accepted_dt, &
                                spack, spack_2d, atol, rtol, get_non_zeros, &
                                N_DYN_CHEM, split_method, cp, cpor, p00, na_extra3d, extra3d, &
                                jphoto, att)

      use mem_grid, only: time
      implicit none

      integer, intent(IN) :: m1
      integer, intent(IN) :: m2
      integer, intent(IN) :: m3

      ! mem_grid
      real, intent(IN) :: dtlt

      ! mem_basic
      real, intent(IN) :: pp(m1, m2, m3)
      real, intent(IN) :: pi0(m1, m2, m3)
      real, intent(IN) :: theta(m1, m2, m3)
      real, intent(IN) :: rv(m1, m2, m3)
      real, intent(IN) :: dn0(m1, m2, m3)

      ! mem_radiate
      real, intent(IN) :: cosz(m2, m3)

      ! mem_micro
      real, intent(IN) :: rcp(m1, m2, m3)

      ! chem1_list
      integer, intent(IN) :: nspecies
      integer, intent(IN) :: nr_photo
      integer, intent(IN) :: nr
      real, intent(IN) :: weight(nspecies)
      character(LEN=10), intent(IN) :: PhotojMethod

      ! mem_chem1
      integer, intent(IN) :: maxnspecies
      integer, intent(IN) :: nspecies_chem_transported
      integer, intent(IN) :: nspecies_chem_no_transported
      integer, intent(IN) :: transp_chem_index(maxnspecies)
      integer, intent(IN) :: no_transp_chem_index(maxnspecies)
      character(LEN=20), intent(IN) :: split_method
      integer, intent(IN) :: N_DYN_CHEM

      ! mem_chem1
      type(chem1_vars), intent(INOUT) :: chem1_g(nspecies)

      ! spack_utils
      integer, intent(IN)    :: nob
      integer, intent(IN)    :: maxblock_size
      integer, intent(IN)    :: block_end(:)    !ngrids
      integer, intent(IN)    :: indexk(:, :)    !block_size,nob(ngr))
      integer, intent(IN)    :: indexi(:, :)    !block_size,nob(ngr))
      integer, intent(IN)    :: indexj(:, :)    !block_size,nob(ngr))
      integer, intent(IN)    :: kij_index(:, :) !(maxblock_size,nob(ngr))
      double precision, intent(INOUT) :: last_accepted_dt (:) !(nob(ngr))
!    INTEGER         ,   INTENT(IN)    :: block_end(nob)
!    INTEGER         ,   INTENT(IN)    :: indexk(maxblock_size,nob)
!    INTEGER         ,   INTENT(IN)    :: indexi(maxblock_size,nob)
!    INTEGER         ,   INTENT(IN)    :: indexj(maxblock_size,nob)
!    INTEGER         ,   INTENT(IN)    :: kij_index(maxblock_size,nob)
!    DOUBLE PRECISION,   INTENT(INOUT) :: last_accepted_dt(nob)

      ! mem_spack
      type(spack_type), intent(INOUT) :: spack(1)
!    TYPE(spack_type_2d), INTENT(INOUT) :: spack_2d(maxblock_size,1)
      type(spack_type_2d), intent(INOUT) :: spack_2d(:, :)
      double precision, intent(IN)    :: atol(nspecies)
      double precision, intent(IN)    :: rtol(nspecies)

      logical, intent(INOUT) :: get_non_zeros

      ! rconstants
      real, intent(IN) :: cp
      real, intent(IN) :: cpor
      real, intent(IN) :: p00

      ! extra
      integer, intent(IN)    :: na_extra3d
      type(ext3d), intent(INOUT) :: extra3d(na_extra3d) ! Output JNO2

      ! FastJX
      double precision, intent(IN) :: jphoto(nr_photo, m1, m2, m3)

      ! uv_atten
      real, intent(IN) :: att(m2, m3)

!  INTEGER,PARAMETER :: block_size=2*32
!  INTEGER,PARAMETER :: block_size=4
!  INTEGER,PARAMETER :: block_size=1

      double precision, parameter         :: pmar = 28.96d0
      double precision, parameter         :: Threshold1 = 0.0d0
      double precision, parameter         :: Threshold2 = -1.d-1

      double precision, parameter         :: Igamma = 0.5d0
      double precision, parameter         :: c43 = 4.d0/3.d0
      double precision, parameter         :: c83 = 8.d0/3.d0
      double precision, parameter         :: c56 = 5.d0/6.d0
      double precision, parameter         :: c16 = 1.d0/6.d0
      double precision, parameter         :: c112 = 1.d0/12.d0
      double precision dble_dtlt_i, fxc, Igamma_dtstep, dt_chem, dt_chem_i
      integer :: i, ijk, n, j, k, ispc, Ji, Jj, k_, i_, j_, kij_, kij

      !integer,parameter ::  NOB_MEM=1   ! if 1 : alloc just one block (uses less memory)
      !                                  ! if 0 : alloc all blocks     (uses more memory)

      !- default number of allocatable blocks (re-use scratch arrays spack()% )
      integer, parameter :: inob = 1

      !- parameters for dynamic timestep control
      integer all_accepted
      double precision :: dt_min, dt_max, dt_actual, dt_new
      double precision :: time_f, time_c
      double precision, parameter :: FacMin = 0.2d0 ! lower bound on step decrease factor (default=0.2)
      double precision, parameter :: FacMax = 6.0d0 ! upper bound on step increase factor (default=6)
      double precision, parameter :: FacRej = 0.1d0 ! step decrease factor after multiple rejections
      double precision, parameter :: FacSafe = 0.9d0 ! step by which the new step is slightly smaller
      ! than the predicted value  (default=0.9)
      double precision, parameter :: uround = 1.d-15, elo = 3.0d0, Roundoff = 1.d-8
      double precision Fac, Tol, err1, max_err1

      double precision, allocatable :: DLmat(:, :, :)
      integer, allocatable :: ipos(:)
      integer, allocatable :: jpos(:)

      integer :: NumberOfNonZeros, nz
      integer :: offset, offsetDLmat, offsetnz
      integer :: blocksize, sizeOfMatrix, maxnonzeros, blocknonzeros
      real :: start, finish
      real :: elapsed_time_solver, elapsed_time, elapsed_time_alloc, elapsed_time_dealloc, elapsed_time_copy

      integer(KIND=8) :: matrix_Id
      integer :: error
      integer, parameter :: complex_t = 0
      integer(KIND=8), allocatable, dimension(:) :: element
      integer(kind=8), external :: sfCreate_solve
      integer(kind=8), external :: sfGetElement
      integer, external :: sfFactor
      character(len=24) :: fileName ! chem_in_xxxxxx.bin
      integer :: iunit, sc_t_size, sc_t_DynSiz

      !---------------------------------------------------------
      !Introdução da escrita dos parâmetros de entrada e saída
      !---------------------------------------------------------
      ! criar o arquivo de dados
      write (filename, fmt='("chem_inp_",I6.6,"-",I4.4,".bin")') int(time),mynum
      open (newunit=iunit, file=filename, status="replace", form="unformatted", action="write")

      write(iunit)  m1
      write(iunit)  m2
      write(iunit)  m3
      write(iunit)  nspecies
      write(iunit)  nr_photo
      write(iunit)  nr
      write(iunit)  maxnspecies
      write(iunit)  nspecies_chem_transported
      write(iunit)  nspecies_chem_no_transported
      write(iunit)  nob
      write(iunit)  maxblock_size
      write(iunit)  na_extra3d
      write(60+mynum,*) mynum," - na_extra3d: ",na_extra3d
      do ispc=1,nspecies
        sc_t_size   = size(chem1_g(ispc)%sc_t)
        sc_t_DynSiz = size(chem1_g(ispc)%sc_t_dyn)
        write(iunit)  sc_t_size  
        write(iunit)  sc_t_DynSiz
      end do
      write(iunit)  dtlt
      write(iunit)  pp  !(m1, m2, m3)
      write(iunit)  pi0 !(m1, m2, m3)
      write(iunit)  theta !(m1, m2, m3)
      write(iunit)  rv !(m1, m2, m3)
      write(iunit)  dn0 !(m1, m2, m3)
      write(iunit)  cosz !(m2, m3)
      write(iunit)  rcp !(m1, m2, m3)
      write(iunit)  weight !(nspecies)
      write(iunit)  PhotojMethod
      write(iunit)  transp_chem_index !(maxnspecies)
      write(iunit)  no_transp_chem_index !(maxnspecies)
      write(iunit)  split_method
      write(iunit)  N_DYN_CHEM
      write(iunit)  block_end!(:)nob
      write(iunit)  indexk!(:, :)(maxblock_size,nob)
      write(iunit)  indexi!(:, :)(maxblock_size,nob)
      write(iunit)  indexj!(:, :)(maxblock_size,nob)
      write(iunit)  kij_index!(:, :)(maxblock_size,nob)
      write(iunit)  atol !(nspecies)
      write(iunit)  rtol !(nspecies)
      write(iunit)  cp
      write(iunit)  cpor
      write(iunit)  p00
      write(iunit)  jphoto !(nr_photo, m1, m2, m3)
      write(iunit)  att !(m2, m3)
      !!!!
      print *,mynum,nspecies,sc_t_size,sc_t_DynSiz
      do ispc=1,nspecies
        write(iunit)  chem1_g(ispc)%sc_p
        !write(iunit,*)  chem1_g(ispc)%sc_pp
        !write(iunit,*)  chem1_g(ispc)%sc_pf
        !write(iunit,*)  chem1_g(ispc)%sc_dd
        !write(iunit,*)  chem1_g(ispc)%sc_wd
        write(iunit)  chem1_g(ispc)%sc_t
        write(iunit)  chem1_g(ispc)%sc_t_dyn
        !print *,mynum,ispc,size(chem1_g(ispc)%sc_t),size(chem1_g(ispc)%sc_t_dyn)
      end do
      write(iunit)  spack(1)%DLdrdc
      write(iunit)  spack(1)%sc_p_new
      !write(iunit,*)  spack(1)%sc_p_4
      write(iunit)  spack(1)%DLr	 
      !write(iunit,*)  spack(1)%DLr3	 
      write(iunit)  spack(1)%jphoto 
      write(iunit)  spack(1)%rk    
      write(iunit)  spack(1)%w    
      write(iunit)  spack(1)%sc_p
      write(iunit)  spack(1)%temp	
      write(iunit)  spack(1)%press 
      write(iunit)  spack(1)%cosz	  
      write(iunit)  spack(1)%att	  
      write(iunit)  spack(1)%vapp 
      write(iunit)  spack(1)%volmol 
      write(iunit)  spack(1)%volmol_i 
      write(iunit)  spack(1)%xlw 
      write(iunit)  spack(1)%err    
      write(60+mynum,*) mynum," - spack_2d: size",size(spack_2d,1),size(spack_2d,2)
      do i = 1, nob
        do ijk = 1, block_end(i)
            write(60+mynum,*) "i, ijk, size(DLmat): ",i,ijk,size(spack_2d(ijk,inob)%DLmat,1),size(spack_2d(ijk,inob)%DLmat,2)
            write(iunit)  spack_2d(ijk,inob)%DLmat 
            write(iunit)  spack_2d(ijk,inob)%DLb1	
            write(iunit)  spack_2d(ijk,inob)%DLb2	
            write(iunit)  spack_2d(ijk,inob)%DLb3	
            write(iunit)  spack_2d(ijk,inob)%DLb4	
            write(iunit)  spack_2d(ijk,inob)%DLk1  
            write(iunit)  spack_2d(ijk,inob)%DLk2  
            write(iunit)  spack_2d(ijk,inob)%DLk3  
            write(iunit)  spack_2d(ijk,inob)%DLk4 
        end do
      end do 
      write(iunit)  last_accepted_dt !(:) !IO !(nob)
      write(iunit)  get_non_zeros !IO
      do i=1,na_extra3d
        write(iunit)  extra3d(i)%d3 !(na_extra3d) ! Output JNO2 !IO
      end do

      close(iunit)

!!  if(.not. spack_alloc) then
!!  !- Determining number of blocks and masks, allocating spack type
!!   call allocIndex(block_size,N_DYN_CHEM)
!!   call alloc_spack(CHEMISTRY)
!!  end if

!  IF(.NOT. get_non_zeros) then
      !- Get always zero elements of Jacobian matrix to save
      !- computational time at solver:
      ! non-opt
      !get_non_zeros = .true.
      !call Prepare(nspecies)

      !  opt
   !!call get_always_zeros(nr_photo,nr,nspecies,get_non_zeros,spack(1)%rk(1,1:nr), &
      !                          spack(1)%jphoto(1,1:nr_photo),p00)

      call get_number_nonzeros(nr_photo, nr, nspecies, spack(1)%rk(1, 1:nr), &
                               spack(1)%jphoto(1, 1:nr_photo), p00, maxnonzeros)

      sizeOfMatrix = nspecies

      allocate (ipos(maxnonzeros)); ipos = 1
      allocate (jpos(maxnonzeros)); jpos = 1
!  ENDIF

      do i = 1, nob !- loop over all blocks

         !- copying structure from input to internal

         do ijk = 1, block_end(i) !index_g%block_end(i)

            kij_ = kij_index(ijk, i) !index_g%kij_index(ijk,i)!- kij brams tendency index
            k_ = indexk(ijk, i) ! index_g%indexk(ijk,i) ! k brams index
            i_ = indexi(ijk, i) ! index_g%indexi(ijk,i) ! i brams index
            j_ = indexj(ijk, i) ! index_g%indexj(ijk,i) ! j brams index

            !- Exner function divided by Cp
            spack(inob)%press(ijk) = (pp(k_, i_, j_) + pi0(k_, i_, j_))/cp

            !- Air temperature (K)
            spack(inob)%temp(ijk) = theta(k_, i_, j_)*spack(inob)%press(ijk)

            !- transform from Exner function to pressure
            spack(inob)%press(ijk) = spack(inob)%press(ijk)**cpor*p00

            !- Water vapor
            spack(inob)%vapp(ijk) =  rv(k_,i_,j_)!&
            !                  * dn0(k_,i_,j_)*6.02D20/18.D0

            !- volmol e' usado para converter de ppbm para molecule/cm^3 (fator fxc incluido)
            ! spack(inob)%volmol(ijk)=(spack(inob)%press(ijk))/(8.314e0*spack(inob)%temp(ijk))
            ! fxc = (6.02e23*1e-15*spack(inob)%volmol(ijk))*pmar
            !
            spack(inob)%volmol(ijk) = (6.02d23*1d-15*pmar)*(spack(inob)%press(ijk))/(8.314d0*spack(inob)%temp(ijk))

            !- inverse of volmol to get back to ppbm at the end of routine
            spack(inob)%volmol_I(ijk) = 1.0d0/spack(inob)%volmol(ijk)

            !- get liquid water content
            spack(inob)%xlw(ijk) = rcp(k_, i_, j_)*dn0(k_, i_, j_)*1.d-3

            !- convert from brams chem (ppbm) arrays to spack (molec/cm3)
            !- no transported species section
            do ispc = 1, nspecies_chem_no_transported

               !- map the species to NO transported ones
               n = no_transp_chem_index(ispc)

               !- initialize no-transported species (don't need to convert, because these
               !- species are already saved using molecule/cm^3 units)
               spack(inob)%sc_p(ijk, n) = chem1_g(n)%sc_p(k_, i_, j_)
!       spack(inob)%sc_p_4(ijk,n) = chem1_g(n)%sc_p(k_,i_,j_)
            end do

         end do

         !- convert from brams chem (ppbm) arrays to spack (molec/cm3)
         !- transported species section
         if (split_method .eq. 'PARALLEL' .and. N_DYN_CHEM .gt. 1) then
            do ijk = 1, block_end(i) !index_g%block_end(i)

               kij_ = kij_index(ijk, i) !index_g%kij_index(ijk,i)!- kij brams tendency index
               k_ = indexk(ijk, i) !index_g%indexk(ijk,i) ! k brams index
               i_ = indexi(ijk, i) !index_g%indexi(ijk,i) ! i brams index
               j_ = indexj(ijk, i) !index_g%indexj(ijk,i) ! j brams index

               do ispc = 1, nspecies_chem_transported

                  !- map the species to transported ones
                  n = transp_chem_index(ispc)

                  !- conversion from ppbm to molecule/cm^3
                  !- get back the concentration at the begin of the chemistry timestep integration:

                  spack(inob)%sc_p(ijk, n) = (chem1_g(n)%sc_p(k_, i_, j_) - &  ! updated mixing ratio
                                              chem1_g(n)%sc_t_dyn(kij_)*N_DYN_CHEM*dtlt) &  ! accumulated tendency
                                             *spack(inob)%volmol(ijk)/weight(n)
                  !- for testing
                  !spack(inob)%sc_p(ijk,n)   = max(0.,spack(inob)%sc_p(ijk,n))
                  !spack(inob)%sc_p_4(ijk,n) =(chem1_g(n)%sc_p(k_,i_,j_) - beta*chem1_g(n)%sc_t_dyn(kij_)*N_DYN_CHEM*dtlt)& !dtlt \E9 o dinamico
                  !                                * spack(inob)%volmol(ijk)/weight(n)
               end do
            end do
         else

            do ijk = 1, block_end(i) !index_g%block_end(i)

               kij_ = kij_index(ijk, i) !index_g%kij_index(ijk,i)!- kij brams tendency index
               k_ = indexk(ijk, i) !index_g%indexk(ijk,i) ! k brams index
               i_ = indexi(ijk, i) !index_g%indexi(ijk,i) ! i brams index
               j_ = indexj(ijk, i) !index_g%indexj(ijk,i) ! j brams index

               do ispc = 1, nspecies_chem_transported

                  !- map the species to transported ones
                  n = transp_chem_index(ispc)

                  !- conversion from ppbm to molecule/cm^3
                  spack(inob)%sc_p(ijk, n) = chem1_g(n)%sc_p(k_, i_, j_)*spack(inob)%volmol(ijk)/weight(n)
                  spack(inob)%sc_p(ijk, n) = max(0.d0, spack(inob)%sc_p(ijk, n))

               end do
            end do
         end if

         !- Photolysis section
         if (trim(PhotojMethod) .eq. 'FAST-JX' .or. trim(PhotojMethod) .eq. 'FAST-TUV') then

            do ijk = 1, block_end(i) !index_g%block_end(i)

               k_ = indexk(ijk, i) !index_g%indexk(ijk,i) ! k brams index
               i_ = indexi(ijk, i) !index_g%indexi(ijk,i) ! i brams index
               j_ = indexj(ijk, i) !index_g%indexj(ijk,i) ! j brams index

               do n = 1, nr_photo
                  spack(inob)%jphoto(ijk, n) = jphoto(n, k_, i_, j_) !fast_JX_g%jphoto(n,k_,i_,j_)
               end do

            end do

         elseif (trim(PhotojMethod) .eq. 'LUT') then

            do ijk = 1, block_end(i) !index_g%block_end(i)

               i_ = indexi(ijk, i) !index_g%indexi(ijk,i) ! i brams index
               j_ = indexj(ijk, i) !index_g%indexj(ijk,i) ! j brams index

               !- UV  attenuation un function AOT
               spack(inob)%att(ijk) = att(i_, j_) !uv_atten_g%att(i_,j_)

               !- get zenital angle (for LUT PhotojMethod)
               spack(inob)%cosz(ijk) = cosz(i_, j_)

            end do
         end if

!- compute kinetical and photochemical reactions (array: spack(inob)%rk)
         call kinetic(nr_photo, spack(inob)%Jphoto &
                      , spack(inob)%rk &
                      , spack(inob)%temp &
                      , spack(inob)%vapp &
                      , spack(inob)%Press &
                      , spack(inob)%cosz &
                      , spack(inob)%att &
                      , 1 &
                      , block_end(i) & !index_g%block_end(i), &
                      , maxblock_size, nr)

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>     !tmp
!
!  !- srf: to use these arrays, be sure NA_EXTRAD3D RAMSIN is at least 8
!    DO ijk=1,block_end(i) !index_g%block_end(i)
!       k_=indexk(ijk,i) !index_g%indexk(ijk,i) ! k brams index
!       i_=indexi(ijk,i) !index_g%indexi(ijk,i) ! i brams index
         !      j_=indexj(ijk,i) !index_g%indexj(ijk,i) ! j brams index

         !  extra3d(7)%d3(k_,i_,j_) =spack(inob)%rk(ijk,1)
         !extra3d(8)%d3(k_,i_,j_) =spack(inob)%rk(ijk,2)
         !extra3d(9)%d3(k_,i_,j_) =spack(inob)%rk(ijk,3)
         !extra3d(10)%d3(k_,i_,j_)=spack(inob)%rk(ijk,4)
         !extra3d(11)%d3(k_,i_,j_)=spack(inob)%rk(ijk,5)
         !extra3d(12)%d3(k_,i_,j_)=spack(inob)%rk(ijk,6)
         !extra3d(13)%d3(k_,i_,j_)=spack(inob)%rk(ijk,7)
         !extra3d(14)%d3(k_,i_,j_)=spack(inob)%rk(ijk,8)
         !extra3d(15)%d3(k_,i_,j_)=spack(inob)%rk(ijk,9)
         !extra3d(16)%d3(k_,i_,j_)=spack(inob)%rk(ijk,10)
         !extra3d(17)%d3(k_,i_,j_)=spack(inob)%rk(ijk,11)
         !extra3d(18)%d3(k_,i_,j_)=spack(inob)%rk(ijk,12)
         !extra3d(19)%d3(k_,i_,j_)=spack(inob)%rk(ijk,13)
         !extra3d(20)%d3(k_,i_,j_)=spack(inob)%rk(ijk,14)
         !extra3d(21)%d3(k_,i_,j_)=spack(inob)%rk(ijk,16)
         !extra3d(22)%d3(k_,i_,j_)=spack(inob)%rk(ijk,17)
         !extra3d(23)%d3(k_,i_,j_)=spack(inob)%rk(ijk,18)
         !extra3d(24)%d3(k_,i_,j_)=spack(inob)%rk(ijk,19)
         !extra3d(25)%d3(k_,i_,j_)=spack(inob)%rk(ijk,20)
         !extra3d(26)%d3(k_,i_,j_)=spack(inob)%rk(ijk,21)
         ! extra3d(27)%d3(k_,i_,j_)=spack(inob)%rk(ijk,22)
         !  extra3d(28)%d3(k_,i_,j_)=spack(inob)%rk(ijk,23)

         ! ENDDO
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>     !tmp

!- ROSENBROCK METHOD ----------------------------------------------------------------------
!- kml: Inicio do equivalente a roschem
!- srf: extending to RODAS 3 (4 stages, order 3)

         dt_chem = last_accepted_dt(i) !index_g%last_accepted_dt(i) ! dble(dtlt)
         dt_min = max(1.0d0, 1.d-2*dble(dtlt*N_DYN_CHEM))
         dt_max = dble(dtlt*N_DYN_CHEM)
         dt_new = 0.0d0

         time_c = 0.0d0
         time_f = dble(dtlt*N_DYN_CHEM)

         run_until_integr_ends: do while (time_c + Roundoff .lt. time_f)

!-   Compute the Jacobian (DLRDC).
            call jacdchemdc(spack(inob)%sc_p &
                            !     CALL jacdchemdc (spack(inob)%sc_p_4  &
                            , spack(inob)%rk &
                            , spack(inob)%DLdrdc & ! Jacobian matrix
                            , nspecies, 1 &
                            , block_end(i) & !index_g%block_end(i), &
                            , maxblock_size, nr)

!-   Compute chemical net production terms (actually, P-L) at initial time (array spack(inob)%DLr)
            call fexchem(spack(inob)%sc_p &
                         !     CALL fexchem (spack(inob)%sc_p_4 &
                         , spack(inob)%rk &
                         , spack(inob)%DLr & !production term
                         , nspecies, 1 &
                         , block_end(i) & !index_g%block_end(i), &
                         , maxblock_size, nr)

            do Ji = 1, nspecies
               do ijk = 1, block_end(i) !index_g%block_end(i)
                  spack_2d(ijk, inob)%DLb1(Ji) = spack(inob)%DLr(ijk, Ji)
               end do
            end do

!-   compute matrix (1/(Igamma*dt)  - Jacobian)  where Jacobian = DLdrdc
!-   fill DLMAT with non-diagonal (fixed in the timestep) Jacobian
            do Jj = 1, nspecies
               do Ji = 1, nspecies
                  do ijk = 1, block_end(i) !index_g%block_end(i)

                     spack_2d(ijk, inob)%DLmat(Ji, Jj) = -spack(inob)%DLdrdc(ijk, Ji, Jj)

                  end do
               end do
            end do

            UntilAccepted: do

!- fill DLMAT with diagonal Jacobian (which changes in time, because of dt_chem)
               Igamma_dtstep = 1.0d0/(Igamma*dt_chem)

               do Jj = 1, nspecies
                  do ijk = 1, block_end(i) !index_g%block_end(i)
                     spack_2d(ijk, inob)%DLmat(Jj, Jj) = Igamma_dtstep - spack(inob)%DLdrdc(ijk, Jj, Jj)
                  end do
               end do

               if ((i .eq. 1)) then
                  blocknonzeros = 0
                  do Ji = 1, nspecies
                     do Jj = 1, nspecies
                        !IF(spack_2d(1,inob)%DLmat(Ji,Jj)<=0 .AND. spack_2d(1,inob)%DLmat(Ji,Jj)>=(-1)*0) CYCLE
                        if ((spack_2d(1, inob)%DLmat(Ji, Jj) .ne. 0.0d+00)) then
                           blocknonzeros = blocknonzeros + 1
                           ipos(blocknonzeros) = Ji
                           jpos(blocknonzeros) = Jj
                        end if
                     end do
                  end do
                  NumberOfNonZeros = blocknonzeros

                  !create matrix
                  matrix_Id = sfCreate_solve(sizeOfMatrix, complex_t, error)

                  if (allocated(element)) deallocate (element)
                  allocate (element(NumberOfNonZeros))

                  do nz = 1, NumberOfNonZeros
                     Ji = ipos(nz)
                     Jj = jpos(nz)
                     element(nz) = sfGetElement(matrix_Id, Ji, Jj)
                     !element(nz)=sfGetElement(matrix_Id,INT(DLmat(1,nz,1)),INT(DLmat(1,nz,2)))
                  end do

                  call sfZero(matrix_Id)
                  do nz = 1, NumberOfNonZeros
                     Ji = ipos(nz)
                     Jj = jpos(nz)
                     call sfAdd1Real(element(nz), spack_2d(1, inob)%DLmat(Ji, Jj))
                     !CALL sfAdd1Real(element(nz),DLmat(1,nz,3))
                  end do

                  error = sfFactor(matrix_Id)
               end if

               !@LNCC: begin points loop
               do ijk = 1, block_end(i) !index_g%block_end(i)
                  call sfZero(matrix_Id)
                  do nz = 1, NumberOfNonZeros
                     Ji = ipos(nz)
                     Jj = jpos(nz)
                     call sfAdd1Real(element(nz), spack_2d(ijk, inob)%DLmat(Ji, Jj))
                     !CALL sfAdd1Real(element(nz),DLmat(ijk,nz,3))
                  end do
                  error = sfFactor(matrix_Id)

!---------------------------------------------------------------------------------------------------------------------
!    1- First step
!-   Compute DLk1 by Solving (1/(Igamma*dt) - DLRDC) DLk1=DLR

                  !- solver sparse
                  call sfSolve(matrix_Id, spack_2d(ijk, inob)%DLb1, spack_2d(ijk, inob)%DLk1)
!
!---------------------------------------------------------------------------------------------------------------------
!    2- Second step
!    compute   K2 by solving (1/0.5 h JAC)K2 = 4/h * K1 +  F(Yn)
!    compute DLK2 by solving (1/Igama*dt - DLRDC)DLK2 = 4/h * DLK1 +  DLb1
                  dt_chem_i = 1.0d0/dt_chem
                  do Ji = 1, nspecies
                     spack_2d(ijk, inob)%DLb2(Ji) = (4.0d0*dt_chem_i)*spack_2d(ijk, inob)%DLk1(Ji) + &
                                                    spack_2d(ijk, inob)%DLb1(Ji)
                  end do
                  call sfSolve(matrix_Id, spack_2d(ijk, inob)%DLb2, spack_2d(ijk, inob)%DLk2)

!---------------------------------------------------------------------------------------------------------------------
!    3- Third step
!    a) update concentrations

                  !dt_chem_i = 1.0d0/dt_chem
                  do Ji = 1, nspecies
                     spack(inob)%sc_p_new(ijk, Ji) = spack(inob)%sc_p(ijk, Ji) + 2.0d0*spack_2d(ijk, inob)%DLk1(Ji)

                     if (spack(inob)%sc_p_new(ijk, Ji) .lt. Threshold1) then
                        spack(inob)%sc_p_new(ijk, Ji) = Threshold1
                        spack_2d(ijk, inob)%DLk1(Ji) = 0.5d0*(spack(inob)%sc_p_new(ijk, Ji) - spack(inob)%sc_p(ijk, Ji))
                     end if
                  end do
!
!    b) update the net production term (DLr= P-L = F(Y3)) at this stage with the first-order
!       approximation with the new concentration
!
                  call fexchem(spack(inob)%sc_p_new &
                               , spack(inob)%rk &
                               , spack(inob)%DLr &
                               , nspecies, ijk &
                               , ijk & !index_g%block_end(i), &
                               , maxblock_size, nr)

!    c) compute   K3 by solving (1 /(0.5 h) - JAC)K3 =   F(Y3) + 0.5 (K1-K2)
                  do Ji = 1, nspecies
                     spack_2d(ijk, inob)%DLb3(Ji) = spack(inob)%DLr(ijk, Ji) + dt_chem_i* &
                                                    (spack_2d(ijk, inob)%DLk1(Ji) - spack_2d(ijk, inob)%DLk2(Ji))
                  end do
                  call sfSolve(matrix_Id, spack_2d(ijk, inob)%DLb3, spack_2d(ijk, inob)%DLk3)

!---------------------------------------------------------------------------------------------------------------------
!    4- Fourth step
!    a) update concentrations
!       Y4 = Yn + 2 * k1 +  K3
                  dt_chem_i = 1.0d0/dt_chem
                  do Ji = 1, nspecies
                     spack(inob)%sc_p_new(ijk, Ji) = spack(inob)%sc_p(ijk, Ji) + &
                                                     2.0d0*spack_2d(ijk, inob)%DLk1(Ji) + &
                                                     spack_2d(ijk, inob)%DLk3(Ji)

                     if (spack(inob)%sc_p_new(ijk, Ji) .lt. Threshold1) then
                        spack(inob)%sc_p_new(ijk, Ji) = Threshold1
                        spack_2d(ijk, inob)%DLk3(Ji) = (spack(inob)%sc_p_new(ijk, Ji) - spack(inob)%sc_p(ijk, Ji)) &
                                                       - 2.0d0*spack_2d(ijk, inob)%DLk1(Ji)
                     end if
                  end do
!    b) update the net production term (DLr= P-L = F(Y4) ) at this stage with the 3rd-order
!       approximation with the new concentration
!
                  call fexchem(spack(inob)%sc_p_new & ! Y4
                               , spack(inob)%rk &
                               , spack(inob)%DLr & ! F(Y4)
                               , nspecies, ijk &
                               , ijk & !index_g%block_end(i)
                               , maxblock_size, nr)

!    c) compute   K4 by solving (1/(0.5 h)- JAC)K4 =    F(Y3) + K1/h -K2/h -8/3 K3/h
                  do Ji = 1, nspecies
                     spack_2d(ijk, inob)%DLb4(Ji) = spack(inob)%DLr(ijk, Ji) + dt_chem_i* & ! F(Y4)
                                                    (spack_2d(ijk, inob)%DLk1(Ji) &
                                                     - spack_2d(ijk, inob)%DLk2(Ji) &
                                                     - c83*spack_2d(ijk, inob)%DLk3(Ji))
                  end do
                  call sfSolve(matrix_Id, spack_2d(ijk, inob)%DLb4, spack_2d(ijk, inob)%DLk4)

               end do
               !@LNCC: end points loop

!---------------------------------------------------------------------------------------------------------------------
!   - the solution
               dt_chem_i = 1.0d0/dt_chem
               do Ji = 1, nspecies
                  do ijk = 1, block_end(i) !index_g%block_end(i)

                     spack(inob)%sc_p_new(ijk, Ji) = spack(inob)%sc_p(ijk, Ji) + &
                                                     2.0d0*spack_2d(ijk, inob)%DLk1(Ji) &
                                                     + spack_2d(ijk, inob)%DLk3(Ji) &
                                                     + spack_2d(ijk, inob)%DLk4(Ji)
                     spack(inob)%sc_p_new(ijk, Ji) = max(spack(inob)%sc_p_new(ijk, Ji), Threshold1)

                     ! IF (spack(inob)%sc_p_new(ijk,Ji) .LT. Threshold1) THEN
                     !    spack(inob)%sc_p_new(ijk,Ji) = Threshold1
                     !    spack_2d(ijk,inob)%DLk3(Ji) = (spack(inob)%sc_p_new(ijk,Ji) - spack(inob)%sc_p(ijk,Ji)) * dt_chem_i
                     ! ENDIF
                     !- keep non-negative values for concentrations
                     !   spack(inob)%sc_p_new(ijk,Ji) = max(spack(inob)%sc_p_new(ijk,Ji),Threshold1)
                     !- update DLK2, in case sc_p_new is changed by the statement above
                     !   spack_2d(ijk,inob)%DLk2(Ji) = 2.0D0*( (spack(inob)%sc_p_new(ijk,Ji)-spack(inob)%sc_p(ijk,Ji)) &
                     !                                         * dt_chem_i - 1.5D0 * spack_2d(ijk,inob)%DLk1(Ji) )
                  end do
               end do

!
!-  Compute the error estimation : spack(inob)%err
               do ijk = 1, block_end(i) !index_g%block_end(i)

                  spack(inob)%err(ijk) = 0.0d0
                  do Ji = 1, nspecies

                     Tol = ATol(Ji) + RTol(Ji)*DMAX1(DABS(spack(inob)%sc_p(ijk, Ji)), DABS(spack(inob)%sc_p_new(ijk, Ji)))

                     err1 = spack_2d(ijk, inob)%DLk4(Ji)

                     spack(inob)%err(ijk) = spack(inob)%err(ijk) + (err1/Tol)**2.0d0

                  end do

                  spack(inob)%err(ijk) = DMAX1(uround, DSQRT(spack(inob)%err(ijk)/nspecies))
               end do

               all_accepted = 1
               do ijk = 1, block_end(i) !index_g%block_end(i)
                  if (spack(inob)%err(ijk) - Roundoff .gt. 1.0d0) then
                     all_accepted = 0; exit
                  end if
               end do

               !- find the maximum error occurred
               max_err1 = maxval(spack(inob)%err(1:block_end(i))) !index_g%block_end(i) ) )

               !- use it to determine the new time step for all block elements
               !- new step size is bounded by FacMin <= Hnew/H <= FacMax
               Fac = min(FacMax, max(FacMin, FacSafe/max_err1**(1.0d0/elo)))

               !- possible new timestep
               dt_new = dt_chem*Fac
               dt_new = max(dt_min, min(dt_max, dt_new))

               !- to reset the timestep resizing in function of the error estimation, use the statements below:
               ! all_accepted = 1; dt_new=dt_max

               if (all_accepted .eq. 0 .and. dt_new .gt. dt_min) then  ! current solution is not accepted

                  !- resize the timestep and try again
                  dt_chem = dt_new

               else    ! current solution is     accepted

                  !- go ahead, updating spack(inob)%sc_p with the solution (spack(inob)%sc_p_new)
                  !- next time
                  time_c = time_c + dt_chem

                  !- next timestep (dt_new but limited by the time_f-time_c, the rest of time integration interval)
                  dt_chem = min(dt_new, time_f - time_c)

                  !- save the accepted timestep for the next integration interval
                  if (time_c .lt. time_f) last_accepted_dt(i) = dt_new !index_g%last_accepted_dt(i) =    dt_new

                  !- pointer (does not work yet)
                  ! spack(inob)%sc_p=>spack(inob)%sc_p_new     ! POINTER
                  !- copy
                  do Ji = 1, nspecies
                     do ijk = 1, block_end(i) !index_g%block_end(i)
                        spack(inob)%sc_p(ijk, Ji) = spack(inob)%sc_p_new(ijk, Ji)
                        ! spack(inob)%sc_p_4(ijk,Ji) = spack(inob)%sc_p_new(ijk,Ji)
                     end do
                  end do

                  exit UntilAccepted

               end if

            end do UntilAccepted

         end do run_until_integr_ends ! time-spliting

         !--------------------------------------------------------------------------------------------
         !- Restoring species tendencies OR updated mixing ratios from internal to brams structure

         !- transported species section

         if (split_method .eq. 'PARALLEL' .and. N_DYN_CHEM .gt. 1) then

            do ijk = 1, block_end(i) !index_g%block_end(i)

               kij_ = kij_index(ijk, i) !index_g%kij_index(ijk,i)
               k_ = indexk(ijk, i) !index_g%indexk(ijk,i)
               i_ = indexi(ijk, i) !index_g%indexi(ijk,i)
               j_ = indexj(ijk, i) !index_g%indexj(ijk,i)

               do ispc = 1, nspecies_chem_transported

                  !- map the species to transported ones
                  n = transp_chem_index(ispc)

                  !- include the chemical tendency at total tendency (convert to unit: ppbm/s)
                  chem1_g(n)%sc_p(k_, i_, j_) = chem1_g(n)%sc_t_dyn(kij_)*N_DYN_CHEM*dtlt + &
                                                spack(inob)%sc_p(ijk, n)*weight(n)*spack(inob)%volmol_I(ijk)

                  chem1_g(n)%sc_p(k_, i_, j_) = max(0., chem1_g(n)%sc_p(k_, i_, j_))

               end do
            end do

         elseif (split_method .eq. 'PARALLEL' .and. N_DYN_CHEM .eq. 1) then

            dble_dtlt_i = 1.0d0/dble(dtlt)
            do ijk = 1, block_end(i) !index_g%block_end(i)

               kij_ = kij_index(ijk, i) !index_g%kij_index(ijk,i)
               k_ = indexk(ijk, i) !index_g%indexk(ijk,i)
               i_ = indexi(ijk, i) !index_g%indexi(ijk,i)
               j_ = indexj(ijk, i) !index_g%indexj(ijk,i)

               do ispc = 1, nspecies_chem_transported

                  !- map the species to transported ones
                  n = transp_chem_index(ispc)

                  !- include the chemical tendency at total tendency (convert to unit: ppbm/s)
!           chem1_g(n)%sc_t(kij_) =                    +  &! use this for update only chemistry (No dyn/emissions)
                  chem1_g(n)%sc_t(kij_) = chem1_g(n)%sc_t(kij_) + &! previous tendency
                                          (spack(inob)%sc_p(ijk, n)*weight(n)*spack(inob)%volmol_I(ijk) - &! new mixing ratio
                                           chem1_g(n)%sc_p(k_, i_, j_)) &! old mixing ratio
                                          *dble_dtlt_i                                                        ! inverse of timestep
               end do
            end do

         else

            do ijk = 1, block_end(i) !index_g%block_end(i)

               kij_ = kij_index(ijk, i) !index_g%kij_index(ijk,i)
               k_ = indexk(ijk, i) !index_g%indexk(ijk,i)
               i_ = indexi(ijk, i) !index_g%indexi(ijk,i)
               j_ = indexj(ijk, i) !index_g%indexj(ijk,i)

               do ispc = 1, nspecies_chem_transported

                  !- map the species to transported ones
                  n = transp_chem_index(ispc)

                  chem1_g(n)%sc_p(k_, i_, j_) = spack(inob)%sc_p(ijk, n)*weight(n)*spack(inob)%volmol_I(ijk)
                  chem1_g(n)%sc_p(k_, i_, j_) = max(0., chem1_g(n)%sc_p(k_, i_, j_))
               end do
            end do

         end if

         !- no transported species section
         do ijk = 1, block_end(i) !index_g%block_end(i)

            kij_ = kij_index(ijk, i) !index_g%kij_index(ijk,i)
            k_ = indexk(ijk, i) !index_g%indexk(ijk,i)
            i_ = indexi(ijk, i) !index_g%indexi(ijk,i)
            j_ = indexj(ijk, i) !index_g%indexj(ijk,i)
            do ispc = 1, nspecies_chem_no_transported

               !- map the species to no transported ones
               n = no_transp_chem_index(ispc)

               !- save no-transported species (keep current unit : molec/cm3)
               chem1_g(n)%sc_p(k_, i_, j_) = max(0., real(spack(inob)%sc_p(ijk, n)))
            end do

         end do

      end do ! enddo loop over all blocks


      write (filename, fmt='("chem_out_",I6.6,"-",I4.4,".bin")') int(time),mynum
      open (newunit=iunit, file=filename, status="replace", form="unformatted", action="write")
      do ispc=1,nspecies
        write(iunit)  chem1_g(ispc)%sc_p
        !write(iunit,*)  chem1_g(ispc)%sc_pp
        !write(iunit,*)  chem1_g(ispc)%sc_pf
        !write(iunit,*)  chem1_g(ispc)%sc_dd
        !write(iunit,*)  chem1_g(ispc)%sc_wd
        write(iunit)  chem1_g(ispc)%sc_t
        write(iunit)  chem1_g(ispc)%sc_t_dyn
      end do
      write(iunit)  spack(1)%DLdrdc
      write(iunit)  spack(1)%sc_p_new
      !write(iunit,*)  spack(1)%sc_p_4
      write(iunit)  spack(1)%DLr	 
      !write(iunit,*)  spack(1)%DLr3	  
      write(iunit)  spack(1)%jphoto 
      write(iunit)  spack(1)%rk    
      write(iunit)  spack(1)%w    
      write(iunit)  spack(1)%sc_p
      write(iunit)  spack(1)%temp	
      write(iunit)  spack(1)%press 
      write(iunit)  spack(1)%cosz	  
      write(iunit)  spack(1)%att	  
      write(iunit)  spack(1)%vapp 
      write(iunit)  spack(1)%volmol 
      write(iunit)  spack(1)%volmol_i 
      write(iunit)  spack(1)%xlw 
      write(iunit)  spack(1)%err    
      do i = 1, nob
        do ijk = 1, block_end(i)
            write(iunit)  spack_2d(ijk,inob)%DLmat 
            write(iunit)  spack_2d(ijk,inob)%DLb1	
            write(iunit)  spack_2d(ijk,inob)%DLb2	
            write(iunit)  spack_2d(ijk,inob)%DLb3	
            write(iunit)  spack_2d(ijk,inob)%DLb4	
            write(iunit)  spack_2d(ijk,inob)%DLk1  
            write(iunit)  spack_2d(ijk,inob)%DLk2  
            write(iunit)  spack_2d(ijk,inob)%DLk3  
            write(iunit)  spack_2d(ijk,inob)%DLk4 
        end do
      end do 
      write(iunit)  last_accepted_dt !(:) !IO !(nob)
      write(iunit)  get_non_zeros !IO
      do i=1,na_extra3d
        write(iunit)  extra3d(i)%d3 !(na_extra3d) ! Output JNO2 !IO
      end do

      close(iunit)



!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!DO n=1,nspecies
!print*,'max1=',spc_name(n),maxval(chem1_g(n)%sc_p(:,:,:)),maxloc(chem1_g(n)%sc_p(:,:,:))
!print*,'max2=',n,spc_name(n),maxval(spack(:)%sc_p(:,n)),maxloc(spack(:)%sc_p(:,n))
!call flush(6)
!END DO
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

   end subroutine chem_rodas3_dyndt

!---------------------------------------------------------------------------------------------------
end module mod_chem_spack_rodas3_dyndt
