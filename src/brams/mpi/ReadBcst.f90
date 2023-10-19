module ReadBcst

  use mem_grid, only: &
       ngrids, nnxp, nnyp, nnzp, nzs, nzg, npatch, &
       time, iyear1, imonth1, idate1, itime1, &
       GlobalSizes

  use node_mod, only:      &
       mchnum, mynum, master_num, &
       nmachs, nodemxp, nodemyp, &
       nodeia, nodeiz, nodeja, nodejz, nodeibcon, &
       nodei0, nodej0

  use an_header, only: &
       IOFileDS,       &
       FieldWriteStoreInfo

  use mem_aerad, only: &
       nwave

  use mem_basic, only:  &
       basic_g

  use mem_turb,  only:   &
       turb_g, idiffk, xkhkm

  use var_tables, only: &
       var_tables_r, &
       num_var,      &
       vtab_r

  use ParLib, only: &
       parf_bcast, &
       parf_minloc, &
       parf_allreduce_max, &
       parf_GatherAllChunks


  implicit none

  private
  public :: ReadStoreOwnChunk
  public :: ReadStoreFullFieldAndOwnChunk
  public :: Broadcast
  public :: ProcWithMin
  public :: MaxCFLOverall
  public :: LocalSizesAndDisp
  public :: PreProcAndGather
  public :: RearrangeAndDump
  public :: DumpFullField
  public :: DumpVTabEntry
  public :: DumpVarTables
  public :: gatherData

!--(DMK-CCATT-INI)-----------------------------------------------------
  public :: storeOwnChunk_3D
!--(DMK-CCATT-FIM)-----------------------------------------------------


  interface ReadStoreOwnChunk
     module procedure ReadStoreOwnChunk_2D, ReadStoreOwnChunk_3D
  end interface

  interface Broadcast
     module procedure Broadcast_I, Broadcast_I1D, &
          Broadcast_R, Broadcast_R1D, Broadcast_R2D, &
          Broadcast_C, Broadcast_C1D
  end interface

  interface gatherData 
     module procedure gatherData2d, gatherData3d, gatherData4d
  end interface

  integer, parameter :: idim_type_min=2
  integer, parameter :: idim_type_max=7
  logical, parameter :: dumpLocal=.false.
  include "i8.h"
  include "UseVfm.h"
contains



  subroutine ReadStoreOwnChunk_2D(grid, fUnit, toStore, fieldName)

    use mem_grid, only:  &
         runtype

    integer, intent(in) :: grid
    integer, intent(in) :: fUnit
    real, pointer :: toStore(:,:)
    character(len=*), intent(in) :: fieldName


    integer :: ldimx, ldimy, lin
    integer :: ia, iz, ja, jz
    integer :: ierr
    character(len=8) :: c0, c1
    character(len=*), parameter :: h="**(ReadStoreOwnChunk_2D)**"

    ! check allocated memory
    if (runtype(1:9)=='MAKEVFILE') then
       ldimx = nnxp(grid)
       ldimy = nnyp(grid)
    else
       ldimx = nodemxp(mynum,grid)
       ldimy = nodemyp(mynum,grid)
    endif
    if (.not. associated(toStore)) then
       call fatal_error(h//" will store at not associated pointer var "//trim(fieldName))
    end if
    lin = size(toStore,1)
    if (ldimx /= lin) then
       write(c0,"(i8)") ldimx
       write(c1,"(i8)") lin
       call fatal_error(h//trim(fieldName)//" x dimension ("//trim(adjustl(c1))//&
            ") differs from required ("//trim(adjustl(c0))//")")
    end if
    lin = size(toStore,2)
    if (ldimy /= lin) then
       write(c0,"(i8)") ldimy
       write(c1,"(i8)") lin
       call fatal_error(h//trim(fieldName)//" y dimension ("//trim(adjustl(c1))//&
            ") differs from required ("//trim(adjustl(c0))//")")
    end if
    ia = nodei0(mynum,grid)+1
    iz = nodei0(mynum,grid)+nodemxp(mynum,grid)
    ja = nodej0(mynum,grid)+1
    jz = nodej0(mynum,grid)+nodemyp(mynum,grid)

    ! allocate scratch area for full domain

    call AllocReadStoreOwnChunk_2D(fUnit, toStore, fieldName, &
         nnxp(grid), nnyp(grid), ldimx, ldimy, ia, iz, ja, jz, &
         mchnum, master_num, runtype)

  end subroutine ReadStoreOwnChunk_2D




  subroutine AllocReadStoreOwnChunk_2D(fUnit, toStore, fieldName, &
       nnxp, nnyp, ldimx, ldimy, ia, iz, ja, jz, &
       mchnum, master_num, runtype)

    integer, intent(in) :: fUnit
    real, pointer :: toStore(:,:)
    character(len=*), intent(in) :: fieldName
    integer, intent(in) :: nnxp
    integer, intent(in) :: nnyp
    integer, intent(in) :: ldimx
    integer, intent(in) :: ldimy
    integer, intent(in) :: ia
    integer, intent(in) :: iz
    integer, intent(in) :: ja
    integer, intent(in) :: jz
    integer, intent(in) :: mchnum
    integer, intent(in) :: master_num
    character(len=*), intent(in) :: runtype


    real :: fullGrid(nnxp,nnyp)
    integer :: lin
    integer :: ierr
    character(len=8) :: c0, c1
    character(len=*), parameter :: h="**(AllocReadStoreOwnChunk_2D)**"
    integer :: i,j,n

    fullGrid=0.0

    !'master process opens file and reads first data into full domain scratch

    if (mchnum == master_num) then
       if (useVfm) then
          call vfirec(fUnit,fullGrid(1,1),nnxp*nnyp,'LIN')
       else
          read(fUnit) fullGrid
       end if
    end if
    
    ! broadcast full domain scratch; 
    !local chunk is extracted and stored at desired variable
    if (runtype(1:9)/='MAKEVFILE') then
       call parf_bcast(fullGrid, int(nnxp,i8), int(nnyp,i8), &
            master_num)
    endif
    
    call mk_2_buff(fullGrid(1,1), toStore(1,1), &
         nnxp, nnyp, ldimx, ldimy, ia, iz, ja, jz)

    if (dumpLocal) then
       write(*,"(a,i4,a,i4,a,i4,a,i4,a,i4,a,i4,a)") &
            h//trim(fieldName)//"(   1:",ldimx, ",   1:",ldimy, &
            ")=fullGrid(",ia,":",iz,",",ja,":",jz,")"
    end if
  end subroutine AllocReadStoreOwnChunk_2D


  subroutine ReadStoreOwnChunk_3D(grid, fUnit, toStore, nz, fieldName)
    use mem_grid, only:  &
         runtype

    integer, intent(in) :: grid
    integer, intent(in) :: fUnit
    integer, intent(in) :: nz
    real, pointer :: toStore(:,:,:)
    character(len=*), intent(in) :: fieldName


    integer :: ldimx, ldimy, lin
    integer :: ia, iz, ja, jz, ka, kz
    integer :: ierr
    character(len=8) :: c0, c1
    character(len=*), parameter :: h="**(ReadStoreOwnChunk_3D)**"

    ! check allocated memory

    if (.not. associated(toStore)) then
       call fatal_error(h//" will store at not associated pointer var "//trim(fieldName))
    end if
    lin = size(toStore,1)
    if (nz /= lin) then
       write(c0,"(i8)") nz
       write(c1,"(i8)") lin
       call fatal_error(h//trim(fieldName)//" z dimension ("//trim(adjustl(c1))//&
            ") differs from required ("//trim(adjustl(c0))//")")
    end if
    ldimx = nodemxp(mynum,grid)
    lin = size(toStore,2)
    if (ldimx /= lin) then
       write(c0,"(i8)") ldimx
       write(c1,"(i8)") lin
       call fatal_error(h//trim(fieldName)//" x dimension ("//trim(adjustl(c1))//&
            ") differs from required ("//trim(adjustl(c0))//")")
    end if
    ldimy = nodemyp(mynum,grid)
    lin = size(toStore,3)
    if (ldimy /= lin) then
       write(c0,"(i8)") ldimy
       write(c1,"(i8)") lin
       call fatal_error(h//trim(fieldName)//" y dimension ("//trim(adjustl(c1))//&
            ") differs from required ("//trim(adjustl(c0))//")")
    end if
    ia = nodei0(mynum,grid)+1
    iz = nodei0(mynum,grid)+nodemxp(mynum,grid)
    ja = nodej0(mynum,grid)+1
    jz = nodej0(mynum,grid)+nodemyp(mynum,grid)

    call AllocReadStoreOwnChunk_3D(fUnit, toStore, fieldName, &
       nz, nnxp(grid), nnyp(grid), ldimx, ldimy, ia, iz, ja, jz, &
       mchnum, master_num, runtype)

    if (dumpLocal) then
       write(*,"(a,i4,a,i4,a,i4,a,i4,a,i4,a,i4,a,i4,a,i4,a)") &
            h//trim(fieldName)//"(   1: ",nz, ",   1:",ldimx, ",   1:",ldimy, &
            ")=fullGrid(   1: ",nz,",",ia,":",iz,",",ja,":",jz,")"
    end if
  end subroutine ReadStoreOwnChunk_3D





  subroutine AllocReadStoreOwnChunk_3D(fUnit, toStore, fieldName, &
       nz, nnxp, nnyp, ldimx, ldimy, ia, iz, ja, jz, &
       mchnum, master_num, runtype)

    integer, intent(in) :: fUnit
    real, pointer :: toStore(:,:,:)
    character(len=*), intent(in) :: fieldName
    integer, intent(in) :: nz
    integer, intent(in) :: nnxp
    integer, intent(in) :: nnyp
    integer, intent(in) :: ldimx
    integer, intent(in) :: ldimy
    integer, intent(in) :: ia
    integer, intent(in) :: iz
    integer, intent(in) :: ja
    integer, intent(in) :: jz
    integer, intent(in) :: mchnum
    integer, intent(in) :: master_num
    character(len=*), intent(in) :: runtype

    real :: fullGrid(nz,nnxp,nnyp)
    integer :: ierr
    character(len=8) :: c0, c1
    character(len=*), parameter :: h="**(AllocReadStoreOwnChunk_3D)**"


    ! master process opens file and reads first data into full domain scratch

    if (mchnum == master_num) then
       if (useVfm) then
          call vfirec(fUnit,fullGrid(1,1,1),nz*nnxp*nnyp,'LIN')
       else
          read(fUnit) fullGrid
       end if
    end if

    ! broadcast full domain scratch; 
    ! local chunk is extracted and stored at desired variable

    call parf_bcast(fullGrid, int(nz,i8), int(nnxp,i8), &
         int(nnyp,i8), master_num)

    call mk_3_buff(fullGrid(1,1,1), toStore(1,1,1), &
         nz, nnxp, nnyp, nz, ldimx, ldimy, ia, iz, ja, jz)

    if (dumpLocal) then
       write(*,"(a,i4,a,i4,a,i4,a,i4,a,i4,a,i4,a,i4,a,i4,a)") &
            h//trim(fieldName)//"(   1: ",nz, ",   1:",ldimx, ",   1:",ldimy, &
            ")=fullGrid(   1: ",nz,",",ia,":",iz,",",ja,":",jz,")"
    end if
  end subroutine AllocReadStoreOwnChunk_3D


  subroutine ReadStoreFullFieldAndOwnChunk(grid, fUnit, full, toStore, fieldName)
    integer, intent(in) :: grid
    integer, intent(in) :: fUnit
    real, pointer :: full(:,:)
    real, pointer :: toStore(:,:)
    character(len=*), intent(in) :: fieldName


    integer :: ldimx, ldimy, lin
    integer :: ia, iz, ja, jz
    integer :: ierr
    character(len=8) :: c0, c1
    character(len=*), parameter :: h="**(ReadStoreFullFieldAndOwnChunk)**"

    ! check full field allocated memory

    !print *,'LFR-DBG Inside read Bor: ', grid,fUnit,size(full,1),size(full,2),size(toStore,1),size(toStore,2),fieldName
    !print *,'LFR-DBG Inside read Nodemxp: ',mynum,grid,nodemxp(mynum,grid); call flush(6)

    if (.not. associated(full)) then
       call fatal_error(h//" will store at not associated pointer var "//trim(fieldName))
    end if
    ldimx = nnxp(grid)
    lin = size(full,1)
    if (ldimx /= lin) then
       write(c0,"(i8)") ldimx
       write(c1,"(i8)") lin
       call fatal_error(h//" full field of "//trim(fieldName)//" x dimension ("//trim(adjustl(c1))//&
            ") differs from required ("//trim(adjustl(c0))//")")
    end if
    ldimy = nnyp(grid)
    lin = size(full,2)
    if (ldimy /= lin) then
       write(c0,"(i8)") ldimy
       write(c1,"(i8)") lin
       call fatal_error(h//" full field of "//trim(fieldName)//" y dimension ("//trim(adjustl(c1))//&
            ") differs from required ("//trim(adjustl(c0))//")")
    end if

    ! check local chunk allocated memory

    if (.not. associated(toStore)) then
       call fatal_error(h//" will store at not associated pointer var "//trim(fieldName))
    end if
    ldimx = nodemxp(mynum,grid)
    lin = size(toStore,1)
    if (ldimx /= lin) then
       write(c0,"(i8)") ldimx
       write(c1,"(i8)") lin
       call fatal_error(h//trim(fieldName)//" x dimension ("//trim(adjustl(c1))//&
            ") differs from required ("//trim(adjustl(c0))//")")
    end if
    ldimy = nodemyp(mynum,grid)
    lin = size(toStore,2)
    if (ldimy /= lin) then
       write(c0,"(i8)") ldimy
       write(c1,"(i8)") lin
       call fatal_error(h//trim(fieldName)//" y dimension ("//trim(adjustl(c1))//&
            ") differs from required ("//trim(adjustl(c0))//")")
    end if
    ia = nodei0(mynum,grid)+1
    iz = nodei0(mynum,grid)+nodemxp(mynum,grid)
    ja = nodej0(mynum,grid)+1
    jz = nodej0(mynum,grid)+nodemyp(mynum,grid)

    ! master process opens file and reads first data into full domain
!print *, 'LFR-DBG->',fieldName,'Reading... vfirec',mchnum 
    if (mchnum == master_num) then
       if (useVfm) then
          call vfirec(fUnit,full(1,1),nnxp(grid)*nnyp(grid),'LIN')
       else
          read(fUnit) full
       end if
    end if
!print *,'FieldName, Max e min lido: ',fieldName,maxval(full),minval(full),mchnum,master_num
    ! broadcast full domain; 
    ! local chunk is extracted and stored at desired variable

    call parf_bcast(full, int(nnxp(grid),i8), int(nnyp(grid),i8), &
         master_num)

    call mk_2_buff(full(1,1), toStore(1,1), &
         nnxp(grid), nnyp(grid), ldimx, ldimy, ia, iz, ja, jz)

    if (dumpLocal) then
       write(*,"(a,i4,a,i4,a,i4,a,i4,a,i4,a,i4,a)") &
            h//trim(fieldName)//"(   1:",ldimx, ",   1:",ldimy, &
            ")=full(",ia,":",iz,",",ja,":",jz,")"
    end if
  end subroutine ReadStoreFullFieldAndOwnChunk




  subroutine Broadcast_I(dataBcst, root, dataName)
    integer, intent(inout) :: dataBcst
    integer, intent(in) :: root
    character(len=*), intent(in) :: dataName
    integer :: ierr
    integer :: IntArr(1)
    character(len=8) :: c0, c1
    character(len=*), parameter :: h="**(Broadcast_I)**"

    if (dumpLocal) then
       write(c0, "(i8)") root
       write(c1, "(i8)") dataBcst
       write(*,"(a)") h//" process "//trim(adjustl(c0))//&
            " before broadcast "//trim(dataName)//" with value "//&
            trim(adjustl(c1))
    endif

    IntArr(1) = dataBcst
    call parf_bcast(IntArr, 1_i8, root)
    dataBcst = IntArr(1)

    if (dumpLocal) then
       write(c0, "(i8)") root
       write(c1, "(i8)") dataBcst
       write(*,"(a)") h//" process "//trim(adjustl(c0))//&
            " broadcast "//trim(dataName)//" with value "//&
            trim(adjustl(c1))
    end if
  end subroutine Broadcast_I






  subroutine Broadcast_I1D(dataBcst, root, dataName)
    integer, intent(inout) :: dataBcst(:)
    integer, intent(in) :: root
    character(len=*), intent(in) :: dataName
    integer :: ierr
    character(len=8) :: c0, c1
    character(len=*), parameter :: h="**(Broadcast_I1D)**"

    if (dumpLocal) then
       write(c0, "(i8)") root
       write(c1, "(i8)") size(dataBcst)
       write(*,"(a)") h//" process "//trim(adjustl(c0))//&
            " before broadcast "//trim(dataName)//" with size "//&
            trim(adjustl(c1))
       call flush(6)
    endif

    call parf_bcast(dataBcst, int(size(dataBcst),i8), root)

    if (dumpLocal) then
       write(c0, "(i8)") root
       write(*,"(a)") h//" process "//trim(adjustl(c0))//&
            " broadcast "//trim(dataName)
    end if
  end subroutine Broadcast_I1D






  subroutine Broadcast_R(dataBcst, root, dataName)
    real, intent(inout) :: dataBcst
    integer, intent(in) :: root
    character(len=*), intent(in) :: dataName
    integer :: ierr
    real :: RealArr
    character(len=8) :: c0, c1
    character(len=*), parameter :: h="**(Broadcast_R)**"

    RealArr = dataBcst
    call parf_bcast(RealArr, root)
    dataBcst = RealArr

    if (dumpLocal) then
       write(c0, "(i8)") root
       write(c1, "(f8.4)") dataBcst
       write(*,"(a)") h//" process "//trim(adjustl(c0))//&
            " broadcast "//trim(dataName)//" with value "//&
            trim(adjustl(c1))
    end if
  end subroutine Broadcast_R






  subroutine Broadcast_R1D(dataBcst, root, dataName)
    real, intent(inout)          :: dataBcst(:)
    integer, intent(in)          :: root
    character(len=*), intent(in) :: dataName
    integer :: ierr
    character(len=8) :: c0, c1
    character(len=*), parameter :: h="**(Broadcast_R1D)**"

    ! broadcast dataBcst

    call parf_bcast(dataBcst, int(size(dataBcst),i8), root)

    if (dumpLocal) then
       write(c0, "(i8)") root
       write(*,"(a)") h//" process "//trim(adjustl(c0))//&
            " broadcast "//trim(dataName)
    end if
  end subroutine Broadcast_R1D






  subroutine Broadcast_R2D(dataBcst, root, dataName)
    real, intent(inout)          :: dataBcst(:,:)
    integer, intent(in)          :: root
    character(len=*), intent(in) :: dataName
    integer :: ierr
    character(len=8) :: c0, c1
    character(len=*), parameter :: h="**(Broadcast_R2D)**"

    ! broadcast dataBcst

    call parf_bcast(dataBcst, int(size(dataBcst,1),i8), &
         int(size(dataBcst,2),i8), root)

    if (dumpLocal) then
       write(c0, "(i8)") root
       write(*,"(a)") h//" process "//trim(adjustl(c0))//&
            " broadcast "//trim(dataName)
    end if
  end subroutine Broadcast_R2D






  subroutine Broadcast_C(dataBcst, root, dataName)
    character(len=*), intent(inout) :: dataBcst
    integer, intent(in) :: root
    character(len=*), intent(in) :: dataName
    integer :: ierr
    character(len=8) :: c0
    character(len=*), parameter :: h="**(Broadcast_C)**"

    call parf_bcast(dataBcst, int(len(dataBcst),i8), root)

    if (dumpLocal) then
       write(c0, "(i8)") root
       write(*,"(a)") h//" process "//trim(adjustl(c0))//&
            " broadcast "//trim(dataName)
    end if
  end subroutine Broadcast_C






  subroutine Broadcast_C1D(dataBcst, root, dataName)
    character(len=*), intent(inout) :: dataBcst(:)
    integer, intent(in) :: root
    character(len=*), intent(in) :: dataName
    integer :: ierr
    character(len=8) :: c0
    character(len=*), parameter :: h="**(Broadcast_C1D)**"

    call parf_bcast(dataBcst, int(len(dataBcst),i8), int(size(dataBcst),i8), root)

    if (dumpLocal) then
       write(c0, "(i8)") root
       write(*,"(a)") h//" process "//trim(adjustl(c0))//&
            " broadcast "//trim(dataName)
    end if
  end subroutine Broadcast_C1D





  subroutine ProcWithMin(val, rank)
    real,    intent(inout) :: val
    integer, intent(out)   :: rank

    integer :: ierr
    real :: bufin(2)
    real :: bufout(2)

    bufin(:) = (/ val, real(mchnum) /)

    call parf_minloc(bufin, bufout)

    val = bufout(1)
    rank = NINT(bufout(2))
  end subroutine ProcWithMin






  subroutine MaxCFLOverall(cflxy, cflz)
    real, intent(inout) :: cflxy(:)
    real, intent(inout) :: cflz(:)

    integer :: sizeCflxy
    integer :: sizeCflz
    integer :: sizeVec
    integer :: ierr
    real, allocatable :: vecIn(:)
    real, allocatable :: vecOut(:)
    character(len=8) :: c0, c1
    character(len=*), parameter :: h="**(MaxCFLOverall)**"

    sizeCflxy = size(cflxy)
    sizeCflz  = size(cflz)

    if (sizeCflxy /= sizeCflz) then
       write(c0,"(i8)") sizeCflxy
       write(c1,"(i8)") sizeCflz
       call fatal_error(h//" sizes of cflxy ("//trim(adjustl(c0))//&
            ") and cflz ("//trim(adjustl(c1))//" disagree")
    end if

    sizeVec=sizeCflxy+sizeCflz

    allocate(vecOut(sizeVec), stat=ierr)
    if (ierr /= 0) then
       write(c0,"(i8)") sizeVec
       write(c1,"(i8)") ierr
       call fatal_error(h//" allocate vecOut ("//trim(adjustl(c0))//&
            ") failed with stat="//trim(adjustl(c1)))
    end if

    allocate(vecIn(sizeVec), stat=ierr)
    if (ierr /= 0) then
       write(c0,"(i8)") sizeVec
       write(c1,"(i8)") ierr
       call fatal_error(h//" allocate vecIn ("//trim(adjustl(c0))//&
            ") failed with stat="//trim(adjustl(c1)))
    end if

    vecIn(1:sizeCflxy) = cflxy
    vecIn(sizeCflxy+1:sizeCflxy+sizeCflz) = cflz

    call parf_allreduce_max(vecIn, vecOut, int(sizeVec,i8))

    cflxy = vecOut(1:sizeCflxy)
    cflz  = vecOut(sizeCflxy+1:sizeCflxy+sizeCflz)

    deallocate(vecOut, stat=ierr)
    if (ierr /= 0) then
       write(c1,"(i8)") ierr
       call fatal_error(h//" deallocate vecOut failed with stat="//trim(adjustl(c1)))
    end if

    deallocate(vecIn, stat=ierr)
    if (ierr /= 0) then
       write(c1,"(i8)") ierr
       call fatal_error(h//" deallocate vecIn failed with stat="//trim(adjustl(c1)))
    end if
  end subroutine MaxCFLOverall



  ! LocalSizesAndDisp : sizes and displacements required by MPI_GATHERV,
  !                      indexed by process id (1:nmachs) and by variable 
  !                      idim_type


  subroutine LocalSizesAndDisp (ngrid, il1, ir2, jb1, jt2, localSize, disp)

    integer, intent(in ) :: ngrid            ! which grid to use
    integer, intent(out) :: il1(nmachs)
    integer, intent(out) :: ir2(nmachs)
    integer, intent(out) :: jb1(nmachs)
    integer, intent(out) :: jt2(nmachs)
    integer, intent(out) :: localSize(nmachs, idim_type_min:idim_type_max)
    integer, intent(out) :: disp(nmachs, idim_type_min:idim_type_max)

    integer :: proc
    integer :: idim_type
    character(len=8) :: c0
    character(len=*), parameter :: h="**(LocalSizesAndDisp)**" 

    ! field size at each process and dimensionality

    do proc = 1, nmachs
       localSize(proc,2) = nodemxp(proc,ngrid)*nodemyp(proc,ngrid)
       localSize(proc,3) = nnzp(ngrid)*nodemxp(proc,ngrid)*nodemyp(proc,ngrid)
       localSize(proc,4) = nzg*nodemxp(proc,ngrid)*nodemyp(proc,ngrid)*npatch
       localSize(proc,5) = nzs*nodemxp(proc,ngrid)*nodemyp(proc,ngrid)*npatch
       localSize(proc,6) = nodemxp(proc,ngrid)*nodemyp(proc,ngrid)*npatch
       localSize(proc,7) = nodemxp(proc,ngrid)*nodemyp(proc,ngrid)*nwave
    end do

    ! displacement from base address of gathered 1D array 
    ! where to start gathering field portion from each process

    do idim_type = idim_type_min, idim_type_max
       disp(1,idim_type) = 0   ! first process at base address
       do proc = 2, nmachs   ! next process at (proc-1) position + localSize(proc-1)
          disp(proc,idim_type) = &
               disp(proc-1,idim_type)+localSize(proc-1,idim_type)
       end do
    end do

    ! gathered 1D array will be unpacked; in this process,
    ! eliminate parallel ghost zones, marking start and end
    ! indices of internal zones

    do proc = 1, nmachs
       if (btest(nodeibcon(proc,ngrid),0)) then
          il1(proc) = nodeia(proc,ngrid) - 1
       else
          il1(proc) = nodeia(proc,ngrid)
       end if
       if (btest(nodeibcon(proc,ngrid),1)) then
          ir2(proc) = nodeiz(proc,ngrid) + 1
       else
          ir2(proc) = nodeiz(proc,ngrid)
       end if
       if (btest(nodeibcon(proc,ngrid),2)) then
          jb1(proc) = nodeja(proc,ngrid) - 1
       else
          jb1(proc) = nodeja(proc,ngrid)
       end if
       if (btest(nodeibcon(proc,ngrid),3)) then
          jt2(proc) = nodejz(proc,ngrid) + 1
       else
          jt2(proc) = nodejz(proc,ngrid)
       end if
    end do
  end subroutine LocalSizesAndDisp



  ! PreProcAndGather: returns field gathered from all processes (at master process)
  !                   pre-process field before gathering if required.



  subroutine PreProcAndGather(preProc, ngrid, idim_type, varn, &
       il1, ir2, jb1, jt2, localSize, disp, thisChunkSize, LocalChunk, &
       sizeGathered, gathered, sizeFullField, FullField)

    logical,          intent(in   ) :: preProc
    integer,          intent(in   ) :: ngrid
    integer,          intent(in   ) :: idim_type
    character(len=*), intent(inout) :: varn
    integer,          intent(in   ) :: il1(nmachs)
    integer,          intent(in   ) :: ir2(nmachs)
    integer,          intent(in   ) :: jb1(nmachs)
    integer,          intent(in   ) :: jt2(nmachs)
    integer,          intent(in   ) :: localSize(nmachs, idim_type_min:idim_type_max)
    integer,          intent(in   ) :: disp(nmachs, idim_type_min:idim_type_max)
    integer,          intent(in   ) :: thisChunkSize
    real,             intent(inout) :: LocalChunk(thisChunkSize)
    integer,          intent(in   ) :: sizeGathered
    real,             intent(out  ) :: gathered(sizeGathered)    !scratch
    integer,          intent(in   ) :: sizeFullField
    real,             intent(out  ) :: FullField(sizeFullField)

    character(len=len(varn)) :: varnOut
    integer :: ierr
    character(len=7)            :: cProc
    character(len=8)            :: c0, c1, c2
    character(len=*), parameter :: h="**(PreProcAndGather)**" 


    write(cProc,"(a5,i2)") " Proc",mchnum

    ! if requiring pre-processing

    if (preProc) then

       if (dumpLocal) then
          write(c0,"(i8)") thisChunkSize
          write(*,"(a)") h//cProc//"  with thisChunkSize="//trim(adjustl(c0))
       end if

       ! pre-process LocalChunk before gathering

       call PreProcForOutput(ngrid, varn, thisChunkSize, LocalChunk, &
            LocalChunk, varnOut)       
       varn = varnOut


    end if

    ! gather field scaterred over slaves; master gets full unpacked field

    call GatherOneFullField(ngrid, idim_type, thisChunkSize, &
         LocalChunk, localSize, disp, il1, ir2, jb1, jt2, &
         sizeGathered, gathered, sizeFullField, FullField)
  end subroutine PreProcAndGather






  subroutine GatherOneFullField(ngrid, idim_type, thisChunkSize, &
       LocalChunk, localSize, disp, il1, ir2, jb1, jt2, &
       sizeGathered, gathered, sizeFullField, FullField)
    integer, intent(in ) :: ngrid
    integer, intent(in ) :: idim_type
    integer, intent(in ) :: thisChunkSize
    real,    intent(in ) :: LocalChunk(thisChunkSize)
    integer, intent(in ) :: localSize(nmachs, idim_type_min:idim_type_max)
    integer, intent(in ) :: disp(nmachs, idim_type_min:idim_type_max)
    integer, intent(in ) :: il1(nmachs)
    integer, intent(in ) :: ir2(nmachs)
    integer, intent(in ) :: jb1(nmachs)
    integer, intent(in ) :: jt2(nmachs)
    integer, intent(in ) :: sizeGathered
    real,    intent(out) :: gathered(sizeGathered)    !scratch
    integer, intent(in ) :: sizeFullField
    real,    intent(out) :: FullField(sizeFullField)

    integer :: ierr
    character(len=8) :: c0, c1
    character(len=*), parameter :: h="**(GatherOneFullField)**"

    ! gather field (with ghost zones and wrong order) at master_num

    if (dumpLocal) then
       write(c0,"(i8)") thisChunkSize
       write(*,"(a)") h//" will gather with local size "//trim(adjustl(c0))
       call flush(6)
    end if
    call GatherAllChunks (LocalChunk, thisChunkSize, idim_type, &
         localSize, disp, gathered, sizeGathered)
    if (dumpLocal) then
       write(*,"(a)") h//" done gathering"
    end if

    ! master_num unpacks fields (removes unnecessary ghost zones and positions entries)

    if (mchnum == master_num) then
       if (dumpLocal) then
          write(*,"(a)") h//" master will RemoveGhost"
       end if
       call RemoveGhost (ngrid, idim_type, sizeGathered, gathered, &
            il1, ir2, jb1, jt2, disp, localSize, sizeFullField, FullField)
       if (dumpLocal) then
          write(*,"(a)") h//" done RemoveGhost"
       end if
    end if
  end subroutine GatherOneFullField






  subroutine GatherAllChunks (LocalChunk, thisChunkSize, idim_type, &
       localSize, disp, gathered, sizeGathered)
    integer, intent(in ) :: thisChunkSize
    real,    intent(in ) :: LocalChunk(thisChunkSize)
    integer, intent(in ) :: idim_type
    integer, intent(in ) :: localSize(nmachs, idim_type_min:idim_type_max)
    integer, intent(in ) :: disp(nmachs, idim_type_min:idim_type_max)
    integer, intent(in ) :: sizeGathered
    real,    intent(out) :: gathered(sizeGathered)

    integer :: ierr
    character(len=8) :: c0, c1
    character(len=*), parameter :: h="**(GatherAllChunks)**" 

    if (dumpLocal) then
       write(c0,"(i8)") thisChunkSize
       write(c1,"(i8)") sizeGathered
       write(*,"(a)") h//" thisChunkSize="//trim(adjustl(c0))//&
            "; sizeGathered="//trim(adjustl(c1))
       call flush(6)
    end if

    ! gather a field

    call parf_GatherAllChunks(LocalChunk, thisChunkSize, idim_type, &
         localSize, disp, gathered, sizeGathered, master_num, nmachs)
    if (dumpLocal) then
       write(*,"(a)") h//" done"
       call flush(6)
    end if

  end subroutine GatherAllChunks







  subroutine RemoveGhost (ngrid, idim_type, sizeGathered, gathered, &
       il1, ir2, jb1, jt2, disp, localSize, sizeFullField, FullField)
    integer, intent(in)  :: ngrid
    integer, intent(in)  :: idim_type
    integer, intent(in)  :: sizeGathered
    real,    intent(in)  :: gathered(sizeGathered)
    integer, intent(in)  :: il1(nmachs)
    integer, intent(in)  :: ir2(nmachs)
    integer, intent(in)  :: jb1(nmachs)
    integer, intent(in)  :: jt2(nmachs)
    integer, intent(in)  :: disp(nmachs, idim_type_min:idim_type_max)
    integer, intent(in)  :: localSize(nmachs, idim_type_min:idim_type_max)
    integer, intent(in)  :: sizeFullField
    real,    intent(out) :: FullField(sizeFullField)

    integer :: proc
    character(len=8) :: c0
    character(len=*), parameter :: h="**(RemoveGhost)**"

    ! unpack gathered field, removing unnecessary ghost zones and
    ! placing entries at correct field positions

    select case (idim_type)
    case (2)
       do proc = 1, nmachs
          if (localSize(proc,idim_type)/=0) then
             call ex_2_buff(FullField, gathered(disp(proc,idim_type)+1), &
                  nnxp(ngrid), nnyp(ngrid), &
                  nodemxp(proc,ngrid), nodemyp(proc,ngrid), &
                  nodei0(proc,ngrid), nodej0(proc,ngrid), &
                  il1(proc), ir2(proc), jb1(proc), jt2(proc))
          end if
       end do
    case (3)
       do proc = 1, nmachs
          if (localSize(proc,idim_type)/=0) then
             call ex_3_buff(FullField, gathered(disp(proc,idim_type)+1), &
                  nnzp(ngrid), nnxp(ngrid), nnyp(ngrid), &
                  nnzp(ngrid), nodemxp(proc,ngrid), nodemyp(proc,ngrid), &
                  nodei0(proc,ngrid), nodej0(proc,ngrid), &
                  il1(proc), ir2(proc), jb1(proc), jt2(proc))
          end if
       end do
    case (4)
       do proc = 1, nmachs
          if (localSize(proc,idim_type)/=0) then
             call ex_4_buff(FullField, gathered(disp(proc,idim_type)+1), &
                  nzg, nnxp(ngrid), nnyp(ngrid), npatch, &
                  nzg, nodemxp(proc,ngrid), nodemyp(proc,ngrid), npatch, &
                  nodei0(proc,ngrid), nodej0(proc,ngrid), &
                  il1(proc), ir2(proc), jb1(proc), jt2(proc))
          end if
       end do
    case (5)
       do proc = 1, nmachs
          if (localSize(proc,idim_type)/=0) then
             call ex_4_buff(FullField, gathered(disp(proc,idim_type)+1), &
                  nzs, nnxp(ngrid), nnyp(ngrid), npatch, &
                  nzs, nodemxp(proc,ngrid), nodemyp(proc,ngrid), npatch, &
                  nodei0(proc,ngrid), nodej0(proc,ngrid), &
                  il1(proc), ir2(proc), jb1(proc), jt2(proc))
          end if
       end do
    case (6)
       do proc = 1, nmachs
          if (localSize(proc,idim_type)/=0) then
             call ex_2p_buff(FullField, gathered(disp(proc,idim_type)+1), &
                  nnxp(ngrid), nnyp(ngrid), npatch, &
                  nodemxp(proc,ngrid), nodemyp(proc,ngrid), npatch, &
                  nodei0(proc,ngrid), nodej0(proc,ngrid), &
                  il1(proc), ir2(proc), jb1(proc), jt2(proc))
          end if
       end do
    case (7)
       do proc = 1, nmachs
          if (localSize(proc,idim_type)/=0) then
             call ex_buff_carma(FullField, gathered(disp(proc,idim_type)+1), &
                  nnxp(ngrid), nnyp(ngrid), nwave, &
                  nodemxp(proc,ngrid), nodemyp(proc,ngrid), nwave, &
                  nodei0(proc,ngrid), nodej0(proc,ngrid), &
                  il1(proc), ir2(proc), jb1(proc), jt2(proc))
          end if
       end do
    case default
       write(c0,"(i8)") idim_type
       call fatal_error(h//" invoked with unknown idim_type ("//&
            trim(adjustl(c0))//")")
    end select
  end subroutine RemoveGhost



  ! RearrangeAndDump: dump field into all required output files;
  !                   rearrange field prior to output if required



  subroutine RearrangeAndDump (ngrid, rearran, &
       histFlag, instFlag, liteFlag, meanFlag, &
       histFileDS, instFileDS, liteFileDS, meanFileDS, &
       varn, idim_type, sizeFullField, FullField, Rear)

    integer,          intent(in   ) :: ngrid
    logical,          intent(in   ) :: rearran
    logical,          intent(in   ) :: histFlag
    logical,          intent(in   ) :: instFlag
    logical,          intent(in   ) :: liteFlag
    logical,          intent(in   ) :: meanFlag
    type(IOFileDS),   intent(inout) :: histFileDS
    type(IOFileDS),   intent(inout) :: instFileDS
    type(IOFileDS),   intent(inout) :: liteFileDS
    type(IOFileDS),   intent(inout) :: meanFileDS
    character(len=*), intent(in   ) :: varn
    integer,          intent(in   ) :: idim_type
    integer,          intent(in   ) :: sizeFullField
    real,             intent(in   ) :: FullField(sizeFullField)
    real,             intent(out  ) :: Rear(sizeFullField)

    character(len=*), parameter :: h="**(RearrangeAndDump)**"

    if (dumpLocal) then
       write(*,"(4(a,l1))") h//" enter field "//trim(varn)//":"//&
            " histFlag=",histFlag, &
            " instFlag=",instFlag, &
            " liteFlag=",liteFlag, &
            " meanFlag=",meanFlag
    end if

    ! if field to be rearranged, rearrange and dump

    if (rearran) then
       call RearrangeForOutput(nnxp(ngrid), nnyp(ngrid), nnzp(ngrid), &
            nzg, nzs, npatch, idim_type, FullField, Rear)

       ! dump rearranged field

       if (instFlag) then
          call FieldWriteStoreInfo(Rear, sizeFullField, &
               varn, idim_type, ngrid, instFileDS)
       end if
       if (liteFlag) then
          call FieldWriteStoreInfo(Rear, sizeFullField, &
               varn, idim_type, ngrid, liteFileDS)
       end if
       if (meanFlag) then
          call FieldWriteStoreInfo(Rear, sizeFullField, &
               varn, idim_type, ngrid, meanFileDS)
       end if
       if (histFlag) then
          call FieldWriteStoreInfo(Rear, sizeFullField, &
               varn, idim_type, ngrid, histFileDS)
       end if
    else

       ! dump original field

       if (instFlag) then
          call FieldWriteStoreInfo(FullField, sizeFullField, &
               varn, idim_type, ngrid, instFileDS)
       end if
       if (liteFlag) then
          call FieldWriteStoreInfo(FullField, sizeFullField, &
               varn, idim_type, ngrid, liteFileDS)
       end if
       if (meanFlag) then
          call FieldWriteStoreInfo(FullField, sizeFullField, &
               varn, idim_type, ngrid, meanFileDS)
       end if
       if (histFlag) then
          call FieldWriteStoreInfo(FullField, sizeFullField, &
               varn, idim_type, ngrid, histFileDS)
       end if
    end if
  end subroutine RearrangeAndDump



  ! DumpFullField: master_process gathers all local chunks of a field and dumps
  !                resulting full field at a file. If dumpGathered is selected,
  !                also dumps gathered field (that includes all ghost zones,
  !                excluded on full field)




  subroutine DumpFullField(ngrid, LocalChunk, name, dumpGathered)

    integer,          intent(in) :: ngrid               ! current grid number
    real,             pointer    :: LocalChunk(:,:,:)   ! field at this process    
    character(len=*), intent(in) :: name    ! field name (file name is <name>.<grid>.<nmachs>
    logical,          intent(in) :: dumpGathered        ! dumps gathered field after full field

    integer, parameter :: idim_type=3
    integer :: il1(nmachs)
    integer :: ir2(nmachs)
    integer :: jb1(nmachs)
    integer :: jt2(nmachs)
    integer :: localSize(nmachs, idim_type_min:idim_type_max)
    integer :: disp(nmachs, idim_type_min:idim_type_max)
    integer :: sizeFullField(idim_type_min:idim_type_max)
    character(len=4) :: cProcs
    character(len=256) :: line
    integer :: proc
    real, pointer :: var_p
    character(len=8) :: c0, c1
    character(len=*), parameter :: h="**(DumpFullField)**"

    if (dumpLocal) then
       write(c0,"(i8)") ngrid
       if (dumpGathered) then
          write(*,"(a)") h//" only dump FullField "//trim(name)//&
               " for grid "//trim(adjustl(c0))
       else
          write(*,"(a)") h//" dump FullField and gathered for field "//trim(name)//&
               " at grid "//trim(adjustl(c0))
       end if
    end if

    call LocalSizesAndDisp (ngrid, il1, ir2, jb1, jt2, localSize, disp)

    if (dumpLocal) then
       do proc = 1, nmachs
          write(line, "(5(a6,i8.8),a12,i8.8,a7,i8.8)") &
               " proc=", proc, &
               "; il1=", il1(proc), &
               "; ir2=", ir2(proc), &
               "; jb1=", jb1(proc), &
               "; jt2=", jt2(proc), &
               "; localSize=", localSize(proc,idim_type), &
               "; disp=", disp(proc,idim_type)
          write(*,"(a)") h//trim(line)
       end do
    end if

    call GlobalSizes(ngrid, nmachs, nwave, sizeFullField)

    var_p => LocalChunk(1,1,1)

    call DumpVTabEntry(ngrid, name, idim_type, var_p, &
         sizeFullField, il1, ir2, jb1, jt2, localSize, disp, &
         name, dumpGathered)

  end subroutine DumpFullField





  ! DumpVTabEntry: Gather and dump full field of one var_table entry into a file;
  !                If requested, also dumps all local chunks with ghost zones.




  subroutine DumpVTabEntry(ngrid, name, idim_type, var_p, &
       sizeFullField, il1, ir2, jb1, jt2, localSize, disp, &
       fPrefix, dumpGathered, unitOpened)

    
    integer,          intent(in) :: ngrid
    character(len=16), intent(in) :: name
    integer,          intent(in) :: idim_type
    real, pointer                :: var_p
    integer,          intent(in) :: sizeFullField(idim_type_min:idim_type_max)
    integer,          intent(in) :: il1(nmachs)
    integer,          intent(in) :: ir2(nmachs)
    integer,          intent(in) :: jb1(nmachs)
    integer,          intent(in) :: jt2(nmachs)
    integer,          intent(in) :: localSize(nmachs,idim_type_min:idim_type_max)
    integer,          intent(in) :: disp(nmachs,idim_type_min:idim_type_max)
    character(len=*), intent(in) :: fPrefix
    logical,          intent(in) :: dumpGathered
    integer, optional, intent(in) :: unitOpened

    character(len=2)            :: cGrid
    character(len=4)            :: cProc
    character(len=10)           :: cProcHeader
    character(len=8)            :: c0, c1, c2
    character(len=*), parameter :: h="**(DumpVTabEntry)**" 

    integer, parameter :: unitLow=10
    integer, parameter :: unitHigh=99
    integer :: sizeGathered
    integer :: sizeLocalChunk
    integer :: iunit
    integer :: ierr
    integer :: proc

    logical :: op

    real, allocatable :: LocalChunk(:)
    real, allocatable :: Gathered(:)
    real, allocatable :: FullField(:)

    ! consistency

    if (idim_type < idim_type_min .or. idim_type > idim_type_max) then
       write(c0,"(i8)") idim_type
       call fatal_error(h//" unknown idim_type="//trim(adjustl(c0)))
    end if

    ! names for dumping

    write(cProc,"(i4.4)")  mchnum
    cProcHeader = " Proc "//cProc
    write(cProc,"(i4.4)")  nmachs
    write(cGrid,"(a1,i1)") "g", ngrid

    ! find available Fortran unit

    if (present(unitOpened)) then
       iunit = unitOpened
    else if (mchnum == master_num) then
       do iunit = unitLow, unitHigh
          inquire (unit=iunit, opened=op)
          if (.not. op) exit
       end do
       if (iunit > unitHigh) then
          call fatal_error(h//" Fortran i/o units exausted")
       end if
    end if

    ! space to copy vtab entry

    sizeLocalChunk=localSize(mynum,idim_type)
    allocate(LocalChunk(sizeLocalChunk), stat=ierr)
    if (ierr /= 0) then
       write(c0,"(i8)") sizeLocalChunk
       write(c1,"(i8)") ierr
       call fatal_error(h//" allocate localSize("// &
            trim(adjustl(c0))//") failed with stat="// &
            trim(adjustl(c1)))
    else if (dumpLocal) then
       write(c0,"(i8)") sizeLocalChunk
       write(*,"(a)") h//cProcHeader//" allocated localSize("//trim(adjustl(c0))//")"
       call flush(6)
    end if

    ! space for gathered full field

    sizeGathered = disp(nmachs,idim_type) + localSize(nmachs,idim_type)
    allocate(Gathered(sizeGathered), stat=ierr)
    if (ierr /= 0) then
       write(c0,"(i8)") sizeGathered
       write(c1,"(i8)") ierr
       call fatal_error(h//" allocate Gathered("// &
            trim(adjustl(c0))//") failed with stat="// &
            trim(adjustl(c1)))
    else if (dumpLocal) then
       write(c0,"(i8)") sizeGathered
       write(*,"(a)") h//cProcHeader//" allocated Gathered("//trim(adjustl(c0))//")"
       call flush(6)
    end if

    ! space for unpacked field

    allocate(FullField(sizeFullField(idim_type)), stat=ierr)
    if (ierr /= 0) then
       write(c0,"(i8)") sizeFullField(idim_type)
       write(c1,"(i8)") ierr
       call fatal_error(h//" allocate FullField("// &
            trim(adjustl(c0))//") failed with stat="// &
            trim(adjustl(c1)))
    else if (dumpLocal) then
       write(c0,"(i8)") sizeFullField(idim_type)
       write(*,"(a)") h//cProcHeader//" allocated FullField("//trim(adjustl(c0))//")"
       call flush(6)
    end if

    ! copy vtab entry into a 1D array

    call CopyLocalChunk(var_p, LocalChunk, sizeLocalChunk)

    ! gather local chunks into a single array at master process (on "Gathered") 
    ! and eliminate ghost zones of the gathered array (on "FullField")

    call GatherOneFullField(ngrid, idim_type, sizeLocalChunk, &
         LocalChunk, localSize, disp, il1, ir2, jb1, jt2, &
         sizeGathered, Gathered, sizeFullField(idim_type), FullField)

    ! only master dumps file

    if (mchnum == master_num) then

       ! open file

       if (.not. present(unitOpened)) then
          open(iunit, file=fPrefix(1:len_trim(fPrefix))//"."//cGrid//"."//cProc, status="replace", &
               action="write", form="unformatted", iostat=ierr)
          if (ierr /= 0) then
             write(c0,"(i8)") ierr
             call fatal_error(h//" opening "//trim(fPrefix)//"."//cGrid//"."//cProc//&
                  " failed with iostat="//trim(adjustl(c0)))
          end if
       end if

       ! dump name, idim_type and dimensions

       write(iunit) name
       write(iunit) idim_type
       select case (idim_type)
       case (2)
          write(iunit) nnxp(ngrid), nnyp(ngrid)
       case (3)
          write(iunit) nnzp(ngrid), nnxp(ngrid), nnyp(ngrid)
       case (4)
          write(iunit) nzg, nnxp(ngrid), nnyp(ngrid), npatch
       case (5)
          write(iunit) nzs, nnxp(ngrid), nnyp(ngrid), npatch
       case (6)
          write(iunit) nnxp(ngrid), nnyp(ngrid), npatch
       case (7)
          write(iunit) nnxp(ngrid), nnyp(ngrid), nwave
       case default
          write(c0,"(i8)") idim_type
          call fatal_error(h//" unknown idim_type="//trim(adjustl(c0)))
       end select

       ! dump Full Field

       write(iunit) FullField

       ! dump flag for dumpGathered

       write(iunit) dumpGathered

       ! dump Gathered (include ghost zones)

       if (dumpGathered) then
          write(iunit) nmachs
          do proc = 1, nmachs
             write(iunit) localSize(proc,idim_type), nnzp(ngrid), &
                  nodemxp(proc,ngrid), nodemyp(proc,ngrid), &
                  nodei0(proc,ngrid), nodej0(proc,ngrid)
             write(iunit) Gathered(disp(proc,idim_type)+1:&
                  disp(proc,idim_type)+localSize(proc,idim_type))
          end do
       end if

       ! close file

       if (.not. present(unitOpened)) then
          close(iunit)
       end if
    end if

    ! deallocated space

    deallocate(LocalChunk, stat=ierr)
    if (ierr /= 0) then
       write(c1,"(i8)") ierr
       call fatal_error(h//" deallocate localSize failed with stat="// &
            trim(adjustl(c1)))
    end if

    deallocate(Gathered, stat=ierr)
    if (ierr /= 0) then
       write(c1,"(i8)") ierr
       call fatal_error(h//" deallocate Gathered failed with stat="// &
            trim(adjustl(c1)))
    end if

    ! maximum space for unpacked field

    deallocate(FullField, stat=ierr)
    if (ierr /= 0) then
       write(c1,"(i8)") ierr
       call fatal_error(h//" deallocate FullField failed with stat="// &
            trim(adjustl(c1)))
    end if
  end subroutine DumpVTabEntry



  ! DumpVarTables: dump all var_table entries for debugging



  subroutine DumpVarTables(fPrefix,dumpGathered)

    include "files.h"

    character(len=*), intent(in) :: fPrefix  ! file prefix
    logical, intent(in) :: dumpGathered

    integer :: ng
    integer :: nv
    integer :: global
    integer :: sizeFullField(idim_type_min:idim_type_max)
    integer :: il1(nmachs)
    integer :: ir2(nmachs)
    integer :: jb1(nmachs)
    integer :: jt2(nmachs)
    integer :: localSize(nmachs,idim_type_min:idim_type_max)
    integer :: disp(nmachs,idim_type_min:idim_type_max)

    integer :: ierr
    integer :: unitOpened
    integer, parameter :: unitLow=10
    integer, parameter :: unitHigh=99
    logical :: op
    character(len=2)   :: cGrid
    character(len=4)   :: cProc
    character(len=8)   :: c0, c1
    character(len=*), parameter :: h="**(DumpVarTables)**"

    ! find available Fortran unit

    do unitOpened = unitLow, unitHigh
       inquire (unit=unitOpened, opened=op)
       if (.not. op) exit
    end do
    if (unitOpened > unitHigh) then
       call fatal_error(h//" Fortran i/o units exausted")
    end if

    ! visit all grids

    do ng = 1, ngrids

       call GlobalSizes(ng, nmachs, nwave, sizeFullField)
       call LocalSizesAndDisp (ng, il1, ir2, jb1, jt2, localSize, disp)

       ! build file name

       write(cProc,"(i4.4)")  nmachs
       write(cGrid,"(a1,i1)") "g", ng

       ! one process opens file

       if (mchnum == master_num) then
          open(unitOpened, file=fPrefix(1:len_trim(fPrefix))//"."//cGrid//"."//cProc, status="replace", &
               action="write", form="unformatted", iostat=ierr)
          if (ierr /= 0) then
             write(c0,"(i8)") ierr
             call fatal_error(h//" opening "//trim(fPrefix)//"."//cGrid//"."//cProc//&
                  " failed with iostat="//trim(adjustl(c0)))
          end if
       end if

       ! visit all vtables entries

       do nv = 1, num_var(ng)
          call DumpVTabEntry(ng, vtab_r(nv,ng)%name, &
               vtab_r(nv,ng)%idim_type, vtab_r(nv,ng)%var_p, &
               sizeFullField, il1, ir2, jb1, jt2, localSize, disp, &
               fPrefix, dumpGathered, unitOpened)
       end do

       ! one process closes file

       if (mchnum == master_num) then
          close(unitOpened)
       end if
    end do
  end subroutine DumpVarTables



  ! Recreating Global Information (Gathering data)
  SUBROUTINE gatherData2D(idim_type, varn, ifm, nnxp, nnyp, &
       nmachs, mchnum, mynum, master_num,                   &
       localData2D, globalData2D)

    USE mem_grid, ONLY : &
         GlobalSizes         ! Subroutine

    USE ParLib, ONLY: &
         parf_bcast ! Subroutine

    IMPLICIT NONE
    INCLUDE "i8.h"
    ! Arguments:
    INTEGER, INTENT(IN)           :: idim_type, ifm, nnxp, nnyp, &
         nmachs, mchnum, mynum, master_num
    CHARACTER(LEN=16), INTENT(IN) :: varn
    REAL, INTENT(IN)              :: localData2D(:,:)
    REAL, INTENT(OUT)             :: globalData2D(:,:)
    ! Local Variables:
    CHARACTER(LEN=16)  :: localVarn
    INTEGER            :: ierr
    INTEGER, PARAMETER :: idim_type_min = 2
    INTEGER, PARAMETER :: idim_type_max = 7
    INTEGER            :: il1(nmachs)
    INTEGER            :: ir2(nmachs)
    INTEGER            :: jb1(nmachs)
    INTEGER            :: jt2(nmachs)
    INTEGER            :: localSize(nmachs,idim_type_min:idim_type_max)
    INTEGER            :: disp(nmachs,idim_type_min:idim_type_max)
    INTEGER            :: maxLocalSize
    INTEGER            :: sizeGathered(idim_type_min:idim_type_max)
    INTEGER            :: maxSizeGathered
    INTEGER            :: sizeFullField(idim_type_min:idim_type_max)
    INTEGER            :: maxsizeFullField
    INTEGER            :: globalSize(idim_type_min:idim_type_max)
    REAL, ALLOCATABLE  :: localChunk(:)
    REAL, ALLOCATABLE  :: gathered(:)
    REAL, ALLOCATABLE  :: fullField(:)

    ! Recreating Global information about Soil Water
    ! grid dependent, field independent constants for gather and unpacking
    ! as a function of idim_type
    CALL LocalSizesAndDisp(ifm, il1, ir2, jb1, jt2, localSize, disp)
    maxLocalSize = MAXVAL(localSize(mynum,:))
    ALLOCATE(localChunk(maxLocalSize), stat=ierr)
    IF (ierr/=0) THEN
       CALL fatal_error("Error allocating localChunk (gatherData)")
    ENDIF
    CALL CopyLocalChunk(localData2D(1,1), localChunk, &
         LocalSize(mynum,idim_type))
    sizeGathered(:) = disp(nmachs,:) + localSize(nmachs,:)
    maxSizeGathered = MAXVAL(sizeGathered)
    ALLOCATE(gathered(maxSizeGathered), stat=ierr)
    IF (ierr/=0) THEN
       CALL fatal_error("Error allocating gathered (gatherData)")
    ENDIF
    ! grid dependent field sizes as a function of idim_type
    CALL GlobalSizes(ifm, nmachs, nwave, globalSize)
    IF (mchnum==master_num) THEN
       sizeFullField(:) = globalSize(:)
    ELSE
       sizeFullField(:) = 1
    END IF
    maxSizeFullField = MAXVAL(sizeFullField)
    ALLOCATE(fullField(sizeFullField(idim_type)), stat=ierr)
    IF (ierr/=0) THEN
       CALL fatal_error("Error allocating fullField (gatherData)")
    ENDIF
    localVarn = trim(varn)
    CALL PreProcAndGather(.FALSE., ifm, idim_type, localVarn, &
         il1, ir2, jb1, jt2, localSize, disp,                 &
         localSize(mynum,idim_type), LocalChunk,              &
         sizeGathered(idim_type), gathered,                   &
         sizeFullField(idim_type), fullField                  )
    
    IF (mchnum==master_num) THEN
       globalData2D = RESHAPE(fullField, (/nnxp, nnyp/))
    ENDIF

    call parf_bcast(globalData2D, int(nnxp,i8), int(nnyp,i8), master_num)

    DEALLOCATE(fullField)
    DEALLOCATE(gathered)
    DEALLOCATE(localChunk)
    
  END SUBROUTINE gatherData2D



  ! Recreating Global Information (Gathering data)
  SUBROUTINE gatherData3D(idim_type, varn, ifm, nnzp, nnxp, nnyp, &
       nmachs, mchnum, mynum, master_num,                         &
       localData3D, globalData3D)

    USE mem_grid, ONLY : &
         GlobalSizes         ! Subroutine

    USE ParLib, ONLY: &
         parf_bcast ! Subroutine

    IMPLICIT NONE
    INCLUDE "i8.h"
    ! Arguments:
    INTEGER, INTENT(IN)           :: idim_type, ifm, nnzp, nnxp, nnyp, &
         nmachs, mchnum, mynum, master_num
    CHARACTER(LEN=16), INTENT(IN) :: varn
    REAL, INTENT(IN)              :: localData3D(:,:,:)
    REAL, INTENT(OUT)             :: globalData3D(:,:,:)
    ! Local Variables:
    CHARACTER(LEN=16)  :: localVarn
    INTEGER            :: ierr
    INTEGER, PARAMETER :: idim_type_min = 2
    INTEGER, PARAMETER :: idim_type_max = 7
    INTEGER            :: il1(nmachs)
    INTEGER            :: ir2(nmachs)
    INTEGER            :: jb1(nmachs)
    INTEGER            :: jt2(nmachs)
    INTEGER            :: localSize(nmachs,idim_type_min:idim_type_max)
    INTEGER            :: disp(nmachs,idim_type_min:idim_type_max)
    INTEGER            :: maxLocalSize
    INTEGER            :: sizeGathered(idim_type_min:idim_type_max)
    INTEGER            :: maxSizeGathered
    INTEGER            :: sizeFullField(idim_type_min:idim_type_max)
    INTEGER            :: maxsizeFullField
    INTEGER            :: globalSize(idim_type_min:idim_type_max)
    REAL, ALLOCATABLE  :: localChunk(:)
    REAL, ALLOCATABLE  :: gathered(:)
    REAL, ALLOCATABLE  :: fullField(:)

    ! Recreating Global information about Soil Water
    ! grid dependent, field independent constants for gather and unpacking
    ! as a function of idim_type
    CALL LocalSizesAndDisp(ifm, il1, ir2, jb1, jt2, localSize, disp)
    maxLocalSize = MAXVAL(localSize(mynum,:))
    ALLOCATE(localChunk(maxLocalSize), stat=ierr)
    IF (ierr/=0) THEN
       CALL fatal_error("Error allocating localChunk (gatherData)")
    ENDIF
    CALL CopyLocalChunk(localData3D(1,1,1), localChunk, &
         LocalSize(mynum,idim_type))
    sizeGathered(:) = disp(nmachs,:) + localSize(nmachs,:)
    maxSizeGathered = MAXVAL(sizeGathered)
    ALLOCATE(gathered(maxSizeGathered), stat=ierr)
    IF (ierr/=0) THEN
       CALL fatal_error("Error allocating gathered (gatherData)")
    ENDIF
    ! grid dependent field sizes as a function of idim_type
    CALL GlobalSizes(ifm, nmachs, nwave, globalSize)
    IF (mchnum==master_num) THEN
       sizeFullField(:) = globalSize(:)
    ELSE
       sizeFullField(:) = 1
    END IF
    maxSizeFullField = MAXVAL(sizeFullField)
    ALLOCATE(fullField(sizeFullField(idim_type)), stat=ierr)
    IF (ierr/=0) THEN
       CALL fatal_error("Error allocating fullField (gatherData)")
    ENDIF
    localVarn = trim(varn)
    CALL PreProcAndGather(.FALSE., ifm, idim_type, localVarn, &
         il1, ir2, jb1, jt2, localSize, disp,                 &
         localSize(mynum,idim_type), LocalChunk,              &
         sizeGathered(idim_type), gathered,                   &
         sizeFullField(idim_type), fullField                  )
    
    IF (mchnum==master_num) THEN
       globalData3D = RESHAPE(fullField, (/nnzp, nnxp, nnyp/))
    ENDIF

    call parf_bcast(globalData3D, &
         int(nnzp,i8), int(nnxp,i8), int(nnyp,i8), master_num)

    DEALLOCATE(fullField)
    DEALLOCATE(gathered)
    DEALLOCATE(localChunk)
    
  END SUBROUTINE gatherData3D



  ! Recreating Global Information (Gathering data)
  SUBROUTINE gatherData4D(idim_type, varn, ifm, mzg, nnxp, nnyp, npat, &
       nmachs, mchnum, mynum, master_num,                              &
       localData4D, globalData4D)

    USE mem_grid, ONLY : &
         GlobalSizes         ! Subroutine

    USE ParLib, ONLY: &
         parf_bcast ! Subroutine

    IMPLICIT NONE
    INCLUDE "i8.h"
    ! Arguments:
    INTEGER, INTENT(IN)           :: idim_type, ifm, mzg, nnxp, nnyp, npat, &
         nmachs, mchnum, mynum, master_num
    CHARACTER(LEN=16), INTENT(IN) :: varn
    REAL, INTENT(IN)              :: localData4D(:,:,:,:)
    REAL, INTENT(OUT)             :: globalData4D(:,:,:,:)
    ! Local Variables:
    CHARACTER(LEN=16)  :: localVarn
    INTEGER            :: ierr
    INTEGER, PARAMETER :: idim_type_min = 2
    INTEGER, PARAMETER :: idim_type_max = 7
    INTEGER            :: il1(nmachs)
    INTEGER            :: ir2(nmachs)
    INTEGER            :: jb1(nmachs)
    INTEGER            :: jt2(nmachs)
    INTEGER            :: localSize(nmachs,idim_type_min:idim_type_max)
    INTEGER            :: disp(nmachs,idim_type_min:idim_type_max)
    INTEGER            :: maxLocalSize
    INTEGER            :: sizeGathered(idim_type_min:idim_type_max)
    INTEGER            :: maxSizeGathered
    INTEGER            :: sizeFullField(idim_type_min:idim_type_max)
    INTEGER            :: maxsizeFullField
    INTEGER            :: globalSize(idim_type_min:idim_type_max)
    REAL, ALLOCATABLE  :: localChunk(:)
    REAL, ALLOCATABLE  :: gathered(:)
    REAL, ALLOCATABLE  :: fullField(:)

    ! Recreating Global information about Soil Water
    ! grid dependent, field independent constants for gather and unpacking
    ! as a function of idim_type
    CALL LocalSizesAndDisp(ifm, il1, ir2, jb1, jt2, localSize, disp)
    maxLocalSize = MAXVAL(localSize(mynum,:))
    ALLOCATE(localChunk(maxLocalSize), stat=ierr)
    IF (ierr/=0) THEN
       CALL fatal_error("Error allocating localChunk (gatherData)")
    ENDIF
    CALL CopyLocalChunk(localData4D(1,1,1,1), localChunk, &
         LocalSize(mynum,idim_type))
    sizeGathered(:) = disp(nmachs,:) + localSize(nmachs,:)
    maxSizeGathered = MAXVAL(sizeGathered)
    ALLOCATE(gathered(maxSizeGathered), stat=ierr)
    IF (ierr/=0) THEN
       CALL fatal_error("Error allocating gathered (gatherData)")
    ENDIF
    ! grid dependent field sizes as a function of idim_type
    CALL GlobalSizes(ifm, nmachs, nwave, globalSize)
    IF (mchnum==master_num) THEN
       sizeFullField(:) = globalSize(:)
    ELSE
       sizeFullField(:) = 1
    END IF
    maxSizeFullField = MAXVAL(sizeFullField)
    ALLOCATE(fullField(sizeFullField(idim_type)), stat=ierr)
    IF (ierr/=0) THEN
       CALL fatal_error("Error allocating fullField (gatherData)")
    ENDIF
    localVarn = trim(varn)
    CALL PreProcAndGather(.FALSE., ifm, idim_type, localVarn, &
         il1, ir2, jb1, jt2, localSize, disp,                 &
         localSize(mynum,idim_type), LocalChunk,              &
         sizeGathered(idim_type), gathered,                   &
         sizeFullField(idim_type), fullField                  )
    
    IF (mchnum==master_num) THEN
       globalData4D = RESHAPE(fullField, (/mzg, nnxp, nnyp, npat/))
    ENDIF

    call parf_bcast(globalData4D, &
         int(mzg,i8), int(nnxp,i8), int(nnyp,i8), int(npat,i8), master_num)

    DEALLOCATE(fullField)
    DEALLOCATE(gathered)
    DEALLOCATE(localChunk)
    
  END SUBROUTINE gatherData4D

!--(DMK-CCATT-INI)-----------------------------------------------------
  !temporary function 
 subroutine storeOwnChunk_3D(grid, fullGrid, toStore, nz, nx, ny, fieldName)
   
   use mem_grid, only:  &
       runtype
       
   use node_mod, only:	&
       mchnum, 		&
       mynum, 		&
       master_num, 	&
       nmachs, 		&
       nodemxp, 	&
       nodemyp, 	&
       nodeia, 		&
       nodeiz, 		&
       nodeja, 		&
       nodejz, 		&
       nodeibcon, 	&
       nodei0, 		&
       nodej0     
  
  use mem_grid, only:   &
      nnzp,			&
      nnxp,			&
      nnyp
 
  use ParLib, only: &
      parf_bcast
 
  include "i8.h"

  integer, intent(in) 			:: grid
  integer, intent(in) 			:: nz
  integer, intent(in)			:: nx
  integer, intent(in)			:: ny
  real,  pointer, dimension(:,:,:)	:: toStore
  real,  dimension(nz,nx,ny)		:: fullGrid
  character(len=*), intent(in)		:: fieldName

    integer 		:: ldimx
    integer		:: ldimy
    integer		:: lin
    integer 		:: ia
    integer		:: iz
    integer		:: ja
    integer		:: jz
    integer		:: ka
    integer		:: kz
    integer 		:: ierr
    character(len=8)	:: c0
    character(len=8)	:: c1
    
    character(len=*), parameter :: h="**(storeOwnChunk_3D)**"

    ! check allocated memory

    if (.not. associated(toStore)) then
       call fatal_error(h//" will store at not associated pointer var "//trim(fieldName))
    end if
    
    lin = size(toStore,1)
    if (nz /= lin) then
       write(c0,"(i8)") nz
       write(c1,"(i8)") lin
       call fatal_error(h//trim(fieldName)//" z dimension ("//trim(adjustl(c1))//&
            ") differs from required ("//trim(adjustl(c0))//")")
    end if
    
    ldimx = nodemxp(mynum,grid)
    lin = size(toStore,2)
    if (ldimx /= lin) then
       write(c0,"(i8)") ldimx
       write(c1,"(i8)") lin
       call fatal_error(h//trim(fieldName)//" x dimension ("//trim(adjustl(c1))//&
            ") differs from required ("//trim(adjustl(c0))//")")
    end if
    
    ldimy = nodemyp(mynum,grid)
    lin = size(toStore,3)
    if (ldimy /= lin) then
       write(c0,"(i8)") ldimy
       write(c1,"(i8)") lin
       call fatal_error(h//trim(fieldName)//" y dimension ("//trim(adjustl(c1))//&
            ") differs from required ("//trim(adjustl(c0))//")")
    end if
    
    ia = nodei0(mynum,grid)+1
    iz = nodei0(mynum,grid)+nodemxp(mynum,grid)
    ja = nodej0(mynum,grid)+1
    jz = nodej0(mynum,grid)+nodemyp(mynum,grid)

    call mk_3_buff(fullGrid(1,1,1), toStore(1,1,1), &
                   nz, nnxp, nnyp, nz, ldimx, ldimy, ia, iz, ja, jz)



  end subroutine storeOwnChunk_3D
!--(DMK-CCATT-FIM)-----------------------------------------------------

end module ReadBcst
