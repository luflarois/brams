module modTimeLineFRN
   !! Módulo para tratar das saídas da timeline para o projeto Furnas
   !!
   !! @note
   !!
   !! **Project**: BRAMS-FURNAS
   !! **Author(s)**: Rodrigues, L.F. [LFR]
   !! **e-mail**: <mailto:luiz.rodrigues@inpe.br>
   !! **Date**:  26Abril2022 17:26
   !!
   !! **Full description**:
   !! Módulo para tratar das saídas da timeline para o projeto Furnas.
   !! Este módulo contem rotinas que identificam as variáveis para os pontos selecionados
   !! no arquivo "pontos.csv" e escreve a linha temporal com os valores para os pontos.
   !! São escritas as variáveis definidas em uma lista.
   !!
   !! @endnote
   !!
   !! @warning
   !!
   !!  [](https://www.gnu.org/graphics/gplv3-127x51.png'')
   !!
   !!     This program is free software: you can redistribute it and/or modify
   !!     it under the terms of the GNU General Public License as published by
   !!     the  Free  Software  Foundation, either version 3 of the License, or
   !!     (at your option) any later version.
   !!
   !!     This program is distributed in the hope that it  will be useful, but
   !!     WITHOUT  ANY  WARRANTY;  without  even  the   implied   warranty  of
   !!     MERCHANTABILITY or FITNESS FOR A  PARTICULAR PURPOSE.  See  the, GNU
   !!     GNU General Public License for more details.
   !!
   !!     You should have received a copy  of the GNU General  Public  License
   !!     along with this program.  If not, see <https://www.gnu.org/licenses/>.
   !!
   !! @endwarning

   use dump
   implicit none
   include 'constants.f90'
   character(len=*), parameter :: sourceName = 'modTimeLineFRN.f90' ! Nome do arquivo fonte
   character(len=*), parameter :: moduleName = 'modTimeLIneFRN' ! Nome do módulo

   integer, parameter :: maxSites = 128
   !! Número máximo de sites permitidos
   integer, parameter :: maxVars = 32
   integer, parameter :: maxVarIn = 256

   private
   public :: writeTimeLineFRN, readSites, createSitesFile, isTimeTograds

   interface writeTimeLineFRN
    module procedure writeTimeLineFRN_2D
    module procedure writeTimeLineFRN_3D
   end interface

   type t_sit
    character(len=32) :: nome
    real :: lat 
    real :: lon
    integer :: xpos
    integer :: ypos
    integer :: fileNumber
    integer :: localXpos
    integer :: localYpos
    real, pointer, dimension(:,:,:) :: varValues ! (var,time,z)
   end type
   !! Tipo com informações dos sites
   type(t_sit) :: sites(maxSites)
   !! Descrição dos sites lidos e suas posições
   integer :: nSites
   !! Número total de sites (obtidos do arquivo de entrada)
   integer :: count_var
   !! Numero de variáveis por estação
   character(len=32), allocatable :: estvar(:)

   type vvar
      character(len=32)  :: varName
      character(len=256) :: varDesc
      character(len=16)  :: varUnit
   end type vvar
   type(vvar) :: varin(maxVarin)
   integer :: realVarIn
   real :: timeToOutput 

contains

   subroutine writeTimeLineFRN_3d(fieldName,OutputField,mxp,myp,mzp,time,oneNamelistFile)
      !! Escreve a time line da variável 3D recebida do Pós
      !!
      !! @note
      !!
      !! **Project**: BRAMS-Furnas
      !! **Author(s)**: Rodrigues, L.F. [LFR]
      !! **e-mail**: <mailto:luiz.rodrigues@inpe.br>
      !! **Date**:  26Abril2022 17:31
      !!
      !! **Full description**:
      !! Escreve a time line da variável 3D recebida do Pós
      !!
      !! @endnote
      !!
      !! @warning
      !!
      !!  [](https://www.gnu.org/graphics/gplv3-127x51.png'')
      !!
      !!     Under the terms of the GNU General Public version 3
      !!
      !! @endwarning
   
      USE node_mod, only: &
      ia, iz, ja, jz, mynum

      use ModNamelistFile, only: namelistFile

      implicit none
      !Parameters:
      character(len=*), parameter :: procedureName = 'writeTimeLineFRN_3d' ! Nome da subrotina
   
      !Variables (input, output, inout)
      type(NamelistFile), pointer :: oneNamelistFile
      character(len=*), intent(in) :: fieldName
      !! Nome do campo a ser avaliado
      integer, intent(in) :: mxp,myp,mzp
      real, intent(in) :: OutputField(mxp,myp,mzp)
      real, intent(in) :: time
   
      !Local variables:
      integer :: i,j,k,n,sout
      integer :: vpos
   
      !Code:
      vpos = checkVar(fieldName)
      if( vpos == 0) return

      !Code:
      do n = 1,nsites  
         if(sites(n)%localXpos>0 .and. sites(n)%localYpos>0) then !Existe algum ponto preenchido
            sout = writeVar3D(n,OutputField,mxp,myp,mzp,oneNamelistFile%inplevs,vpos,int(time))
         endif
      enddo


   end subroutine writeTimeLineFRN_3d

   function checkVar(varName) result(varPos)
      !! Verifica se a variável passada está no array de variáveis
      !!
      !! @note
      !!
      !! **Project**: BRAMS_FURNAS
      !! **Author(s)**: Rodrigues, L.F. [LFR]
      !! **e-mail**: <mailto:luiz.rodrigues@inpe.br>
      !! **Date**:  29Abril2022 16:12
      !!
      !! **Full description**:
      !! Verifica se a variável passada está no array de variáveis
      !!
      !! @endnote
      !!
      !! @warning
      !!
      !!  [](https://www.gnu.org/graphics/gplv3-127x51.png'')
      !!
      !!     Under the terms of the GNU General Public version 3
      !!
      !! @endwarning
   
      implicit none
      !Parameters:
      character(len=*), parameter :: procedureName = 'checkVar' ! Nome da função
   
      !Variables (input):
      character(len=*), intent(in) :: varName
      integer :: varPos
   
      !Local variables:
      integer :: i

      varPos = 0
   
      !Code:
      !print *,count_var
      do i = 1,count_var
         !print *,i!,to_lower(trim(varName))
         if(trim(varName)==trim(estvar(i))) then
            varPos = i
            return
         endif
      enddo
   
   end function checkVar

   subroutine writeTimeLineFRN_2d(fieldName,OutputField,mxp,myp,time)
      !! Escreve a time line da variável 2D recebida do Pós
      !!
      !! @note
      !!
      !! **Project**: BRAMS-Furnas
      !! **Author(s)**: Rodrigues, L.F. [LFR]
      !! **e-mail**: <mailto:luiz.rodrigues@inpe.br>
      !! **Date**:  26Abril2022 17:31
      !!
      !! **Full description**:
      !! Escreve a time line da variável 2D recebida do Pós
      !!
      !! @endnote
      !!
      !! @warning
      !!
      !!  [](https://www.gnu.org/graphics/gplv3-127x51.png'')
      !!
      !!     Under the terms of the GNU General Public version 3
      !!
      !! @endwarning
   
      USE node_mod, only: &
      ia, iz, ja, jz, mynum

      implicit none
      !Parameters:
      character(len=*), parameter :: procedureName = 'writeTimeLineFRN_2d' ! Nome da subrotina
   
      !Variables (input, output, inout)
      character(len=*), intent(in) :: fieldName
      !! Nome do campo a ser avaliado
      integer, intent(in) :: mxp,myp
      real, intent(in) :: OutputField(mxp,myp)
      real, intent(in) :: time

      !Local variables:
      integer :: i,j,k,n,sout
      integer :: vpos
   
      !Code:
      vpos = checkVar(fieldName)
      if( vpos == 0) return

      do n = 1,nsites  
         if(sites(n)%localXpos>0 .and. sites(n)%localYpos>0) then !Existe algum ponto preenchido
            sout = writeVar2D(n,OutputField,mxp,myp,vpos,int(time))
         endif
      enddo


   end subroutine writeTimeLineFRN_2d

 !=============================================================================================
 subroutine readSites(mchnum,master_num,oneNamelistFile)
    !! Lê o arquivo de sites de localização de pontos
    !!
    !! @note
    !!
    !! **Project**  : BRAMS/FURNAS
    !! **Author(s)**: Luiz Flávio Rodrigues [lufla]
    !! **e-mail**   : <luiz.rodrigues@inpe.br>
    !! **Date**: 13 December 2021 (Monday)
    !!
    !! **Full Description**: 
    !! Lê o arquivo de sites de localização de pontos
    !!
    !! @endnote
    !!
    !! @warning
    !! ![](https://licensebuttons.net/l/by-sa/3.0/88x31.png"")
    !! Now is under CC Attribution-ShareAlike 4.0 International, please see:
    !! &copy; <https://creativecommons.org/licenses/by-sa/4.0/>
    !! @endwarning
    !!
     
    !Uses 
   use mem_grid, only: &
      oneGlobalGridData, ngrids, nnxp, nnyp, &
      iyear1,        & ! (IN)
      imonth1,       & ! (IN)
      idate1,        & ! (IN)
      itime1,        &  ! (IN)
      timmax,        &
      dtlong,        &
      timeunit

   USE ReadBcst, ONLY: &
      Broadcast

   USE node_mod, only: &
      ia, iz, ja, jz, &
      nodei0, nodej0, &
      mynum

   use ModNamelistFile, only: namelistFile
 
    implicit none
    
    !Parameters
    character(len=*),parameter :: procedureName='**readSites**' !Name of this procedure

    type(NamelistFile), pointer :: oneNamelistFile
    integer, intent(in) :: mchnum
    !! Número do processador atual
    integer, intent(in) :: master_num
    !! Número do processador mestre

    !Local Variables
    integer :: fileNum
    character(len=256) :: linha
    character(len=256) :: campos4(4)
    character(len=256) :: campos3(3)
    integer :: count, rsp, i, j, n
    real :: val, dist, d1, d2
    character(len=32) :: nome
    real :: lat(1)
    real :: lon(1)
    integer :: xpos
    integer :: ypos
    integer :: fileNumber
    integer :: localXpos
    integer :: localYpos
    character(len=8) :: varname
    integer :: totsec
 

    !Code area
    if(mchnum == master_num) then

      if(.not. fileExist("estvars.dat")) iErrNumber = dumpMessage(c_tty,c_yes,'','',c_fatal &
          ,'Arquivo de variáveis, estvars.dat, não encontrado, verifique!')

      open(unit = 87, file = "estvars.dat", status = 'old', action = 'read')
      count_var = 0
      do
         count_var = count_var+1
         read(87,fmt='(A)', END=30) linha
      enddo

30    count_var = count_var-1
      allocate(estvar(count_var))
      
      rewind(87)

      print *,'=== variables to stations -  CSV File ==='
      do i=1, count_var
         read(87,fmt='(A)', END=35) estvar(i)
         print *,i,':',estvar(i)
      enddo
      print *,'========================================='
35    close(unit=87)


      if(.not. fileExist("estacoes.csv")) iErrNumber = dumpMessage(c_tty,c_yes,'','',c_fatal &
      ,'Arquivo de estações, estacoes.csv, não encontrado, verifique!')

      print *,'===============  Stations ==============='
       open(unit = 87, file = "estacoes.csv", status = 'old', action = 'read')
       !Pula a linha de cabeçalho
       read(87,*)
       count = 0
       print *,'Limite Lon min,max :', oneGlobalGridData(1)%global_glon(1,1),',' &
              ,oneGlobalGridData(1)%global_glon(nnxp(1),1)
       print *,'Limite Lat min,max :', oneGlobalGridData(1)%global_glat(1,1),',' &
               ,oneGlobalGridData(1)%global_glat(1,nnyp(1))

       do 
          count = count+1
          read(87,fmt='(A)',END=50) linha
          !Faz um parse separando os campos
          campos3 = parseCsv(linha,3,';')
          sites(count)%nome = trim(campos3(1))
          read(campos3(2),*) sites(count)%lat
          read(campos3(3),*) sites(count)%lon
          !
          sites(count)%xpos = 0
          sites(count)%ypos = 0

          !Checa se o site está dentro do domínio do modelo
          print *,'Nome da estacao:',trim(sites(count)%nome)
          print *,'Lon do site    :', sites(count)%lon
          print *,'Lat do site    :', sites(count)%lat

          if(sites(count)%lon<=oneGlobalGridData(1)%global_glon(1,1) &
            .or. sites(count)%lon>=oneGlobalGridData(1)%global_glon(nnxp(1),1)) &
            iErrNumber = dumpMessage(c_tty,c_yes,'','',c_fatal &
          ,'Longitude do ponto fora do domínio do Modelo!',sites(count)%lon,"F6.2")

          if(sites(count)%lat<=oneGlobalGridData(1)%global_glat(1,1) &
          .or. sites(count)%lat>=oneGlobalGridData(1)%global_glat(1,nnyp(1))) &
          iErrNumber = dumpMessage(c_tty,c_yes,'','',c_fatal &
          ,'Latitude do ponto fora do domínio do Modelo!',sites(count)%lat,"F6.2")

          !Procura o ponto em longitude do modelo mais próximo a lon lida 
          !guarda esse ponto em xpos
          dist = 1000.0
          do i = 2,nnxp(1)
             d1 = abs(sites(count)%lon-oneGlobalGridData(1)%global_glon(i-1,1))
             d2 = abs(sites(count)%lon-oneGlobalGridData(1)%global_glon(i,1))
             if(d1 < dist .and. d1 < d2) then 
               sites(count)%xpos = i-1
               dist = d1
             elseif(d2 < dist .and. d2 < d1) then 
               sites(count)%xpos = i-1
               dist = d2
             endif
          enddo
          if(sites(count)%xpos == 0) iErrNumber = dumpMessage(c_tty,c_yes,'','',c_fatal &
          ,'Longitude do ponto fora do domínio do Modelo!',sites(count)%lon,"F6.2")    
          !Procura o ponto em latitude do modelo mais próximo a lat lida 
          !guarda esse ponto em ypos
          dist = 1000.0
          do i = 2,nnyp(1)
             d1 = abs(sites(count)%lat-oneGlobalGridData(1)%global_glat(1,i-1))
             d2 = abs(sites(count)%lat-oneGlobalGridData(1)%global_glat(1,i))
             if(d1 < dist .and. d1 < d2) then 
               sites(count)%ypos = i-1
               dist = d1
             elseif(d2 < dist .and. d2 < d1) then 
               sites(count)%ypos = i-1
               dist = d2
             endif
          enddo
          if(sites(count)%ypos == 0) iErrNumber = dumpMessage(c_tty,c_yes,'','',c_fatal &
          ,'Latitude do ponto fora do domínio do Modelo!',sites(count)%lat,"F6.2")     
       enddo
50     nSites = count-1

      close(unit = 87)
    endif
    
    CALL Broadcast(count_var, master_num, "count_var")
    if(mchnum /= master_num) allocate(estvar(count_var))
    do i=1,count_var
       varname = estvar(i)
       CALL Broadcast(varname, master_num, "varname")
       estvar(i) = varname
    end do

    CALL Broadcast(nsites, master_num, "nsites")
    do n=1,nsites
      nome = sites(n)%nome
      lat  = sites(n)%lat
      lon  = sites(n)%lon
      xpos = sites(n)%xpos
      ypos = sites(n)%ypos
      fileNumber =sites(n)%fileNumber
      CALL Broadcast(nome, master_num, "nome")
      CALL Broadcast(lat, master_num, "lat")
      !print *,'%%%%:',mchnum,lat
      CALL Broadcast(lon, master_num, "lon")
      CALL Broadcast(xpos, master_num, "xpos")
      CALL Broadcast(ypos, master_num, "ypos")
      CALL Broadcast(fileNumber, master_num, "fileNumber")
      sites(n)%nome = nome
      sites(n)%lat  = lat(1)
      sites(n)%lon  = lon(1)
      !print *,'Lat=',mchnum,lat,sites(n)%lat; call flush(6)
      !print *,'Nom=',mchnum,sites(n)%nome; call flush(6)
      sites(n)%xpos = xpos
      !print *,'Xpos=',mchnum,sites(n)%xpos; call flush(6)
      sites(n)%ypos = ypos
      sites(n)%fileNumber = fileNumber

      !Procura pelos pontos no processador local
      do i=ia,iz
         if(nodei0(mynum,1)+i == sites(n)%xpos) then 
            do j=ja,jz
               if(nodej0(mynum,1)+j == sites(n)%ypos) then 
                  !Encontrou os pontos xpos e ypos nesse processador
                  !guardar no site
                  sites(n)%localXpos = i
                  sites(n)%localypos = j
                  print *,'Localizacao    pos X:',n, sites(n)%xpos,'=',sites(n)%localxpos, ' no proc ',mynum
                  print *,'Localizacao    pos Y:',n, sites(n)%ypos,'=',sites(n)%localypos, ' no proc ',mynum
                  print *,'---------------------------------------------------------------'                  
                  !rsp = createSitesFile(n,oneNamelistFile%inplevs)
               endif
            enddo
         endif
      enddo 

      if(sites(n)%localXpos>0 .and. sites(n)%localypos>0) then
         allocate(sites(n)%varValues(count_var,0:int(timmax),oneNamelistFile%inplevs))
      endif

    enddo

 end subroutine readSites 

 function writeVar3D(siteNumber,OutputField,mxp,myp,mzp,levels,vpos,time) result(sout)
    !! Escreve um campo 3D no arquivo correspondente
    !!
    !! @note
    !!
    !! **Project**: BRAMS_FURNAS
    !! **Author(s)**: Rodrigues, L.F. [LFR]
    !! **e-mail**: <mailto:luiz.rodrigues@inpe.br>
    !! **Date**:  28Abril2022 19:34
    !!
    !! **Full description**:
    !! Escreve um campo 3D no arquivo correspondente
    !!
    !! @endnote
    !!
    !! @warning
    !!
    !!  [](https://www.gnu.org/graphics/gplv3-127x51.png'')
    !!
    !!     Under the terms of the GNU General Public version 3
    !!
    !! @endwarning
 
    implicit none
    !Parameters:
    character(len=*), parameter :: procedureName = 'writeVar3D' ! Nome da função

    !Variables (input):
    integer, intent(in) :: siteNumber
    integer, intent(in) :: mxp,myp,mzp,levels
    real, intent(in) :: OutputField( mxp,myp,mzp)
    integer, intent(in) :: vpos
    integer, intent(in) :: time

   
    !Local variables:
    integer :: sout,lev
    character(len=2) :: clev
 
    !Code:
    do lev = 1,levels
      sites(siteNumber)%varValues(vPos,time,lev) = OutputField(sites(siteNumber)%localXpos,sites(siteNumber)%localYpos,lev)
    enddo

    sout=1
 
 end function writeVar3D

 function writeVar2D(siteNumber,OutputField,mxp,myp,vpos,time) result(sout)
   !! Escreve um campo 2D no arquivo correspondente
   !!
   !! @note
   !!
   !! **Project**: BRAMS_FURNAS
   !! **Author(s)**: Rodrigues, L.F. [LFR]
   !! **e-mail**: <mailto:luiz.rodrigues@inpe.br>
   !! **Date**:  28Abril2022 19:34
   !!
   !! **Full description**:
   !! Escreve um campo 2D no arquivo correspondente
   !!
   !! @endnote
   !!
   !! @warning
   !!
   !!  [](https://www.gnu.org/graphics/gplv3-127x51.png'')
   !!
   !!     Under the terms of the GNU General Public version 3
   !!
   !! @endwarning

   implicit none
   !Parameters:
   character(len=*), parameter :: procedureName = 'writeVar2D' ! Nome da função

   !Variables (input):
   integer, intent(in) :: siteNumber
   integer, intent(in) :: mxp,myp
   real, intent(in) :: OutputField(mxp,myp)
   integer, intent(in) :: time

   integer, intent(in) :: vpos

   integer :: sout

   !Local variables:

   !Code:
   sites(siteNumber)%varValues(vPos,time,1) = OutputField(sites(siteNumber)%localXpos,sites(siteNumber)%localYpos)

   sout=1

end function writeVar2D


 function createSitesFile(oneNamelistFile) result(isok)
    !! Cria os arquivos csv para saídas das timelines
    !!
    !! @note
    !!
    !! **Project**: BRAMS_Furnas
    !! **Author(s)**: Rodrigues, L.F. [LFR]
    !! **e-mail**: <mailto:luiz.rodrigues@inpe.br>
    !! **Date**:  26Abril2022 17:53
    !!
    !! **Full description**:
    !! Cria os arquivos csv para saídas das timelines
    !!
    !! @endnote
    !!
    !! @warning
    !!
    !!  [](https://www.gnu.org/graphics/gplv3-127x51.png'')
    !!
    !!     Under the terms of the GNU General Public version 3
    !!
    !! @endwarning
   use mem_grid, ONLY: &
      zm, grid_g

   use mem_grid, ONLY: &
      iyear1,        & ! (IN)
      imonth1,       & ! (IN)
      idate1,        & ! (IN)
      itime1,        &  ! (IN)
      timmax,        &
      dtlong,        &
      timeunit

   use io_params         , only: frqanl

   use ModNamelistFile, only: namelistFile

   use ModDateUtils, only: &
      date_add_to_dble

   USE node_mod, only: &
      ia, iz, ja, jz, &
      nodei0, nodej0, &
      mynum

    implicit none
    !Parameters:
    character(len=*), parameter :: procedureName = 'createSitesFile' ! Nome da função
    character(len=*), parameter :: prefixFormat="(A19,','"
    character(len=*), parameter :: sufixFormat ="(F14.4,','),F14.4)"

    character(len=*), parameter :: prefixFormat2="(I4.4,'-',I2.2,'-',I2.2,'T',I2.2,':',I2.2,':',I2.2,','"
    character(len=*), parameter :: sufixFormat2 ="(F14.4,','),F14.4)"

    integer, parameter :: fileNumber = 90
 
    !Variables (input):
    !integer, intent(in) :: nfile
    !integer, intent(in) :: nlevels
     integer :: isok
     type(NamelistFile), pointer :: oneNamelistFile
 
    !Local variables:
    integer :: fn
    character(len=2) :: clev
    integer :: sout
    integer :: l,n,v,seconds,totsec,secc,i
    character(len=10) :: cdate !"YYYYMMDDHH" 
    integer :: iyy,imm,idd,ihh,hh,mm,ss
    integer :: count
    character(len=256) :: linha
    character(len=256) :: campos4(4)
    character(len=16) :: vardim
    character(len=256) :: vardesc

    !Code:
    if(.not. fileExist("variables.csv")) iErrNumber = dumpMessage(c_tty,c_yes,'','',c_fatal &
    ,'Arquivo de variáveis, variables.csv, não encontrado, verifique!')      

    !Lê as variaveis do arquivo CSV para obter unidade e descrição
    open(unit = 87, file = "variables.csv", status = 'old', action = 'read')
    count = 0
    !print *,'=== Campos de variables.csv ==='      
    do 
       count = count+1
       read(87,fmt='(A)',END=40) linha
       !Faz um parse separando os campos
       campos4 = parseCsv(linha,4,',')
       varin(count)%varName = campos4(2)
       varin(count)%varDesc = campos4(3)
       varin(count)%varUnit = campos4(4)
!       print *,count,':',varin(count)%varName,':',varin(count)%varUnit
    enddo
40    close(unit=87)
    realVarIn = count-1


    write(clev,fmt="(I2.2)") oneNamelistFile%inplevs
    write(cdate,fmt='(I4,I2.2,I2.2,I2.2)') iyear1,imonth1,idate1,itime1

    !do n=1,nsites
    !  print *,mynum,n,sites(n)%nome
    !  print *,mynum,sites(n)%lon,sites(n)%lat
    !  print *,mynum,sites(n)%localXpos,sites(n)%localYpos
    !enddo

    do n=1,nsites
      if(sites(n)%localXpos<=0 .or. sites(n)%localYpos<=0) cycle
      print *,'*** Escrevendo para site ',n,trim(sites(n)%nome),sites(n)%lon,sites(n)%lat
      do v=1,count_var
         seconds = 0
         open(unit = fileNumber, file = trim(sites(n)%nome)//'_'//trim(estvar(v))//'_'//cdate//'.csv' &
              ,status = "replace", action = "write")
         !Obtém a unidade da variável a ser escrita
         vardim = getUnitInVarList(estvar(v))
         vardesc = getDescInVarList(estvar(v))

         !Escreve um cabeçalho com informações no arquivo CSV
         write(fileNumber,fmt='(A,A,A,F9.5,A,F9.5,A,A)') '# Local: ',trim(sites(n)%nome) &
               ,', Lat=',sites(n)%lat,', Lon= ',sites(n)%lon
         write(fileNumber,fmt='(A,A,A,A,A,A)') '# Sinal: ',trim(estvar(v)),' ,Unidade: ',trim(vardim) &
            ,' ,Desc: ',trim(vardesc)
         
         write(fileNumber,fmt=prefixFormat//clev//sufixFormat) 'Tempo' &
               ,(zm(l)*grid_g(1)%rtgt(sites(n)%localxpos,sites(n)%localypos),l=1,oneNamelistFile%inplevs)
         
         !escreve os tempos e os valores da variável para cada nível
         do seconds = 0,int(timmax),int(frqanl)
            secc = seconds !-int(frqanl)+1

            call date_add_to_dble(iyear1,imonth1,idate1,itime1,dble(secc),'s' &
            ,iyy,imm,idd,ihh)
            hh=int(ihh/10000)
            mm=int((ihh-hh*10000)/100)
            ss=int(ihh-hh*10000-mm*100)            

            write(fileNumber,fmt=prefixFormat2//clev//sufixFormat2) iyy,imm,idd,hh,mm,ss &
                  ,sites(n)%varValues(v,int(secc),1:oneNamelistFile%inplevs)
         enddo

         close(unit = fileNumber)
      enddo
   enddo
   isok = 0

 end function createSitesFile

 function closeSitesFile() result(stt)
    !! Close the sites files
    !!
    !! @note
    !!
    !! **Project**: BRAMS-FURNAS
    !! **Author(s)**: Rodrigues, L.F. [LFR]
    !! **e-mail**: <mailto:luiz.rodrigues@inpe.br>
    !! **Date**:  09Maio2022 12:51
    !!
    !! **Full description**:
    !! Close the sites files
    !!
    !! @endnote
    !!
    !! @warning
    !!
    !!  [](https://www.gnu.org/graphics/gplv3-127x51.png'')
    !!
    !!     Under the terms of the GNU General Public version 3
    !!
    !! @endwarning
 
    implicit none
    !Parameters:
    character(len=*), parameter :: procedureName = 'closeSitesFile' ! Nome da função
 
    !Variables (input):
    integer :: stt
 
    !Local variables:
 
    !Code:
    !close(unit = sites(nfile)%fileNumber)

    stt = 0
 
 end function closeSitesFile

 !=============================================================================================
 function parseCsv(linha,nCampos,separador) result(campos)
    !! Parse uma linha de csv
    !!
    !! @note
    !!
    !! **Project**  : BRAMS/FURNAS
    !! **Author(s)**: Luiz Flávio Rodrigues [lufla]
    !! **e-mail**   : <luiz.rodrigues@inpe.br>
    !! **Date**: 13 December 2021 (Monday)
    !!
    !! **Full Description**: 
    !! Parse uma linha de csv
    !!
    !! @endnote
    !!
    !! @warning
    !! ![](https://licensebuttons.net/l/by-sa/3.0/88x31.png"")
    !! Now is under CC Attribution-ShareAlike 4.0 International, please see:
    !! &copy; <https://creativecommons.org/licenses/by-sa/4.0/>
    !! @endwarning
    !!
     
    !Uses
 
    implicit none
    
    !Parameters
    character(len=*),parameter :: procedureName='**parseCsv**' !Name of this procedure
 
    !Variables intent(in/out)
    character(len=*), intent(in) :: linha
    !! linha com os campos separados pelo separador
    integer, intent(in) :: nCAmpos
    !! Total de campos
    character, intent(in) :: separador
    !! caracter que separa os campos

    !Local Variables
    character(len=256) :: campos(nCampos), resto
    integer :: beg,ipos,i,sizeLin
 
    !Code area
    resto = trim(linha)
    beg = 1
    do i = 1,nCampos
         sizeLin=len(resto)
         !print *,'resto Ini:',resto,sizelin
         ipos = index(resto,separador)
         !print *,ipos
         campos(i)=resto(beg:ipos-1)
         resto=trim(resto(ipos+1:sizeLin))
         !print *,'resto Fin:',resto
         if(i == nCampos) then
             campos(i)=resto
         endif
     enddo
    
 end function parseCsv 

 !=============================================================================================
 logical function fileExist(fileName)
 !! Check if fileName exist
 !!
 !! @note
 !! ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
 !!
 !! **Brief**: Check if fileName exist
 !!
 !! **Documentation**: <http://brams.cptec.inpe.br/documentation/>
 !!
 !! **Author(s)**: Luiz Flavio Rodrigues **&#9993;**<luiz.rodrigues@inpe.br>
 !!
 !! **Date**: 26 August 2020 (Wednesday)
 !! @endnote
 !!
 !! @changes
 !! &#9744; <br/>
 !! @endchanges
 !! @bug
 !!
 !!@endbug
 !!
 !!@todo
 !!  &#9744; <br/>
 !! @endtodo
 !!
 !! @warning
 !! Now is under CC-GPL License, please see
 !! &copy; <https://creativecommons.org/licenses/GPL/2.0/legalcode.pt>
 !! @endwarning
 !!
 
 !Use area
 use dump

 implicit none

 include "constants.f90"
 character(len=*),parameter :: procedureName='**fileExist**' !Name of this procedure
 !
 !Local Parameters

 !Input/Output variables
 character(len=*), intent(in) :: fileName

 !Local variables

 !Code

 inquire(file=fileName(1:len_trim(fileName)),exist=fileExist)


end function fileExist 

    !=============================================================================================
function to_lower(strIn) result(strOut)
   !! Convert case to lower case
   !!
   !! @note
   !! ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
   !!
   !! **Brief**: Convert case to lower case
   !!
   !! **Documentation**: <http://brams.cptec.inpe.br/documentation/>
   !!
   !! **Author(s)**: Luiz Flavio Rodrigues **&#9993;**<luiz.rodrigues@inpe.br>
   !!
   !! **Date**: 26 August 2020 (Wednesday)
   !! @endnote
   !!
   !! @changes
   !! &#9744; <br/>
   !! @endchanges
   !! @bug
   !!
   !!@endbug
   !!
   !!@todo
   !!  &#9744; <br/>
   !! @endtodo
   !!
   !! @warning
   !! Now is under CC-GPL License, please see
   !! &copy; <https://creativecommons.org/licenses/GPL/2.0/legalcode.pt>
   !! @endwarning
   !!
   
   !Use area
   use dump

   implicit none

   include "constants.f90"
   character(len=*),parameter :: procedureName='**to_lower**' !Name of this procedure
   !
   !Local Parameters

   !Input/Output variables
   character(*), intent(in) :: strIn
   !! String to be converted
   character(len=len(strIn)) :: strOut
   !! Return string converted

   !Local variables
   integer :: i
   !Code
   do i = 1, len(strIn)
       select case(strIn(i:i))
       case("A":"Z")
          strOut(i:i) = achar(iachar(strIn(i:i))+32)
       end select
    end do

end function to_lower

function getUnitInVarList(varName) result(varU)
   !! retorna a unidade da variável verificada no variables.csv
   !!
   !! @note
   !!
   !! **Project**: BRAMS-Furnras
   !! **Author(s)**: Rodrigues, L.F. [LFR]
   !! **e-mail**: <mailto:luiz.rodrigues@inpe.br>
   !! **Date**:  15Julho2022 13:46
   !!
   !! **Full description**:
   !! retorna a unidade da variável verificada no variables.csv
   !!
   !! @endnote
   !!
   !! @warning
   !!
   !!  [](https://www.gnu.org/graphics/gplv3-127x51.png'')
   !!
   !!     Under the terms of the GNU General Public version 3
   !!
   !! @endwarning
   use ModPostTypes, only: &
      all_post_variables

   implicit none
   !Parameters:
   character(len=*), parameter :: procedureName = 'getUnitInVarList' ! Nome da função

   !Variables (input):
   character(len=32), intent(in) :: varName
   character(len=16) :: varU

   !Local variables:
   integer :: i
   character(len=32) :: vnam

   !Code:
   varU = ''
   do i = 1, realVarIn
      vnam = varin(i)%varName(3:len_trim(varin(i)%varName)-1)
      !print *,'Checando: #',trim(varName),'#',trim(vnam),':',trim(varName) == trim(vnam)
      if(trim(varName) == trim(vnam)) then
         varU = varin(i)%varUnit
         exit
      endif
   enddo


end function getUnitInVarList

function getDescInVarList(varName) result(varU)
   !! retorna a descricao da variável verificada no variables.csv
   !!
   !! @note
   !!
   !! **Project**: BRAMS-Furnras
   !! **Author(s)**: Rodrigues, L.F. [LFR]
   !! **e-mail**: <mailto:luiz.rodrigues@inpe.br>
   !! **Date**:  15Julho2022 13:46
   !!
   !! **Full description**:
   !! retorna a descricao da variável verificada no variables.csv
   !!
   !! @endnote
   !!
   !! @warning
   !!
   !!  [](https://www.gnu.org/graphics/gplv3-127x51.png'')
   !!
   !!     Under the terms of the GNU General Public version 3
   !!
   !! @endwarning
   use ModPostTypes, only: &
      all_post_variables

   implicit none
   !Parameters:
   character(len=*), parameter :: procedureName = 'getUnitInVarList' ! Nome da função

   !Variables (input):
   character(len=32), intent(in) :: varName
   character(len=256) :: varU

   !Local variables:
   integer :: i
   character(len=32) :: vnam

   !Code:
   varU = ''
   do i = 1, realVarIn
      vnam = varin(i)%varName(3:len_trim(varin(i)%varName)-1)
      !print *,'Checando: #',trim(varName),'#',trim(vnam),':',trim(varName) == trim(vnam)
      if(trim(varName) == trim(vnam)) then
         varU = varin(i)%varDesc
         exit
      endif
   enddo


end function getDescInVarList

function isTimeTograds() result(istime)
   !! Retorna se é tempo de escrever arquivo grads qdo IPOS=10
   !!
   !! @note
   !!
   !! **Project**: BRAMS-Furnas
   !! **Author(s)**: Rodrigues, L.F. [LFR]
   !! **e-mail**: <mailto:luiz.rodrigues@inpe.br>
   !! **Date**:  19Julho2022 08:27
   !!
   !! **Full description**:
   !! Retorna se é tempo de escrever arquivo grads qdo IPOS=10
   !!
   !! @endnote
   !!
   !! @warning
   !!
   !!  [](https://www.gnu.org/graphics/gplv3-127x51.png'')
   !!
   !!     Under the terms of the GNU General Public version 3
   !!
   !! @endwarning
   use mem_grid, only: &
      time, &
      dtlongn, &
      timmax

   use meteogram, only: &
      meteogramFreq

   use io_params, only : & ! 
      ipos

   implicit none
   !Parameters:
   character(len=*), parameter :: procedureName = 'isTimeTograds' ! Nome da função

   !Variables (input):
   logical :: istime

   !Local variables:

   !Code:
   istime = .false.
   if (IPOS/=10) then 
      return
   elseif (IPOS==10 .and.  (mod(time,meteogramFreq)<dtlongn(1)    .or.  &
      time>=timmax - 0.01*dtlongn(1))) then
      istime = .true.
   endif

end function isTimeTograds

end module modTimeLineFRN
