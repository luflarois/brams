!############################# Change Log ##################################
! Interface entre BRAMS e JULESsfclyr_jules.f90sfclyr_jules.f90sfclyr_jules.f90
! 

SUBROUTINE sfclyr_jules(mzp,mxp,myp,iaI,izI,jaI,jzI,jdim,julesFile,i0,j0)

  !--- Modulos do BRAMS ---

   use io_params, only : frqanl,afilout,frqhis,hfilin

   use node_mod, only: MYNUM
   
   USE leaf_coms, ONLY : slcpd,slmsts,gzotheta 
 
   USE mem_leaf      
 
   USE mem_jules      
 
   USE mem_basic,     ONLY: basic_g
 
   USE mem_grid ,     ONLY : npatch, nzg, grid_g, time, dtlong, iyear1,imonth1, idate1 &
                             ,istp,itime1, timmax,zt,ngrid,runtype,nnxp,nnyp                   
   USE rconstants,    ONLY : cpi,cp,alvl,vonk,g,p00,cpor
 
   USE mem_radiate,   ONLY : radiate_g 
 
   USE mem_turb,      ONLY : turb_g
   
   USE mem_cuparm,    ONLY: cuparm_g, nnqparm
 
   USE micphys,       ONLY: level
 
   USE mem_micro,     ONLY: micro_g
  
   USE chem1_list,    ONLY: CO2,chemical_mechanism
  
   USE mem_chem1,     ONLY: chem1_g, chemistry


  !--- Modulos do JULES ---
   USE io_constants, ONLY: max_file_name_len
  
   USE mem_brams_jules, ONLY: mynumB,glatB,glonB,nxB,nyB,land_fracB,precipB,swdownB &
                       ,lwdownB,diff_radB,tempB,upsB,vpsB,pstarB,qB,fracB       &
                       ,timestepB,main_run_startB,main_run_endB,ntimestepB,sm_levelsB &
                       ,dzsoilB,sthuB,tsoilB,tstarB,viesALL,output_periodB,dir_run_idB &
                       ,dump_periodB,runtypeB,hfilinB,z1_uvB,z1_tqB

   USE gridbox_mean_mod, ONLY: surftiles_to_gbm

   USE jules_surface_types_mod,       ONLY : npft,ntype

   USE csigma, ONLY :  sbcon
  
   USE model_time_mod, ONLY : end_of_run
  
   USE ancil_info, ONLY : land_pts,row_length
                            
   USE jules_fields_mod, ONLY: progs,psparms,ainfo,trifctltype
   
   !USE p_s_parms, ONLY : sthu_soilt,sthf_soilt

   !USE prognostics, ONLY :  tstar_surft, lai_pft

   USE sf_diags_mod, ONLY: sf_diag

   USE datetime_mod, ONLY: datetime,datetime_to_string,datetime_advance,datetime_diff

   !USE trifctl, ONLY :  NPP_gb,   &                   ! Net primary productivity (kg C/m2/s)
   !                     RESP_S_soilt, &                  ! Soil respiration (kg C/m2/s)
   !                     RESP_P_gb, &                  ! Plant respiration (kg C/m2/s)
   !                     GPP_gb                        ! Gross primary productivity (kg C/m2/s)

   USE gridmean_fluxes, ONLY: fqw_1_ij,         &  ! latent heat flux 
                              ftl_1_ij,         &  ! fluxo de calor sensivel do JULES
                              taux_1_ij,        &  !   W'ly component of surface wind stress (N/sq m)
                              tauy_1_ij            !   S'ly component of surface wind stress (N/sq m)

   USE fluxes, ONLY : emis_surft, land_albedo_ij
  
   !USE aero,          ONLY : co2_3d_ij  

   use digitalFilter, only:
      pref_vies
   

   IMPLICIT NONE

   character(len=*), intent(in) :: julesFile
  
   INTEGER               :: nsoil, fase,ia,iz,ja,jz,hh,mm,a,nt_vies_aux
   INTEGER, PARAMETER :: fat_dtlong=1  ! > 1 para nao executar o JULES em todos os timestep do BRAMS
  
   INTEGER, INTENT(IN) :: mzp,mxp,myp,jdim,iaI,izI,jaI,jzI,i0,j0
  
   INTEGER            :: ng,i ,j,v,l,n,ip,k,ntimestep, itime1_seg,maq,file_size
   INTEGER,SAVE :: nt_vies

   LOGICAL              :: there, start=.TRUE.

   REAL               :: zoverl,cx,wtol,psin,piv,dec &
                        ,tempk,fracliq,zts2,ths2,ustar2,tstar2,rstar2

   REAL, DIMENSION(:), ALLOCATABLE :: rlongupJ
    
   CHARACTER (LEN=2)  :: tB_str    

   REAL, DIMENSION(mxp,myp) :: pcpgl

   LOGICAL, SAVE , ALLOCATABLE :: run_jules(:)
  
   REAL, ALLOCATABLE :: dens2(:,:),viesG(:,:),viesM(:,:,:),aux1(:,:),aux2(:,:)
  
   TYPE(datetime) :: dti,dtf  ! Placeholder for a datetime in a calculation

   DATA wtol/1e-20/

   CHARACTER(LEN=max_file_name_len) :: nml_dir  ! Directory containing namelists

   CHARACTER(LEN=256) :: aux,file_vies!,pref_vies
   
   CHARACTER(LEN=4), PARAMETER :: vars(5)=(/'UVEL','VVEL','TEMP','UMID','PNMM'/)
   
   LOGICAL :: rem_vies


!Levar essa variavel para o RAMSIMN
!pref_vies='NONE' ! prefix to vies files - 'NONE' or '' to no remove vies      
!pref_vies='/home/demerval/intercomparacao_dk/Oper_DSM/vies_emq/vies_ERA5+OBS_' ! 'NONE' to no remove vies  !Levar essa variavel para o RAMSIMN

IF (LEN(TRIM(pref_vies)) .lt. 5) THEN  
   rem_vies=.false.
ELSE
   rem_vies=.true.
ENDIF

nml_dir=trim(julesFile)

INQUIRE(FILE=trim(nml_dir)//'/drive.nml',EXIST=there)

IF (.not. there) THEN
   PRINT*;PRINT*;PRINT*;PRINT*, 'Not found:  '//trim(nml_dir)//'/drive.nml'
   PRINT*;PRINT*, 'Check JULESIN variable in RAMSIN'
   STOP 'STOP in: sfclyr_jules.f90 (a)'
ENDIF

ia=iaI-1
iz=izI+1
ja=jaI-1
jz=jzI+1

mynumB=MYNUM
ng=ngrid
nxB=iz-ia+1
nyB=jz-ja+1

IF (.not. ALLOCATED(swdownB)) THEN
   ALLOCATE( swdownB(nxB,nyB),lwdownB(nxB,nyB),diff_radB(nxB,nyB),precipB(nxB,nyB), &
             tempB(nxB,nyB),upsB(nxB,nyB),vpsB(nxB,nyB),pstarB(nxB,nyB),qB(nxB,nyB), &
             tstarB(nxB,nyB) )
ENDIF

!--- Precipitacao total ---
pcpgl(:,:)=0.
IF (nnqparm(ng) > 0 .and. level >= 3) THEN
   pcpgl(:,:)=cuparm_g(ng)%conprr(:,:) + micro_g(ng)%pcpg(:,:)
ELSEIF(nnqparm(ng) == 0 .and. level >= 3) THEN
   pcpgl(:,:)=micro_g(ng)%pcpg(:,:)
ENDIF


!--- Remocao de VIES ---{
IF (start .and. rem_vies) THEN
   DO v=1,5  ! 1=UVEL 2=VVEL  3=TEMP  4=UMID  5=PNMM
      file_vies=trim(pref_vies)//vars(v)//'.gra'
      INQUIRE(FILE=trim(file_vies),EXIST=there)
      IF (there) THEN
         INQUIRE(FILE=file_vies, SIZE=file_size)
         nt_vies_aux=file_size/nnxp(ng)/nnyp(ng)/4
         IF ( file_size/1./nnxp(ng)/nnyp(ng) .ne. real(file_size/nnxp(ng)/nnyp(ng)) ) THEN
            PRINT*, 'ERROR: Vies file not compatible: '//trim(file_vies) 
            PRINT*, nnxp(ng),nnyp(ng),file_size/1./nnxp(ng)/nnyp(ng),real(file_size/nnxp(ng)/nnyp(ng))
            STOP 'STOP in: sfclyr_jules.f90 (b)'
         ENDIF

         IF (.not. ALLOCATED(viesALL)) THEN
            nt_vies=nt_vies_aux
            ALLOCATE(viesALL(nxB,nyB,nt_vies,5))
            viesALL(:,:,:,:)=0.
         ELSE
            IF (nt_vies_aux /= nt_vies) THEN
               PRINT*, 'ERROR: Arquivos de vies com numero de tempos distintos',nt_vies_aux,nt_vies
               STOP 'STOP in: sfclyr_jules.f90 (c)'
            ENDIF
         ENDIF

         ALLOCATE(viesG(nnxp(ng),nnyp(ng)))
         INQUIRE(IOLENGTH=maq) piv
         
         IF (mynum==1) PRINT*, 'Reading: '//TRIM(file_vies)
         OPEN(77,FILE=TRIM(file_vies),FORM='unformatted', ACCESS='direct',RECL=nnxp(ng)*nnyp(ng)*maq)
         DO i=1,nt_vies
            READ(77,REC=i) viesG
            viesALL(:,:,i,v)=viesG(i0+1:i0+nxB,j0+1:j0+nyB)
         ENDDO  

         DEALLOCATE (viesG)
         CLOSE (77)
      ELSE
         PRINT*, 'Not found: '//trim(file_vies)
         STOP 'STOP in: sfclyr_jules.f90 (d)'
      ENDIF
   ENDDO
ENDIF

IF (.not. ALLOCATED(viesM)) THEN
   ALLOCATE(viesM(nxB,nyB,5))  ! 1=UVEL 2=VVEL  3=TEMP  4=UMID  5=PNMM
   viesM(:,:,:)=0.
ENDIF

IF (rem_vies) THEN
   IF (.not. ALLOCATED(viesALL)) THEN
      PRINT*, 'pref_vies not is NONE, but not found any vies file with name: '//trim(pref_vies)//'????.gra'
      STOP 'STOP in: sfclyr_jules.f90 (e)'
   ENDIF
 
   ALLOCATE(aux1(nxB,nyB),aux2(nxB,nyB))

   a=int(time/3600.+1.0)
   IF (a+1>nt_vies) a=nt_vies-1
   dec=time-(a-1)*3600.
   !--- Vento ---
   viesM(:,:,1)=(dec*viesALL(:,:,a+1,1) - dec*viesALL(:,:,a,1) + 3600.*viesALL(:,:,a,1))/3600.
   viesM(:,:,2)=(dec*viesALL(:,:,a+1,2) - dec*viesALL(:,:,a,2) + 3600.*viesALL(:,:,a,2))/3600.

   !--- Temperarura ---
   viesM(:,:,3)=(dec*viesALL(:,:,a+1,3) - dec*viesALL(:,:,a,3) + 3600.*viesALL(:,:,a,3))/3600.

   !--- Umidade ---
   viesM(:,:,4)=(dec*viesALL(:,:,a+1,4) - dec*viesALL(:,:,a,4) + 3600.*viesALL(:,:,a,4))/3600. !vies do Td [C]
   aux1(:,:)=p00 * ((basic_g(ng)%pi0(2,ia:iz,ja:jz) + basic_g(ng)%pp(2,ia:iz,ja:jz)) * cpi) ** cpor  ! Press [Pa]
   aux2(:,:)=6.122*exp(17.67*(viesM(:,:,4))/(viesM(:,:,4)+243.5))  ! vies do "e" em [mb]
   viesM(:,:,4)=0.622*aux2(:,:)/(aux1(:,:)-0.38*aux2(:,:))  ! vies de "q" [kg/kg]

   !--- Pressao ---
   viesM(:,:,5)=(dec*viesALL(:,:,a+1,5) - dec*viesALL(:,:,a,5) + 3600.*viesALL(:,:,a,5))/3600.

   DEALLOCATE(aux1,aux2)
ENDIF
!--- Remocao de VIES ---}

runtypeB=runtype
hfilinB=hfilin

swdownB(:,:)   =radiate_g(ng)%rshort(ia:iz, ja:jz) 
lwdownB(:,:)   =radiate_g(ng)%rlong(ia:iz, ja:jz)

!--- TMP ate resolver o problema de borda da radiacao ---{
swdownB(ia,:)=swdownB(ia+1,:)
swdownB(iz,:)=swdownB(iz-1,:)
swdownB(:,ja)=swdownB(:,ja+1)
swdownB(:,jz)=swdownB(:,jz-1)

lwdownB(ia,:)=lwdownB(ia+1,:)
lwdownB(iz,:)=lwdownB(iz-1,:)
lwdownB(:,ja)=lwdownB(:,ja+1)
lwdownB(:,jz)=lwdownB(:,jz-1)
!------------------------------------------------}


! 1=UVEL 2=VVEL  3=TEMP  4=UMID  5=PNMM
precipB(:,:)   =pcpgl(ia:iz, ja:jz)
diff_radB(:,:) =0.0 !radiate_g(ng)%rshortdif(ia:iz, ja:jz)  !TMP
upsB(:,:)      =basic_g(ng)%up(2,ia:iz, ja:jz) -viesM(:,:,1)
vpsB(:,:)      =basic_g(ng)%vp(2,ia:iz, ja:jz) -viesM(:,:,2)
tempB(:,:)     =0.5* (basic_g(ng)%theta(1,ia:iz,ja:jz) + basic_g(ng)%theta(2,ia:iz,ja:jz)) * &
                .5 * cpi * (basic_g(ng)%pi0(1,ia:iz,ja:jz) + basic_g(ng)%pi0(2,ia:iz,ja:jz)   &
                + basic_g(ng)%pp(1,ia:iz,ja:jz) + basic_g(ng)%pp(2,ia:iz,ja:jz)) - viesM(:,:,3)

qB(:,:)        =basic_g(ng)%rv(2,ia:iz, ja:jz)/(1+basic_g(ng)%rv(2,ia:iz, ja:jz))-viesM(:,:,4)
pstarB(:,:)    =p00 * ((basic_g(ng)%pi0(2,ia:iz, ja:jz) + basic_g(ng)%pp(2,ia:iz, ja:jz)) * cpi) ** cpor - viesM(:,:,5)*100.

IF (start) THEN
   start=.false.

   dump_periodB=nint(frqhis)
   output_periodB=nint(frqanl)
   dir_run_idB=trim(afilout)
   
   z1_uvB=zt(2) * grid_g(ng)%rtgt(ia,ja)  !Utilizando a altura do primeiro ponto (a variacao eh inferior a 1 metro)
   z1_tqB=2.0

   sm_levelsB=nzg
   IF (.not. ALLOCATED(dzsoilB)) ALLOCATE(dzsoilB(sm_levelsB))
   IF (.not. ALLOCATED(sthuB)) ALLOCATE(sthuB(sm_levelsB,nxB,nyB))
   IF (.not. ALLOCATED(tsoilB)) ALLOCATE(tsoilB(sm_levelsB,nxB,nyB))
  
   dzsoilB(1)=-1.0*slz(sm_levelsB)
   DO k=2,sm_levelsB
       dzsoilB(k)=-1.0*slz(sm_levelsB-k+1) + slz(sm_levelsB-k+2)
   ENDDO

   CALL sfcdata
   
   !--- ACOPLANDO - BRAMS p/ JULES ------{
   DO j=ja,jz
      DO i=ia,iz
         DO k=1,sm_levelsB
            nsoil = nint(leaf_g(ng)%soil_text(k,i,j,2))
            sthuB(sm_levelsB-k+1,i,j)=leaf_g(ng)%soil_water(k,i,j,2)/slmsts(nsoil)

            CALL qwtk(leaf_g(ng)%soil_energy(k,i,j,2),     &
                 leaf_g(ng)%soil_water(k,i,j,2)*1.e3, slcpd(nsoil),tempk,fracliq)
            tsoilB(sm_levelsB-k+1,i,j)=tempk
         ENDDO
      ENDDO
   ENDDO
      
   !--- Utilizando a media da temperatura do primeiro (abaixo do solo) e segundo nivel sigma ---
   tstarB(:,:) = 0.5* (basic_g(ng)%theta(1,ia:iz,ja:jz) + basic_g(ng)%theta(2,ia:iz,ja:jz)) * &
                .5 * cpi * (basic_g(ng)%pi0(1,ia:iz,ja:jz) + basic_g(ng)%pi0(2,ia:iz,ja:jz)   &
                + basic_g(ng)%pp(1,ia:iz,ja:jz) + basic_g(ng)%pp(2,ia:iz,ja:jz))
   !-------------------------------------} 

   ALLOCATE(glatB(nxB,nyB),glonB(nxB,nyB),land_fracB(nxB,nyB))
   glatB=grid_g(ng)%glat(ia:iz,ja:jz)
   glonB=grid_g(ng)%glon(ia:iz,ja:jz)

      land_fracB=1.0  ! calculando em todos os pontos, pois senao dah problema no history, pois 
                      ! o patch_area estah mudando ao longo do tempo.

   timestepB=nint(dtlong*fat_dtlong)  !DSM foi colocado em sfclyr_jules.f90 um fator multiplicando o dtlong para nao chamar o JULES a cada timeStep do BRAMS
   dti%year=iyear1
   dti%month=imonth1
   dti%day=idate1
   itime1_seg=(itime1/100)*3600 + mod(itime1,100)*60
   dti%time=itime1_seg

   dtf = datetime_advance(dti, nint(timmax))

   if (trim(runtypeB)=='HISTORY') then
      aux=trim(hfilinB(index(hfilinB,'-head.txt',BACK = .TRUE.)-17:index(hfilinB,'-head.txt',BACK = .TRUE.)-1))
      read(aux(1:4),*) dti%year
      read(aux(6:7),*) dti%month
      read(aux(9:10),*) dti%day
      
      read(aux(12:13),*) hh
      read(aux(14:15),*) mm
      dti%time=hh*3600+mm*60
   endif

   main_run_startB=datetime_to_string(dti)
   
   main_run_endB=datetime_to_string(dtf)

   ntimestep=datetime_diff(dtf, dti)/dtlong+1
   ALLOCATE (run_jules(ntimestep))

   ntimestepB=1
   run_jules=.FALSE.
   run_jules(1)=.TRUE.
   DO i=1+fat_dtlong,ntimestep,fat_dtlong
      run_jules(i)=.TRUE.
      ntimestepB=ntimestepB+1
   ENDDO
 
   if(mynum==1) then
      open(unit=66,file='jules.log',status='old',position='append',action='write')
      write(unit=66,fmt='(A)') '---- Inicio da FASE-1 ---'
      close(unit=66)
   endif
   fase=1; call jules_subroutine(nml_dir,fase)  ! faz leitura do namelist (le o ntype)

   !--- Converte o vegetacao do BRAMS (Leaf3) para a do JULES ---{
   ALLOCATE(fracB(nxB,nyB,ntype))
   CALL frac_from_leaf(ntype,nxB,nyB,npatch   &
                        ,leaf_g(ng)%patch_area(ia:iz,ja:jz,:)  &
                        ,leaf_g(ng)%leaf_class(ia:iz,ja:jz,:)  &
                        ,fracB)
   !-------------------------------------------------------------}

   if(mynum==1) then
      open(unit=66,file='jules.log',status='old',position='append',action='write')
      write(unit=66,fmt='(A)') '---- Inicio da FASE-2 ---'
      close(unit=66)
   endif
   fase=2; call jules_subroutine(nml_dir,fase)  ! finaliza a inicializacao do JULES

   IF (ALLOCATED(dzsoilB)) DEALLOCATE(dzsoilB)
   IF (ALLOCATED(sthuB)) DEALLOCATE(sthuB)
   IF (ALLOCATED(tsoilB)) DEALLOCATE(tsoilB)

ENDIF

IF (.not. run_jules(istp)) RETURN

if(mynum==1 .and. time==0) then
   open(unit=66,file='jules.log',status='old',position='append',action='write')
   write(unit=66,fmt='(A)') '---- Inicio da FASE-3 ---'
   close(unit=66)
endif
!--- ACOPLANDO - BRAMS p/ JULES ------{
!DSM    IF ( (chemical_mechanism .eq. 'RELACS_TUV' .or. chemical_mechanism .eq. 'CO2') &
!DSM          .and. CHEMISTRY == 0) THEN
!DSM       aerotype%co2_3d_ij(1:nxB, 1:nyB) = chem1_g(CO2,ng)%sc_p (2, ia:iz, ja:jz) * 1.E-9
!DSM    ELSE
!DSM       aerotype%co2_3d_ij = 384.
!DSM    END IF
!-------------------------------------}

fase=3; call jules_subroutine(nml_dir,fase) !--- para cada timstep do BRAMS

!--- ACOPLANDO - JULES p/ BRAMS ------{
!--- Acoplando o albedt ---{
radiate_g(ng)%albedt(:,:)=(land_albedo_ij(:,:,1)+land_albedo_ij(:,:,2)+land_albedo_ij(:,:,3)+land_albedo_ij(:,:,4))/4

ALLOCATE(rlongupJ(land_pts))
rlongupJ(:)=sbcon * surftiles_to_gbm(emis_surft * progs%tstar_surft**4, ainfo)
DO l=1,land_pts
   j = ( ainfo%land_index(l)-1 ) / row_length + 1
   i = ainfo%land_index(l) - ( j-1 ) * row_length
   IF (i<1 .or. i>nxB .or. j<1 .or. j>nyB .or. i+ia-1>iz .or. j+ja-1>jz) THEN
      PRINT*, "ERRO... conversao incorreta de l para i,j -> l,i,j=",l,i,j
      STOP 'STOP in: sfclyr_jules.f90 (f)'
   ENDIF
   !--- Acoplando o rlongup ---{
   radiate_g(ng)%rlongup(i+ia-1,j+ja-1) = rlongupJ(l)

   !--- Acoplando a umidade do solo ---
   DO k=1,sm_levelsB
      nsoil = nint(leaf_g(ng)%soil_text(k,i,j,2))
         leaf_g(ng)%soil_water(k,i+ia-1,j+ja-1,2:npatch)=psparms%sthu_soilt(l,1,sm_levelsB-k+1)*slmsts(nsoil) + &
                                                         psparms%sthf_soilt(l,1,sm_levelsB-k+1)*slmsts(nsoil)
   ENDDO

   !--- Acoplando o lai ---
   leaf_g(ng)%veg_lai(i+ia-1,j+ja-1,2:npatch) = SUM(progs%lai_pft(l,:) * ainfo%frac_surft(l,1:npft))

ENDDO
DEALLOCATE(rlongupJ)

!--- Acoplando fluxo de calor latente ---
turb_g(ng)%sflux_r(ia:iz,ja:jz) = sf_diag%latent_heat(:,:)/alvl

!--- Acoplando fluxo de calor sensivel ---
turb_g(ng)%sflux_t(ia:iz,ja:jz) = ftl_1_ij(:,:)/cp

ALLOCATE(dens2(nxB,nyB))
!--- Acoplando fluxo sflux_u e sflux_v ---
dens2(:,:) = (basic_g(ng)%dn0(1,ia:iz,ja:jz) + basic_g(ng)%dn0(2,ia:iz,ja:jz)) * .5
turb_g(ng)%sflux_u(ia:iz,ja:jz) = -1*taux_1_ij(:,:)/dens2(:,:)
turb_g(ng)%sflux_v(ia:iz,ja:jz) = -1*tauy_1_ij(:,:)/dens2(:,:)

DO j=ja,jz
   DO i=ia,iz
      !--- Acoplando ustar, tstar, rstar ---
      ustar2 = max(0.000001,sqrt( sqrt( (turb_g(ng)%sflux_u(i+ia-1,j+ja-1))**2 + (turb_g(ng)%sflux_v(i+ia-1,j+ja-1))**2 )))
      tstar2 = ftl_1_ij(i,j)/(dens2(i,j) * cp * ustar2)
      rstar2 = sf_diag%latent_heat(i,j)/(dens2(i,j) * alvl * ustar2)
      leaf_g(ng)%ustar(i+ia-1,j+ja-1,1:npatch)=ustar2
      leaf_g(ng)%tstar(i+ia-1,j+ja-1,1:npatch)=tstar2
      leaf_g(ng)%rstar(i+ia-1,j+ja-1,1:npatch)=rstar2

      !--- Acoplando fluxo sflux_w ---
      zts2 = zt(2) * grid_g(ng)%rtgt(i+ia-1,j+ja-1)
      ths2 = basic_g(ng)%theta(2,i+ia-1,j+ja-1)
      gzotheta = g * zts2 / ths2
      zoverl = gzotheta * vonk * tstar2 / (ustar2 * ustar2)
      IF (zoverl < 0.) THEN
         cx = zoverl * sqrt(sqrt(1. - 15. * zoverl))
      ELSE
         cx = zoverl / (1.0 + 4.7 * zoverl)
      ENDIF      
      psin = sqrt((1.-2.86 * cx) / (1. + cx * (-5.39 + cx * 6.998 )))
      ip=2
      turb_g(ng)%sflux_w(i+ia-1,j+ja-1) = (0.27 * max(6.25 * (1. - cx) * &
        psin,wtol)- 1.18 * cx * psin) * ustar2 * ustar2 * leaf_g(ng)%patch_area(i+ia-1,j+ja-1,ip)    
   ENDDO
ENDDO
DEALLOCATE(dens2)
!-------------------------------------}

!--- ESCREVENDO AS VARIAVEIS DO JULES NO OUTPUT DO BRAMS ---{
jules_g(ng)%u10mj(ia:iz,ja:jz)=sf_diag%u10m(:,:)
jules_g(ng)%v10mj(ia:iz,ja:jz)=sf_diag%v10m(:,:)
jules_g(ng)%t2mj(ia:iz,ja:jz)=sf_diag%t1p5m(:,:)
jules_g(ng)%rv2mj(ia:iz,ja:jz)=sf_diag%q1p5m(:,:)/(1-sf_diag%q1p5m(:,:))


 DO l=1,land_pts
   j = ( ainfo%land_index(l)-1 ) / row_length + 1
   i = ainfo%land_index(l) - ( j-1 ) * row_length


   IF (i<1 .or. i>nxB .or. j<1 .or. j>nyB .or. i+ia-1>iz .or. j+ja-1>jz) THEN
      PRINT*, "ERRO... conversao incorreta de l para i,j -> l,i,j=",l,i,j
      STOP 'STOP in: sfclyr_jules.f90 (g)'
   ENDIF

   jules_g(ng)%gpp(i+ia-1,j+ja-1)=trifctltype%gpp_gb(l)
   jules_g(ng)%resp_p(i+ia-1,j+ja-1)=trifctltype%resp_p_gb(l)
   jules_g(ng)%npp(i+ia-1,j+ja-1)=trifctltype%npp_gb(l)
   jules_g(ng)%resp_s(i+ia-1,j+ja-1)=SUM(sum(trifctltype%resp_s_soilt(l,1,:,:), 2), 1)
 ENDDO
!-----------------------------------------------------------}

!IF ( end_of_run ) THEN  !esta funcao eh do JULES - pode tambem utilizar o ultimo timestep do BRAMS
!LFR
if(time>=timmax) then
   if(mynum==1) then 
      open(unit=66,file='jules.log',status='old',position='append',action='write')
      write(unit=66,fmt='(A)') '---- Inicio da FASE-4 ---'
      close(unit=66)
   endif


   fase=4; call jules_subroutine(nml_dir,fase) ! apos o ultimo timestep do BRAMS
ENDIF 
 
end subroutine sfclyr_jules


!{DSM
!--------------------------------------------------------------------------------------------------
!--- Encontra a fracao de vegetacao a partir do mapa da fracao de vegetacao definida pelo leaf3 ---
!--------------------------------------------------------------------------------------------------
SUBROUTINE frac_from_leaf( ntype,nx,ny,npatch,patch_area,leaf_class,frac )
   IMPLICIT NONE
   INTEGER, INTENT(IN)  :: ntype,nx,ny,npatch
   REAL, INTENT(IN)     :: patch_area(nx,ny,npatch), leaf_class(nx,ny,npatch)
   REAL, INTENT(OUT)    :: frac(nx,ny,ntype)

   INTEGER              :: j,i,tJ,n
   CHARACTER (LEN=80)   :: veg(ntype)
   CHARACTER (LEN=2)    :: tB_str          
   
   !--- Convertendo o tipo de vegetacao do BRAMS para o JULES ---
   CALL brams2jules(veg,ntype)
  
   if (ntype == 9 ) then
      frac=0.
                         
      DO j=1,ny
         DO i=1,nx
            IF (i > nx .or. j > ny) THEN
               PRINT*, 'ERRO!!! i > nx ou j > ny - i, nx, j, ny =',i, nx, j, ny
               STOP 'STOP in: sfclyr_jules.f90 (h)'
            ENDIF

            DO n=2,npatch
               WRITE(tB_str,'(i2.2)') nint(leaf_class(i,j,n))
               
               !--- Encontrando o tipo correspondente ao JULES ---
               DO tJ=1,ntype+1    !--- o +1 eh apenas para checar se foi encontrado (condicao abaixo)
                  IF (index(veg(tJ),tB_str)/=0) exit
               ENDDO
               
               !--- Checando se encontrou um indice valido para a vegetacao do JULES ---
               IF (tJ>ntype) THEN
                  PRINT*, 'ERRO!!! Nao foi encontrado uma correspondencia entre BRAMS e JULES'
                  STOP 'STOP in: sfclyr_jules.f90 (i)'
               ENDIF
               
               frac(i,j,tJ)=frac(i,j,tJ) + max(0.,patch_area(i,j,n))

            ENDDO  !--- DO n=1,npatch
               
            frac(i,j,7)=frac(i,j,7)+max(0.,patch_area(i,j,1))
            
            WHERE ( frac(i,j,:) <= 1.00004E-06 ) frac(i,j,:) = 1.00004E-06  ! valor minimo que o JULES trabalha

            frac(:,:,9)=0.0 !--- Mas ice fraction deve iniciar com zero.
            
            n=1
            DO WHILE (sum(frac(i,j,:)) > 1.0)
               IF (frac(i,j,n)>0.01) frac(i,j,n)=frac(i,j,n)-0.01                             !delta=0.01
               n=n+1
               if (n>8) n=1     
            ENDDO
            
            !--- Para garantir que a fracao total seja igual a 1 ---
            frac(i,j,1)=1.-( frac(i,j,2)+frac(i,j,3)+frac(i,j,4)     &
                           +frac(i,j,5)+frac(i,j,6)+frac(i,j,7)+frac(i,j,8)+frac(i,j,9) )
                                       
         ENDDO  !i=1,nx
      ENDDO  !j=1,ny
   elseif (ntype == 13 ) then
      frac=0.
                         
      DO j=1,ny
         DO i=1,nx
            IF (i > nx .or. j > ny) THEN
               PRINT*, 'ERRO!!! i > nx ou j > ny - i, nx, j, ny =',i, nx, j, ny
               STOP 'STOP in: sfclyr_jules.f90 (j)'
            ENDIF

            DO n=2,npatch
               WRITE(tB_str,'(i2.2)') nint(leaf_class(i,j,n))
               
               !--- Encontrando o tipo correspondente ao JULES ---
               DO tJ=1,ntype+1    !--- o +1 eh apenas para checar se foi encontrado (condicao abaixo)
                  IF (index(veg(tJ),tB_str)/=0) exit
               ENDDO
               
               !--- Checando se encontrou um indice valido para a vegetacao do JULES ---
               IF (tJ>ntype) THEN
                  PRINT*, 'ERRO!!! Nao foi encontrado uma correspondencia entre BRAMS e JULES'
                  STOP 'STOP in: sfclyr_jules.f90 (k)'
               ENDIF
               
               frac(i,j,tJ)=frac(i,j,tJ) + max(0.,patch_area(i,j,n))

            ENDDO  !--- DO n=1,npatch
               
            frac(i,j,11)=frac(i,j,11)+max(0.,patch_area(i,j,1))
            
            WHERE ( frac(i,j,:) <= 1.00004E-06 ) frac(i,j,:) = 1.00004E-06  ! valor minimo que o JULES trabalha

            frac(:,:,13)=0.0 !--- Mas ice fraction deve iniciar com zero.
            
            n=1
            DO WHILE (sum(frac(i,j,:)) > 1.0)
               IF (frac(i,j,n)>0.01) frac(i,j,n)=frac(i,j,n)-0.01                             !delta=0.01
               n=n+1
               if (n>12) n=1     
            ENDDO
            
            !--- Para garantir que a fracao total seja igual a 1 ---
            frac(i,j,1)=1.-( frac(i,j,2)+frac(i,j,3)+frac(i,j,4)     &
                           +frac(i,j,5)+frac(i,j,6)+frac(i,j,7)+frac(i,j,8)  &
                           +frac(i,j,9)+frac(i,j,10)+frac(i,j,11)+frac(i,j,12)+frac(i,j,13) )
                                       
         ENDDO  !i=1,nx
      ENDDO  !j=1,ny
   else
      PRINT*, "ntype=",ntype
      PRINT*, "ATENCAO... ntype <> 9 ou de 13, Deve-se ajustar a subrotina brams2jules em init_frac2brams.f90"
      STOP 'STOP in: sfclyr_jules.f90 (l)'
   endif
  
END SUBROUTINE frac_from_leaf

!------------------------------------------------------------------------------
!--- Convertendo o tipo de vegetacao do BRAMS para o JULES ---
!-------------------------------------------------------------------------------
SUBROUTINE brams2jules(veg,ntype)
   IMPLICIT NONE
   INTEGER, INTENT(in)              :: ntype
   CHARACTER (len=80), INTENT(out)  :: veg(ntype)
   
   if (ntype == 9 ) then
      veg(1)='06 07 20'       !tJ=1 => BT=broadleaf trees
      veg(2)='04 05 14'       !tJ=2 => NT=needleleaf trees
      veg(3)='15 08 16 17'    !tJ=3 => C3G=C3 (temperate) grass
      veg(4)='09'             !tJ=4 => C4G=C4 (tropical) grass
      veg(5)='11 12 13 18'    !tJ=5 => shrub
      veg(6)='19 21'          !tJ=6 => urban
      veg(7)='00 01'          !tJ=7 => lake=inland water
      veg(8)='03 10'          !tJ=8 => soil=bare soil
      veg(9)='02'             !tJ=9 => ice
   elseif (ntype == 13 ) then
      veg(1)='07 20'         !tJ=1 => BT="Broadleaf"  tropical 
      veg(2)='14'            !tJ=2 => Broadleaf temperada
      veg(3)='06'            !tJ=3 => Broadleaf deciduas
      veg(4)='04'            !tJ=4 => Needle-leaf" "evergreen
      veg(5)='05'            !tJ=5 => "Needle-leaf" deciduas
      veg(6)='15 08 11 16 17'!tJ=6 => Gramineas - C3
      veg(7)='09'            !tJ=7 => Gramineas - C4
      veg(8)='12 18'         !tJ=8 => shrub - Cerrado "evergreen"
      veg(9)='13'            !tJ=9 => shrub - Cerrado deciduo
      veg(10)='19 21'        !tJ=10 => urban
      veg(11)='00 01'        !tJ=11 => lake=inland water
      veg(12)='03 10'        !tJ=12 => bare soil - Solo nu
      veg(13)='02'           !tJ=13 => ice - Gelo
   else
      PRINT*, "ntype=",ntype
      PRINT*, "ATENCAO... ntype <> 9 ou de 13, Deve-se ajustar a subrotina brams2jules em init_frac2brams.f90"
      STOP  'STOP in: sfclyr_jules.f90 (m)'
   endif
END SUBROUTINE brams2jules

!DSM}
