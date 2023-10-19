!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

!--(DMK-CCATT-INI)----------------------------------------------------------------
subroutine chem_ISAN_file_inv (iyear1,imonth1,idate1,itime1,timmax,CHEM_ASSIM, CHEMISTRY)
!--(DMK-CCATT-OLD)----------------------------------------------------------------
!subroutine ISAN_file_inv (iyear1,imonth1,idate1,itime1,timmax)
!--(DMK-CCATT-OLD)----------------------------------------------------------------
  use ModDateUtils
  use isan_coms, only: &
       IGRIDFL, &
       IAPR, &
       GUESS1ST, &
       IARAWI, &
       IASRFCE, &
       ISAN_INC, &
       NPDATES, &
       IPROC_DATES, &
       I1ST_FLG, &
       IUPA_FLG, &
       ISFC_FLG, &
       IPROC_FLAG, &
       IPROC_NAMES, &
       MAXISFILES, &
       icFileType, &
       icPrefix

  use mem_grid, only: &
       time     
  use dump, only: &
    dumpMessage  

  implicit none

  include "files.h" 
  include "constants.f90"
  include "UseVfm.h"
  
  integer, intent (IN) :: CHEM_ASSIM, CHEMISTRY
  integer :: iyear1,imonth1,idate1,itime1
  real :: timmax

  character(len=*),parameter :: h='**(chem_ISAN_file_inv)**'

  ! Local variables ----------------------------------------------------------
  !
  !  nfgfiles = Number of First Guess FILES   (defined value = zero)
  !  nupfiles = Number of Upper Air FILES     (defined value = zero)
  !  nsffiles = Number of SurFace input FILES (defined value = zero)
  !
  ! -------------------------------------------------------------------------
  integer :: nfgfiles,nc,nf,lnf,nn,ndates,nupfiles,nsffiles,isan_err_flag
  integer :: inyear,inmonth,indate,inhour
  integer :: iyearf,imonthf,idatef,ihourf
  integer :: iyear2,imonth2,idate2,ihour2,ihour1
  real :: tinc

  character(len=14) :: itotdate_fg(maxisfiles),itotdate_up(maxisfiles)  &
       ,itotdate_sf(maxisfiles),itotdates(4*maxisfiles),idate_end
  character(len=f_name_length), dimension(maxisfiles) :: fnames_fg, &
       fnames_up, fnames_sf
  character(len=f_name_length) :: rams_filelist_arg
  logical :: there,fileExist
  integer :: localTime
  character(len=f_name_length) :: sVarName 
  integer :: iyears,imonths,idates,ihours
  integer :: indice

  ! Define variables -------------------------------------------------------
  nfgfiles = 0
  nupfiles = 0
  nsffiles = 0
  ! ------------------------------------------------------------------------

  !          Go through first guess, upper air, surface input files
  !            and make inventory

  if(igridfl.ne.0) then
     if(iapr(1:1).ne.' '.and.iapr(1:1).ne.char(0) ) then

!--(DMK-CCATT-INI)----------------------------------------------------------------
     print*,'----> performing the data inventory'

        nfgfiles=-1
        nc=len_trim(iapr)
        indice = 1

!!$ if(guess1st.eq.'PRESS') then
!!$ rams_filelist_arg = iapr(1:nc)//'????-??-??-????'
!!$ call RAMS_filelist(fnames_fg, rams_filelist_arg, nfgfiles)

        nfgfiles = ((timmax/3600) / (isan_inc/100)) + 1    
        
        
        call date_add_to(iyear1,imonth1,idate1,itime1*100  &
             ,time,'s',iyears,imonths,idates,ihours)

        localTime = itime1
        indice = 1  

        do nf=1,nfgfiles

           call date_add_to(iyears, imonths, idates, localTime*100,  &
                0., 's', iyears, imonths, idates, ihours)

!--(DMK-CCATT-INI)----------------------------------------------------------------
           if(guess1st.eq.'PRESS')  then 
	    if(CHEM_ASSIM == 1 .and. CHEMISTRY >= 0)  then 
	        write(sVarName,100) iapr(1:len_trim(iapr)),'-',iyears,'-',imonths,'-',idates,'-',ihours/100   
                if (useVfm) then
                   sVarName = trim(sVarName) // '.vfm'
                else
                   sVarName = trim(sVarName) // '.bin'
                end if
	    else 
	        write(sVarName,101) iapr(1:len_trim(iapr)),iyears,'-',imonths,'-',idates,'-',ihours/100   
	    endif
	   endif
	   100 format(a,a1,i4.4,a1,i2.2,a1,i2.2,a1,i4.4)
	   101 format(a,   i4.4,a1,i2.2,a1,i2.2,a1,i4.4)


           if(guess1st.eq.'RAMS')   then 
	      write(sVarName,102) iapr(1:len_trim(iapr)),'-A-',iyears,'-',imonths,'-',idates,'-',ihours   
              sVarName = sVarName // "-head.txt"
	      102 format(a,a3,i4.4,a1,i2.2,a1,i2.2,a1,i6.6)
           endif


!--(DMK-CCATT-OLD)----------------------------------------------------------------
!           write(sVarName,100) iapr(1:len_trim(iapr)),iyears,'-',imonths,'-',idates,'-',ihours/100   
!100        format(a,i4.4,a1,i2.2,a1,i2.2,a1,i4.4)
!--(DMK-CCATT-FIM)----------------------------------------------------------------
!           if(guess1st.eq.'RAMS') then
!              sVarName = sVarName // "-head.txt"
!--(DMK-CCATT-INI)----------------------------------------------------------------
!           else
!srf - introduzindo diferenças entre assimilaçao de dp + chem e dp meteo somente 
!	     if(CHEM_ASSIM == 1 .and. CHEMISTRY >= 0) then 
!	       sVarName = TRIM(sVarName) // '.vfm'
!	     else
!	       sVarName = TRIM(sVarName)
!	     endif
!srf
!--(DMK-CCATT-FIM)----------------------------------------------------------------	      
!           endif


!--(DMK-CCATT-INI)----------------------------------------------------------------
           print*,'sVarName=',TRIM(sVarName)
!--(DMK-CCATT-FIM)----------------------------------------------------------------

           inquire(file=sVarName(1:len_trim(sVarName)),exist=there)

           if (there) then
              fnames_fg(indice) = trim(sVarName)
              indice = indice + 1
           endif

!--srf change LE 18000 to LE 2100 for  3/3 ivars
           if(localTime .LE. 2100)then
              localTime = localTime + isan_inc
           else
              localTime = 000000
              localTime = localTime + isan_inc
           endif
        enddo

        nfgfiles = indice - 1

!!$ endif

!!$ if(guess1st.eq.'RAMS') then
!!$ rams_filelist_arg = iapr(1:nc)//'*-head.txt'
!!$ call RAMS_filelist(fnames_fg, rams_filelist_arg, nfgfiles)
!!$ endif

        !      do nf=1,nfgfiles
        !         print*,'reading fg:',nc,nfgfiles,fnames_fg(nf)
        !      enddo

        if(nfgfiles.gt.maxisfiles) then
           print*,'too many first guess files'
           call fatal_error('lots_of_first_guess')
        endif

        do nf=1,nfgfiles

           lnf=len_trim(fnames_fg(nf))

           if(guess1st.eq.'PRESS') then
              read(fnames_fg(nf)(lnf-18:lnf-4),20) inyear,inmonth,indate,inhour
           end if
!--(DMK-CCATT-FIM)----------------------------------------------------------------

           ! form of a-A-2000-07-01-060000-head.txt
           if(guess1st.eq.'RAMS')  then
              read(fnames_fg(nf)(lnf-25:lnf-10),20) inyear,inmonth,indate,inhour
           end if

20         format(i4,1x,i2,1x,i2,1x,i4)

           call date_make_big(inyear,inmonth,indate,inhour*100,itotdate_fg(nf))

        enddo

        call RAMS_dintsort(nfgfiles,itotdate_fg,fnames_fg)

        !      do nf=1,nfgfiles
        !         print*,'fg files:',nf,itotdate_fg(nf),fnames_fg(nf)
        !      enddo

     endif
  endif


  if(igridfl.ne.3) then
     if(iarawi(1:1).ne.' '.and.iarawi(1:1).ne.char(0) ) then

        nupfiles=-1
        nc=len_trim(iarawi)
!!$         rams_filelist_arg = iarawi(1:nc)//'????-??-??-????'
!!$         call RAMS_filelist(fnames_up, rams_filelist_arg, nupfiles)

        nupfiles = ((timmax/3600) / (isan_inc/100)) + 1    

        call date_add_to(iyear1,imonth1,idate1,itime1*100  &
             ,time,'s',iyears,imonths,idates,ihours)

        localTime = itime1
        indice = 1  

        do nf=1,nupfiles

           call date_add_to(iyears, imonths, idates, localTime*100,  &
                0., 's', iyears, imonths, idates, ihours)

           write(sVarName,400) iarawi(1:len_trim(iarawi)),iyears,'-',imonths,'-',idates,'-',ihours/100   
400        format(a,i4.4,a1,i2.2,a1,i2.2,a1,i4.4)

           inquire(file=sVarName(1:len_trim(sVarName)),exist=there)

           if (there) then
              fnames_up(indice) = trim(sVarName)
              indice = indice + 1
           endif

           if(localTime .LE. 1800)then
              localTime = localTime + isan_inc
           else
              localTime = 000000
              localTime = localTime + isan_inc
           endif
        enddo

        nupfiles = indice - 1

        if(nupfiles.gt.maxisfiles) then
           print*,'too many upper air files'
           call fatal_error('lots_of_upper_air')
        endif

        do nf=1,nupfiles
           lnf=len_trim(fnames_up(nf))
           read(fnames_up(nf)(lnf-14:lnf),20) inyear,inmonth,indate,inhour

           call date_make_big(inyear,inmonth,indate,inhour*100,itotdate_up(nf))

        enddo

        call RAMS_dintsort(nupfiles,itotdate_up,fnames_up)

        !      do nf=1,nupfiles
        !         print*,'up files:',nf,itotdate_up(nf),fnames_up(nf)
        !      enddo

     endif


     if(iasrfce(1:1).ne.' '.and.iasrfce(1:1).ne.char(0) ) then

        nsffiles=-1
        nc=len_trim(iasrfce)
!!$        rams_filelist_arg = iasrfce(1:nc)//'????-??-??-????'
!!$        call RAMS_filelist(fnames_sf, rams_filelist_arg, nsffiles)

        nsffiles = ((timmax/3600) / (isan_inc/100)) + 1    

        call date_add_to(iyear1,imonth1,idate1,itime1*100  &
             ,time,'s',iyears,imonths,idates,ihours)

        localTime = itime1
        indice = 1  

        do nf=1,nsffiles

           call date_add_to(iyears, imonths, idates, localTime*100,  &
                0., 's', iyears, imonths, idates, ihours)

           write(sVarName,300) iasrfce(1:len_trim(iasrfce)),iyears,'-',imonths,'-',idates,'-',ihours/100   
300        format(a,i4.4,a1,i2.2,a1,i2.2,a1,i4.4)

           inquire(file=sVarName(1:len_trim(sVarName)),exist=there)

           if (there) then
              fnames_sf(indice) = trim(sVarName)
              indice = indice + 1
           endif

           if(localTime .LE. 1800)then
              localTime = localTime + isan_inc
           else
              localTime = 000000
              localTime = localTime + isan_inc
           endif
        enddo

        nsffiles = indice - 1

        if(nsffiles.gt.maxisfiles) then
           print*,'too many surface air files'
           call fatal_error('lots_of_surface')
        endif

        do nf=1,nsffiles
           lnf=len_trim(fnames_sf(nf))
           read(fnames_sf(nf)(lnf-14:lnf),20) inyear,inmonth,indate,inhour

           call date_make_big(inyear,inmonth,indate,inhour*100,itotdate_sf(nf))
        enddo

        call RAMS_dintsort(nsffiles,itotdate_sf,fnames_sf)

        !      do nf=1,nsffiles
        !         print*,'sf files:',nf,itotdate_sf(nf),fnames_sf(nf)
        !      enddo

     endif

  endif

  ! put dates in order, removing duplicates

  call RAMS_sort_dint3(nfgfiles,itotdate_fg  &
       ,nupfiles,itotdate_up,nsffiles,itotdate_sf  &
       ,ndates,itotdates)

  call RAMS_unique_dint(ndates,itotdates)

  !print*,'dates:',ndates
  !do nn=1,ndates
  !   print*,'dates:',itotdates(nn)
  !enddo

  !  start printing section
  !--------------------------------------------------------------

  print*,' '
  print*,' '
  print*,' '
  print*,'---------------------------------------------------------'
  write (*,fmt='(A,I1,A,I3,A)')'-- ISAN Input File Date Inventory for icFileType = ',icFileType,', with ',ndates,' dates.'
  print*,'---------------------------------------------------------'
  do nn=1,ndates
    print*,'---- Date:',itotdates(nn)

    if(icFileType==1) then
      call makeGrib2fileName(trim(icPrefix),itime1,isan_inc,nn &
            ,iproc_names(nn,5),iproc_names(nn,6))
      print*,'---- First guess grib2 file:'  &
                ,iproc_names(nn,5)(1:len_trim(iproc_names(nn,5)))
    elseif(icFileType==2) then
      iproc_names(nn,5)=trim(icPrefix)
      print*,'---- First guess netCDF file:'  &
                ,iproc_names(nn,5)(1:len_trim(iproc_names(nn,5)))
    elseif(icFileType==3) then
      call makeGeosfileName(trim(icPrefix),isan_inc,nn &
            ,iproc_names(nn,5),iproc_names(nn,6),itotdates(nn),iyear1,imonth1,idate1,itime1)
    elseif(icFileType==4) then
      call makeGradsfileName(trim(icPrefix),isan_inc,nn &
            ,iproc_names(nn,5),iproc_names(nn,6),itotdates(nn),iyear1,imonth1,idate1,itime1)
    else
      do nf=1,nfgfiles
        if(itotdates(nn).eq.itotdate_fg(nf)) then
           print*,'---- First guess file:'  &
                ,fnames_fg(nf)(1:len_trim(fnames_fg(nf)))
        endif
      enddo
    endif

    do nf=1,nupfiles
        if(itotdates(nn).eq.itotdate_up(nf)) then
           print*,'---- Upper air   file:'  &
                ,fnames_up(nf)(1:len_trim(fnames_up(nf)))
        endif
     enddo
     do nf=1,nsffiles
        if(itotdates(nn).eq.itotdate_sf(nf)) then
           print*,'---- Surface obs file:'  &
                ,fnames_sf(nf)(1:len_trim(fnames_sf(nf)))
        endif
     enddo
     print*,'------------------------------------------------------'
  enddo
  !stop 666
  !--------------------------------------------------------------

  ! Find dates we are going to process

  !print*,'start dates:',timmax,iyear1,imonth1,idate1,itime1,isan_inc

  !Find end date
  call date_add_to(iyear1,imonth1,idate1,itime1*100  &
       ,timmax,'s',iyearf,imonthf,idatef,ihourf)
  call date_make_big(iyearf,imonthf,idatef,ihourf,idate_end)
  !print*,'end dates:',iyearf,imonthf,idatef,ihourf,idate_end

  ihour1=itime1*100
  tinc= (isan_inc/100) * 60.  + mod(isan_inc,100)
  npdates = 1
  call date_add_to (iyear1,imonth1,idate1,ihour1  &
       ,tinc*(npdates-1),'m',iyear2,imonth2,idate2,ihour2)
  call date_make_big(iyear2,imonth2,idate2,ihour2  &
       ,iproc_dates(npdates))
  do while (iproc_dates(npdates) .lt. idate_end)
     npdates = npdates + 1
     
!--(DMK-CCATT-INI)---------------------------------------------------------     
!   print*,npdates,iproc_dates(npdates),maxisfiles,tinc
!--(DMK-CCATT-FIM)---------------------------------------------------------     

     call date_add_to (iyear1,imonth1,idate1,ihour1  &
          ,tinc*(npdates-1),'m',iyear2,imonth2,idate2,ihour2)
     call date_make_big(iyear2,imonth2,idate2,ihour2  &
          ,iproc_dates(npdates))
  end do


  !   We have the dates to process. Find if files exist and set overall
  !     "go ahead" flag. Put filenames in iproc_names array if they will be used.

  !      iproc_flag array info:
  !      1) Process?     0=no, 1=yes
  !      2) First guess? 0=no file, 1=exists, 2=exists, don't use, 3=interpolate
  !      3) Upper air?   0=no file, 1=exists, 2=exists, don't use
  !      4) Surface?     0=no file, 1=exists, 2=exists, don't use

  isan_err_flag=0

  print*,' '
  print*,' '
  print*,' '
  print*,'---------------------------------------------------------'
  print*,'-----------  ISAN Processing Information    -------------'
  print*,'---------------------------------------------------------'
  print*,'----    Flags:  IGRIDFL =',igridfl,'  I1ST_FLG=',i1st_flg
  print*,'----    Flags:  IUPA_FLG=',iupa_flg,'  ISFC_FLG=',isfc_flg
  print*,'---------------------------------------------------------'

  do nn=1,npdates
     iproc_flag(nn,1)=0
     iproc_flag(nn,2)=0
     iproc_flag(nn,3)=0
     iproc_flag(nn,4)=0

     print*,'------------------------------------------------------'
     print*,'---- Date:', iproc_dates(nn)

     iproc_flag(nn,1)=1
     
     if(icFileType==1) then
      print *,itime1; call flush(6)
      call makeGrib2fileName(trim(icPrefix),itime1,isan_inc,nn &
            ,iproc_names(nn,5),iproc_names(nn,6))
      print*,'---- First guess grib2 (NCEP) file:'  &
                ,iproc_names(nn,5)(1:len_trim(iproc_names(nn,5)))
     elseif(icFileType==2) then
        iproc_names(nn,5)=trim(icPrefix)
       print*,'---- First guess nc4 (NASA) file:'  &
                ,iproc_names(nn,5)(1:len_trim(iproc_names(nn,5)))       
     elseif(icFileType==3) then
        call makeGeosfileName(trim(icPrefix),isan_inc,nn &
            ,iproc_names(nn,5),iproc_names(nn,6),iproc_dates(nn),iyear1,imonth1,idate1,itime1)
     elseif(icFileType==4) then
        call makeGradsfileName(trim(icPrefix),isan_inc,nn &
            ,iproc_names(nn,5),iproc_names(nn,6),iproc_dates(nn),iyear1,imonth1,idate1,itime1)
     endif

    if(icFileType==1) then

        inquire(file=iproc_names(nn,5)(1:len_trim(iproc_names(nn,5))), exist=fileExist )
        if(.not. fileExist)  iErrNumber=dumpMessage(c_tty,c_yes,h,'468' &
              ,c_fatal,'File '//iproc_names(nn,5)(1:len_trim(iproc_names(nn,5))) &
              //' not found. Please, verify and solve it!')
        write(*,fmt='(A)') '---- First guess file' &
              //iproc_names(nn,5)(1:len_trim(iproc_names(nn,5)))//' exists.'
        iproc_flag(nn,2)=1
    elseif(icFileType==2) then
        inquire(file=iproc_names(nn,5)(1:len_trim(iproc_names(nn,5))), exist=fileExist )
        if(.not. fileExist)  iErrNumber=dumpMessage(c_tty,c_yes,h,'468' &
              ,c_fatal,'File '//iproc_names(nn,5)(1:len_trim(iproc_names(nn,5))) &
              //' not found. Please, verify and solve it!')
        write(*,fmt='(A)') '---- First guess file' &
              //iproc_names(nn,5)(1:len_trim(iproc_names(nn,5)))//' exists.'
        iproc_flag(nn,2)=1
    elseif(icFileType==3) then
        inquire(file=iproc_names(nn,5)(1:len_trim(iproc_names(nn,5))), exist=fileExist )
        if(.not. fileExist)  iErrNumber=dumpMessage(c_tty,c_yes,h,'468' &
              ,c_fatal,'File '//iproc_names(nn,5)(1:len_trim(iproc_names(nn,5))) &
              //' not found. Please, verify and solve it!')
        write(*,fmt='(A)') '---- First guess file' &
              //iproc_names(nn,5)(1:len_trim(iproc_names(nn,5)))//' exists.'
        iproc_flag(nn,2)=1
    elseif(icFileType==4) then
        write(*,fmt='(A)') 'Checking for '//iproc_names(nn,5)(1:len_trim(iproc_names(nn,5)))
        inquire(file=iproc_names(nn,5)(1:len_trim(iproc_names(nn,5))), exist=fileExist )
        if(.not. fileExist)  iErrNumber=dumpMessage(c_tty,c_yes,h,'468' &
              ,c_fatal,'File '//iproc_names(nn,5)(1:len_trim(iproc_names(nn,5))) &
              //' not found. Please, verify and solve it!')
        write(*,fmt='(A)') '---- First guess file' &
              //iproc_names(nn,5)(1:len_trim(iproc_names(nn,5)))//' exists.'
        iproc_flag(nn,2)=1
    else

      ! first guess
      iproc_flag(nn,2)=0
      do nf=1,nfgfiles
         if(iproc_dates(nn).eq.itotdate_fg(nf)) then
            print*,'---- First guess file exists.'
            if(igridfl.eq.0) then
               iproc_flag(nn,2)=2
            else
               iproc_flag(nn,2)=1
            endif
            iproc_names(nn,1)=fnames_fg(nf)
            goto 71
         endif
      enddo 
      ! file doesn't exist
      if(i1st_flg.eq.1) then
         print*,'---- First guess file does not exist.'  &
              ,' Will not process this time.'
         iproc_flag(nn,1)=0
      elseif(i1st_flg.eq.2) then
         isan_err_flag=1
         print*,'---- First guess file does not exist.'  &
              ,' Will STOP run.'
         iproc_flag(nn,1)=0
      elseif(i1st_flg.eq.3) then
         print*,'---- First guess file does not exist.'  &
              ,' Will attempt interpolation.'
         iproc_flag(nn,2)=3
      endif 

    endif
71   continue

     ! upper air
     iproc_flag(nn,3)=0
     do nf=1,nupfiles
        if(iproc_dates(nn).eq.itotdate_up(nf)) then
           print*,'---- Upper air file exists.'
           if(igridfl.eq.3) then
              iproc_flag(nn,3)=2
           else
              iproc_flag(nn,3)=1
           endif
           iproc_names(nn,2)=fnames_up(nf)
           goto 72
        endif
     enddo

     ! file doesn't exist
     if(iupa_flg.eq.1) then
        print*,'---- Upper air file does not exist.'  &
             ,' Will not process this time.'
        iproc_flag(nn,1)=0
     elseif(iupa_flg.eq.2) then
        isan_err_flag=1
        print*,'---- Upper air file does not exist.',' Will STOP run.'
        iproc_flag(nn,1)=0
     elseif(iupa_flg.eq.3) then
        print*,'---- Upper air file does not exist.',' Will try running without.'
        iproc_flag(nn,3)=3
     endif
72   continue

     ! surface
     iproc_flag(nn,4)=0
     do nf=1,nsffiles
        if(iproc_dates(nn).eq.itotdate_sf(nf)) then
           print*,'---- Surface obs file exists.'
           if(igridfl.eq.3) then
              iproc_flag(nn,4)=2
           else
              iproc_flag(nn,4)=1
           endif
           iproc_names(nn,3)=fnames_sf(nf)
           goto 73
        endif
     enddo

     ! file doesn't exist
     if(isfc_flg.eq.1) then
        print*,'---- Surface file does not exist.',' Will not process this time.'
        iproc_flag(nn,1)=0
     elseif(isfc_flg.eq.2) then
        isan_err_flag=1
        print*,'---- Surface file does not exist.',' Will STOP run.'
        iproc_flag(nn,1)=0
     elseif(isfc_flg.eq.3) then
        print*,'---- Surface file does not exist.',' Will try running without.'
        iproc_flag(nn,4)=3
     endif
73   continue

  enddo
  print*,'---------------------------------------------------------'
  if (isan_err_flag.ne.0) then
     print*,'ISAN run stopping because of errors!'
     print*,'See previous output listing.'
     call fatal_error('isan_file_errors')
  endif

!--(DMK-CCATT-INI)----------------------------------------------------------------
end subroutine chem_ISAN_file_inv
!--(DMK-CCATT-OLD)----------------------------------------------------------------
!end subroutine ISAN_file_inv
!--(DMK-CCATT-FIM)----------------------------------------------------------------
subroutine makeGrib2fileName(prefix,itime1,isan_inc,nn,gribFileName,invFileName)
    !Header_vai_aqui!
    implicit none

    include "constants.f90"

    !Parameters (constants)
    character(len=*),parameter :: h='**(makeGrib2fileName)**'

    ! Input/Output variables
    character(len=*), intent(in) :: prefix 
    !#
    integer,intent(in)    :: itime1
    !#
    integer,intent(in)    :: isan_inc
    !#
    integer,intent(in)    :: nn
    !#
    character(len=*),intent(out) :: gribFileName
    !#
    character(len=*),intent(out) :: invFileName
    !#
    !Local variables
    character(len=3) :: ctime

    !Code
    !print *,'In: ',itime1,nn,isan_inc
    write(ctime,fmt='(I3.3)') int(itime1/100)+(nn-1)*isan_inc/100
    gribFileName=trim(prefix)//ctime
    invFileName=trim(prefix)//ctime//'.inv'

    
end subroutine makeGrib2fileName

subroutine makeGradsfileName(prefix,isan_inc,nn,gradsFileName,invFileName,iproc_date,iyear1,imonth1,idate1,itime1)
    !Header_vai_aqui!
    implicit none

    include "constants.f90"

    !Parameters (constants)
    character(len=*),parameter :: h='**(makeGradsfileName)**'

    ! Input/Output variables
    character(len=*), intent(in) :: prefix 
    !#
    character(len=*), intent(in) :: iproc_date 
    !#
    integer,intent(in)    :: itime1
    !#
    integer,intent(in)    :: iyear1
    !#
    integer,intent(in)    :: imonth1
    !#
    integer,intent(in)    :: idate1
    !#    
    integer,intent(in)    :: isan_inc
    !#
    integer,intent(in)    :: nn
    !#
    character(len=*),intent(out) :: gradsFileName
    !#
    character(len=*),intent(out) :: invFileName
    !#
    !Local variables
    character(len=2) :: ctime,cdate,cmonth
    character(len=4) :: cyear!,cntime


    !Code
    !print *,'In: nn,isan_inc,iproc_date,iyear1,imonth1,idate1,itime1= ',nn,isan_inc,iproc_date,iyear1,imonth1,idate1,itime1
    write(ctime,fmt='(I2.2)') int(itime1/100)
    write(cyear,fmt='(I4.4)') iyear1
    write(cmonth,fmt='(I2.2)') imonth1
    write(cdate,fmt='(I2.2)') idate1
    gradsFileName=trim(prefix)//iproc_date(1:10)//'.gra'
    invFileName=trim(prefix)//ctime//'.inv'
    
end subroutine makeGradsfileName


subroutine makeGeosfileName(prefix,isan_inc,nn,geosFileName,invFileName,iproc_date,iyear1,imonth1,idate1,itime1)
    !Header_vai_aqui!
    implicit none

    include "constants.f90"

    !Parameters (constants)
    character(len=*),parameter :: h='**(makeGeosfileName)**'

    ! Input/Output variables
    character(len=*), intent(in) :: prefix 
    !#
    character(len=*), intent(in) :: iproc_date 
    !#
    integer,intent(in)    :: itime1
    !#
    integer,intent(in)    :: iyear1
    !#
    integer,intent(in)    :: imonth1
    !#
    integer,intent(in)    :: idate1
    !#    
    integer,intent(in)    :: isan_inc
    !#
    integer,intent(in)    :: nn
    !#
    character(len=*),intent(out) :: geosFileName
    !#
    character(len=*),intent(out) :: invFileName
    !#
    !Local variables
    character(len=2) :: ctime,cdate,cmonth
    character(len=4) :: cyear!,cntime


    !Code
    print *,'In: nn,isan_inc,iproc_date,iyear1,imonth1,idate1,itime1= ',nn,isan_inc,iproc_date,iyear1,imonth1,idate1,itime1
    ! GEOS.fp.fcst.inst3_3d_asm_Cp.20200213_00+20200214_0300.V01.nc4
    !write(cntime,fmt='(I4.4)') int(itime1)+(nn-1)*isan_inc
    write(ctime,fmt='(I2.2)') int(itime1/100)
    write(cyear,fmt='(I4.4)') iyear1
    write(cmonth,fmt='(I2.2)') imonth1
    write(cdate,fmt='(I2.2)') idate1
    geosFileName=trim(prefix)//cyear//cmonth//cdate//'_'//ctime//'+'//iproc_date(1:8)//'_'//iproc_date(9:12)//'.V01.nc4'
    invFileName=trim(prefix)//ctime//'.inv'
    
end subroutine makeGeosfileName
