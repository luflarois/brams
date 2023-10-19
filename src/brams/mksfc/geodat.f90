!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine geodat(n2, n3, datr, hfn, ofn, vt2da, vt2db, ngr, vnam)

  use mem_grid, only: &
       ngrid, nnxp, nnyp, deltaxn, deltayn, xtn, ytn, platn, plonn

  use io_params, only: &
       itopsflg, toptenh, toptwvl, iz0flg

  implicit none
  
  include "files.h"

  ! Arguments:
  integer, intent(IN)            :: n2, n3
  real, intent(INOUT)            :: vt2da(*), vt2db(*), datr(n2,n3)
  character(len=*), intent(IN) :: hfn, ofn
  character(len=03), intent(IN)  :: vnam
  ! Local Variables:
  integer :: ngr

  
  character(len=f_name_length) :: title
  integer :: lb, iblksizo, no, isbego, iwbego, iodim, mof, niq, njq, np
  real    :: offlat, offlon, deltallo, deltaxq, deltayq, deltaxp, deltayp, erad
  integer :: ierr
  real, allocatable :: dato(:)

  LB = len_trim(HFN)
  if (LB<=0) then
     print *, '==================================================='
     print *, '|  Problem in GEODAT, Input data prefix incorrect !'
     print *, '|  Grid :', ngrid
     print *, '|  File prefix:', trim(HFN)
     print *, '==================================================='
     call fatal_error('GEODAT-file')
  endif

  if ( (vnam(1:2)=='TO' .or. vnam(1:2)=='ZO') .and.  &
       (ITOPSFLG(NGR)==1 .and. TOPTENH(NGR)>1.)) then
     print *, '==================================================='
     print *, '|  Problem in GEODAT, Silhouette weight too high !'
     print *, '|  Grid :', NGR
     print *, '|  Weight (range TOPTENH=0-1):', TOPTENH(NGR)
     print *, '==================================================='
     call fatal_error('GEODAT')
  endif

  !     configure grid specs for raw data to rams grid (R) transfer

  !     raw data grid (O)
  if (vnam=='TOD' .or. vnam=='ZOD') then  !using dted data
     iblksizo =   5 ! 5 degree squares
     no       = 600 ! 30" data over a 5 degree square
     isbego   = -90
     iwbego   =   0
     offlat   =   0.
     offlon   =   0.
     DELTALLO = FLOAT(IBLKSIZO)/FLOAT(NO)
  else
     TITLE    = HFN(1:LB)//'HEADER'
     LB       = len_trim(TITLE)
     call rams_f_open(29, title(1:lb), 'FORMATTED', 'OLD', 'READ', 0)
     read (29,*) IBLKSIZO, NO, ISBEGO, IWBEGO, offlat, offlon
     close (29)
     DELTALLO = FLOAT(IBLKSIZO)/FLOAT(NO-1)
  endif

  iodim = max(100000, 4*no*no)
  MOF   = IODIM/(NO*NO)

  allocate(dato(iodim+mof+mof), STAT=ierr)
  if (ierr/=0) call fatal_error("Error allocating dato")

  !     temp grid (Q) - smoothing only applied to topo
  if (vnam(1:2)=='TO') then
     DELTAXQ = 0.5*abs(TOPTWVL(NGR))*DELTAXN(NGR)
     DELTAYQ = 0.5*abs(TOPTWVL(NGR))*DELTAYN(NGR)
  else
     DELTAXQ = DELTAXN(NGR)
     DELTAYQ = DELTAYN(NGR)
  endif
  NIQ = int(FLOAT(NNXP(NGR)-1)*DELTAXN(NGR)/DELTAXQ) + 4
  NJQ=int(FLOAT(NNYP(NGR)-1)*DELTAYN(NGR)/DELTAYQ)+4

  !     interpollated raw data grid (P)
  NP      = min(10, max(1,int(DELTAXQ/(DELTALLO*111000.))))
  DELTAXP = DELTAXQ/FLOAT(NP)
  DELTAYP = DELTAYQ/FLOAT(NP)

  call SFCOPQR(NO, MOF, NP, NIQ, NJQ, N2, N3, XTN(1,NGR), YTN(1,NGR), &
       platn(ngr), plonn(ngr),                                        &
       ERAD, DELTALLO, DELTAXP, DELTAYP, DELTAXQ, DELTAYQ, IBLKSIZO,  &
       ISBEGO, IWBEGO, DATO(1), VT2DA, VT2DB, DATR,                   &
       OFN, offlat, offlon, VNAM, NGR, itopsflg(ngr), iz0flg(ngr))

  deallocate(dato)

end subroutine geodat

!     ******************************************************************

subroutine sfcopqr(no, mof, np, niq, njq, n2, n3, xt, yt, platn, plonn, &
     erad, deltallo, deltaxp, deltayp, deltaxq, deltayq, iblksizo,      &
     isbego, iwbego, dato, datp, datq, datr,                            &
     ofn, offlat, offlon, vnam, ngr, itopsflg, iz0flg)

  use teb_spm_start, only: TEB_SPM

  implicit none
  ! Arguments:
  integer, intent(IN)            :: no, mof, np, niq, njq, n2, n3
  real, intent(IN)               :: xt(n2), yt(n3)
  real, intent(IN)               :: platn, plonn, erad, deltallo, &
       deltaxp, deltayp, deltaxq, deltayq
  integer, intent(IN)            :: iblksizo, isbego, iwbego
  real, intent(INOUT)            :: dato(no,no,mof), datp(np,np), datq(niq,njq), &
       datr(n2,n3)
  character(len=*), intent(IN) :: ofn
  real, intent(IN)               :: offlat, offlon
  character(len=003), intent(IN) :: vnam
  integer, intent(IN)            :: ngr, itopsflg, iz0flg
  ! Local Variables:

  include "files.h"
  
  character(len=f_name_length) :: title3
  character(len=003) :: title1
  character(len=004) :: title2
  logical            :: l1 
  integer, parameter :: maxmiss=1000
  character(len=256) :: fnmiss(maxmiss)
  real, allocatable  :: sdq(:,:), shaq(:,:), sdr(:,:), datre(:,:)
  real, allocatable  :: iso(:), iwo(:)
  real :: xcentr, ycentr, glatp, glonp, rio_full, rjo_full, &
       xq, yq, xp, yp, wio1, wio2, wjo1, wjo2, sha, rha, rh2, sh, rh, &
       xq1, yq1, xr, yr, rval, diff, difflcl
  integer :: nmiss, nono, nofr, iof, iq, jq, ip, jp, iwoc, isoc, &
       io1, io2, jo1, jo2, lb, nn, isocpt, isocpo, iwocpo, iwocph, iwocpt, &
       io_full, jo_full, iofr, jofr, ir, jr, is, js, i, j
  integer            :: ierr

  allocate (sdq(niq,njq), STAT=ierr)
  if (ierr/=0) call fatal_error("Error allocating sdq")
  allocate (shaq(niq,njq), STAT=ierr)
  if (ierr/=0) call fatal_error("Error allocating shaq")
  allocate (sdr(n2,n3), STAT=ierr)
  if (ierr/=0) call fatal_error("Error allocating sdr")
  allocate (datre(n2,n3), STAT=ierr)
  if (ierr/=0) call fatal_error("Error allocating datre")
  allocate (iso(mof), STAT=ierr)
  if (ierr/=0) call fatal_error("Error allocating iso")
  allocate (iwo(mof), STAT=ierr)
  if (ierr/=0) call fatal_error("Error allocating iwo")
  print *, '###'
  nmiss  = 0

  nono   = no*no
  XCENTR = 0.5*(XT(1)+XT(N2))
  YCENTR = 0.5*(YT(1)+YT(N3))
  NOFR   = 0

  do IOF=1,MOF
     ISO(IOF) = 0
     IWO(IOF) = 0
  enddo
  print *, 'Inside sfcopqr ',vnam,NJQ,NIQ,NP
  do JQ=1,NJQ
     do IQ=1,NIQ
        XQ = (FLOAT(IQ)-0.5*FLOAT(NIQ+1))*DELTAXQ + XCENTR
        YQ = (FLOAT(JQ)-0.5*FLOAT(NJQ+1))*DELTAYQ + YCENTR
        do JP=1,NP
           do IP=1,NP
              XP = XQ + (FLOAT(IP)-0.5*FLOAT(NP+1))*DELTAXP
              YP = YQ + (FLOAT(JP)-0.5*FLOAT(NP+1))*DELTAYP

              call xy_ll(GLATP, GLONP, platn, plonn, xp, yp)

              glatp = max(-89.9999, min(89.9999, (glatp - offlat)))

!--(DMK-CCATT-INI)-----------------------------------------------------
            !srf-rams60 mod
            glonp = max(-179.999,min(179.999,glonp - offlon))
!--(DMK-CCATT-ORI)-----------------------------------------------------
!            glonp = glonp - offlon
!--(DMK-CCATT-END)-----------------------------------------------------

              if (glonp>= 180.) glonp = glonp - 360.
              if (glonp<=-180.) glonp = glonp + 360.

              rio_full = (glonp - float(iwbego))/deltallo
              rjo_full = (glatp - float(isbego))/deltallo

              io_full  = int(rio_full)
              jo_full  = int(rjo_full)

              iwoc     = (io_full/(no-1))*iblksizo + iwbego
              isoc     = (jo_full/(no-1))*iblksizo + isbego

              wio2     = rio_full - float(io_full)
              wjo2     = rjo_full - float(jo_full)

              wio1     = 1. - wio2
              wjo1     = 1. - wjo2

              io1      = mod(io_full,no-1) + 1
              jo1      = mod(jo_full,no-1) + 1

              io2      = io1 + 1
              jo2      = jo1 + 1

              do IOFR=1,NOFR
                 JOFR = IOFR
                 if (ISO(IOFR)==ISOC .and. IWO(IOFR)==IWOC) GO TO 10
              enddo

              !                  not using dted data
              if (vnam/='TOD' .and. vnam/='ZOD') then
                 ISOCPT = abs(ISOC)/10
                 ISOCPO = abs(ISOC) - ISOCPT*10
                 IWOCPH = abs(IWOC)/100
                 IWOCPT = (abs(IWOC)-IWOCPH*100)/10
                 IWOCPO = abs(IWOC) - IWOCPH*100 - IWOCPT*10
                 if (ISOC>=0) then
                    write(TITLE1,'(2I1,A1)') ISOCPT, ISOCPO, 'N'
                 else
                    write(TITLE1,'(2I1,A1)') ISOCPT, ISOCPO, 'S'
                 endif
                 if (IWOC>=0) then
                    write(TITLE2,'(3I1,A1)') IWOCPH, IWOCPT, IWOCPO, 'E'
                 else
                    write(TITLE2,'(3I1,A1)') IWOCPH, IWOCPT, IWOCPO, 'W'
                 endif
                 LB     = len_trim(OFN)
                 TITLE3 = OFN(1:LB)//TITLE1//TITLE2
!!$                 print *, "DEBUG-ALF:LB,trim(OFN)=",LB,trim(OFN) 
                 LB     = len_trim(TITLE3)
!!$                 print *, "DEBUG-ALF:LB,trim(TITLE3)=",LB,trim(TITLE3) 
                 inquire(FILE=TITLE3(1:LB), EXIST=L1)
                 if(L1) print *, "DEBUG-ALF:EXIST?=", TITLE3(1:LB),L1
                 if (.not.L1) then
                    do nn=1,nmiss
                       if (TITLE3(1:LB)==fnmiss(nn)) goto 302
                    enddo
                    nmiss         = nmiss + 1
                    fnmiss(nmiss) = TITLE3(1:LB)
302                 continue
                    DATP(IP,JP)   = 0.
                    goto 20
                 endif
              endif
              
              if (NOFR>=MOF) then
                 do IOF=1,MOF
                    ISO(IOF) = 0
                    IWO(IOF) = 0
                 enddo
                 NOFR = 0
              endif
              NOFR = NOFR + 1
              JOFR = NOFR
 
              !                 using dted data
              if (vnam=='TOD' .or. vnam=='ZOD') then
                 call dted(no, ofn, isoc, iwoc, dato(1,1,nofr))
              else
                 call rams_f_open(29, TITLE3(1:LB), 'FORMATTED', 'OLD', 'READ', 0)
!!$                 call rams_f_open(29, trim(TITLE3), 'FORMATTED', 'OLD', &
!!$                      'READ', 0)
                 call VFIREC(29, DATO(1,1,NOFR), NONO, 'LIN')
                 close(29)
              endif

              ISO(NOFR) = ISOC
              IWO(NOFR) = IWOC

10            continue

!!$              ! DEBUG-ALF
!!$              if (io1<1 .or. jo1<1 .or. jofr<1 .or. io2<1 .or. jo2<1) then
!!$                 print *, "DEBUG-ALF:ip,jp,io1,io2,jo1,jo2,jofr=", &
!!$                      ip,jp,io1,io2,jo1,jo2,jofr
!!$                 print *, "DEBUG-ALF:GLONP,rio_full,io_full,iwbego,deltallo,no=", &
!!$                      GLONP,rio_full,io_full,iwbego,deltallo,no
!!$                 call flush(6)
!!$              endif
!!$              !

              datp(ip,jp) = wio1*(wjo1*dato(io1,jo1,jofr) + &
                   wjo2*dato(io1,jo2,jofr))               + &
                   wio2*(wjo1*dato(io2,jo1,jofr)          + &
                   wjo2*dato(io2,jo2,jofr))

20            continue
           enddo
        enddo

        !           std dev for envelope orog and topo based zo schemes
        SHA = 0.
        RHA = 0.
        RH2 = 0.
        do JP=1,NP
           SH = 0.
           RH = 0.
           do IP=1,NP
              SH  = max(SH, DATP(IP,JP))
              RH  = RH  + DATP(IP,JP)
              RH2 = RH2 + DATP(IP,JP)**2
           enddo
           SHA = SHA + SH/(2.*FLOAT(NP))
           RHA = RHA + RH
        enddo
        DATQ(IQ,JQ) = RHA/FLOAT(NP*NP)
        SDQ(IQ,JQ)  = sqrt(max(0., RH2/NP**2-DATQ(IQ,JQ)**2))
        do IP=1,NP
           SH = 0.
           do JP=1,NP
              SH = max(SH, DATP(IP,JP))
           enddo
           SHA = SHA + SH/(2.*FLOAT(NP))
        enddo
        SHAQ(IQ,JQ) = SHA

     enddo
     !         print*,'finished sfcopqr row jq = ',jq
  enddo

  !     envelope and zo schemes

  if ( (vnam=='TOP' .or. vnam=='TOD') .and. &
       (ITOPSFLG==2 .or. ITOPSFLG==3) .and. &
       (NP*NP)<8)                           &
       print*,'Warning - trying to calc a std dev for: ', NP*NP, ' points'
  if ( (vnam=='ZOT' .or. vnam=='ZOD') .and. &
       IZ0FLG==1 .and. (NP*NP)<8)           &
       print*,'Warning - trying to calc a std dev for: ', NP*NP, ' points'

  if (vnam=='TOP' .or. vnam=='TOD')  &
       call TOPOQ(NIQ, NJQ, DELTAXQ, DELTAYQ, DATQ, SDQ, SHAQ, &
       DATRE, NGR, N2, N3)
  if (vnam=='ZOT' .or. vnam=='ZOD')  &
       call ZOQ(NIQ, NJQ, DATQ, SDQ, NGR)

  XQ1 = (1.-0.5*FLOAT(NIQ+1))*DELTAXQ + XCENTR
  YQ1 = (1.-0.5*FLOAT(NJQ+1))*DELTAYQ + YCENTR
  do JR=1,N3
     do IR=1,N2
        XR = (XT(IR)-XQ1)/DELTAXQ + 1.
        YR = (YT(JR)-YQ1)/DELTAYQ + 1.
        call GDTOST(DATQ, NIQ, NJQ, XR, YR, RVAL)

        !           envelope orog and zo schemes

        if (vnam=='ZOT' .or. vnam=='ZOD') then
           DATR(IR,JR) = max(0., RVAL)
           !               print*,'z0r',IR,JR,DATR(IR,JR)
        else
           if (TEB_SPM==1) then
              if (vnam/='FUS') then
                 DATR(IR,JR) = max(0., RVAL)
              else
                 datr(ir,jr) = rval
              endif
           else
              DATR(IR,JR) = max(0., RVAL)
           endif
        endif

        call GDTOST(SDQ, NIQ, NJQ, XR, YR, RVAL)
        SDR(IR,JR) = max(0., RVAL)
     enddo
  enddo

  if (nmiss>0) then
     print *, '-----------------------------------------------------'
     print *, 'Input physiographical data file processing: (sfcopqr)'
     print *, '-----------------------------------------------------'
     print *, '  Input data blocks not found (data assumed to be zero):'
     do nn=1,nmiss
        print *, nn, trim(fnmiss(nn))
     enddo
     print *, '-----------------------------------------------------'
  endif

  !     check to find the largest change in topo height

  if (vnam=='TOP' .or. vnam=='TOD') then
     diff     =    0.
     difflcl  =    0.
     is       = -999
     js       = -999
     do j=2,n3-1
        do i=2,n2-1
           !print*,'1=',difflcl,max(difflcl,abs(datr(i,j)-datr(i-1,j)))
           difflcl = max(difflcl, abs(datr(i,j)-datr(i-1,j)))
           !print*,'2=',difflcl,max(difflcl,abs(datr(i,j)-datr(i+1,j)))
           difflcl = max(difflcl, abs(datr(i,j)-datr(i+1,j)))
           !print*,'3=',difflcl,max(difflcl,abs(datr(i,j)-datr(i,j-1)))
           difflcl = max(difflcl, abs(datr(i,j)-datr(i,j-1)))
           !print*,'4=',difflcl,max(difflcl,abs(datr(i,j)-datr(i,j+1)))
           difflcl = max(difflcl, abs(datr(i,j)-datr(i,j+1)))
           if (abs(diff-difflcl)>1.) then
              is = i
              js = j
           endif
           diff = max(diff, difflcl)
        enddo
     enddo
     write (6,100) ' Max d(topo) on grid @i,j=', ngr, is, js, diff
100  format(a,3i4,f8.1)
  endif

  deallocate(SDQ, SHAQ, SDR, DATRE)
  deallocate(ISO, IWO)

end subroutine sfcopqr

!**********************************************************************

subroutine topoq(niq, njq, deltaxq, deltayq, datq, sdq, shaq, datre,  &
     ngr, n2, n3)

  use io_params, only: &
       itopsflg, toptenh, toptwvl

  implicit none
  ! Arguments:
  integer, intent(IN) :: niq, njq, n2, n3, ngr
  real, intent(IN)    :: deltaxq, deltayq
  real, intent(INOUT) :: datq(niq,njq)
  real, intent(IN)    :: sdq(niq,njq), shaq(niq,njq)
  real, intent(INOUT) :: datre(n2,n3)
  ! Local Variables:
  integer :: iq, jq, jmin, imin, ire, jre, imax, jmax
  real    :: rad, count, total, remax, remin, average

  !     orographic schemes

  if (ITOPSFLG(ngr)<0) then                         ! No orography
     do jq=1,njq
        do iq=1,niq
           datq(iq,jq) = 0.
           !            print*,'None',iq,jq,datq(iq,jq)
        enddo
     enddo
     print *,'No orography'
  elseif (ITOPSFLG(ngr)<0) then                     ! Average
     print *, 'No orography enhancement applied'
  elseif (ITOPSFLG(ngr)==1) then                    ! Silhouette
     do jq=1,njq
        do iq=1,niq
           datq(iq,jq) = SHAQ(IQ,JQ)*toptenh(ngr) + &
                DATQ(IQ,JQ)*(1.-toptenh(ngr))
           !                  print*,'Silhouette',iq,jq,datq(iq,jq)
        enddo
     enddo
     print *, 'Silhouette Orography applied with'
     print *, 'weighting = ', toptenh(ngr)
  elseif (ITOPSFLG(ngr)==2) then                    ! Envelope
     do jq=1,njq
        do iq=1,niq
           datq(iq,jq) = datq(iq,jq) + toptenh(ngr)*sdq(iq,jq)
           !                  print*,'EO',iq,jq,datq(iq,jq)
        enddo
     enddo
     print *, 'Envelope Orography applied with'
     print *, 'enhancement = ', toptenh(ngr), ' x std dev'
  else if (ITOPSFLG(ngr)>=3) then                   ! Reflected Envelope

     !        the radius we want to search for the current pts relative
     !        height should correspond well to half the filtering wavelength
     !        used on the topo (toptwvl)

     Rad = abs(toptwvl(ngr))/2
     do jq=1,njq
        do iq=1,niq
           datre(iq,jq) = datq(iq,jq)
        enddo
     enddo
     do jq=1,njq
        do iq=1,niq
           count = 0.
           total = 0.
           remax = datre(iq,jq)
           remin = datre(iq,jq)
           jmin  = jq - nint(Rad)
           imin  = iq - nint(Rad)
           jmax  = jq + nint(Rad)
           imax  = iq + nint(Rad)
           do jre=max(1,jmin),min(njq,jmax)
              do ire=max(1,imin),min(niq,imax)
                 if ((float(iq-ire)**2 + float(jq-jre)**2)<=Rad**2) then
                    count = count + 1.
                    total = total + datre(ire,jre)
                    remax = max(remax, datre(ire,jre))
                    remin = min(remin, datre(ire,jre))
                 endif
              enddo
           enddo
           average = total/count
           if (remax/=remin) &
                datq(iq,jq) = datre(iq,jq) + (datre(iq,jq)-average)/  &
                ((remax-remin)/2)*toptenh(ngr)*sdq(iq,jq)
           !               print*,'REO',iq,jq,datre(iq,jq),sdq(iq,jq),datq(iq,jq)
           !               print*,'avg,n',average,count,remax,remin
        enddo
     enddo
     print *, 'Reflected Envelope Orography applied with'
     print *, 'enhancement = ', toptenh(ngr), ' x std dev'
     print *, 'and search radius (grid points) = ', Rad
  endif

end subroutine topoq

!**********************************************************************

subroutine ZOQ(NIQ, NJQ, DATQ, SDQ, NGR)

  use io_params, only: &
       itopsflg, iz0flg, z0fact, z0max

  use mem_leaf, only: &
       zrough

  implicit none
  ! Arguments:
  integer, intent(IN) :: NIQ, NJQ, NGR
  real, intent(INOUT) :: DATQ(NIQ,NJQ)
  real, intent(IN)    :: SDQ(NIQ,NJQ)
  ! Local Variables
  integer :: iq, jq

  !     topo base roughness length.

  do jq=1,njq
     do iq=1,niq
        if (ITOPSFLG(ngr)<0) then  ! No orography
           datq(iq,jq) = zrough
        elseif (iz0flg(ngr)==1) then
           datq(iq,jq) = min(z0fact*sdq(iq,jq), z0max(NGR))
        else
           datq(iq,jq) = zrough
        endif
        !            print*,'z0',iq,jq,datq(iq,jq)
     enddo
  enddo
  if (ITOPSFLG(ngr)<0) then  ! No orography
     print *, 'No orography'
!!$  else
!!$     print *, 'Subgrid terrain roughness applied with'
!!$     print *, 'factor  = ', z0fact, ' x std dev'
!!$     print *, 'maximum = ', z0max(NGR)
  endif

end subroutine ZOQ

!**********************************************************************

subroutine dted(no, pathname, lat, lon, dato)

  !     Let's try and bypass all the bookeeping and just read the file
  !     Note that the latitude bands are 5 degrees less than those
  !     specified in the original code from Sarma since we are not
  !     considering the max latitude but rather the start latitude.

  implicit none
  ! Arguments:
  integer, intent(IN)            :: no, lat, lon
  real, intent(OUT)              :: dato(no,no)
  character(len=256), intent(IN) :: pathname
  ! Local Variables:
!!$  character(len=80) :: fname
  integer :: ifact, notfnd

  ifact = 6
  if (lat>=0) then
     if (lat<=75)  ifact = 4
     if (lat<=70)  ifact = 3
     if (lat<=65)  ifact = 2
     if (lat<=45)  ifact = 1
  else
     ifact = 1
     if (lat<=-45) ifact = 2
     if (lat<=-65) ifact = 3
     if (lat<=-70) ifact = 4
     if (lat<=-75) ifact = 6
  end if
  call dtedint(no, ifact, lon, lat, notfnd, pathname, dato)

end subroutine dted

! ----------------------------------------------------------------------

subroutine dtedint(no, iwres, lon, lat, notfnd, pathname, dato)

  !     This routine (call by readdtdt) determines the filename
  !     of the dted terrain file, based on the latitude and longitude,
  !     and calls the c routine that reads the file.
  !     uses shell script to read a data set that has been compressed
  !     Most recent development Paul Boris.
  !     new input data holding arrays:

  implicit none
  ! Arguments:
  integer, intent(IN)            :: no
  integer, intent(INOUT)         :: iwres
  integer, intent(IN)            :: lon, lat
  integer, intent(INOUT)         :: notfnd
  real, intent(OUT)              :: dato(no,no)
  character(len=256), intent(IN) :: pathname
  ! Local Variables:
  real :: readin2(360000)
  character(len=256) :: newname1
  character(len=005) :: fmtstr
  character(len=003) :: degns,degew
  character(len=012) :: dtedfile
  character(len=004) :: subdir
  character(len=256) ::  extdted !, callname
  real :: dtedwk(600,600)
  integer :: no_blanks, i, num_chrs, ifile, lc, isave, ndtx, ndty, &
       ik, jk, ibtre, j, id
  real :: rvaln, wt
  real, external :: readdted1

  no_blanks = 1
  do i = 1, len(pathname)
     if (pathname(i:i)/=' ') no_blanks = i
  enddo

  num_chrs = no_blanks + 17
  write (fmtstr,'(a2,i2,a1)') '(a',num_chrs,')'

  write (degew,'(i3)') abs(lon)
  if (degew(1:1)==' ') degew(1:1) = '0'
  if (degew(2:2)==' ') degew(2:2) = '0'
  write (degns,'(i3)') abs(lat)
  if (degns(1:1)==' ') degns(1:1) = '0'
  if (degns(2:2)==' ') degns(2:2) = '0'
  if (lon<0) then
     dtedfile(5:5) = 'w'
  else
     dtedfile(5:5) = 'e'
  endif
  dtedfile(6:8) = degew(1:3)
  if (lat<0) then
     dtedfile(1:1) = 's'
  else
     dtedfile(1:1) = 'n'
  endif
  dtedfile(2:4) = degns(1:3)

  !     Change from filenames .mvk.gz <---> .mgz
  !      dtedfile(9:12) = '.mgz'
  dtedfile(9:12) = '.mvk'

  !     open the dted data file and read the terrain data:

  !      print*, 'OPENING DTED DATA FILE: ',dtedfile
  subdir   = '/'//dtedfile(1:1)//dtedfile(5:5)//'/'
  newname1 = pathname(1:no_blanks)//subdir

  ifile = 42
  open(ifile, file='namefils.out')
  write (ifile, fmt='(a12)') dtedfile
  close (ifile)

  !     Call all the decompress stuff from this Fortran code.
  lc = len_trim(newname1)

  !     Change from filenames .mvk.gz <---> .mgz
  !      extdted=newname1(1:lc)//dtedfile
  extdted = newname1(1:lc)//dtedfile//'.gz'

  lc = len_trim(extdted)
  write (newname1,1002) extdted(1:lc)
1002 format('cp ', a, ' tmp.gz')
  print *, newname1(1:index(newname1,'tmp.gz')+5)
  call system(newname1)
  call system('chmod ugo+rw tmp.gz')
  call system('gzip -d tmp.gz')
  write (newname1,1003) dtedfile
1003 format('mv tmp ',a)
  call system(newname1)

!!$  call azero(no*no,dato)
  dato = 0.

  isave = iwres
  rvaln = READDTED1(iwres, readin2)
  if (iwres==-100000) then
     notfnd = 1
     return
  endif
  iwres = isave

  ndtx = 600/iwres
  ndty = 600

  do jk=1,ndty
     do ik=1,ndtx
        ibtre         = (jk-1)*ndtx + ik
        dtedwk(ik,jk) = readin2(ibtre)
     enddo
  enddo

  write (newname1,1004) dtedfile
1004 format('rm -f ', a)
  call system(newname1)

  do j=1,no
     do i=1,no-iwres+1
        id        = (i+iwres-1)/iwres
        wt        = mod(real(i+iwres-1), real(iwres))/real(iwres)
        dato(i,j) = (1.-wt)*dtedwk(id,j) + wt*dtedwk(id+1,j)
     enddo
  enddo

end subroutine dtedint

!**************************************************************************

subroutine geodat_var(n2, n3, datr, hfn, ofn, vt2da, vt2db, ngr, vnam)

  use mem_grid, only: &
       ngrid, xmn, ymn, XTN, YTN, platn, plonn,nnxp, nnyp, deltaxn, deltayn

  use io_params, only: &
       ITOPSFLG, TOPTENH, iz0flg, toptwvl

  implicit none
  include "files.h"

  ! Arguments:
  integer, intent(IN)            :: n2, n3, ngr
  real, intent(INOUT)            :: vt2da(*), vt2db(*), datr(n2,n3)
  ! Debug by Lamosa
  ! change len to *
  ! old:
  !character(len=f_name_length), intent(IN) :: hfn, ofn
  ! new:
  character(len=*), intent(IN) :: hfn, ofn

  character(len=03), intent(IN)  :: vnam

  ! Local Variables:
  !integer :: ngr
  character(len=f_name_length) :: title
  integer :: lb, iblksizo, no, isbego, iwbego, iodim, mof, niq, njq, np, i, j
  real :: offlat, offlon, deltallo, deltaxq, deltayq, deltaxp, deltayp, erad
  real, allocatable:: dato(:)
  !-srf
  real    :: TOPTWVL_2d(n2,n3)
  integer :: ierr

  print *, '====================================================='
  if(VNAM(1:2)=='TO') then
     print *, 'starting topography on grid:', NGR
  elseif(VNAM(1:2)=='ZO') then
     print *, 'starting surface roughness on grid:', NGR
  else
     print *, 'starting '//vnam//' data on grid:', NGR
  endif

  LB = len_trim(HFN)
  if (LB<=0) then
     print *, '==================================================='
     print *, '|  Problem in GEODAT, Input data prefix incorrect !'
     print *, '|  Grid :', ngrid
     print *, '|  File prefix:', trim(HFN)
     print *, '==================================================='
     call fatal_error('GEODAT-file')
  endif

  if ((vnam(1:2)=='TO' .or. vnam(1:2)=='ZO') .and. &
       (ITOPSFLG(NGR)==1 .and. TOPTENH(NGR)>1.)) then
     print *, '==================================================='
     print *, '|  Problem in GEODAT, Silhouette weight too high !'
     print *, '|  Grid :', NGR
     print *, '|  Weight (range TOPTENH=0-1):', TOPTENH(NGR)
     print *, '==================================================='
     call fatal_error('GEODAT')
  endif

  !     configure grid specs for raw data to rams grid (R) transfer

  !     raw data grid (O)
  if (vnam=='TOD' .or. vnam=='ZOD') then  !using dted data
     iblksizo =   5 ! 5 degree squares
     no       = 600     ! 30" data over a 5 degree square
     isbego   = -90
     iwbego   =   0
     offlat   =   0.
     offlon   =   0.
     DELTALLO = FLOAT(IBLKSIZO)/FLOAT(NO)
  else
     TITLE = HFN(1:LB)//'HEADER'
     LB    = len_trim(TITLE)
     call rams_f_open(29, title(1:lb), 'FORMATTED', 'OLD', 'READ', 0)
     READ(29,*) IBLKSIZO, NO, ISBEGO, IWBEGO, offlat, offlon
     CLOSE(29)
     DELTALLO = FLOAT(IBLKSIZO)/FLOAT(NO-1)
  endif

  iodim = max(100000, 4*no*no)
  MOF   = IODIM/(NO*NO)

  allocate(dato(iodim+mof+mof), STAT=ierr)
  if (ierr/=0) call fatal_error("Error allocating dato")
  
  
  !- 30jun2012srf -opt
  ! - part 1  - homegeneous smoothing
     !     temp grid (Q) - smoothing only applied to topo
  if (vnam(1:2)=='TO') then
     DELTAXQ = 0.5*abs(TOPTWVL(NGR))*DELTAXN(NGR)
     DELTAYQ = 0.5*abs(TOPTWVL(NGR))*DELTAYN(NGR)
  else
     DELTAXQ = DELTAXN(NGR)
     DELTAYQ = DELTAYN(NGR)
  endif
  NIQ = int(FLOAT(NNXP(NGR)-1)*DELTAXN(NGR)/DELTAXQ) + 4
  NJQ=int(FLOAT(NNYP(NGR)-1)*DELTAYN(NGR)/DELTAYQ)+4

  !     interpollated raw data grid (P)
  NP      = min(10, max(1,int(DELTAXQ/(DELTALLO*111000.))))
  DELTAXP = DELTAXQ/FLOAT(NP)
  DELTAYP = DELTAYQ/FLOAT(NP)

  call SFCOPQR(NO, MOF, NP, NIQ, NJQ, N2, N3, XTN(1,NGR), YTN(1,NGR), &
       platn(ngr), plonn(ngr),                                        &
       ERAD, DELTALLO, DELTAXP, DELTAYP, DELTAXQ, DELTAYQ, IBLKSIZO,  &
       ISBEGO, IWBEGO, DATO(1), VT2DA, VT2DB, DATR,                   &
       OFN, offlat, offlon, VNAM, NGR, itopsflg(ngr), iz0flg(ngr))

  ! - part 2  - variable smoothing

  !-srf:  get TOPTWVL with grid dependence
  call get_TOPTWVL(ngr, n2, n3, TOPTWVL_2d)
  !-srf end

  do j=2,n3-1
     do i=2,n2-1

        if(TOPTWVL_2d(i,j) == abs(TOPTWVL(NGR)) ) cycle

        !     temp grid (Q) - smoothing only applied to topo
        if (vnam(1:2)=='TO') then

           !-srf : orig with TOPTWVL cte qq i,j
           !DELTAXQ=0.5*TOPTWVL(NGR) * ( xmn(i,ngr)-xmn(i-1,ngr) )!*DELTAXN(NGR)
           !DELTAYQ=0.5*TOPTWVL(NGR) * ( ymn(j,ngr)-ymn(j-1,ngr) ) 
           
           !-srf : TOPTWVL f(i,j)
           DELTAXQ = 0.5*TOPTWVL_2d(i,j)*(xmn(i,ngr) - xmn(i-1,ngr)) !*DELTAXN(NGR)
           DELTAYQ = 0.5*TOPTWVL_2d(i,j)*(ymn(j,ngr) - ymn(j-1,ngr)) 
           !print*,'geodat1=',TOPTWVL_2d(i,j), DELTAYQ
        else
           !srf -orig
           !DELTAXQ=DELTAXN(NGR)
           !DELTAYQ=DELTAYN(NGR)
           DELTAXQ = (xmn(i,ngr) - xmn(i-1,ngr))
           DELTAYQ = (ymn(j,ngr) - ymn(j-1,ngr))
        endif

        !-versao1
        !NIQ=INT(FLOAT(NNXP(NGR)-1)*( xmn(i,ngr)-xmn(i-1,ngr) )/DELTAXQ)+4
        !NJQ=INT(FLOAT(NNYP(NGR)-1)*( ymn(j,ngr)-ymn(j-1,ngr) )/DELTAYQ)+4
        !-versao opt
        NIQ = INT((xmn(i,ngr) - xmn(i-1,ngr))/DELTAXQ) + 4
        NJQ = INT((ymn(j,ngr) - ymn(j-1,ngr))/DELTAYQ) + 4

        !     interpollated raw data grid (P)
        NP      = MIN(10, MAX(1, INT(DELTAXQ/(DELTALLO*111000.))))
        DELTAXP = DELTAXQ/FLOAT(NP)
        DELTAYP = DELTAYQ/FLOAT(NP)

        CALL SFCOPQR_VAR(NO, MOF, NP, NIQ, NJQ, N2, N3, XTN(1,NGR), YTN(1,NGR), &
             platn(ngr), plonn(ngr),                                            &
             ERAD, DELTALLO, DELTAXP, DELTAYP, DELTAXQ, DELTAYQ, IBLKSIZO,      &
             ISBEGO, IWBEGO, DATO(1), VT2DA, VT2DB, DATR,                       &
             OFN, offlat, offlon, VNAM, NGR, itopsflg(ngr), iz0flg(ngr), i, j)

     enddo
  enddo

  deallocate(dato)

  !srf-
  DATR(1,:)  = DATR(2,:)
  DATR(N2,:) = DATR(N2-1,:)
  DATR(:,1)  = DATR(:,2)
  DATR(:,N3) = DATR(:,N3-1)

END subroutine geodat_var

!**********************************************************************

subroutine sfcopqr_var(no, mof, np, niq, njq, n2, n3, xt, yt, platn, plonn, &
     erad, deltallo, deltaxp, deltayp, deltaxq, deltayq, iblksizo,          &
     isbego, iwbego, dato, datp, datq, datr,                                &
     ofn, offlat, offlon, vnam, ngr, itopsflg, iz0flg, iii, jjj)

  implicit none
  ! Arguments:
  integer, intent(IN)            :: no, mof, np, niq, njq, n2, n3
  real, intent(IN)               :: xt(n2), yt(n3)
  real, intent(IN)               :: platn, plonn, erad, deltallo, &
       deltaxp, deltayp, deltaxq, deltayq
  integer, intent(IN)            :: iblksizo, isbego, iwbego
  real, intent(INOUT)            :: dato(no,no,mof), datp(np,np), datq(niq,njq), &
       datr(n2,n3)
  character(len=256), intent(IN) :: ofn
  real, intent(IN)               :: offlat, offlon
  character(len=03), intent(IN)  :: vnam
  integer, intent(IN)            :: ngr, itopsflg, iz0flg
  integer, intent(IN)            :: iii, jjj
  ! Local Variables:

  include "files.h"

  character(len=f_name_length) :: title3
  character(len=003) :: title1
  character(len=004) :: title2
  logical            :: l1
  integer, parameter :: maxmiss=1000
  character(len=256) :: fnmiss(maxmiss)
  real, allocatable  :: sdq(:,:), shaq(:,:), sdr(:,:), datre(:,:)
  real, allocatable  :: iso(:), iwo(:)
  real               :: xcentr, ycentr, glatp, glonp, rio_full, rjo_full,  &
       xq, yq, xp, yp, wio1, wio2, wjo1, wjo2, sha, rha, rh2, sh, rh, &
       xq1, yq1, xr, yr, rval, diff, difflcl
  integer            :: nmiss, nono, nofr, iof, iq, jq, ip, jp, iwoc, isoc, &
       io1, io2, jo1, jo2, &
       lb, nn, isocpt, isocpo, iwocpo, iwocph, iwocpt, io_full, jo_full, &
       iofr, jofr, ir, jr, is, js, i, j
  integer            :: ierr

  allocate (sdq(niq,njq), STAT=ierr)
  if (ierr/=0) call fatal_error("Error allocating sdq")
  allocate (shaq(niq,njq), STAT=ierr)
  if (ierr/=0) call fatal_error("Error allocating shaq")
  allocate (sdr(n2,n3), STAT=ierr)
  if (ierr/=0) call fatal_error("Error allocating sdr")
  allocate (datre(n2,n3), STAT=ierr)
  if (ierr/=0) call fatal_error("Error allocating datre")
  allocate (iso(mof), STAT=ierr)
  if (ierr/=0) call fatal_error("Error allocating iso")
  allocate (iwo(mof), STAT=ierr)
  if (ierr/=0) call fatal_error("Error allocating iwo")

  nmiss  = 0

  nono   = no*no
  !- para versao 1
  !XCENTR=0.5*(XT(1)+XT(N2))   
  !YCENTR=0.5*(YT(1)+YT(N3))   
  !- para versao OPT
  XCENTR = 0.5*(XT(III)+XT(III+1)) 
  YCENTR = 0.5*(YT(JJJ)+YT(JJJ+1)) 

  NOFR   = 0

  DO IOF=1,MOF
     ISO(IOF) = 0
     IWO(IOF) = 0
  ENDDO

  DO JQ=1,NJQ
     DO IQ=1,NIQ
        XQ = (FLOAT(IQ) - 0.5*FLOAT(NIQ+1))*DELTAXQ + XCENTR
        YQ = (FLOAT(JQ) - 0.5*FLOAT(NJQ+1))*DELTAYQ + YCENTR
        DO JP=1,NP
           DO IP=1,NP
              XP = XQ + (FLOAT(IP) - 0.5*FLOAT(NP+1))*DELTAXP
              YP = YQ + (FLOAT(JP) - 0.5*FLOAT(NP+1))*DELTAYP

              call xy_ll(GLATP, GLONP, platn, plonn, xp, yp)

              glatp = max(-89.9999, min(89.9999, glatp - offlat))
              !srf-rams60 mod
              !            glonp = glonp - offlon
              glonp = max(-179.999, min(179.999, glonp - offlon))
              !
              if (glonp>= 180.) glonp = glonp - 360.
              if (glonp<=-180.) glonp = glonp + 360.

              rio_full = (glonp - float(iwbego))/deltallo
              rjo_full = (glatp - float(isbego))/deltallo

              io_full  = int(rio_full)
              jo_full  = int(rjo_full)

              iwoc     = (io_full/(no-1))*iblksizo + iwbego
              isoc     = (jo_full/(no-1))*iblksizo + isbego

              wio2     = rio_full - float(io_full)
              wjo2     = rjo_full - float(jo_full)

              wio1     = 1. - wio2
              wjo1     = 1. - wjo2

              io1      = mod(io_full,no-1) + 1
              jo1      = mod(jo_full,no-1) + 1

              io2      = io1 + 1
              jo2      = jo1 + 1

              DO IOFR=1,NOFR
                 JOFR = IOFR
                 IF (ISO(IOFR)==ISOC .AND. IWO(IOFR)==IWOC) GO TO 10
              ENDDO

              !                  not using dted data
              if (vnam/='TOD' .and. vnam/='ZOD') then
                 ISOCPT = ABS(ISOC)/10
                 ISOCPO = ABS(ISOC) - ISOCPT*10
                 IWOCPH = ABS(IWOC)/100
                 IWOCPT = (ABS(IWOC) - IWOCPH*100)/10
                 IWOCPO = ABS(IWOC) - IWOCPH*100 - IWOCPT*10
                 IF (ISOC>=0) THEN
                    WRITE(TITLE1,'(2I1,A1)') ISOCPT, ISOCPO, 'N'
                 ELSE
                    WRITE(TITLE1,'(2I1,A1)') ISOCPT, ISOCPO, 'S'
                 ENDIF
                 IF (IWOC>=0) THEN
                    WRITE(TITLE2,'(3I1,A1)') IWOCPH, IWOCPT, IWOCPO, 'E'
                 ELSE
                    WRITE(TITLE2,'(3I1,A1)') IWOCPH, IWOCPT, IWOCPO, 'W'
                 ENDIF
                 LB     = len_trim(OFN)
                 TITLE3 = OFN(1:LB)//TITLE1//TITLE2
                 LB     = len_trim(TITLE3)
                 INQUIRE(FILE=TITLE3(1:LB), EXIST=L1)

                 IF (.NOT.L1) THEN
                    do nn=1,nmiss
                       if (TITLE3(1:LB)==fnmiss(nn)) goto 302
                    enddo
                    nmiss         = nmiss + 1
                    fnmiss(nmiss) = TITLE3(1:LB)
302                 continue
                    DATP(IP,JP)   = 0.
                    GOTO 20
                 ENDIF
              ENDIF

              IF (NOFR>=MOF) THEN
                 DO IOF=1,MOF
                    ISO(IOF) = 0
                    IWO(IOF) = 0
                 ENDDO
                 NOFR = 0
              ENDIF
              NOFR = NOFR + 1
              JOFR = NOFR

              !                 using dted data
              if (vnam=='TOD' .or. vnam=='ZOD') then
                 call dted(no, ofn, isoc, iwoc, dato(1,1,nofr))
              else
                 call rams_f_open  &
                      (29, TITLE3(1:LB), 'FORMATTED', 'OLD', 'READ', 0)
                 CALL VFIREC(29, DATO(1,1,NOFR), NONO, 'LIN')
                 CLOSE(29)
              endif

              ISO(NOFR) = ISOC
              IWO(NOFR) = IWOC

10            CONTINUE

!!$              ! DEBUG-ALF
!!$              if (io1<1 .or. jo1<1 .or. jofr<1 .or. io2<1 .or. jo2<1) then
!!$                 print *, "DEBUG-ALF2:ip,jp,io1,io2,jo1,jo2,jofr=", &
!!$                      ip,jp,io1,io2,jo1,jo2,jofr
!!$                 print *, "DEBUG-ALF2:GLONP,rio_full,io_full,iwbego,deltallo,no=", &
!!$                      GLONP,rio_full,io_full,iwbego,deltallo,no
!!$                 call flush(6)
!!$              endif
!!$              !

              datp(ip,jp) = wio1*(wjo1*dato(io1,jo1,jofr) + &
                   wjo2*dato(io1,jo2,jofr))               + &
                   wio2*(wjo1*dato(io2,jo1,jofr)          + &
                   wjo2*dato(io2,jo2,jofr))

20            CONTINUE
           ENDDO
        ENDDO

        !           std dev for envelope orog and topo based zo schemes
        SHA = 0.
        RHA = 0.
        RH2 = 0.
        DO JP=1,NP
           SH = 0.
           RH = 0.
           DO IP=1,NP
              SH  = MAX(SH,DATP(IP,JP))
              RH  = RH  + DATP(IP,JP)
              RH2 = RH2 + DATP(IP,JP)**2
           ENDDO
           SHA = SHA + SH/(2.*FLOAT(NP))
           RHA = RHA + RH
        ENDDO
        DATQ(IQ,JQ) = RHA/FLOAT(NP*NP)
        SDQ(IQ,JQ)  = SQRT(max(0., RH2/NP**2 - DATQ(IQ,JQ)**2))
        DO IP=1,NP
           SH = 0.
           DO JP=1,NP
              SH = MAX(SH, DATP(IP,JP))
           ENDDO
           SHA = SHA + SH/(2.*FLOAT(NP))
        ENDDO
        SHAQ(IQ,JQ) = SHA

     ENDDO
     !         print*,'finished sfcopqr row jq = ',jq
  ENDDO

  !     envelope and zo schemes

  if ((vnam=='TOP' .or. vnam=='TOD') .and. (ITOPSFLG==2 .or. ITOPSFLG==3) .and. &
       NP*NP<8) &
       print*,'Warning - trying to calc a std dev for: ', NP*NP, ' points'
  if ((vnam=='ZOT' .or. vnam=='ZOD') .and. IZ0FLG==1 .and. NP*NP<8) &
       print*,'Warning - trying to calc a std dev for: ', NP*NP, ' points'

  if (vnam=='TOP' .or. vnam=='TOD') &
       CALL TOPOQ(NIQ, NJQ, DELTAXQ, DELTAYQ, DATQ, SDQ, SHAQ, DATRE, NGR, N2, N3)
  if (vnam=='ZOT' .or. vnam=='ZOD')  &
       CALL ZOQ(NIQ, NJQ, DATQ, SDQ, NGR)

  XQ1 = (1. - 0.5*FLOAT(NIQ+1))*DELTAXQ + XCENTR
  YQ1 = (1. - 0.5*FLOAT(NJQ+1))*DELTAYQ + YCENTR

  !-srf ::: telesc
  !DO JR=1,N3
  !   DO IR=1,N2
  DO JR=jjj,jjj
     DO IR=iii,iii
        XR = (XT(IR)-XQ1)/DELTAXQ + 1.
        YR = (YT(JR)-YQ1)/DELTAYQ + 1.
        CALL GDTOST(DATQ, NIQ, NJQ, XR, YR, RVAL)

        !           envelope orog and zo schemes

        if (vnam=='ZOT' .or. vnam=='ZOD') then
           DATR(IR,JR) = MAX(0., RVAL)
           !               print*,'z0r',IR,JR,DATR(IR,JR)
        else
           DATR(IR,JR) = MAX(0., RVAL)
        endif

        CALL GDTOST(SDQ, NIQ, NJQ, XR, YR, RVAL)
        SDR(IR,JR) = MAX(0., RVAL)
     ENDDO
  ENDDO
  
  if (nmiss>0) then
     print *, '-----------------------------------------------------'
     print *, 'Input physiographical data file processing:'
     print *, '-----------------------------------------------------'
     print *, '  Input data blocks not found (data assumed to be zero):'
     do nn=1,nmiss
        print *, fnmiss(nn)
     enddo
     print *, '-----------------------------------------------------'
  endif

  !srf ELIMINAR PARA telescopica
  go to 222

  !     check to find the largest change in topo height

  if (vnam=='TOP' .or. vnam=='TOD') then
     diff    =    0.
     difflcl =    0.
     is      = -999
     js      = -999
     do j=2,n3-1
        do i=2,n2-1
           !print*,'1=',difflcl,max(difflcl,abs(datr(i,j)-datr(i-1,j)))
           difflcl = max(difflcl, abs(datr(i,j)-datr(i-1,j)))
           !print*,'2=',difflcl,max(difflcl,abs(datr(i,j)-datr(i+1,j)))
           difflcl = max(difflcl, abs(datr(i,j)-datr(i+1,j)))
           !print*,'3=',difflcl,max(difflcl,abs(datr(i,j)-datr(i,j-1)))
           difflcl = max(difflcl, abs(datr(i,j)-datr(i,j-1)))
           !print*,'4=',difflcl,max(difflcl,abs(datr(i,j)-datr(i,j+1)))
           difflcl = max(difflcl, abs(datr(i,j)-datr(i,j+1)))
           if (abs(diff-difflcl)>1.) then
              is = i
              js = j
           endif
           diff = max(diff,difflcl)
        enddo
     enddo
     write(6,100) ' Max d(topo) on grid @i,j=', ngr, is, js, diff
100  format(a,3i4,f8.1)
  endif

  !srf
222 continue

  deallocate(SDQ, SHAQ, SDR, DATRE)
  deallocate(ISO, IWO)

END subroutine sfcopqr_var

!**********************************************************************

subroutine get_TOPTWVL(ngr, n2, n3, TOPTWVL_2d)

  use mem_grid, only: &
       platn, plonn, xmn, ymn

  use io_params, only: &
       TOPTWVL

  implicit none
  ! Arguments:
  integer, intent(IN) :: n2, n3, ngr
  real, intent(OUT)   :: TOPTWVL_2d(n2,n3)
  ! Local Variables:
  real    :: rlat(n2,n3), rlon(n2,n3)
  integer :: i, j    

  !print*,n2,n3,platn(ngr),plonn(ngr),xmn(1,ngr),ymn(1,ngr)
  !-calculate lat, lon of each grid box T-points
  do j=1,n3
     do i=1,n2
        call xy_ll(rlat(i,j), rlon(i,j), platn(ngr), plonn(ngr), &
             xmn(i,ngr), ymn(j,ngr))
     enddo
  enddo

  do j=1,n3
     do i=1,n2
        TOPTWVL_2d(i,j) = abs(TOPTWVL(NGR))

        if (rlat(i,j)<-15.) then
           if (rlon(i,j)<-60.) TOPTWVL_2d(i,j) = 15.
        endif

        if (rlat(i,j)>=-15. .and. rlat(i,j)<-9.) then
           if (rlon(i,j)<-62.) TOPTWVL_2d(i,j) = 15.
        endif

        if (rlat(i,j)>=-9. .and. rlat(i,j)<-1.) then
           if (rlon(i,j)<-70.) TOPTWVL_2d(i,j) = 15.
        endif

        if (rlat(i,j)>=-1.) then
           if (rlon(i,j)<=-57 .and. rlon(i,j)>-67.) TOPTWVL_2d(i,j) = 11.
        endif

        if (rlat(i,j)>=-1.) then
           if (rlon(i,j)<=-67.) TOPTWVL_2d(i,j) = 15.
        endif

        if (rlat(i,j)<=-10.) then
           if (rlon(i,j)>=40.)  TOPTWVL_2d(i,j) = 11.
        endif
        !print*,i,j,TOPTWVL_2d(i,j)
     enddo
  enddo

end subroutine get_TOPTWVL

!**********************************************************************
!--------------------------------------------------------------------------
subroutine geodat_var_opt(n2, n3, datr, hfn, ofn, vt2da, vt2db, ngr, vnam)

  use mem_grid, only: &
       ngrid, nnxp, nnyp, deltaxn, deltayn, xtn, ytn, platn, plonn

  use io_params, only: &
       itopsflg, toptenh, toptwvl, iz0flg

  implicit none
  
  include "files.h"

  ! Arguments:
  integer, intent(IN)            :: n2, n3
  real, intent(INOUT)            :: vt2da(*), vt2db(*), datr(n2,n3)
  character(len=*), intent(IN) :: hfn, ofn
  character(len=03), intent(IN)  :: vnam
  ! Local Variables:
  integer :: ngr

  
  character(len=f_name_length) :: title
  integer :: lb, iblksizo, no, isbego, iwbego, iodim, mof, niq, njq, np
  real    :: offlat, offlon, deltallo, deltaxq, deltayq, deltaxp, deltayp, erad
  integer :: ierr
  real, allocatable :: dato(:)

  integer, parameter :: max_number_toptwvl=20
  
!--srf 4.3 setting : 10 - 25 - 25
!--srf 5.1 test  1 : 10 - 15 - 15
  
  real, dimension(2:max_number_toptwvl) :: TOPTWVLX=(/&
  0., & ! 2
  0., & ! 3
  0., & ! 4
  0., & ! 5
  10., & ! 6
  0., & ! 7
  0., & ! 8
  0., & ! 9
  0., & ! 10
  15.,& ! 11   ! 15 , 25 
  0., & ! 12
  0., & ! 13
  0., & ! 14
  0., & ! 15
  0., & ! 16
  0., & ! 17
  0., & ! 18
  0., & ! 19 
  15. & ! 19  ! 15, 25
  /)
    
  real, allocatable :: datrx(:,:,:)
  integer itpx
  
  TOPTWVLX(nint(abs(TOPTWVL(NGR))))=abs(TOPTWVL(NGR))
  print*, 'TOPTWVLX=',TOPTWVLX

 
  LB = len_trim(HFN)
  if (LB<=0) then
     print *, '==================================================='
     print *, '|  Problem in GEODAT, Input data prefix incorrect !'
     print *, '|  Grid :', ngrid
     print *, '|  File prefix:', trim(HFN)
     print *, '==================================================='
     call fatal_error('GEODAT-file')
  endif

  if ( (vnam(1:2)=='TO' .or. vnam(1:2)=='ZO') .and.  &
       (ITOPSFLG(NGR)==1 .and. TOPTENH(NGR)>1.)) then
     print *, '==================================================='
     print *, '|  Problem in GEODAT, Silhouette weight too high !'
     print *, '|  Grid :', NGR
     print *, '|  Weight (range TOPTENH=0-1):', TOPTENH(NGR)
     print *, '==================================================='
     call fatal_error('GEODAT')
  endif

  !     configure grid specs for raw data to rams grid (R) transfer

  !     raw data grid (O)
  if (vnam=='TOD' .or. vnam=='ZOD') then  !using dted data
     iblksizo =   5 ! 5 degree squares
     no       = 600 ! 30" data over a 5 degree square
     isbego   = -90
     iwbego   =   0
     offlat   =   0.
     offlon   =   0.
     DELTALLO = FLOAT(IBLKSIZO)/FLOAT(NO)
  else
     TITLE    = HFN(1:LB)//'HEADER'
     LB       = len_trim(TITLE)
     call rams_f_open(29, title(1:lb), 'FORMATTED', 'OLD', 'READ', 0)
     read (29,*) IBLKSIZO, NO, ISBEGO, IWBEGO, offlat, offlon
     close (29)
     DELTALLO = FLOAT(IBLKSIZO)/FLOAT(NO-1)
  endif

  iodim = max(100000, 4*no*no)
  MOF   = IODIM/(NO*NO)

  allocate(dato(iodim+mof+mof), STAT=ierr)
  if (ierr/=0) call fatal_error("Error allocating dato")
  
  IF (vnam(1:2)=='TO') then
     allocate(datrx(max_number_toptwvl,n2,n3), STAT=ierr)
     if (ierr/=0) call fatal_error("Error allocating datrx")
     datrx= 0.0
     do itpx=2,max_number_toptwvl
   	   if(TOPTWVLX(itpx) == 0.0) CYCLE
   	   print*,'gerando topo para TOPTWVL=',itpx,TOPTWVLX(itpx)
   	   !	 temp grid (Q) - smoothing only applied to topo
   	      !DELTAXQ = 0.5*abs(TOPTWVL(NGR))*DELTAXN(NGR)
   	      !DELTAYQ = 0.5*abs(TOPTWVL(NGR))*DELTAYN(NGR)
   	      DELTAXQ = 0.5*abs(TOPTWVLX(itpx))*DELTAXN(NGR)
   	      DELTAYQ = 0.5*abs(TOPTWVLX(itpx))*DELTAYN(NGR)
   	   NIQ = int(FLOAT(NNXP(NGR)-1)*DELTAXN(NGR)/DELTAXQ) + 4
   	   NJQ=int(FLOAT(NNYP(NGR)-1)*DELTAYN(NGR)/DELTAYQ)+4

   	   !	 interpollated raw data grid (P)
   	   NP	   = min(10, max(1,int(DELTAXQ/(DELTALLO*111000.))))
   	   DELTAXP = DELTAXQ/FLOAT(NP)
   	   DELTAYP = DELTAYQ/FLOAT(NP)
         print*,'#'
   	   call SFCOPQR(NO, MOF, NP, NIQ, NJQ, N2, N3, XTN(1,NGR), YTN(1,NGR), &
   		platn(ngr), plonn(ngr), 				       &
   		ERAD, DELTALLO, DELTAXP, DELTAYP, DELTAXQ, DELTAYQ, IBLKSIZO,  &
   		ISBEGO, IWBEGO, DATO(1), VT2DA, VT2DB, DATRX(ITPX,:,:), 	       &
   		OFN, offlat, offlon, VNAM, NGR, itopsflg(ngr), iz0flg(ngr))
         print*,'##'
     enddo
     print*,'topo para TOPTWVL ',max_number_toptwvl, ' gerado ...'
     call datrx2datr(ngr,n2,n3,max_number_toptwvl,TOPTWVLX,DATRX,DATR)
     deallocate(datrx)

 ELSE
        DELTAXQ = DELTAXN(NGR)
 	DELTAYQ = DELTAYN(NGR)
    	NIQ = int(FLOAT(NNXP(NGR)-1)*DELTAXN(NGR)/DELTAXQ) + 4
    	NJQ=int(FLOAT(NNYP(NGR)-1)*DELTAYN(NGR)/DELTAYQ)+4

    	!     interpollated raw data grid (P)
    	NP	= min(10, max(1,int(DELTAXQ/(DELTALLO*111000.))))
    	DELTAXP = DELTAXQ/FLOAT(NP)
    	DELTAYP = DELTAYQ/FLOAT(NP)

    	call SFCOPQR(NO, MOF, NP, NIQ, NJQ, N2, N3, XTN(1,NGR), YTN(1,NGR), &
    	     platn(ngr), plonn(ngr),					    &
    	     ERAD, DELTALLO, DELTAXP, DELTAYP, DELTAXQ, DELTAYQ, IBLKSIZO,  &
    	     ISBEGO, IWBEGO, DATO(1), VT2DA, VT2DB, DATR,		    &
    	     OFN, offlat, offlon, VNAM, NGR, itopsflg(ngr), iz0flg(ngr))
 ENDIF  

 deallocate(dato)

end subroutine geodat_var_opt
!-----------------------------------------------------------------------------
subroutine datrx2datr(ngr,n2,n3,max_number_toptwvl,TOPTWVLX,DATRX,DATR)

  use mem_grid, only: &
       platn, plonn, xmn, ymn

  use io_params, only: &
       TOPTWVL

  implicit none
  ! Arguments:
  integer, intent(IN) :: n2, n3, ngr,max_number_toptwvl
  real, intent(IN)   :: DATRX(max_number_toptwvl,n2,n3)
  real, intent(IN)   :: TOPTWVLX(max_number_toptwvl)
  
  real, intent(OUT)   :: DATR(n2,n3)
  ! Local Variables:
  real    :: rlat(n2,n3), rlon(n2,n3),TOPTWVL_2d(n2,n3)
  integer :: i, j    

  !print*,n2,n3,platn(ngr),plonn(ngr),xmn(1,ngr),ymn(1,ngr)
  !-calculate lat, lon of each grid box T-points
  do j=1,n3
     do i=1,n2
        call xy_ll(rlat(i,j), rlon(i,j), platn(ngr), plonn(ngr), &
             xmn(i,ngr), ymn(j,ngr))
     enddo
  enddo

!ORIG 
!GO TO 100
!  do j=1,n3
!     do i=1,n2
!        TOPTWVL_2d(i,j) = abs(TOPTWVL(NGR))!
!
!        if (rlat(i,j)<-15.) then
!           if (rlon(i,j)<-60.) TOPTWVL_2d(i,j) = 15.
!        endif
!
!	 if (rlat(i,j)>=-15. .and. rlat(i,j)<-9.) then
!	    if (rlon(i,j)<-62.) TOPTWVL_2d(i,j) = 15.
!	 endif
!
!	 if (rlat(i,j)>=-9. .and. rlat(i,j)<-1.) then
!	    if (rlon(i,j)<-70.) TOPTWVL_2d(i,j) = 15.
!	 endif
!
!	 if (rlat(i,j)>=-1.) then
!	    if (rlon(i,j)<=-57 .and. rlon(i,j)>-67.) TOPTWVL_2d(i,j) = 11.
!	 endif
!
!	 if (rlat(i,j)>=-1.) then
!	    if (rlon(i,j)<=-67.) TOPTWVL_2d(i,j) = 15.
!	 endif
!
!	 if (rlat(i,j)<=-10.) then
!	    if (rlon(i,j)>=40.)  TOPTWVL_2d(i,j) = 11.
!	 endif
!	 !print*,i,j,TOPTWVL_2d(i,j)
!     enddo
!  enddo
!GO TO 200
!
!
!100  CONTINUE
  do j=1,n3
     do i=1,n2
        TOPTWVL_2d(i,j) = abs(TOPTWVL(NGR))
	
	if (rlat(i,j) <-1.) then
           if (rlon(i,j)<-62.) TOPTWVL_2d(i,j) = 20.!15.
        endif

        if (rlat(i,j)>=-1.) then
           if (rlon(i,j)<-68.) TOPTWVL_2d(i,j) = 20.!15.
        endif

        if (rlat(i,j)>=-1) then
           if (rlon(i,j)>=-68. .AND. rlon(i,j) <= -58.) TOPTWVL_2d(i,j) = 11.
        endif


        if (rlat(i,j)<=-30. .AND. rlat(i,j) <= -17.) then
           if (rlon(i,j)>=60.)  TOPTWVL_2d(i,j) = 6.
        endif
     enddo
  enddo

!200 CONTINUE
 if(nint(maxval(TOPTWVL_2d))>max_number_toptwvl) &
 stop 'maxval(TOPTWVL_2d must be <= max_number_toptwvl'
 
  do j=1,n3
     do i=1,n2
     datr(i,j)=datrx(nint(TOPTWVL_2d(i,j)),i,j)
     enddo
  enddo

  print*,'min max topo=',minval(datr),maxval(datr)


end subroutine datrx2datr
