!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine isenio (inout,iun,n1,n2)

  use isan_coms

  implicit none
  include "i8.h"
  include "UseVfm.h"
  integer :: iun,n1,n2
  character(len=*) :: inout
  integer :: nlt,nx3,ny3,ninn,l !npts
  integer(kind=8) :: npts

  if (useVfm) then

     ! vfm coded file

     if(inout.eq.'IN') THEN
        read(iun,920) iyy,imm,idd,ihh,nx3,ny3,ninn,(levth(l),l=1,ninn)
920     format(7i4,(13i6))
        if(nx3.ne.n1.or.ny3.ne.n2.or.ninn.ne.nisn) then
           print*,'Isentropic stage grid dimensions do not match'
           print*,'   configuration file on read !'
           print*,' File dimens - ',nx3,ny3,ninn
           print*,' Run  dimens - ',n1,n2,nisn
           stop 'IO3-2'
        endif

        npts=n1*n2
        do nlt=1,nisn
           call vfirec(iun,pi_u(1,1,nlt),npts,'LIN')
           call vmissr(pi_u(:,:,nlt),npts,1e30,-998.)
           call vfirec(iun,pi_v(1,1,nlt),npts,'LIN')
           call vmissr(pi_v(:,:,nlt),npts,1e30,-998.)
           call vfirec(iun,pi_s(1,1,nlt),npts,'LIN')
           call vmissr(pi_p(:,:,nlt),npts,1e30,-.5)
           call vfirec(iun,pi_p(1,1,nlt),npts,'LIN')
           call vmissr(pi_s(:,:,nlt),npts,1e30,-.5)
           call vfirec(iun,pi_r(1,1,nlt),npts,'LIN')
           call vmissr(pi_r(:,:,nlt),npts,1e30,-.5)
        enddo

        call vfirec(iun,rs_u,npts,'LIN')
        call vmissr(rs_u,npts,1e30,-998.)
        call vfirec(iun,rs_v,npts,'LIN')
        call vmissr(rs_v,npts,1e30,-998.)
        call vfirec(iun,rs_p,npts,'LIN')
        call vmissr(rs_p,npts,1e30,-.5)
        call vfirec(iun,rs_t,npts,'LIN')
        call vmissr(rs_t,npts,1e30,-.5)
        call vfirec(iun,rs_r,npts,'LIN')
        call vmissr(rs_r,npts,1e30,-.5)
        call vfirec(iun,rs_s,npts,'LIN')
        call vmissr(rs_s,npts,1e30,-.5)
        call vfirec(iun,rs_top,npts,'LIN')
        call vmissr(rs_top,npts,1e30,-.5)
        call vfirec(iun,rs_qual,npts,'LIN')
        call vmissr(rs_qual,npts,1e30,-.5)

        call vfirec(iun,rs_slp,npts,'LIN')
        call vmissr(rs_slp,npts,1e30,-.5)
        call vfirec(iun,rs_sfp,npts,'LIN')
        call vmissr(rs_sfp,npts,1e30,-.5)
        call vfirec(iun,rs_sft,npts,'LIN')
        call vmissr(rs_sft,npts,1e30,-.5)
        call vfirec(iun,rs_snow,npts,'LIN')
        call vmissr(rs_snow,npts,1e30,-.5)
        call vfirec(iun,rs_sst,npts,'LIN')
        call vmissr(rs_sst,npts,1e30,-.5)

        print 201,' *****  Isentropic file input *****************'  &
             ,iyear,imonth,idate,ihour,n1,n2,nisn  &
             ,(levth(l),l=1,nisn)
201     format(//,a,//  &
             ,' *',7X,' Date (year,month,day,hour)  - ',4I5,/  &
             ,' *',7X,' Number of X,Y points        - ',2I5,/  &
             ,' *',7X,' Number of isentropic levels - ',I5,/  &
             ,' *',7X,' Isentropic levels (K)       - '/,(32X,8I5))
        print '(a)',' **********************************************'

     endif

     if(inout.eq.'out') then

        write(iun,920) iyear,imonth,idate,ihour,n1,n2,nisn  &
             ,(levth(l),l=1,nisn)

        npts=n1*n2
        do nlt=1,nisn
           call vmissw(pi_u(:,:,nlt),npts,pi_scra(:,:,nlt),1E30,-999.)
           call vforec(iun,pi_scra,npts,18,pi_scrb,'LIN')
           call vmissw(pi_v(:,:,nlt),npts,pi_scra(:,:,nlt),1E30,-999.)
           call vforec(iun,pi_scra,npts,18,pi_scrb,'LIN')
           call vmissw(pi_p(:,:,nlt),npts,pi_scra(:,:,nlt),1E30,-1.)
           call vforec(iun,pi_scra,npts,18,pi_scrb,'LIN')
           call vmissw(pi_s(:,:,nlt),npts,pi_scra(:,:,nlt),1E30,-1.)
           call vforec(iun,pi_scra,npts,18,pi_scrb,'LIN')
           call vmissw(pi_r(:,:,nlt),npts,pi_scra(:,:,nlt),1E30,-1.)
           call vforec(iun,pi_scra,npts,18,pi_scrb,'LIN')
        ENDDO

        call vmissw(rs_u,npts,pi_scra(:,:,1),1E30,-999.)
        call vforec(iun,pi_scra,npts,18,pi_scrb,'LIN')
        call vmissw(rs_v,npts,pi_scra(:,:,1),1E30,-999.)
        call vforec(iun,pi_scra,npts,18,pi_scrb,'LIN')
        call vmissw(rs_p,npts,pi_scra(:,:,1),1E30,-1.)
        call vforec(iun,pi_scra,npts,18,pi_scrb,'LIN')
        call vmissw(rs_t,npts,pi_scra(:,:,1),1E30,-1.)
        call vforec(iun,pi_scra,npts,18,pi_scrb,'LIN')
        call vmissw(rs_r,npts,pi_scra(:,:,1),1E30,-1.)
        call vforec(iun,pi_scra,npts,18,pi_scrb,'LIN')
        call vmissw(rs_s,npts,pi_scra(:,:,1),1E30,-1.)
        call vforec(iun,pi_scra,npts,18,pi_scrb,'LIN')
        call vmissw(rs_top,npts,pi_scra(:,:,1),1E30,-1.)
        call vforec(iun,pi_scra,npts,18,pi_scrb,'LIN')
        call vmissw(rs_qual,npts,pi_scra(:,:,1),1E30,-1.)
        call vforec(iun,pi_scra,npts,18,pi_scrb,'LIN')

        call vmissw(rs_slp,npts,pi_scra(:,:,1),1E30,-1.)
        call vforec(iun,pi_scra,npts,18,pi_scrb,'LIN')
        call vmissw(rs_sfp,npts,pi_scra(:,:,1),1E30,-1.)
        call vforec(iun,pi_scra,npts,18,pi_scrb,'LIN')
        call vmissw(rs_sft,npts,pi_scra(:,:,1),1E30,-1.)
        call vforec(iun,pi_scra,npts,18,pi_scrb,'LIN')
        call vmissw(rs_snow,npts,pi_scra(:,:,1),1E30,-1.)
        call vforec(iun,pi_scra,npts,18,pi_scrb,'LIN')
        call vmissw(rs_sst,npts,pi_scra(:,:,1),1E30,-1.)
        call vforec(iun,pi_scra,npts,18,pi_scrb,'LIN')

        print 201,' *****  Isentropic file written *************'  &
             ,iyear,imonth,idate,ihour,n1,n2,nisn  &
             ,(levth(l),l=1,nisn)

        print 303,igridfl,gobsep,gobrad
303     format(/,  &
             ' Grid flag (IGRIDFL)               -',I4,/  &
             ,' Grid-obs separation in degrees    -',F5.2,/  &
             ,' Grid-obs radius influence degrees -',F5.2)

     endif

  else

     ! binary file

     if(inout.eq.'IN') THEN

        read (iun) iyy,imm,idd,ihh,nx3,ny3,ninn,(levth(l),l=1,ninn)
        if(nx3.ne.n1.or.ny3.ne.n2.or.ninn.ne.nisn) then
           print*,'Isentropic stage grid dimensions do not match'
           print*,'   configuration file on read !'
           print*,' File dimens - ',nx3,ny3,ninn
           print*,' Run  dimens - ',n1,n2,nisn
           stop 'IO3-2'
        endif

        npts=n1*n2
        do nlt=1,nisn
           read (iun) pi_u(:,:,nlt)
           call vmissr(pi_u(:,:,nlt),npts,1e30,-998.)
           read (iun) pi_v(:,:,nlt)
           call vmissr(pi_v(:,:,nlt),npts,1e30,-998.)
           read (iun) pi_p(:,:,nlt)
           call vmissr(pi_p(:,:,nlt),npts,1e30,-.5)
           read (iun) pi_s(:,:,nlt)
           call vmissr(pi_s(:,:,nlt),npts,1e30,-.5)
           read (iun) pi_r(:,:,nlt)
           call vmissr(pi_r(:,:,nlt),npts,1e30,-.5)
        enddo

        read (iun) rs_u(:,:)
        call vmissr(rs_u,npts,1e30,-998.)
        read (iun) rs_v(:,:)
        call vmissr(rs_v,npts,1e30,-998.)
        read (iun) rs_p(:,:)
        call vmissr(rs_p,npts,1e30,-.5)
        read (iun) rs_t(:,:)
        call vmissr(rs_t,npts,1e30,-.5)
        read (iun) rs_r(:,:)
        call vmissr(rs_r,npts,1e30,-.5)
        read (iun) rs_s(:,:)
        call vmissr(rs_s,npts,1e30,-.5)
        read (iun) rs_top(:,:)
        call vmissr(rs_top,npts,1e30,-.5)
        read (iun) rs_qual(:,:)
        call vmissr(rs_qual,npts,1e30,-.5)

        read (iun) rs_slp(:,:)
        call vmissr(rs_slp,npts,1e30,-.5)
        read (iun) rs_sfp(:,:)
        call vmissr(rs_sfp,npts,1e30,-.5)
        read (iun) rs_sft(:,:)
        call vmissr(rs_sft,npts,1e30,-.5)
        read (iun) rs_snow(:,:)
        call vmissr(rs_snow,npts,1e30,-.5)
        read (iun) rs_sst(:,:)
        call vmissr(rs_sst,npts,1e30,-.5)

        print 201,' *****  Isentropic file input *****************'  &
             ,iyear,imonth,idate,ihour,n1,n2,nisn  &
             ,(levth(l),l=1,nisn)
        print '(a)',' **********************************************'

     endif

     if(inout.eq.'out') then

        write (iun) iyear,imonth,idate,ihour,n1,n2,nisn,(levth(l),l=1,nisn)

        npts=n1*n2
        do nlt=1,nisn
           call vmissw(pi_u(:,:,nlt),npts,pi_scra(:,:,nlt),1E30,-999.)
           write (iun) pi_scra(:,:,nlt)
           call vmissw(pi_v(:,:,nlt),npts,pi_scra(:,:,nlt),1E30,-999.)
           write (iun) pi_scra(:,:,nlt)
           call vmissw(pi_p(:,:,nlt),npts,pi_scra(:,:,nlt),1E30,-1.)
           write (iun) pi_scra(:,:,nlt)
           call vmissw(pi_s(:,:,nlt),npts,pi_scra(:,:,nlt),1E30,-1.)
           write (iun) pi_scra(:,:,nlt)
           call vmissw(pi_r(:,:,nlt),npts,pi_scra(:,:,nlt),1E30,-1.)
           write (iun) pi_scra(:,:,nlt)
        ENDDO

        call vmissw(rs_u,npts,pi_scra(:,:,1),1E30,-999.)
        write (iun) pi_scra(:,:,nlt)
        call vmissw(rs_v,npts,pi_scra(:,:,1),1E30,-999.)
        write (iun) pi_scra(:,:,nlt)
        call vmissw(rs_p,npts,pi_scra(:,:,1),1E30,-1.)
        write (iun) pi_scra(:,:,nlt)
        call vmissw(rs_t,npts,pi_scra(:,:,1),1E30,-1.)
        write (iun) pi_scra(:,:,nlt)
        call vmissw(rs_r,npts,pi_scra(:,:,1),1E30,-1.)
        write (iun) pi_scra(:,:,nlt)
        call vmissw(rs_s,npts,pi_scra(:,:,1),1E30,-1.)
        write (iun) pi_scra(:,:,nlt)
        call vmissw(rs_top,npts,pi_scra(:,:,1),1E30,-1.)
        write (iun) pi_scra(:,:,nlt)
        call vmissw(rs_qual,npts,pi_scra(:,:,1),1E30,-1.)
        write (iun) pi_scra(:,:,nlt)

        call vmissw(rs_slp,npts,pi_scra(:,:,1),1E30,-1.)
        write (iun) pi_scra(:,:,nlt)
        call vmissw(rs_sfp,npts,pi_scra(:,:,1),1E30,-1.)
        write (iun) pi_scra(:,:,nlt)
        call vmissw(rs_sft,npts,pi_scra(:,:,1),1E30,-1.)
        write (iun) pi_scra(:,:,nlt)
        call vmissw(rs_snow,npts,pi_scra(:,:,1),1E30,-1.)
        write (iun) pi_scra(:,:,nlt)
        call vmissw(rs_sst,npts,pi_scra(:,:,1),1E30,-1.)
        write (iun) pi_scra(:,:,nlt)

        print 201,' *****  Isentropic file written *************'  &
             ,iyear,imonth,idate,ihour,n1,n2,nisn  &
             ,(levth(l),l=1,nisn)

        print 303,igridfl,gobsep,gobrad

     endif




  end if
  return
end subroutine isenio

!***************************************************************************

subroutine sigzio (inout,iun,n1,n2)

  use isan_coms

  implicit none
  include "i8.h"
  include "UseVfm.h"
  integer :: iun,n1,n2
  character(len=*) :: inout

  integer :: nlt,l,ninn,nx3,ny3

  integer(kind=i8) :: npts


    if (useVfm) then

       ! vfm coded file

       if(inout.eq.'IN') then
          read(iun,920) iyy,imm,idd,ihh,nx3,ny3,ninn  &
               ,(sigz(l),l=1,ninn)
920       format(7i4,(9f8.2))
          if(nx3.ne.n1.or.ny3.ne.n2.or.ninn.ne.nsigz)then
             print*,'Sigma-z grid dimensions do not match'
             print*,'   input data on read !'
             print*,' File  dimensions - ',nx3,ny3,ninn
             print*,' Input dimensions - ',n1,n2,nsigz
             stop 'iO3-2'
          endif

          npts=n1*n2
          do nlt=1,nsigz
             call vfirec(iun,ps_u(1,1,nlt),npts,'LIN')
             call vmissr(ps_u(:,:,nlt),npts,1e30,-998.)
             call vfirec(iun,ps_v(1,1,nlt),npts,'LIN')
             call vmissr(ps_v(:,:,nlt),npts,1e30,-998.)
             call vfirec(iun,ps_p(1,1,nlt),npts,'LIN')
             call vmissr(ps_p(:,:,nlt),npts,1e30,-.5)
             call vfirec(iun,ps_t(1,1,nlt),npts,'LIN')
             call vmissr(ps_t(:,:,nlt),npts,1e30,-.5)
             call vfirec(iun,ps_r(1,1,nlt),npts,'LIN')
             call vmissr(ps_r(:,:,nlt),npts,1e30,-.5)
          enddo

          print 201,' *****  Sigma-z file input *****************'  &
               ,iyear,imonth,idate,ihour,n1,n2,nsigz  &
               ,(sigz(l),l=1,nsigz)
201       format(//,a,//  &
               ,' *',7X,' Date (year,month,day,hour)  - ',4I5,/  &
               ,' *',7X,' Number of X,Y points        - ',2I5,/  &
               ,' *',7X,' Number of sigma-z levels    - ',I5,/  &
               ,' *',7X,' Sigma-z levels (m)          - '/,(32X,7F8.1))
          print '(a)',' **********************************************'

       else if(inout.eq.'OUT') then
          write(iun,920) iyear,imonth,idate,ihour,n1,n2,nsigz  &
               ,(sigz(l),l=1,nsigz)

          npts=n1*n2
          do nlt=1,nsigz
             call vmissw(ps_u(:,:,nlt),npts,ps_scra(:,:,nlt),1E30,-999.)
             call vforec(iun,ps_scra,npts,18,ps_scrb,'LIN')
             call vmissw(ps_v(:,:,nlt),npts,ps_scra(:,:,nlt),1E30,-999.)
             call vforec(iun,ps_scra,npts,18,ps_scrb,'LIN')
             call vmissw(ps_p(:,:,nlt),npts,ps_scra(:,:,nlt),1E30,-1.)
             call vforec(iun,ps_scra,npts,18,ps_scrb,'LIN')
             call vmissw(ps_t(:,:,nlt),npts,ps_scra(:,:,nlt),1E30,-1.)
             call vforec(iun,ps_scra,npts,18,ps_scrb,'LIN')
             call vmissw(ps_r(:,:,nlt),npts,ps_scra(:,:,nlt),1E30,-1.)
             call vforec(iun,ps_scra,npts,18,ps_scrb,'LIN')
          enddo

          print 201,' *****  Sigma-z file written *************'  &
               ,iyear,imonth,idate,ihour,n1,n2,nsigz   &
               ,(sigz(l),l=1,nsigz)

       endif

    else

       ! binary file

       if(inout.eq.'IN') then
          read(iun) iyy,imm,idd,ihh,nx3,ny3,ninn  &
               ,(sigz(l),l=1,ninn)
          if(nx3.ne.n1.or.ny3.ne.n2.or.ninn.ne.nsigz)then
             print*,'Sigma-z grid dimensions do not match'
             print*,'   input data on read !'
             print*,' File  dimensions - ',nx3,ny3,ninn
             print*,' Input dimensions - ',n1,n2,nsigz
             stop 'iO3-2'
          endif

          npts=n1*n2
          do nlt=1,nsigz
             read(iun) ps_u(:,:,nlt)
             call vmissr(ps_u(:,:,nlt),npts,1e30,-998.)
             read(iun) ps_v(:,:,nlt)
             call vmissr(ps_v(:,:,nlt),npts,1e30,-998.)
             read(iun) ps_p(:,:,nlt)
             call vmissr(ps_p(:,:,nlt),npts,1e30,-.5)
             read(iun) ps_t(:,:,nlt)
             call vfirec(iun,ps_t(1,1,nlt),npts,'LIN')
             read(iun) ps_r(:,:,nlt)
             call vmissr(ps_r(:,:,nlt),npts,1e30,-.5)
          enddo

          print 201,' *****  Sigma-z file input *****************'  &
               ,iyear,imonth,idate,ihour,n1,n2,nsigz  &
               ,(sigz(l),l=1,nsigz)
          print '(a)',' **********************************************'

       else if(inout.eq.'OUT') then
          write(iun) iyear,imonth,idate,ihour,n1,n2,nsigz  &
               ,(sigz(l),l=1,nsigz)

          npts=n1*n2
          do nlt=1,nsigz
             call vmissw(ps_u(:,:,nlt),npts,ps_scra(:,:,nlt),1E30,-999.)
             write(iun) ps_scra(:,:,nlt)
             call vmissw(ps_v(:,:,nlt),npts,ps_scra(:,:,nlt),1E30,-999.)
             write(iun) ps_scra(:,:,nlt)
             call vmissw(ps_p(:,:,nlt),npts,ps_scra(:,:,nlt),1E30,-1.)
             write(iun) ps_scra(:,:,nlt)
             call vmissw(ps_t(:,:,nlt),npts,ps_scra(:,:,nlt),1E30,-1.)
             write(iun) ps_scra(:,:,nlt)
             call vmissw(ps_r(:,:,nlt),npts,ps_scra(:,:,nlt),1E30,-1.)
             write(iun) ps_scra(:,:,nlt)
          enddo

          print 201,' *****  Sigma-z file written *************'  &
               ,iyear,imonth,idate,ihour,n1,n2,nsigz   &
               ,(sigz(l),l=1,nsigz)

       endif
    end if

  return
end subroutine sigzio

!***************************************************************************

subroutine vmissw (af,n,as,fm,fx)

  implicit none
  include "i8.h"
  integer(kind=i8), intent(in) :: n
  real, intent(in)             :: af(n), fm, fx
  real, intent(inout)          :: as(n)
  integer(kind=i8)             :: i

  do i=1,n
     as(i)=af(i)
     if(af(i)>=fm) as(i)=fx
  enddo

  return
end subroutine vmissw

!***************************************************************************

subroutine vmissr (af,n,fm,fx)

  implicit none
  include "i8.h"
  integer(kind=i8), intent(in) :: n
  real, intent(inout)          :: af(n)
  real, intent(in)             :: fm, fx
  integer(kind=i8)             :: i

  do i=1,n
     if(af(i)<=fx) af(i)=fm
  enddo

  return
end subroutine vmissr

