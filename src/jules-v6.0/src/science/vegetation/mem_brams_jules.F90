!############################# Change Log ##################################
! 5.0.0
!
! Demerval S. Moreira - 27/Ago/2012
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


Module mem_brams_jules
   implicit none
   integer:: mynumB, nxB,nyB,ntB,timestepB,ntimestepB,sm_levelsB,output_periodB,dump_periodB
   real, allocatable :: glatB(:,:), glonB(:,:), land_fracB(:,:), dzsoilB(:), &
                        precipB(:,:), swdownB(:,:), lwdownB(:,:), &
                        diff_radB(:,:), tempB(:,:), upsB(:,:), vpsB(:,:), &
                        pstarB(:,:), qB(:,:), fracB(:,:,:), sthuB(:,:,:),tsoilB(:,:,:), &
                        tstarB(:,:),viesALL(:,:,:,:),z1_uvB,z1_tqB
   CHARACTER(len=19) :: main_run_startB,main_run_endB
   CHARACTER(len=256) :: dir_run_idB,hfilinB
   CHARACTER(len=16) :: runtypeB
   
   
End Module mem_brams_jules
