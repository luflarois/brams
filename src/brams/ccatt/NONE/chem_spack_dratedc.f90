  
 MODULE mod_chem_spack_dratedc
  
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: dratedc ! subroutine
 CONTAINS
  
   SUBROUTINE dratedc(rk,y,dw,ngas,ijkbeg,ijkend,maxblock_size,nr)
 
!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!
!     This routine computes the derivative of reaction  rates.
!     This routine is automatically generated by SPACK.
!     Mechanism: ../Mechanism/CB07   
!     Species: ../Mechanism/ciCB07 
!
!------------------------------------------------------------------------
!
!     -- INPUT VARIABLES
!
!     RK: kinetic rates.
!     Y: chemical concentrations.
!
!     -- INPUT/OUTPUT VARIABLES
!
!     -- OUTPUT VARIABLES
!
!     DW: derivative of reaction rates wrt Y.
!
!------------------------------------------------------------------------
!
!     -- REMARKS
!
!------------------------------------------------------------------------
!
!     -- MODIFICATIONS
!
!------------------------------------------------------------------------
!
!     -- AUTHOR(S)
!
!     SPACK.
!
!------------------------------------------------------------------------
 
      IMPLICIT NONE
 
 
 
     INTEGER	       , INTENT(IN)  :: ngas                      
     INTEGER	       , INTENT(IN)  :: ijkbeg			  
     INTEGER	       , INTENT(IN)  :: ijkend			  
     INTEGER	       , INTENT(IN)  :: maxblock_size		  
     INTEGER	       , INTENT(IN)  :: nr			  
     DOUBLE PRECISION , INTENT(IN)  :: rk(maxblock_size,nr)	  
     DOUBLE PRECISION , INTENT(IN)  :: y(maxblock_size,NGAS) 	  
     DOUBLE PRECISION , INTENT(OUT) :: dw(maxblock_size,nr,NGAS) 
 
   END SUBROUTINE dratedc
 
  END MODULE mod_chem_spack_dratedc
 