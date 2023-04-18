module potential

use util;

IMPLICIT NONE;

CONTAINS


!%%%%%%%%%%%%%%%%%%%%%%%%%%%   Hydrogen

!%%%%%%%%%   Calculation of the atomic potential
complex*16  FUNCTION pot(x)

  IMPLICIT NONE;

  complex*16, intent(in) :: x;

  real*8 :: Z,A,B,C,F,G,Rx;

   Z = 1.0d0;

   A = 1.51d0; !! This is He


   pot = -Z/sqrt(A*A + x*x);


end function pot


!%%%%%%%%%   Calculation of the gradient along z of potential
complex*16  FUNCTION gradzpot(x)


 IMPLICIT NONE;

 complex*16, intent(in) :: x;

 real*8 :: Z,A,B,C,F,G,Rx;

  Z = 1.0d0;

  A = 1.51d0; !! This is He


  gradzpot = Z*x/(A*A + x*x)**(1.5d0);

end function gradzpot





end module potential
