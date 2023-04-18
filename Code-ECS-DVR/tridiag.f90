module tridiag

use util;

IMPLICIT NONE;

CONTAINS


!%%%%%%%%%   Calculation of the norm of a vector expanded on l
real*8  FUNCTION norm(psi,r,Nr,Lmax)

  IMPLICIT NONE;

  integer, intent(in) :: Nr,Lmax;
  Complex*16,dimension(:,:),Allocatable :: psi;
  real*8,dimension(:),Allocatable :: r;

  integer :: i,l;
  real*8 :: h;
  h = r(2)-r(1);
  norm = 0.0d0;
  do i = 1 , Nr
    do l = 0 , Lmax
       !norm = norm + 0.5*real(psi(l,i)*dconjg(psi(l,i)) + psi(l,i+1)*dconjg(psi(l,i+1)))*(r(i+1)-r(i));
       norm = norm + real(psi(l,i)*dconjg(psi(l,i)))*h;
    enddo
  enddo

  ! Here there is no need to multiply by r^2 since psi = phi/r
  norm = dsqrt(norm);

end function norm



!%%%%%%%%%   Calculation of the norm of a vector
real*8  FUNCTION norm_vect(psi,r,Nr)

  IMPLICIT NONE;

  integer, intent(in) :: Nr;
  Complex*16,dimension(:),Allocatable :: psi;
  real*8,dimension(:),Allocatable :: r;

  integer :: i;
  real*8 :: h;
  h = r(2)-r(1);

  norm_vect = 0.d0;
  do i = 1 , Nr
    !norm_vect = norm_vect + 0.5*(psi(i)*dconjg(psi(i)) + psi(i+1)*dconjg(psi(i+1)))*(r(i+1)-r(i));
    norm_vect = norm_vect + psi(i)*dconjg(psi(i))*h;
  enddo

  ! Here there is no need to multiply by r^2 since psi = phi/r
  norm_vect = dsqrt(norm_vect);

end function norm_vect



!%%%%%%%%%   Coupling coefficients
real*8  FUNCTION ang_c(l,m)

  IMPLICIT NONE;

  integer, intent(in) :: l,m;

  ang_c = ((real(l,8)+1.0d0)*(real(l,8)+1.0d0)-real(m,8)*real(m,8))/((2.0d0*real(l,8)+1.0d0)*(2.0d0*real(l,8)+3.0d0));
  ang_c = sqrt(ang_c);

end function ang_c




!%%%%%%%%%   subroutine to normalise psi
SUBROUTINE normalize(psi,r,Nr,Lmax)

  IMPLICIT NONE;
  integer, intent(in) :: Nr,Lmax;
  Complex*16,dimension(:,:),Allocatable :: psi;
  real*8,dimension(:),Allocatable,intent(in) :: r;

  real*8 :: norme;

  norme = norm(psi,r,Nr,Lmax);

  psi = psi/norme; ! Here it is a matrix


END SUBROUTINE normalize




SUBROUTINE tridiag_gauss(a,b,c,r,u,Nmax)
! This routine solves the system M*u = r where M is a tri-diagonal matrix filled up
! with elements a(2:Nmax) (linf), b(1:Nmax)(diag) and c(1:Nmax-1) (lsup)

  INTEGER j;
  INTEGER , intent(in) :: Nmax;
  complex*16,dimension(:), intent(in) :: a,b,c,r;
  complex*16,dimension(:) :: u;
  complex*16 :: bet;
  complex*16,dimension(Nmax) :: gam;


  if((real(b(1)).eq.0).and.(aimag(b(1)).eq.0)) pause 'tridiag : rewrite the equation -> b(1) = 0.0';

  bet = b(1);
  u(1) = r(1)/bet;

  gam(:) = 0.0d0;

  do j = 2 , Nmax                        ! Forward substitution
    gam(j) = c(j-1)/bet;
    bet = b(j) - a(j)*gam(j);
    if((real(bet).eq.0).and.(aimag(bet).eq.0)) pause 'tridiag failed';
    u(j) = (r(j) - a(j)*u(j-1))/bet;
  enddo

  do j = Nmax-1 , 1 , -1                        ! Backward substitution
    u(j) = u(j)-gam(j+1)*u(j+1);
  enddo

END SUBROUTINE tridiag_gauss




!      SUBROUTINE solve_tridiag(a,b,c,v,x,n)     !http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
!      implicit none
!!      a - sub-diagonal (means it is the diagonal below the main diagonal)
!!      b - the main diagonal
!!      c - sup-diagonal (means it is the diagonal above the main diagonal)
!!      v - right part
!!      x - the answer
!!      n - number of equations

!        integer,intent(in) :: n
!        real(8),dimension(n),intent(in) :: a,b,c,v
!        real(8),dimension(n),intent(out) :: x
!        real(8),dimension(n) :: bp,vp
!        real(8) :: m
!        integer i

!! Make copies of the b and v variables so that they are unaltered by this sub
!        bp(1) = b(1)
!        vp(1) = v(1)

!        !The first pass (setting coefficients):
!firstpass: do i = 2,n
!         m = a(i)/bp(i-1)
!         bp(i) = b(i) - m*c(i-1)
!         vp(i) = v(i) - m*vp(i-1)
!        end do firstpass

!         x(n) = vp(n)/bp(n)
!        !The second pass (back-substition)
!backsub:do i = n-1, 1, -1
!          x(i) = (vp(i) - c(i)*x(i+1))/bp(i)
!        end do backsub

!    end SUBROUTINE solve_tridiag



!SUBROUTINE tridiag_gauss(a,b,c,r,u,Nmax)
!! This routine solves the system M*u = r where M is a tri-diagonal matrix filled up
!! with elements a(2:Nmax) (linf), b(1:Nmax)(diag) and c(1:Nmax-1) (lsup)
!
!  INTEGER j;
!  INTEGER , intent(in) :: Nmax;
!  complex*16,dimension(:), intent(in) :: a,b,c,r;
!  complex*16,dimension(:) :: u;
!  complex*16 :: bet;
!  complex*16,dimension(Nmax) :: gam;
!
!
!  if((real(b(1)).eq.0).and.(aimag(b(1)).eq.0)) pause 'tridiag : rewrite the equation -> b(1) = 0.0';
!
!  bet = b(1);
!  u(1) = r(1)/bet;
!
!  gam(:) = 0.0d0;
!
!  do j = 2 , Nmax                        ! Forward substitution
!    gam(j) = c(j-1)/bet;
!    bet = b(j) - a(j)*gam(j);
!    if((real(bet).eq.0).and.(aimag(bet).eq.0)) pause 'tridiag failed';
!    u(j) = (r(j) - a(j)*u(j-1))/bet;
!  enddo
!
!  do j = Nmax-1 , 1 , -1                        ! Backward substitution
!    u(j) = u(j)-gam(j+1)*u(j+1);
!  enddo
!
!END SUBROUTINE tridiag_gauss






end module tridiag
