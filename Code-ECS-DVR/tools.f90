!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! THIS MODULE CONTAINS THE FOLLOWING FUNCTIONS :
!!! Init_data() :
!!! Init_matrix() :
!!! plgndr(l,m,x) :
!!! Y_sp(l,m,theta,phi) :
!!! angle_2pi(x,y) :
!!! plot_projectxy(xmin,xmax,Nx,ymin,ymax,Ny,r,Nr,Lmax,psi) :
!!! plot_projectxz(xmin,xmax,Nx,zmin,zmax,Nz,r,Nr,Lmax,psi) :
!!! index_radius(radius,r,Nr) :


module tools

use util; use tridiag; use potential; use coulombwave;

IMPLICIT NONE;

CONTAINS




!%%%%%%%%%%  Data initialisation
SUBROUTINE Init_data()

  Integer :: st;
  !character(1) :: tab;
  !tab(1) = char(9);

  real*8 :: t_tot;

  write(*,*)

  open(UNIT=1,FILE="param.txt",FORM="FORMATTED", ACCESS="SEQUENTIAL",status='old',action='read');

  write(*,*) "Calculation of the TDSE 1D in the velocity gauge using DVR basis";
  write(*,*) "(c) Harmonic Group, CELIA, CNRS";
  write(*,*)
  write(*,*) "Parameters of the calculation extracted from param.txt";
  write(*,*)


  READ(UNIT=1,fmt=*,IOSTAT=st) E_guess;
  write(*,10) E_guess; !Here char(9) in horizontal tab
10  format ("Initial energy to calculate the exact one (E_guess) :",T60,f13.6," a.u.");

  READ(UNIT=1,fmt=*,IOSTAT=st) xmax;
  write(*,40) xmax;
40  format ("Size of the box (xmax) : ",T60,f13.6," a.u.");

  READ(UNIT=1,fmt=*,IOSTAT=st) Nseq;
  write(*,50) Nseq;
50  format ("Number of points defining the sequence from 0 to Rmax (Nseq) : ",T60,i5);

  READ(UNIT=1,fmt=*,IOSTAT=st) NDVR;
  write(*,51) NDVR;
51  format ("Order of the polynome defined in the sequence (NDVR) : ",T60,i5);

  READ(UNIT=1,fmt=*,IOSTAT=st) x_ecs;
  write(*,20) x_ecs;
20  format ("Position at which the ECS starts : ",T60,f13.6," a.u.");

  READ(UNIT=1,fmt=*,IOSTAT=st) theta_ecs;
  write(*,30) theta_ecs;
30  format ("Theta of the ECS : ",T60,f13.6);

  READ(UNIT=1,fmt=*,IOSTAT=st) Nexp;
  write(*,60) Nexp;
60  format ("Number of points of the expansion of the grid (Nexp): ",T60,i5);

  READ(UNIT=1,fmt=*,IOSTAT=st) omega;
  write(*,70) omega;
70  format ("Frequency of the laser (omega) : ",T60,f13.6," a.u.");

  READ(UNIT=1,fmt=*,IOSTAT=st) omega2;
  write(*,170) omega2;
170  format ("Frequency of the second pulse (omega2) : ",T60,f13.6," a.u.");

  READ(UNIT=1,fmt=*,IOSTAT=st) F0;
  write(*,80) F0;
80  format ("Amplitude of the field (F0) : ",T60,f13.6," a.u.");

  READ(UNIT=1,fmt=*,IOSTAT=st) F02;
  write(*,180) F02;
180  format ("Amplitude of the second field (F02) : ",T60,f13.6," a.u.");

  READ(UNIT=1,fmt=*,IOSTAT=st) dt;
  write(*,100) dt;
100  format ("Time step : ",T60,f13.6," a.u.");

  READ(UNIT=1,fmt=*,IOSTAT=st) Nc;
  write(*,101) Nc;
101  format ("Number of cycles in the pulse : ",T60,i5);

  READ(UNIT=1,fmt=*,IOSTAT=st) Nc2;
  write(*,201) Nc2;
201  format ("Number of cycles of the second pulse : ",T60,i5);

  !Nt = Nt*Nc;


  READ(UNIT=1,fmt=*,IOSTAT=st) phi_CEP;
  write(*,102) phi_CEP;
102  format ("Value od the CEP : ",T60,f12.6," in units of period for a sin² envelope");
  phi_CEP = (phi_CEP*(2.0d0*pi/omega))*omega;


  READ(UNIT=1,fmt=*,IOSTAT=st) phi_CEP2;
  write(*,202) phi_CEP2;
202  format ("Value of the CEP of the second field: ",T60,f13.6," in units of period for a sin² envelope");
  phi_CEP2 = (phi_CEP2*(2.0d0*pi/omega2))*omega2;


  READ(UNIT=1,fmt=*,IOSTAT=st) param_env;
  if ( param_env == 0) then
   write(*,*) "The envelope is 1";
  endif
  if ( param_env == 1) then
   write(*,*) "The envelope is a flat-top";
  endif
  if ( param_env == 2) then
   write(*,*) "The envelope is a sin²";
  endif


  READ(UNIT=1,fmt=*,IOSTAT=st) param_env2;
  if ( param_env2 == 0) then
   write(*,*) "The envelope of the second pulse is 1";
  endif
  if ( param_env2 == 1) then
   write(*,*) "The envelope of the second pulse is a flat-top";
  endif
  if ( param_env2 == 2) then
   write(*,*) "The envelope of the second pulse is a sin²";
  endif


  READ(UNIT=1,fmt=*,IOSTAT=st) t_on;
  write(*,103) t_on;
103  format ("Turn-on duration : ",T60,f13.6," in unit of period");

  READ(UNIT=1,fmt=*,IOSTAT=st) t_on2;
  write(*,203) t_on2;
203  format ("Turn-on duration of the second pulse : ",T60,f13.6," in unit of period");

  READ(UNIT=1,fmt=*,IOSTAT=st) t_off;
  write(*,104) t_off;
104  format ("Turn-off duration : ",T60,f13.6," in unit of period");

  READ(UNIT=1,fmt=*,IOSTAT=st) t_off2;
  write(*,204) t_off2;
204  format ("Turn-off duration of the second pulse : ",T60,f13.6," in unit of period");


  READ(UNIT=1,fmt=*,IOSTAT=st) t_delay;
  write(*,205) t_delay;
205  format ("Delay between the second and the first pulse : ",T60,f13.6," a.u.");

  t_tot = t_delay + 0.0d0*pi*(dfloat(Nc2)/omega2) + 2.0d0*pi*(dfloat(Nc)/omega);
  !t_tot = t_delay + 1.0d0*pi*(dfloat(Nc2)/omega2) + 1.0d0*pi*(dfloat(Nc)/omega);

  Nt = int(t_tot/dt);


  READ(UNIT=1,fmt=*,IOSTAT=st) E_start;
  write(*,105) E_start;
105  format ("Starting energy of the electron spectrum : ",T60,f13.6," a.u.");

  READ(UNIT=1,fmt=*,IOSTAT=st) E_step;
  write(*,106) E_step;
106  format ("Step of energy of the electron spectrum : ",T60,f13.6," a.u.");

  READ(UNIT=1,fmt=*,IOSTAT=st) dE;
  write(*,107) dE;
107  format ("Energy resolution of the electron spectrum : ",T60,f13.6," a.u.");

  READ(UNIT=1,fmt=*,IOSTAT=st) NE;
  write(*,108) NE;
108  format ("Number of points of the electron spectrum : ",T60,i5);

  READ(UNIT=1,fmt=*,IOSTAT=st) Nlanczos;
  write(*,120) Nlanczos;
120  format ("Size of the Lanczos basis : ",T60,i5);


  close(UNIT=1);

  write(*,*)


END SUBROUTINE Init_data




!%%%%%%%%%%  Data initialisation
SUBROUTINE Init_matrix()

  Integer :: i,j,k,l,m,index_PHP;
  real*8,dimension(:),Allocatable :: x_i,w_i;
  real*8 :: eps;
  complex*16 :: somme;

  Complex*16 :: exp_somme,ALPHA,BETA,inter_xi,inter_wi;

  Integer :: LDA,LWORK, KD, LDAB, INFO;

  CHARACTER :: JOBZ,UPLO,TRANSA,TRANSB;

  complex*16,dimension(:,:),Allocatable :: Tinter,D1_inter;


  COMPLEX*16,dimension(:,:),Allocatable :: H0_mat_format,eig_vect;
  COMPLEX*16,dimension(:,:),Allocatable :: mat_inter,H0_mat_diag,eig_vect_complex;
  COMPLEX*16,dimension(:),Allocatable :: x_inter,WORK,RWORK;


  Allocate(x_i(1:NDVR)); Allocate(w_i(1:NDVR));

  Allocate(Tinter(1:(Nseq)*(NDVR-1)+1,1:(Nseq)*(NDVR-1)+1),D1_inter(1:(Nseq)*(NDVR-1)+1,1:(Nseq)*(NDVR-1)+1))

  allocate(mat_inter(1:Nx,1:Nx))


  LDA = Nx;
  JOBZ = 'V'; !'N'
  UPLO = 'L';
  INFO = 0 ;
  LWORK = 2*LDA-1;
  KD = NDVR-1;
  LDAB = Nx
  Allocate(H0_mat_format(1:LDAB,1:Nx),x_inter(1:Nseq*(NDVR-1)),WORK(1:Nx*Nx),RWORK(1:3*LDA-2),eig_vect(1:Nx,1:Nx));
  Allocate(H0_mat_diag(1:Nx,1:Nx),eig_vect_complex(1:Nx,1:Nx));





  xseq(0) = -xmax;
  do i = 1 , Nseq
     if (i.lt.200) then
       xseq(i) = xseq(i-1) + 2.0d0*xmax/real(NSeq,8);
     else
       xseq(i) = xseq(i-1) + 2.0d0*xmax/real(NSeq,8);
     endif

  enddo

  write(*,*)
  write(*,*) "dx = ",xseq(2)-xseq(1)," in the inner region"
  write(*,*) "dx = ", xseq(301)-xseq(300)," in the outer region"
  write(*,*)



  eps = 1.0d-16;
  CALL lglnodes(x_i,w_i,NDVR,eps)


  do i = 1 , Nseq
     do j = 1 , NDVR

      xi_seq(i,j) = xseq(i-1) + (xseq(i)-xseq(i-1))*(x_i(j)+1.0d0)/2.0d0;
      wi_seq(i,j) = (xseq(i)-xseq(i-1))*w_i(j)/2.0d0;

     enddo
  enddo





            do i = 1 , Nseq
               do j = 1 , NDVR

                 inter_xi = xi_seq(i,j)
                 if(abs(inter_xi).ge.x_ecs) then

                   inter_wi = exp(xi*theta_ecs)*wi_seq(i,j);
                   if(abs(abs(inter_xi)-x_ecs).lt.(1.0d-8)) then !! This is the particular case of the Rc being at the edge of a sequence

                     if((real(inter_xi).lt.0 .AND. (j.EQ.NDVR)) .OR. (real(inter_xi).gt.0 .AND. (j.EQ.1))) then
                       inter_wi = exp(xi*theta_ecs)*wi_seq(i,j);
                     else
                       inter_wi = wi_seq(i,j);
                     endif

                   ENDIF



                  if(real(inter_xi).gt.0) then
                     inter_xi = x_ecs + exp(xi*theta_ecs)*(inter_xi - x_ecs);
                   else
                     inter_xi = -x_ecs + exp(xi*theta_ecs)*(inter_xi + x_ecs);
                   endif

                 else
                   inter_xi = xi_seq(i,j)
                   inter_wi = wi_seq(i,j)
                 endif

                 xi_seq(i,j) = inter_xi
                 wi_seq(i,j) = inter_wi

               enddo
            enddo






  DO i = 1 , Nseq

     DO j = 1 , NDVR

         DO k = 1 , NDVR
            somme = 0.0d0;
            DO l = 1 , NDVR
               somme = somme + wi_seq(i,l)*der_lagrange(i,j,l)*der_lagrange(i,k,l)
            ENDDO

            T_reduced(i,j,k) = somme;
            D1_reduced(i,j,k) = wi_seq(i,j)*der_lagrange(i,k,j);

         ENDDO

     ENDDO


  ENDDO







! write(*,*) "TEST ESSAI"
! DO i = 1 , NSeq
!   write(*,*) " i = ",i
!   write(*,*)
!   write(*,*) xi_seq(i,:)
!   write(*,*)
!   write(*,*) wi_seq(i,:)
! ENDDO


!write(*,*) "TESt Kinetic MATRIX"
!write(*,*) "Nseq = 278"
!write(*,*) xi_seq(278,:)
!write(*,*) "Nseq = 223"
!write(*,*) xi_seq(223,:)
!write(*,*) "Nseq = 5"
!write(*,*) wi_seq(278,:)
!write(*,*) "Nseq = 496"
!write(*,*) wi_seq(223,:)
!write(*,*) "Nseq = 277"
!write(*,*) T_reduced(277,8,9)
!write(*,*) "Nseq = 223"
!write(*,*) T_reduced(222,3,4)
!write(*,*) "Full"
!DO i = 1 , NDVR
!  DO j = 1 , NDVR
!   write(*,*) i,j, T_reduced(277,i,j)
!  ENDDO
! ENDDO
!
! write(*,*) "Full"
! DO i = 1 , NDVR
!   DO j = 1 , NDVR
!    write(*,*) i,j, T_reduced(222,i,j)
!   ENDDO
!  ENDDO



k = 2; m = 1;
do i = 1 , Nx


   x(i) = xi_seq(m,k)
   weight(i) = wi_seq(m,k)
   x_rb(i) = xseq(m-1) + (xseq(m)-xseq(m-1))*(x_i(k)+1.0d0)/2.0d0;

    if (k .eq. (NDVR))then
        if (m .eq. Nseq) then
         weight(i) = wi_seq(m,k)
         write(*,*) "WRONG USE !!!!! ----- !!!!"
        else
         weight(i) = wi_seq(m,NDVR)+wi_seq(m+1,1)
        endif
    endif
    !weight_rb(i) = weight(i);

    if( k .lt. NDVR) then
        k = k + 1;
    else
        m = m +1; k = 2;
    endif

    !write(*,*) Nx,i,m,k,weight(i),wi_seq(m,k-1)

enddo


! do m = 1 , Nseq
!  do  k = 1, NDVR-1
!
!   i = (m-1)*(NDVR-1) + k;
!
!   x_inter(i) = xi_seq(m,k)
!   weight(i) = wi_seq(m,k)
!   x_rb(i) = xseq(m-1) + (xseq(m)-xseq(m-1))*(xi(k)+1.0d0)/2.0d0;
!
!    if (k .eq. (NDVR-1))then
!         weight(i) = wi_seq(m,NDVR)+wi_seq(m+1,1)
!    else
!      weight(i) = wi_seq(m,k)
!    endif
!
!    !write(*,*) i,x_inter(i)
!
!enddo
!enddo



 Tinter(:,:) = (0.0d0,0.0d0); D1_inter(:,:) = (0.0d0,0.0d0);
 do i = 2 , Nseq - 1

    do j = 1 , NDVR

       do k = 1 , NDVR

            if ((j .NE. 1 ) .AND. (k .NE. 1) .AND. (j .NE. NDVR) .and. (k .NE. NDVR)) then
                Tinter((i-1)*(NDVR-1)+j,(i-1)*(NDVR-1) + k) = 0.5d0*T_reduced(i,j,k)/sqrt(wi_seq(i,j)*wi_seq(i,k))
                D1_inter((i-1)*(NDVR-1)+j,(i-1)*(NDVR-1)+k) = D1_reduced(i,j,k)/sqrt(wi_seq(i,j)*wi_seq(i,k))
            endif

            if ((j .EQ. 1) .AND. ( k .EQ. 1)) then
                Tinter((i-1)*(NDVR-1)+j,(i-1)*(NDVR-1)+ k) = &
                0.5d0*(T_reduced(i,1,1)+T_reduced(i-1,NDVR,NDVR))/(wi_seq(i,1) + wi_seq(i-1,NDVR))
                D1_inter((i-1)*(NDVR-1)+j,(i-1)*(NDVR-1)+k) = &
                (D1_reduced(i,1,1)+D1_reduced(i-1,NDVR,NDVR))/(wi_seq(i,1) + wi_seq(i-1,NDVR))
            endif

            if ((j .EQ. NDVR) .AND. ( k .EQ. NDVR)) then
                Tinter((i-1)*(NDVR-1)+j,(i-1)*(NDVR-1)+ k) = &
                0.5d0*(T_reduced(i,NDVR,NDVR)+T_reduced(i+1,1,1))/(wi_seq(i,NDVR) + wi_seq(i+1,1))
                D1_inter((i-1)*(NDVR-1)+j,(i-1)*(NDVR-1)+k) = &
                (D1_reduced(i,NDVR,NDVR)+D1_reduced(i+1,1,1))/(wi_seq(i,NDVR) + wi_seq(i+1,1))
            endif

            if ((j .EQ. 1) .AND. ( k .NE. 1) .AND. (k .NE. NDVR)) then
               Tinter((i-1)*(NDVR-1)+j,(i-1)*(NDVR-1)+ k) = &
               0.5d0*(T_reduced(i,1,k))/sqrt((wi_seq(i,k))*(wi_seq(i,1) + wi_seq(i-1,NDVR)))
               D1_inter((i-1)*(NDVR-1)+j,(i-1)*(NDVR-1)+k) = &
               (D1_reduced(i,1,k))/sqrt((wi_seq(i,k))*(wi_seq(i,1) + wi_seq(i-1,NDVR)))
            endif

            if ((j .EQ. 1) .AND. ( k .EQ. NDVR)) then
                Tinter((i-1)*(NDVR-1)+j,(i-1)*(NDVR-1)+ k) = &
                0.5d0*(T_reduced(i,1,NDVR))/sqrt((wi_seq(i,NDVR)+wi_seq(i+1,1))*(wi_seq(i,1) + wi_seq(i-1,NDVR)))
                D1_inter((i-1)*(NDVR-1)+j,(i-1)*(NDVR-1)+k) = &
                (D1_reduced(i,1,NDVR))/sqrt((wi_seq(i,NDVR)+wi_seq(i+1,1))*(wi_seq(i,1) + wi_seq(i-1,NDVR)))
            endif

            if ((j .NE. 1) .AND. ( k .EQ. 1) .AND. (j .NE. NDVR)) then
                Tinter((i-1)*(NDVR-1)+j,(i-1)*(NDVR-1)+ k) = &
                0.5d0*(T_reduced(i,j,1))/sqrt((wi_seq(i,j))*(wi_seq(i,1) + wi_seq(i-1,NDVR)))
                D1_inter((i-1)*(NDVR-1)+j,(i-1)*(NDVR-1)+k) = &
                (D1_reduced(i,j,1))/sqrt((wi_seq(i,j))*(wi_seq(i,1) + wi_seq(i-1,NDVR)))
            endif


            if ((j .EQ. NDVR) .and. ( k .eq. 1)) then
                Tinter((i-1)*(NDVR-1)+j,(i-1)*(NDVR-1)+ k) = &
                0.5d0*(T_reduced(i,NDVR,1))/sqrt((wi_seq(i,NDVR)+wi_seq(i+1,1))*(wi_seq(i,1) + wi_seq(i-1,NDVR)))
                D1_inter((i-1)*(NDVR-1)+j,(i-1)*(NDVR-1)+k) = &
                (D1_reduced(i,NDVR,1))/sqrt((wi_seq(i,NDVR)+wi_seq(i+1,1))*(wi_seq(i,1) + wi_seq(i-1,NDVR)))
            endif

            if ((j .eq. NDVR) .and. (k .ne. NDVR) .and. ( k .ne. 1)) then
                Tinter((i-1)*(NDVR-1)+j,(i-1)*(NDVR-1)+ k) = &
                0.5d0*(T_reduced(i,NDVR,k))/sqrt((wi_seq(i,NDVR)+wi_seq(i+1,1))*wi_seq(i,k))
                D1_inter((i-1)*(NDVR-1)+j,(i-1)*(NDVR-1)+k) = &
                (D1_reduced(i,NDVR,k))/sqrt((wi_seq(i,NDVR)+wi_seq(i+1,1))*wi_seq(i,k))
            endif

            if ((j .ne. 1) .and. (j .ne. NDVR) .and. ( k .eq. NDVR)) then
                Tinter((i-1)*(NDVR-1)+j,(i-1)*(NDVR-1)+ k) = &
                0.5d0*(T_reduced(i,j,NDVR))/sqrt(wi_seq(i,j)*(wi_seq(i+1,1) + wi_seq(i,NDVR)))
                D1_inter((i-1)*(NDVR-1)+j,(i-1)*(NDVR-1)+k) = &
                (D1_reduced(i,j,NDVR))/sqrt(wi_seq(i,j)*(wi_seq(i+1,1) + wi_seq(i,NDVR)))
            endif


       enddo

    enddo

 enddo

 !!! Case for i = 1

     do j = 1 , NDVR

       do k = 1 , NDVR

            if ((j .NE. 1 ) .AND. (k .NE. 1) .AND. (j .NE. NDVR) .and. (k .NE. NDVR)) then
                Tinter(j,k) = 0.5d0*T_reduced(1,j,k)/sqrt(wi_seq(1,j)*wi_seq(1,k))
                D1_inter(j,k) = D1_reduced(1,j,k)/sqrt(wi_seq(1,j)*wi_seq(1,k))
            endif

            if ((j .EQ. 1) .AND. ( k .EQ. 1)) then
                Tinter(j,k) = &
                0.5d0*(T_reduced(1,1,1))/(wi_seq(1,1))
                D1_inter(j,k) = &
                (D1_reduced(1,1,1))/(wi_seq(1,1))
            endif

            if ((j .EQ. NDVR) .AND. ( k .EQ. NDVR)) then
                Tinter(j,k) = &
                0.5d0*(T_reduced(1,NDVR,NDVR)+T_reduced(2,1,1))/(wi_seq(1,NDVR) + wi_seq(2,1))
                D1_inter(j,k) = &
                (D1_reduced(1,NDVR,NDVR)+D1_reduced(2,1,1))/(wi_seq(1,NDVR) + wi_seq(2,1))
            endif

            if ((j .EQ. 1) .AND. ( k .NE. 1) .AND. (k .NE. NDVR)) then
               Tinter(j,k) = &
               0.5d0*(T_reduced(1,1,k))/sqrt((wi_seq(1,k))*(wi_seq(1,1)))
               D1_inter(j,k) = &
               (D1_reduced(1,1,k))/sqrt((wi_seq(1,k))*(wi_seq(1,1)))
            endif

            if ((j .EQ. 1) .AND. ( k .EQ. NDVR)) then
                Tinter(j,k) = &
                0.5d0*(T_reduced(1,1,NDVR))/sqrt((wi_seq(1,NDVR)+wi_seq(2,1))*(wi_seq(1,1)))
                D1_inter(j,k) = &
                (D1_reduced(1,1,NDVR))/sqrt((wi_seq(1,NDVR)+wi_seq(2,1))*(wi_seq(1,1)))
            endif

            if ((j .NE. 1) .AND. ( k .EQ. 1) .AND. (j .NE. NDVR)) then
                Tinter(j,k) = &
                0.5d0*(T_reduced(1,j,1))/sqrt((wi_seq(1,j))*(wi_seq(1,1)))
                D1_inter(j,k) = &
                (D1_reduced(1,j,1))/sqrt((wi_seq(1,j))*(wi_seq(1,1)))
            endif

            if ((j .EQ. NDVR) .and. ( k .eq. 1)) then
                Tinter(j,k) = &
                0.5d0*(T_reduced(1,NDVR,1))/sqrt((wi_seq(1,NDVR)+wi_seq(2,1))*(wi_seq(1,1)))
                D1_inter(j,k) = &
                (D1_reduced(1,NDVR,1))/sqrt((wi_seq(1,NDVR)+wi_seq(2,1))*(wi_seq(1,1)))
            endif

            if ((j .eq. NDVR) .and. (k .ne. NDVR) .and. ( k .ne. 1)) then
                Tinter(j,k) = &
                0.5d0*(T_reduced(1,NDVR,k))/sqrt((wi_seq(1,NDVR)+wi_seq(2,1))*wi_seq(1,k))
                D1_inter(j,k) = &
                (D1_reduced(1,NDVR,k))/sqrt((wi_seq(1,NDVR)+wi_seq(2,1))*wi_seq(1,k))
            endif

            if ((j .ne. 1) .and. (j .ne. NDVR) .and. ( k .eq. NDVR)) then
                Tinter(j,k) = &
                0.5d0*(T_reduced(1,j,NDVR))/sqrt(wi_seq(1,j)*(wi_seq(2,1) + wi_seq(1,NDVR)))
                D1_inter(j,k) = &
                (D1_reduced(1,j,NDVR))/sqrt(wi_seq(1,j)*(wi_seq(2,1) + wi_seq(1,NDVR)))
            endif


       enddo

    enddo




 !!! case for i = Nseq

    do j = 1 , NDVR

       do k = 1 , NDVR

            if ((j .NE. 1 ) .AND. (k .NE. 1) .AND. (j .NE. NDVR) .and. (k .NE. NDVR)) then
                Tinter((Nseq-1)*(NDVR-1)+j,(Nseq-1)*(NDVR-1)+ k) = 0.5d0*T_reduced(Nseq,j,k)/sqrt(wi_seq(Nseq,j)*wi_seq(Nseq,k))
                D1_inter((Nseq-1)*(NDVR-1)+j,(Nseq-1)*(NDVR-1)+ k) = D1_reduced(Nseq,j,k)/sqrt(wi_seq(Nseq,j)*wi_seq(Nseq,k))
            endif

            if ((j .EQ. 1) .AND. ( k .EQ. 1)) then
                Tinter((Nseq-1)*(NDVR-1)+j,(Nseq-1)*(NDVR-1)+ k) = &
                0.5d0*(T_reduced(Nseq,1,1)+T_reduced(Nseq-1,NDVR,NDVR))/(wi_seq(Nseq,1) + wi_seq(Nseq-1,NDVR))
                D1_inter((Nseq-1)*(NDVR-1)+j,(Nseq-1)*(NDVR-1)+k) = &
                (D1_reduced(Nseq,1,1)+D1_reduced(Nseq-1,NDVR,NDVR))/(wi_seq(Nseq,1) + wi_seq(Nseq-1,NDVR))
            endif

            if ((j .EQ. NDVR) .AND. ( k .EQ. NDVR)) then
                Tinter((Nseq-1)*(NDVR-1)+j,(Nseq-1)*(NDVR-1)+ k) = &
                0.5d0*(T_reduced(Nseq,NDVR,NDVR))/(wi_seq(Nseq,NDVR))
                D1_inter((Nseq-1)*(NDVR-1)+j,(Nseq-1)*(NDVR-1)+k) = &
                (D1_reduced(Nseq,NDVR,NDVR))/(wi_seq(Nseq,NDVR))
            endif

            if ((j .EQ. 1) .AND. ( k .NE. 1) .AND. (k .NE. NDVR)) then
               Tinter((Nseq-1)*(NDVR-1)+j,(Nseq-1)*(NDVR-1)+ k) = &
               0.5d0*(T_reduced(Nseq,1,k))/sqrt((wi_seq(Nseq,k))*(wi_seq(Nseq,1) + wi_seq(Nseq-1,NDVR)))
               D1_inter((Nseq-1)*(NDVR-1)+j,(Nseq-1)*(NDVR-1)+k) = &
               (D1_reduced(Nseq,1,k))/sqrt((wi_seq(Nseq,k))*(wi_seq(Nseq,1) + wi_seq(Nseq-1,NDVR)))
            endif

            if ((j .EQ. 1) .AND. ( k .EQ. NDVR)) then
                Tinter((Nseq-1)*(NDVR-1)+j,(Nseq-1)*(NDVR-1)+ k) = &
                0.5d0*(T_reduced(Nseq,1,NDVR))/sqrt((wi_seq(Nseq,NDVR))*(wi_seq(Nseq,1) + wi_seq(Nseq-1,NDVR)))
                D1_inter((Nseq-1)*(NDVR-1)+j,(Nseq-1)*(NDVR-1)+k) = &
                (D1_reduced(Nseq,1,NDVR))/sqrt((wi_seq(Nseq,NDVR))*(wi_seq(Nseq,1) + wi_seq(Nseq-1,NDVR)))
            endif

            if ((j .NE. 1) .AND. ( k .EQ. 1) .AND. (j .NE. NDVR)) then
                Tinter((Nseq-1)*(NDVR-1)+j,(Nseq-1)*(NDVR-1)+ k) = &
                0.5d0*(T_reduced(Nseq,j,1))/sqrt((wi_seq(Nseq,j))*(wi_seq(Nseq,1) + wi_seq(Nseq-1,NDVR)))
                D1_inter((Nseq-1)*(NDVR-1)+j,(Nseq-1)*(NDVR-1)+k) = &
                (D1_reduced(Nseq,j,1))/sqrt((wi_seq(Nseq,j))*(wi_seq(Nseq,1) + wi_seq(Nseq-1,NDVR)))
            endif


            if ((j .EQ. NDVR) .and. ( k .eq. 1)) then
                Tinter((Nseq-1)*(NDVR-1)+j,(Nseq-1)*(NDVR-1)+ k) = &
                0.5d0*(T_reduced(Nseq,NDVR,1))/sqrt((wi_seq(Nseq,NDVR))*(wi_seq(Nseq,1) + wi_seq(Nseq-1,NDVR)))
                D1_inter((Nseq-1)*(NDVR-1)+j,(Nseq-1)*(NDVR-1)+k) = &
                (D1_reduced(Nseq,NDVR,1))/sqrt((wi_seq(Nseq,NDVR))*(wi_seq(Nseq,1) + wi_seq(Nseq-1,NDVR)))
            endif

            if ((j .eq. NDVR) .and. (k .ne. NDVR) .and. ( k .ne. 1)) then
                Tinter((Nseq-1)*(NDVR-1)+j,(Nseq-1)*(NDVR-1)+ k) = &
                0.5d0*(T_reduced(Nseq,NDVR,k))/sqrt((wi_seq(Nseq,NDVR))*wi_seq(Nseq,k))
                D1_inter((Nseq-1)*(NDVR-1)+j,(Nseq-1)*(NDVR-1)+k) = &
                (D1_reduced(Nseq,NDVR,k))/sqrt((wi_seq(Nseq,NDVR))*wi_seq(Nseq,k))
            endif

            if ((j .ne. 1) .and. (j .ne. NDVR) .and. ( k .eq. NDVR)) then
                Tinter((Nseq-1)*(NDVR-1)+j,(Nseq-1)*(NDVR-1)+ k) = &
                0.5d0*(T_reduced(Nseq,j,NDVR))/sqrt(wi_seq(Nseq,j)*(wi_seq(Nseq,NDVR)))
                D1_inter((i-1)*(NDVR-1)+j,(i-1)*(NDVR-1)+k) = &
                (D1_reduced(Nseq,j,NDVR))/sqrt(wi_seq(Nseq,j)*(wi_seq(Nseq,NDVR)))
            endif


       enddo

    enddo



   do i = 1 , Nx
      do j = 1,Nx
          T_mat(i,j) = Tinter(i+1,j+1)
          D1_mat(i,j) = D1_inter(i+1,j+1)
      enddo
      !x(i) = x_inter(i+1);
   enddo



 !! USed in E_initialise
  H0_mat(:,:) = T_mat(:,:);


  do i = 1 , Nx
     Pot_mat(i) = pot(x(i));
     H0_mat(i,i) = H0_mat(i,i) + Pot_mat(i);
  enddo

 write(*,*)
 !write(*,*) H0_mat(5,5)
 write(*,*) "Test matrix symmetry, Nx = ",Nx
 write(*,*)

 open(UNIT=1,FILE="Hamiltonian_re.dat",FORM="FORMATTED",action='write');
 open(UNIT=2,FILE="Hamiltonian_im.dat",FORM="FORMATTED",action='write');

 do i = 1 , Nx

     do j = 1,Nx
         if((abs(H0_mat(i,j)-H0_mat(Nx-i+1,Nx-j+1)).gt.1e-8).OR.(abs(H0_mat(i,j)-H0_mat(j,i)).gt.1e-8)) then
           write(*,*) "ERROR SYMMETRY (H0_mat): ", i,j
           write(*,*) H0_mat(i,j)
           write(*,*) H0_mat(Nx-i+1,Nx-j+1)
           write(*,*) H0_mat(j,i)

         endif

         if(j.ne.Nx) then
          write(UNIT=1,fmt='(XE13.6)',advance='no') real(H0_mat(i,j));
          write(UNIT=2,fmt='(XE13.6)',advance='no') aimag(H0_mat(i,j));
         else
          write(UNIT=1,fmt='(XE13.6)') real(H0_mat(i,j));
          write(UNIT=2,fmt='(XE13.6)') aimag(H0_mat(i,j));
        endif

     enddo

      if((abs(x(i)+x(Nx+1-i)).gt.1e-8).OR.(abs(weight(i)-weight(Nx+1-i)).gt.1e-8)) then
        write(*,*) "ERROR SYMMETRY (w,x): ", i
        write(*,*) x(i),xi_seq(1,2),xi_seq(Nseq,NDVR-1)
        write(*,*) x(Nx+1-i)
        write(*,*) weight(i),wi_seq(1,2),wi_seq(Nseq,NDVR-1)
        write(*,*) weight(Nx+1-i)
      endif

    enddo

 close(1)
 close(2)
 deallocate(x_i,w_i,D1_inter,Tinter,H0_mat_format,x_inter)





END SUBROUTINE Init_matrix







!%%%%%% Zero et poids de Gauss-Lobatto
SUBROUTINE lglnodes(xi,wi,NDVR,eps)

  integer, intent(in) :: NDVR;
  real*8,  intent(in) :: eps;
  real*8, dimension(:),Allocatable :: xi,wi,xold;
  real*8, dimension(:,:),Allocatable :: P;

  integer :: i,j,k;
  real*8 :: epss,test; !fact,pll,pmm,pmmp1,somx2,coef,n1,n2;
  Allocate(xold(1:NDVR)); Allocate(P(1:NDVR,0:NDVR));


    wi(:) = 0.0d0;
    xi(:) = 0.0d0;
    xold(:) = 0.0d0;

    ! The Legendre Vandermonde Matrix
    P(:,:) = 0.0d0;

    epss = eps

    ! Use the Chebyshev-Gauss-Lobatto nodes as the first guess
    do i = 1 , NDVR
        xi(i) = -cos(acos(-1.0)*real(i-1,8)/real(NDVR-1,8))
    enddo

    ! Compute P using the recursion relation
    ! Compute its first and second derivatives and
    ! update x using the Newton-Raphson method.

    xold(:) = 2.0

    do i = 0 , 10
        xold(:) = xi(:)

        P(:,0) = 1.0
        P(:,1) = xi(:)

        do k = 2 , NDVR
           do j = 1 , NDVR
            P(j,k) = ( real(2*k-1,8)*xi(j)*P(j,k-1) - real(k-1,8)*P(j,k-2) ) / real(k,8)
           enddo
        enddo


        do j = 1 , NDVR
          xi(j) = xold(j) - ( xi(j)*P(j,NDVR-1) - P(j,NDVR-2) )/( real(NDVR,8)*P(j,NDVR-1))
        enddo

        test = 0.0d0;
        do j = 1 , NDVR
           if((xi(j)-xold(j))**(2.0d0) > test) then
              test = (xi(j)-xold(j))**(2.0d0)
           endif

        enddo
        write(*,*) i,test;
        if (test < epss ) then
            go to 10
        endif

    enddo


10  do i = 1 , NDVR
      wi(i) = 2.0d0 / ( real(NDVR*(NDVR-1),8)*(P(i,NDVR-1)*P(i,NDVR-1)));
    enddo



    deallocate(P,xold)


END SUBROUTINE



!%%%%%% Legendre polynome
complex*16 FUNCTION lagrange_pol(i,j,x)
   integer, intent(in) :: i,j;
   complex*16, intent(in) :: x;

   integer :: k;

   lagrange_pol = 0.0d0

    if (x.ne.xi_seq(i,j)) then
        lagrange_pol = 1.0d0
        do k = 1 , NDVR
            if (k .ne. j) then
                lagrange_pol = lagrange_pol*(x-xi_seq(i,k))/(xi_seq(i,j)-xi_seq(i,k))
            endif
        enddo

    else

        lagrange_pol = 1.0d0;

    endif


END FUNCTION




!%%%%%% Legendre polynome
complex*16 FUNCTION der_lagrange(i,j,k)
   integer, intent(in) :: i,j,k;
   integer :: m;

   der_lagrange = 0.0d0;

    if (k.eq.j) then
        der_lagrange = 0.0d0
        do m = 1 , NDVR
            if (m.ne.j) then
                der_lagrange = der_lagrange + 1./(xi_seq(i,j)-xi_seq(i,m))
            endif
        enddo
    else
        der_lagrange = 1.0d0
        do m = 1 , NDVR
            if ((m.ne.j).AND.(m.ne.k)) then
                der_lagrange = der_lagrange*(xi_seq(i,k)-xi_seq(i,m))/(xi_seq(i,j)-xi_seq(i,m));
            endif
        enddo
        der_lagrange = der_lagrange/(xi_seq(i,j) - xi_seq(i,k));
    endif



END FUNCTION





!%%%%%% Legendre polynome
real*8 FUNCTION plgndr(l,m,x)

  integer, intent(in) :: m,l;
  real*8, intent(in) :: x;
  integer i,ll,mm;
  real*8 fact,pll,pmm,pmmp1,somx2,coef,n1,n2;

  mm = m; coef = 1.0;
  if(abs(mm).gt.l.or.abs(x).gt.1) pause 'bad argument in plgndr'
  ! case m < 0
  if(m.lt.0) then
    coef = (-1.0d0)**float(mm);
    mm = -mm;

    n1 = 0.0d0;
    do i = 1 , l-mm
      n1 = n1+ dlog(dfloat(i));
    enddo
    n2 = 0.0d0;
    do i = 1 , l+mm
      n2 = n2*dlog(dfloat(i));
    enddo
    coef = coef*dexp(n1-n2);

  endif

  pmm = 1.0d0;
  if (mm.gt.0) then
     somx2 = sqrt((1.0d0-x)*(1.0d0+x));
     fact = 1.0d0;
     do i = 1,mm
        pmm = -pmm*fact*somx2;
	fact = fact + 2.0d0;
     enddo
  endif

  if (l == mm) then
	plgndr=pmm*coef;
  else
	pmmp1=x*(2*mm+1)*pmm
	if (l == mm+1) then
		plgndr=pmmp1*coef;
	else
 	  do ll=mm+2,l
	    pll=(x*(2*ll-1)*pmmp1-(ll+mm-1)*pmm)/(ll-mm);
	    pmm=pmmp1;
	    pmmp1=pll;
	  end do
	plgndr=pll*coef;
	endif
  endif

  !write(*,*) "legendre", plgndr

END FUNCTION



!%%%%%% Spherical harmonics
complex*16 FUNCTION Y_sp(l,m,theta,phi)

  integer, intent(in) :: m,l;
  real*8, intent(in) :: theta,phi;
  integer i;
  real*8 coef,x,n1,n2;

    coef = dfloat(2*l+1)/(4.0d0*pi);
    n1 = 0.0d0;
    do i = 1 , l-m
      n1 = n1 + dlog(dfloat(i));
    enddo
    n2 = 0.0d0;
    do i = 1 , l+m
      n2 = n2 + dlog(dfloat(i));
    enddo
    coef = coef*dexp(n1-n2); coef = dsqrt(coef);

    x = dcos(theta);
    Y_sp = coef*plgndr(l,m,x)*(dcos(dfloat(m)*phi)+xi*dsin(dfloat(m)*phi));


END FUNCTION


!%%%%%%%%%%% Calculate the angle from 0 to 2pi
real*8 FUNCTION angle_2pi(x,y)
   real*8, intent(in) :: x,y;

   if(x.eq.0) then
     angle_2pi = 0.5*pi;
       if(y.lt.0) then
         angle_2pi = 1.5*pi;
       endif
   else
     angle_2pi = atan(abs(y/x));

     if(x.lt.0.and.y.lt.0) then
       angle_2pi = pi+angle_2pi;
     endif
     if(x.lt.0.and.0.lt.y) then
       angle_2pi = pi-angle_2pi;
     endif
     if(0.lt.x.and.y.lt.0) then
       angle_2pi = 2.0*pi-angle_2pi;
     endif
   endif

   if(y.eq.0.and.x.lt.0) then
       angle_2pi = pi;
   endif
   if(y.eq.0.and.0.le.x) then
       angle_2pi = 0.0d0;
   endif

END FUNCTION


!%%%%%%%%%%% Calculate the angle from 0 to 2pi
INTEGER FUNCTION index_radius(radius,x,Nx)
   real*8, intent(in) :: radius;
   real*8,dimension(:),Allocatable :: x;
   integer, intent(in) :: Nx;

   if (radius.lt.x(Nx)) then
     index_radius = 1;
     do while (x(index_radius).lt.radius)
       index_radius = index_radius+1;
     enddo
     index_radius = index_radius - 1;
   else
     index_radius = -1;
   endif


END FUNCTION







!%%%%%%%%%%% Calculate the energy associated to H of the wave-function psi
complex*16 FUNCTION E_calculation(psi,x,Nx)
   complex*16,dimension(:),Allocatable, intent(in) :: x;
   integer, intent(in) :: Nx;
   Complex*16,dimension(:), Allocatable, intent(in) :: psi;

   complex*16 :: E_coul,E_cin,E_cent,coeff;
   integer :: i,l,k,n,index_set;
   Complex*16,dimension(:),allocatable :: u;

   Allocate(u(1:Nx));


   E_coul = 0.0d0; E_cin = 0.0d0;

  do i = 1,Nx
      ! Calculation of the coulomb component
	   E_coul = E_coul + psi(i)*dconjg(psi(i))*pot(x(i));
  enddo

        !write(*,*) "l",l,E_coul

   	 ! Calculation of the cinetic component

      DO i = 1 , NDVR-1 ! Because of the boundary at x = -xmax
       u(i) = 0.0d0*xi;
       do k = 1 , NDVR-1
         u(i) = u(i) +  T_mat(i,k)*psi(k)
         !write(*,*) "k","u",u(i),T_mat(i,k),psi(l,k)
       ENDDO
       !write(*,*) "index",n,i,i,Nr,u(i)
      ENDDO
      index_set = NDVR-2


      DO n = 2 , NSeq-1

        do k = 2 , NDVR
          u(index_set+1) = u(index_set+1) +  T_mat(index_set+1,index_set+k)*psi(index_set+k)
          !write(*,*) "k","u",index_set+1,k,u(index_set+1),T_mat(index_set+1,index_set+k),psi(l,index_set+k)
        ENDDO
        !write(*,*) "index",n,i,index_set+1,Nr,u(index_set+1)

          DO i = 2 , NDVR
            u(index_set+i) = 0.0d0*xi;
            do k = 1 , NDVR
              u(index_set+i) = u(index_set+i) +  T_mat(index_set+i,index_set+k)*psi(index_set+k)
            ENDDO
            !write(*,*) "index",n,i,index_set+i,Nr,u(index_set+i)
          ENDDO

          index_set = index_set + NDVR-1;

     ENDDO

     if(Nseq .eq.  (Nseq)*(NDVR-1)-1) then ! Boundary et x = xmax

       do k = 2 , NDVR-1
         u(index_set+1) = u(index_set+1) +  T_mat(index_set+1,index_set+k)*psi(index_set+k)
         !write(*,*) "k","u",index_set+1,k,u(index_set+1),T_mat(index_set+1,index_set+k),psi(l,index_set+k)
       ENDDO
       !write(*,*) "index",n,i,index_set+1,Nr,u(index_set+1)

       DO i = 2 , NDVR-1
         u(index_set+i) = 0.0d0*xi;
         do k = 1 , NDVR
           u(index_set+i) = u(index_set+i) +  T_mat(index_set+i,index_set+k)*psi(index_set+k)
         ENDDO
         !write(*,*) "index",n,i,index_set+i,Nr,u(index_set+i)
       ENDDO

     else

       do k = 2 , NDVR
         u(index_set+1) = u(index_set+1) +  T_mat(index_set+1,index_set+k)*psi(index_set+k)
         !write(*,*) "k","u",index_set+1,k,u(index_set+1),T_mat(index_set+1,index_set+k),psi(l,index_set+k)
       ENDDO
       !write(*,*) "index",n,i,index_set+1,Nr,u(index_set+1)

       DO i = 2 , NDVR
         u(index_set+i) = 0.0d0*xi;
         do k = 1 , NDVR
           u(index_set+i) = u(index_set+i) +  T_mat(index_set+i,index_set+k)*psi(index_set+k)
         ENDDO
         !write(*,*) "index",n,i,index_set+i,Nr,u(index_set+i)
       ENDDO

     endif

         do i = 1,Nx
	         E_cin = E_cin + u(i)*dconjg(psi(i));
         enddo



   E_calculation = E_cin+E_coul;

   write(*,*) "energies decomposition, Viriel th: ",E_cin,E_coul,E_calculation

   deallocate(u);

END FUNCTION E_calculation





complex*16 FUNCTION E_initialise(psi,x,Nx,E_guess,CV)
   complex*16,dimension(:),Allocatable, intent(in) :: x;
   complex*16,dimension(:),Allocatable :: psi;
   integer, intent(in) :: Nx;
   real*8, intent(in) :: E_guess,CV;

   !real*8 :: h,E,E_old,test,norme,delta2_1,u_old;
   !integer :: i;

   complex*16,dimension(:,:),Allocatable ::H0_mat_diag,eig_vect,VL,H0_mat_format
   complex*16, dimension(:),Allocatable :: WORK,eig;
   real*8, dimension(:),Allocatable :: RWORK;

   CHARACTER :: JOBZ,UPLO,JOBVR,JOBVL;
   INTEGER :: i,j,k,LDA,LDB,LWORK,INFO,N_init,LDZ,LDVL,LDVR,Nlanczos_INIT,N_CV,KD,LDAB;

   REAL*8 :: EPS,test_CV,coeff_lanc,beta,coef,error;

   complex*16, dimension(:,:),Allocatable :: Z,Lanczos_basis;
   complex*16, dimension(:),Allocatable :: diag_Hess,off_diag_Hess,psi_inter,psi_inter_1,psi_inter_2,psi_1(:);

   CHARACTER*8 :: method;


   method = 'diago'; !lanczos !diago

   !! Init
   !! H0_mat has been defined in init_matrix()

  if(method .EQ. 'diago') then


    write(*,*) "Diagonalization for the GS"
    N_init = 1;
    LDA = Nx;
    JOBVR = 'V'
    JOBVL = 'N'; !'N'
    UPLO = 'L';
    INFO = 0 ;
    LWORK = 2*Nx;
    KD = NDVR-1;
    LDVL = Nx;
    LDVR = Nx;


    Allocate(H0_mat_format(1:Nx,1:Nx),eig(1:Nx),WORK(1:2*Nx),RWORK(1:2*Nx),eig_vect(1:Nx,1:Nx),VL(1:Nx,1:Nx));


    DO i = 1 , Nx
      !DO j = i , Nx !min(Nx,i+KD)
      DO j = 1 , Nx !min(Nx,i+KD)

         !H0_mat_format(1+j-i,i) = H0_mat(j,i)
         H0_mat_format(i,j) = H0_mat(i,j)

      ENDDO
   ENDDO

   !!!!!!!###################################
   !! Diagonalisation method
   CALL ZGEEV( JOBVL, JOBVR, Nx, H0_mat_format, LDA, eig, VL, LDVR,  eig_vect, LDVR, WORK, LWORK, RWORK, INFO )

   error = 1.0;
   Do i = 1 , Nx

     if(abs(eig(i) - (E_guess -xi*0.2d0)).lt.error) then
       N_init = i;
       error = abs(eig(i) - (E_guess -xi*0.2d0))
     endif

   ENDDO



   psi(:) = eig_vect(:,N_init);
   E_initialise = eig(N_init);
   !!!!!!!###################################

   open(UNIT=1,FILE="eigen_energies.dat",FORM="FORMATTED",action='write');
   DO i = 1 , Nx

     write(UNIT=1,fmt='(XE13.6,XE13.6)') real(eig(i)),aimag(eig(i));

     write(*,*) i,x(i),weight(i),theta_ecs

   ENDDO

   close(1)





   !deallocate(H0_mat_diag,eig,WORK,RWORK);
   deallocate(H0_mat_format,eig,eig_vect,WORK,RWORK);

 else if(method .EQ. 'lanczos') then



!   write(*,*) "Lanczos algorithm for the GS calculation"
!
!   Nlanczos_INIT = 15
!   EPS = 1.0d0-12
!
!   JOBZ = 'V'; !! For eigenvector 'V' and eigenvalue only 'N'
!   JOBVR = 'V';
!   JOBVL = 'N';
!   LDZ = Nlanczos_INIT; LDA = Nlanczos_INIT;
!   LDVL = Nlanczos_INIT; LDVR = Nlanczos_INIT;
!
!
!
!
!   allocate(diag_Hess(1:Nlanczos_INIT),off_diag_Hess(1:Nlanczos_INIT-1),WORK(1:2*Nlanczos_INIT),Z(1:Nlanczos_INIT,1:Nlanczos_INIT));
!
!   allocate(Lanczos_basis(1:Nx,1:Nlanczos_INIT),psi_inter(1:Nx),psi_inter_1(1:Nx),psi_inter_2(1:Nx),psi_1(1:Nx));
!
!
!   do i = 1 , Nx
!      psi_1(i) = exp(-x(i)*x(i));
!   enddo
!
!  N_CV = 0; test_CV = 1.0d0;
!  DO WHILE (test_CV .gt. (1.0e-16))
!
!  	Lanczos_basis(1:Nx,1) = psi_1(1:Nx)/sqrt(scalar_product_vect(psi_1,psi_1,x,Nx));
!    psi_inter_1(1:Nx) = Lanczos_basis(1:Nx,1);
!
!
!           call apply_H0_vect(psi_inter,psi_inter_1,x,Nx);
!           diag_Hess(1) = scalar_product_vect(psi_inter_1,psi_inter,x,Nx);
!   	       psi_inter_2(1:Nx) = psi_inter(1:Nx) - diag_Hess(1)*Lanczos_basis(1:Nx,1);
!           off_diag_Hess(1) = sqrt(scalar_product_vect(psi_inter_2,psi_inter_2,x,Nx));
!
!           if(off_diag_Hess(1).le.eps) then
!   		       Lanczos_basis(1:Nx,2) = 0.0d0;
!              write(*,*) "00" ,1,off_diag_Hess(1);
!           else
!             Lanczos_basis(1:Nx,2) = psi_inter_2(1:Nx)/off_diag_Hess(1);
!                     !write(*,*) "01" ,k,test_norm;
!   	     endif
!
!   	     psi_inter_1(1:Nx) = Lanczos_basis(1:Nx,2);
!
!           do k = 2 , Nlanczos_INIT-1
!
!             call apply_H0_vect(psi_inter,psi_inter_1,x,Nx);
!             diag_Hess(k) = scalar_product_vect(psi_inter_1,psi_inter,x,Nx);
!
!   		     psi_inter_2(1:Nx) = psi_inter(1:Nx) - diag_Hess(k)*Lanczos_basis(1:Nx,k)&
!   		     & - off_diag_Hess(k-1)*Lanczos_basis(1:Nx,k-1);
!           off_diag_Hess(k) = sqrt(scalar_product_vect(psi_inter_2,psi_inter_2,x,Nx));
!
!   !write(*,*) "Norme pour off_diag_hess" ,k,off_diag_Hess(k);
!
!                   if(off_diag_Hess(k).le.eps) then
!   		              Lanczos_basis(1:Nx,k+1) = 0.0d0;
!                     write(*,*) "00" ,k,off_diag_Hess(k);
!                   else
!                     Lanczos_basis(1:Nx,k+1) = psi_inter_2(1:Nx)/off_diag_Hess(k);
!                     !write(*,*) "01" ,k,test_norm;
!   		            endif
!
!                   psi_inter_1(1:Nx) = Lanczos_basis(1:Nx,k+1);
!
!           enddo
!
!
!           call apply_H0_vect(psi_inter,psi_inter_1,x,Nx);
!           diag_Hess(Nlanczos_INIT) = scalar_product_vect(psi_inter_1,psi_inter,x,Nx);
!
!
!           CALL DSTEV( JOBZ, Nlanczos_INIT, diag_Hess, off_diag_Hess, Z, LDZ, WORK, INFO);
!
!           psi(1:Nx) = 0.0d0*xi;
!           do j = 1 , Nlanczos_INIT
!
!           	coeff_lanc = 0.0d0*xi;
!           	do k = 1 , Nlanczos_INIT
!              		coeff_lanc = coeff_lanc + Z(1,k)*Z(j,k)/(E_guess - diag_Hess(k));
!           	enddo
!
!                   psi(1:Nx) = psi(1:Nx) + coeff_lanc*Lanczos_basis(1:Nx,j);
!
!           enddo
!
!
!
!           psi(1:Nx) = psi(1:Nx)/sqrt(scalar_product_vect(psi,psi,x,Nx));
!
!           test_CV = 0.0d0
!           do j = 1 , Nx
!             if(abs(abs(psi(j)) - abs(psi_1(j))) .gt. test_CV) then
!                 test_CV = abs(abs(psi(j)) - abs(psi_1(j)));
!             endif
!
!           ENDDO
!
!           psi_1(:) = psi(:); N_CV = N_CV + 1;
!
!           !E_initialise = E_calculation_vector(psi,r,Nr,l_init)
!
!         !write(*,*) "i : ",i
!         !write(*,*) diag_Hess
!         !write(*,*) "Energy in fucntion : ",E_initialise
!         !write(*,*) "wave fucntion 1 : ",psi(13)
!         !write(*,*) "wave fucntion 1500 : ",psi(1500)
!         !write(*,*) "wave fucntion 3000 : ",psi(3000)
!         !write(*,*) "CV : ",test_CV, N_CV
!         !write(*,*)
!
!   ENDDO
!
!   DO i = 1 , Nx
!
!      if(abs(psi(i)) .LT. (1.0e-14)) then
!        psi(i) = 0.0d0;
!      ENDIF
!
!   ENDDO

!   E_initialise = E_calculation_vector(psi,x,Nx)
!
!   deallocate(diag_Hess,off_diag_Hess,WORK,Z,Lanczos_basis,psi_inter,psi_inter_1,psi_inter_2,psi_1);

   ENDIF



   open(UNIT=1,FILE="psi_gs.dat",FORM="FORMATTED",action='write');

       do i = 1,Nx
         !write(UNIT=1,fmt='(XE13.6,XE13.6,XE13.6)') real(x_rb(i)),REAL(psi(i)),AIMAG(psi(i));
         write(UNIT=1,fmt='(XE13.6,XE13.6,XE13.6)') real(x_rb(i)),REAL(psi(i)/sqrt(weight(i))),AIMAG(psi(i)/sqrt(weight(i)));
         !write(*,*) i,REAL(psi(i)/sqrt(weight(i))),AIMAG(psi(i)/sqrt(weight(i)))
       enddo

   close(1);




END FUNCTION E_initialise




!! Calculates the scalar product between two kets
FUNCTION scalar_product_vect(psi_1,psi_2,x,Nx) result(somme)

  IMPLICIT NONE;

  integer, intent(in) :: Nx;
  real*8,dimension(:),Allocatable, intent(in) :: psi_1,psi_2;
  real*8,dimension(:),Allocatable, intent(in) :: x;
  real*8 :: somme;
  integer :: i;

  somme = 0.0d0;
  do i = 1 , Nx
    somme = somme + psi_1(i)*psi_2(i);
  enddo


END FUNCTION scalar_product_vect

!%%%%%%%%%%% Calculate the energy associated to H of the wave-function
!%%%%%%%%%%% psi assumed to be projected on one spherical harmonic
complex*16 FUNCTION E_calculation_vector(psi,x,Nx)

  complex*16,dimension(:),Allocatable, intent(in) :: x;
  integer, intent(in) :: Nx;
  complex*16,dimension(:), Allocatable, intent(in) :: psi;

  complex*16 :: E_coul,E_cin,E_cent,coeff;
  integer :: i,l,k,n,index_set;
  complex*16,dimension(:),allocatable :: u;

  Allocate(u(1:Nx));


  E_coul = 0.0d0; E_cin = 0.0d0;


       do i = 1,Nx
          ! Calculation of the coulomb component
         E_coul = E_coul + psi(i)*psi(i)*pot(x(i));
       enddo

       !write(*,*) "l",l,E_coul

    ! Calculation of the cinetic component

     DO i = 1 , NDVR-1 ! Because of the boundary at x = -xmax
      u(i) = 0.0d0*xi;
      do k = 1 , NDVR-1
        u(i) = u(i) +  T_mat(i,k)*psi(k)
        !write(*,*) "k","u",u(i),T_mat(i,k),psi(l,k)
      ENDDO
      !write(*,*) "index",n,i,i,Nr,u(i)
     ENDDO
     index_set = NDVR-2


     DO n = 2 , NSeq-1

       do k = 2 , NDVR
         u(index_set+1) = u(index_set+1) +  T_mat(index_set+1,index_set+k)*psi(index_set+k)
         !write(*,*) "k","u",index_set+1,k,u(index_set+1),T_mat(index_set+1,index_set+k),psi(l,index_set+k)
       ENDDO
       !write(*,*) "index",n,i,index_set+1,Nr,u(index_set+1)

         DO i = 2 , NDVR
           u(index_set+i) = 0.0d0*xi;
           do k = 1 , NDVR
             u(index_set+i) = u(index_set+i) +  T_mat(index_set+i,index_set+k)*psi(index_set+k)
           ENDDO
           !write(*,*) "index",n,i,index_set+i,Nr,u(index_set+i)
         ENDDO

         index_set = index_set + NDVR-1;

    ENDDO

    if(Nx .eq.  (Nseq)*(NDVR-1)-1) then ! Boundary et x = xmax

      do k = 2 , NDVR-1
        u(index_set+1) = u(index_set+1) +  T_mat(index_set+1,index_set+k)*psi(index_set+k)
        !write(*,*) "k","u",index_set+1,k,u(index_set+1),T_mat(index_set+1,index_set+k),psi(l,index_set+k)
      ENDDO
      !write(*,*) "index",n,i,index_set+1,Nr,u(index_set+1)

      DO i = 2 , NDVR-1
        u(index_set+i) = 0.0d0*xi;
        do k = 1 , NDVR
          u(index_set+i) = u(index_set+i) +  T_mat(index_set+i,index_set+k)*psi(index_set+k)
        ENDDO
        !write(*,*) "index",n,i,index_set+i,Nr,u(index_set+i)
      ENDDO

    else

      do k = 2 , NDVR
        u(index_set+1) = u(index_set+1) +  T_mat(index_set+1,index_set+k)*psi(index_set+k)
        !write(*,*) "k","u",index_set+1,k,u(index_set+1),T_mat(index_set+1,index_set+k),psi(l,index_set+k)
      ENDDO
      !write(*,*) "index",n,i,index_set+1,Nr,u(index_set+1)

      DO i = 2 , NDVR
        u(index_set+i) = 0.0d0*xi;
        do k = 1 , NDVR
          u(index_set+i) = u(index_set+i) +  T_mat(index_set+i,index_set+k)*psi(index_set+k)
        ENDDO
        !write(*,*) "index",n,i,index_set+i,Nr,u(index_set+i)
      ENDDO

    endif

        do i = 1,Nx
          E_cin = E_cin + real(u(i)*psi(i),8);
        enddo

  E_calculation_vector = E_cin+E_coul;

  write(*,*) "energies decomposition, Viriel th: ",E_cin,E_coul,E_calculation_vector

  deallocate(u);

END FUNCTION E_calculation_vector






subroutine apply_H0_vect(psi_out,psi_in,x,Nx)

  IMPLICIT NONE;

  integer, intent(in) :: Nx;
  complex*16,dimension(:),Allocatable, intent(in) :: psi_in;
  complex*16,dimension(:),Allocatable :: psi_out,u;
  complex*16,dimension(:), Allocatable, intent(in) :: x;

  integer :: i,l,k,j,n,index_set;

  Allocate(u(1:Nx));



          do i = 1,Nx
             ! Calculation of the coulomb component
  	        psi_out(i) =  psi_in(i)*pot(x(i));
          enddo

          !write(*,*) "l",l,E_coul

     	 ! Calculation of the cinetic component

        DO i = 1 , NDVR-1 ! Because of the boundary at x = -xmax
         u(i) = 0.0d0;
         do k = 1 , NDVR-1
           u(i) = u(i) +  T_mat(i,k)*psi_in(k)
           !write(*,*) "k","u",u(i),T_mat(i,k),psi(l,k)
         ENDDO

        ENDDO
        index_set = NDVR-2

        DO n = 2 , NSeq-1

            do k = 2 , NDVR
              u(index_set+1) = u(index_set+1) +  T_mat(index_set+1,index_set+k)*psi_in(index_set+k)
              !write(*,*) "k","u",index_set+1,k,u(index_set+1),T_mat(index_set+1,index_set+k),psi(l,index_set+k)
            ENDDO
            !write(*,*) "index",n,i,index_set+1,Nr,u(index_set+1)

            DO i = 2 , NDVR
              u(index_set+i) = 0.0d0*xi;
              do k = 1 , NDVR
                u(index_set+i) = u(index_set+i) +  T_mat(index_set+i,index_set+k)*psi_in(index_set+k)
              ENDDO
              !write(*,*) "index",n,i,index_set+i,Nr,u(index_set+i)
            ENDDO

            index_set = index_set + NDVR-1;

       ENDDO

       if(Nseq .eq.  (Nseq)*(NDVR-1)-1) then ! Boundary et x = xmax

         do k = 2 , NDVR-1
           u(index_set+1) = u(index_set+1) +  T_mat(index_set+1,index_set+k)*psi_in(index_set+k)
           !write(*,*) "k","u",index_set+1,k,u(index_set+1),T_mat(index_set+1,index_set+k),psi(l,index_set+k)
         ENDDO
         !write(*,*) "index",n,i,index_set+1,Nr,u(index_set+1)

         DO i = 2 , NDVR-1
           u(index_set+i) = 0.0d0*xi;
           do k = 1 , NDVR
             u(index_set+i) = u(index_set+i) +  T_mat(index_set+i,index_set+k)*psi_in(index_set+k)
           ENDDO
           !write(*,*) "index",n,i,index_set+i,Nr,u(index_set+i)
         ENDDO

       else

         do k = 2 , NDVR
           u(index_set+1) = u(index_set+1) +  T_mat(index_set+1,index_set+k)*psi_in(index_set+k)
           !write(*,*) "k","u",index_set+1,k,u(index_set+1),T_mat(index_set+1,index_set+k),psi(l,index_set+k)
         ENDDO
         !write(*,*) "index",n,i,index_set+1,Nr,u(index_set+1)

         DO i = 2 , NDVR
           u(index_set+i) = 0.0d0*xi;
           do k = 1 , NDVR
             u(index_set+i) = u(index_set+i) +  T_mat(index_set+i,index_set+k)*psi_in(index_set+k)
           ENDDO
           !write(*,*) "index",n,i,index_set+i,Nr,u(index_set+i)
         ENDDO

       endif

      psi_out(:) = psi_out(:) + u(:);



     deallocate(u);

END SUBROUTINE apply_H0_vect




INTEGER FUNCTION fact(h)
   integer,intent(in) :: h;
   integer m;

   fact = 1;
   do m = 1,h
     fact = fact*m;
   enddo

end FUNCTION fact



REAL*8 FUNCTION bessin(n)
  INTEGER :: n,j,IACC,m
  REAL*8  :: bi,bim,bip,tox,bessi0,BIGNO,BIGNI,bessi,bess0;

  if(n .eq. 0) then
     bessin = 1.2660658777520083356d0
  endif
  if(n .eq. 1) then
     bessin = 0.56515910399248502721d0
  endif
  if(n .ge. 2) then

    IACC=40
    BIGNO=1.0e10
    BIGNI=1.0e-10
    bess0 = 1.2660658777520083356d0;

    tox=2.0d0; bip=0.0d0; bi=1.0d0
    bessi=0.d0
    m = 2*((n+int(sqrt(float(IACC*n))))) !Downward  recurrence from evenm.
    do j = m, 1, -1 !MakeIACClarger to increase accuracy.
      bim = bip+real(j,8)*tox*bi !The  downward recurrence.
      bip = bi
      bi = bim
      if (abs(bi) .gt. BIGNO) then
        bessi = bessi*BIGNI
        bi = bi*BIGNI
        bip = bip*BIGNI
      endif
      if (j .eq. n) then
         bessi=bip
      endif
    enddo
    bessin = bess0*bessi/bi

  endif




END FUNCTION bessin




end module tools
