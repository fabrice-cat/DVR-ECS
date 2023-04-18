module prop

use util; use tools; use tridiag; use potential;

IMPLICIT NONE;

CONTAINS



!%%%%%%%%%%%%%%  Definition of the electric field
real*8 FUNCTION E_field(F0,omega,t,phi_CEP,nc,t_on,t_off,param_env)
  real*8,intent(in) :: F0,omega,t,phi_CEP,t_on,t_off;
  integer,intent(in) :: Nc,param_env;

  real*8 :: h;

  h = 1.0d-10;


  E_field = A_field(F0,omega,t+h,phi_CEP,nc,t_on,t_off,param_env)-A_field(F0,omega,t,phi_CEP,nc,t_on,t_off,param_env);
  E_field = -E_field/h;



end function E_field


!%%%%%%%%%%%%%%  Definition of the potential vector
real*8 FUNCTION A_field(F0,omega,t,phi_CEP,nc,t_on,t_off,param_env)
  real*8,intent(in) :: F0,omega,t,phi_CEP,t_on,t_off;
  integer,intent(in) :: Nc,param_env;

  real*8 :: omegap,ts,tf,amp,coef,period;

  period = 2.0d0*pi/omega;

  omegap = omega/dfloat(Nc);

   SELECT CASE (param_env)
      CASE (0)
         A_field = (F0/omega)*dsin(omega*t);

         if(t.gt.(dfloat(nc)*period).or.(t.lt.0)) then
            A_field = 0.0d0;
         endif

      CASE (1)
         if (t.lt.(t_on*2.0d0*pi/omega)) then
             A_field = (F0/omega)*(sin(0.5d0*pi*t/(t_on*2.0d0*pi/omega)))**(4.0);
         else
              if (((float(Nc)-t_off)*2.0d0*pi/omega).lt.t) then
                  amp = (F0/omega)*dcos((t_off-t_on)*2.0d0*pi);
                  ts = (float(Nc)-t_off)*2.0d0*pi/omega;
                  tf = nc*2.0d0*pi/omega;
                  A_field = amp*(dsin(0.5d0*pi*(t-tf)/(ts-tf))**(4.0));
              else
                  A_field = (F0/omega)*dcos(omega*(t-t_on*2.0d0*pi/omega));
              endif
         endif
         if(t.gt.(dfloat(nc)*period).or.(t.lt.0)) then
            A_field = 0.0d0;
         endif

      CASE (2)
         A_field = (F0/omega)*dsin(0.50d0*omegap*t)*dsin(0.50d0*omegap*t)*dsin(omega*t+phi_CEP);

         if(t.gt.(dfloat(nc)*period).or.(t.le.0)) then
            A_field = 0.0d0;
         endif

   END SELECT


end function A_field


!%%%%%%%%%%%%%%  Definition of the potential vector
real*8 FUNCTION A_field_XUV(F0,omega,t,phi_CEP,nc,t_on,t_off,param_env)
  real*8,intent(in) :: omega,t,phi_CEP,t_on,t_off;
  integer,intent(in) :: Nc,param_env;

  real*8 :: F0,omegap,ts,tf,amp,coef,period;

  period = 2.0d0*pi/omega;

  omegap = omega/dfloat(Nc);

  F0 = 1.0d-5;

  A_field_XUV = F0*(cos(11.0d0*omega*t)+cos(13.0d0*omega*t)+cos(15.0d0*omega*t)+cos(17.0d0*omega*t)+cos(19.0d0*omega*t))/omega;
  A_field_XUV = A_field_XUV*dsin(0.50d0*omegap*t)*dsin(0.50d0*omegap*t)


  if(t.gt.(dfloat(nc)*period).or.(t.lt.0)) then
       A_field_XUV = 0.0d0;
  endif

end function A_field_XUV




subroutine propagation(F0,omega,phi_CEP,Nc,Nt,t_on,t_off,t_on2,t_off2,param_env,param_env2,t_delay,psi,Nx,Nexp,x,Nlanczos)
  real*8,intent(in) :: F0,omega,phi_CEP,t_on,t_off,t_on2,t_off2;
  Complex*16,dimension(:), Allocatable :: psi;
  Complex*16,dimension(:,:), Allocatable :: Lanczos_basis;
  Complex*16,dimension(:), Allocatable :: x;
  Integer,intent(in) :: Nc,param_env,param_env2,Nt,Nexp,Nlanczos;
  real*8,intent(in) :: t_delay;
  Integer :: Nx,test,N_CV;
  CHARACTER :: JOBZ,JOBVL,JOBVR;


  Integer :: i,j,k,l,m,number_cycle,ref,ind,ind_abs,l_finish,INFO,LDZ,coeff_taylor,Nb_taylor,index_prop,N_int_field,k_min,Dk,Nstop;
  real*8 :: t_start,t_stop,coeff,lambda,alpha,beta,eps,dt_test;
  real*8 :: A_champ,c1,A_champ1,c2,A_champ2,c3,A_champ3,E_champ,t,t1,t2,t3,tc,dipole,test_norm,intermediaire,test_CV,theta_CAP,dt_ref;
  Complex*16,dimension(:), Allocatable :: mem_psi,psi_inter_1,psi_inter_2,psi_inter,psi_test;
  real*8,dimension(:,:), Allocatable :: Z;
  real*8,dimension(:), Allocatable :: diag_Hess,off_diag_Hess,WORK,vect_inter;
  Complex*16 :: absorbtion,proj,init,zptf,coeff_lanc,ref_psi,ortho1,ortho2,ortho3,coeff_lanc1,coeff_lanc2,VCAP;

  INTEGER :: LDA,LDVL,LDVR,LWORK;
  real*8,dimension(:,:), Allocatable :: VL,VR,Hessenberg;
  real*8,dimension(:), Allocatable :: WR,WI;

  COMPLEX*16 :: somme_inter;

  open(UNIT=1,FILE="dipole.dat",FORM="FORMATTED",action='write');
  open(UNIT=2,FILE="test_cv.dat",FORM="FORMATTED",action='write');


  write(*,*)
  write(*,*) "Propagation procedure ... "
  write(*,*)

  allocate(psi_test(1:Nx),psi_inter_1(1:Nx),psi_inter_2(1:Nx),psi_inter(1:Nx));
  allocate(vect_inter(1:Nx));


  !! Init
  JOBZ = 'V'; !! For eigenvector 'V' and eigenvalue only 'N'
  JOBVR = 'V';
  JOBVL = 'N';
  LDZ = Nlanczos; LDA = Nlanczos;
  LDVL = Nlanczos; LDVR = Nlanczos;
  !LWORK = 4*Nlanczos;


  allocate(diag_Hess(1:Nlanczos),off_diag_Hess(1:Nlanczos-1),WORK(1:2*Nlanczos),Z(1:Nlanczos,1:Nlanczos));

  allocate(Hessenberg(1:Nlanczos,1:Nlanczos),VL(1:Nlanczos,1:Nlanczos),VR(1:Nlanczos,1:Nlanczos),WR(1:Nlanczos),WI(1:Nlanczos));


  allocate(Lanczos_basis(1:Nx,1:Nlanczos));

  t = 0.0d0; tc = 2.0d0*pi/omega; number_cycle = 1; call CPU_TIME(t_start);
  write(*,'(a)',advance='no') "Cycle number : "; write(*,'(I5)',advance='no') number_cycle;
  write(*,'(a)',advance='no') "   Size of the box : ";write(*,'(I5)') Nx;
  ref = int(dfloat(Nt)/(20.0d0*dfloat(nc))); ind = ref;


  !dt = ((2.0d0*pi/omega2)*dfloat(Nc2)+t_delay)/dfloat(Nt);

  !dt = (2.0d0*pi/omega)*dfloat(Nc)/dFloat(Nt);


  eps = 1.0d-16; ref_psi = psi(1); N_CV = 1; dt_ref = dt;


  N_int_field = 10

  Lanczos_basis(1:Nx,1:Nlanczos) = 0.0d0*xi;

!  write(*,*) "Taylor propagation"

  !!!!!! INIT WITH TAYLOR
Nb_taylor = 0;
!  coeff_taylor = 10; Nb_taylor = 10;
!  do i = 1,Nb_taylor*coeff_taylor

   !if(i.gt.Nt) then
   ! A_champ = 0.0d0; E_champ = 0.0d0;
   !else

    !A_champ = A_field(F0,omega,t,phi_CEP,Nc,t_on,t_off,param_env) + A_field(F02,omega2,t-t_delay,phi_CEP2,Nc2,t_on2,t_off2,param_env2);
    !E_champ = E_field(F0,omega,t,phi_CEP,Nc,t_on,t_off,param_env);

!     A_champ = A_field(F0,omega,t,phi_CEP,Nc,t_on,t_off,param_env) + 0.0d0*A_field_XUV(F02,omega2,t-t_delay,phi_CEP2,Nc2,t_on2,t_off2,param_env2);

   !endif

!    psi_inter_1(abs(m_init):lmax,1:Nr) = psi(abs(m_init):lmax,1:Nr);
!    do k = 6 , 1 , -1
!      call apply_H(psi_inter,psi_inter_1,A_champ,r,Lmax,Nr,m_init);
!      psi_inter_1(abs(m_init):lmax,1:Nr) = psi(abs(m_init):lmax,1:Nr) - xi*(dt/real(coeff_taylor,8))*psi_inter(abs(m_init):lmax,1:Nr)/real(k,8);
      !psi(:,:) = psi(:,:) + (-xi*dt*0.001d0)**(dfloat(k))*psi_inter(:,:)/fact(k);
      !psi_inter_1(:,:) = psi_inter(:,:);
!    enddo

!    psi(abs(m_init):lmax,1:Nr) = psi_inter_1(abs(m_init):lmax,1:Nr);

    !psi(:,:) = dconjg(psi(:,:));

!    t = t + dt/real(coeff_taylor,8);

!    write(*,*) i,Nb_taylor*coeff_taylor,A_champ!,psi(0,1);
!    write(*,*) "norm after TAYLOR :", real(scalar_product(psi,psi,Lmax,r,Nr,m_init));

!  enddo

!  ind = Nb_taylor;

  write(*,*) "Regular propagation"

  !write(*,*)

  write(*,*) A_champ;
  write(*,*) psi(1);
  write(*,*) ref_psi*exp(-xi*t*(-0.5d0));

  index_prop = 0;
!go to 1000


!  do k = 1 , Nr
!        CALL RANDOM_NUMBER(intermediaire);
!  	psi(0,k) = psi(0,k) + intermediaire*1.0d-10;
!  enddo



  !!!!!! TRUE PROPAGATION


!  psi(:,:) = 0.0d0*xi;
!
!  do j = 1 , Nr
!
!
!   psi(10,j) = (6.0d0+9.0d0*r(j)*r(j)-18.0d0*r(j)**(3.0d0))*exp(-0.3d0*r(j))*(1.0d0-exp(-10.d0*(r(j)-0.2d0*dfloat(Nr))**(2.0d0)));
!
!  enddo
!
!  psi(:,:) = psi(:,:)/dsqrt(real(scalar_product(psi,psi,Lmax,r,Nr,m_init)));
!
!  write(*,*) "test matlab"
!  write(*,*) psi(0,1)



  psi_test(:)  = psi(:);
  dt_test = dt/10.0d0;

  ! Loop on the time
  do i = Nb_taylor+1,Nt

   !write(*,*)
   !write(*,*) i,Nt

   if(i.gt.Nt) then
    A_champ = 0.0d0; E_champ = 0.0d0;
   else

    !A_champ = A_field(F0,omega,t,phi_CEP,Nc,t_on,t_off,param_env) + A_field(F02,omega2,t-t_delay,phi_CEP2,Nc2,t_on2,t_off2,param_env2);
    !E_champ = E_field(F0,omega,t,phi_CEP,Nc,t_on,t_off,param_env);

    !if(index_prop .eq. 10) then
     A_champ = A_field(F0,omega,t+dt*0.5d0,phi_CEP,Nc,t_on,t_off,param_env)!+ 0.0d0*A_field_XUV(F02,omega2,t-t_delay,phi_CEP2,Nc2,t_on2,t_off2,param_env2);
     !index_prop = 0
    !else
    !  index_prop = index_prop +1;
    !endif
     t1 = t;

     t2 = t +dt;

     t3 = t +dt*(0.5d0+sqrt(3.0d0/20.0d0));

     A_champ1 = A_field(F0,omega,t1,phi_CEP,Nc,t_on,t_off,param_env)

     A_champ2 = A_field(F0,omega,t2,phi_CEP,Nc,t_on,t_off,param_env)

     !A_champ3 = A_field(F0,omega,t3,phi_CEP,Nc,t_on,t_off,param_env)

     !c1 = 37.0d0/240.0d0+10.0d0*sqrt(5.0d0/3.0d0)/87.0d0; c2 = -1.0d0/30.0d0; c3 = 37.0d0/240.0d0-10.0d0*sqrt(5.0d0/3.0d0)/87.0d0

     E_champ = A_field(F0,omega,t+dt,phi_CEP,Nc,t_on,t_off,param_env) + 0.0d0*A_field_XUV(F02,omega2,t-t_delay+dt,phi_CEP2,Nc2,t_on2,t_off2,param_env2);

     E_champ = -(E_champ - A_champ)/dt;

   endif

    !!write(*,*) param_env,param_env2,t_delay,A_champ,A_field(F0,omega,t,phi_CEP,Nc,t_on,t_off,param_env),A_field(F02,omega2,t-t_delay,phi_CEP2,Nc2,t_on2,t_off2,param_env2);

    !absorbtion = 1.0d-15; !psi(0,165);




      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     !                            Propagation procedure


        !write(*,*) "mod de psi", psi(0,1)*dconjg(psi(0,1))*h;
        !write(*,*) "norm", scalar_product(psi,psi,Lmax,r,Nr,m_init);





166  	Hessenberg(1:Nlanczos,1:Nlanczos) = 0.0d0;

!      psi_inter(:,:) = psi(:,:);
!      psi(:,:) = (0.0d0,0.0d0);
!      Nstop = 150;
!      do l = abs(m_init) , Lmax
!
!          do j = 1 , Nr
!
!
!             k_min = int(max(1,j-Nstop))
!             Dk = int(min(Nr,j+Nstop))
!
!              somme_inter = (0.0d0,0.0d0)
!              DO k = 1,Nr !k_min , Dk
!                somme_inter = somme_inter + expm_mat(l,j,k)*psi_inter(l,k)
!              ENDDO
!              psi(l,j) = somme_inter;
!
!            !psi(l,j) = exp(-xi*dt*0.25*real(l*(l+1),8)/(r(j)*r(j)))*psi(l,j)
!
!          enddo
!      enddo

      !psi(:,:) = psi(:,:)/sqrt(real(scalar_product(psi,psi,Lmax,r,Nr,m_init)));

        !alpha = -0.5d0
        !!! Create the Lanczos basis

    !write(*,*) "test wave function"

    !somme_inter = (0.0d0,0.0d0)
    !DO k = 1 , Nr
    !  somme_inter = somme_inter + psi_inter(5,k)*psi_inter(5,k)
    !ENDDO

    !write(*,*) somme_inter

    !write(*,*) psi_inter(5,1)

    !write(*,*) psi(5,1000)


    !somme_inter = (0.0d0,0.0d0)
    !DO k = 1 , Nr
    !  somme_inter = somme_inter + psi(0,k)*dconjg(psi(0,k))
    !ENDDO

    !write(*,*) somme_inter

    !PAUSE


        666        DO m = 1 , N_CV


        psi_inter(:) = psi(:);
        !psi(:,:) = (0.0d0,0.0d0);
        Nstop = 150;


                A_champ = A_field(F0,omega,t+dt*0.5d0 + real(m-1,8)*dt,phi_CEP,Nc,t_on,t_off,param_env)
        	      Lanczos_basis(1:Nx,1) = psi_inter(1:Nx)/sqrt(real(scalar_product(psi_inter,psi_inter,x,Nx)));
                psi_inter_1(1:Nx) = Lanczos_basis(1:Nx,1);

                !call apply_OMEGA(psi_inter,psi_inter_1,c1,A_champ1,c2,A_champ2,c3,A_champ3,r,Lmax,Nr,m_init,dt);
                call apply_H(psi_inter,psi_inter_1,A_champ,x,Nx);
                !call apply_HCF(psi_inter,psi_inter_1,A_champ,r,Lmax,Nr,m_init,l_cut);
                !call apply_Apz(psi_inter,psi_inter_1,A_champ,r,Lmax,Nr,m_init,l_cut);
                !psi_inter(:,:) = psi_inter(:,:) - alpha*psi_inter_1(:,:);
                diag_Hess(1) = real(scalar_product(psi_inter_1,psi_inter,x,Nx));
        	      psi_inter_2(1:Nx) = psi_inter(1:Nx) - diag_Hess(1)*Lanczos_basis(1:Nx,1);
                off_diag_Hess(1) = sqrt(real(scalar_product(psi_inter_2,psi_inter_2,x,Nx)));
                !Lanczos_basis(abs(m_init):Lmax,1:Nr,2) = psi_inter_2(abs(m_init):Lmax,1:Nr)/off_diag_Hess(1);

                if(off_diag_Hess(1).le.eps) then
        		       Lanczos_basis(1:Nx,2) = 0.0d0*xi;
                   write(*,*) "00" ,1,off_diag_Hess(1);
                else
                  Lanczos_basis(1:Nx,2) = psi_inter_2(1:Nx)/off_diag_Hess(1);
                          !write(*,*) "01" ,k,test_norm;
        	     endif

        	     psi_inter_1(1:Nx) = Lanczos_basis(1:Nx,2);

                do k = 2 , Nlanczos-1
                	!call apply_OMEGA(psi_inter,psi_inter_1,c1,A_champ1,c2,A_champ2,c3,A_champ3,r,Lmax,Nr,m_init,dt);
                  call apply_H(psi_inter,psi_inter_1,A_champ,x,Nx);
                  !call apply_HCF(psi_inter,psi_inter_1,A_champ,x,Nx);
                  !call apply_Apz(psi_inter,psi_inter_1,A_champ,r,Lmax,Nr,m_init,l_cut);
                  !psi_inter(:,:) = psi_inter(:,:) - alpha*psi_inter_1(:,:);
        !psi_inter_2(abs(m_init):Lmax,1:Nr) = psi_inter(abs(m_init):Lmax,1:Nr) - off_diag_Hess(k-1)*Lanczos_basis(abs(m_init):Lmax,1:Nr,k-1)
                	diag_Hess(k) = real(scalar_product(psi_inter_1,psi_inter,x,Nx));

        		     psi_inter_2(1:Nx) = psi_inter(1:Nx) - diag_Hess(k)*Lanczos_basis(1:Nx,k)&
        		     & - off_diag_Hess(k-1)*Lanczos_basis(1:Nx,k-1);
                off_diag_Hess(k) = sqrt(real(scalar_product(psi_inter_2,psi_inter_2,x,Nx)));

        !write(*,*) "Norme pour off_diag_hess" ,k,off_diag_Hess(k);

                        if(off_diag_Hess(k).le.eps) then
        		              Lanczos_basis(1:Nx,k+1) = 0.0d0*xi;
                          write(*,*) "00" ,k,off_diag_Hess(k);
                        else
                          Lanczos_basis(1:Nx,k+1) = psi_inter_2(1:Nx)/off_diag_Hess(k);
                          !write(*,*) "01" ,k,test_norm;
        		            endif


                        psi_inter_1(1:Nx) = Lanczos_basis(1:Nx,k+1);

                enddo

        	!psi_inter_1(abs(m_init):Lmax,1:Nr) = Lanczos_basis(abs(m_init):Lmax,1:Nr,Nlanczos)

                !call apply_OMEGA(psi_inter,psi_inter_1,c1,A_champ1,c2,A_champ2,c3,A_champ3,r,Lmax,Nr,m_init,dt);
                call apply_H(psi_inter,psi_inter_1,A_champ,x,Nx);
                !call apply_HCF(psi_inter,psi_inter_1,A_champ,r,Lmax,Nr,m_init,l_cut);
                !call apply_Apz(psi_inter,psi_inter_1,A_champ,r,Lmax,Nr,m_init,l_cut);
                !psi_inter(:,:) = psi_inter(:,:) - alpha*psi_inter_1(:,:);
                !psi_inter_2(abs(m_init):Lmax,1:Nr) = psi_inter(abs(m_init):Lmax,1:Nr) - off_diag_Hess(Nlanczos-1)*Lanczos_basis(abs(m_init):Lmax,1:Nr,Nlanczos-1)
                diag_Hess(Nlanczos) = real(scalar_product(psi_inter_1,psi_inter,x,Nx));


              !psi_inter_1(abs(m_init):Lmax,1:Nr) = Lanczos_basis(abs(m_init):Lmax,1:Nr,1);
              !psi_inter(abs(m_init):Lmax,1:Nr) = Lanczos_basis(abs(m_init):Lmax,1:Nr,2);

              ! write(*,*) "ortho", real(scalar_product(psi_inter,psi_inter_1,Lmax,r,Nr,m_init));


               test_CV = 1.0d0;
               do k = 1 , Nlanczos-1
        	       test_CV = test_CV*off_diag_Hess(k)*dt/dfloat(k);
               enddo


               if(test_CV .gt. (1.0d-20)) then
        6661      N_CV = N_CV*2; dt = dt/2.0d0;

                  test_CV = 1.0d0;
                  do k = 1 , Nlanczos-1
                    test_CV = test_CV*off_diag_Hess(k)*dt/dfloat(k);
                  enddo
                  if( test_CV .gt. (1.0d-20) ) then
                    GOTO 6661
                  endif

                  GO TO 666;
               endif



               !write(*,*) i,diag_Hess(:)
               !write(*,*) off_diag_Hess
               !write(*,*)

                CALL DSTEV( JOBZ, Nlanczos, diag_Hess, off_diag_Hess, Z, LDZ, WORK, INFO);

        !        write(*,*) "Ritz 2",Nr,Lmax,Nlanczos
        !        do j = 1 , Nlanczos
        !		      write(*,*) diag_Hess(j),Z(2,j);
        !		      write(*,*) INFO,WR(j),Z(j,3);
        !        enddo

                !write(*,*) i,diag_Hess(:)
                !write(*,*)


           	    test_norm = 0.0d0;
                psi(1:Nx) = 0.0d0*xi;
                do j = 1 , Nlanczos


                	coeff_lanc = 0.0d0*xi;
                	do k = 1 , Nlanczos
                      !coeff_lanc = coeff_lanc + (exp(-xi*dt*diag_Hess(k)-0.6d0*dt)-1.0d0)*Z(1,k)*Z(j,k)/(xi*diag_Hess(k)+0.6d0);
                   		coeff_lanc = coeff_lanc + exp(-xi*dt*(diag_Hess(k)))*Z(1,k)*Z(j,k);
                      !coeff_lanc = coeff_lanc + (1.0d0/(1.0d0+xi*dt*(diag_Hess(k))/2.0d0))*Z(1,k)*Z(j,k);
                	enddo



                        !write(*,*) real(scalar_product(psi_inter,psi,Lmax,r,Nr,m_init));

                        !test_norm = scalar_product(psi_inter,psi,Lmax,r,Nr,m_init);

                        psi(1:Nx) = psi(1:Nx) + coeff_lanc*Lanczos_basis(1:Nx,j);


                        test_norm = test_norm + real(coeff_lanc*conjg(coeff_lanc));


                enddo



              ENDDO



        !      test_CV = error_test(psi,psi_test,Lmax,Nr,m_init)
              write(*,*) "time",i,Nt
              write(*,*) "error : ",test_CV,N_CV,dt



        !      write(UNIT=2,fmt='(XE13.6,XE13.6,XE13.6,XE13.6,XE13.6,XE13.6)') t+dt,real(psi(1,50)),aimag(psi(1,50)),real(psi_test(1,50)),aimag(psi_test(1,50)),test_CV;





       N_CV = 1; dt = dt_ref;

!       psi_inter(:,:) = psi(:,:);
!       psi(:,:) = (0.0d0,0.0d0);
!       do l = abs(m_init)  , Lmax
!           do j = 1 , Nr
!
!
!              k_min = int(max(1,j-Nstop))
!              Dk = int(min(Nr,j+Nstop))
!              somme_inter = (0.0d0,0.0d0)
!
!               DO k = 1,Nr !k_min ,  Dk
!                 somme_inter = somme_inter + expm_mat(l,j,k)*psi_inter(l,k)
!               ENDDO
!               psi(l,j) = somme_inter;
!
!               !psi(l,j) = exp(-xi*dt*0.25*real(l*(l+1),8)/(r(j)*r(j)))*psi(l,j)
!
!           enddo
!       enddo


       write(*,*) "REAL norm :", sqrt(real(scalar_product(psi,psi,x,Nx)));




     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ( t < tc) then

     ! index for the number of cycle
     if (ind.eq.ref) then
      call print_symbol();
      ind = 1;
     else
       ind = ind + 1;
     endif

    else

     tc = tc + 2.0d0*pi/omega;
     number_cycle = number_cycle + 1; Nx = Nx + Nexp;

     ! expand the grid and psi
!    deallocate(r); allocate(r(1:Nr));
!       do j = 1,Nr
!         r(j) = dfloat(j)*h;
!       enddo
!     allocate(mem_psi(abs(m_init):Lmax,1:Nr)); mem_psi(abs(m_init):lmax,1:Nr) = 0.0d0; mem_psi(abs(m_init):lmax,1:Nr-Nexp) = psi(abs(m_init):lmax,1:Nr-Nexp);
!     deallocate(psi); allocate(psi(abs(m_init):Lmax,1:Nr)); psi = mem_psi; deallocate(mem_psi);
!     vect_inter = psi0; deallocate(psi0); allocate(psi0(1:Nr)); psi0 = 0.0d0; psi0(1:Nr-Nexp) = vect_inter(1:Nr-Nexp);

!     deallocate(Lanczos_basis,psi_inter_1,psi_inter_2,psi_inter);
!     allocate(Lanczos_basis(abs(m_init):Lmax,1:Nr,NLanczos)); allocate(psi_inter_1(abs(m_init):Lmax,1:Nr),psi_inter_2(abs(m_init):Lmax,1:Nr),psi_inter(abs(m_init):Lmax,1:Nr))


!     deallocate(vect_inter); allocate(vect_inter(1:Nr));


     call CPU_TIME(t_stop);
     write(*,*)
     write(*,*) "Duration of the calculation for the cycle: ",t_stop-t_start,"sec";
     call CPU_TIME(t_start);

     write(*,*) "norme de psi",psi(1),test_norm;
     write(*,*) "norme de psi",scalar_product(psi,psi,x,Nx);
     write(*,*) "Convergence test1",test_CV
     proj = 0.0d0;
     do j = 1  , Nx
        proj = proj + real(psi(j)*conjg(psi(j)),8)
     enddo
     write(*,'(a)',advance='no') "Population of lmax : "; write(*,'(XE13.6)') real(proj,8);
    ! test_CV = 1.0d0;
     !do k = 1 , Nlanczos-1
      !test_CV = test_CV*off_diag_Hess(k)*dt/real(k,8);
     !enddo

     !write(*,*) "Convergence test2",test_CV

     write(*,*); write(*,*);
     write(*,'(a)',advance='no') "Cycle number : "; write(*,'(I5)',advance='no') number_cycle;
     write(*,'(a)',advance='no') "   Size of the box : ";write(*,'(I5)') Nx;
     write(*,'(a)',advance='no') "xmax : "; write(*,'(XE13.6)') x(Nx);




    endif




    ! Apply the attenuation factor at the boundary r = rmax ; absorbing wall
    !ind_abs = 20;
    !do j = Nr-ind_abs , Nr
  !    absorbtion = exp(-1.0*(r(Nr-ind_abs)-r(j))**(2.0));
  !    psi(:,j) =  absorbtion*psi(:,j);
  !  enddo
    !call normalize(psi,r,Nr,Lmax);

    ! Apply the attenuation factor at the boundary r = rmax ; absorbing wall
    !ind_abs = 20; theta_CAP = 0.2d0;
    !do j = Nr-ind_abs , Nr
      !VCAP = -25.0d0*xi*(1.0d0-cos(pi*(r(j)-r(Nr-ind_abs))/(2.0d0*(r(Nr)-r(Nr-ind_abs)))));
      !absorbtion = exp(-xi*VCAP);
    !  absorbtion = cos(pi*(r(j)-r(Nr-ind_abs))/(2.0d0*(r(Nr)-r(Nr-ind_abs))))**(0.15d0);
    !  psi(:,j) =  absorbtion*psi(:,j);
    !enddo
    !call normalize(psi,r,Nr,Lmax);

    ! Calculate the dipole  <z..>(t)
    dipole = 0.0d0;
      do j = 1 , Nx
        dipole = dipole + real(psi(j)*conjg(psi(j))*gradzpot(x(j)),8);
      enddo







    ! Project on the GS
    proj = 0.0d0*(1+xi);
    do j = 1 , Nx
      proj = proj + psi0(j)*psi(j);
    enddo

    t = t + dt; time = t;

    write(UNIT=1,fmt='(XE13.6,XE13.6,XE13.6,XE13.6,XE13.6,XE13.6)') t,A_champ,E_champ,real(proj*conjg(proj)),dipole,test_CV;
    !write(UNIT=1,fmt='(XE12.6,XE12.6)') t,A_champ;



  enddo





  call CPU_TIME(t_stop);
  write(*,*)
  write(*,*) "Duration of the calculation for the cycle: ",t_stop-t_start,"sec";




  ! Calculate the dipole  <z.>(tf) This is used for the calculation of <z.>(omega)
!  zptf = 0.0d0*xi;
!  do l = 0,Lmax-1
!   proj = 0.0d0*xi;
 !  do j = 1 , Nr-1

     ! Calculate the term <conj(phi_(l+1))phi_l/r>
 !     proj = proj + (psi(l,j)*dconjg(psi(l+1,j))*r(j)**(-1.0d0)+psi(l,j+1)*dconjg(psi(l+1,j+1))*r(j+1)**(-1.0d0))*h*0.5d0;
 !  enddo

     ! Calculate the term <der(phi_(l+1))/dr*conj(phi_l)/r>
          !! psi_inter is delta1*psi_(l+1)
!       vect_inter(1) = lambda*psi(l+1,1) + psi(l+1,2);
!       do k = 2,Nr-1
!           vect_inter(k) = psi(l+1,k+1) - psi(l+1,k-1);
!       enddo
!       vect_inter(Nr) = - psi(l+1,Nr-1);
!       vect_inter(:) = vect_inter(:)*(2.0d0*h)**(-1.0d0);

         !! This is M1
!       diag(1) = (4.0d0+lambda)/6.0d0; sup(1) = 1.0d0/6.0d0;
!       do j = 2 , Nr
!   	inf(j) = 1.0d0/6.0d0;
!        sup(j) = 1.0d0/6.0d0;
!        diag(j) = 4.0d0/6.0d0;
!       enddo
!
!       call tridiag_gauss(inf,diag,sup,psi_inter1,psi_inter2,Nr); ! psi_inter2 = derivative of phi_(l+1)
!
!     ! add the term <der(phi_(l+1))/dr*conj(phi_l)>
!    do j = 1 , Nr-1
!       proj = proj + (psi_inter2(j)*dconjg(psi(l,j))+psi_inter2(j+1)*dconjg(psi(l,j+1)))*h*0.5d0;
!    enddo


 !  zptf = zptf + 2.0d0*aimag(proj)*ang_c(l,m_init);
 ! enddo



1000     write(*,*) "norme de psi",psi(1),scalar_product(psi,psi,x,Nx);
         write(*,*) "norme de psi",test_norm;
         write(*,*) "Convergence test",test_CV


  !open(UNIT=2,FILE="zptf.dat",FORM="FORMATTED",action='write');

  !write(UNIT=2,fmt='(XE13.6,XE13.6)') real(zptf),aimag(zptf);



 close(1); close(2);

 deallocate(psi_inter_1); deallocate(psi_inter_2); deallocate(psi_inter);
 deallocate(diag_Hess); deallocate(off_diag_Hess); deallocate(Lanczos_basis);

end subroutine propagation



!!!!!!!!!!!!!!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!!!!!!!!!!!!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   END OF FUNCTION PROPAGATION
!!!!!!!!!!!!!!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




subroutine apply_H0(psi_out,psi_in,x,Nx)

  IMPLICIT NONE;

  integer, intent(in) :: Nx;
  complex*16,dimension(:),Allocatable, intent(in) :: psi_in;
  complex*16,dimension(:),Allocatable :: psi_out;
  complex*16,dimension(:),Allocatable :: u;
  complex*16,dimension(:), Allocatable, intent(in) :: x;

  integer :: i,l,k,j,n,index_set;

  Allocate(u(1:Nx));


          do i = 1,Nx
             ! Calculation of the coulomb component
  	        psi_out(i) = psi_in(i)*pot(x(i));
          enddo

          !write(*,*) "l",l,E_coul

     	 ! Calculation of the cinetic component

        DO i = 1 , NDVR-1 ! Because of the boundary at r = 0
         u(i) = 0.0d0*xi;
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

       if(Nseq .eq.  (Nseq)*(NDVR-1)-1) then ! Boundary et r = Rmax

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

END SUBROUTINE apply_H0





subroutine apply_Dr(psi_out,psi_in,x,Nx)

  IMPLICIT NONE;

  integer, intent(in) :: Nx;
  complex*16,dimension(:),Allocatable, intent(in) :: psi_in;
  complex*16,dimension(:),Allocatable :: psi_out;
  complex*16,dimension(:),Allocatable :: u;
  complex*16,dimension(:), Allocatable, intent(in) :: x;

  integer :: i,l,k,j,n,index_set;

  complex*16 :: proj;

  Allocate(u(1:Nx));



      psi_out(:) = 0.0d0*xi;


        DO i = 1 , NDVR-1 ! Because of the boundary at r = 0
         u(i) = 0.0d0*xi;
         do k = 1 , NDVR-1
           u(i) = u(i) +  D1_mat(i,k)*psi_in(k)
           !write(*,*) "k","u",u(i),T_mat(i,k),psi(l,k)
         ENDDO

        ENDDO
        index_set = NDVR-2


        DO n = 2 , NSeq-1

          do k = 2 , NDVR
            u(index_set+1) = u(index_set+1) +  D1_mat(index_set+1,index_set+k)*psi_in(index_set+k)
            !write(*,*) "k","u",index_set+1,k,u(index_set+1),T_mat(index_set+1,index_set+k),psi(l,index_set+k)
          ENDDO
          !write(*,*) "index",n,i,index_set+1,Nr,u(index_set+1)

            DO i = 2 , NDVR
              u(index_set+i) = 0.0d0*xi;
              do k = 1 , NDVR
                u(index_set+i) = u(index_set+i) +  D1_mat(index_set+i,index_set+k)*psi_in(index_set+k)
              ENDDO
              !write(*,*) "index",n,i,index_set+i,Nr,u(index_set+i)
            ENDDO

            index_set = index_set + NDVR-1;

       ENDDO

       if(Nx .eq.  (Nseq)*(NDVR-1)-1) then ! Boundary et r = Rmax

         do k = 2 , NDVR-1
           u(index_set+1) = u(index_set+1) +  D1_mat(index_set+1,index_set+k)*psi_in(index_set+k)
           !write(*,*) "k","u",index_set+1,k,u(index_set+1),T_mat(index_set+1,index_set+k),psi(l,index_set+k)
         ENDDO
         !write(*,*) "index",n,i,index_set+1,Nr,u(index_set+1)

         DO i = 2 , NDVR-1
           u(index_set+i) = 0.0d0*xi;
           do k = 1 , NDVR
             u(index_set+i) = u(index_set+i) +  D1_mat(index_set+i,index_set+k)*psi_in(index_set+k)
           ENDDO
           !write(*,*) "index",n,i,index_set+i,Nr,u(index_set+i)
         ENDDO

       else

         do k = 2 , NDVR
           u(index_set+1) = u(index_set+1) +  D1_mat(index_set+1,index_set+k)*psi_in(index_set+k)
           !write(*,*) "k","u",index_set+1,k,u(index_set+1),T_mat(index_set+1,index_set+k),psi(l,index_set+k)
         ENDDO
         !write(*,*) "index",n,i,index_set+1,Nr,u(index_set+1)

         DO i = 2 , NDVR
           u(index_set+i) = 0.0d0*xi;
           do k = 1 , NDVR
             u(index_set+i) = u(index_set+i) +  D1_mat(index_set+i,index_set+k)*psi_in(index_set+k)
           ENDDO
           !write(*,*) "index",n,i,index_set+i,Nr,u(index_set+i)
         ENDDO

       endif


       psi_out(:) = u(:);



     deallocate(u);

END SUBROUTINE apply_Dr





!! Applies the Hamiltonian on the ket psi_inter_in and the result is psi_inter_out
subroutine apply_H(psi_out,psi_in,A_champ,x,Nx)

  IMPLICIT NONE;

  integer, intent(in) :: Nx;
  complex*16,dimension(:),Allocatable, intent(in) :: psi_in;
  complex*16,dimension(:),Allocatable :: psi_out;
  complex*16,dimension(:), Allocatable, intent(in) :: x;
  real*8, intent(in) :: A_champ;
  integer :: i,l,k,j,Nstop,k_min,Dk;
  Complex*16 :: somme_inter;

  psi_out(1:Nx) = 0.0d0*xi;

  psi_inter(1:Nx) = 0.0d0*xi

  !!! Apply Hmix
         !! This computes a vector containing the derivative in r of psi_in for each l
         CALL apply_Dr(psi_out,psi_in,x,Nx)


  !! Apply H0


        psi_inter(1:Nx) = 0.0d0*xi
        CALL apply_H0(psi_inter,psi_in,x,Nx)

         do j = 1 , Nx
            psi_out(j) = -xi*A_champ*psi_out(j) + psi_inter(j)
         enddo




end subroutine apply_H



!! Calculates the scalar product between two kets
FUNCTION scalar_product(psi_1,psi_2,x,Nx) result(somme)

  IMPLICIT NONE;

  integer, intent(in) :: Nx;
  complex*16,dimension(:),Allocatable, intent(in) :: psi_1,psi_2;
  complex*16,dimension(:),Allocatable, intent(in) :: x;
  complex*16 :: somme;
  integer :: i,l;


	       somme = 0.0d0*xi;
          do i = 1 , Nx
            !scalar_product = scalar_product + 0.5d0*real(psi_2(l,i)*dconjg(psi_1(l,i)) + psi_2(l,i+1)*dconjg(psi_1(l,i+1)))*(r(i+1)-r(i));
            somme = somme + conjg(psi_1(i))*psi_2(i);
          enddo


END FUNCTION scalar_product




FUNCTION error_test(psi_1,psi_2,Nx) result(error)
  IMPLICIT NONE;

  integer, intent(in) :: Nx;
  complex*16,dimension(:),Allocatable, intent(in) :: psi_1,psi_2;
  real*8 :: error,test;
  integer :: i,l;


	       error = 0.0d0;

          do i = 1 , Nx

            test = sqrt(real(conjg(psi_1(i)-psi_2(i))*(psi_1(i)-psi_2(i)),8));
            if(test .gt. error) then
              error = test;
            endif

          enddo


END FUNCTION error_test


end module prop
