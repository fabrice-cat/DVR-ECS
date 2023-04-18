
!%%%%%%%%%%%%%%%%%  Program  %%%%%%%%%%%%%%
PROGRAM TDSE_3D

  use util; use tools; use tridiag; use prop; use potential;

  IMPLICIT NONE

  complex*16 :: ortho,somme,energy;

  real*8 :: norme;
  integer :: i,l,j;

  Complex*16,dimension(:),Allocatable :: psi_out;

  !%%%% Open the file param.txt and read the parameters of the calculation
  call Init_data();

  !%%%% Allocate matrices
  Allocate(xseq(0:Nseq)); Allocate(xi_seq(1:Nseq,1:NDVR)); Allocate(wi_seq(1:Nseq,1:NDVR));

  Allocate(T_reduced(1:Nseq,1:NDVR,1:NDVR),D1_reduced(1:Nseq,1:NDVR,1:NDVR));



  Nx = (Nseq)*(NDVR-1)+1
  Nx = Nx - 2; !! -2 for 0 at the two boundaries and -1 for the boundary at r = 0

  Allocate(T_mat(1:Nx,1:Nx),H0_mat(1:Nx,1:Nx),D1_mat(1:Nx,1:Nx),Pot_mat(1:Nx));
  Allocate(psi0(1:Nx),x(1:Nx),weight(1:Nx),x_rb(1:Nx),weight_rb(1:Nx));

  Allocate(psi(1:Nx),psi_out(1:Nx),psi_inter(1:Nx),psi_inter2(1:Nx));
  Allocate(psi_inter_loc(1:Nx),psi_inter2_loc(1:Nx),psi_inter3_loc(1:Nx));

!  Allocate(psiexc(1:Nr));

  call Init_matrix();

  CV = 1.0d-30;

  ! write(*,*) "symmetry test, Nx",Nx
  ! write(*,*) "(10,10) - (4490,4490)",H0_mat(10,10),H0_mat(4490,4490)
  ! write(*,*)
  ! write(*,*) "x(2500) - x(2000)", x(2500), x(2000)
  ! write(*,*)
  ! write(*,*) "Pot(2500) - Pot(1999)", Pot_mat(2500), Pot_mat(2000)
  ! write(*,*)
  ! write(*,*) "(2500,2500) - (2000,2000)",T_mat(2500,2500),T_mat(2000,2000)
  ! write(*,*)
  ! write(*,*) "(2500,2501) - (2000,1999)",T_mat(2500,2501),T_mat(2000,1999)
  ! write(*,*)
  ! write(*,*) "(10,10) - (4490,4490)",T_mat(10,10),T_mat(4490,4490)
  ! write(*,*)
  ! write(*,*) "(10,11) - (4490,4489)",T_mat(10,11),T_mat(4490,4489)
  ! write(*,*)
  ! write(*,*) "(10,9) - (4490,4491)",T_mat(10,9),T_mat(4490,4491)
  ! write(*,*)
  ! write(*,*) "x(900) - x(3600)", x(900), x(3600)
  ! write(*,*)
  ! write(*,*) "(900,900) - (3600,3600)",T_mat(900,900),T_mat(3600,3600)
  ! write(*,*)
  ! write(*,*) "(900,901) - (3600,3599)",T_mat(900,901),T_mat(3600,3599)
  ! write(*,*)
  ! write(*,*) "(900,899) - (3600,3601)",T_mat(900,899),T_mat(3600,3601)
  ! write(*,*)

energy = E_initialise(psi0,x,Nx,E_guess,CV)

 write(*,*) energy

!! INit the wave function
 psi(:) = 0.0d0*xi;
 do i = 1 , Nx
    psi(i) = psi0(i);
enddo

 write(*,*) "wave fucntion",psi(5),psi0(5)

 somme = 0.0d0
 Do i = 1 ,Nx
   somme = somme + psi(i)*psi(i);
 ENDDO
  write(*,*) "Norme direct",somme;
  somme = 0.0d0
  Do i = 1 ,Nx
    somme = somme + dconjg(psi(i))*psi(i);
  ENDDO
  write(*,*) "Norme conjugates",somme;

  write(*,*) "with apply H0";
  energy = E_calculation(psi,x,Nx)
  write(*,*) energy;



  !PAUSE
  CALL apply_Dr(psi_out,psi,x,Nx)


  open(UNIT=1,FILE="derpsi_gs.dat",FORM="FORMATTED",action='write');

      do i = 1,Nx
        !write(UNIT=1,fmt='(XE13.6,XE13.6,XE13.6)') r(i),real(Xin(i)/sqrt(weight(i))),dimag(Xin(i)/sqrt(weight(i)));
        write(UNIT=1,fmt='(XE13.6,XE13.6,XE13.6)') x(i),real(psi_out(i)/sqrt(weight(i))),0.0d0;
      enddo

  close(1);

!  write(*,*) energy,E_guess,100.0d0*abs((energy-E_guess)/E_guess),l_init,m_init;


!  psi(abs(m_init):lmax,1:Nr) = 0.0d0*xi;
!  psi(l_init,1:Nr) = psi0(1:Nr);


 !call propagation(F0,omega,phi_CEP,Nc,Nt,t_on,t_off,t_on2,t_off2,param_env,param_env2,t_delay,psi,Nx,Nexp,x,Nlanczos);


 open(UNIT=1,FILE="psi_end_mod.dat",FORM="FORMATTED",action='write');

  do i = 1 , Nx

        somme =  psi(i)*dconjg(psi(i));

     write(UNIT=1,fmt='(XE13.6,XE13.6)') x(i),somme/weight(i);
  enddo


 ! close(1)

 ! write(*,*) "Save the wavefunction ..."
 ! open(UNIT=1,FILE="psi_end_re.dat",FORM="FORMATTED",action='write');
 ! open(UNIT=2,FILE="psi_end_im.dat",FORM="FORMATTED",action='write');
 ! open(UNIT=3,FILE="r_num.dat",FORM="FORMATTED",action='write');

 ! do i = 1,Nr
 !
 !    do l = abs(m_init),Lmax
 !
 !      if(l.eq.Lmax) then
 !       write(UNIT=1,fmt='(XE13.6)') real(psi(l,i));
 !       write(UNIT=2,fmt='(XE13.6)') dimag(psi(l,i));
 !      else
 !       write(UNIT=1,fmt='(XE13.6)',advance='no') real(psi(l,i));
 !       write(UNIT=2,fmt='(XE13.6)',advance='no') dimag(psi(l,i));
 !      endif
 !
 !    enddo
 !
 !    write(UNIT=3,fmt='(XE13.6)') r(i);
 ! enddo


 ! close(1); close(2);close(3)


  write(*,*) "Photoelectron calculation ..."

  	 ortho = 0.0d0*xi;
     do i = 1 , Nx
            ortho = ortho + psi0(i)*psi(i);
     enddo


  psi(1:Nx) = psi(1:Nx) - ortho*psi0(1:Nx);

!  call window_analysis(psi,r,Nr,Lmax,m_init,zc,E_start,E_step,dE,NE,theta_start,d_theta);
!zc = 1.0d0
!call phoelectron_spectrum_CW(psi,r,Nr,lmax,m_init,zc,E_step,NE,theta_start,d_theta)


deallocate(xseq,wi_seq,xi_seq);





END PROGRAM TDSE_3D
