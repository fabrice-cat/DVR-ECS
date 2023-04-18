module util
INTEGER:: Nlanczos;
real*8 :: x_ecs,theta_ecs,xmax,dx,E_guess,omega,F0,CV,zc,t,phi_CEP,t_on,t_off,dt;
INTEGER :: Nseq,NDVR,Nx;
real*8 :: F02,t_on2,t_off2,phi_CEP2,omega2,t_delay;
real*8 :: E_start,E_step,dE,theta_start,d_theta;
integer :: NE,Nc,Nt,Nexp,param_env;
integer :: Nc2,param_env2;
real*8 :: theta,phi,time,E_cut_sup,E_cut_inf;
Complex*16,PARAMETER :: xi=(0.0d0,1.0d0);
real*8,PARAMETER :: pi = acos(-1.0d0);
real*8 :: EI_H = 13.6d0;
Complex*16,dimension(:),Allocatable :: psi0,psi,psi_inter,psi_inter2,psi_inter_loc,psi_inter2_loc,psi_inter3_loc;
INTEGER,dimension(:),Allocatable :: max_index_PHP;
complex*16,dimension(:),Allocatable :: x,xseq,Pot_mat,weight,x_rb,weight_rb;
complex*16,dimension(:,:),Allocatable :: xi_seq,wi_seq;
complex*16,dimension(:,:,:),Allocatable :: T_reduced,D1_reduced,proj_PHP;
complex*16,dimension(:,:,:),Allocatable :: expm_mat;
complex*16,dimension(:,:),Allocatable :: T_mat,H0_mat,D1_mat;
end module util
