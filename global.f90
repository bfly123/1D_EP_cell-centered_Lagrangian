module global        
implicit none
integer           nv                                !�����������
integer           JX                                !x���������� 
integer           kind_problem                                !x���������� 
double precision dLx                               !x����������򳤶�  
double precision TT                                !����ʱ��
double precision SF                                !����ʱ��
double precision, allocatable::U(:,:)            !ͨ��Ut=F
double precision, allocatable::X(:)
end module

module init_tran   !��ʼ�߽�ֵ
double precision rho1                              !���߽�����
double precision rho2                              !�Ҳ�߽�����
double precision p1
double precision p2
double precision u1
double precision u2
double precision sxx2
double precision sxx1
end module

module global_cont
implicit none
double precision::Y0
double precision::rho0
double precision::gamma0
double precision::miu
double precision::a0
double precision::pi
double precision::s0
end module

module mymkl
include 'mkl.fi' 
!include  'lapack.f90'
end module
