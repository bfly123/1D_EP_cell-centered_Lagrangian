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
double precision::Y0=3.d8
double precision::rho0=2785
double precision::gamma0=2.d0
double precision::miu=2.76d10
double precision::a0=5328
double precision::pi=6*dasin(0.5d0)
double precision::s0=1.338
end module

