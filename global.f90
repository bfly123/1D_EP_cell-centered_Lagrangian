module global        
implicit none
integer           nv                                !�����������
integer           JX                                !x���������� 
integer           kind_problem                                !x���������� 
double precision dLx                               !x����������򳤶�  
double precision TT                                !����ʱ��
double precision,allocatable::U(:,:)            !ͨ��Ut=F
double precision,allocatable::X(:)
end module

module init_tran   !��ʼ�߽�ֵ
double precision rou1                              !���߽�����
double precision rou2                              !�Ҳ�߽�����
double precision p1
double precision p2
double precision u1
double precision u2
end module

module global_cont
implicit none
double precision::Y0=9.d7
double precision::rho0=8930
double precision::gamma0=2.d0
double precision::miu=4.5d10
double precision::a0=3940
double precision::pi=6*dasin(0.5d0)
double precision::s0=1.49
end module

