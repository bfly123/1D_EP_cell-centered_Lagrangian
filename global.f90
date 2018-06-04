module global        
implicit none

integer           M1                                !x����������
integer           nv                                !�����������
integer           JX                                !x���������� 
integer           nk                                !ʶ��������
integer           kind_problem                      !��������
integer           kind_split
double precision dLx                               !x����������򳤶�  
double precision dx                                !x����Ԫ�ߴ�
double precision SF                                !CFL��
double precision TT                                !����ʱ��
double precision t                                 !������ʱ����
double precision dt                                !ʱ�䲽��


double precision,allocatable::U(:,:)            !ͨ��Ut=F
double precision,allocatable::f(:,:)            !
double precision,allocatable::fp(:,:)           !x������ͨ��
double precision,allocatable::fm(:,:)           !x����
double precision,allocatable::x(:)

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
double precision::Y0=1.4d0
double precision::rho0=1.4d0
double precision::gamma0=1.4d0
double precision::miu=1.4d0
double precision::a0=1.4d0
double precision::pi=6*dasin(0.5d0)
end module
