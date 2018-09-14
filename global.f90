module global        
implicit none
integer           nv                                !虚拟网格点数
integer           JX                                !x方向网格数 
integer           kind_problem                                !x方向网格数 
double precision dLx                               !x方向计算区域长度  
double precision TT                                !计算时间
double precision SF                                !计算时间
double precision, allocatable::U(:,:)            !通量Ut=F
double precision, allocatable::X(:)
end module

module init_tran   !初始边界值
double precision rho1                              !左侧边界条件
double precision rho2                              !右侧边界条件
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
