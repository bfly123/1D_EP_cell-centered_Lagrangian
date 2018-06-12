module global        
implicit none
integer           nv                                !虚拟网格点数
integer           JX                                !x方向网格数 
integer           kind_problem                                !x方向网格数 
double precision dLx                               !x方向计算区域长度  
double precision TT                                !计算时间
double precision,allocatable::U(:,:)            !通量Ut=F
double precision,allocatable::X(:)
end module

module init_tran   !初始边界值
double precision rou1                              !左侧边界条件
double precision rou2                              !右侧边界条件
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

