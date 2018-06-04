module global        
implicit none

integer           M1                                !x方向网格数
integer           nv                                !虚拟网格点数
integer           JX                                !x方向网格数 
integer           nk                                !识别器类型
integer           kind_problem                      !问题类型
integer           kind_split
double precision dLx                               !x方向计算区域长度  
double precision dx                                !x方向单元尺寸
double precision SF                                !CFL数
double precision TT                                !计算时间
double precision t                                 !计算中时间标记
double precision dt                                !时间步长


double precision,allocatable::U(:,:)            !通量Ut=F
double precision,allocatable::f(:,:)            !
double precision,allocatable::fp(:,:)           !x方向正通量
double precision,allocatable::fm(:,:)           !x方向负
double precision,allocatable::x(:)

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
double precision::Y0=1.4d0
double precision::rho0=1.4d0
double precision::gamma0=1.4d0
double precision::miu=1.4d0
double precision::a0=1.4d0
double precision::pi=6*dasin(0.5d0)
end module
