subroutine u_f(u4,u5)
use constant
double precision u4(0:2)
double precision u5(0:2)
double precision uu,uv,c,p,rho,e

    rho=u4(0)
    uu=u4(1)/rho
    e=u4(2)
    p=(GAMA-1)*(e-0.5*rho*(uu*uu+uv*uv))

    u5(0)=rho*uu
    u5(1)=rho*uu*uu+p
    u5(2)=(e+p)*uu

end
