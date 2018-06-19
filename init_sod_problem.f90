subroutine init_sod_problem
use init_tran
use global_cont
use global
implicit none
integer i,j,k
double precision f_eta,dx,p

    nv=4
    jx=800
    dlx=5.0d-2
	dx=dlx/jx
    allocate(U(-nv:jx+nv,0:3))
    allocate(X(-nv:jx+nv))
    
    u=0


    TT=1.d-6
    TT=0
	SF=0.45d0

    rho2=2785d0
    u2=0
    p2=1.d-6
	sxx2=0
    
    rho1=2785d0
    u1=800.d0
    p1=1.d-6
    sxx1=0

    do i=-nv,jx+nv
	x(i)=(i-0.5d0)*dx 
        if(x(i).ge.5.d-3)then !X(i-1/2)
            U(i,0)=rho2
            U(i,1)=rho2*u2
            U(i,2)=(p2-rho0*a0**2*f_eta(rho2))/(rho0*gamma0)*rho2+0.5d0*rho2*u2*u2
			p = (u(i,2)/rho2- 0.5d0*u2*u2)*rho0*gamma0+rho0*a0**2*f_eta(rho2)
			U(i,3)=0
        else
            U(i,0)=rho1
            U(i,1)=rho1*u1
            !U(i,2)=(p1-rho0*a0**2*f_eta(rho1))/(rho0*gamma0)*rho1+0.5d0*rho1*u1*u1
			!p = (u(i,2)/rho1- 0.5d0*u1*u1)*rho0*gamma0+rho0*a0**2*f_eta(rho1)
            p=1.0d-6/5570.d0+320000.d0
			p = (p- 320000.d0)*5570.d0
			U(i,3)=0
        endif
		write(*,*)p,rho0*gamma0,0.5d0*u1*u1,rho0*a0*a0*f_eta(rho1)
    enddo
end 
