subroutine init_Piston_problem
use init_tran
use global_cont
use global
implicit none
integer i,j,k
double precision f_eta,dx,p

		Y0=9.d7
		rho0=8930
		gamma0=2.d0
		miu=4.5d10
		a0=3940
		pi=6*dasin(0.5d0)
		s0=1.49

    nv=4
    jx=800
    dlx=1.0
	dx=dlx/jx
    allocate(U(-nv:jx+nv,0:3))
    allocate(X(-nv:jx+nv))
    
    u=0


    TT=1.5d-4
	SF=0.05d0

    rho2=8930d0
    u2=20.d0
    p2=1.d5
	sxx2=0
    
    rho1=8930d0
    u1=0.d0
    p1=1.d5
    sxx1=0

    do i=-nv,jx+nv
	!if (i>10)then
			x(i)=i*dx 
	        U(i,0)=rho1
            U(i,1)=rho1*u1
            U(i,2)=(p1-rho0*a0**2*f_eta(rho1))/(rho0*gamma0)*rho1+0.5d0*rho1*u1*u1
			U(i,3)=0
	!	else 
	!ji!		x(i)=i*dx 
       !     U(i,0)=rho2
       !j     U(i,1)=rho2*u2
        !    U(i,2)=(p2-rho0*a0**2*f_eta(rho2))/(rho0*gamma0)*rho2+0.5d0*rho2*u2*u2
		!	U(i,3)=0
		!endif
    enddo
end 
