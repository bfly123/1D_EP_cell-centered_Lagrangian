subroutine init_accuracy_test
use init_tran
use global_cont
use global
implicit none
integer i,j,k
double precision f_eta,dx,p,a,s,e0,b
double precision ue(0:3)

		Y0=9.d10
		rho0=8930
		gamma0=2.d0
		miu=4.5d10
		a0=3940
		pi=6*dasin(0.5d0)
		s0=1.49

		s= 6.d7
		e0= 0.d-14
		a = 1.d4
		b=0.1

    nv=4
    jx=400
    dlx=1.0
	dx=dlx/jx
    allocate(U(-nv:jx+nv,0:3))
    allocate(X(-nv:jx+nv))
    
    u=0


    TT=0.01d0
	SF=0.1d0

    do i=-nv,jx+nv
	!if (i>10)then
			x(i)=i*dx 
			u(i,0)= rho0*(1-b*sin(2*pi*x(i)))
			u(i,1) = u(i,0)* a 
			u(i,2)= (e0+0.5d0*a**2)*u(i,0)
			u(i,3)= s*sin(2*pi*x(i))
	enddo
end 


