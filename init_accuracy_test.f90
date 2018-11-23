subroutine init_accuracy_test
use init_tran
use global_cont
use global
implicit none
integer i,j,k
double precision f_eta,dx,p,a,s,e0,b
double precision ue(0:3)

		Y0=0
		rho0=2
		gamma0=2.d0
		miu=4.5d5
		a0=394
		pi=6*dasin(0.5d0)
		s0=1.49

		s= 6.d2
		e0= 1.d-1
		a = 1.d3
		b=0.2

    nv=6
    jx=2000
    dlx=1.0
	dx=dlx/jx
    allocate(U(-nv:jx+nv,0:3))
    allocate(Uo(-nv:jx+nv,0:3))
    allocate(X(-nv:jx+nv))
    
    u=0

    TT=1.d0
	SF=0.6d0

	do i=-nv,jx+nv
		x(i)=(i+0.5)*dx 
	enddo

    do i=-nv,jx+nv
	!if (i>10)then
			uo(i,0)= rho0*(1+b*sin(pi*(x(i)+x(i-1))))
			uo(i,1) =   1+b*sin(pi*(x(i)+x(i-1)))
			uo(i,2)= 1! -s*sin(pi*(x(i)+x(i-1))) + 1.5d0*rho0*a0**2 !( 1-b*sin(pi*(x(i)+x(i-1))))
		!	(e0+0.5d0*a**2)*u(i,0)
			uo(i,3)=10!  s*sin(pi*(x(i)+x(i-1)))
			call trans_ue_to_u(uo(i,:),u(i,:))
	enddo
!call output1(uo)
!	pause
end 


