subroutine init_entropy_wave
use init_tran
use global_cont
use global
implicit none
integer i,j,k
double precision f_eta,dx,p,a,s,e0,b
double precision ue(0:3)

		Y0=  9.d10
		rho0=8930
		gamma0=2.d0
		miu=4.5d10
		a0=3940
		pi=6*dasin(0.5d0)
		s0=1.49

		s= 6.d4
		e0= 0.1
		a = 1.d4
		b=0.2

    nv=4
    jx=200
    dlx=1.0
	dx=dlx/jx
    allocate(U(-nv:jx+nv,0:3))
    allocate(Uo(-nv:jx+nv,0:3))
    allocate(X(-nv:jx+nv))
    
    u=0

    TT=0.1d0
	SF=0.5d0

	do i=-nv,jx+nv
		x(i)=(i+0.5d0)*dx 
	enddo
    do i=-nv+1,jx+nv-1
	!if (i>10)then
		uo(i,0)= rho0*(1-b*sin(2*pi*(x(i)+x(i-1))/2))
		uo(i,1)= a
		uo(i,2)=1.1d0*a0**2*rho0
		uo(i,3)=s *sin(pi*(x(i)+x(i-1)))
		call trans_ue_to_u(uo(i,:),u(i,:))
	enddo
end 
