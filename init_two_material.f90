subroutine init_two_material
use init_tran
use global_cont
use global
implicit none
integer i,j,k
double precision f_eta,dx,p


    nv=5
    jx=200
    dlx=5.0d0
	dx=dlx/jx
    allocate(U(-nv:jx+nv,0:3))
    allocate(Uo(-nv:jx+nv,0:3))
    allocate(X(-nv:jx+nv))
    
    u=0


    !TT=5.d-6
    TT=2.d-3
    !TT=1.d-6
	SF=0.2d0
!	SF=0.1d0

    rho1=893d0
    u2=0
    p1=1.d-6
	sxx1=0
    
    rho2=278.5d0
    u1=90.0d0
    !u1=50.d0
    p2=1.d-6
    sxx2=0

	  inter = int(2.5d0/dx)
    do i=-nv,jx+nv
	x(i)=(i-0.5d0)*dx 
        if(i.ge.inter)then !X(i-1/2)
			call state_choose(2)
			uo(i,0)=rho2
			uo(i,1)=u2
			uo(i,2)=p2
			uo(i,3)=0
			call trans_ue_to_u(uo(i,:),u(i,:))
        else
			call state_choose(1)
			uo(i,0)=rho1
			uo(i,1)=u1
			uo(i,2)=p1
			uo(i,3)=0
			call trans_ue_to_u(uo(i,:),u(i,:))
        endif
    enddo

end 
