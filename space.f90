subroutine space(u4,f4,u_half)
use global
use global_cont
implicit none

integer i,j,k
double precision u4(-nv:jx+nv,0:3)
double precision ul(-nv:jx+nv,0:3)
double precision ulo(-nv:jx+nv,0:3)
double precision uo(-nv:jx+nv,0:3)
double precision ur(-nv:jx+nv,0:3)
double precision uro(-nv:jx+nv,0:3)
double precision h(-nv:jx+nv,0:3)
double precision f4(-nv:jx+nv,0:3)
double precision u_half(-nv:jx+nv)
double precision  rho,uu,p,f_eta
integer kind_split

kind_split =2

select case(kind_split) 
case(1)
	!call  upwind(nv,jx,u4(:,i),ul(:,i),ur(:,i))
		do i=-nv,jx+nv
			rho=u4(i,0)
        uu=U(i,1)/rho
	  p = (u(i,2)/rho- 0.5*uu**2)*rho0*gamma0+rho0*a0**2*f_eta(rho)
	  uo(i,0)=rho
	  uo(i,1)=uu
	  uo(i,2)=p
	  uo(i,3)=u4(i,3)
	enddo
	do i=0,3
	call  WENO5_new(nv,jx,uo(:,i),ulo(:,i),uro(:,i))
	enddo

case(2)

	do i=0,3
	call  WENO5_new(nv,jx,u(:,i),ul(:,i),ur(:,i))
	enddo

	call  HLLC_EP(nv,jx,u4(:,1)/u4(:,0),ul,ur,h,u_half)
case (3)
		call  LF_splitting(u4,ul,ur)
		do i =0,3
			call WENO5_new_LF(nv,jx,ul(:,i),ur(:,i),h(:,i))
		enddo
		u_half(:) = u4(:,1)/u4(:,0)

	endselect


    do k=0,3
     do i=-nv+2,jx+nv-2
        f4(i,k)=h(i,k)-h(i-1,k)
     enddo
    enddo

    end subroutine
 
 
