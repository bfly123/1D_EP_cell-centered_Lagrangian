subroutine space(u4,x4,dt)
use global
use constant
implicit none

integer i,j,k
double precision u4(-nv:jx+nv,0:3)
double precision x4(-nv:jx+nv)
double precision ul(-nv:jx+nv,0:3)
double precision ur(-nv:jx+nv,0:3)
double precision h(-nv:jx+nv,0:3)
integer ind(-nv:jx+nv,0:2)
    
    ind=0
    call SW_splitting(u4)
	call  upwind(nv,jx,u4,ul,ur)
	call  HLLC_EP(nv,jx,ul,ur,h)
    call Output1(h)  
    do k=0,3
     do i=0,jx
		dx=x(i)-x(i-1)
        f(i,k)=-1.d0/dx*(h(i,k)-h(i-1,k))
     enddo
     
    enddo
  

    end subroutine
 
 
