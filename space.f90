subroutine space(u4,f4,u_half)
use global
implicit none

integer i,j,k
double precision u4(-nv:jx+nv,0:3)
double precision ul(-nv:jx+nv,0:3)
double precision ur(-nv:jx+nv,0:3)
double precision h(-nv:jx+nv,0:3)
double precision f4(-nv:jx+nv,0:3)
double precision u_half(-nv:jx+nv)
    

    call SW_splitting(u4)
	call  upwind(nv,jx,u4(,0:3),ul(,0:3),ur(,0:3))
	call  HLLC_EP(nv,jx,ul,ur,h,u_half)
    call Output1(h)  

    do k=0,3
     do i=0,jx
        f(i,k)=h(i,k)-h(i-1,k)
     enddo
     
    enddo
  

    end subroutine
 
 
