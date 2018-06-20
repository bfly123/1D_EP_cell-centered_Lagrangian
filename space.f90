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
    
	do i =0,3	
	!call  upwind(nv,jx,u4(:,i),ul(:,i),ur(:,i))
	call  WENO5_new(nv,jx,u4(:,i),ul(:,i),ur(:,i))
	enddo


	call  HLLC_EP(nv,jx,u4(:,1)/u4(:,0),ul,ur,h,u_half)
    !call Output1(h)  

    do k=0,3
     do i=-nv+2,jx+nv-2
        f4(i,k)=h(i,k)-h(i-1,k)
     enddo
    enddo

    end subroutine
 
 
