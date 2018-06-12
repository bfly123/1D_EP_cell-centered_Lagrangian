      subroutine CFL(dt)
		  use global
      implicit none
      integer i
	  double precision dxmin
      
	  dxmin=x(1)-x(0)
      do i=1,Jx
		dx=x(i)-x(i-1)
		if( dx.le.dxmin ) then
			dxmin = dx
		endif
	enddo
	    
     dt=SF*dxmin

      end
