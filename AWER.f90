 subroutine ADER(dt)
	  using global
	  implicit none
	  integer i,j,k

   call bound(u)
	do i=-nv+1,nv+jx
		dx1(i)= x(i)-x(i-1)
    	udx(i,0:2)=u(i,0:2)*dx1(i)
    	udx(i,3)=u(i,3)
		urho(i)=u(i,0)
	enddo

	  call space_ADER(U,f,u_half)

	  call Gauss(f,dt)

		x1(:)=x(:)+dt*u_half(:)
	   udx(:,0:2) = udx(:,0:2)+f(:,0:2)

	  end








