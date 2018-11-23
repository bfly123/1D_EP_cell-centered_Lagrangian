 subroutine R_K(t,dt)
	 use global
 implicit none
 integer i,j,k
 double precision udx1(-nv:jx+nv,0:3)
 double precision udx2(-nv:jx+nv,0:3)
 double precision udx3(-nv:jx+nv,0:3)
 double precision udx4(-nv:jx+nv,0:3)

 double precision x1(-nv:jx+nv)
 double precision x2(-nv:jx+nv)
 double precision x3(-nv:jx+nv)
 double precision x4(-nv:jx+nv)

 double precision sxx1(-nv:jx+nv)
 double precision sxx2(-nv:jx+nv)
 double precision sxx3(-nv:jx+nv)
 double precision sxx4(-nv:jx+nv)

 double precision udx(-nv:jx+nv,0:3)
 double precision src(-nv:jx+nv,0:3)
 double precision urho(-nv:jx+nv)
 double precision f(-nv:jx+nv,0:3)
 double precision u_half(-nv:jx+nv)
 double precision dx1(-nv:jx+nv)
 double precision dx,dt,dxmin,fgamma,t
 integer time_order

  time_order=3

   call bound(u)
	do i=-nv+1,nv+jx
		dx1(i)= x(i)-x(i-1)
    	udx(i,0:2)=u(i,0:2)*dx1(i)
    	udx(i,3)=u(i,3)
		urho(i)=u(i,0)
	enddo
  select case(time_order)
case(5)
	  !*************1**********
	call space(U,f,u_half)

	!x1(:)=x(:)+dt*u_half(:)
	x1(:)=x(:)+dt*u_half(:)

	do i=-nv+1,nv+jx
    !	udx1(i,3)=udx(i,3)-dt*f(i,3)/dx1(i)
		dx1(i)= x1(i)-x1(i-1)
    	udx1(i,0:2)=udx(i,0:2)-dt*f(i,0:2)
    	!u(i,3)=udx1(i,3)
		u(i,0:2) = udx1(i,0:2)/dx1(i)
		call rho_sxx(urho(i),u(i,0),u(i,3))
	!	udx1(i,3)=u(i,3)*dx
		urho(i)=u(i,0)
	enddo
!	x(:)=x1(:)

	  !*************2**********
    call bound(u)
	call space(U,f,u_half)
	x2(:)=3.d0/4*x(:)+1.d0/4*x1(:)+1.d0/4*dt*u_half(:)

	do i=-nv+1,nv+jx
    	!udx2(i,3)=3.d0/4*udx(i,3)+1.d0/4*udx1(i,3)-1.d0/4*dt*f(i,3)/dx1(i)
		dx1(i)= x2(i)-x2(i-1)
    	udx2(i,0:2)=3.d0/4*udx(i,0:2)+1.d0/4*udx1(i,0:2)-1.d0/4*dt*f(i,0:2)
    	!u(i,3)=udx2(i,3)
		u(i,0:2) = udx2(i,0:2)/dx1(i)
		call rho_sxx(urho(i),u(i,0),u(i,3))
		urho(i)=u(i,0)
		!u(i,3) = fgamma(u(i,3))
	enddo
!
!	  !*************3**********

    call bound(u)
	call space(U,f,u_half)
		x(:)=1.d0/3*x(:)+2.d0/3*x2(:)+2.d0/3*dt*u_half(:)
	do i=-nv+1,nv+jx
    	!udx(i,3)=1.d0/3*udx(i,3)+2.d0/3*udx2(i,3)-2.d0/3*dt*f(i,3)/dx1(i)
		dx1(i)= x(i)-x(i-1)
    	udx(i,0:2)=1.d0/3*udx(i,0:2)+2.d0/3*udx2(i,0:2)-2.d0/3*dt*f(i,0:2) -dt*dx1(i)*src(i,0:2)
    	!u(i,3)=udx(i,3)
		u(i,0:2) = udx(i,0:2)/dx1(i)
		call rho_sxx(urho(i),u(i,0),u(i,3))
	!urho(i)=u(i,0)
		!u(i,3) = fgamma(u(i,3))
	enddo
!
  case(3)

call source1(t,src)
	  !*************1**********
	call space(U,f,u_half)
	!x1(:)=x(:)+dt*u_half(:)
!	call bound_v(u_half)
	x1(:)=x(:)+dt*u_half(:)

	do i=-nv+1,nv+jx
    	udx1(i,3)=udx(i,3)-dt*f(i,3)/dx1(i)
		dx1(i)= x1(i)-x1(i-1)
    	udx1(i,0:2)=udx(i,0:2)-dt*f(i,0:2) +dt*src(i,0:2)
    	u(i,3)=udx1(i,3)
		u(i,0:2) = udx1(i,0:2)/dx1(i)

		u(i,3) = fgamma(u(i,3))
		call trans_u_to_ue(u(i,:),uo(i,:))
		!call trans_u_to_ue()
	!	call rho_sxx(urho(i),u(i,0),u(i,3))
	!	urho(i)=u(i,0)
	!	udx1(i,3)=u(i,3)*dx
	enddo
!	x(:)=x1(:)

	  !*************2**********

call source1(t+dt,src)
    call bound(u)
	call space(U,f,u_half)
!	call bound_v(u_half)
	x2(:)=3.d0/4*x(:)+1.d0/4*x1(:)+1.d0/4*dt*u_half(:)

	do i=-nv+1,nv+jx
    	udx2(i,3)=3.d0/4*udx(i,3)+1.d0/4*udx1(i,3)-1.d0/4*dt*f(i,3)/dx1(i)
		dx1(i)= x2(i)-x2(i-1)
    	udx2(i,0:2)=3.d0/4*udx(i,0:2)+1.d0/4*udx1(i,0:2)-1.d0/4*dt*f(i,0:2)  +1.d0/4*dt*src(i,0:2)
    	u(i,3)=udx2(i,3)
		u(i,0:2) = udx2(i,0:2)/dx1(i)
		!call rho_sxx(urho(i),u(i,0),u(i,3))
		!urho(i)=u(i,0)
		u(i,3) = fgamma(u(i,3))
		call trans_u_to_ue(u(i,:),uo(i,:))
	enddo
!

	call source1(t+5.d0/4*dt,src)

	!write(*,*) src
	!pause
!	  !*************3**********
    call bound(u)
	call space(U,f,u_half)
!	call bound_v(u_half)
		x(:)=1.d0/3*x(:)+2.d0/3*x2(:)+2.d0/3*dt*u_half(:)
	do i=-nv+1,nv+jx
    	udx(i,3)=1.d0/3*udx(i,3)+2.d0/3*udx2(i,3)-2.d0/3*dt*f(i,3)/dx1(i)
		dx1(i)= x(i)-x(i-1)
    	udx(i,0:2)=1.d0/3*udx(i,0:2)+2.d0/3*udx2(i,0:2)-2.d0/3*dt*f(i,0:2) +2.d0/3*dt*src(i,0:2)
    !	udx(i,0:2)=udx(i,0:2)+dt*dx1(i)*src(i,0:2)
    	u(i,3)=udx(i,3)
		u(i,0:2) = udx(i,0:2)/dx1(i)
		!call rho_sxx(urho(i),u(i,0),u(i,3))
		!urho(i)=u(i,0)
		u(i,3) = fgamma(u(i,3))

		call trans_u_to_ue(u(i,:),uo(i,:))
	enddo
!
case(1)

call bound(u)
call space(U,f,u_half)
call source1(t,src)
	x(:)=x(:)+dt*u_half(:)
	do i=-nv+1,nv+jx
    	udx(i,3)=udx(i,3)-dt*f(i,3)/dx1(i)
		dx1(i)= x(i)-x(i-1)
    	udx(i,0:2)=udx(i,0:2)-dt*f(i,0:2) +src(i,0:2)*dt
    	u(i,3)=udx1(i,3)
		u(i,0:2) = udx(i,0:2)/dx1(i)
		u(i,3) = fgamma(u(i,3))
	enddo
!o
case(4)
 !*************1**********
	call space(U,f,u_half)
	x1(:)=x(:)+1.d0/2*dt*u_half(:)
	do i=-nv+1,nv+jx
		dx= x1(i)-x1(i-1)
    	udx1(i,:)=udx(i,:)-1.d0/2*dt*f(i,:)
		u(i,:) = udx1(i,:)/dx
		u(i,3) = fgamma(u(i,3))
		!udx1(i,3)=u(i,3)*dx
	enddo

	  !*************2**********
    call bound(u)
	call space(U,f,u_half)
	x2(:)=x(:)+1.d0/2*dt*u_half(:)

	do i=-nv+1,nv+jx
		dx= x2(i)-x2(i-1)
    	udx2(i,:)=udx(i,:)-1.d0/2*dt*f(i,:)
		u(i,:) = udx2(i,:)/dx
		u(i,3) = fgamma(u(i,3))
		udx2(i,3)=u(i,3)*dx
	enddo

	  !*************3**********
    call bound(u)
	call space(U,f,u_half)
	x3(:)=x(:)+dt*u_half(:)
	do i=-nv+1,nv+jx
		dx= x(i)-x(i-1)
    	udx3(i,:)=udx(i,:)-dt*f(i,:)
		u(i,:) = udx3(i,:)/dx
		u(i,3) = fgamma(u(i,3))
		udx3(i,3)=u(i,3)*dx
	enddo

	  !*************4**********
  call bound(u)
	call space(U,f,u_half)
	x(:)=1.d0/3*(-x(:)+x1(:)+2*x2(:)+x3(:))+1.d0/6*dt*u_half(:)
	do i=-nv+1,nv+jx
		dx= x(i)-x(i-1)
    	udx(i,:)=1.d0/3*(-udx(i,:)+udx1(i,:)+2*udx2(i,:)+udx3(i,:))-1.d0/6*dt*f(i,:)
		u(i,:) = udx(i,:)/dx
		u(i,3) = fgamma(u(i,3))
	enddo

end select

	call bound(U)


 !   call space(u)
 !   ut1=u+dt*f
 !   
 !   call bound(ut1)
 !   call space(ut1)
 !       ut2=3.d0/4*u+1.d0/4*ut1+1.d0/4*dt*f
 ! 
	!call bound(ut2)
	!call space(ut2)
 !
 !      u=1.d0/3*u+2.d0/3*ut2+2.d0/3*dt*f
 !   call bound(u)      


 !call bound(u)
 !call space(u)
 !ut1=u+1.d0/2*dt*f
 !call bound(ut1)
 !call space(ut1)
 !ut2=u+1.d0/2*dt*f
 !call bound(ut2)
 !call space(ut2)
 !ut3=u+dt*f
 !call bound(ut3)
 !call space(ut3)
 !u=1.d0/3*(-u+ut1+2*ut2+ut3)+1.d0/6*dt*f
 !call bound(u)
endsubroutine
 

subroutine rho_sxx(rho1,rho2,sxx)
	use global_cont
	implicit none 
	double precision rho1,rho2,sxx,fgamma

	sxx=-4.d0/3*miu*log(rho2/rho1)+sxx
	sxx=fgamma(sxx)
end



