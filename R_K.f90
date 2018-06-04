 subroutine R_K
 use global
 use constant
 implicit none
 integer i,j,k
 double precision ut2(-nv:jx+nv,0:2)
 double precision ut3(-nv:jx+nv,0:2)
 double precision ut1(-nv:jx+nv,0:2)
 double precision ut4(-nv:jx+nv,0:2)
   
 !   call bound(u)
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

    call bound(u)

	x(:)=x(:)+dt*h_G(:,1)/h_G(i,0)
    call space(u,x,dt)
    u=u+dt*f
    call bound(u)
endsubroutine
     
