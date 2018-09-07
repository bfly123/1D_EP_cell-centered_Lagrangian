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
double precision AL(-nv:jx+nv,0:3,0:3)
double precision AR(-nv:jx+nv,0:3,0:3)
double precision f4(-nv:jx+nv,0:3)
double precision u_half(-nv:jx+nv)
double precision  rho,uu,p,f_eta
integer kind_split

kind_split =4

select case(kind_split) 
case(1)
		do i=-nv,jx+nv
			rho=u4(i,0)
        uu=U4(i,1)/rho
	  p = (u4(i,2)/rho- 0.5*uu**2)*rho0*gamma0+rho0*a0**2*f_eta(rho)
	  uo(i,0)=rho
	  uo(i,1)=uu
	  uo(i,2)=p
	  uo(i,3)=u4(i,3)
	enddo
	do i=0,3
!call  WENO5_new(nv,jx,uo(:,i),ulo(:,i),uro(:,i))

	!call  WENO3_new(nv,jx,uo(:,i),ulo(:,i),uro(:,i))
!	call  WENO3LIU_new(nv,jx,uo(:,i),ulo(:,i),uro(:,i))
	call  upwind(nv,jx,uo(:,i),ulo(:,i),uro(:,i))
	enddo
	do i=-nv,jx+nv
	 rho= ulo(i,0)
	  uu=ulo(i,1)
	  p=ulo(i,2)
	  ul(i,0)=rho
	  ul(i,1)=rho*uu
      Ul(i,2)=(p-rho0*a0**2*f_eta(rho))/(rho0*gamma0)*rho+0.5d0*rho*uu*uu
	 ul(i,3)= ulo(i,3)
	
	  rho= uro(i,0)
	  uu=uro(i,1)
	  p=uro(i,2)
	  ur(i,0)=rho
	  ur(i,1)=rho*uu
      Ur(i,2)=(p-rho0*a0**2*f_eta(rho))/(rho0*gamma0)*rho+0.5d0*rho*uu*uu
	 ur(i,3)= uro(i,3)
	 enddo

!call output1(ul)	
!pause
	call  HLLC_EP(nv,jx,u4(:,1)/u4(:,0),ul,ur,h,u_half)

case(2)

	do i=0,3
	call  WENO5_new(nv,jx,u4(:,i),ul(:,i),ur(:,i))
	!call  WENO3_new(nv,jx,u4(:,i),ul(:,i),ur(:,i))
!	call  WENO3_new_change(nv,jx,x,u4(:,i),ul(:,i),ur(:,i))
!	call  WENO3LIU_new(nv,jx,u4(:,i),ul(:,i),ur(:,i))
	!call upwind3(nv,jx,u4(:,i),ul(:,i),ur(:,i))
	enddo
!
!call output1(ul)	
	  !read(*,*)i
	call  HLLC_EP(nv,jx,u4(:,1)/u4(:,0),ul,ur,h,u_half)
case (3)
		call  LF_splitting(u4,ul,ur)
		do i =0,3
			call WENO5_new_LF(nv,jx,ul(:,i),ur(:,i),h(:,i))
		enddo
		u_half(:) = u4(:,1)/u4(:,0)
	case(4)
		do i =-nv,jx+nv
		call eigen_var(u4(i,:),AR(i,:,:),AL(i,:,:))
!		write(*,*) AR(i,2,3)
!		pause
		enddo
!	do i=-nv,jx+nv
!		uo(i,0:3)=matmul(u4(i,0:3),transpose(AL(i,0:3,0:3)))
!	enddo
	uo=0
	do i=0,3
		do j=0,3
	uo(:,i)=uo(:,i)+AL(:,i,j)*u4(:,j)
	enddo
	enddo
!	 write(*,*) u4(0,0:3)
!
!u4=0
!	do i=0,3
!		do j=0,3
!	u4(:,i)=u4(:,i)+Ar(:,j,i)*uo(:,j)
!	enddo
!	enddo
!	 write(*,*) u4(0,0:3)
!	 pause

	do i=0,3
	!call  WENO5_new(nv,jx,uo(:,i),ulo(:,i),uro(:,i))
	call  upwind(nv,jx,uo(:,i),ulo(:,i),uro(:,i))
!	call  WENO3_new(nv,jx,uo(:,i),ulo(:,i),uro(:,i))
	enddo
!	do i=-nv,jx+nv
!	ul(i,0:3)=matmul(ulo(i,0:3),transpose(AR(i,0:3,0:3)))
!	ur(i,0:3)=matmul(uro(i,0:3),transpose(AR(i,0:3,0:3)))
!	enddo
	ul=0
	ur=0
	do i=0,3
		do j=0,3
	ul(:,i)=ul(:,i)+AR(:,i,j)*ulo(:,j)
	ur(:,i)=ur(:,i)+AR(:,i,j)*uro(:,j)
	enddo
	enddo
call output1(ul)	
!pause
!	  read(*,*)i
	!call  HLLC_EP(nv,jx,u4(:,1)/u4(:,0),ul,ur,h,u_half)
	call  HLLC_EP(nv,jx,u4(:,1)/u4(:,0),ul,ur,h,u_half)
!	u_half(:)=u4(:,1)/u4(:,0)
	case(5)

	do i=0,3
	!call  WENO5_new(nv,jx,uo(:,i),ulo(:,i),uro(:,i))
	call  upwind(nv,jx,u(:,i),ul(:,i),ur(:,i))
!	call  WENO3_new(nv,jx,uo(:,i),ulo(:,i),uro(:,i))
	enddo
	do i = -nv,jx+nv
		call eigen_var(ul(i,:),AR(i,:,:),AL(i,:,:))
		ul(i,0:3)=matmul(ul(i,0:3),transpose(AL(i,0:3,0:3)))
	!	ul(i,0:3)=matmul(ul(i,0:3),transpose(AR(i,0:3,0:3)))
		call eigen_var(ur(i,:),AR(i,:,:),AL(i,:,:))
		ur(i,0:3)=matmul(ur(i,0:3),transpose(AL(i,0:3,0:3)))
	!	ur(i,0:3)=matmul(ur(i,0:3),transpose(AR(i,0:3,0:3)))
		call eigen_var(u4(i,:),AR(i,:,:),AL(i,:,:))
		ul(i,0:3)=matmul(ul(i,0:3),transpose(AR(i,0:3,0:3)))
		ur(i,0:3)=matmul(ur(i,0:3),transpose(AR(i,0:3,0:3)))
	enddo

!call output1(ul)	
!pause
!	  read(*,*)i
	call  HLLC_EP(nv,jx,u4(:,1)/u4(:,0),ul,ur,h,u_half)
!	u_half(:)=u4(:,1)/u4(:,0)
endselect


    do k=0,3
     do i=-nv+2,jx+nv-2
        f4(i,k)=h(i,k)-h(i-1,k)
     enddo
    enddo

    end subroutine
 
 
