subroutine space(u4,f4,u_half)
use global
use global_cont
implicit none

integer i,j,k
double precision u4(-nv:jx+nv,0:3)
double precision ul(-nv:jx+nv,0:3)
double precision ulx(-nv:jx+nv,0:3)
double precision ulo(-nv:jx+nv,0:3)
double precision uo(-nv:jx+nv,0:3)
double precision ux(-nv:jx+nv,0:3)
double precision ur(-nv:jx+nv,0:3)
double precision uro(-nv:jx+nv,0:3)
double precision urx(-nv:jx+nv,0:3)
double precision h(-nv:jx+nv,0:3)
double precision AL(-nv:jx+nv,0:3,0:3)
double precision AR(-nv:jx+nv,0:3,0:3)
double precision f4(-nv:jx+nv,0:3)
double precision u_half(-nv:jx+nv)
double precision  rho,uu,p,f_eta
integer kind_split

kind_split =2

select case(kind_split) 
case(1)
		do i=-nv,jx+nv
			call trans_u_to_ue(u4(i,:),uo(i,:))
			call eigen_var_OR(uo(i,:),Ar(i,:,:))
			call eigen_var_OL(Ar(i,:,:),AL(i,:,:))
			ux(i,0:3) = matmul(AL(i,0:3,0:3),uo(i,0:3))
		enddo
 		 
!	ux=0
!	do i=0,3
!		do j=0,3
!	ux(:,i)=ux(:,i)+AL(:,i,j)*uo(:,j)
!	enddo
!	enddo
!	ux=uo
   
	do i=0,3
call  WENO3_new(nv,jx,ux(:,i),ulx(:,i),urx(:,i))
	!call  WENO3_origin(nv,jx,ux(:,i),ulx(:,i),urx(:,i))
!	call  WENO3_new_change(nv,jx,x,ux(:,i),ulx(:,i),urx(:,i))
!	call  WENO5_new(nv,jx,ux(:,i),ulx(:,i),urx(:,i))
	!call upwind(nv,jx,ux(:,i),ulx(:,i),urx(:,i))
!	call upwind(nv,jx,ux(:,i),ulx(:,i),urx(:,i))
	enddo
	
ulo=0
uro=0
	do i=0,3
		do j=0,3
	ulo(:,i)=ulo(:,i)+Ar(:,i,j)*ulx(:,j)
	uro(:,i)=uro(:,i)+Ar(:,i,j)*urx(:,j)
	enddo
	enddo
   
	do i=-nv,jx+nv
		call trans_ue_to_u(ulo(i,:),ul(i,:))
		call trans_ue_to_u(uro(i,:),ur(i,:))
	enddo
!call output1(ul)	
!pause
!call  HLLC_EP(nv,jx,u4(:,1)/u4(:,0),ul,ur,h,u_half)
call  HLLC_EP_new(nv,jx,u4(:,1)/u4(:,0),ul,ur,h,u_half)
case(2)

	do i=0,3
	!call  WENO5_new(nv,jx,u4(:,i),ul(:,i),ur(:,i))

	!call  WENO3_new_change(nv,jx,x,u4(:,i),ul(:,i),ur(:,i))
	!call  WENO3_new(nv,jx,u4(:,i),ul(:,i),ur(:,i))
	!call  WENO3_new_change(nv,jx,x,u4(:,i),ul(:,i),ur(:,i))
!	call  WENO3LIU_new(nv,jx,u4(:,i),ul(:,i),ur(:,i))
	call upwind(nv,jx,u4(:,i),ul(:,i),ur(:,i))
	enddo
!
!call output1(ul)	
!pause
	  !read(*,*)i
	!call  HLLC_EPM(nv,jx,u4(:,1)/u4(:,0),ul,ur,h,u_half)
	call  HLLC_EP_new(nv,jx,u4(:,1)/u4(:,0),ul,ur,h,u_half)
!	call  HLLC_EP(nv,jx,u4(:,1)/u4(:,0),ul,ur,h,u_half)
case (3)
		call  LF_splitting(u4,ul,ur)
		do i =0,3
			call WENO5_new_LF(nv,jx,ul(:,i),ur(:,i),h(:,i))
		enddo
		u_half(:) = u4(:,1)/u4(:,0)

	case(4)
		do i =-nv,jx+nv
		call eigen_var_R(u4(i,:),AR(i,:,:))
		call eigen_var_L(AR(i,:,:),AL(i,:,:))
		enddo
		 
		uo=0
	do i=0,3
		do j=0,3
	uo(:,i)=uo(:,i)+AL(:,i,j)*u4(:,j)
	enddo
	enddo

	do i=0,3
	!call  WENO5_new(nv,jx,uo(:,i),ulo(:,i),uro(:,i))
!	call  upwind(nv,jx,uo(:,i),ulo(:,i),uro(:,i))
	call  WENO3_new(nv,jx,uo(:,i),ulo(:,i),uro(:,i))
	enddo

	ul=0
	ur=0
	do i=0,3
		do j=0,3
	ul(:,i)=ul(:,i)+AR(:,i,j)*ulo(:,j)
	ur(:,i)=ur(:,i)+AR(:,i,j)*uro(:,j)
	enddo
	enddo

!call output1(ul)	
!pause
!	  read(*,*)i
call  HLLC_EP(nv,jx,u4(:,1)/u4(:,0),ul,ur,h,u_half)
!	call  HLLC_EPM(nv,jx,u4(:,1)/u4(:,0),ul,ur,h,u_half)

!	call  HLLC_EP_new(nv,jx,u4(:,1)/u4(:,0),ul,ur,h,u_half)
endselect


    do k=0,3
     do i=-nv+2,jx+nv
        f4(i,k)=h(i,k)-h(i-1,k)
     enddo
    enddo

    end subroutine
 
 
