 subroutine ADER(t,dt)
	  use global
	  use global_cont
	  implicit none
	  integer i,j,k
	  double precision dx1(-nv:jx+nv)
	  double precision uL(-nv:jx+nv,0:3)
	  double precision uR(-nv:jx+nv,0:3)
	  double precision f(-nv:jx+nv,0:3)
	  double precision udx(-nv:jx+nv,0:3)

	  double precision u1(-nv:jx+nv,0:3)
	  double precision uug(-nv:jx+nv)
	  double precision puR_px(-nv:jx+nv,0:3,3)
	  double precision puL_px(-nv:jx+nv,0:3,3)
	  double precision pu_px(-nv:jx+nv,0:3,3)
	  double precision du_dt(-nv:jx+nv,0:3,3)
	  double precision ux(-nv:jx+nv,0:3)
	  double precision ulx(-nv:jx+nv,0:3)
	  double precision ulo(-nv:jx+nv,0:3)
	  double precision uro(-nv:jx+nv,0:3)
	  double precision urx(-nv:jx+nv,0:3)
	  double precision uo(-nv:jx+nv,0:3)
	  double precision AL(-nv:jx+nv,0:3,0:3)
	  double precision AR(-nv:jx+nv,0:3,0:3)
	  double precision src(-nv:jx+nv,0:3)
	  double precision ue(0:3)
	  double precision dt,fgamma,t

    call bound(u)
	do i=-nv+1,nv+jx
		dx1(i)= x(i)-x(i-1)
    	udx(i,0:2)=u(i,0:2)*dx1(i)
	enddo

		do i=-nv,jx+nv
			call trans_u_to_ue(u(i,:),uo(i,:))
			call eigen_var_OR(uo(i,:),Ar(i,:,:))
			call eigen_var_OL(Ar(i,:,:),AL(i,:,:))
			ux(i,0:3) = matmul(AL(i,0:3,0:3),uo(i,0:3))
		enddo



	do i = 0,3
!	call subcell_WENO3(nv,jx,dx1(2),u(:,i),uL(:,i),uR(:,i),pUL_px(:,i,:),pUR_px(:,i,:))

	!	call  WENO3_new(nv,jx,U(:,i), uL(:,i),uR(:,i))
call  WENO3_new(nv,jx,ux(:,i),ulx(:,i),urx(:,i))
!	call  upwind(nv,jx,U(:,i), uL(:,i),uR(:,i))
	enddo

	do i=-nv,jx+nv
		ulo(i,:)= matmul(AR(i,:,:),ulx(i,:))
		uro(i,:)= matmul(AR(i,:,:),urx(i,:))
	enddo
   
	do i=-nv,jx+nv
		call trans_ue_to_u(ulo(i,:),ul(i,:))
		call trans_ue_to_u(uro(i,:),ur(i,:))
	enddo

!src = 0

call Riemann_solver(nv,jx,U(:,1)/U(:,0),uL,uR,pUL_px,pUR_px,U1,pu_px(:,:,1:2))

!call  HLLC_EP(nv,jx,u(:,1)/u(:,0),ul,ur,F,uug)

call material_derivative(U1,pU_px,dU_dt)  

!F=F*dt
!uug=uug*dt
call Gauss(U1,t,dt,dU_dt,F,uug,src)
!	do i = -nv,jx+nv
!	call trans_u_to_ue(U1(i,:),ue(:))
!	call trans_ue_to_FLagrangian(ue(:),F(i,:))
!	F=F*dt
!	enddo
!	uug(:)= U1(:,1)/u1(:,0)*dt

!uug=uug*dt

	!call output1(F)

	x=x+ uug 

	do i=-nv+1,jx+nv-1
		U(i,3) = U(i,3) - (F(i,3)-F(i-1,3))/dx1(i) + src(i,3)
		udx(i,0:2) = udx(i,0:2) - F(i,0:2)+F(i-1,0:2)+src(i,0:2)*dx1(i)
		u(i,3) = fgamma(u(i,3))
	enddo

	do i =-nv,jx+nv
		dx1(i)=x(i)- x(i-1)
	enddo

	do i=-nv+1,nv+jx
		U(i,0:2) = Udx(i,0:2)/dx1(i)
	enddo
	call bound(u)
	!call output

	!pause
 end








