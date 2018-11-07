 subroutine ADER(dt)
	  use global
	  use global_cont
	  implicit none
	  integer i,j,k
	  double precision dx1(-nv:jx+nv)
	  double precision uL(-nv:jx+nv,0:3)
	  double precision uR(-nv:jx+nv,0:3)

	  double precision u1(-nv:jx+nv,0:3)
	  double precision puR_px(-nv:jx+nv,0:3)
	  double precision puL_px(-nv:jx+nv,0:3)
	  double precision pu_px(-nv:jx+nv,0:3)
	  double precision du_dt(-nv:jx+nv,0:3,3)

   call bound(u)
	do i=-nv+1,nv+jx
		dx1(i)= x(i)-x(i-1)
    	udx(i,0:2)=u(i,0:2)*dx1(i)
	enddo

	do i = 0,3
		call subcell_WENO3(nv,jx,dx1(1),u(:,i),uL(:,i),uR(:,i),pUL_px(:,i,:),pUR_px(:,i,:))
	enddo

	call Reimann_solver(uL,uR,pUL_px,pUR_px,U1,pu_px(:,:,1:2))

	call material_derivative(U1,pU_px,dU_dt)  


	call trans_U_to_ue(U1,ue)
	call trans_ue_to_FLagrangian(U1,F)

	call Gauss(U1,dU_dt,F,uug)

	x=x+ uug 

	U(:,3) = U(:,3) + f(:,3)/dx1(:)

	udx(:,0:2) = udx(:,0:2)+f(:,0:2)

	do i =-nv,jx+nv
	dx1(i)=x(i)- x(i-1)
	enddo
	U(:,0:2) = Udx(:,0:2)/dx1(:)
 end








