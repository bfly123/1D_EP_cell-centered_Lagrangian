subroutine Gauss(U1,t,dt,dU_dt,F,uug,src)
	  use global_cont
	  use global
	  implicit none
	  double precision U1(-nv:jx+nv,0:3)
	  double precision Ue(0:3)
	  double precision U2(-nv:jx+nv,0:3)
	  double precision du_dt(-nv:jx+nv,0:3,3)
	  double precision F(-nv:jx+nv,0:3)
	  double precision F1(-nv:jx+nv,0:3)
	  double precision F2(-nv:jx+nv,0:3)
	  double precision src(-nv:jx+nv,0:3)
	  double precision src1(-nv:jx+nv,0:3)
	  double precision src2(-nv:jx+nv,0:3)
	  double precision uug(-nv:jx+nv)
	  double precision uu1(-nv:jx+nv)
	  double precision uu2(-nv:jx+nv)
	  double precision dt,t1,t2,t
	  integer i

	  t1=dt*(3.d0-sqrt(3.0))/6

	  call source1(t+t1,src1)

	  U2 = U1    !+ dU_dt(:,:,1)*t1   ! + dU_dt(:,:,2)*t1**2/2  
	  uu1(:)= U2(:,1)/U2(:,0)

	  do i=-nv,jx+nv
		  call trans_U_to_ue(u2(i,0:3),ue(0:3))
		  call trans_ue_to_FLagrangian(ue(0:3),F1(i,0:3))
	  enddo

	  t2=dt*(3.d0+sqrt(3.0))/6

	  call source1(t+t2,src2)

	  U2 = U1   ! + dU_dt(:,:,1)*t2  ! + dU_dt(:,:,2)*t2**2/2  
	  uu2(:)= U2(:,1)/U2(:,0)

	  do i=-nv,jx+nv
		call trans_U_to_ue(u2(i,0:3),ue(0:3))
		call trans_ue_to_FLagrangian(ue(0:3),F2(i,0:3))
	  enddo

	  F= dt*(F1+F2)/2
	  uug =dt*(uu1+uu2)/2
	  src=dt*(src1+src2)/2
	!call output1(F2)
 end









