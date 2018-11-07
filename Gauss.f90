subroutine Gauss(U1,dU_dt,F,uug)
	  use global_cont
	  implicit none
	  double precision U1(-nv:jx+nv,0:3)
	  double precision du_dt(-nv:jx+nv,0:3,3)
	  double precision F(-nv:jx+nv,0:3)
	  double precision uug(-nv:jx+nv)
	  double precision uu1(-nv:jx+nv)
	  double precision uu2(-nv:jx+nv)
	  double precision dt,t1,t2

	  t1=dt*(3.d0-sqrt(3))/6

	  U2 = U1+ dU_dt(:,:,1)*t1+ dU_dt(:,:,2)*t1**2/2  
	  uu1(:)= U2(:,1)/U2(:,0)

	  call trans_U_to_ue(u2,ue)
	  call trans_ue_to_FLagrangian(ue,F1)

	  t1=dt*(3.d0-sqrt(3))/6

	  U2 = U1+ dU_dt(:,:,1)*t2+ dU_dt(:,:,2)*t2**2/2  
	  uu2(:)= U2(:,1)/U2(:,0)

	  call trans_U_to_ue(u2,ue)
	  call trans_ue_to_FLagrangian(ue,F2)

	  F= dt*(F1+F2)/2
	  uug =dt*(uu1+uu2)/2
 end









