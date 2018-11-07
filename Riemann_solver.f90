subroutine  Reimann_solver(uL,uR,pUL_px,pUR_px,U1,pu_px)
	!  use global
	  use global_cont
	  implicit none
	  integer i,j,k
	  double precision  tmp,nv,jx

	  !double precision u(-nv:jx+nv)
	  double precision ul(-nv:jx+nv,0:3)
	  double precision ur(-nv:jx+nv,0:3)
	  double precision h(-nv:jx+nv,0:3)
	  double precision u1(-nv:jx+nv,0:3)
	  double precision puL_px(-nv:jx+nv,0:3,2)
	  double precision puR_px(-nv:jx+nv,0:3,2)
	  double precision pu_px(-nv:jx+nv,0:3,2)
	 
	  double precision ue(0:3)
	  double precision A(0:3,0:3)
	  double precision rhoL,rhoL_star,uuL,sxxL,pL,cL,sL,sxxL_star,pL_star,sigma_star
	  double precision rhoR,rhoR_star,uuR,sxxR,pR,cR,sR,sxxR_star,pR_star
	  double precision s_star,sigmaL,sigmaR


	  !********************HLLC_EPM**************************
do i =-nv, jx+nv-1
	 call trans_u_to_ue(ul(i,:),ue(:))
	  rhoL= ue(0)
	  uuL = ue(1)
	  pL =  ue(2)
	  sxxL= ue(3)
	  sigmaL=-pL+sxxL
	 call sound(ue,cL)
	 call trans_u_to_ue(uR(i,:),ue(:))
	  rhoR= ue(0)
	  uuR = ue(1)
	  pR =  ue(2)
	  sxxR= ue(3)
	  sigmaR=-pR+sxxR

	  call sound(ue,cR)
	  
	  sL=min(uuL-cL,uuR-cR)
	  sR=max(uuL+cL,uuR+cR)
     s_star = (sigmaL-sigmaR+rhoL*uuL*(sL-uuL)-rhoR*uuR*(sR-uuR))/(rhoL*(sL-uuL)-rhoR*(sR-uuR))
	 rhoL_star=rhoL*(uuL-sL)/(s_star-sL)
	 rhoR_star=rhoR*(uuR-sR)/(s_star-sR)

	 tmp = 3.d0/4/miu*log(rhoL_star/rhoL)+sxxL

	 if (dabs(sxxL).ge.2.d0/3*Y0) then
		 sxxL_star=sxxL
		else if (dabs(tmp).ge.2.d0/3*Y0) then
		sxxL_star=sxxL
	else
		sxxL_star = tmp
	endif

	 tmp = 3.d0/4/miu*log(rhoR_star/rhoR)+sxxR
 if (dabs(sxxR).ge.2.d0/3*Y0) then
		 sxxR_star=sxxR
		else if (dabs(tmp).ge.2.d0/3*Y0) then
		sxxR_star=sxxR
	else
		sxxR_star = tmp
	endif
	sigma_star = sigmaL-rhoL*(sL-uuL)*(s_star-uuL)
	pL_star=sxxL_star-sigma_star
	pR_star=sxxR_star-sigma_star

 if (S_star.ge.u(i))then
	    ue(0)= rhoL_star
		ue(1)= s_star
		ue(2)= pL_star
	    ue(3)= sxxL_star
		else 
	
	    ue(0)= rhoR_star
		ue(1)= s_star
		ue(2)= pR_star
	    ue(3)= sxxR_star
	endif
 call trans_ue_to_u(ue,u1(i,:))


 !*************solve pu_px with a HLL****************
 call trans_ue_A(U1(i,:), ue, A)

 do j=1,2
	pu_px(i,:,j) = (sR*puR_px(i,:,j) -sL*puL_px(i,:,j)+ &
				 matmul(A,puL_px(i,:,j)) - matmul(A,puR_px(i,:,j)))/(sR-sL)
 enddo

enddo

end











