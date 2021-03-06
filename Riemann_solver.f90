!subroutine  Riemann_solver(nv,jx,U,uL,uR,pUL_px,pUR_px,U1,pu_px)
!!subroutine  Riemann_solver(nv,jx,U,uL,uR,h,u_half)
!	!  use global
!	  use global_cont
!	  implicit none
!	  double precision  tmp,tmp1
!	  integer nv,i,jx,j
!	  double precision u(-nv:jx+nv)
!	  double precision u1(-nv:jx+nv,0:3)
!	  double precision ul(-nv:jx+nv,0:3)
!	  double precision ur(-nv:jx+nv,0:3)
!	  double precision h(-nv:jx+nv,0:3)
!	  double precision u_half(-nv:jx+nv)
!	  double precision ue(0:3)
!	  double precision A(0:3,0:3)
!
!	  double precision rhoL,rhoL_star,uuL,sxxL,pL,cL,sL,sxxL_star,pL_star,sigma_star
!	  double precision rhoR,rhoR_star,uuR,sxxR,pR,cR,sR,sxxR_star,pR_star,eR,eL
!	  double precision cap_rhoR,cap_uuR,cap_sxxR,cap_pR,cap_sigmaR
!	  double precision cap_rhoL,cap_uuL,cap_sxxL,cap_pL,cap_sigmaL
!	  double precision s_star,sigmaL,sigmaR,c1,c2,f_eta
!	  double precision pu_px(-nv:jx+nv,0:3,2)
!	  double precision puL_px(-nv:jx+nv,0:3,2)
!	  double precision puR_px(-nv:jx+nv,0:3,2)
!
!do i =-nv, jx+nv-1
!!	 call trans_u_to_ue(ul(i,:),ue(:))
!	  rhoL= uL(i,0)
!	  uuL = uL(i,1)
!	  pL =  uL(i,2)
!	  sxxL= uL(i,3)
!	  sigmaL=-pL+sxxL
!	 !call state_p_to_ei(rhoL,pL,eL)
!	 call sound(uL(i.:),cL)
!
!	 !call trans_u_to_ue(uR(i,:),ue(:))
!	  rhoR= uR(i,0)
!	  uuR = uR(i,1)
!	  pR =  uR(i,2)
!	  sxxR= uR(i,3)
!	  sigmaR=-pR+sxxR
!	  !call state_p_to_ei(rhoR,pR,eR)
!	  call sound(uR(i,:),cR)
!	  
!	  sL=min(uuL-cL,uuR-cR)
!	  sR=max(uuL+cL,uuR+cR)
!     s_star = (sigmaL-sigmaR+rhoL*uuL*(sL-uuL)-rhoR*uuR*(sR-uuR))/(rhoL*(sL-uuL)-rhoR*(sR-uuR))
!	 rhoL_star=rhoL*(uuL-sL)/(s_star-sL)
!	 rhoR_star=rhoR*(uuR-sR)/(s_star-sR)
!
!! if (abs(tmp).gt.2.d0/3*Y0)then
!!		write(*,*)i, abs(tmp)-2.d0/3*Y0
!!	pause
!!	endif
!
!!	 if(abs(tmp).ge.2.d0/3*Y0)then
!!		 write(*,*)i, tmp,"**********"
!!		 pause
!!	endif
!
!	 tmp =-4.d0/3*miu*log(rhoL_star/rhoL)+sxxL
!
!	 !if (dabs(sxxL).lt.2.d0/3*Y0.and.abs(tmp).ge.2.d0/3*Y0) then
!	 if (abs(tmp).ge.2.d0/3*Y0) then
!!		 write(*,*) i,"left"
!!		 write(*,*) rhoL,uuL,pL,sxxL
!!		 write(*,*) i,"star"
!!		 write(*,*) rhoL_star,rhoR_star,s_star
!
!		if(tmp.ge.2.d0/3*Y0)then
!		cap_sxxL=2.d0/3*Y0
!		cap_rhoL=rhoL*exp(-Y0/miu/2+3.d0/4*sxxL/miu)
!		else if(tmp.le.-2.d0/3*Y0)then
!			cap_sxxL=-2.d0/3*Y0
!			cap_rhoL= rhoL*exp(Y0/miu/2+3.d0/4*sxxL/miu)
!		endif
!		!tmp1= rhoL*cap_rhoL/(cap_rhoL-rhoL) 
!		tmp1=1/rhoL-1/cap_rhoL
!		c1=1/rho0/gamma0
!		c2=a0**2/gamma0
!		cap_pL=(2*(c2*f_eta(cap_rhoL)+eL)-tmp1*(sigmaL+cap_sxxL))/(2*c1-tmp1)
!		cap_sigmaL=-cap_pL+cap_sxxL
!		if(rhoL_star.ge.rhoL)then
!			 cap_uuL=uuL-sqrt((sigmaL-cap_sigmaL)*tmp1)
!		else
!			 cap_uuL=uuL+sqrt((sigmaL-cap_sigmaL)*tmp1)
!		endif
!		rhoL=cap_rhoL
!		uuL=cap_uuL
!		pL=cap_pL
!		sxxL=cap_sxxL
!		sigmaL=cap_sigmaL
!		ue(0)=rhoL
!		ue(1)=uuL
!		ue(2)=pL
!		ue(3)=sxxL
!	call sound(ue,cL)
!	!write(*,*)ue(:)
!!	pause
!	endif
!	
!
!!	 if (abs(tmp).gt.2.d0/3*Y0)then
!!		write(*,*)i, abs(tmp)-2.d0/3*Y0
!!	pause
!!	endif
!
!
!	 tmp =-4.d0/3*miu*log(rhoR_star/rhoR)+sxxR
!	 !if (dabs(sxxR).lt.2.d0/3*Y0.and.abs(tmp).ge.2.d0/3*Y0) then
!	 if (abs(tmp).ge.2.d0/3*Y0) then
!	!	 write(*,*) i,"******"
!	!	 write(*,*) rhoR,uuR,pR,sxxR
!	!	 pause
!		 if(tmp.ge.2.d0/3*Y0)then
!			cap_sxxR= 2.d0/3*Y0
!			cap_rhoR= rhoR*exp(-Y0/miu/2+3.d0/4*sxxR/miu)
!		else if(tmp.le.-2.d0/3*Y0)then
!			cap_sxxR=-2.d0/3*Y0
!			cap_rhoR= rhoR*exp(Y0/miu/2+3.d0/4*sxxR/miu)
!		endif
!		!tmp1= rhoR*cap_rhoR/(cap_rhoR-rhoR) 
!		tmp1= 1/rhoR-1/cap_rhoR 
!		c1=1/rho0/gamma0
!		c2=a0**2/gamma0
!		cap_pR=(2*(c2*f_eta(cap_rhoR)+eR)-tmp1*(sigmaR+cap_sxxR))/(2*c1-tmp1)
!		cap_sigmaR=-cap_pR+cap_sxxR
!		if(rhoR_star.ge.rhoR)then
!			cap_uuR=uuR+sqrt((sigmaR-cap_sigmaR)*tmp1)
!		else
!			cap_uuR=uuR-sqrt((sigmaR-cap_sigmaR)*tmp1)
!		endif
!!		 write(*,*)i, cap_pR,pR 
!!		 pause
!		rhoR=cap_rhoR
!		uuR=cap_uuR
!		pR=cap_pR
!		sxxR=cap_sxxR
!		sigmaR=cap_sigmaR
!		ue(0)=rhoR
!		ue(1)=uuR
!		ue(2)=pR
!		ue(3)=sxxR
!		call sound(ue,cR)
!!	write(*,*)ue(:)
!		 !pause
!!		 j=i
!!	pause
!endif
!
!	
!	sL=min(uuL-cL,uuR-cR)
!	sR=max(uuL+cL,uuR+cR)
!
!    s_star = (sigmaL-sigmaR+rhoL*uuL*(sL-uuL)-rhoR*uuR*(sR-uuR))/(rhoL*(sL-uuL)-rhoR*(sR-uuR))
!
!	rhoL_star=rhoL*(uuL-sL)/(s_star-sL)
!	rhoR_star=rhoR*(uuR-sR)/(s_star-sR)
!
!
!	tmp=-4.d0/3*miu*log(rhoL_star/rhoL)+sxxL
! if (tmp.ge.2.d0/3*Y0) then
!	 sxxL_star=2.d0/3*Y0
! else if(tmp.le.-2.d0/3*Y0)then
!	 sxxL_star=-2.d0/3*Y0
! else
!	 sxxL_star=tmp
!endif
!
!	tmp= -4.d0/3*miu*log(rhoR_star/rhoR)+sxxR
! if (tmp.ge.2.d0/3*Y0) then
!	 sxxR_star=2.d0/3*Y0
! else if(tmp.le.-2.d0/3*Y0)then
!	 sxxR_star=-2.d0/3*Y0
! else
!	 sxxR_star=tmp
!endif
!
!	sigma_star = sigmaL-rhoL*(sL-uuL)*(s_star-uuL)
!	pL_star=sxxL_star-sigma_star
!	pR_star=sxxR_star-sigma_star
!!	if(i==j)then
!!		write(*,*)j,pL_star,pR_star,sxxL_star,sxxR_star,s_star
!!		pause
!!	endif
!
!
!!	if (S_star.ge.u(i))then
!!		h(i,0)=0
!!		h(i,1)=pL_star-sxxL_star
!!		h(i,2)=(pL_star-sxxL_star)*s_star
!!		h(i,3)=-4.d0*miu/3*s_star
!!		u_half(i)=s_star
!!	else 
!!		h(i,0)=0
!!		h(i,1)=pR_star-sxxR_star
!!		h(i,2)=(pR_star-sxxR_star)*s_star
!!		h(i,3)=-4.d0*miu/3*s_star
!!		u_half(i)=s_star
!!	endif
!
! if (S_star.ge.u(i))then
!	    ue(0)= rhoL_star
!		ue(1)= s_star
!		ue(2)= pL_star
!	    ue(3)= sxxL_star
!		else 
!	
!	    ue(0)= rhoR_star
!		ue(1)= s_star
!		ue(2)= pR_star
!	    ue(3)= sxxR_star
!	endif
! call trans_ue_to_u(ue,u1(i,:))
!! call trans_ue_to_FLagrangian(ue,h(i,:))
!!	u_half(i) = s_star
!
!
! !*************solve pu_px with a HLL****************
! call trans_ue_Ae(ue, A)
!
! do j=1,2
!	pu_px(i,:,j) = (sR*puR_px(i,:,j) -sL*puL_px(i,:,j)+ &
!				 matmul(A,puL_px(i,:,j))- matmul(A,puR_px(i,:,j)))/(sR-sL)
! enddo
!
!enddo
!
!end

subroutine  HLLC_EP_NEW_origin(nv,jx,u,uL,uR,u1)
!subroutine  Riemann_solver(nv,jx,U,uL,uR,h,u_half)
	!  use global
	  use global_cont
	  implicit none
	  double precision  tmp,tmp1
	  integer nv,i,jx,j
	  double precision u(-nv:jx+nv)
	  double precision u1(-nv:jx+nv,0:3)
	  double precision ul(-nv:jx+nv,0:3)
	  double precision ur(-nv:jx+nv,0:3)
	  double precision h(-nv:jx+nv,0:3)
	  double precision u_half(-nv:jx+nv)
	  double precision ue(0:3)
	  double precision A(0:3,0:3)

	  double precision rhoL,rhoL_star,uuL,sxxL,pL,cL,sL,sxxL_star,pL_star,sigma_star
	  double precision rhoR,rhoR_star,uuR,sxxR,pR,cR,sR,sxxR_star,pR_star,eR,eL
	  double precision cap_rhoR,cap_uuR,cap_sxxR,cap_pR,cap_sigmaR
	  double precision cap_rhoL,cap_uuL,cap_sxxL,cap_pL,cap_sigmaL
	  double precision s_star,sigmaL,sigmaR,c1,c2,f_eta
	  double precision pu_px(-nv:jx+nv,0:3,2)
	  double precision puL_px(-nv:jx+nv,0:3,2)
	  double precision puR_px(-nv:jx+nv,0:3,2)

do i =-nv, jx+nv-1
!	 call trans_u_to_ue(ul(i,:),ue(:))
	  rhoL= uL(i,0)
	  uuL = uL(i,1)
	  pL =  uL(i,2)
	  sxxL= uL(i,3)
	  sigmaL=-pL+sxxL
	 call state_p_to_ei(rhoL,pL,eL)
	 call sound(uL(i,:),cL)

	 !call trans_u_to_ue(uR(i,:),ue(:))
	  rhoR= uR(i,0)
	  uuR = uR(i,1)
	  pR =  uR(i,2)
	  sxxR= uR(i,3)
	  sigmaR=-pR+sxxR
	  call state_p_to_ei(rhoR,pR,eR)
	  call sound(uR(i,:),cR)
	  
	  sL=min(uuL-cL,uuR-cR)
	  sR=max(uuL+cL,uuR+cR)
     s_star = (sigmaL-sigmaR+rhoL*uuL*(sL-uuL)-rhoR*uuR*(sR-uuR))/(rhoL*(sL-uuL)-rhoR*(sR-uuR))
	 rhoL_star=rhoL*(uuL-sL)/(s_star-sL)
	 rhoR_star=rhoR*(uuR-sR)/(s_star-sR)

! if (abs(tmp).gt.2.d0/3*Y0)then
!		write(*,*)i, abs(tmp)-2.d0/3*Y0
!	pause
!	endif

!	 if(abs(tmp).ge.2.d0/3*Y0)then
!		 write(*,*)i, tmp,"**********"
!		 pause
!	endif

	 tmp =-4.d0/3*miu*log(rhoL_star/rhoL)+sxxL

	 if (dabs(sxxL).lt.2.d0/3*Y0.and.abs(tmp).ge.2.d0/3*Y0) then
	! if (abs(tmp).ge.2.d0/3*Y0) then
!		 write(*,*) i,"left"
!		 write(*,*) rhoL,uuL,pL,sxxL
!		 write(*,*) i,"star"
!		 write(*,*) rhoL_star,rhoR_star,s_star

		if(tmp.ge.2.d0/3*Y0)then
		cap_sxxL=2.d0/3*Y0
		cap_rhoL=rhoL*exp(-Y0/miu/2+3.d0/4*sxxL/miu)
		else if(tmp.le.-2.d0/3*Y0)then
			cap_sxxL=-2.d0/3*Y0
			cap_rhoL= rhoL*exp(Y0/miu/2+3.d0/4*sxxL/miu)
		endif
		!tmp1= rhoL*cap_rhoL/(cap_rhoL-rhoL) 
		tmp1=1/rhoL-1/cap_rhoL
		c1=1/rho0/gamma0
		c2=a0**2/gamma0
		cap_pL=(2*(c2*f_eta(cap_rhoL)+eL)-tmp1*(sigmaL+cap_sxxL))/(2*c1-tmp1)
		cap_sigmaL=-cap_pL+cap_sxxL
		if(rhoL_star.ge.rhoL)then
			 cap_uuL=uuL-sqrt((sigmaL-cap_sigmaL)*tmp1)
		else
			 cap_uuL=uuL+sqrt((sigmaL-cap_sigmaL)*tmp1)
		endif
		rhoL=cap_rhoL
		uuL=cap_uuL
		pL=cap_pL
		sxxL=cap_sxxL
		sigmaL=cap_sigmaL
		ue(0)=rhoL
		ue(1)=uuL
		ue(2)=pL
		ue(3)=sxxL
	call sound(ue,cL)
	!write(*,*)ue(:)
!	pause
	endif
	

!	 if (abs(tmp).gt.2.d0/3*Y0)then
!		write(*,*)i, abs(tmp)-2.d0/3*Y0
!	pause
!	endif


	 tmp =-4.d0/3*miu*log(rhoR_star/rhoR)+sxxR
	if (dabs(sxxR).lt.2.d0/3*Y0.and.abs(tmp).ge.2.d0/3*Y0) then
	! if (abs(tmp).ge.2.d0/3*Y0) then
	!	 write(*,*) i,"******"
	!	 write(*,*) rhoR,uuR,pR,sxxR
	!	 pause
		 if(tmp.ge.2.d0/3*Y0)then
			cap_sxxR= 2.d0/3*Y0
			cap_rhoR= rhoR*exp(-Y0/miu/2+3.d0/4*sxxR/miu)
		else if(tmp.le.-2.d0/3*Y0)then
			cap_sxxR=-2.d0/3*Y0
			cap_rhoR= rhoR*exp(Y0/miu/2+3.d0/4*sxxR/miu)
		endif
		!tmp1= rhoR*cap_rhoR/(cap_rhoR-rhoR) 
		tmp1= 1/rhoR-1/cap_rhoR 
		c1=1/rho0/gamma0
		c2=a0**2/gamma0
		cap_pR=(2*(c2*f_eta(cap_rhoR)+eR)-tmp1*(sigmaR+cap_sxxR))/(2*c1-tmp1)
		cap_sigmaR=-cap_pR+cap_sxxR
		if(rhoR_star.ge.rhoR)then
			cap_uuR=uuR+sqrt((sigmaR-cap_sigmaR)*tmp1)
		else
			cap_uuR=uuR-sqrt((sigmaR-cap_sigmaR)*tmp1)
		endif
!		 write(*,*)i, cap_pR,pR 
!		 pause
		rhoR=cap_rhoR
		uuR=cap_uuR
		pR=cap_pR
		sxxR=cap_sxxR
		sigmaR=cap_sigmaR
		ue(0)=rhoR
		ue(1)=uuR
		ue(2)=pR
		ue(3)=sxxR
		call sound(ue,cR)
!	write(*,*)ue(:)
		 !pause
!		 j=i
!	pause
endif

	
	sL=min(uuL-cL,uuR-cR)
	sR=max(uuL+cL,uuR+cR)

    s_star = (sigmaL-sigmaR+rhoL*uuL*(sL-uuL)-rhoR*uuR*(sR-uuR))/(rhoL*(sL-uuL)-rhoR*(sR-uuR))

	rhoL_star=rhoL*(uuL-sL)/(s_star-sL)
	rhoR_star=rhoR*(uuR-sR)/(s_star-sR)


	tmp=-4.d0/3*miu*log(rhoL_star/rhoL)+sxxL
 if (tmp.ge.2.d0/3*Y0) then
	 sxxL_star=2.d0/3*Y0
 else if(tmp.le.-2.d0/3*Y0)then
	 sxxL_star=-2.d0/3*Y0
 else
	 sxxL_star=tmp
endif

	tmp= -4.d0/3*miu*log(rhoR_star/rhoR)+sxxR
 if (tmp.ge.2.d0/3*Y0) then
	 sxxR_star=2.d0/3*Y0
 else if(tmp.le.-2.d0/3*Y0)then
	 sxxR_star=-2.d0/3*Y0
 else
	 sxxR_star=tmp
endif

	sigma_star = sigmaL-rhoL*(sL-uuL)*(s_star-uuL)
	pL_star=sxxL_star-sigma_star
	pR_star=sxxR_star-sigma_star
!	if(i==j)then
!		write(*,*)j,pL_star,pR_star,sxxL_star,sxxR_star,s_star
!		pause
!	endif


!	if (S_star.ge.u(i))then
!		h(i,0)=0
!		h(i,1)=pL_star-sxxL_star
!		h(i,2)=(pL_star-sxxL_star)*s_star
!		h(i,3)=-4.d0*miu/3*s_star
!		u_half(i)=s_star
!	else 
!		h(i,0)=0
!		h(i,1)=pR_star-sxxR_star
!		h(i,2)=(pR_star-sxxR_star)*s_star
!		h(i,3)=-4.d0*miu/3*s_star
!		u_half(i)=s_star
!	endif

 if (S_star.ge.u(i))then
	    u1(i,0)= rhoL_star
		u1(i,1)= s_star
		u1(i,2)= pL_star
	    u1(i,3)= sxxL_star
		else 
	
	    u1(i,0)= rhoR_star
		u1(i,1)= s_star
		u1(i,2)= pR_star
	    u1(i,3)= sxxR_star
	endif

enddo

end




















