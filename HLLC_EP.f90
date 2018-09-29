subroutine HLLC_EP(nv,jx,u,ul,ur,h,u_half)
	  use global_cont
	  implicit none
	  double precision  f_eta,fgamma 
	  integer nv,i,jx,j
	  double precision u(-nv:jx+nv)
	  double precision ul(-nv:jx+nv,0:3)
	  double precision ur(-nv:jx+nv,0:3)
	  double precision h(-nv:jx+nv,0:3)
	  double precision u_half(-nv:jx+nv)
	  double precision ue(0:3),FL(0:3),FR(0:3),u_hll(0:3)

	  double precision rhoL,uuL,sxxL,pL,cL,sL,sxxL_bar,sxxL_star,sxL_star,pL_star,sigmaxL_star
	  double precision rhoR,uuR,sxxR,pR,cR,sR,sxxR_bar,sxxR_star,sxR_star,pR_star,sigmaxR_star
	  double precision s_barStar,feta,feta1,feta_eta,a_squre,s_star,f_eta_eta

do i =-nv, jx+nv-1

	 call trans_u_to_ue(ul(i,:),ue(:))
	  rhoL= ue(0)
	  uuL = ue(1)
	  pL =  ue(2)
	  sxxL= ue(3)
	 call sound(ue,cL)

	 call trans_u_to_ue(uR(i,:),ue(:))
	  rhoR= ue(0)
	  uuR = ue(1)
	  pR =  ue(2)
	  sxxR= ue(3)
	  call sound(ue,cR)
	  
	  sL=min(uuL-cL,uuR-cR)
	  sR=max(uuL+cL,uuR+cR)
      s_barStar = (sxxR-sxxL)/(rhoL*(sL-uuL)-rhoR*(sR-uuR))
      sxxL_bar=sxxL+rhoL*(sL-uuL)*s_barStar
	sxxR_bar=sxxR+rhoR*(sR-uuR)*s_barStar
	sxxL_star=Fgamma(sxxL_bar)
	sxxR_star=Fgamma(sxxR_bar)
	
	sxL_star=(sxxL_star-sxxL)/(rhoL*(sL-uuL))
	sxR_star=(sxxR_star-sxxR)/(rhoR*(sR-uuR))
	
	s_star=(pR-pL+rhoL*(uuL-sxL_star)*(sL-uuL)-rhoR*(uuR-sxR_star)*(sR-uuR))/(rhoL*(sL-uuL)-rhoR*(sR-uuR))

	!u_half(i)=s_star
	pL_star=pL+rhoL*(sL-uuL)*(s_star+sxL_star-uuL)
	pR_star=pR+rhoR*(sR-uuR)*(s_star+sxR_star-uuR)
	
!	sigmaxL_star=sxxL_star-PL_star
!	sigmaxR_star=sxxR_star-PR_star
	
!	if (sL.ge.u(i))then
!		h(i,0)=0
!		h(i,1)=pL-sxxL
!		h(i,2)=(pL-sxxL)*uuL
!		h(i,3)=-4*miu/3*uuL
!		u_half(i)=uuL
!	!if (s_star.le.0)then
!	!if (s_star.ge.))then
!else if (S_star.ge.u(i))then
!		h(i,0)=0
!		h(i,1)=pL_star-sxxL_star
!		h(i,2)=(pL_star-sxxL_star)*s_star
!		h(i,3)=-4*miu/3*s_star
!		u_half(i)=s_star
!	else if (SR.ge.u(i))then
!		h(i,0)=0
!		h(i,1)=pR_star-sxxR_star
!		h(i,2)=(pR_star-sxxR_star)*s_star
!		h(i,3)=-4*miu/3*s_star
!		u_half(i)=s_star
!	else
!		h(i,0)=0
!		h(i,1)=pR-sxxR
!		h(i,2)=(pR-sxxR)*uuR
!		h(i,3)=-4*miu/3*uuR
!		u_half(i)=uuR
!	endif
 if (S_star.ge.u(i))then
		h(i,0)=0
		h(i,1)=pL_star-sxxL_star
		h(i,2)=(pL_star-sxxL_star)*s_star
		h(i,3)=-4.d0*miu/3*s_star
		u_half(i)=s_star
	else 
		h(i,0)=0
		h(i,1)=pR_star-sxxR_star
		h(i,2)=(pR_star-sxxR_star)*s_star
		h(i,3)=-4.d0*miu/3*s_star
		u_half(i)=s_star
	endif

! if (S_star.le.u(i))then
!		h(i,0)=0
!		h(i,1)=pL-sxxL
!		h(i,2)=(pL-sxxL)*uuL
!		h(i,3)=-4.d0*miu/3*s_star
!		u_half(i)=s_star
!	else 
!		h(i,0)=0
!		h(i,1)=pR-sxxR
!		h(i,2)=(pR-sxxR)*uuR
!		h(i,3)=-4.d0*miu/3*s_star
!		u_half(i)=s_star
!	endif

enddo

endsubroutine

subroutine HLL_EP(nv,jx,u,ul,ur,h,u_half)
	  use global_cont
	  implicit none
	  double precision  f_eta,fgamma 
	  integer nv,i,jx
	  double precision u(-nv:jx+nv)
	  double precision ul(-nv:jx+nv,0:3)
	  double precision ur(-nv:jx+nv,0:3)
	  double precision h(-nv:jx+nv,0:3)
	  double precision u_half(-nv:jx+nv)
	  double precision ue(0:3),FL(0:3),FR(0:3),u_hll(0:3)

	  double precision rhoL,uuL,sxxL,pL,cL,sL
	  double precision rhoR,uuR,sxxR,pR,cR,sR
	  double precision rho,uu,sxx,p

do i =-nv, jx+nv-1

	 call trans_u_to_ue(ul(i,:),ue(:))
	  rhoL= ue(0)
	  uuL = ue(1)
	  pL =  ue(2)
	  sxxL= ue(3)
	 call sound(ue,cL)
	 call trans_ue_to_Feuler(ue,FL)

	 call trans_u_to_ue(uR(i,:),ue(:))
	  rhoR= ue(0)
	  uuR = ue(1)
	  pR =  ue(2)
	  sxxR= ue(3)
	  call sound(ue,cR)
	 call trans_ue_to_Feuler(ue,FR)

	  sL=min(uuL-cL,uuR-cR)
	  sR=max(uuL+cL,uuR+cR)

	  u_hll(:)=(sr*uR(i,:)-sL*uL(i,:)+FL(:)-FR(:))/(sR-sL)

	  call trans_u_to_ue(u_hll,ue)
	  rho=ue(0)
	  uu=ue(1)
	  p=ue(2)
	  sxx=ue(3)
	h(i,0)=0
	h(i,1)=p-sxx
	h(i,2)=(p-sxx)*uu
	h(i,3)=-4.d0*miu/3*uu
	u_half(i)=uu

! if (S_star.le.u(i))then
!		h(i,0)=0
!		h(i,1)=pL-sxxL
!		h(i,2)=(pL-sxxL)*uuL
!		h(i,3)=-4.d0*miu/3*s_star
!		u_half(i)=s_star
!	else 
!		h(i,0)=0
!		h(i,1)=pR-sxxR
!		h(i,2)=(pR-sxxR)*uuR
!		h(i,3)=-4.d0*miu/3*s_star
!		u_half(i)=s_star
!	endif

enddo

endsubroutine

subroutine HLLC_EPM(nv,jx,u,ul,ur,h,u_half)
	  use global_cont
	  implicit none
	  double precision  tmp
	  integer nv,i,jx
	  double precision u(-nv:jx+nv)
	  double precision ul(-nv:jx+nv,0:3)
	  double precision ur(-nv:jx+nv,0:3)
	  double precision h(-nv:jx+nv,0:3)
	  double precision u_half(-nv:jx+nv)
	  double precision ue(0:3)

	  double precision rhoL,rhoL_star,uuL,sxxL,pL,cL,sL,sxxL_star,pL_star,sigma_star
	  double precision rhoR,rhoR_star,uuR,sxxR,pR,cR,sR,sxxR_star,pR_star
	  double precision s_star,sigmaL,sigmaR

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
		h(i,0)=0
		h(i,1)=pL_star-sxxL_star
		h(i,2)=(pL_star-sxxL_star)*s_star
		h(i,3)=-4.d0*miu/3*s_star
		u_half(i)=s_star
	else 
		h(i,0)=0
		h(i,1)=pR_star-sxxR_star
		h(i,2)=(pR_star-sxxR_star)*s_star
		h(i,3)=-4.d0*miu/3*s_star
		u_half(i)=s_star
	endif

! if (S_star.le.u(i))then
!		h(i,0)=0
!		h(i,1)=pL-sxxL
!		h(i,2)=(pL-sxxL)*uuL
!		h(i,3)=-4.d0*miu/3*s_star
!		u_half(i)=s_star
!	else 
!		h(i,0)=0
!		h(i,1)=pR-sxxR
!		h(i,2)=(pR-sxxR)*uuR
!		h(i,3)=-4.d0*miu/3*s_star
!		u_half(i)=s_star
!	endif

enddo
end

subroutine HLLC_EP_new(nv,jx,u,ul,ur,h,u_half)
	  use global_cont
	  implicit none
	  double precision  tmp,tmp1
	  integer nv,i,jx,j
	  double precision u(-nv:jx+nv)
	  double precision ul(-nv:jx+nv,0:3)
	  double precision ur(-nv:jx+nv,0:3)
	  double precision h(-nv:jx+nv,0:3)
	  double precision u_half(-nv:jx+nv)
	  double precision ue(0:3)

	  double precision rhoL,rhoL_star,uuL,sxxL,pL,cL,sL,sxxL_star,pL_star,sigma_star
	  double precision rhoR,rhoR_star,uuR,sxxR,pR,cR,sR,sxxR_star,pR_star,eR,eL
	  double precision cap_rhoR,cap_uuR,cap_sxxR,cap_pR,cap_sigmaR
	  double precision cap_rhoL,cap_uuL,cap_sxxL,cap_pL,cap_sigmaL
	  double precision s_star,sigmaL,sigmaR,c1,c2,f_eta

do i =-nv, jx+nv-1
	 call trans_u_to_ue(ul(i,:),ue(:))
	  rhoL= ue(0)
	  uuL = ue(1)
	  pL =  ue(2)
	  sxxL= ue(3)
	  sigmaL=-pL+sxxL
	 call state_p_to_ei(rhoL,pL,eL)
	 call sound(ue,cL)

	 call trans_u_to_ue(uR(i,:),ue(:))
	  rhoR= ue(0)
	  uuR = ue(1)
	  pR =  ue(2)
	  sxxR= ue(3)
	  sigmaR=-pR+sxxR
	  call state_p_to_ei(rhoR,pR,eR)
	  call sound(ue,cR)
	  
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

	 !if (dabs(sxxL).lt.2.d0/3*Y0.and.abs(tmp).ge.2.d0/3*Y0) then
	 if (abs(tmp).ge.2.d0/3*Y0) then
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
	 !if (dabs(sxxR).lt.2.d0/3*Y0.and.abs(tmp).ge.2.d0/3*Y0) then
	 if (abs(tmp).ge.2.d0/3*Y0) then
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

	if (S_star.ge.u(i))then
		h(i,0)=0
		h(i,1)=pL_star-sxxL_star
		h(i,2)=(pL_star-sxxL_star)*s_star
		h(i,3)=-4.d0*miu/3*s_star
		u_half(i)=s_star
	else 
		h(i,0)=0
		h(i,1)=pR_star-sxxR_star
		h(i,2)=(pR_star-sxxR_star)*s_star
		h(i,3)=-4.d0*miu/3*s_star
		u_half(i)=s_star
	endif

! if (S_star.le.u(i))then
!		h(i,0)=0
!		h(i,1)=pL-sxxL
!		h(i,2)=(pL-sxxL)*uuL
!		h(i,3)=-4.d0*miu/3*s_star
!		u_half(i)=s_star
!	else 
!		h(i,0)=0
!		h(i,1)=pR-sxxR
!		h(i,2)=(pR-sxxR)*uuR
!		h(i,3)=-4.d0*miu/3*s_star
!		u_half(i)=s_star
!	endif

enddo

endsubroutine


