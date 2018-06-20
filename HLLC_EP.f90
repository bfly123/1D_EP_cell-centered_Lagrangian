subroutine HLLC_EP(nv,jx,u,ul,ur,h,u_half)
	  use global_cont
	  implicit none
	  double precision  f_eta,fgamma 
	  integer nv,i,jx
	  double precision u(-nv:jx+nv)
	  double precision ul(-nv:jx+nv,0:3)
	  double precision ur(-nv:jx+nv,0:3)
	  double precision h(-nv:jx+nv,0:3)
	  double precision u_half(-nv:jx+nv)

	  double precision rhoL,uuL,sxxL,pL,cL,sL,sxxL_bar,sxxL_star,sxL_star,pL_star,sigmaxL_star
	  double precision rhoR,uuR,sxxR,pR,cR,sR,sxxR_bar,sxxR_star,sxR_star,pR_star,sigmaxR_star
	  double precision s_barStar,feta,feta1,feta_eta,a_squre,s_star,f_eta_eta

do i =-nv, jx+nv-1
	  rhoL= ul(i,0)
	  uuL = ul(i,1)/rhoL
	  sxxL= ul(i,3)
	  feta = f_eta(rhoL)
	  pL = (ul(i,2)/rhoL- 0.5*uuL**2)*rho0*gamma0+rho0*a0**2*feta
	  feta_eta=f_eta_eta(rhoL)
	  
	  !write(*,*)i, rhoL, feta,feta1,feta_eta 


	  a_squre=a0**2 *feta_eta + pL/rhoL**2 *rho0 *gamma0
	  cL=sqrt(a_squre-rho0/rhoL**2*gamma0*sxxL+4.d0/3*miu/rhoL)

	  rhoR= ur(i,0)
	  uuR = ur(i,1)/rhoR
	  feta = f_eta(rhoR)
	  pR = (ur(i,2)/rhoR-0.5*uuR**2)*rho0*gamma0+rho0*a0**2*feta
	  sxxR= ur(i,3)

	  feta_eta=f_eta_eta(rhoR)
	  a_squre=a0 **2 *feta_eta + pR/rhoR**2 *rho0 *gamma0
	  cR=sqrt(a_squre-rho0/rhoR**2*sxxR*gamma0+4.d0/3*miu/rhoR)

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

	u_half(i)=s_star
	pL_star=pL+rhoL*(sL-uuL)*(s_star+sxL_star-uuL)
	pR_star=pR+rhoR*(sR-uuR)*(s_star+sxR_star-uuR)
	
!	sigmaxL_star=sxxL_star-PL_star
!	sigmaxR_star=sxxR_star-PR_star
	
	if (s_star.ge.u(i))then
	!if (s_star.le.0)then
	!if (s_star.ge.))then
		h(i,0)=0
		h(i,1)=pL_star-sxxL_star
		h(i,2)=(pL_star-sxxL_star)*s_star
	else
		h(i,0)=0
		h(i,1)=pR_star-sxxR_star
		h(i,2)=(pR_star-sxxR_star)*s_star
	endif
	h(i,3)=-4*miu/3*s_star
enddo

endsubroutine

