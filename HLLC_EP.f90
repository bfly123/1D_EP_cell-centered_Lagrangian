subroutine HLLC_EP(nv,jx,u,ul,ur,h)
	  implicit none
	  integer nv
	  integer jx
	  double precision u(-nv:jx+nv,0:3)
	  double precision ul(-nv:jx+nv,0:3)
	  double precision ur(-nv:jx+nv,0:3)
	  double precision h(-nv:jx+nv,0:3)
	  double precision hp(-nv:jx+nv,0:3)
	  double precision hm(-nv:jx+nv,0:3)


do i =-1, jx
	  rhoL= ul(i,0)
	  uL = ul(i,1)/rhoL
	  sxxL= ul(i,3)
	  feta = f_eta(rhoL)
	  pL = (ul(i,2)- 0.5*uL**2)*rho0*gamma0+rho0*a0**2*feta
	  feta1 = f_eta(rhoL*(1+0.01))
	  feta_eta=(feta1-feta)/(0.01*rhoL)*rho0

	  a_squre=a_0 **2 *feta_eta + pL/rhoL**2 *rho0 *gamma0
	  cL=sqrt(a_squre-rho0/rhoL**2*sxxL+4.d0/3*miu/rhoL)

	  rhoR= ur(i,0)
	  uR = ur(i,1)/rhoR
	  feta = f_eta(rhoR)
	  feta1 = f_eta(rhoR*(1+0.01))
	  pR = (ur(i,2)-0.5*uR**2)*rho0*gamma0+rho0*a0**2*feta
	  sxxR= ur(i,3)

	  feta_eta=(feta1-feta)/(0.01*rhoR)*rho0
	  a_squre=a_0 **2 *feta_eta + pR/rhoR**2 *rho0 *gamma0
	  cR=sqrt(a_squre-rho0/rhoR**2*sxxR+4.d0/3*miu/rhoR)
	  sL=min(uL-cL,uR-cR)
	  sR=min(uL+cL,uR+cR)
     s_barStar = (sxxR-sxxL)/(rhoL*(sL-uL)-rhoR*(sR-uR))
	sxxL_bar=sxxL+rhoL*(sL-uL)*s_barStar
	sxxR_bar=sxxR+rhoR*(sR-uR)*s_barStar
	sxxL_star=Fgamma(sxxL_bar)
	sxxR_star=Fgamma(sxxR_bar)
	
	sxL_star=(sxxL_star-sxxL)/(rhoL*(sL-uL))
	sxR_star=(sxxR_star-sxxR)/(rhoR*(sR-uR))
	
	s_star=(pR-pL+rhoL*(uL-sxL_star)*(sL-uL)-rhoR*(uR-sxR_star)*(sR-uR))/(rhoL*(sL-uL)-rhoR*(sR-uR))
	
	pL_star=pL+rhoL*(sL-uL)*(s_star+sxL_star-uL)
	pR_star=pR+rhoR*(sR-uR)*(s_star+sxR_star-uR)
	
	sigmaxL_star=sxxL_star-PL_star
	sigmaxR_star=sxxR_star-PR_star
	
	
	if (s_star.le.0)then

	h_G( : ) = hL(i,:)
	else 
	h_G(: ) = hR(i,:)
	endif

	hR(i,0)=0
	hR(i,1)=pR_star-sxxR_star
	hR(i,2)=(pR_star-sxxR_star)*uR

	hL(i,0)=0
	hL(i,1)=pL_star-sxxL_star
	hL(i,2)=(pL_star-sxxL_star)*uL

	if(h_G(i,1)/h_g(i,0).le.s_star) then
		h(i,0:2)= hL(i,0:2)
	else
		h(i,0:2)=hR(i,0:2)
	endif
		h(i,3)=4*miu/3*h_G(i,1)/h_G(i,0)
enddo

double precision function f_eta(rho)
	  implicit none
	  use global_cont
	  double precision:: rho,eta

	  eta=rho/rho0
	  f_eta= (eta -1)*(eta-gamma0*(eta-1)/2)/(eta-s*(eta-1))^2

	  end function

double precision function fgamma(sxx)
	implicit none
	  use global_cont
	  double precision sxx
	  if (sxx.ge.2.d0/3*Y0) then
		  fgamma = 2.d0/3*Y0
	 else
		  fgamma=sxx
	endif
end function



	  

