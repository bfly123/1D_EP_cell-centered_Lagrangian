subroutine LF_splitting(u4,fp,fm)
use global
use global_cont
implicit none
integer i,j,k
double precision uu,uv,rho,p,cc,u_max,u_max_tmp,sxx,f_eta,f_eta_eta,a_squre
double precision uc1,uc2,vc1,vc2,vvc1,vvc2,vv,w2,feta_eta
double precision u4(-nv:jx+nv,0:3)
double precision f(-nv:jx+nv,0:3)
double precision fp(-nv:jx+nv,0:3)
double precision fm(-nv:jx+nv,0:3)
integer nx1,ny1
!c El 为特征值,其中El(:,1)为x方向的特征值（4个， u, u, u+c, u-c)
      
	  u_max=0
      do i=-nv,jx+nv
	  uu=u4(i,1)/u4(i,0)
	  rho=u4(i,0)
	  p = (u4(i,2)/rho- 0.5*uu**2)*rho0*gamma0+rho0*a0**2*f_eta(rho)
	  sxx=u4(i,3)
	  feta_eta=f_eta_eta(rho)
	  a_squre=a0**2 *feta_eta + p/rho**2 *rho0 *gamma0
	  cc=sqrt(a_squre-rho0/rho**2*gamma0*sxx+4.d0/3*miu/rho)
	  
      u_max_tmp=abs(uu)+cc
	  if (u_max_tmp.ge.u_max)then
		  u_max=u_max_tmp
	endif
	  enddo

	  do i=-nv,jx+nv
        
	  uu=u4(i,1)/u4(i,0)
	  rho=u4(i,0)
	  p = (u4(i,2)/rho- 0.5*uu**2)*rho0*gamma0+rho0*a0**2*f_eta(rho)
	  sxx=u4(i,3)
	 

     f(i,0)=0
     f(i,1)=p-sxx
     f(i,2)=(p-sxx)*uu
     f(i,3)=-4.d0/3*miu*uu
    
      fp(i,0)=0
      fm(i,0)=0
      fp(i,1)=1.d0/2*f(i,1)
      fm(i,1)=1.d0/2*f(i,1)
      fp(i,2)=1.d0/2*(f(i,2)+uu*u4(i,2))
      fm(i,2)=1.d0/2*(f(i,2)-uu*u4(i,2))
      fp(i,3)=1.d0/2*f(i,3)
      fm(i,3)=1.d0/2*f(i,3)
     enddo
    
   
       end

