      subroutine CFL(dt)
		  use global
		  use global_cont
      implicit none
      integer i
	  double precision dxmin, dt,dx
	  double precision rho,feta,uu,sxx,p,feta1,feta_eta,a_squre,c,f_eta,f_eta_eta
      
	  dxmin=x(1)-x(0)

      do i=0,Jx
	  dx=x(i)-x(i-1)

	  rho= uo(i,0)
	  !feta = f_eta(rho)
	  uu = uo(i,1)

	  sxx=uo(i,3)
	  p = uo(i,2) !(u(i,2)/rho- 0.5*uu**2)*rho0*gamma0+rho0*a0**2*feta

	  call sound(uo(i,:),c)
	  !feta_eta=f_eta_eta(rho)
	  !a_squre=a0**2 *feta_eta + p/rho**2 *rho0 *gamma0
	  !c=sqrt(a_squre-rho0/rho**2*gamma0*sxx+4.d0/3*miu/rho)
		if( dx.le.dxmin ) then
			dxmin = dx
		endif
	enddo
	  
     dt=SF*dxmin/c
      end
