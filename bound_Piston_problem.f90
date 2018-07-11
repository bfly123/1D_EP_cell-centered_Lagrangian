subroutine bound_Piston_problem(u4)
 use global
    use init_tran
      use global_cont
      implicit none
      integer i,j,k
      double precision u4(-nv:jx+nv,0:3) 
      double precision  f_eta,rho,uu,p,sxx

!�ұ߽�
	DO I=1,nv
		rho=u(jx-i,0)
		uu=u(jx-i,1)/rho
	  u(jx+i,0)=rho
	  u(jx+i,1)=-uu*rho
      U(jx+i,2)= u(jx-i,2)
	  u(jx+i,3)=u(jx-i,3)
	enddo
         
!��߽�  
	   DO I=-nv,-1
	   rho=u(-i,0)
	   uu=u(-i,1)/rho
	   p=(u(-i,2)/rho-0.5*uu**2)*(rho0*gamma0)+rho0*a0**2*f_eta(rho)

	  u(i,0)=rho
	  u(i,1)=u2*rho
      U(i,2)=(p-rho0*a0**2*f_eta(rho))/(rho0*gamma0)*rho+0.5d0*rho*u2*u2
	  u(i,3)=u(-i,3)

	   enddo
    
	END
	
	