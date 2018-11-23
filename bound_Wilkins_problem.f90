subroutine bound_Wilkins_problem(u4)
 use global
    use init_tran
      use global_cont
      implicit none
      integer i,j,k
      double precision u4(-nv:jx+nv,0:3) 
      double precision  f_eta,rho,uu,p,sxx,dx

!ÓÒ±ß½ç
	DO I=1,nv
	   !dx=x(jx)-x(jx-1)
	   !x(jx+i)=x(jx)+i*dx
        U(jx+i,:)=U(jx-i,:) 
        U(jx,1)=0

		Uo(jx+i,:)=Uo(jx-i,:)
		Uo(jx,1)=0
	enddo
         
!×ó±ß½ç  
	   DO I=-nv,-1
	   dx=x(1)-x(0)
		!x(i)=x(0)-i*dx 
	   rho=rho0
	   uu=u(0,1)/rho
	  u(i,0)=rho
	  u(i,1)=uu*rho
      U(i,2)=(-rho0*a0**2*f_eta(rho))/(rho0*gamma0)*rho0+0.5d0*rho*uu*uu
	  u(i,3)=0

	  uo(i,0)=rho0
	  uo(i,1)=uu
	  uo(i,2)=0
	  uo(i,3)=0

	   enddo
    
	END
	
	
