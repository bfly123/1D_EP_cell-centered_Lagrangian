subroutine bound_sod_problem(u4)
 use global
    use init_tran
      use global_cont
      implicit none
      integer i,j,k
      double precision u4(-nv:jx+nv,0:3) 
      double precision  f_eta,rho,uu,p,sxx

!ÓÒ±ß½ç
	DO I=1,nv
        U(jx+i,:)=U(jx-i,:) 
        U(jx,1)=0
	enddo
         
!×ó±ß½ç  
	   DO I=-nv,-1
	   rho=u(0,0)
	   uu=u(0,1)/rho
	  u(i,0)=rho
	  u(i,1)=uu*rho
      U(i,2)=(p1-rho0*a0**2*f_eta(rho))/(rho0*gamma0)*rho0+0.5d0*rho*uu*uu
	  u(i,3)=0

	   

	   enddo
    
	END
	
	
