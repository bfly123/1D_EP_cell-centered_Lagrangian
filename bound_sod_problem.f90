subroutine bound_sod_problem(u4)
 use global
    use init_tran
      use global_cont
      implicit none
      integer i,j,k
      double precision u4(-nv:jx+nv,0:3) 
      double precision  f_eta

!ÓÒ±ß½ç
	DO I=1,nv
        U(jx+i,:)=U(jx-i,:) 
        U(jx,1)=0
	enddo
         
!×ó±ß½ç  
	   DO I=-nv,-1
	   U4(i,:)=U4(0,:)
	!		U(i,0)=rho1
    !        U(i,1)=rho1*u1
    !        U(i,2)=(p1-rho0*a0**2*f_eta(rho1))/(rho0*gamma0)*rho1+0.5*rho1*u1**2
	!		U(i,3)=0
      
	   enddo
    
	END
	
	
