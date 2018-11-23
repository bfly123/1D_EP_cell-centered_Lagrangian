subroutine bound_Piston_problem(u4)
 use global
    use init_tran
      use global_cont
      implicit none
      integer i,j,k
      double precision u4(-nv:jx+nv,0:3) 
      double precision ue(0:3) 
      double precision  f_eta,rho,uu,p,sxx

!ÓÒ±ß½ç
	DO I=1,nv
	  uo(jx+i,0)=uo(jx-i,0)
	  uo(jx+i,1)=-uo(jx-i,1)
      Uo(jx+i,2)= uo(jx-i,2)
	  uo(jx+i,3)=u(jx-i,3)

	  call trans_Ue_to_u(uo(jx+i,:),u(jx+i,:))

	enddo
         
!×ó±ß½ç  
	   DO I=-nv,-1
!		rho=u(0,0)
!	  uu=u(0,1)/rho
!	  p = (u(0,2)/rho- 0.5*uu**2)*rho0*gamma0+rho0*a0**2*f_eta(rho)
!	  uu=u2
	  !u(i,:)=u(1,:)
!	  u(i,1)=uu*rho
!      U(i,2)=(p-rho0*a0**2*f_eta(rho))/(rho0*gamma0)*rho+0.5d0*rho*uu*uu
!	  u(i,3)=0 !u(0,3)
		  uo(i,:)=uo(0,:)
		  uo(i,1)=u2
	  call trans_ue_to_u(uo(i,:),u(i,:))
	  !u(i,3)= 0 !u(0,3)

   	   enddo
    
	END
	
	
