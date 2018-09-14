subroutine eigen_var(u,AL,Ar)
use global_cont
	  implicit none 
double precision u(0:3),ue(0:3),u1(0:3)
double precision AR(0:3,0:3),A(0:3,0:3),AL(0:3,0:3)
double precision AR1(4,4),AL1(4,4),A1(4,4)
double precision  rho,uu,p,sxx,f_eta_eta,cc,a2, gamma,b1,h,phi2,f_eta,s2,c2,ei,ee,ee1,p_d
integer i,j,info,k,m,n

     call trans_u_to_ue(u,ue) 
	  rho=ue(0)
      uu=Ue(1)
	  p= ue(2)
	  sxx=ue(3)

	  
	  a2=a0**2 *f_eta_eta(rho) 
	  phi2=-4.d0*miu/3.d0/rho
	  gamma = gamma0*rho0/rho
	  ee=Ue(2)-uu**2/2
	  
	  A(0,0)=0 
	  A(1,0)=1
	  A(2,0)=0
	  A(3,0)=0

	  A(0,1)=gamma*(uu**2/2-ee)+a2-uu**2 
	  A(1,1)= uu*(2-gamma)
	  A(2,1)=gamma
	  A(3,1)=-1
	  
	  A(0,2)=-uu*(ee*gamma-uu**2/2*gamma+uu**2/2)+a2*uu
	  A(1,2)=uu**2*(0.5d0-gamma)
	  A(2,2)=uu*(1+gamma)
	  A(3,2)=-uu


	  A(0,3)= -phi2*uu
	  A(1,3)= phi2
	  A(2,3)= 0 
	  A(3,3)= uu

	  do i=0,3
		do j=0,3
			A1(i+1,j+1)=A(i,j)
		enddo
	enddo
n=4
	call MKL_AL_AR(A1,AR1,n)
pause
    call reverse(AR1,AL1,4)
    
	  do i=0,3
		do j=0,3
			AR(i,j)=AR1(i+1,j+1)
			AL(i,j)=AL1(i+1,j+1)
		enddo
		enddo
		write(*,*) matmul(AR,AL)
 end



subroutine eigen_var_R(u,Ar)
use global_cont
	  implicit none 
double precision u(0:3),ue(0:3),u1(0:3)
double precision AR(0:3,0:3)
double precision  rho,uu,p,sxx,f_eta_eta,cc,a2, gamma1,b1,h,phi2,f_eta,s2,c2,ei,ee,ee1,p_d
integer i,j,info,k,m

     call trans_u_to_ue(u,ue) 
	  rho=ue(0)
      uu=Ue(1)
	  p= ue(2)
	  sxx=ue(3)

	  ee=u(2)/rho !E
	  ei=ee-0.5*uu**2 !e
	  ee1=ei-0.5d0*uu**2
	  gamma1=gamma0*rho0/rho

	  a2=a0**2 *f_eta_eta(rho) 
	  s2=a2+ (p-sxx)/rho**2*rho0*gamma0
	  c2=s2+4.d0/3*miu/rho
	  cc=dsqrt(c2)
	  b1=a2-gamma1 *ee 
	  h=ee + (p - sxx)/rho

	  phi2=-4.d0*miu/3.d0/rho

	  Ar(0,0)= 1.d0/b1
	  Ar(1,0)= uu/b1
	  Ar(2,0)= 0.d0
	  Ar(3,0)= 1.d0

      Ar(0,1)= -gamma1/b1
	  Ar(1,1)= -gamma1*uu/b1
	  Ar(2,1)= 1.d0
	  Ar(3,1)= 0.d0

      Ar(0,2)= 1.d0/(s2-c2)
	  Ar(1,2)= 1.d0/(s2-c2) *(uu-cc)
	  Ar(2,2)= 1.d0/(s2-c2)*(h-uu*cc)
	  Ar(3,2)= 1.d0

      Ar(0,3)= 1.d0/(s2-c2)
	  Ar(1,3)= 1.d0/(s2-c2)*(uu+cc)
	  Ar(2,3)= 1.d0/(s2-c2)*(h+uu*cc)
	  Ar(3,3)= 1.d0
	  end

subroutine eigen_var_L(AR,AL)
use global_cont
	  implicit none 
double precision u(0:3),ue(0:3),u1(0:3)
double precision AR(0:3,0:3)
double precision A(0:3,0:3)
double precision AR1(4,4)
double precision AL(0:3,0:3)
double precision AL1(4,4)
integer i,j

     !call eigen_var_OR(u,AR)
	 do i=0,3
		do j=0,3
		AR1(i+1,j+1)=AR(i,j)
		enddo
	enddo

	 call reverse(AR1,AL1,4)
	do i=0,3
		do j=0,3
		AL(i,j)=AL1(i+1,j+1)
		enddo
	enddo
end



subroutine eigen_var_OR(u,Ar)
use global_cont
	  implicit none 
double precision u(0:3)
double precision AR(0:3,0:3)
double precision  rho,uu,p,sxx,cc,S2,phi2

	  rho=u(0)
      uu=U(1)
	  p= u(2)
	  sxx=u(3)
	  call sound(u,cc)

	  phi2=-4.d0*miu/3.d0/rho
	  S2=cc**2+phi2

	  Ar(0,0)= 0
	  Ar(1,0)= 0
	  Ar(2,0)= 1.d0
	  Ar(3,0)= 1.d0

      Ar(0,1)= 1.d0
	  Ar(1,1)= 0
	  Ar(2,1)= 0
	  Ar(3,1)= 0

      Ar(0,2)= 1.d0/phi2
	  Ar(1,2)= -cc/rho /phi2
	  Ar(2,2)=  s2/phi2
	  Ar(3,2)= 1.d0

      Ar(0,3)= 1.d0/phi2
	  Ar(1,3)= cc/rho /phi2
	  Ar(2,3)=  s2/phi2
	  Ar(3,3)= 1.d0
	  end

subroutine eigen_var_OL(AR,AL)
use global_cont
	  implicit none 
double precision u(0:3),ue(0:3),u1(0:3)
double precision AR(0:3,0:3)
double precision A(0:3,0:3)
double precision AR1(4,4)
double precision AL(0:3,0:3)
double precision AL1(4,4)
integer i,j

     !call eigen_var_OR(u,AR)
	 do i=0,3
		do j=0,3
		AR1(i+1,j+1)=AR(i,j)
		enddo
	enddo

	 call reverse(AR1,AL1,4)
	do i=0,3
		do j=0,3
		AL(i,j)=AL1(i+1,j+1)
		enddo
	enddo
end



