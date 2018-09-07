subroutine eigen_var(u,Ar,AL)
use global_cont
	  implicit none 
double precision u(0:3),ue(0:3),u1(0:3)
double precision AR(0:3,0:3)
double precision A(0:3,0:3)
double precision AR1(4,4)
double precision AL(0:3,0:3)
double precision AL1(4,4)
double precision  rho,uu,p,sxx,f_eta_eta,cc,a2, gamma1,b1,h,phi2,f_eta,s2,c2,ei,ee,ee1
integer i,j,info,k,m


	  rho=u(0)
      uu=U(1)/rho
	  ee=u(2)/rho !E
	  ei=ee-0.5*uu**2 !e
	  ee1=ei-0.5d0*uu**2
	  p = ei*rho0*gamma0+rho0*a0**2*f_eta(rho)

	  sxx=u(3)
	  a2=a0**2 *f_eta_eta(rho)  !+ p/rho**2 *rho0 *gamma0
	  s2=a2+ p/rho**2 *rho0 *gamma0-rho0/rho**2*gamma0*sxx
	  c2=s2+4.d0/3*miu/rho
	  !cc=sqrt(a2-rho0/rho**2*gamma0*sxx+4.d0/3*miu/rho)
	  cc=sqrt(c2)
	  gamma1=gamma0*rho0/rho

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

      Ar(0,2)= 1.d0/phi2
	  Ar(1,2)= 1.d0/phi2 *(uu-cc)
	  Ar(2,2)= 1.d0/phi2*(h-uu*cc)
	  Ar(3,2)= 1.d0

      Ar(0,3)= 1.d0/phi2
	  Ar(1,3)= 1.d0/phi2 *(uu+cc)
	  Ar(2,3)= 1.d0/phi2*(h+uu*cc)
	  Ar(3,3)= 1.d0

!	  do i=1,4
!	  do j =1,4
!		Ar1(i,j)=Ar(i-1,j-1)
!enddo
!enddo

		AL(0,0)= (1.d0-s2/c2)*(a2-gamma1*ee1)
                         
		AL(1,0)= (4.d0*a2**2-4*a2*(s2+2.d0*ei*gamma1)     &
		            +gamma1*(4.d0*ei*s2+4.d0*ei**2.d0*gamma1-uu**2*(-4*c2+2.d0*s2+gamma1*uu**2)))/(4.d0*c2*gamma1)
		         
		AL(2,0)= -0.5d0*(1-s2/c2)*(a2+cc*uu- ee1*gamma1)
		AL(3,0)=  0.5d0*(1-s2/c2)*(-a2+cc*uu+ ee1*gamma1)
                         
		AL(0,1)= -(1.d0-s2/c2)*uu*gamma1
		AL(1,1)= (-a2-c2+s2+ee*gamma1)*uu/c2
		AL(2,1)=  0.5d0*(1.d0-s2/c2)*(cc+uu*gamma1)
		AL(3,1)= -0.5d0*(1.d0-s2/c2)*(cc-uu*gamma1)
                         
		AL(0,2)=  (1.d0-s2/c2)*gamma1
		AL(1,2)= (a2+c2-s2-ee*gamma1)/c2
		AL(2,2)= -0.5d0*(1.d0-s2/c2)*gamma1
		AL(3,2)= -0.5d0*(1.d0-s2/c2)*gamma1
                         
		AL(0,3)=  s2/c2
		AL(1,3)= (-a2+s2+ee*gamma1)/c2/gamma1
		AL(2,3)=  0.5d0*(1.d0-s2/c2)
		AL(3,3)=  0.5d0*(1.d0-s2/c2)

!call reverse(Ar1,am1,4)
A=0
!	  do i=0,3
!	  do j =0,3
!	  do k=0,3
!		 A(i,j)= A(i,j)+AL(i,k)*AR(k,j) 
!		 enddo
!		enddo
!		enddo
!write(*,*) A(0:3,0)
!write(*,*) A(0:3,1)
!write(*,*) A(0:3,2)
!write(*,*) A(0:3,3)
!!enddo
!write(*,*) matmul(AR,AL)
!write(*,*)AR(2,3)
!pause

	end



