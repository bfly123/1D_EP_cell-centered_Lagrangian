subroutine eigen_var(u,ue,Ar)
use global_cont
	  implicit none 
double precision u(0:3),ue(0:3),u1(0:3)
double precision AR(0:3,0:3)
double precision A(0:3,0:3)
double precision AR1(4,4)
double precision AM(0:3,0:3)
double precision AM1(4,4)
double precision  rho,uu,p,sxx,f_eta_eta,cc,a_squre, gamma1,b1,h,phi2,f_eta
integer i,j,info,k,m



	  rho=u(0)
      uu=U(1)/rho
	  p = (u(2)/rho- 0.5*uu**2)*rho0*gamma0+rho0*a0**2*f_eta(rho)
	  sxx=u(3)
	  a_squre=a0**2 *f_eta_eta(rho) + p/rho**2 *rho0 *gamma0
	  cc=sqrt(a_squre-rho0/rho**2*gamma0*sxx+4.d0/3*miu/rho)

	  gamma1=gamma0*rho0/rho

	  b1=a0**2*f_eta_eta(rho)-gamma1 * u(2)/rho
	  h=u(2)/rho + (p - sxx)/rho

	  phi2=-4.d0*miu/3/rho


	  Ar(0,0)= 1/b1
	  Ar(1,0)= uu/b1
	  Ar(2,0)= 0
	  Ar(3,0)= 1

      Ar(0,1)= -gamma1/b1
	  Ar(1,1)= -gamma1*uu/b1
	  Ar(2,1)= 1
	  Ar(3,1)= 0

      Ar(0,2)= 1/phi2
	  Ar(1,2)= 1/phi2 *(uu-cc)
	  Ar(2,2)= 1/phi2*(h-uu*cc)
	  Ar(3,2)= 1

      Ar(0,3)= 1/phi2
	  Ar(1,3)= 1/phi2 *(uu+cc)
	  Ar(2,3)= 1/phi2*(h+uu*cc)
	  Ar(3,3)= 1

	  do i=1,4
	  do j =1,4
		Ar1(i,j)=Ar(i-1,j-1)
enddo
enddo
call reverse(Ar1,am1,4)

	  do i=1,4
	  do j =1,4
		Am(i-1,j-1)=Am1(i,j)
		enddo
		enddo
!ue=0	
!do i=0,3
!do j=0,3
!ue(i)=ue(i)+Am(i,j)*u(j)
!enddo
!enddo
ue=matmul(u,transpose(Am))
!write(*,*) matmul(u,transpose(Am))-ue
!read(*,*)i

	end



