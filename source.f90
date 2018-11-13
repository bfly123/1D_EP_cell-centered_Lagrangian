subroutine Source1(t,src)
	  use global
	  use global_cont
	  implicit none
	  double precision eta, a,b, s, s21, s2, s3,t,xt,f_eta_eta,rho_x,sxx_x
	  double precision src(-nv:jx+nv,0:3)
	  integer i

	  src=0
	  select case(kind_problem)

	  case(4)

	  a=1.d4
	  b=0.1d0
	  s=6.d7

	  do i =-nv,jx+nv

	 
	  eta=1-b*sin(2*pi*(x(i)-a*t))
	  rho_x=-2*pi*rho0*b*cos(2*pi*(x(i)-a*t))
	  sxx_x=2*pi*s*cos(2*pi*(x(i)-a*t))
	  s2=a0**2*f_eta_eta(eta*rho0)*rho_x-sxx_x
	  s3=a*s2


	 ! s2=-(2*pi*b*rho0*a0**2*f_eta_eta(eta)+2*pi*s)* cos(2*pi*(x(i)-a*t))
	 !! s3=s2*a

	 ! s21=(1-gamma0/2)*(eta-1)/(eta-s0*(eta-1))**2-&
	 !   	2.d0*(1-s0)*(eta-0.5d0*gamma0*(eta-1))*(eta-1)/(eta-s0*(eta-1))**3&
	 !   	+(eta-0.5d0*gamma0*(eta-1))/(eta-s0*(eta-1))**2
	 ! s2=s21*(-2*b*pi*cos(2*pi*(x(i)-a*t))) -2*pi*s*sin(2*pi*(x(i)-a*t))
	 ! s3=a*s2
	  src(i,0)=0
	  src(i,1)=s2
	  src(i,2)=s3
	  src(i,3)=0
	  
	  enddo
  endselect

!  write(*,*) src(:,1)
!  pause

	  end
















