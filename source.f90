subroutine Source1(t,src)
	  use global
	  use global_cont
	  implicit none
	  double precision eta, a,b, s, s20, s2, s3,t,xt,f_eta_eta,rho_x,sxx_x,xx
	  double precision src(-nv:jx+nv,0:3)
	  double precision  dx,x1,x2,x3,s21,s22,s23
	  integer i

	  src=0
	  select case(kind_problem)

	  case(4)

	  a=1.d3
	  b=0.2d0
	  s=6.d4

	 do i =-nv,jx+nv

	 dx = x(i)-x(i-1)
	 x1=(x(i)+x(i-1))/2

	  eta=1-b*sin(2*pi*(x1-a*t))
	  rho_x=-2*pi*rho0*b*cos(2*pi*(x1-a*t))
	  sxx_x=2*pi*s*cos(2*pi*(x1-a*t))
	  s21=a0**2*f_eta_eta(eta*rho0)*rho_x-sxx_x
	  s21=-2*pi*s*cos(2*pi*(x1-a*t))-1.5d0*rho0*a0**2*b*(2*pi*cos(x1-a*t))

	  x2=x1+sqrt(0.6)*(dx/2)

	  eta=1-b*sin(2*pi*(x2-a*t))
	  rho_x=-2*pi*rho0*b*cos(2*pi*(x2-a*t))
	  sxx_x=2*pi*s*cos(2*pi*(x2-a*t))
	  s22=a0**2*f_eta_eta(eta*rho0)*rho_x-sxx_x
	  s22=-2*pi*s*cos(2*pi*(x2-a*t))-1.5d0*rho0*a0**2*b*(2*pi*cos(x2-a *t))

	  x3=x1- sqrt(0.6)*(dx/2)

	 eta=1-b*sin(2*pi*(x3-a*t))
	 rho_x=-2*pi*rho0*b*cos(2*pi*(x3-a*t))
	  sxx_x=2*pi*s*cos(2*pi*(x3-a*t))
	  s22=-2*pi*s*cos(2*pi*(x3-a*t))-1.5d0*rho0*a0**2*b*(2*pi*cos(x3-a*t))
	  !s23=a0**2*f_eta_eta(eta*rho0)*rho_x-sxx_x

	 ! s23=-2*pi*s*cos(2*pi*(x3-a*t))-1.5d0*rho0*a0**2*b*(2*pi*cos(x3-a*t))

	  s2=(s21*8.d0/9+s22*5.d0/9+s23*5.d0/9)*dx/2
	 s2= s21*dx

	  s3=s2*a

	 !@x1=(x(i)+x(i-1))/2
	 !eta=1-b*sin(2*pi*(x1-a*t))
	  !s21=(1-gamma0/2)*(eta-1)/(eta-s0*(eta-1))**2-&
	  !  	2.d0*(1-s0)*(eta-0.5d0*gamma0*(eta-1))*(eta-1)/(eta-s0*(eta-1))**3&
	  !  	+(eta-0.5d0*gamma0*(eta-1))/(eta-s0*(eta-1))**2
	  !s2=s21*(-2*b*pi*cos(2*pi*(x1-a*t))) -2*pi*s*sin(2*pi*(x1-a*t))
	  !s2=-2*pi*s*cos(2*pi*(x1-a*t))
	  !s3=a*s2

	  src(i,0)=0
	  src(i,1)=s2 !*dx
	  src(i,2)=s3!*dx
	  src(i,3)=0
	  
	  enddo
  endselect

!  call output1(src)
!
!  pause
!
src=0
!  write(*,*) src(:,1)
!  pause

	  end
















