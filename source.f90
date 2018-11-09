subroutine source(t,src)
	  use global
	  use global_cont
	  implicit none

	  double precision scr(-nv:jx+nv,0:3)

	  scr=0
	  select case(kind_problem)

	  case(4)

	  a=1.d4
	  b=0.2d0
	  

	  eta=1-b*sin(2*pi*(x-a*t))
	  s21=(1-gamma0/2)*(eta-1)/(eta-s0*(eta-1))**2-&
			2.d0*(1-s0)*(eta-0.5d0*gamma0*(eta-1))*(eta-1)/(eta-s0*(eta-1))**3&
			+(eta-0.5d0*gamma0*(eta-1))/(eta-s0*(eta-1))**2
	  s2=s21*(-2*b*pi*cos(2*pi*(x(i)-a*t))) -2*pi*s0*sin(2*pi*(x-a*t))
	  s3=a*s2








