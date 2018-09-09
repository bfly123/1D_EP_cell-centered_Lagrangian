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

subroutine eigen_var_L(u,AL)
use global_cont
	  implicit none 
double precision u(0:3),ue(0:3),u1(0:3)
double precision AR(0:3,0:3)
double precision A(0:3,0:3)
double precision AR1(4,4)
double precision AL(0:3,0:3)
double precision AL1(4,4)
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
	  gamma1= gamma0*rho0/rho

	  a2=a0**2 *f_eta_eta(rho) 
	  phi2=-4.d0*miu/3.d0/rho
	  s2=a2+ (p-sxx)/rho**2*rho0*gamma0
	  c2=s2+4.d0/3*miu/rho
	  cc=dsqrt(c2)

	  b1=a2-gamma1 *ee 
	  h=ee + (p - sxx)/rho

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

	end


