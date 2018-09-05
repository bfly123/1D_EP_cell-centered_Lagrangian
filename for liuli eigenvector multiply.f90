     
	 call multiply2x1(el,uu,ul,nn)

	 call multiply2x1(er,uu,ul,nn)


!   EE: TOTAL ENERGY PER VOLUME OF ONE CELL
subroutine LEFT_EIGENVECTORS_1D_ep(ei,PP,DD,uu,sxx,eigen,nser_eos,npa_eos,para_eos,shearm)
 use exact_riemann_solvers_realeos
IMPLICIT NONE
INTEGER KIND
PARAMETER(KIND=8)       
REAL(KIND),intent(out)  :: eigen(4,4)
REAL(KIND),intent(in)   :: ei,PP,DD,uu,sxx,shearm
integer,intent(in)      :: nser_eos,npa_eos
REAL(KIND),intent(in)   :: para_eos(20)
REAL(KIND) d0,a0,gamma0,gamma,ee,hh,dcc,epsn,p_d,p_e,a2,b1,cc,ee1,s2,c2,s

        d0     = para_eos(1)
        a0     = para_eos(2)
        gamma0 = para_eos(3)
        s      = para_eos(4)

        gamma = d0*gamma0/dd

        EE = ei+uu**2/2

		HH = ee+(pp-sxx)/dd      

        call derivative_eos(pp,dd,ei,nser_eos,para_eos,npa_eos,p_d,p_e) 
        a2 = p_d

		b1 = a2-gamma*ee

		EE1= ei-uu**2/2

		s2  = p_d+(pp-sxx)/dd**2*p_e
		c2  = s2+4.d0/3.d0*shearm/dd
		cc  = dsqrt(c2)

		eigen(1,1)= (1-s2/c2)*(a2-gamma*ee1)
                         
		eigen(2,1)= (4*a2**2-4*a2*(s2+2*ei*gamma)     &
		            +gamma*(4*ei*s2+4*ei**2*gamma-uu**2*(-4*c2+2*s2+gamma*uu**2)))/(4*c2*gamma)
		         
		eigen(3,1)= -0.5d0*(1-s2/c2)*(a2+cc*uu- ee1*gamma)
		eigen(4,1)=  0.5d0*(1-s2/c2)*(-a2+cc*uu+ ee1*gamma)
                         
		eigen(1,2)= -(1-s2/c2)*uu*gamma
		eigen(2,2)= (-a2-c2+s2+ee*gamma)*uu/c2
		eigen(3,2)=  0.5d0*(1-s2/c2)*(cc+uu*gamma)
		eigen(4,2)= -0.5d0*(1-s2/c2)*(cc-uu*gamma)
                         
		eigen(1,3)=  (1-s2/c2)*gamma
		eigen(2,3)= (a2+c2-s2-ee*gamma)/c2
		eigen(3,3)= -0.5d0*(1-s2/c2)*gamma
		eigen(4,3)= -0.5d0*(1-s2/c2)*gamma
                         
		eigen(1,4)=  s2/c2
		eigen(2,4)= (-a2+s2+ee*gamma)/c2/gamma
		eigen(3,4)=  0.5d0*(1-s2/c2)
		eigen(4,4)=  0.5d0*(1-s2/c2)

	return
	end

!   EE: TOTAL ENERGY PER VOLUME OF ONE CELL
subroutine RIGHT_EIGENVECTORS_1D_ep(ei,PP,DD,uu,sxx,eigen,nser_eos,npa_eos,para_eos,shearm)
 use exact_riemann_solvers_realeos
IMPLICIT NONE
INTEGER KIND
PARAMETER(KIND=8)       
REAL(KIND),intent(out)  :: eigen(4,4)
REAL(KIND),intent(in)   :: ei,PP,DD,uu,sxx
integer,intent(in)      :: nser_eos,npa_eos
REAL(KIND),intent(in)   :: para_eos(20)
REAL(KIND) d0,a0,gamma0,gamma,ee,hh,dcc,shearm,epsn,p_d,p_e,a2,b1,s2,c2,cc,s


        d0     = para_eos(1)
        a0     = para_eos(2)
        gamma0 = para_eos(3)
        s      = para_eos(4)

        gamma = d0*gamma0/dd

        EE = ei+uu**2/2
		HH = ee+(pp-sxx)/dd      


        call derivative_eos(pp,dd,ei,nser_eos,para_eos,npa_eos,p_d,p_e) 
        a2 = p_d

		s2  = p_d+(pp-sxx)/dd**2*p_e
		c2  = s2+4.d0/3.d0*shearm/dd
		cc  = dsqrt(c2)

		b1 = a2-gamma*EE

		dcc = 4.d0/3.d0*shearm/dd

		eigen(1,1)=  1.d0/b1
		eigen(2,1)=  uu/b1 
		eigen(3,1)=  0.d0
		eigen(4,1)=  1.d0

		eigen(1,2)=  -gamma/b1
		eigen(2,2)=  -gamma*uu /b1
		eigen(3,2)=  1.d0
		eigen(4,2)=  0.d0

		eigen(1,3)=1.d0/(s2-c2)
		eigen(2,3)=(cc-uu)/(c2-s2)
		eigen(3,3)=-(hh-uu*cc)/(c2-s2)
		eigen(4,3)=1.d0

		eigen(1,4)=1.d0/(s2-c2)
		eigen(2,4)=(uu+cc)/(s2-c2)
		eigen(3,4)=(hh+uu*cc)/(s2-c2)
		eigen(4,4)=1.d0

	return
	end

subroutine derivative_eos(pres,dens,enin,nser_eos,para_eos,npa_eos,p_d,p_e) 
!   calculate pressure's derivative to density
!   p_d is the derivative of p to density
!   p_e is the derivative of p to intel energy

include "precision.h"
dimension para_eos(20)

!       用于状态方程mie-gruieisen的检验
!       状态方程为
!       xx=dens/d0 
!       p=d0*a0^2*f(xx) + d0*gamma0*enin
!      f(xx) = (xx-1)*[xx-1/2*gamma0*(xx-1)]/[xx-s*(xx-1)]^2
         d0      = para_eos(1)
         a0      = para_eos(2)
         gamma0  = para_eos(3)
         s       = para_eos(4)
         xx      = dens/d0
		 xx_d    = 1.d0/d0

        p_d = (-2*a0**2*d0*(-1 + dens/d0)*(dens/d0 - ((-1 + dens/d0)*gamma0)/2.)*(1/d0 - s/d0))/(dens/d0 - (-1 + dens/d0)*s)**3 +      &
              (a0**2*d0*(-1 + dens/d0)*(1/d0 - gamma0/(2.*d0)))/(dens/d0 - (-1 + dens/d0)*s)**2 +                                      &
              (a0**2*(dens/d0 - ((-1 + dens/d0)*gamma0)/2.))/(dens/d0 - (-1 + dens/d0)*s)**2
	        
		 p_e     = d0 * gamma0

	return
	endsubroutine derivative_eos

subroutine multiply2x1(e1,e2,e3,nn)
	integer kind,nn,n,m
     parameter(kind=8)
	REAL(KIND),dimension(nn,nn) ::  e1
	REAL(KIND),dimension(nn)    ::  e2,e3

	do n=1,nn
        e3(n)=0.d0
        do k=1,nn
	      e3(n)=e3(n)+e1(n,k)*e2(k)
	    end do
	end do

	return
	end
