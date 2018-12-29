		subroutine state_e_to_p(ei,rho,p)
		use global_cont
		implicit none 
		double precision rho,p,f_eta,ei
		!mie-Gruneisen equation  :  p = rho0 a0^2 f(eta) +rho0 Gamma0 ei
		p=rho0*a0**2*f_eta(rho) + rho0*Gamma0*ei
		end

		subroutine  state_de_to_dp(dei,rho,drho,dp) !p=p(e,rho)
		use global_cont
		implicit none 
		double precision rho,p,f_eta_eta,dei,dp,drho
		!mie-Gruneisen equation  :  p = rho0 a0^2 f(eta) +rho0 Gamma0 ei
		dp=a0**2*f_eta_eta(rho)*drho+rho0*Gamma0*dei
		end

		subroutine state_p_to_ei(rho,p,ei)
		use global_cont
		implicit none 
		double precision rho,p,f_eta,ei
		!mie-Gruneisen equation  :  p = rho0 a0^2 f(eta) +rho0 Gamma0 ei
		ei=(p-rho0*a0**2*f_eta(rho))/(rho0*Gamma0)
		end


	function f_eta(rho)
	  use global_cont
	  implicit none
	  double precision:: rho,eta,f_eta

	  eta=rho/rho0
	  f_eta= (eta -1.d0)*(eta-gamma0*(eta-1.d0)/2.d0)/(eta-s0*(eta-1))**2

	  end function

	 function f_eta_eta(rho)
	 use global_cont
	 implicit none
	 double precision rho,f_eta_eta,eta
	  eta=rho/rho0
	  f_eta_eta= (eta+(s0-gamma0)*(eta-1))/(eta-s0*(eta-1))**3
	  end function

	 function f_eta_eta_eta(rho)
	 use global_cont
	 implicit none
	 double precision rho,f_eta_eta_eta,eta
	  eta=rho/rho0
	  f_eta_eta_eta=((1+s0-gamma0)*(eta-s0*(eta-1))-3.d0*(1-s0)*(eta+(s0-gamma0)*(eta-1)))/(eta-s0*(eta-1))**4
	  end function

	 subroutine state_choose(i)
	 use global_cont
	 implicit none
	 integer i

	 select case(i)

	 case(2)

	Y0=3.d8
	rho0=2785
	gamma0=2.d0
	miu=2.76d10
	a0=5328
	pi=6*dasin(0.5d0)
	s0=1.338d0

	case(1) 

		Y0=9.d7
		rho0=8930
		gamma0=2.d0
		miu=4.5d10
		a0=3940
		pi=6*dasin(0.5d0)
		s0=1.49d0

		endselect
end





