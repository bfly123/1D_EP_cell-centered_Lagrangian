		subroutine state_e_to_p(ei,rho,p)
		use global_cont
		implicit none 
		double precision rho,p,f_eta,ei

		!mie-Gruneisen equation  :  p = rho0 a0^2 f(eta) +rho0 Gamma0 ei
		p=rho0*a0**2*f_eta(rho) + rho0*Gamma0*ei
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


