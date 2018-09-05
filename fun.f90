function f_eta(rho)
	  use global_cont
	  implicit none
	  double precision:: rho,eta,f_eta

	  eta=rho/rho0
	  f_eta= (eta -1.d0)*(eta-gamma0*(eta-1.d0)/2.d0)/(eta-s0*(eta-1))**2

	  end function

function fgamma(sxx)
	  use global_cont
	implicit none
	  double precision sxx,fgamma
	  if (abs(sxx).le.2.d0/3*Y0) then
		  fgamma =sxx 
	  else if (sxx.ge. 2.d0/3*Y0)  then
		  fgamma=2.d0/3*Y0
		else 
		  fgamma=-2.d0/3*Y0
	endif
end function
 function f_eta_eta(rho)
	 use global_cont
	 implicit none
	 double precision rho,f_eta_eta,eta
	  eta=rho/rho0
	  f_eta_eta= (eta+(s0-gamma0)*(eta-1))/(eta-s0*(eta-1))**3
	  end function



       
