subroutine trans_u_to_ue(u,ue)
	  !use global_cont
	  implicit none
	  double precision u(0:3)
	  double precision ue(0:3)
	  double precision f_eta,ei

		ue(0)=u(0)
		ue(1)=u(1)/u(0)
		ei=u(2)/u(0)-0.5d0*ue(1)**2
		call state_e_to_p(ei,u(0), ue(2)) !p=p(e,rho)
		ue(3)=u(3)
		end


subroutine trans_ue_to_u(ue,u)
	  use global_cont
	  implicit none
	  double precision u(0:3)
	  double precision ue(0:3)
	  double precision f_eta,ei

		u(0)=ue(0)
		u(1)=ue(1)*ue(0)
		call state_p_to_ei(u(0),ue(2),ei)
		u(2) = (ei+0.5d0*ue(1)**2)*ue(0)
		u(3)=ue(3)
		end

subroutine trans_ue_to_FLagrangian(ue,F)
	  use global_cont
	  implicit none
	  double precision ue(0:3)
	  double precision F(0:3)
		F(0)=0
		F(1)=ue(2)-ue(3)
		F(2) = F(1)*ue(1) 
		F(3)=-4.d0*miu/3*ue(3)
		end
subroutine trans_ue_to_Feuler(ue,F)
	  use global_cont
	  implicit none
	  double precision ue(0:3)
	  double precision F(0:3)
	  double precision f_eta,ei,ee
		F(0)=ue(0)*ue(1)
		F(1)=ue(0)*ue(1)**2+ue(2)-ue(3)
		call state_p_to_ei(ue(0),ue(2),ei)
		ee=ei+0.5d0*ue(1)**2
		F(2) = (ue(0)*ee+ue(2)-ue(3))*ue(1)
		F(3)=-4.d0*miu/3*ue(3)
		end


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
		ei=p-rho0*a0**2*f_eta(rho)/(rho0*Gamma0)
		end





	  

