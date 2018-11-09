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

subroutine trans_du_to_due(u,du,due)
	  !use global_cont
	  implicit none
	  double precision u(0:3)
	  double precision du(0:3)
	  double precision ue(0:3)
	  double precision due(0:3)
	  double precision f_eta,ei,dei

		call trans_u_to_ue(u,ue)
		due(0)=du(0)
		due(1)=du(1)/u(0)-u(1)*du(0)/u(0)**2
		ei=u(2)/u(0)-0.5d0*ue(1)**2
		dei=du(2)/u(0)- u(2)/u(0)*du(0)-ue(1)*due(1)
		call state_de_to_dp(dei,u(0),du(0),due(2)) !p=p(e,rho)
		due(3)=du(3)
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
		F(3)=-4.d0*miu/3*ue(1)
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

		subroutine sound(ue,c)
		use global_cont
		implicit none 
		double precision ue(0:3)
		double precision c,rho,uu,p,sxx,a_squre,f_eta_eta

		rho=ue(0)
		uu=ue(1)
		p=ue(2)
		sxx=ue(3)

	  a_squre=a0**2 *f_eta_eta(rho) + p/rho**2*rho0 *gamma0
	  c=sqrt(a_squre-rho0/rho**2*gamma0*sxx+4.d0/3*miu/rho)
	  end






	  

