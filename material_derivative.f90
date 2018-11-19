!subroutine material_derivative(U1,pU_px,du_dt)
!use global_cont
!use global
!implicit none
!double precision U1(-nv:jx+nv,0:3)
!double precision Uo(-nv:jx+nv,0:3)
!double precision Ue(0:3)
!double precision dU_dt(-nv:jx+nv,0:3,3)
!double precision pU_px(-nv:jx+nv,0:3,3)
!double precision pUe_px(0:3)
!double precision pUe_pt(0:3)
!double precision pU_pt(0:3,3)
!double precision pU2_pxt(0:3)
!double precision A(0:3,0:3)
!double precision A_pA_px(0:3,0:3)
!double precision A2(0:3,0:3)
!double precision pA_px(0:3,0:3)
!double precision pA_pt(0:3,0:3)
!double precision puu_pt,puu_px,uu,f_eta_eta
!integer i,j,k
!
!
!!call trans_px_to_pt(pu_px,pu_pt)
!!call time_derivative(U_half,pU_pt,pU_px)
!do i=-nv,jx+nv
!
!	call trans_u_to_Ue(U1(i,0:3),ue)
!	call trans_ue_A(U1(i,0:3),ue,A)
!
!!*pU_pt=- A*pU_px
!	Pu_pt(0:3,1)=- matmul(A,pU_px(i,0:3,1))
!
!	call  trans_du_to_due(U1(i,0:3),pu_pt(0:3,1),pue_pt)  !**********\partial U/\parttial t to \partial U_e/\partial t
!	call  trans_du_to_due(U1(i,0:3),pu_px(i,0:3,1),pue_px) 
!	
!	call trans_ue_dA(U1(i,0:3),Ue,pU_pt(0:3,1),pUe_pt,pA_pt)
!	call trans_ue_dA(U1(i,0:3),Ue,pU_px(i,0:3,1),pUe_px,pA_px)
!	
!	A_pA_px=matmul(A,pA_px)
!
!	A2=matmul(A,A)
!
!	pU2_pxt(0:3)=- matmul(pA_px,pU_px(i,0:3,1))- matmul(A,pU_px(i,0:3,2))
!	pU_pt(0:3,2)=- matmul(pA_pt,pU_px(i,0:3,1))- matmul(A,pU2_pxt)
!
!	uu=Ue(1)
!
!	puu_pt=pUe_pt(1)
!	puu_px=pUe_px(1)
!	
!	!dU_dt(i,0:3,1) = pU_pt(0:3,1)+uu*pU_px(i,0:3,1)
!	dU_dt(i,0:3,1) = pU_pt(0:3,1) !+uu*pU_px(i,0:3,1)
!	dU_dt(i,0:3,2) = pU_pt(0:3,2) !+puu_pt*pU_px(i,0:3,1)+2*uu*pU2_pxt(0:3)+uu*puu_px*pU_px(i,0:3,1)+uu**2*pU_px(i,0:3,2)
!
!enddo
!
!end

subroutine material_derivative_try(Ue,pUe,du)
	use global_cont
	use global
	implicit none
	double precision U1(-nv:jx+nv,0:3)
	double precision Ue(-nv:jx+nv,0:3)
	double precision pUe(-nv:jx+nv,0:3,2)
	double precision dU(-nv:jx+nv,0:3,2)
	double precision rho,uu,p,sxx,E,puu,pp,psxx,ei
	double precision p2rho,p2uu,p2p,p2sxx
	double precision drho,duu,dp,dsxx,dE,dei
	double precision dprho,dpuu,dpP,dpsxx,dpE,dpei
	integer i,j,k

	do i = -nv,jx+nv
	rho = ue(i,0)
	uu = ue(i,1)
	p = ue(i,2)
	sxx = ue(i,3)

	call  state_p_to_ei(rho,p ,ei)
	E=ei+uu**2/2
!	call trans_ue_to_u(ue(i,:),u1(i,:))

	prho = pue(i,0,1)
	puu  = pue(i,1,1)
	pp   = pue(i,2,1)
	psxx = pue(i,3,1)

	p2rho = pue(i,0,2)
	p2uu  = pue(i,1,2)
	p2p   = pue(i,2,2)
	p2sxx = pue(i,3,2)


	dU(i,0,1) = -rho*puu
	dU(i,1,1) = -rho*uu*puu-pp+psxx
	dU(i,2,1) = -rho*E*puu-(pp-psxx)*uu-(p-sxx)*puu
	dU(i,3,1) = 4.d0/3*miu*puu


	drho = -rho*puu
	duu  = -1/rho*(pP-psxx)
	dE   = -uu/rho*(pP-psxx)-1/rho*puu*(p-sxx)
	dei  = dE - 2*uu*duu
	dp   = a0**2*f_eta_eta(rho)*drho + rho0*gamma0*dei
	dsxx = 4.d0/3*miu*puu


	dprho = -prho*puu-rho*p2uu
	dpuu  = prho/rho**2*(pP-psxx)-1/rho*(p2P-p2sxx)
	dpE   = -puu/rho*(pP-psxx)+uu/rho**2*prho*(pP-psxx)-uu/rho*(p2P-p2sxx)&
			+1/rho**2*prho*puu*(p-sxx)-1/rho*p2uu*(p-sxx)-1/rho*puu*(pP-psxx)
	dpei  = dpE - 2puu*duu-2*uu*dpuu
	dpP   = a0**2*f_eta_eta_eta(rho)/rho0*prho*drho+a0**2*f_eta_eta(rho)*dprho&
			+rho0*gamma0*dpei
	dpSxx = 4.d0/3*miu*dpuu


	dU(i,0,2) = -drho*dpuu -rho*dpuu
	dU(i,1,2) = -drho*uu*puu-rho*duu*puu-rho*uu*dpuu-dpp+dpsxx
	dU(i,2,2) = -drho*E*puu-rho*dE*puu-rho*E*dpuu- duu*(pP-psxx)&
				-uu*(dpP-dpsxx)-dpuu*(p-sxx)-puu*(dp-dsxx)
	dU(i,3,2) =  4.d0/3*miu*dpuu
	enddo
	end










subroutine trans_ue_A(u,ue,A)
	use global_cont
	implicit none
	double precision A(0:3,0:3)
	double precision u(0:3),ue(0:3)
	double precision rho,uu,p,sigma_x,sxx,ei,p_rho,gamma,E,f_eta_eta

	E=u(2)/u(0)
	rho=ue(0)
	uu=ue(1)
	p=ue(2)
	sxx=ue(3)
	sigma_x=-p+sxx
	gamma=gamma0*rho0/rho
	ei=e-uu*uu/2
	p_rho=a0**2*f_eta_eta(rho)

	A(0,0)=0
	A(1,0)=-uu**2+p_rho+gamma*(uu**2/2-ei)
	A(2,0)=(gamma*(uu**2/2-ei)-ei+sigma_x/rho +p_rho)*uu
	A(3,0)= 4.d0/3*miu*uu/rho

	A(0,1)=1
	A(1,1)=uu*(2-gamma)
	A(2,1)=-gamma*uu**2-sigma_x/rho +ei
	A(3,1)= -4.d0/3*miu/rho

	A(0,2)=0
	A(1,2)=gamma
	A(2,2)=(1+gamma)*uu
	A(3,2)=0

	A(0,3)=0
	A(1,3)=-1
	A(2,3)=-uu
	A(3,3)=uu
	end

subroutine trans_ue_Ae(ue,Ae)
	use global_cont
	implicit none
	double precision Ae(0:3,0:3)
	double precision u(0:3),ue(0:3)
	double precision rho,uu,p,sigma_x,sxx,c

	rho=ue(0)
	uu=ue(1)
	p=ue(2)
	sxx=ue(3)
	call sound(ue,c)

	Ae(0,0) = uu 
	Ae(0,1) = rho
	Ae(0,2) = 0
	Ae(0,3) = 0


	Ae(1,0) = 0 
	Ae(1,1) = uu
	Ae(1,2) = 1.d0/rho
	Ae(1,3) = -1.d0/rho 


	Ae(2,0) = 0 
	Ae(2,1) = rho*c**2-rho0/rho*gamma0*sxx
	Ae(2,2) = uu
	Ae(2,3) = 0

	Ae(3,0) = 0 
	Ae(3,1) = -4.d0/3*miu
	Ae(3,2) = 0
	Ae(3,3) = uu

	end








subroutine trans_ue_dA(u,ue,du,due,dA)
use global_cont
implicit none
double precision dA(0:3,0:3)
double precision ue(0:3),due(0:3),u(0:3),du(0:3)
double precision rho,uu,p,sigma_x,sxx,ei,p_rho,gamma
double precision drho,duu,dp,dsigma,dsxx,dei,dgamma,dE,E
	double precision f_eta_eta,f_eta_eta_eta


	rho=ue(0)
	uu=ue(1)
	p=ue(2)
	sxx=ue(3)


	E=u(2)/u(0)

	dE=du(2)/u(0)-u(2)/u(0)**2*du(0)

	ei=E-0.5*uu**2 !e
	gamma=gamma0*rho0/rho
	sigma_x=-p+sxx
	p_rho=a0**2*f_eta_eta(rho)

	drho=due(0)
	duu=due(1)
	dp=due(2)
	dsxx=due(3)

	dsigma=-dp+dsxx
	dei=dE-uu*duu
	dGamma = -gamma0*rho0*drho/rho**2


	p_rho=a0**2*f_eta_eta(rho)

	dA(0,0)= 0
	dA(0,1) =0
	dA(0,2)= 0
	dA(0,3)=0

	dA(1,0)= -2*uu*duu +a0**2 *f_eta_eta_eta(rho)*drho/rho0+&
			dGamma*(uu*uu/2-ei) +gamma*(uu*duu-dei)
	dA(1,1)=  duu*(2-gamma) - uu*dgamma
	dA(1,2)= dgamma
	dA(1,3)=0

	dA(2,0) = (gamma*(uu*uu/2-ei)-ei+sigma_x/rho+p_rho)*duu +&
			(dgamma*(uu**2/2-ei)+gamma*(uu*duu-dei)-dei+dsigma/rho-&
			sigma_x*drho/rho**2+a0**2*f_eta_eta_eta(rho)*drho/rho0)*uu
	dA(2,1)= -dgamma *uu*uu -2*uu*gamma*duu - dsigma/rho&
			+sigma_x*drho/rho**2+dei
	dA(2,2)=(1+gamma)*duu+dgamma*uu
	dA(2,3)=-duu

	dA(3,0)= 4.d0/3*miu*duu/rho - 4.d0/3*miu*uu*drho/rho**2
	dA(3,1)= 4.d0/3*miu*drho/rho**2
	dA(3,2)= 0
	dA(3,3)= duu

	end






	



