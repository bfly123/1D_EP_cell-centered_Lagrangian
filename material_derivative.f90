subroutine material_derivative(U1,pU_px,du_dt)
use global_cont
use global
implicit none
double precision U1(-nv:jx+nv,0:3)
double precision Ue(0:3)
double precision dU_dt(-nv:jx+nv,0:3,3)
double precision pU_px(-nv:jx+nv,0:3,3)
double precision pUe_px(0:3)
double precision pUe_pt(0:3)
double precision pU_pt(0:3,3)
double precision pU2_pxt(0:3)
double precision A(0:3,0:3)
double precision A_pA_px(0:3,0:3)
double precision A2(0:3,0:3)
double precision pA_px(0:3,0:3)
double precision pA_pt(0:3,0:3)
double precision puu_pt,puu_px,uu,f_eta_eta
integer i,j,k


!call trans_px_to_pt(pu_px,pu_pt)
!call time_derivative(U_half,pU_pt,pU_px)
do i=-nv,jx+nv

	call trans_u_to_Ue(U1(i,0:3),ue)
	call trans_ue_A(U1(i,0:3),ue,A)

!*pU_pt=- A*pU_px
	Pu_pt(0:3,1)=- matmul(A,pU_px(i,0:3,1))

	call  trans_du_to_due(U1(i,0:3),pu_pt(0:3,1),pue_pt)  !**********\partial U/\parttial t to \partial U_e/\partial t
	call  trans_du_to_due(U1(i,0:3),pu_px(i,0:3,1),pue_px) 
	
	call trans_ue_dA(U1(i,0:3),Ue,pUe_pt,pA_pt)
	call trans_ue_dA(U1(i,0:3),Ue,pUe_px,pA_px)
	
	A_pA_px=matmul(A,pA_px)

	A2=matmul(A,A)

	pU2_pxt(0:3)=- matmul(pA_px,pU_px(i,0:3,1))- matmul(A,pU_px(i,0:3,2))
	pU_pt(0:3,2)=- matmul(pA_pt,pU_px(i,0:3,1))- matmul(A,pU2_pxt)

	uu=Ue(1)

	puu_pt=pUe_pt(1)
	puu_px=pUe_px(1)
	
	dU_dt(i,0:3,1) = pU_pt(0:3,1)+uu*pU_px(i,0:3,1)
	dU_dt(i,0:3,2) = pU_pt(0:3,2)+puu_pt*pU_px(i,0:3,1)+2*uu*pU2_pxt(0:3)+uu*puu_px*pU_px(i,0:3,1)+uu**2*pU_px(i,0:3,2)

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


subroutine trans_ue_dA(u,ue,due,dA)
use global_cont
implicit none
double precision dA(0:3,0:3)
double precision ue(0:3),due(0:3),u(0:3)
double precision rho,uu,p,sigma_x,sxx,ei,p_rho,gamma
double precision drho,duu,dp,dsigma,dsxx,dei,dgamma,dE,E
	double precision f_eta_eta,f_eta_eta_eta


	rho=ue(0)
	uu=ue(1)
	p=ue(2)
	sxx=ue(3)


	E=u(2)/u(0)

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






	



