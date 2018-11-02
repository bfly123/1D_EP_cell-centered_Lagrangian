subroutine material_derivative(U_half,dU_dt,pU_px)
use global_cont
use global
implicit none
double precision U_half(-nv:jx+nv,0:3)
double precision dU_dt(-nv:jx+nv,0:3,3)
double precision pU_px(-nv:jx+nv,0:3,3)
double precision pUe_px(-nv:jx+nv,0:3,3)
double precision pU_pt(-nv:jx+nv,0:3,3)
double precision pU2_pxt(-nv:jx+nv,0:3,3)
double precision pUe_pt(-nv:jx+nv,0:3,3)


call subcell(u,pu_px)
!call trans_px_to_pt(pu_px,pu_pt)
!call time_derivative(U_half,pU_pt,pU_px)

call  trans_du_to_due(u_half,pu_pt,pue_pt) 
call  trans_du_to_due(u_half,pu_px,pue_px) 

do i=-nv,jx+nv
call trans_ue_A(u,ue,A)
call trans_ue_dA(u,ue,pue_pt,pA_pt)
call trans_ue_dA(u,ue,pue_px,pA_px)
	do j=0,3
		pu_pt(i,j,1)=0
		do k=0,3
			pU_pt(i,j,1)=pU_pt(i,j,1)- A(j,k)*Pu_px(i,k,1)
		enddo
		pu_pt(i,j,2)=0
		A_PA_px=A*pA_px
		A2=A*A
		pU2_pxt=-PA_px*pU_px-A*pU_px(2)
		do k=0,3
			pU_pt(i,j,2)=pU_pt(i,j,2)- pA_pt(j,k)*Pu_px(i,k,1)+A_pA_px(j,k)*pu_px+A2(j,k)*pU_px(i,k,2)
		enddo
	enddo

do i=-nv,jx+nv
	uu=u(i,1)/u(i,0)
	duu_dt=pU_pt(i,1)/pU_pt(i,0)
	duu_dx=pU_px(i,1)/pU_px(i,0)
	do j =0,3
		du_dt(i,j,1) = pu_pt(i,j,1)+uu*pu_px(i,j,1)
		du_dt(i,j,2) = pu_pt(i,j,2)+2*uu*pu_px(i,j,1)
	enddo
enddo




call _solve_pA()

do i=-nv,jx+nv
	do j =0,3
		uu=u_half(i,1)/u_half(i,0)
		du_dt(i,j,2) = pu_pt(i,j,2)+pue_pt(i,1)*pu_px(i,j,1)&
	+2*uu*p2U_pxpt+uu*pu_pt(i,j,1)+uu*pue_px(i,1)*pU_px(i,j,1)
	+uu**2*pu_px(i,j,2)
	enddo
enddo



do i=-nv,jx+nv
	do j =0,3
		uu=u_half(i,1)/u_half(i,0)
		du_dt(i,j,2) = pu_pt(i,j,2)+u_half*pu_px(i,j,1)
	enddo
enddo


end


	subroutine trans_px_to_pt(u,pU_px,pu_pt)
	use global_cont
	implicit none
	double precision pu_px(0:3,3)
	double precision u(0:3)
	double precision ue(0:3)
	double precision pu_pt(0:3,3)
	double precision A(0:3,0:3)
	double precision pA_pt(0:3,0:3)
	double precision pA_px(0:3,0:3)

	call trans_u_to_ue(u,ue) 
	call trans_ue_A(u(3))


end


subroutine trans_ue_A(u,ue,A)
	use global_cont
	implicit none
	double precision A(0:3,0:3)
	double precision u(0:3),ue(0:3)
	double precision rho,uu,p,sigma_x,sxx,ei,p_rho,gamma,E

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
double precision drho,duu,dp,dsigma,dsxx,dei,dgamma
	double precision f_eta_eta,f_eta_eta_eta


	rho=ue(0)
	uu=ue(1)
	p=ue(2)
	sxx=ue(3)


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

	dA(2,0) = (gamma*(uu*uu/2-ei)-ei+sigma_x/rho+p_rho)*du +&
			(dgamma*(uu**2/2-ei)+gamma*(uu*duu-dei)-dei+dsigma/rho-&
			sigma_x*drho/rho**2+a0**2*f_eta_eta_eta(rho)*drho/rho0)*uu
	dA(2,1)= -dgamma *uu*uu -2*uu*gamma - dsigma/rho&
			+sigma_x*drho/rho**2+dei
	dA(2,2)=(1+gamma)*duu+dgamma*uu
	dA(2,3)=-duu

	dA(3,0)= 4.d0/3*miu*duu/rho - 4.d0/3*miu*uu*drho/rho**2
	dA(3,1)= 4.d0/3*miu*drho/rho**2
	dA(3,2)= 0
	dA(3,3)= duu

	end







	



