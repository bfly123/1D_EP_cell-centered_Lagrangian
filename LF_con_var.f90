subroutine LF_con_var(u4,h)
    use global
    use constant
    implicit none
    
    integer i,j,k,ij
    double precision u4(-nv:jx+nv,0:2)
    double precision h(-nv:jx+nv,0:2)
    double precision F_L(0:2)
    double precision F_R(0:2)
    double precision uL(-nv:jx+nv,0:2) 
    double precision uR(-nv:jx+nv,0:2)
    double precision cc,rho,uu,p,uv,a,tmp1
    double precision rho_L,rho_R,uu_L,uu_R,p_L,p_R,c_L,c_R
    double precision u_roe,c_roe,S_L,s_R,S_star,s_min,s_max,Chi_L,Chi_R
    double precision U_Star_L(0:2)
    double precision U_star_R(0:2)
    
    
     tmp1=real(1.d0/(gama-1.d0))
    
    
  !  do i=-nv,jx+nv
  !          rho=U4(i,0)
  !          uu=U4(i,1)/rho
  !          
  !          z(n_spec)=1
  !         do ij=1,n_spec-1
  !             z(ij)=u(i,ij+2)
  !             z(n_spec)=Z(n_spec)-Z(ij)
  !          enddo
  !      
  !          a=0
  !          do ij=1,n_spec
  !             a=a+q0(ij)*z(ij)  
  !          enddo
  !          
  !          p=(GAMA-1.d0)*(U4(i,2)-0.5d0*U4(i,0)*(uu*uu)-rho*a)    
  !          u_m(i,0)=rho
  !          u_m(i,1)=uu
  !          u_m(i,2)=p
  !          
  !      enddo
  !!      do i=-1,jx+1
  !!      do j=-1,jy+1
  !!          f(i,j,:)=u_m(i,j,:)-u_m(i,j-1,:)
  !!      enddo
  !!      enddo
  !!call Output1(f)
  
!Çóur ul h f1
        
        do k=0,2  
           call weno5_new(nv,jx,u4(:,k),uL(:,k),uR(:,k)) 
           !call upwind(nv,jx,fp(:,j,k),fm(:,j,k),h(:,j,k))  
    !      call weno5(nv,jx,fp(:,j,k),fm(:,j,k),h(:,j,k))
        enddo
       do i=-1,jx+1
            rho_L=uL(i,0)
            uu_L=uL(i,1)/rho_L 
            p_L=(GAMA-1.d0)*(UL(i,2)-0.5d0*rho_L*(uu_L*uu_L))
            
            F_L(0)=rho_L*uu_L
            F_L(1)=rho_L*uu_L*uu_L+p_L
            F_L(2)=uu_L*(uL(i,2)+p_L)
            C_L=sqrt(gama*p_L/rho_L)
            
             
            rho_R=uR(i,0)
            uu_R=uR(i,1)/rho_R 
            p_R=(GAMA-1.d0)*(UR(i,2)-0.5d0*rho_R*(uu_R*uu_R))
            
            F_R(0)=rho_R*uu_R
            F_R(1)=rho_R*uu_R*uu_R+p_R
            F_R(2)=uu_R*(uR(i,2)+p_R)
            C_R=sqrt(gama*p_R/rho_R)

           !!LF 
            h(i,:)=0.5d0*(F_L(:)+F_R(:)+max(0.d0,(c_L+abs(uu_L)),(c_R+abs(uu_R)))*(uL(i,:)-UR(i,:)))!max((c_L+abs(uu_L)),(c_R+abs(uu_R)))*(u_L(:)-U_R(:)))
           !*************************!HLLC
           !u_roe=(sqrt(rho_L)*uu_L+sqrt(rho_R)*uu_R)/(sqrt(rho_L)+sqrt(rho_R))
           !c_roe=(sqrt(rho_L)*c_L+sqrt(rho_R)*c_R)/(sqrt(rho_L)+sqrt(rho_R))
           !S_L=min((u_roe-c_roe),uu_L-c_L)
           !S_R=max((u_roe+c_roe),uu_R+c_R)
           !S_star=(p_R-p_L+rho_L*uu_L*(s_L-uu_L)-rho_R*uu_R*(S_R-uu_R))/(rho_L*(S_L-uu_L)-rho_R*(S_R-uu_R))
           !S_min=min(0.d0,s_L)
           !S_max=max(0.d0,S_R)
           !Chi_L=(s_L-uu_L)/(s_L-s_star)
           !chi_R=(s_R-uu_R)/(s_R-s_star)
           !
           !u_star_L(0)=Chi_L*rho_L
           !u_star_L(1)=Chi_L*rho_L*S_star
           !u_star_L(2)=Chi_L*(uL(i,2)+(s_star-uu_L)*(rho_L*S_star+p_L/(s_L-uu_L)))
           !
           !u_star_R(0)=Chi_R*rho_R
           !u_star_R(1)=Chi_R*rho_R*S_star
           !u_star_R(2)=Chi_R*(uR(i,2)+(s_star-uu_R)*(rho_R*S_star+p_R/(s_R-uu_R)))
           !h(i,:)=(1+sign(1.0d0,s_star))/2*(F_L(:)+s_min*(u_star_L(:)-uL(i,:)))&
           !      +(1-sign(1.0d0,s_star))/2*(f_R(:)+s_max*(u_star_R(:)-uR(i,:)))
       
    enddo
    
end