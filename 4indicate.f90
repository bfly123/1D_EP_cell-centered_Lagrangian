 subroutine indicate_MP(ind,nv,n,lf,rf,dx)
 implicit none 
double precision Q_M
double precision u11(3,-nv:n+nv)
double precision dx
integer nv,n,i
integer ind(-nv:n+nv)
 double precision lf(-nv:n+nv),&
          RF(-nv:n+nv),H(-nv:n+nv)
double precision p_aa(6)
double precision p_d_m4x(-nv:n+nv)
double precision  p_d(-nv:n+nv)
double precision p_u_MD,p_u_UL, p_u_LC,p_u_half, l_u_min,l_u_max,q_median
double precision ss,u_p,u_m,m_u_half
        
        ss=1.d-20
        H(:)=lf(:)+rf(:)

        do i=-nv+1,n+nv-1
            u11(1,i)=1.d0/24*(h(i+1)+22*h(i)+h(i-1))
            u11(2,i)=1.d0/2*(h(i+1)-h(i-1))  
            u11(3,i)=1.d0/2*(h(i+1)-2*h(i)+h(i-1)) 
       enddo  
     
     do i=-2,n+2
        p_d(i)=u11(1,i+1)-2*u11(1,i)+u11(1,i-1)
      enddo
      do i=-1,n+1
         
          p_aa(1)=4*p_d(i)-p_d(i+1)
          p_aa(2)=4*p_d(i+1)-p_d(i)
          p_aa(3)=p_d(i)
          p_aa(4)=p_d(i+1)
          p_aa(5)=p_d(i-1)
          p_aa(6)=p_d(i+2)
          
          p_d_m4x(i)=q_m(p_aa,6)
      enddo
    do i=-2,n+2 
       u_p=u11(1,i+1)-6*u11(2,i+1)+30*u11(3,i+1)
       u_m=u11(1,i)+6*u11(2,i)+30*u11(3,i)   
       p_u_MD=1.d0/2*(u11(1,i)+u11(1,i+1)-p_d_M4X(i)) 
       p_u_UL=u11(1,i)+2*(u11(1,i)-u11(1,i-1))
       p_u_LC=u11(1,i)+1.d0/2*(u11(1,i)-u11(1,i-1))+4.d0/3*p_d_m4x(i-1)
       
       l_u_min=dmax1(dmin1(u11(1,i),u11(1,i+1),p_u_MD)&
                 ,dmin1(u11(1,i),p_u_UL,P_u_LC))
       l_u_max=dmin1(dmax1(u11(1,i),u11(1,i+1),p_u_MD)&
                 ,dmax1(u11(1,i),p_u_UL,P_u_LC))
     
        p_u_half=q_median(u_p,L_U_min,L_U_max) 
        m_u_half=q_median(u_m,L_U_min,L_U_max) 
 
          if(u_p.ne.p_u_half.or.u_m.ne.m_u_half)then 
             ind(i)=0
           endif
           
        enddo

 end