subroutine indicate_MMP(ind,nv,n,lf,rf,dx)
implicit none
double precision u11(3,-nv:n+nv)

double precision dx,ss
integer nv,n,i
integer ind(-nv:n+nv)
 double precision lf(-nv:n+nv),&
          RF(-nv:n+nv),H(-nv:n+nv)
 double precision u_mod,u_p_half,u_m_half
 double precision du_min,dmin_u,p_xy,p_phi
  
      ss=1.e-10

     H(:)=lf(:)+rf(:)

    do i=-nv+1,n+nv-1
       u11(1,i)=1.d0/24*(h(i+1)+22*h(i)+h(i-1))
       u11(2,i)=1.d0/2*(h(i+1)-h(i-1))  
       u11(3,i)=1.d0/2*(h(i+1)-2*h(i)+h(i-1)) 
   enddo  

do i=-2,n+2
  u_p_half=u11(1,i)-6*u11(2,i)+30*u11(3,i)
  u_m_half=u11(1,i)+6*u11(2,i)+30*u11(3,i)
  du_min=u11(1,i)-dmin1(u11(1,i-1),u11(1,i),u11(1,i+1))
  dmin_u=u11(1,i)-dmin1(u_p_half,u_m_half)
  p_xy=du_min/(dmin_u+1.e-12)
  p_phi=dmin1(1.,p_xy)
  
  if(abs(p_phi-1).ge.1.e-12)then
    ind(i)=0
   endif
 enddo
 
 end