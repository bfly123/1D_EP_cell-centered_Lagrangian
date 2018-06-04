subroutine indicate_KXRCF(ind,nv,n,lf,rf,dx)

implicit none 
double precision P_M,Q_M
double precision u11(3,-nv:n+nv)

double precision dx
integer nv,n,i
integer ind(-nv:n+nv)
 double precision lf(-nv:n+nv),&
          RF(-nv:n+nv),H(-nv:n+nv)
 double precision p_U_H,pp_a,pp_b,pp_c,p_u_abs


      H(:)=lf(:)+rf(:)
      do i=-nv+1,n+nv-1

            u11(1,i)=1.d0/24*(h(i+1)+22*h(i)+h(i-1))
            u11(2,i)=1.d0/2*(h(i+1)-h(i-1))  
            u11(3,i)=1.d0/2*(h(i+1)-2*h(i)+h(i-1)) 
       
        enddo  
 
  do I=-2,N+2 
     
     p_u_h=abs(u11(1,i)-1.d0/2*u11(2,i)+1.d0/6*u11(3,i)&
           -u11(1,i-1)-1.d0/2*u11(2,i-1)-1.d0/6*u11(3,i-1))
      pp_a=u11(1,i)-1.d0/12*u11(3,i)
      pp_b=u11(2,i)
      pp_c=u11(3,i)
      
      p_u_abs=sqrt(pp_a**2+1.d0/12*pp_b**2+1.d0/80*pp_c**2+1.d0/6*pp_a*pp_c)
            
     if(p_u_h/p_u_abs/(dx*0.5)**1.5.gt.1.d0)then

        ind(i)=0 
     endif
    
      p_u_h=abs(u11(1,i)+1.d0/2*u11(2,i)+1.d0/6*u11(3,i)&
            -u11(1,i+1)+1.d0/2*u11(2,i+1)-1.d0/6*u11(3,i+1))
      pp_a=u11(1,i)-1.d0/12*u11(3,i)
      pp_b=u11(2,i)
      pp_c=u11(3,i)
      p_u_abs=sqrt(pp_a**2+1.d0/12*pp_b**2+1.d0/80*pp_c**2+1.d0/6*pp_a*pp_c)
            
        if(p_u_h/p_u_abs/(dx*0.5)**1.5.gt.1.d0)then
           ind(i)=0 
        endif
   enddo
   
  end