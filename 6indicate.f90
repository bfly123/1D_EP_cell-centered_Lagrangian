subroutine indicate_BDF(ind,nv,n,lf,rf,dx)

implicit none 
double precision P_M,Q_M
double precision u11(3,-nv:n+nv)

double precision dx,ss
integer nv,n,i
integer ind(-nv:n+nv)
 double precision lf(-nv:n+nv),&
          RF(-nv:n+nv),H(-nv:n+nv)
 double precision p_aaa(3)
 double precision u_mod
      ss=1.e-6

      H(:)=lf(:)+rf(:)

    do i=-nv+1,n+nv-1
       u11(1,i)=1.d0/24*(h(i+1)+22*h(i)+h(i-1))
       u11(2,i)=1.d0/2*(h(i+1)-h(i-1))  
       u11(3,i)=1.d0/2*(h(i+1)-2*h(i)+h(i-1)) 
   enddo  
  
   do I=-2,N+2 
   
      p_aaa(1)=3*u11(3,i)
      p_aaa(2)=u11(2,i+1)-u11(2,i)
      p_aaa(3)=u11(2,i)-u11(2,i-1)
      u_mod=1.d0/3*q_M(p_aaa,3)
      
      
      
          if(abs(u_mod-u11(3,i)).ge.ss)then
             ind(i)=0
          endif
    enddo
    
    end
        