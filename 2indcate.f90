!TVB Ê¶±ğÆ÷
subroutine indicate_TVB(ind,nv,n,lf,rf,dx)
implicit none 
double precision P_M
double precision u11(3,-nv:n+nv)
double precision p_a(3)
double precision u0,u2,u01,u21,dx
integer nv,n,i
integer ind(-nv:n+nv)
 double precision lf(-nv:n+nv),RF(-nv:n+nv),h(-nv:n+nv)
       
        H(:)=lf(:)+rf(:)

      do i=-nv+1,n+nv-1

            u11(1,i)=1.d0/24*(h(i+1)+22*h(i)+h(i-1))
            u11(2,i)=1.d0/2*(h(i+1)-h(i-1))  
            u11(3,i)=1.d0/2*(h(i+1)-2*h(i)+h(i-1)) 
       
        enddo  
        p_M=0  
        do i=-2,n+2
           if(P_M.le.abs(u11(1,i+1)+u11(1,i-1)-2*u11(1,i)))then 
              P_M=abs(u11(1,i+1)+u11(1,i-1)-2*u11(1,i)) 
           endif
       enddo     
        P_M=2.d0/3*P_M   
    
    
     do I=-2,N+2 
           u0=1.d0/2*u11(2,i)+u11(3,i)*1.d0/6
           u2=1.d0/2*u11(2,i+1)-1.d0/6*u11(3,i+1)
          
            p_a(1)=u0
            p_a(2)=u11(1,i+1)-u11(1,i)
            p_a(3)=u11(1,i)-u11(1,i-1)  
        
         if(abs(u0).lt.p_M*dx**2)then
            u01=u0
         elseif(p_a(1)*p_a(2)>0.and.p_a(1)*p_a(3)>0)then
         u01=dmin1(abs(p_a(1)),abs(p_a(2)),&
             abs(p_a(3)))*p_a(1)/abs(p_a(1))
!            u01=q_m(p_a,3)
         else 
         u01=0.
         endif
             
            p_a(1)=u2
            p_a(2)=u11(1,i+1)-u11(1,i)
            p_a(3)=u11(1,i)-u11(1,i-1)
                    
           if(abs(u2).lt.p_M*dx**2)then
            u21=u2
          elseif(p_a(1)*p_a(2)>0.and.p_a(1)*p_a(3)>0)then
         u21=dmin1(abs(p_a(1)),abs(p_a(2)),&
             abs(p_a(3)))*p_a(1)/abs(p_a(1)) 
            
           ! u21=q_m(p_a,3)
           else
           u21=0.
         endif

     if ((u01.ne.u0).or.(u21.ne.u2))then 
         
         ind(i)=0 
     endif
  enddo
  end