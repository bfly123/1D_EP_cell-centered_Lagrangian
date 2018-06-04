subroutine indicate_XS(ind,nv,n,lf,rf,dx)
implicit none 
double precision P_u_Min,p_u_max,p_fhi,p_b,p_gama
double precision p_a(-nv:n+nv)
double precision ss,dx
integer nv,n,i
integer ind(-nv:n+nv)
 double precision lf(-nv:n+nv),&
          RF(-nv:n+nv),H(-nv:n+nv)


     ss=1.0e-6
     H(:)=lf(:)+rf(:)

     
      do i=-nv+1,n+nv-1
          P_a(i)=(h(i)-h(i-1)) **2+ss
       enddo
       
         p_u_max=0
         p_u_min=0
       do i=-nv+1,n+nv-1
         if(p_u_max<=h(i))then
         
          P_u_max=h(i)
         endif
         if(P_u_min>=h(i))then
         P_U_Min=h(i)
         endif
       enddo
     do I=-1,N+1 
       P_gama=(P_U_MAx-P_u_min)**2/p_a(i)
       p_b=(P_a(i)/p_a(i-1)+p_a(i+1)/p_a(i+2))**2
       
       p_fhi=p_b/(p_b+p_gama) 
       
      if(p_fhi>=dx**2)then
      ind(i+1)=0
      ind(i)=0
      
      endif
     enddo
 end