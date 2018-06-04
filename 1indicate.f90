subroutine indicate_ATV(ind,nv,n,lf,rf,dx)
implicit none 

double precision dx,TV
integer nv,n,i
integer ind(-nv:n+nv)
 double precision lf(-nv:n+nv),&
          RF(-nv:n+nv),H(-nv:n+nv)
 
 
   
   H(:)=lf(:)+rf(:)
   TV=0
  do i=-nv+1,n+nv-1
     
  TV=TV+abs(h(i)-h(i-1))
  enddo
  TV=TV/(n+2*nv-3)
   
 do I=-2,N+2 
   if(abs(h(i)-h(i-1)).gt.0.3*TV)then
   ind(i)=0
   ind(i+1)=0
   endif
 enddo
 
 end
   