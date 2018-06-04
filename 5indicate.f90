subroutine indicate_MR(ind,nv,n,lf,rf,dx)
implicit none 
double precision dx
integer nv,n,i
integer ind(-nv:n+nv)
double precision lf(-nv:n+nv),&
          RF(-nv:n+nv),H(-nv:n+nv)
          
double precision d

h=lf+rf

do I=-2,N+2
   d=h(i)-1.d0/2*(h(i-1)+h(i+1))
   if(abs(d).ge.0.25*dx)then
     ind(i)=0
    endif
enddo

end
