subroutine indicate_shen(ind,nv,n,lf,rf,dx)

implicit none 

double precision p_a(3)
double precision u0,u2,u01,u21,dx
integer nv,n,i
integer ind(-nv:n+nv)
 double precision lf(-nv:n+nv),&
          RF(-nv:n+nv),H(-nv:n+nv)
double precision is0,is1,is2, tao5,dmin
double precision ss
   H(:)=lf(:)+rf(:)
   
   do I=-1,N
          IS0=13.d0*(h(I-1)-2.*h(I)+h(I+1))**2/12.d0&
                          +(h(I-1)-4.*h(I)+3.*h(I+1))**2/4.d0
          IS1=13.d0*(h(I-1)-2.*h(I)+h(I+1))**2/12.d0&
                          +(h(I-1)-h(I+1))**2/4.d0 
          IS2=13.d0*(h(I+1)-2.*h(I+2)+h(I+3))**2/12.d0&
                          +(3.*h(I+1)-4.*h(I+2)+h(I+3))**2/4.d0
          tao5=abs(is0-is2)
          dmin=dmin1(is0,is1,is2)
      if (tao5.gt.dmin) then  
          ind(i-1)=0
          ind(i)=0
        endif
enddo


end
      