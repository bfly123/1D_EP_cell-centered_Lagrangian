SUBROUTINE WENO(nv,jx,fp,fm,h,ind)
implicit none
integer nv,jx,il,iu,n,i
double precision fp(-nv:jx+nv),&
          fm(-nv:jx+nv),h(-nv:jx+nv)
 double precision lf(-nv+1:jx+nv+1),&
          RF(-nv+1:jx+nv+1),H1(-nv+1:jx+nv+1)&
            ,HR(-nv+1:jx+nv+1)
integer ind(-nv:jx+nv),ind1(-nv+1:jx+nv+1)
il= -nv+1
iu=jx+nv+1
n=jx+1
do i=-nv+1,jx+nv+1
   lf(i)=fp(i-1)
   rf(i)=fm(i-1)
   ind1(i)=ind(i-1)
enddo                    
call WENO_NEW(IL,IU,N,LF,RF,H1,HR,ind1)
! call WENONEW(IL,IU,N,LF,RF,H1,HR)
 
do i=-nv,jx+nv
    h(i)=h1(i+1)+hr(i+1)
enddo


end


