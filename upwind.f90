subroutine upwind(nv,N,LF,RF,H)
implicit none
integer nv ,i,j,k,n
double precision lf(-nv:nv+n)
double precision rf(-nv:nv+n)
double precision h(-nv:nv+n)

do i=-1,n+1
  h(i)=lf(i)+rf(i+1)
  enddo
  
end


