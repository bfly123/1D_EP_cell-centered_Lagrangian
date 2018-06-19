subroutine upwind(nv,N,u,ul,ur)
implicit none
integer nv ,i,j,k,n
double precision ul(-nv:nv+n)
double precision ur(-nv:nv+n)
double precision u(-nv:nv+n)

do i=-nv,n+nv+1
ul(i)= u(i)
ur(i)= u(i+1)

  enddo
  
end


