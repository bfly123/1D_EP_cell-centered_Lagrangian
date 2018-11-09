subroutine subcell_WENO3(nv,jx,dx,u,uL,uR,pUL,pUR)
	  !use global_cont
	  implicit none
	  integer nv,jx,i
	  double precision u(-nv:jx+nv)
	  double precision uL(-nv:jx+nv)
	  double precision uR(-nv:jx+nv)
	  double precision puL(-nv:jx+nv,3)
	  double precision puR(-nv:jx+nv,3)
	  double precision  IS0,IS1,c0,c1,a0,a1,b0,b1,u1,ss,u2,d0,d1,d2,dx
	
	  ss=1.d-8

	  !***UL
	 do i=-nv+1,jx+nv-1

	  IS0= (U(i)-U(i-1))**2
	  IS1= (U(i+1)-u(i))**2

	  c0= 2.d0/3
	  c1= 1.d0/3
	  a0= c0/(ss+IS0)**2
	  a1= c1/(ss+IS1)**2

	  b0=a0/(a0+a1)
	  b1=a1/(a0+a1)

	  u1= b0*(u(i)+u(i-1))/2+b1*(3*u(i)-u(i+1))/2

	  c0= 1.d0/3
	  c1= 2.d0/3
	  a0= c0/(ss+IS0)**2
	  a1= c1/(ss+IS1)**2

	  b0=a0/(a0+a1)
	  b1=a1/(a0+a1)

	  u2= b0*(3*u(i)-u(i-1))/2+b1*(u(i)+u(i+1))/2

	  d0=u1
	  d1= 6*u(i)-2*u2-4*u1
	  d2=3*(u1+u2-2*u(i))

	  uL(i)= u2
	  puL(i,1)= (d1+2*d2)/dx
	  puL(i,2)= 2*d2/dx**2  

	  uR(i-1)=u1
	  puR(i-1,1) = d1/dx
	  puR(i-1,2) = 2*d2/dx**2

	  enddo
end
