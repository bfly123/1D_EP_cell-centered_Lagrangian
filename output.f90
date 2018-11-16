    subroutine Output
		use global
    use global_cont
    implicit none
    double precision rou,uu,p,f_eta,a,b,s,e0,error1,error2,rhoE,rho
    integer i,j

		s= 6.d7
		e0= 0.d-14
		a = 1.d4
		b=0.1

        open(1,file='result.dat',status='unknown')
      do 80 i=0,Jx       
        rou=U(i,0)
        uu=U(i,1)/rou
	  p = (u(i,2)/rou- 0.5*uu**2)*rho0*gamma0+rho0*a0**2*f_eta(rou)
        write(1,83)(x(i)+x(i+1))/2,rou,uu,p,U(i,3)
        !write(1,83)(x(i)+x(i+1))/2,u(i,:)
80    continue
      close(1) 

        open(2,file='error.dat',status='unknown')
        open(3,file='exact.dat',status='unknown')

		error1=0
		error2=0
      do i=0,Jx       

	  rho= rho0*(1-b*sin(2*pi*((x(i)+x(i+1))/2-a*0.01))) 

	  rhoE= rho*(e0+ 0.5d0*a**2)

	  write(2,*) x(i),abs(u(i,0)-rho)/rho, abs(u(i,2)-rhoE)/rhoE
	  write(3,83) (x(i)+x(i+1))/2,rho,a,rhoE, s*sin(2*pi*(x(i)-a*0.01))
	  if( abs(u(i,0)-rho)/rho .ge.error1)then
		  error1= abs(u(i,0)-rho)/rho
	  endif
	  if( abs(u(i,2)-rhoE)/rhoE .ge.error2)then
		  error2= abs(u(i,2)-rhoE)/rhoE
	  endif
	  enddo
	  close(2)
	  close(3)
	  write(*,*)  error1,error2

	  


     
83    format(D20.10,D20.10,D20.10,D20.10,D20.10,D20.10)
      end		


