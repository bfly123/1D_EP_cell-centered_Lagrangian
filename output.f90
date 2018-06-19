    subroutine Output
		use global
    use global_cont
    implicit none
    double precision rou,uu,p,f_eta
    integer i,j

        open(1,file='result.dat',status='unknown')
      do 80 i=0,Jx       
        rou=U(i,0)
        uu=U(i,1)/rou
	  p = (u(i,2)/rou- 0.5*uu**2)*rho0*gamma0+rho0*a0**2*f_eta(rou)
        write(1,83)(x(i)+x(i+1))/2,rou,uu,p,U(i,3)
80    continue
      close(1) 
     
83    format(D20.10,D20.10,D20.10,D20.10,D20.10,D20.10)
      end		


