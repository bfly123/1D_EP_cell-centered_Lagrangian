    subroutine Output
    use global
    use constant
    implicit none
    double precision rou,uu,p
    integer i,j

        open(1,file='result.dat',status='unknown')
      do 80 i=0,Jx       
        rou=U(i,0)
        uu=U(i,1)/rou
        p=(GAMA-1)*(U(i,2)-0.5*U(i,0)*(uu*uu))
        write(1,83)dx*i,rou,uu,p
80    continue
      close(1) 
     
83    format(D20.10,D20.10,D20.10,D20.10,D20.10,D20.10)
      end		



