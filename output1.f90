    subroutine Output1(u4)
    use global
    implicit none
    double precision u4(-nv:jx+NV,0:3)
    double precision rou,uu,p
    integer i,j

        open(1,file='result11.dat',status='unknown')
      do 80 i=0,Jx       
        write(1,83)(x(i)+x(i+1))/2,u4(i,0),u4(i,1),u4(i,2),u4(i,3)
80    continue
      close(1) 
     
83    format(D20.10,D20.10,D20.10,D20.10,D20.10)
      end		



