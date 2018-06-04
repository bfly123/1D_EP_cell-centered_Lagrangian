subroutine init_two_blast_waves
use global
use init_tran
use constant
implicit none
integer i,j,k
double precision rou3,u3,p3

    nk=1
      
      
      
    nv=4
    M1=1000 
    jx=M1

    dlx=1.

    dx=dlx/jx


    allocate(U(-nv:jx+nv,0:2))
    allocate(f(-nv:jx+nv,0:2))
    allocate (fp(-nv:jx+nv,0:2))
    allocate (fm(-nv:jx+nv,0:2))
    
    u=0
    f=0
    fp=0
    fm=0

     SF=0.05
    TT=0.038

    rou1=1.
    u1=0.
    p1=1000

    rou2=1.
    u2=0.0
    p2=100
    
    
   rou3=1.
   u3=0.
   p3=0.01
    
    do i=0,jx
        if(i*dx.ge.0.9)then
            U(i,0)=rou2
            U(i,1)=rou2*u2
            U(i,2)=p2/(GAMA-1)+rou2*(u2*u2)/2
        else if(i*dx.le.0.1)then
            U(i,0)=rou1
            U(i,1)=rou1*u1
            U(i,2)=p1/(GAMA-1)+rou1*(u1*u1)/2
        else
            U(i,0)=rou3
            U(i,1)=rou3*u3
            U(i,2)=p3/(GAMA-1)+rou3*(u3*u3)/2
        endif
    enddo
    
end 