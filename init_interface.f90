subroutine init_interface
use global
use init_tran
use constant
implicit none
integer i,j,k

    nk=1
      
      
      
    nv=4
    M1=401 
    jx=M1

    dlx=2

    dx=dlx/jx


    allocate(U(-nv:jx+nv,0:2))
    allocate(f(-nv:jx+nv,0:2))
    allocate (fp(-nv:jx+nv,0:2))
    allocate (fm(-nv:jx+nv,0:2))
    
    u=0
    f=0
    fp=0
    fm=0

    SF=0.1
    TT=0.28

    rou2=1.0
    u2=1
    p2=1.0/1.4
    
    rou1=10
    u1=1
    p1=1.0/1.4
    
    do i=0,jx
        if(i.ge.int(jx/2))then
            U(i,0)=rou2
            U(i,1)=rou2*u2
            U(i,2)=p2/(GAMA-1)+rou2*(u2*u2)/2
        else
            U(i,0)=rou1
            U(i,1)=rou1*u1
            U(i,2)=p1/(GAMA-1)+rou1*(u1*u1)/2
        endif
    enddo
    
end 
