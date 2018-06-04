subroutine init_lax_problem
use global
use init_tran
use constant
implicit none
integer i,j,k

    nk=10
      
      
      
    nv=4
    M1=400
    jx=M1

    dlx=10

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
    TT=1.3

    rou1=0.445
    u1=0.698
    p1=3.528

    rou2=0.5
    u2=0.0
    p2=0.571
   
    do i=0,jx
        if(i*dx.ge.5)then
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