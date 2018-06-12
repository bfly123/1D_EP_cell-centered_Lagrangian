subroutine init_sod_problem
use init_tran
use global_cont
use global
implicit none
integer i,j,k

    nv=4
    jx=401 
    dlx=2
    dx=dlx/jx

    allocate(U(-nv:jx+nv,0:3))
    allocate(X(-nv:jx+nv))
    
    u=0
    f=0

    TT=0.28

    rou2=0.125
    u2=0
    p2=0.1
    
    rou1=1.0
    u1=0
    p1=1.0
    
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
