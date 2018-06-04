subroutine init_shu_osher_problem 
use global
use init_tran
use constant
implicit none
integer i,j,k

    nk=10
      
     
      
    nv=4
    M1=401 
    jx=M1

    dlx=10

    dx=dlx/jx
    

    allocate(U(-nv:jx+nv,0:2))
    allocate(f(-nv:jx+nv,0:2))
    allocate (fp(-nv:jx+nv,0:2))
    allocate (fm(-nv:jx+nv,0:2))
    allocate(x(-nv:jx+nv))
    u=0
    f=0
    fp=0
    fm=0

    SF=0.1
    TT=1.8
    
    
    do i=0,jx
    
    X(i)=i*dx
    rou2=1.+0.3*SIN(5*X(i))
    u2=0
    p2=1.0

    rou1=3.857143
    u1=2.629369
    p1=10.333333
        if(X(i).le.-4)then
            U(i,0)=rou1
            U(i,1)=rou1*u1
            U(i,2)=p1/(GAMA-1)+rou1*(u1*u1)/2
        else
            U(i,0)=rou2
            U(i,1)=rou2*u2
            U(i,2)=p2/(GAMA-1)+rou2*(u2*u2)/2
        endif
    enddo
    
end 