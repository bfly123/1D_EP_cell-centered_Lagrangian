subroutine init_boundreflect
use global
use init_tran
use constant
implicit none
integer i,j,k

    nk=10
      
     
      
    nv=4
    M1=401 
    jx=M1

    dlx=5

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

    SF=0.1d0
    TT=2.0d0
    
    
    do i=0,jx
    
    X(i)=i*dx
    rou2=1.d0
    u2=0.0d0
    p2=1.0d0

    rou1=27.d0/7+dsin(x(i)/2*pi)
    u1=2*dsqrt(1.4d0)/0.9d0
    p1=31.d0/3
    
    !rou1=3.857143
    !u1=2.629369
    !p1=10.333333
        !if(X(i).ge.2*pi)then
        !    U(i,0)=rou2
        !    U(i,1)=rou2*u2
        !    U(i,2)=p2/(GAMA-1)+rou2*(u2*u2)/2
        !    
        !else if(x(i).ge.pi)then
        !    rou1=3.857143
        !    U(i,0)=rou1
        !    U(i,1)=rou1*u1
        !    U(i,2)=p1/(GAMA-1)+rou1*(u1*u1)/2
        !else
        !    
        !    U(i,0)=rou1
        !    U(i,1)=rou1*u1
        !    U(i,2)=p1/(GAMA-1)+rou1*(u1*u1)/2
        !endif
        
        if (x(i).ge.2.d0)then
            
            U(i,0)=rou2
            U(i,1)=rou2*u2
            U(i,2)=p2/(GAMA-1)+rou2*(u2*u2)/2 
        else 
            !rou1=3.857143
            U(i,0)=rou1
            U(i,1)=rou1*u1
            U(i,2)=p1/(GAMA-1)+rou1*(u1*u1)/2
        endif
        
    enddo
    
end 