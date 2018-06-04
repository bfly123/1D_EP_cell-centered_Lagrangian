subroutine bound_boundreflect(u4)
 use global
    use init_tran
      use constant
      implicit none
      integer i,j,k
      double precision u4(-nv:jx+nv,0:2) 

!ÓÒ±ß½ç
	DO I=1,nv
        U(jx+i,0)=U(jx-i,0) 
        U(jx+i,1)=-u(jx-i,1) 
        U(jx+i,2)=u(jx-i,2)
        !U(jx,1)=0
    enddo
    
    rou1=27.d0/7
    u1=2*dsqrt(1.4d0)/0.9d0
    p1=31.d0/3     
!×ó±ß½ç  
	   DO I=-nv,-1
          U(i,0)=rou1 
          U(i,1)=rou1*u1 
          U(i,2)=p1/(GAMA-1)+rou1*(u1*u1)/2 
          
        enddo 


 
    
	END
	
	
