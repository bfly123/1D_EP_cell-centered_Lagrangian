subroutine bound_lax_problem(u4)
 use global
    use init_tran
      use constant
      implicit none
      integer i,j,k
      double precision u4(-nv:jx+nv,0:2) 

!�ұ߽�
	DO I=1,nv
        U(jx+i,0)=rou2 
        U(jx+i,1)=rou2*u2 
        U(jx+i,2)=p2/(GAMA-1)+rou2*(u2*u2)/2   
	enddo
         
!��߽�  
	   DO I=-nv,-1
          U(i,0)=rou1 
          U(i,1)=rou1*u1 
          U(i,2)=p1/(GAMA-1)+rou1*(u1*u1)/2 
          
        enddo 


 
    
	END
	
	
