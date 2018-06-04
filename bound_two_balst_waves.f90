subroutine bound_two_blast_waves(u4)
 use global
    use init_tran
      use constant
      implicit none
      integer i,j,k
      double precision u4(-nv:jx+nv,0:2) 

!ÓÒ±ß½ç
	DO I=1,nv
        U(jx+i,0)=u(jx-i,0)
        U(jx+i,1)=-u(jx-i,1) 
        U(jx+i,2)=u(jx-i,2)   
	enddo
         
!×ó±ß½ç  
	   DO I=-nv,-1
          U(i,0)=u(-i,0) 
          U(i,1)=-u(-i,1) 
          U(i,2)=u(-i,2) 
          
        enddo 


 
    
	END
	
	
