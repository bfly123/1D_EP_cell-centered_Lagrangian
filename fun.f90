      function sing(x)
      implicit none 
      double precision x,sing,n
    
      !double precision x,y,sing
      if( x>0)then 
      n=1
      else if(x<0)then 
      n=-1
      else if(x==0)then
      n=0
      endif
      sing=n
     ! return 
      end 
      
      
      function q_median(x,y,z)
      IMPLICIT DOUBLE PRECISION(a-h,l,o-z)
    !  double precision x,y,z
      dimension a(2)
      integer sing
      ! double precision q_median
      a(1)=y-x
      a(2)=z-x
      
      q_median=x+q_M(a,2)
      
      end 
      
      function q_M(a,k)
      IMPLICIT none
      integer k,j,i
      double precision  a(k)
      double precision q_M,m
      double precision sing
     
      
      m=sing(a(1))
      do i=1,k
      if(sing(a(i)).ne.m)then
         m=0
      endif 
      enddo
    
      
      q_m=abs(a(1)) 
       
       do i=1,k
          if(q_m.ge.abs(a(i)))then
             q_m=abs(a(i))
           endif
      enddo
      
      q_m=q_m*m
      
      end
       
       