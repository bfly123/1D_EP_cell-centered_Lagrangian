subroutine weno3_new(nv,jx,F,HL,HR)
 IMPLICIT none
 double precision ss,a300,a301,a302,a310,a311,a312,a320,&
                   a321,a322,c30,c31,c32,is0,is1,is2,tao5,&
                   q30,q31,q32,aa0,aa1,aa2,w0,w1,w2
 integer  n,nv,jx,i
 double precision hl(-nv:nv+jx),hr(-nv:nv+jx),f(-nv:nv+jx)
         SS=1E-20
        
        do i=-1,jx+1
             
             Q30=f(i)-f(i-1)
             Q31=f(i+1)-f(i)
      
!!      
                w0=1.d0/(1.d0+2.d0*(q30*q30+ss)/(q31*q31+ss))
                hl(i)= f(i)+(w0*(q30-q31)+q31)/2.d0
         
            Q30=f(i+1)-f(i+2)
             Q31=-Q31
            
                w0=1.d0/(1.d0+2.d0*(q30*q30+ss)/(q31*q31+ss))

                hr(i)= f(i+1)+(w0*(q30-q31)+q31)/2.d0

    enddo
end

subroutine weno3_origin(nv,jx,F,HL,HR)
 IMPLICIT none
 double precision ss,a300,a301,a302,a310,a311,a312,a320,&
                   a321,a322,c30,c31,c32,is0,is1,is2,tao5,&
                   q30,q31,q32,aa0,aa1,aa2,w0,w1,w2,b1,b2,w3,w4
 integer  n,nv,jx,i
 double precision hl(-nv:nv+jx),hr(-nv:nv+jx),f(-nv:nv+jx)
         SS=1E-8
        
        do i=-1,jx+1
             
             Q30=-1.d0/2*f(i-1)+3.d0/2*f(i)
             Q31=1.d0/2*f(i+1)+1.d0/2*f(i)
      
			 b1=(f(i)-f(i-1))**2
			 b2=(f(i+1)-f(i))**2

			 w3=1.d0/3/(ss+b1)
			 w4=2.d0/3/(ss+b2)

			 w1=w3/(w3+w4)
			 w2=w4/(w3+w4)
!!      
                hl(i)= w1*q30+w2*q31
         
            
             Q30=-1.d0/2*f(i+2)+3.d0/2*f(i+1)
             Q31=1.d0/2*f(i+1)+1.d0/2*f(i)

			 b1=(f(i+1)-f(i+2))**2
			 b2=(f(i+1)-f(i))**2

			 w3=1.d0/3/(ss+b1)
			 w4=2.d0/3/(ss+b2)

			 w1=w3/(w3+w4)
			 w2=w4/(w3+w4)

                hr(i)= w1*q30+w2*q31
    enddo
end
subroutine weno3LIU_new(nv,jx,F,HL,HR)
 IMPLICIT none
 double precision ss, w,q1,q2,b1,b2
 integer  n,nv,jx,i
 double precision hl(-nv:nv+jx),hr(-nv:nv+jx),f(-nv:nv+jx)
         SS=1E-6
        
        do i=-1,jx+1
             
             q1=f(i)-f(i-1)
             q2=f(i+1)-f(i)

    		 b1=(abs(q1+q2)-abs(3.2*q1-q2))**2+ss
    		 b2=(abs(q1+q2)-abs(3.2*q2-q1))**2+ss
		!	b1=q1**2+ss
		!	b2=q2**2+ss
!!      
                w=1.d0/(1.d0+2.d0*(b1/b2))
                hl(i)= f(i)+(w*(q1-q2)+q2)/2.d0
         

             q1=f(i+1)-f(i+2)
             q2=-q2

    		 b1=(abs(q1+q2)-abs(3.2*q1-q2))**2+ss
    		 b2=(abs(q1+q2)-abs(3.2*q2-q1))**2+ss
!!      
			!b1=q1**2+ss
			!b2=q2**2+ss

                w=1.d0/(1.d0+2.d0*(b1/b2))
                hr(i)= f(i+1)+(w*(q1-q2)+q2)/2.d0
       
           
    enddo
end

subroutine weno3_new_change(nv,jx,x,F,HL,HR)
 IMPLICIT none
 double precision ss,w0,w1,w2,c11,c12,&
				   c21,c22,e1,e2,b1,b2,q1,q2,w3,w4
 integer  n,i,j,nv,jx
 double precision hl(-nv:nv+jx),hr(-nv:nv+jx),f(-nv:nv+jx),x(-nv:nv+jx),h(-nv:nv+jx)

         SS=1E-4
		do i=-nv+1,jx+nv-1
			h(I)=(x(i)+x(i+1))/2-(x(i-1)+x(i))/2
		enddo
	!	write(*,*) h(:)
		do i=-nv,3
		h(i)= h(5)
		enddo

		do i=-5,nv
		h(jx+i)= h(jx-6)
		enddo
		!pause

        do i=-1,jx+1
		  c11= h(i-1)**2/(h(i-1)**2-h(i-2)*(2*h(i-1)+h(i-2))) !-1/2
		!j  c11=-1.d0/2
		  c12=-h(i-2)*(2*h(i-1)+h(i-2))/(h(i-1)**2-h(i-2)*(2*h(i-1)+h(i-2))) !3/2
		  !c12 =3.d0/2
			c21=h(i)**2/(h(i)**2+h(i-1)**2) !1/2
		!	c21=1.d0/2
			c22=h(i-1)**2/(h(i)**2+h(i-1)**2) !1/2

             q1=c12*f(i)+c11*f(i-1)
             q2=c21*f(i)+c22*f(i+1)

			 b1=h(i)**2/(h(i)**2+h(i-1)**2)*(f(i)-f(i-1))**2
			 b2=h(i-1)**2/(h(i)**2+h(i-1)**2)*(f(i+1)-f(i))**2

			 e2=-c11*(h(i-1)+h(i-2))**3-(c12-c11)*h(i-1)**3 !2
			 !e1=-c11*(h(i-1)+h(i-2))**3-(c12-c11)*h(i-1)**3 !2
			 e1=c22*h(i-1)**3+c21*h(i)**3
			 !e2=2.d0
			 !e1=1.d0

			 w2=e2/(e1+e2)
			 w1=e1/(e1+e2)

                w3=w1/(ss+b1)**2
                w4=w2/(ss+b2)**2
				 w1= w3/(w4+w3)
				 w2= w4/(w4+w3)
                hl(i)= w1*q1+w2*q2

	!	  c11= h(i-1)**2/(h(i-1)**2-h(i-2)*(2*h(i-1)+h(i-2)))
		  c11= h(i)**2  /(h(i)  **2-h(i+1)*(2*h(i)  +h(i+1)))
		  !c11=-1.d0/2
	!	  c12=-h(i-2)*(2*h(i-1)+h(i-2))/(h(i-1)**2-h(i-2)*(2*h(i-1)+h(i-2)))
		  c12=-h(i+1)*(2*h(i)  +h(i+1))/(h(i)  **2-h(i+1)*(2*  h(i)+h(i+1)))
		  !c12=3.d0/2
	!		c21=h(i)**2/(h(i)**2+h(i-1)**2)
			c21=h(i-1)**2/(h(i-1)**2+h(i)**2)
		!	c21=1.d0/2
	!		c22=h(i-1)**2/(h(i)**2+h(i-1)**2)
			c22=h(i)**2/(h(i-1)**2+h(i)**2)
		!	c22=1.d0/2

            q1=c12*f(i+1)+c11*f(i+2)
            q2=c22*f(i)+c21*f(i+1)

			 b1=h(i)**2/(h(i+1)**2+h(i)**2)*(f(i+2)-f(i+1))**2
		!	 b1=0.5*(f(i+2)-f(i+1))**2
			 b2=h(i+1)  **2/(h(i+1)**2+h(i)**2)*(f(i+1)-f(i))**2
		!	 b2=0.5*(f(i+1)-f(i))**2

			 e2=-(c11*(h(i)+h(i+1))**3+(c12-c11)*h(i)**3)
		!	 e2=2.d0
			 e1=c22*h(i-1)**3+c21*h(i)**3
		!	 e1=1.d0
			 !write(*,*) c11,c12,c21,c22,e1,e2
!			 read(*,*) j
			 w2=e2/(e1+e2)
			 w1=e1/(e1+e2)
!!      
                w3=w1/(ss+b1)**2
                w4=w2/(ss+b2)**2

			 w1=w3/(w4+w3)
			 w2=w4/(w4+w3)

                hr(i)= w1*q1+w2*q2
    enddo
	!pause
end

subroutine weno3_new_two_matter(nv,jx,inter,U,UL,UR,k)
      implicit none
	  integer nv,jx,i,inter,k
	  double precision u(-nv:jx+nv)
	  double precision u1(-nv:jx+nv)
	  double precision uL(-nv:jx+nv)
	  double precision uR(-nv:jx+nv)
	  double precision  IS0,IS1,c0,c1,a0,a1,b0,b1,ss,u2,d0,d1,d2,dx
	
	  ss=1.d-10

	  !***UL
	  U1=U
	 do i=-nv+1,jx+nv-1
	 if(i==inter)then
		 select  case(k)
	  case(0,3)
		  U1(i-1) = U(i)
	  endselect
	  else if(i==inter-1) then

		  select  case(k)
	  case(0,3)
		  U1(i+1)= U(i) 
		endselect
		endif

	  IS0= (U1(i)-U1(i-1))**2
	  IS1= (U1(i+1)-u1(i))**2

	  c0= 2.d0/3
	  c1= 1.d0/3
	  a0= c0/(ss+IS0)
	  a1= c1/(ss+IS1)

!	  a0=c0
!	  a1=c1

	  b0=a0/(a0+a1)
	  b1=a1/(a0+a1)

	  uR(i-1)= b0*(u1(i)+u1(i-1))/2+b1*(3*u1(i)-u1(i+1))/2

	  c0= 1.d0/3
	  c1= 2.d0/3

	  a0= c0/(ss+IS0)
	  a1= c1/(ss+IS1)
	  
!	  a0=c0
!	  a1=c1

	  b0=a0/(a0+a1)
	  b1=a1/(a0+a1)

	  uL(i)= b0*(3*u1(i)-u1(i-1))/2+b1*(u1(i)+u1(i+1))/2

	  u1(i-1)=u(i-1)
	  u1(i+1)=u(i+1)
	  enddo
 
end

