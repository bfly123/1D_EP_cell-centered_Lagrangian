subroutine weno3_new(nv,jx,F,HL,HR)
 IMPLICIT none
 double precision ss,a300,a301,a302,a310,a311,a312,a320,&
                   a321,a322,c30,c31,c32,is0,is1,is2,tao5,&
                   q30,q31,q32,aa0,aa1,aa2,w0,w1,w2
 integer  n,nv,jx,i
 double precision hl(-nv:nv+jx),hr(-nv:nv+jx),f(-nv:nv+jx)
         SS=1E-4
        
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
				   c21,c22,e1,e2,b1,b2,q1,q2
 integer  n,nv,jx,i,j
 double precision hl(-nv:nv+jx),hr(-nv:nv+jx),f(-nv:nv+jx),x(-nv:nv+jx),h(-nv:nv+jx)

         SS=1E-4
		do i=-nv+1,jx+nv-1
			h(I)=(x(i)+x(i+1))/2-(x(i-1)+x(i))/2
			write(*,*) h(i)
		enddo

        do i=-1,jx+1
			c11=-h(i-1)/(h(i-2)*(h(i-1)+h(i-2)))
			c12=(h(i-2)**2+2*h(i-1)*h(i-2))/(h(i-1)*h(i-2)*(h(i-1)+h(i-2)))
			c21=(h(i)/h(i-1)/(h(i)+h(i-1)))
			c22=(h(i-1)/h(i)/(h(i)+h(i-1)))

             q1=c12*f(i)+c11*f(i-1)
             q2=c21*f(i)+c22*f(i+1)

			 b1=(f(i)-f(i-1))**2
			 b2=(f(i+1)-f(i))**2

			 e1=c11*h(i)**3 +c11*(h(i-1)+h(i))**3+c12*h(i)**3
			 e2=c22*h(i+1)**3+c21*h(i)**3
			 w2=e1/(e1+e2)
			 w1=1-e1/(e1+e2)
!!      
                w1=w1/(ss+b1)**2
                w2=w2/(ss+b2)**2
                hl(i)= w1*q1+w2*q2

			c11=-h(i)/(h(i+1)*(h(i)+h(i+1)))
			c12=(h(i+1)**2+2*h(i)*h(i+1))/(h(i)*h(i+1)*(h(i)+h(i+1)))
			c21=(h(i-1)/h(i)/(h(i)+h(i-1)))
			c22=(h(i)/h(i-1)/(h(i)+h(i-1)))

            q1=c12*f(i+1)+c11*f(i+2)
            q2=c22*f(i)+c21*f(i+1)

			 b1=(f(i+2)-f(i+1))**2
			 b2=(f(i+1)-f(i))**2

			 e1=c11*h(i+1)**3 +c11*(h(i+2)+h(i+1))**3+c12*h(i+1)**3
			 e2=c22*h(i)**3+c21*h(i+1)**3
			 !write(*,*) c11,c12,c21,c22,e1,e2
			 read(*,*) j
			 w2=e1/(e1+e2)
			 w1=1-e1/(e1+e2)
!!      
                w1=w1/(ss+b1)**2
                w2=w2/(ss+b2)**2
                hr(i)= w1*q1+w2*q2
    enddo
end
