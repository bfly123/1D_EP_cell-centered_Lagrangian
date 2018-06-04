subroutine weno5_new(nv,jx,F,HL,HR)
 IMPLICIT none
 double precision ss,a300,a301,a302,a310,a311,a312,a320,&
                   a321,a322,c30,c31,c32,is0,is1,is2,tao5,&
                   q30,q31,q32,aa0,aa1,aa2,w0,w1,w2
 integer  n,nv,jx,i
 double precision hl(-nv:nv+jx),hr(-nv:nv+jx),f(-nv:nv+jx)
         SS=1E-5
         a300=1.d0/3.d0
         a301=-7.d0/6.d0
         a302=11.d0/6.d0
         a310=-1.d0/6.d0
         a311=5.d0/6.d0
         a312=1.d0/3.d0
         a320=1.d0/3.d0
         a321=5.d0/6.d0
         a322=-1.d0/6.d0
         c30=1.d0/10.d0
         c31=6.d0/10.d0
         c32=3.d0/10.d0 

        do i=-1,jx
          IS0=13.d0*(F(I-2)-2.*F(I-1)+F(I))**2/12.d0&
                          +(F(I-2)-4.*F(I-1)+3.*F(I))**2/4.d0
          IS1=13.d0*(F(I-1)-2.*F(I)+F(I+1))**2/12.d0&
                          +(F(I-1)-F(I+1))**2/4.d0
          IS2=13.d0*(f(I)-2.*F(I+1)+f(I+2))**2/12.d0&
                         +(3.*f(I)-4.*f(I+1)+f(I+2))**2/4.d0
             tao5=abs(is0-is2)
             
             Q30=A300*f(i-2)+A301*f(i-1)+A302*f(i)
             Q31=A310*f(i-1)+A311*f(i)+A312*f(i+1)
             Q32=A320*f(i)+A321*f(i+1)+A322*f(i+2)
            
            aa0=C30/(ss+IS0)**2
            aa1=c31/(ss+IS1)**2
            aa2=c32/(ss+IS2)**2

!                aa0=C30*(1.0d0+(tao5/(ss+IS0))**2)
!                aa1=C31*(1.0d0+(tao5/(ss+IS1))**2)
!                aa2=C32*(1.0d0+(tao5/(ss+IS2))**2)   
!      
                w0=aa0/(aa0+aa1+aa2)
                w1=aa1/(aa0+aa1+aa2)
                w2=aa2/(aa0+aa1+aa2)
                hl(i)=w0*Q30+w1*Q31+w2*Q32 
         
!         aa0=0.1*(1.0d0+(tao5/(ss+IS0)))
!          aa1=0.6*(1.0d0+(tao5/(ss+IS1)))
!          aa2=0.3*(1.0d0+(tao5/(ss+IS2))) 
            !aa0=0.1d0/(ss+IS0)**2
            !aa1=0.6d0/(ss+IS1)**2
            !aa2=0.3d0/(ss+IS2)**2

!     h(i)=1.d0/30*lf(i-2)-13.d0/60*lf(i-1)+47.d0/60*lf(i)+27.d0/60*lf(i+1)-1.d0/20*lf(i+2)
  
         enddo
    do i=-1,jx     
          IS0=13.d0*(f(I+3)-2.*f(I+2)+f(I+1))**2/12.d0&
                          +(f(I+3)-4.*f(I+2)+3.*f(I+1))**2/4.d0
          IS1=13.d0*(f(I+2)-2.*f(I+1)+f(I))**2/12.d0&
                          +(f(I+2)-f(I))**2/4.d0
          IS2=13.d0*(f(I+1)-2.*f(I)+f(I-1))**2/12.d0&
                          +(3.*f(I+1)-4.*f(I)+f(I-1))**2/4.d0
          tao5=abs(is0-is2)
         
            Q30=A300*f(i+3)+A301*f(i+2)+A302*f(i+1)
             Q31=A310*f(i+2)+A311*f(i+1)+A312*f(i)
             Q32=A320*f(i+1)+A321*f(i)+A322*f(i-1)
            
                aa0=C30/(ss+IS0)**2
                aa1=c31/(ss+IS1)**2
                aa2=c32/(ss+IS2)**2

!                aa0=C30*(1.0d0+(tao5/(ss+IS0))**2)
!                aa1=C31*(1.0d0+(tao5/(ss+IS1))**2)
!                aa2=C32*(1.0d0+(tao5/(ss+IS2))**2)          


                 w0=aa0/(aa0+aa1+aa2)
                 w1=aa1/(aa0+aa1+aa2)
                 w2=aa2/(aa0+aa1+aa2)

                 hr(i)=w0*Q30+w1*Q31+w2*Q32
      
!          aa0=0.1*(1.0d0+(tao5/(ss+IS0))**2)
!          aa1=0.6*(1.0d0+(tao5/(ss+IS1))**2)
!          aa2=0.3*(1.0d0+(tao5/(ss+IS2))**2) 
              

    enddo
end        
