subroutine weno5(nv,LF,RF,H)
 IMPLICIT none
 double precision ss,a300,a301,a302,a310,a311,a312,a320,&
                   a321,a322,c30,c31,c32,is0,is1,is2,tao5,&
                   q30,q31,q32,aa0,aa1,aa2,w0,w1,w2
 integer  nv,n,i
 double precision lf(-nv:nv),rf(-nv:nv)
 double precision h,hr
         SS=1E-20
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
         
        i=0 
          IS0=13.d0*(LF(I-2)-2.*LF(I-1)+LF(I))**2/12.d0&
                          +(LF(I-2)-4.*LF(I-1)+3.*LF(I))**2/4.d0
          IS1=13.d0*(LF(I-1)-2.*LF(I)+LF(I+1))**2/12.d0&
                          +(LF(I-1)-LF(I+1))**2/4.d0
          IS2=13.d0*(LF(I)-2.*LF(I+1)+LF(I+2))**2/12.d0&
                         +(3.*LF(I)-4.*LF(I+1)+LF(I+2))**2/4.d0
             tao5=abs(is0-is2)

             Q30=A300*LF(i-2)+A301*LF(i-1)+A302*LF(i)
             Q31=A310*LF(i-1)+A311*LF(i)+A312*LF(i+1)
             Q32=A320*LF(i)+A321*LF(i+1)+A322*LF(i+2)


         aa0=C30/(ss+IS0)**2
          aa1=C31/(ss+IS1)**2
          aa2=C32/(ss+IS2)**2
         
         aa0=0.1*(1.0d0+(tao5/(ss+IS0))**2)
          aa1=0.6*(1.0d0+(tao5/(ss+IS1))**2)
          aa2=0.3*(1.0d0+(tao5/(ss+IS2))**2) 
             
                w0=aa0/(aa0+aa1+aa2)
                w1=aa1/(aa0+aa1+aa2)
                w2=aa2/(aa0+aa1+aa2)
                h=w0*Q30+w1*Q31+w2*Q32
     
  
         
          IS0=13.*(RF(I+3)-2.*RF(I+2)+RF(I+1))**2/12.&
                          +(RF(I+3)-4.*RF(I+2)+3.*RF(I+1))**2/4.
          IS1=13.*(RF(I+2)-2.*RF(I+1)+RF(I))**2/12.&
                          +(RF(I+2)-RF(I))**2/4.
          IS2=13.*(RF(I+1)-2.*RF(I)+RF(I-1))**2/12.&
                          +(3.*RF(I+1)-4.*RF(I)+RF(I-1))**2/4.
          tao5=abs(is0-is2)
         
             Q30=A300*RF(i+3)+A301*RF(i+2)+A302*RF(i+1)
             Q31=A310*RF(i+2)+A311*RF(i+1)+A312*RF(i)
             Q32=A320*RF(i+1)+A321*RF(i)+A322*RF(i-1)


          aa0=C30/(ss+IS0)**2
          aa1=C31/(ss+IS1)**2
          aa2=C32/(ss+IS2)**2
      
          aa0=0.1*(1.0d0+(tao5/(ss+IS0))**2)
          aa1=0.6*(1.0d0+(tao5/(ss+IS1))**2)
          aa2=0.3*(1.0d0+(tao5/(ss+IS2))**2) 
              
                 w0=aa0/(aa0+aa1+aa2)
                 w1=aa1/(aa0+aa1+aa2)
                 w2=aa2/(aa0+aa1+aa2)

                 hr=w0*Q30+w1*Q31+w2*Q32
         h=hr+h
end        
