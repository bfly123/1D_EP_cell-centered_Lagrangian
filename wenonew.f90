      SUBROUTINE WENONEW(nv,n,lf,rf,h,ind)
!      SUBROUTINE WENO(nv,jx,fp,fm,h,ind)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!      implicit none
      INTEGER SCHEME,nv,n
      
      DOUBLE PRECISION LF(-nv:n+nv),RF(-nv:n+nv),H(-nv:n+nv),HR(-nv:n+nv)
      INTEGER          end_p(-1:4000),start_p(-1:4000),ind(-nv:n+nv)
     
      DOUBLE PRECISION IS0,IS1,IS2
      DOUBLE PRECISION a3(4000),b3(4000),&
                       c3(4000),d3(4000),xh(4000)



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
!-------------boundry
       do j=1,2
          if (j.eq.1) i=-1
          if (j.eq.2) i=N

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
                   
!------------------
          aa0=C30*(1.0d0+(tao5/(ss+IS0))**2)
          aa1=C31*(1.0d0+(tao5/(ss+IS1))**2)
          aa2=C32*(1.0d0+(tao5/(ss+IS2))**2) 
                

          w0=aa0/(aa0+aa1+aa2)
          w1=aa1/(aa0+aa1+aa2)
          w2=aa2/(aa0+aa1+aa2)

          h(i)=w0*Q30+w1*Q31+w2*Q32
       end do


       start_p=N+1
       NTC=1
       start_p(NTC)=0
       do I=0,N-1
          IS0=13.d0*(LF(I-2)-2.*LF(I-1)+LF(I))**2/12.d0&
                          +(LF(I-2)-4.*LF(I-1)+3.*LF(I))**2/4.d0
          IS1=13.d0*(LF(I-1)-2.*LF(I)+LF(I+1))**2/12.d0&
                          +(LF(I-1)-LF(I+1))**2/4.d0
          IS2=13.d0*(LF(I)-2.*LF(I+1)+LF(I+2))**2/12.d0&
                         +(3.*LF(I)-4.*LF(I+1)+LF(I+2))**2/4.d0
          tao5=abs(is0-is2)
          tao0=abs(is0-is1)
          tao1=abs(is2-is1)
          dmin=dmin1(is0,is1,is2)
!          if (tao5.gt.dmin) then 
            if(ind(i)==0.or.ind(i+1)==0.or.ind(i-1)==0.or.ind(i-2)==0.or.ind(i+2)==0.or.ind(i+3)==0)then    
             end_p(NTC)=i
             NTC=NTC+1
             start_p(NTC)=i+1

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


                aa0=C30*(1.0d0+(tao5/(ss+IS0))**2)
                aa1=C31*(1.0d0+(tao5/(ss+IS1))**2)
                aa2=C32*(1.0d0+(tao5/(ss+IS2))**2)   
         
                w0=aa0/(aa0+aa1+aa2)
                w1=aa1/(aa0+aa1+aa2)
                w2=aa2/(aa0+aa1+aa2)
                h(i)=w0*Q30+w1*Q31+w2*Q32
!             end if
          end if
       end do            ! end of (do i=5,n-5)
        end_p(NTC)=N
       do k=1,NTC
          m=end_p(k)-start_p(k)
          if(m>=1)then
             do j=1,m
                i=start_p(k)+j-1
                a3(j)=1.0d0/3.0d0
                b3(j)=1.0d0
                c3(j)=a3(j)
                d3(j)=(lf(i+2)+29.d0*lf(i+1)+29.d0*lf(i)+lf(i-1))&
                 /36.d0


             end do
             d3(1)=d3(1)-a3(1)*h(start_p(k)-1)
             d3(m)=d3(m)-c3(m)*h(end_p(k))
             call tridiagsolve(m,a3,b3,c3,d3,xh)
             do j=1,m
                h(start_p(k)+j-1)=xh(j)
             end do
            endif
       end do          ! end of (do k=1,NTC)
!----------------

!c****************************************************************
!c****************************************************************
!c************************************http://www4.ncsu.edu/eos/users/w/white/www/white/ma580/chap2.5.PDF****************************
!c********************Negative PART*******************************
       do j=1,2
          if (j.eq.1) i=-1
          if (j.eq.2) i=N
          IS0=13.*(RF(I+3)-2.*RF(I+2)+RF(I+1))**2/12.&
                          +(RF(I+3)-4.*RF(I+2)+3.*RF(I+1))**2/4.
          IS1=13.*(RF(I+2)-2.*RF(I+1)+RF(I))**2/12.&
                          +(RF(I+2)-RF(I))**2/4.
          IS2=13.*(RF(I+1)-2.*RF(I)+RF(I-1))**2/12.&
                          +(3.*RF(I+1)-4.*RF(I)+RF(I-1))**2/4.


          Q30=A300*RF(i+3)+A301*RF(i+2)+A302*RF(i+1)
          Q31=A310*RF(i+2)+A311*RF(i+1)+A312*RF(i)
          Q32=A320*RF(i+1)+A321*RF(i)+A322*RF(i-1)



            aa0=C30*(1.0d0+(tao5/(ss+IS0))**2)
            aa1=C31*(1.0d0+(tao5/(ss+IS1))**2)
            aa2=C32*(1.0d0+(tao5/(ss+IS2))**2)          

       

            w0=aa0/(aa0+aa1+aa2)
            w1=aa1/(aa0+aa1+aa2)
            w2=aa2/(aa0+aa1+aa2)

            hr(i)=w0*Q30+w1*Q31+w2*Q32
         end do



       start_p=N+1
       NTC=1
       start_p(NTC)=0
       do I=0,N-1
           IS0=13.*(RF(I+3)-2.*RF(I+2)+RF(I+1))**2/12.&
                          +(RF(I+3)-4.*RF(I+2)+3.*RF(I+1))**2/4.
          IS1=13.*(RF(I+2)-2.*RF(I+1)+RF(I))**2/12.&
                          +(RF(I+2)-RF(I))**2/4.
          IS2=13.*(RF(I+1)-2.*RF(I)+RF(I-1))**2/12.&
                          +(3.*RF(I+1)-4.*RF(I)+RF(I-1))**2/4.
          tao5=abs(is0-is2)
          tao0=abs(is0-is1)
          tao1=abs(is2-is1)
          dmin=dmin1(is0,is1,is2)
!          if (tao5.gt.dmin) then   
            if(ind(i+1)==0.or.ind(i+2)==0.or.ind(i)==0.or.ind(i-1)==0.or.ind(i+3)==0.or.ind(i+4)==0)then      
             end_p(NTC)=i
             NTC=NTC+1
             start_p(NTC)=i+1

             Q30=A300*RF(i+3)+A301*RF(i+2)+A302*RF(i+1)
             Q31=A310*RF(i+2)+A311*RF(i+1)+A312*RF(i)
             Q32=A320*RF(i+1)+A321*RF(i)+A322*RF(i-1)


                aa0=C30*(1.0d0+(tao5/(ss+IS0))**2)
                aa1=C31*(1.0d0+(tao5/(ss+IS1))**2)
                aa2=C32*(1.0d0+(tao5/(ss+IS2))**2)          

      

                 w0=aa0/(aa0+aa1+aa2)
                 w1=aa1/(aa0+aa1+aa2)
                 w2=aa2/(aa0+aa1+aa2)

                 hr(i)=w0*Q30+w1*Q31+w2*Q32
!            end if
          end if
       end do            ! end of (do i=5,n-5)
        end_p(NTC)=N
       do k=1,NTC
          m=end_p(k)-start_p(k)
            if(m>=1)then
             do j=1,m
                i=start_p(k)+j-1
                a3(j)=1.0d0/3.0d0
                b3(j)=1.0d0
                c3(j)=a3(j)
                d3(j)=(rf(i+2)+29.d0*rf(i+1)+29.d0*rf(i)+rf(i-1))&
                 /36.d0
      
             end do
             d3(1)=d3(1)-a3(1)*hr(start_p(k)-1)
             d3(m)=d3(m)-c3(m)*hr(end_p(k))
             call tridiagsolve(m,a3,b3,c3,d3,xh)
             do j=1,m
                hr(start_p(k)+j-1)=xh(j)
             end do
           endif
       end do            ! end of (do k=1,NTC)
       
       h=h+hr
 END