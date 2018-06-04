subroutine characterReconstruct(u4,h,ind)
use constant
use global
implicit none

    integer i,j,k,m,n,nx1,ny1,ii,jj,ntc,flag
    double precision u4(-nv:jx+nv,0:2),f1(-nv:jx+nv,0:2),h(-nv:jx+nv,0:2),&
                        u5(-2:3,0:2),xh(10000,0:2)
    double precision a3(10000),b3(1000),c3(10000),&
                      d3(10000,0:2)
    integer start_p(-1:10000),end_p(-1:10000),ind(-nv:jx+nv)

    double precision d1,uu1,uv1,p1,h1,d2,uu2,uv2,p2,h2,&
                      r1,r2,r0,u0,v0,h0,c0,t0,v_max
    double precision r(3,3),l(3,3),a_p(3),a_m(3)
    double precision v_p(-nv:nv,3),h_v(0:2),v_m(-nv:nv,3)
    double precision ss
    ss=1.0e-10
!    if ( kind_problem==6)then
!      nx1=int(0.6/dx+ss)
!      ny1=int(0.2/dy+ss)
!      do i=1,nv
!        do j=1,nv
!           u4(nx1+i,ny1-j,:)=corner(i,j,:,1)
!        enddo
!      enddo
!    endif 

           
        start_p=jx+2
        NTC=1
        start_p(NTC)=-1
        
        do i=-1,jx
            flag=0
!            if(i==-1.or.i==jx)then
!                flag=1
!            else 
                do j=-3+i,i+4
                    if(ind(j)==0)then
                        flag=1
                    endif
                enddo
!            endif
           if (flag==1)then
    !          if(i>=1.and.i<=jx-1)then   
                 end_p(NTC)=i
                 NTC=NTC+1
                 start_p(NTC)=i+1 
    !          endif
                 
                d1=U4(i,0)
                uu1=u4(i,1)/u4(i,0)
                p1=(GAMA-1)*(U4(i,2)-0.5*U4(i,0)*(uu1*uu1))
                H1=gama/(gama-1.d0)*p1/d1+(uu1*uu1)*0.5d0 
                
                
                d2=u4(i+1,0)
                uu2=u4(i+1,1)/u4(i+1,0)
                p2=(GAMA-1)*(U4(i+1,2)-0.5*U4(i+1,0)*(uu2*uu2))
                H2=gama/(gama-1.d0)*p2/d2+(uu2*uu2)*0.5d0 
            
                r1=sqrt(d1) ; r2= sqrt(d2); r0= r1+r2
                u0=(r1*uu1+r2*uu2)/r0
                H0=(r1*H1+r2*H2)/r0                     
                c0=sqrt((gama-1.d0)*(H0-u0*u0*0.5d0))  
        	 
                !ÌØÕ÷¾ØÕóRx
                r(1,1)=1
                r(1,2)=1
                r(1,3)=1
        
                r(2,1)=u0
                r(2,2)=u0+c0
                r(2,3)=u0-c0
         
                r(3,1)=u0*u0/2
                r(3,2)=u0*u0/2+u0*c0+c0*c0/(gama-1)
                r(3,3)=u0*u0/2-u0*c0+c0*c0/(gama-1)
             


                !×óÌØÕ÷¾ØÕó
                t0=0.5/(c0*c0)

                l(1,1)=t0*(2*c0*c0-(gama-1)*u0*u0)
                l(1,2)=t0*(2*u0*(gama-1))
                l(1,3)=-t0*2*(gama-1)

                l(2,1)=t0*(-c0*u0+(gama-1)/2*u0*u0)
                l(2,2)=t0*(c0-(gama-1)*u0)
                l(2,3)=t0*(gama-1)


                l(3,1)=t0*(c0*u0+(gama-1)/2*u0*u0)
                l(3,2)=t0*(-c0-(gama-1)*u0)
                l(3,3)=t0*(gama-1)

            
    !        a_p(1)=0.5*(u0+abs(u0))
    !        a_p(2)=0.5*(u0+abs(u0))
    !        a_p(3)=0.5*(u0-c0+abs(u0-c0))
    !        a_p(4)=0.5*(u0+c0+abs(u0+c0))
    !        
    !        a_m(1)=0.5*(u0-abs(u0))
    !        a_m(2)=0.5*(u0-abs(u0))
    !        a_m(3)=0.5*(u0-c0-abs(u0-c0))
    !        a_m(4)=0.5*(u0+c0-abs(u0+c0))
                v_max=max(abs(u0),abs(u0+c0),abs(u0-c0))  
         
                a_p(1)=0.5*(u0+v_max)
                a_p(2)=0.5*(u0+c0+v_max)
                a_p(3)=0.5*(u0-c0+v_max)

                a_m(1)=0.5*(u0-v_max)
                a_m(2)=0.5*(u0+c0-v_max)
                a_m(3)=0.5*(u0-c0-v_max)            
         !********A=RaL
         !********v=aLu   
                do m=-nv+1,nv
                    do k=1,3
                        v_p(m,k)=0
                        v_m(m,k)=0
                        do n=0,2
                            v_p(m,k)=v_p(m,k)+a_p(k)*l(k,n+1)*u4(i+m,n)
                            v_m(m,k)=v_m(m,k)+a_m(k)*l(k,n+1)*u4(i+m,n)
                        enddo
                    enddo   
                enddo
    !         do m=-nv+1,nv
    !            do k=1,4
    !                v_p(m,k)=0
    !                v_m(m,k)=0
    !                do n=0,3
    !                    v_p(m,k)=v_p(m,k)+l(k,n+1)*fp(i+m,j,n)
    !                    v_m(m,k)=v_m(m,k)+l(k,n+1)*fm(i+m,j,n)
    !                enddo
    !            enddo   
    !        enddo
            
                do k=0,2
                    call weno5(nv,v_p(-nv:nv,k+1),v_m(-nv:nv,k+1),H_v(k))
                enddo
     !f_half=RV_P       
            
                do k=0,2
                    h(i,k)=0
                    do m=0,2
                        h(i,k)=h(i,k)+r(k+1,m+1)*H_v(m)
                    enddo
                enddo
!            else
!                do ii=-2,3
!                    call u_f(u4(ii+i,j,:),u5(ii,:))
!                enddo 
!                
!                do n=0,3
!                    h(i,j,n)=1.d0/60*u5(3,n)-2.d0/15*u5(2,n)+37.d0/60*u5(1,n)+1.d0/60*u5(-2,n)-2.d0/15*u5(-1,n)+37.d0/60*u5(0,n)  
!                enddo
            
            endif       
        enddo !do i=-1,jx
!    enddo !j
      
       end_p(NTC)=jx
       
       do k=1,NTC
          m=end_p(k)-start_p(k)
          if(m>=1)then
            do jj=1,m
                i=start_p(k)+jj-1
                do ii=-1,2
                    call u_f(u4(ii+i,:),u5(ii,:))
                enddo 
                
!                do n=0,3
!                    h(i,j,n)=1.d0/60*u5(3,n)-2.d0/15*u5(2,n)+37.d0/60*u5(1,n)+1.d0/60*u5(-2,n)-2.d0/15*u5(-1,n)+37.d0/60*u5(0,n)  
!                enddo
                a3(jj)=1.0d0/3.0d0
                b3(jj)=1.0d0
                c3(jj)=a3(jj)
                do n=0,2
                      d3(jj,n)=(u5(2,n)+29.d0*u5(1,n)+29.d0*u5(0,n)+u5(-1,n))&
                               /36.d0
                enddo
            end do !jj=1,m
            
            do n=0,2
                d3(1,n)=d3(1,n)-a3(1)*h(start_p(k)-1,n)
                d3(m,n)=d3(m,n)-c3(m)*h(end_p(k),n)
                 call tridiagsolve(m,a3,b3,c3,d3(:,n),xh(:,n))
                do jj=1,m
                   h(start_p(k)+jj-1,n)=xh(jj,n)
!                   h(start_p(k)+jj-1,j,n)=1.d0/60*u5(3,n)-2.d0/15*u5(2,n)+37.d0/60u5(1,n)+1.d0/60*u5(-2,n)-2.d0/15*u5(-1,n)+37.d0/60u5(0,n)
                end do
               
           enddo
        endif
    enddo !k=1,NTC

end
