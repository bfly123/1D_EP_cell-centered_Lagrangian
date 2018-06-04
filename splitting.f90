subroutine SW_splitting(u4)
use global
use constant
implicit none
integer i,j,k
double precision tmp0,tmp1,tmp2,tmp3
double precision uu,rho,p,cc
double precision e1,e2,e3,e1p,e2p,e3p,e1m,e2m,e3m
double precision uc1,uc2,vc1,vc2,vvc1,vvc2,vv,w2
double precision u4(-nv:jx+nv,0:2)
double precision ss
double precision nx1,ny1
!c El 为特征值,其中El(:,1)为x方向的特征值（4个， u, u, u+c, u-c)
        
	  tmp1=2.d0*(gama-1.d0)
	  tmp2=1.d0/(2.d0*gama)
	  tmp3=(3.d0-gama)/(2.d0*(gama-1.d0)) 
      ss=1e-2
            
	  do i=-nv,jx+nv
       
	  uu=u4(i,1)/u4(i,0)
	  rho=u4(i,0)
	  p=(gama-1)*(u4(i,2)-rho*(uu*uu)/2.d0)
	  cc=sqrt(gama*p/rho)
	  E1=uu
      E2=uu-cc
      E3=uu+cc
	  E1P=(E1+abs(E1)+ss)/2.d0
	  E2P=(E2+abs(E2)+ss)/2.d0
	  E3P=(E3+abs(E3)+ss)/2.d0
	  E1M=E1-E1P
	  E2M=E2-E2P
	  E3M=E3-E3P
	  
      tmp0=rho/(2.d0*gama) 
	  uc1=uu-cc
	  uc2=uu+cc
	  vvc1=(uc1*uc1)/2.d0
	  vvc2=(uc2*uc2)/2.d0
	  vv=(gama-1.d0)*(uu*uu)
      W2=tmp3*cc*cc

        fp(i,0)=tmp0*(tmp1*E1P+E2P+E3P)
        fp(i,1)=tmp0*(tmp1*E1P*uu+E2P*uc1+E3P*uc2)
        fp(i,2)=tmp0*(E1P*vv+E2p*vvc1+E3P*vvc2+W2*(E2P+E3P))
        
        fm(i,0)=tmp0*(tmp1*E1M+E2M+E3M)
        fm(i,1)=tmp0*(tmp1*E1M*uu+E2M*uc1+E3M*uc2)
        fm(i,2)=tmp0*(E1M*vv+E2M*vvc1+E3M*vvc2+W2*(E2M+E3M))
     
     enddo

end

subroutine LF_splitting(u4)
use global
use constant
implicit none
integer i,j,k
double precision tmp0,tmp1,tmp2,tmp3
double precision uu,uv,rho,p,cc,u_max,u_max_tmp
double precision e1,e2,e3
double precision uc1,uc2,vc1,vc2,vvc1,vvc2,vv,w2
double precision u4(-nv:jx+nv,0:2)
double precision ss
integer nx1,ny1
!c El 为特征值,其中El(:,1)为x方向的特征值（4个， u, u, u+c, u-c)
        
	  tmp1=2.d0*(gama-1.d0)
	  tmp2=1.d0/(2.d0*gama)
	  tmp3=(3.d0-gama)/(2.d0*(gama-1.d0)) 
      ss=1.e-10
      
	  u_max=0
      do i=-nv,jx+nv
        
	  uu=u4(i,1)/u4(i,0)
	  rho=u4(i,0)
	  p=(gama-1)*(u4(i,2)-rho*(uu*uu)/2.d0)
	  cc=sqrt(gama*p/rho)
	  
	  E1=uu
      E2=uu-cc
      E3=uu+cc
      u_max_tmp=abs(uu)+cc
	  if (u_max_tmp.ge.u_max)then
		  u_max=u_max_tmp
	endif
	  enddo
	  do i=-nv,jx+nv
        
	  uu=u4(i,1)/u4(i,0)
	  rho=u4(i,0)
	  p=(gama-1)*(u4(i,2)-rho*(uu*uu)/2.d0)
	  cc=sqrt(gama*p/rho)
	  
	  E1=uu
      E2=uu-cc
      E3=uu+cc
     ! u_max=abs(uu)+cc

      tmp0=rho/(2.d0*gama) 
	  uc1=uu-cc
	  uc2=uu+cc
	  vvc1=(uc1*uc1)/2.d0
	  vvc2=(uc2*uc2)/2.d0
	  vv=(gama-1.d0)*(uu*uu)
      W2=tmp3*cc*cc

     f(i,0)=tmp0*(tmp1*E1+E2+E3)
     f(i,1)=tmp0*(tmp1*E1*uu+E2*uc1+E3*uc2)
     f(i,2)=tmp0*(E1*vv+E2*vvc1+E3*vvc2+W2*(E2+E3))
    
      fp(i,:)=1.d0/2*(f(i,:)+u_max*U4(i,:))
      fm(i,:)=1.d0/2*(f(i,:)-u_max*U4(i,:))
     enddo
    
   
       end

