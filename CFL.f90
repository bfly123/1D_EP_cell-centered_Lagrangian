      subroutine CFL(u4)
      use constant
      use global
      implicit none
      double precision dmaxvel
      double precision uu
      double precision vv
      double precision vel
      double precision p
      double precision CF
!      double precision dt
      integer i
      integer j
      double precision u4(-nv:jx+nv,0:2)
      
      dmaxvel=1e-10
      do  i=0,Jx
           uu=U4(i,1)/U4(i,0)
           p=(GAMA-1)*(U(i,2)-0.5*U4(i,0)*(uu*uu))  
           vel=dsqrt(GAMA*p/U4(i,0))+dsqrt(uu*uu)
           if(vel.gt.dmaxvel)dmaxvel=vel
        enddo
     
     dt=Sf*dx/dmaxvel
!     dt=cf*dx
      end
