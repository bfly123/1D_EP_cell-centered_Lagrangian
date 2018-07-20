subroutine upwind3(nv,jx,F,HL,HR)
 IMPLICIT none
 double precision ss,a300,a301,a302,a310,a311,a312,a320,&
                   a321,a322,c30,c31,c32,is0,is1,is2,tao5,&
                   q30,q31,q32,aa0,aa1,aa2,w0,w1,w2
 integer  n,nv,jx,i
 double precision hl(-nv:nv+jx),hr(-nv:nv+jx),f(-nv:nv+jx)
         SS=1E-15
		 do i=-1,jx+1
hl(i)=-1.d0/6*f(i-1)+5.d0/6*f(i)+1.d0/3*f(i+1)
hr(i)=-1.d0/6*f(i+2)+5.d0/6*f(i+1)+1.d0/3*f(i)
		 enddo
		 end

