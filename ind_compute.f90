subroutine ind_compute(u4,ind)
use global 
implicit none

    integer i,j,k
    double precision u4(-nv:jx+nv,0:2),f_t(-nv:jx+nv,0:2)
    integer ind(-nv:jx+nv,0:2),ind1(-nv:jx+nv,0:2)
                        
    
    f_t=0
    ind1=1
    ind=1

    do k=0,2
            call indicate_choose(ind(:,k),nv,jx,f_t(:,0),u4(:,k),nk,dx)            
        
    enddo
ind(:,0)=ind(:,0)*ind(:,1)*ind(:,2)
!ind1=0
!ind=0
!    
!    ind(:,:,0)=ind(:,:,0)*ind(:,:,1)*ind(:,:,2)*ind(:,:,3)
!    ind1(:,:,0)=ind1(:,:,0)*ind1(:,:,1)*ind1(:,:,2)*ind1(:,:,3)
    !ind1=0
    !ind=0
!    ind=ind*ind1
    !ind1=ind

    !do i=-nv,jx+nv
    !    do j=0,jy
    !       do k=-nv+j,j+nv
    !          if(ind(i,k,0)==0)then
    !            ind1(i,j,0)=0
    !          endif
    !       enddo
    !    enddo
    !enddo
    !
    !do i=0,jx
    !    do j=-nv,jy+nv
    !       do k=-nv+i,i+nv
    !          if(ind(k,j,0)==0)then
    !            ind1(i,j,0)=0
    !          endif
    !       enddo
    !    enddo
    !enddo
    !ind=ind1          
    ind(-nv:0,0)=0
    ind(jx:jx+nv,0)=0

    !ind=0
    !


end subroutine          
