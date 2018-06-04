subroutine indicate_choose(ind,nv,n,lf,rf,nk,dx)
implicit none 
integer nv,n,nk
integer ind(-nv:n+nv)
 double precision lf(-nv:n+nv),RF(-nv:n+nv)
 double precision dx      
 select case(nk)
    case(1)
        call indicate_ATV(ind,nv,n,lf,rf,dx)
    case(2)  
         call indicate_TVB(ind,nv,n,lf,rf,dx)     
    case(3)
         call indicate_XS(ind,nv,n,lf,rf,dx)
    case(4)
         call indicate_MP(ind,nv,n,lf,rf,dx)
    case(5)
         call indicate_MR(ind,nv,n,lf,rf,dx)
    case(6) 
         call indicate_BDF(ind,nv,n,lf,rf,dx) 
    case(7)
         call indicate_BSB(ind,nv,n,lf,rf,dx) 
    case(8)
         call indicate_MMP(ind,nv,n,lf,rf,dx)     
    case(9)
         call indicate_KXRCF(ind,nv,n,lf,rf,dx)    
    case(10)
         call indicate_SHEN(ind,nv,n,lf,rf,dx)     
 end select  


end