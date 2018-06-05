subroutine timesolve
 use global
implicit none
integer it
   
  
   T=0
   do it=1,1000000000
       

   call cfl(u)
   
      if(t>=tt)then
         exit
      else if(t+dt>tt)then
         dt=tt-t
      endif   
      t=t+dt
      write(*,*)'T=',T,'dt=',dt
      call R_K
      
   enddo
           
end  
   
