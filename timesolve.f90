subroutine timesolve
 use global
implicit none
integer it
double precision dt,t
  
   T=0
  do while(t< tt)
      call CFL(dt)
	  if(t+dt>tt)then
		dt=tt-t
      endif   
      t=t+dt
      write(*,*)'T=',T,'dt=',dt
   enddo
end  
   
