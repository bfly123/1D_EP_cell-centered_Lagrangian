subroutine timesolve
 use global
implicit none
double precision dt,t
integer ti
  
   T=0
   ti=0
  do while(t<tt)
      call CFL(dt)
	  if(t+dt>tt)then
		dt=tt-t
      endif   
      t=t+dt
	  ti=ti+1
	  call R_K(dt)
      write(*,*)'T=',T,'dt=',dt
!stop
   enddo
end  
   
