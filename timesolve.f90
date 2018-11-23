subroutine timesolve
 use global
implicit none
double precision dt,t
integer ti,t1
  
   T=0
   ti=0
  do while(t<tt)
      call CFL(dt)

	  if(t+dt>tt)then
		dt=tt-t
      endif   
	  t1=t
      t=t+dt
      write(*,*)'T=',T,'dt=',dt
	  ti=ti+1
	 call R_K(t1,dt)
!	  call ADER(t1,dt)
!stop
   enddo
end  
   
