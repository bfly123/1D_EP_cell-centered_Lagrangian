subroutine timesolve
 use global
implicit none
double precision dt,t,time_begin,time_end
integer ti,t1,i
  
   T=0
   ti=0
   ! do while(t<tt)
   call cpu_time(time_begin)
   do i = 1,100

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
   call cpu_time(time_end)
      write(*,*) time_end - time_begin ,"Second"
end  
   
