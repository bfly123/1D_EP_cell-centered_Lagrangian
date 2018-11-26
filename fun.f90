
function fgamma(sxx,i,inter)
	  use global_cont
	implicit none
	  double precision sxx,fgamma
	  integer i,inter
	  if(i.ge.inter)then
		  call state_choose(2)
	  else
		  call state_choose(1)
	  endif

	  if (abs(sxx).le.2.d0/3*Y0) then
		  fgamma =sxx 
	  else if (sxx.ge. 2.d0/3*Y0)  then
		  fgamma=2.d0/3*Y0
		else 
		  fgamma=-2.d0/3*Y0
	endif
end function


       
