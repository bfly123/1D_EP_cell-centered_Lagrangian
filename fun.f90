
function fgamma(sxx)
	  use global_cont
	implicit none
	  double precision sxx,fgamma
	  if (abs(sxx).le.2.d0/3*Y0) then
		  fgamma =sxx 
	  else if (sxx.ge. 2.d0/3*Y0)  then
		  fgamma=2.d0/3*Y0
		else 
		  fgamma=-2.d0/3*Y0
	endif
end function


       
