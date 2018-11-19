subroutine  VL_splitting(uo,UL,UR,Uh)
	  implicit none
	  double precision uL(0:3)
	  double precision uR(0:3)
	  double precision uH(0:3)
	  double precision uo(0:3)
	  double precision c,uu

	  uu=uo(1)
	  call sound(uo,c)

	  if(uu.ge.0)then
		uh(0) = uL(0)
		uh(1) = uL(1)
	else
		uh(0) = uR(0)
		uh(1) = ur(1)
	endif

		if(uu-c.ge.0)then
		uh(2) = uL(2)
		else
		uh(2) = ur(2)
		endif

		if(uu+c.ge.0)then
		uh(3) = uL(3)
		else
		uh(3) = uR(3)
		endif

		end
		



