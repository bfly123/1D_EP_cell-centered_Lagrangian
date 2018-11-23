subroutine bound(u4)
use global
implicit none
double precision u4(-nv:jx+nv,0:3)

select case(kind_problem)

    case(1)
        call bound_Wilkins_problem(u4) 
    case(2)
        call bound_Piston_problem(u4)        
    case(4)
        call bound_accuracy_test(u4)
    case(5)
        call bound_accuracy_test(u4)
   !     call bound_boundreflect(u4,X4)
   ! case(6)
   !     call bound_interface(u4,X4)
end select

end

subroutine bound_v(u1)
	use global
	implicit none 
	double precision u1(-nv:jx+nv)
	integer i

	select case(kind_problem)

	case(4,5)
		do i = -nv,-1
		u1(i) = u1(jx+i)
		enddo
		do i=1,nv
		u1(jx+i) = u1(i)
		enddo
	end select
	end



