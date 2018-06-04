subroutine bound(u4)
use global
implicit none
double precision u4(-nv:jx+nv,0:2)

select case(kind_problem)

    case(1)
        call bound_sod_problem(u4) 
    case(2)
        call bound_sod_problem(u4)        
!        call bound_periodic_boundary_condition(u4)
    case(3)
        call bound_lax_problem(u4)
    case(4)
        call bound_two_blast_waves(u4)
    case(5)
        call bound_boundreflect(u4)
    case(6)
        call bound_interface(u4)
end select

end
