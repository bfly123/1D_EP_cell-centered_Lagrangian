subroutine init
use global
implicit none 

!**********kind of problem*************
!**1 double mach reflection 
!     
    kind_problem=6 
    !case 1 sod_problem
    !case 2 shu_osher
    !case 3 lax_problem
    !case 4 two_blast
    
    kind_split=2
      
  select case(kind_problem) 
    
    case(1)
        call init_sod_problem
    case(2)
        call init_shu_osher_problem 
    case(3)
        call init_lax_problem
    case(4)
       call init_two_blast_waves
    case(5)
        call init_boundreflect
    case(6)
        call init_interface !接触间断
  end select   
     
   
    

    
end
       
    
