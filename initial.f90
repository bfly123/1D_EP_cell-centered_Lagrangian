subroutine init
use global
implicit none 

!**********kind of problem*************
!**1 double mach reflection 
!     
    !case 1 Wilkins' problem
    !case 2 Piston problem
    !case 3 lax_problem
    !case 4 two_blast
  kind_problem =6
  select case(kind_problem) 
    
    case(1)
        call init_Wilkins_problem
    case(2)
        call init_Piston_problem 
    case(4)
       call init_accuracy_test
    case(5)
        call init_entropy_wave
    case(6)
        call init_two_material
  end select   
     

    
end
       
    
