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
  kind_problem =2  
      
  select case(kind_problem) 
    
    case(1)
        call init_Wilkins_problem
    case(2)
        call init_Piston_problem 
    !case(3)
    !    call init_lax_problem
    !case(4)
    !   call init_two_blast_waves
    !case(5)
    !    call init_boundreflect
    !case(6)
    !    call init_interface !接触间断
  end select   
     
   
    

    
end
       
    
