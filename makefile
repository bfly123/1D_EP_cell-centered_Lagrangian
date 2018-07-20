Objects = global.o fun.o bound_choose.o  bound_Wilkins_problem.o  bound_Piston_problem.o CFL.o  init_Wilkins_problem.o  init_Piston_problem.o initial.o  Main.o output.o  R_K.o space.o  timesolve.o upwind.o HLLC_EP.o  WENO5_NEW.o WENO5_NEW_LF.o LF_splitting.o eigen_var.o  reverse.o output1.o WENO3_new.o upwind3.o

source = global.f90 fun.f90 bound_choose.f90   bound_Wilkins_problem.f90  bound_Piston_problem.f90  CFL.f90  init_Wilkins_problem.f90 init_Piston_problem.f90  initial.f90  Main.f90 output.f90  R_K.f90 space.f90 timesolve.f90  upwind.f90  HLLC_EP.f90 WENO5_NEW.f90 WENO5_NEW_LF.f90 LF_splitting.f90 eigen_var.f90 reverse.f90 output1.f90WENO3_new.f90 upwind3.f90
Bin = main.exe
F90 = ifort -mkl
${Bin}:${Objects}
	${F90} -o ${Bin} ${Objects}

%.o :%.f90
	${F90} -c -o $@ $<


.PHONY:clean
clean:
	rm -f main.exe $(Objects) *.mod

