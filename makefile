Objects = global.o fun.o bound_choose.o  bound_sod_problem.o CFL.o  init_sod_problem.o  initial.o  Main.o output.o  R_K.o space.o  timesolve.o upwind.o HLLC_EP.o  WENO5_NEW.o

source = global.f90 fun.f90 bound_choose.f90   bound_sod_problem.f90  CFL.f90  init_sod_problem.f90  initial.f90  Main.f90 output.f90  R_K.f90 space.f90 timesolve.f90  upwind.f90  HLLC_EP.f90 WENO5_NEW.f90
Bin = main.exe
F90 = gfortran

${Bin}:${Objects}
	${F90} -o ${Bin} ${Objects}

%.o :%.f90
	${F90} -c -o $@ $<


.PHONY:clean
clean:
	rm -f main.exe $(Objects) *.mod

