Objects = global.o bound_choose.o bound_lax_problem.o  bound_sod_problem.o bound_two_balst_waves.o boundreflect.o CFL.o Characterreconst.o fun.o init_boundreflect.o init_lax_problem.o init_shu.o init_sod_problem.o init_two_blast.o initial.o LF_con_var.o Main.o output.o output1.o R_K.o space.o splitting.o timesolve.o tridiagsolve.o u_f.o upwind.o WENO.o weno5.o WENO5_NEW.o wenonew.o weno_new.o init_interface.o bound_interface.o
source = global.f90 bound_choose.f90 bound_lax_problem.f90  bound_sod_problem.f90 bound_two_balst_waves.f90 boundreflect.f90 CFL.f90 Characterreconst.f90 fun.f90 init_boundreflect.f90 init_lax_problem.f90 init_shu.f90 init_sod_problem.f90 init_two_blast.f90 initial.f90 LF_con_var.f90 Main.f90 output.f90 output1.f90 R_K.f90 space.f90 splitting.f90 timesolve.f90 tridiagsolve.f90 u_f.f90 upwind.f90 WENO.f90 weno5.f90 WENO5_NEW.f90  wenonew.f90 weno_new.f90 init_interface.f90 bound_interface.f90
Bin = main.exe
F90 = gfortran -openmp

${Bin}:${Objects}
	${F90} -o ${Bin} ${Objects}

%.o :%.f90
	${F90} -c -o $@ $<


.PHONY:clean
clean:
	rm -f main.exe $(Objects) *.mod

