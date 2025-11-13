clean:
	rm -f *.o *.mod *.out

all: clean quadde_test mosig_michalski

quadde_module:
	gfortran -c quadde_module.f90

quadde_test: quadde_module
	gfortran -o a.out quadde_test.f90 quadde_module.o && ./a.out

mosig_michalski: quadde_module
	gfortran -o a.out mosig_michalski_PE.f90 quadde_module.o && ./a.out
