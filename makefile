slateclib = ./slatec/lib/libslatec.a

all: clean quadde_test mosig_michalski

clean:
	rm -f *.o *.mod *.out

quadde_module:
	gfortran -fPIC -c quadde_module.f90

quadde_test: quadde_module
	gfortran -o a.out quadde_test.f90 quadde_module.o && ./a.out

mosig_michalski: quadde_module
	gfortran -g -o a.out mosig_michalski_PE.f90 quadde_module.o $(slateclib) && ./a.out

slatec:
	cd slatec && make FC=gfortran && cd ..
