FC=gfortran
FFLAGS=-fPIC
slateclib = ./slatec/lib/libslatec.a

all: clean quadde_test mosig_michalski

clean:
	rm -f *.o *.mod *.out

quadde_module:
	$(FC) $(FFLAGS) -c quadde_module.f90

quadde_test: quadde_module
	$(FC) $(FFLAGS) -o a.out quadde_test.f90 quadde_module.o && ./a.out

mosig_michalski: quadde_module $(slateclib)
	$(FC) $(FFLAGS) -o a.out mosig_michalski_PE.f90 quadde_module.o $(slateclib) && ./a.out

$(slateclib): slatec
	cd slatec && make FC=$(FC) && cd ..
