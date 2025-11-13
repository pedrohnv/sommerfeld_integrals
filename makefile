FC=gfortran
FFLAGS=-fPIC
slateclib = ./slatec/lib/libslatec.a

all: clean quadde_test mosig_michalski_PE_test

clean:
	rm -f *.o *.mod

quadde_module:
	$(FC) $(FFLAGS) -c quadde_module.f90

quadde_test: quadde_module
	$(FC) $(FFLAGS) -o test quadde_test.f90 quadde_module.o && ./test
	rm test

mosig_michalski_PE_module: quadde_module
	$(FC) $(FFLAGS) -c mosig_michalski_PE_module.f90

mosig_michalski_PE_test: quadde_module mosig_michalski_PE_module $(slateclib)
	$(FC) $(FFLAGS) -o test mosig_michalski_PE_test.f90 quadde_module.o mosig_michalski_PE_module.o $(slateclib) && ./test
	rm test

$(slateclib): slatec
	cd slatec && make FC=$(FC) && cd ..
