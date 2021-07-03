#
#                    Define compiler and options for compilation
#FC = ifort       
FC = gfortran
#FFLAGS = -ffast -w -convert big_endian  -Vaxlib
#LDFLAGS = -ffast -w  -convert big_endian  -Vaxlib
FFLAGS = -O2 -w -fconvert=big-endian  
LDFLAGS = -O2 -w  -fconvert=big-endian  

#
BIN = nbody1.o nbodyaux.o integral.o

nbody1 : nbody1.h nbody1.o nbodyaux.o
#                                         Compile nbody1 
	$(FC) $(LDFLAGS) -o  nbody1.exe nbody1.o nbodyaux.o
	size nbody1.exe

equilibrium : equilibrium.o nbody1.h  nbodyaux.o
	$(FC) $(LDFLAGS) -o  equilibrium.exe equilibrium.o nbodyaux.o
clean:                                 # remove object files
	rm  *.exe *.o
#      Implicit rules for compilation
.f.o: 
	$(FC) -c $(FFLAGS) $<


