objects = constants.o readsnap.o readsampoint.o fields.o visco.o
F90     = mpifort
# FFLAGS =  -O3 -fopenmp -ffast-math -fmax-errors=4
FFLAGS = -O3
execname: $(objects)
	$(F90) $(FFLAGS) -o visco $(objects) 
%.o: %.f90
	$(F90) -c $<

clean:
	rm *.mod $(objects)
	rm visco
#	rm ./data/*
	
