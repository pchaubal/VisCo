objects = constants.o readsnap.o readsampoint.o utils.o fields.o visco.o
F90     = mpifort
execname: $(objects)
	$(F90) -o visco $(objects) 
%.o: %.f90
	$(F90) -c $<

clean:
	rm *.mod $(objects)
	rm visco
#	rm ./data/*
	
