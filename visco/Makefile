objects = constants.o readsnap.o fields.o visco.o
F90     = mpifort
execname: $(objects)
	$(F90) -o execute $(objects) 
%.o: %.f90
	$(F90) -c $<

clean:
	rm *.mod $(objects)
	rm execute
#	rm ./data/*
	
