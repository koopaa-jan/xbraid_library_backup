# Variables
CC = mpicc
CFLAGS = -Wall -fPIC
LIBPATHS=
INCLUDES=
LIBS = -lm -lxbraid_dyn -lmpi

# Targets
all: ex-01_dyn ex-heat_equation ex-heat_equation_2D 

ex-01_dyn: ex-01_dyn.o
	$(CC) $(CFLAGS) $(INCLUDES) ex-01_dyn.o -o $@ $(LIBPATHS) $(LIBS)

ex-01_dyn.o: ex-01_dyn.c
	$(CC) -c $(CFLAGS) $(INCLUDES) $< -o $@

ex-heat_equation: ex-heat_equation.o
	$(CC) $(CFLAGS) $(INCLUDES) ex-heat_equation.o -o $@ $(LIBPATHS) -Wl,-rpath,$(LIBPATHS) $(LIBS)

ex-heat_equation.o: ex-heat_equation.c
	$(CC) -c $(CFLAGS) $(INCLUDES) $< -o $@

ex-heat_equation_2D: ex-heat_equation_2D.o
	$(CC) $(CFLAGS) $(INCLUDES) ex-heat_equation_2D.o -o $@ $(LIBPATHS) -Wl,-rpath,$(LIBPATH) $(LIBS)

ex-heat_equation_2D.o: ex-heat_equation_2D.c
	$(CC) -c $(CFLAGS) $(INCLUDES) $< -o $@


clean:
	rm -f *.o *.out.* ex-01_dyn ex-heat_equation ex-heat_equation_2D 