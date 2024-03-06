# Variables
CC = @CC@
CFLAGS = -Wall -fPIC $(AM_CFLAGS)
LIBS = -lm $(AM_LDFLAGS) -lmpi -Wl,-Bstatic -lbraid -Wl,-Bdynamic

# Targets
all: src/build/libxbraid_dyn.so ex-01_dyn ex-heat_equation

src/build/libxbraid_dyn.so: src/build/braid_dyn.o src/build/dmr_dyn.o src/build/uvector_dyn.o
	$(CC) -shared -o $@ $^ $(LIBS)

src/build/braid_dyn.o: src/braid_dyn.c src/braid_dyn.h src/_braid_dyn.h
	$(CC) -c $(CFLAGS) $< -o $@

src/build/dmr_dyn.o: src/dmr_dyn.c src/dmr_dyn.h
	$(CC) -c $(CFLAGS) $< -o $@

src/build/uvector_dyn.o: src/uvector_dyn.c src/_braid_dyn.h
	$(CC) -c $(CFLAGS) $< -o $@

ex-01_dyn: src/ex-01_dyn.o
	$(CC) $(CFLAGS) $(INCLUDES) src/ex-01_dyn.o -o $@ $(LIBS) -L/opt/hpc/build/xbraid_library_dyn/src/build -lxbraid_dyn -Wl,-rpath,/opt/hpc/build/xbraid_library_dyn/src/build

src/ex-01_dyn.o: src/ex-01_dyn.c
	$(CC) -c $(CFLAGS) $(INCLUDES) $< -o $@

ex-heat_equation: src/ex-heat_equation.o
	$(CC) $(CFLAGS) $(INCLUDES) src/ex-heat_equation.o -o $@ $(LIBS) -L/opt/hpc/build/xbraid_library_dyn/src/build -lxbraid_dyn -Wl,-rpath,/opt/hpc/build/xbraid_library_dyn/src/build

src/ex-heat_equation.o: src/ex-heat_equation.c
	$(CC) -c $(CFLAGS) $(INCLUDES) $< -o $@

install:
	install -d $(DESTDIR)$(prefix)/lib
	install -m 755 src/build/libxbraid_dyn.so $(DESTDIR)$(prefix)/lib
	install -d $(DESTDIR)$(prefix)/include
	install -m 644 src/braid_dyn.h $(DESTDIR)$(prefix)/include
	install -m 644 src/_braid_dyn.h $(DESTDIR)$(prefix)/include
	install -m 644 src/dmr_dyn.h $(DESTDIR)$(prefix)/include
	

clean:
	rm -f src/build/*.o src/build/libxbraid_dyn.so src/ex-01_dyn.o ex-01_dyn ex-01.out.* ex-heat_equation ex-h_e.out.*