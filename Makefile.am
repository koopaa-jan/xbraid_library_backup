# Variables
CC = @CC@
CFLAGS = -Wall -fPIC $(AM_CFLAGS)
LIBS = $(AM_LDFLAGS) -lm -lmpi -ldmr -Wl,-Bstatic -lbraid -Wl,-Bdynamic

# Targets
all: src/build/libxbraid_dyn.so  

src/build/libxbraid_dyn.so: src/build/braid_dyn.o src/build/uvector_dyn.o
	$(CC) -shared -o $@ $^ $(LIBS)

src/build/braid_dyn.o: src/braid_dyn.c src/braid_dyn.h src/_braid_dyn.h
	$(CC) -c $(CFLAGS) $< -o $@

src/build/uvector_dyn.o: src/uvector_dyn.c src/_braid_dyn.h
	$(CC) -c $(CFLAGS) $< -o $@


install:
	install -d $(DESTDIR)$(prefix)/lib
	install -m 755 src/build/libxbraid_dyn.so $(DESTDIR)$(prefix)/lib
	install -d $(DESTDIR)$(prefix)/include
	install -m 644 src/braid_dyn.h $(DESTDIR)$(prefix)/include
	install -m 644 src/_braid_dyn.h $(DESTDIR)$(prefix)/include
	

clean:
	rm -f src/build/*.o src/build/libxbraid_dyn.so 