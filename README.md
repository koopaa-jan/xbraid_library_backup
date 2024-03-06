# Dynamic Resource Extensions for Xbraid
This repository contains dynamic resource extensions for the xbraid library 

## Prerequisites
- [Xbraid](https://www.xbraid.org/)
- [Open-MPI with DPP](https://gitlab.inria.fr/dynres/dyn-procs/ompi)

## Installation
The library can be installed with autotools. In this directory:

1. Run autogen
```
./autogen.sh
```
2. Run configure (optionally specify an install prefix and the paths to xbraid and mpi if they are not installed in a standard location)
```
./configure [--prefix=/path/to/install-dir] [--with-xbraid=/path/to/xbraid] [--with-mpi=/path/to/mpi]
```
3. Run make
```
make && make install
```