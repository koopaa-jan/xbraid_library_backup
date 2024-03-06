#!/bin/sh

echo "*******************************************************"
echo "Let's build the dynamic resource extensions for xbraid!"
echo "*******************************************************"
echo ""

autoreconf -i

echo "autogen: SUCCESS! You can now run ./configure [--prefix=/path/to/install-dir] [--with-xbraid=/path/to/xbraid] [--with-mpi=/path/to/mpi]"