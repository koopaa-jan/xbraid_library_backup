#!/bin/sh

AC_INIT([xbraid_dyn], [1.0], [jan.grimm@tum.de])
AM_INIT_AUTOMAKE
AC_PROG_CC

# Add the --with-xbraid option
AC_ARG_WITH([xbraid],
  [AS_HELP_STRING([--with-xbraid],
    [Specify the path to xbraid library and include directory])],
  [XBRAID_PATH="$withval"],
  [XBRAID_PATH=""])

# Check if the xbraid path is provided
if test "x$XBRAID_PATH" != "x"; then
  AC_MSG_NOTICE([Using user-specified xbraid path: $XBRAID_PATH])
else
  AC_MSG_NOTICE([No xbraid path specified, using default])
  # Set a default path or behavior if not specified
  XBRAID_PATH="/usr/local"
fi

# Add the specified xbraid library and include paths to AM_LDFLAGS and AM_CFLAGS
AM_LDFLAGS="-L$XBRAID_PATH $AM_LDFLAGS"
AM_CFLAGS="-I$XBRAID_PATH $AM_CFLAGS"
AC_SUBST([XBRAID_PATH])

# Add the --with-mpi option
AC_ARG_WITH([mpi],
  [AS_HELP_STRING([--with-mpi],
    [Specify the path to mpi library and include directory])],
  [MPI_PATH="$withval"],
  [MPI_PATH=""])

# Check if the xbraid path is provided
if test "x$MPI_PATH" != "x"; then
  AC_MSG_NOTICE([Using user-specified mpi path: $MPI_PATH])
else
  AC_MSG_NOTICE([No mpi path specified, using default])
  # Set a default path or behavior if not specified
  MPI_PATH="/usr/local"
fi

# Add the specified mpi library and include paths to AM_LDFLAGS and AM_CFLAGS
AM_LDFLAGS="-L$MPI_PATH/lib $AM_LDFLAGS"
AM_CFLAGS="-I$MPI_PATH/include $AM_CFLAGS"
AC_SUBST([MPI_PATH])
AC_SUBST([AM_CFLAGS])
AC_SUBST([AM_LDFLAGS])



# Define the installation prefix
AC_PREFIX_DEFAULT([/usr/local])

# Checks for libraries
AC_CHECK_LIB([mpi], [MPI_Session_dyn_v2a_psetop])
AC_CHECK_LIB([xbraid], [braid_SetPrintLevel])

# Output Makefile
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
