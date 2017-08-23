#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.65])
AC_INIT([SaddleSum], [1.5.0], [qmbp@ncbi.nlm.nih.gov])

AC_CONFIG_SRCDIR([RSaddleSum/src/RSaddleSum.c])


# Checks for programs.
AC_PROG_CC([gcc])

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lm':
AC_CHECK_LIB([m], [sqrt])

# Checks for header files.
AC_CHECK_HEADERS([inttypes.h stddef.h stdint.h stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
AC_TYPE_SIZE_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT8_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_STRTOD
AC_CHECK_FUNCS([floor memmove memset pow sqrt strtol strtoul])

AC_CONFIG_FILES([build/Makefile
                 swig/Makefile
                 Makefile
                 ])
AC_OUTPUT
