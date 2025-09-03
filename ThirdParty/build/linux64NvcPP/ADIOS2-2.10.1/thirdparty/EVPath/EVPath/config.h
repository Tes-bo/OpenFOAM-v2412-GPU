#ifndef _CM_CONFIG_H
#define _CM_CONFIG_H

/* config.h.in.  Generated from configure.ac by autoheader.  */

/* Define to 1 if you have the <enet/enet.h> header file. */
/* #undef HAVE_ENET_ENET_H */

/* Define to desired port range Low:High */
#define EVPATH_DEFAULT_PORT_RANGE "ANY"

/* Define to 1 if you have the enet header and libraries */
#define ENET_FOUND

/* Define to 1 to if ZPL ENET transport was built */
/* #undef ZPL_ENET_AVAILABLE */

/* Define to 1 if you have the udt4 header and libraries */
/* #undef UDT4_FOUND */

/* Define to 1 if you have the IB header and libraries */
/* #undef IB_FOUND */

/* Define to 1 if you have the LIBFABRIC header and libraries */
/* #undef LIBFABRIC_FOUND */

/* Define to 1 if you have the df_shm header and libraries */
/* #undef DF_SHM_FOUND */

/* Place where evpath transport modules are installed */
#define EVPATH_MODULE_INSTALL_DIR "/home/sibo/TesboCFD/ThirdParty-v2412/platforms/linux64NvcPP/ADIOS2-2.10.1/lib64/adios2-evpath-modules-2_10"

/* Place where evpath tests are to be installed */
/* #undef EVPATH_TEST_INSTALL_DIR */

/* Place where the ssh binary can be found */
#define SSH_PATH "/usr/bin/ssh"

/* Define to 1 if you have the <dlfcn.h> header file. */
/* #undef HAVE_DLFCN_H */

/* Define if fd_set has element fds_bits */
/* #undef HAVE_FDS_BITS */

/* Define to 1 if you have the `getdomainname' function. */
#define HAVE_GETDOMAINNAME

/* Define to 1 if you have the `getloadavg' function. */
#define HAVE_GETLOADAVG

/* Define to 1 if you have the `gettimeofday' function. */
#define HAVE_GETTIMEOFDAY

/* Define to 1 if you have the `pthread' library (-lpthread). */
/* #undef HAVE_LIBPTHREAD */

/* Uses the Mac OS X style conventions for sysctl */
/* #undef HAVE_MAC_SYSCTL */

/* Define if system uses LINUX-style sysctl */
/* #undef HAVE_LINUX_SYSCTL */

/* Define to 1 if you have the <malloc.h> header file. */
#define HAVE_MALLOC_H

/* Define to 1 if you have the <math.h> header file. */
/* #undef HAVE_MATH_H */

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H

/* Define to 1 if you have the <netdb.h> header file. */
#define HAVE_NETDB_H

/* Define to 1 if you have the <stdarg.h> header file. */
/* #undef HAVE_STDARG_H */

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H

/* Define to 1 if you have the <strings.h> header file. */
/* #undef HAVE_STRINGS_H */

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H

/* Define to 1 if you have the `sysconf' function. */
#define HAVE_SYSCONF

/* Define to 1 if you have the `sysinfo' function. */
#define HAVE_SYSINFO

/* Define to 1 if you have the <sys/epoll.h> header file. */
#define HAVE_SYS_EPOLL_H

/* Define to 1 if you have the <cod.h> header file. */
#define HAVE_COD_H

/* Define to 1 if you have the <sys/select.h> header file. */
#define HAVE_SYS_SELECT_H

/* Define to 1 if you have the <sys/sockio.h> header file. */
/* #undef HAVE_SYS_SOCKIO_H */

/* Define to 1 if you have the <sys/stat.h> header file. */
/* #undef HAVE_SYS_STAT_H */

/* Define to 1 if you have the <sys/sysctl.h> header file. */
/* #undef HAVE_SYS_SYSCTL_H */

/* Define to 1 if you have the <sys/times.h> header file. */
#define HAVE_SYS_TIMES_H

/* Define to 1 if you have the <sys/time.h> header file. */
#define HAVE_SYS_TIME_H

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H

/* Define to 1 if you have the <sys/uio.h> header file. */
#define HAVE_SYS_UIO_H

/* Define to 1 if you have the <sys/un.h> header file. */
#define HAVE_SYS_UN_H

/* Define to 1 if you have the `uname' function. */
#define HAVE_UNAME

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H

/* Define to 1 if you have the <windows.h> header file. */
/* #undef HAVE_WINDOWS_H */

/* Define to 1 if you have the <winsock2.h> header file. */
/* #undef HAVE_WINSOCK2_H */

/* The size of `int', as computed by sizeof. */
#define SIZEOF_INT 4

/* The size of `long', as computed by sizeof. */
#define SIZEOF_LONG 8

/* Define to 1 if you have the `writev' function. */
#define HAVE_WRITEV

/* Define to 1 if you have the `getifaddrs' function. */
#define HAVE_GETIFADDRS

/* Define to 1 if you have the `clock_gettime' function. */
#define HAVE_CLOCK_GETTIME

/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS

/* Define this if Pthreads should be used for running tests */
#define USE_PTHREADS

/* Version number of package */
/* #undef VERSION */

/* Define so that glibc/gnulib argp.h does not typedef error_t. */
/* #undef __error_t_defined */

/* Define to empty if `const' does not conform to ANSI C. */
/* #undef const */

/* Define to a type to use for `error_t' if it is not otherwise available. */
/* #undef error_t */

/* Define to `int' if <sys/types.h> does not define. */
/* #undef pid_t */

/* Define to `int' if <sys/types.h> does not define. */
/* #undef pid_t */

/* Define to the shared module suffix in use */
#define CMAKE_SHARED_MODULE_SUFFIX ".so"

/* Set to 1 to build without dynamic linking  */
/* #undef NO_DYNAMIC_LINKING */

/* Set to 1 if NNTI libraries and include file are found */
/* #undef NNTI_FOUND */

/* Set to 1 if NVML libraries and include file are found */
/* #undef NVML_FOUND */

/* Set to 1 if CM should default to CMSelfFormats */
#define CM_SELF_FORMATS 1

/* Set to 1 if Cmake found MPI for C */
#define MPI_C_FOUND

/* Define if byteorder is bigendian */
/* #undef WORDS_BIGENDIAN */

#define CM_DEFAULT_TRANSPORT "sockets"

#define CM_LIBRARY_PREFIX "adios2_"

#define IPCONFIG_ENVVAR_PREFIX "ADIOS2_"

#endif
