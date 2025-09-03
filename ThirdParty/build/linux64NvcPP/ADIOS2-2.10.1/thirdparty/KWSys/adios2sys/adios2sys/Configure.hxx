/* Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
   file Copyright.txt or https://cmake.org/licensing#kwsys for details.  */
#ifndef adios2sys_Configure_hxx
#define adios2sys_Configure_hxx

/* Include C configuration.  */
#include <adios2sys/Configure.h>

/* Whether wstring is available.  */
#define adios2sys_STL_HAS_WSTRING 1
/* Whether <ext/stdio_filebuf.h> is available. */
#define adios2sys_CXX_HAS_EXT_STDIO_FILEBUF_H                         \
  1
/* Whether the translation map is available or not. */
#define adios2sys_SYSTEMTOOLS_USE_TRANSLATION_MAP                     \
  1

#if defined(__SUNPRO_CC) && __SUNPRO_CC > 0x5130 && defined(__has_attribute)
#  define adios2sys_has_cpp_attribute(x) __has_attribute(x)
#elif defined(__has_cpp_attribute)
#  define adios2sys_has_cpp_attribute(x) __has_cpp_attribute(x)
#else
#  define adios2sys_has_cpp_attribute(x) 0
#endif

#if __cplusplus >= 201103L
#  define adios2sys_NULLPTR nullptr
#else
#  define adios2sys_NULLPTR 0
#endif

#ifndef adios2sys_FALLTHROUGH
#  if __cplusplus >= 201703L &&                                               \
    adios2sys_has_cpp_attribute(fallthrough)
#    define adios2sys_FALLTHROUGH [[fallthrough]]
#  elif __cplusplus >= 201103L &&                                             \
    adios2sys_has_cpp_attribute(gnu::fallthrough)
#    define adios2sys_FALLTHROUGH [[gnu::fallthrough]]
#  elif __cplusplus >= 201103L &&                                             \
    adios2sys_has_cpp_attribute(clang::fallthrough)
#    define adios2sys_FALLTHROUGH [[clang::fallthrough]]
#  endif
#endif
#ifndef adios2sys_FALLTHROUGH
#  define adios2sys_FALLTHROUGH static_cast<void>(0)
#endif

#undef adios2sys_has_cpp_attribute

/* If building a C++ file in kwsys itself, give the source file
   access to the macros without a configured namespace.  */
#if defined(KWSYS_NAMESPACE)
#  if !adios2sys_NAME_IS_KWSYS
#    define kwsys adios2sys
#  endif
#  define KWSYS_NAME_IS_KWSYS adios2sys_NAME_IS_KWSYS
#  define KWSYS_STL_HAS_WSTRING adios2sys_STL_HAS_WSTRING
#  define KWSYS_CXX_HAS_EXT_STDIO_FILEBUF_H                                   \
    adios2sys_CXX_HAS_EXT_STDIO_FILEBUF_H
#  define KWSYS_FALLTHROUGH adios2sys_FALLTHROUGH
#  define KWSYS_SYSTEMTOOLS_USE_TRANSLATION_MAP                               \
    adios2sys_SYSTEMTOOLS_USE_TRANSLATION_MAP
#endif

#endif
