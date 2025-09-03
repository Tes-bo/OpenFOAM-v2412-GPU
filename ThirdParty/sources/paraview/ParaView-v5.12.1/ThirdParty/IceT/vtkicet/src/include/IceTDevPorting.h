/* -*- c -*- *******************************************************/
/*
 * Copyright (C) 2003 Sandia Corporation
 * Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 * the U.S. Government retains certain rights in this software.
 *
 * This source code is released under the New BSD License.
 */

#ifndef __IceTDevPorting_h
#define __IceTDevPorting_h

#include <IceT.h>

#include <stdio.h>
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif

/* Returns the size of the type given by the identifier (ICET_INT, ICET_FLOAT,
   etc.)  in bytes. */
ICET_EXPORT IceTInt icetTypeWidth(IceTEnum type);

#ifdef WIN32
#define strncpy(dest, src, size) strncpy_s(dest, size, src, _TRUNCATE)
#define fdopen _fdopen
#endif

/* A portable version of getenv. This version requires you to give a buffer
   to result the result into. Returns ICET_TRUE if the environment variable
   was found.*/
ICET_EXPORT IceTBoolean icetGetEnv(const char *variable_name,
                                   char *buffer,
                                   IceTSizeType buffer_size);

/* A portable version of putenv. This version takes the variable name and
   value as separate arguments. */
ICET_EXPORT void icetPutEnv(const char *name, const char *value);

/* A portable version of snprintf. The behavior might not be perfectly
   consistent across platforms. */
ICET_EXPORT IceTSizeType icetSnprintf(char *buffer, IceTSizeType size,
                                      const char *format, ...);

#ifdef __cplusplus
}
#endif

#endif /*__IceTDevPorting_h*/
