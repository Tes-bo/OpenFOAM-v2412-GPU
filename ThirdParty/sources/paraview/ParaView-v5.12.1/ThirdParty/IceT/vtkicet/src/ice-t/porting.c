/* -*- c -*- *******************************************************/
/*
 * Copyright (C) 2003 Sandia Corporation
 * Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 * the U.S. Government retains certain rights in this software.
 *
 * This source code is released under the New BSD License.
 */

/* Turn on non-POSIX features of GNU. Many of these are widely available. */
#define _GNU_SOURCE

#include <IceTDevPorting.h>

#include <IceT.h>

#include <IceTDevDiagnostics.h>

#include <stdlib.h>
#include <stdarg.h>

#ifndef WIN32
#include <sys/time.h>
#else
#include <windows.h>
#include <winbase.h>
#endif

#ifndef WIN32
double icetWallTime(void)
{
    static struct timeval start = { 0, 0 };
    struct timeval now;
    struct timeval *tp;

  /* Make the first call to icetWallTime happen at second 0.  This should
     allow for more significant bits in the microseconds. */
    if (start.tv_sec == 0) {
        tp = &start;
    } else {
        tp = &now;
    }

    gettimeofday(tp, NULL);

    return (tp->tv_sec - start.tv_sec) + 0.000001*(double)tp->tv_usec;
}
#else /*WIN32*/
double icetWallTime(void)
{
    static DWORD start = 0;

    if (start == 0) {
        start = GetTickCount();
        return 0.0;
    } else {
        DWORD now = GetTickCount();
        return 0.001*(now-start);
    }
}
#endif /*WIN32*/

IceTInt icetTypeWidth(IceTEnum type)
{
    switch (type) {
      case ICET_BOOLEAN:
          return sizeof(IceTBoolean);
      case ICET_BYTE:
          return sizeof(IceTByte);
      case ICET_SHORT:
          return sizeof(IceTShort);
      case ICET_INT:
          return sizeof(IceTInt);
      case ICET_FLOAT:
          return sizeof(IceTFloat);
      case ICET_DOUBLE:
          return sizeof(IceTDouble);
      case ICET_POINTER:
          return sizeof(IceTVoid *);
      case ICET_VOID:
          return 1;
      case ICET_NULL:
          return 0;
      default:
          icetRaiseError(ICET_INVALID_VALUE, "Bad type identifier %d.", type);
    }

    return 0;
}

#ifndef WIN32
IceTBoolean icetGetEnv(const char *variable_name,
                       char *buffer,
                       IceTSizeType buffer_size)
{
    const char *value = getenv(variable_name);
    if (value != NULL)
    {
        strncpy(buffer, value, buffer_size);
        return ICET_TRUE;
    } else {
        buffer[0] = '\0';
        return ICET_FALSE;
    }
}

void icetPutEnv(const char *name, const char *value)
{
    setenv(name, value, ICET_TRUE);
}
#else /*WIN32*/
IceTBoolean icetGetEnv(const char *variable_name,
                       char *buffer,
                       IceTSizeType buffer_size)
{
    size_t required_size; /* Not really used. */
    errno_t error;

    error = getenv_s(&required_size, buffer, buffer_size, variable_name);
    if ((error == 0) && (required_size > 0))
    {
        return ICET_TRUE;
    } else {
        buffer[0] = '\0';
        return ICET_FALSE;
    }
}

void icetPutEnv(const char *name, const char *value)
{
    _putenv_s(name, value);
}
#endif /*WIN32*/

ICET_EXPORT IceTSizeType icetSnprintf(char *buffer, IceTSizeType size,
                                      const char *format, ...)
{
    va_list format_args;
    IceTSizeType num_written;

    va_start(format_args, format);
#ifdef WIN32
    num_written = _vsnprintf_s(buffer, size, _TRUNCATE, format, format_args);
#else
    num_written = vsnprintf(buffer, size, format, format_args);
#endif
    va_end(format_args);

    return num_written;
}
