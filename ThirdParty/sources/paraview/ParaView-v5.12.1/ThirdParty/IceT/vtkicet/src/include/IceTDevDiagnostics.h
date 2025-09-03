/* -*- c -*- *******************************************************/
/*
 * Copyright (C) 2003 Sandia Corporation
 * Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 * the U.S. Government retains certain rights in this software.
 *
 * This source code is released under the New BSD License.
 */

#ifndef __IceTDevDiagnostics_h
#define __IceTDevDiagnostics_h

#include <IceT.h>

#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif

ICET_EXPORT void icetRaiseDiagnostic(IceTEnum type,
				     IceTBitField level,
                                     const char *file,
                                     int line,
                                     const char *format,
                                     ...);

#define icetRaiseError(type, ...)			\
    icetRaiseDiagnostic(type, ICET_DIAG_ERRORS, __FILE__, __LINE__, __VA_ARGS__)

#define icetRaiseWarning(type, ...)			\
    icetRaiseDiagnostic(type, ICET_DIAG_WARNINGS, __FILE__, __LINE__, __VA_ARGS__)

#ifdef DEBUG
#define icetRaiseDebug(...)				\
    icetRaiseDiagnostic(ICET_NO_ERROR, ICET_DIAG_DEBUG, __FILE__, __LINE__, __VA_ARGS__)
#else /* DEBUG */
#define icetRaiseDebug(...)
#endif /* DEBUG */

ICET_EXPORT void icetDebugBreak(void);

#ifdef __cplusplus
}
#endif

#endif /* __IceTDevDiagnostics_h */
