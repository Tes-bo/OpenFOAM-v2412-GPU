/* -*- c -*- *******************************************************/
/*
 * Copyright (C) 2003 Sandia Corporation
 * Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 * the U.S. Government retains certain rights in this software.
 *
 * This source code is released under the New BSD License.
 */

#include <IceTDevDiagnostics.h>

#include <IceT.h>

#include <IceTDevCommunication.h>
#include <IceTDevContext.h>
#include <IceTDevPorting.h>

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>

#include <signal.h>

#ifndef WIN32
#include <unistd.h>
#endif

#define MAX_MESSAGE_LEN 1024

static IceTEnum currentError = ICET_NO_ERROR;
static IceTEnum currentLevel;

void icetRaiseDiagnostic(IceTEnum type,
                         IceTBitField level,
                         const char *file,
                         int line,
                         const char *format,
                         ...)
{
#define ICET_MESSAGE_SIZE 1024
    static int raisingDiagnostic = 0;
    IceTBitField diagLevel;
    static char full_message[ICET_MESSAGE_SIZE+1];
    IceTSizeType offset;
    int rank;
    va_list format_args;

#define FINISH        raisingDiagnostic = 0; return

    if (raisingDiagnostic) {
        printf("PANIC:%s:%d: diagnostic raised while rasing diagnostic!\n",
               file, line);
        icetStateDump();
#ifdef DEBUG
        icetDebugBreak();
#endif
        return;
    }
    if (icetGetState() == NULL) {
        printf("PANIC:%s:%d: diagnostic raised when no context was current!\n",
               file, line);
#ifdef DEBUG
        icetDebugBreak();
#endif
        return;
    }
    raisingDiagnostic = 1;
    full_message[0] = '\0';
    offset = 0;
#define ADD_TO_MESSAGE(...)                                                 \
    offset += icetSnprintf(full_message + offset,                           \
                      ICET_MESSAGE_SIZE - offset,                           \
                      __VA_ARGS__);                                         \
    if (offset > ICET_MESSAGE_SIZE) {                                       \
        printf("PANIC:%s:%d: Diagnostic message too large!\n", file, line); \
        icetStateDump();                                                    \
        return;                                                             \
    }

    if ((currentError == ICET_NO_ERROR) || (level < currentLevel)) {
        currentError = type;
        currentLevel = level;
    }
    diagLevel = icetUnsafeStateGetInteger(ICET_DIAGNOSTIC_LEVEL)[0];
    if ((diagLevel & level) != level) {
      /* Don't do anything if we are not reporting the raised diagnostic. */
        FINISH;
    }

    rank = icetCommRank();
    if ((diagLevel & ICET_DIAG_ALL_NODES) != 0) {
        /* Reporting on all nodes. */
        ADD_TO_MESSAGE("ICET,%d:", rank);
    } else if (rank == 0) {
        /* Rank 0 is lone reporter. */
        ADD_TO_MESSAGE("ICET:");
    } else {
      /* Not reporting because not rank 0. */
        FINISH;
    }
    /* Add description of diagnostic type. */
    switch (level & 0xFF) {
      case ICET_DIAG_ERRORS:
          ADD_TO_MESSAGE("ERROR:");
          break;
      case ICET_DIAG_WARNINGS:
          ADD_TO_MESSAGE("WARNING:");
          break;
      case ICET_DIAG_DEBUG:
          ADD_TO_MESSAGE("DEBUG:");
          break;
    }
#ifdef DEBUG
    ADD_TO_MESSAGE("%s:%d:", file, line);
#else
    /* shut up warnings */
    (void)file;
    (void)line;
#endif

    ADD_TO_MESSAGE(" %s\n", format);
    va_start(format_args, format);
    vprintf(full_message, format_args);
    va_end(format_args);
    fflush(stdout);

#ifdef DEBUG
    if ((level & 0xFF) == ICET_DIAG_ERRORS) {
        icetDebugBreak();
    }
#endif

    FINISH;
}

IceTEnum icetGetError(void)
{
    IceTEnum error = currentError;
    currentError = ICET_NO_ERROR;
    return error;
}

void icetDiagnostics(IceTBitField mask)
{
    icetStateSetInteger(ICET_DIAGNOSTIC_LEVEL, mask);
}

void icetDebugBreak(void)
{
#if 0
    printf("Waiting for debugger in process %d\n", getpid());
    sleep(100);
#endif
    raise(SIGSEGV);
}
