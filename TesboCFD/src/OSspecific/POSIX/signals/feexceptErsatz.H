// ============================================================================
// Generic handling of floating point exceptions
//
// 2018 Alexey Matveichev
// 2021 Tatsuya Shimizu - ARM64 support
// 2023 Guanyang Xue - improve ARM64 support
// ----------------------------------------------------------------------------
//
// Original Author
//     David N. Williams
//
// License
//     Public Domain
//
// Description
//     Alternate, nondefault handling of IEEE 754 floating-point exceptions
//     in OS X, based on the GNU functions feenableexcept(), fedisableeexcept()
//     and fegetexcept() [in libm]
//
// http://www-personal.umich.edu/~williams/archive/computation/fe-handling-example.c
// ============================================================================

#ifndef Foam_feexceptErsatz_H
#define Foam_feexceptErsatz_H

#ifdef __APPLE__
#include <fenv.h>

// Workaround for Apple Silicon - has SIGILL (illegal instruction) not SIGFPE
#if defined __arm64__
#undef  SIGFPE
#define SIGFPE SIGILL
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline int feenableexcept(unsigned int excepts)
{
    static fenv_t fenv;
    unsigned int new_excepts = excepts & FE_ALL_EXCEPT;
    unsigned int old_excepts;   // previous masks

    if (fegetenv(&fenv))
    {
        return -1;
    }

#if defined __arm64__
    old_excepts = fenv.__fpsr & FE_ALL_EXCEPT;

    fenv.__fpsr |= new_excepts;
    fenv.__fpcr |= (new_excepts << 8);
#else
    old_excepts = fenv.__control & FE_ALL_EXCEPT;

    fenv.__control &= ~new_excepts;
    fenv.__mxcsr   &= ~(new_excepts << 7);
#endif

    return fesetenv(&fenv) ? -1 : old_excepts;
}


inline int fedisableexcept(unsigned int excepts)
{
    static fenv_t fenv;
    unsigned int new_excepts = excepts & FE_ALL_EXCEPT;
    unsigned int old_excepts;   // all previous masks

    if (fegetenv(&fenv))
    {
        return -1;
    }

#if defined __arm64__
    old_excepts = fenv.__fpsr & FE_ALL_EXCEPT;

    fenv.__fpsr &= ~new_excepts;
    fenv.__fpcr &= ~(new_excepts << 8);
#else
    old_excepts = fenv.__control & FE_ALL_EXCEPT;

    fenv.__control |= new_excepts;
    fenv.__mxcsr   |= (new_excepts << 7);
#endif

    return fesetenv(&fenv) ? -1 : old_excepts;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif  // __APPLE__
#endif

// ************************************************************************* //
