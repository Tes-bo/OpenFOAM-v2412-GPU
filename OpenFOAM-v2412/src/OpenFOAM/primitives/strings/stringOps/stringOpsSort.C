/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenCFD Ltd.
-------------------------------------------------------------------------------
Changes for OpenFOAM

   Code cleanup, reduction, elimination of redundant code.
\*---------------------------------------------------------------------------*/

//========================================================================
// Copyright (c) 1998-2010,2011 Free Software Foundation, Inc.
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights to use, copy, modify, merge, publish,
// distribute, distribute with modifications, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
// IN NO EVENT SHALL THE ABOVE COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
// DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
// OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
// THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//
// Except as contained in this notice, the name(s) of the above copyright
// holders shall not be used in advertising or otherwise to promote the
// sale, use or other dealings in this Software without prior written
// authorization.
//========================================================================

//========================================================================
//  Original Author: Jan-Marten Spit <jmspit@euronet.nl>
//========================================================================

#include "stringOpsSort.H"
#include <cctype>
#include <cstring>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Local tweaks/preferences:
//
// [DIGITS_ALWAYS_FIRST] : as per original code
// This results in "file123.txt" sorting before "file.txt", which is
// inconsistent with what 'ls -v' produces
// - normally do not want this (Mark Olesen: Oct-2017)
#undef DIGITS_ALWAYS_FIRST

// [MANUAL_NUMCOMPARE] : handwritten code instead of strncmp
// The original code has a mix of strncmp for equality but handwritten code
// for greater-than/less-than.
//
// -> this does to be unneeded, rely on strncmp() return values
//    (Mark Olesen: Oct-2017)
#undef MANUAL_NUMCOMPARE

// [DEBUG_NATSTRCMP] : debug info to std::cerr (for development purposes only)
#undef DEBUG_NATSTRCMP

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef DEBUG_NATSTRCMP
#include <iostream>

template<class T>
static inline void debugPrint(const char* text, T item1, T item2)
{
    std::cerr << text << ": <" << item1 << "> <" << item2 << ">\n";
}

// Character sequences, end pointer is _inclusive_.
static inline void debugPrint
(
    const char* b1, const char* e1,
    const char* b2, const char* e2
)
{
    std::cerr << "<";
    while (b1 < e1) { std::cerr << *b1++; }
    std::cerr << *e1 << "> ";

    std::cerr << "<";
    while (b2 < e2) { std::cerr << *b2++; }
    std::cerr << *e2 << ">\n";
}
#endif

#ifdef MANUAL_NUMCOMPARE
// Manual comparison of identical length (digit) sequences
static inline int manual_numcompare(const char* s1, const char* s2)
{
    while (*s1 && *s2)
    {
        const int cmp = (*s1 - *s2);
        if (!cmp) return cmp;

        ++s1; ++s2;
    }
    return 0;  // Length check done before in caller.
}
#endif


// ------------------------------------------------------------------------- //

int Foam::stringOps::natstrcmp(const char* s1, const char* s2)
{
    #ifdef DEBUG_NATSTRCMP
    debugPrint("natstrcmp", s1, s2);
    #endif

    // States for state engine
    enum stateType { SCAN, ALPHA, NUMERIC };

    // Number of leading zeroes
    unsigned zeros1 = 0;
    unsigned zeros2 = 0;

    // Pointers to begin/end of integer sequences (without leading zeros)
    const char* numbeg1 = nullptr;
    const char* numbeg2 = nullptr;
    const char* numend1 = nullptr;
    const char* numend2 = nullptr;

    stateType state = SCAN;

    const char* p1 = s1;
    const char* p2 = s2;

    while (*p1 && *p2)
    {
        // Bitmask for digits vs alpha
        // - 0: neither are digits
        // - 1: p1 is the only digit
        // - 2: p2 is the only digit
        // - 3: both p1 and p2 are digits

        const unsigned digitMask =
            ((isdigit(*p2) ? 2:0) | (isdigit(*p1) ? 1:0));

        switch (state)
        {
            case SCAN:
            {
                #ifdef DEBUG_NATSTRCMP
                debugPrint("SCAN", *p1, *p2);
                #endif

                switch (digitMask)
                {
                    case 0:         // (alpha,alpha)
                    {
                        state = ALPHA;

                        if (*p1 == *p2)
                        {
                            ++p1; ++p2;
                        }
                        else
                        {
                            // Lexical compare
                            return (*p1 - *p2);
                        }
                        break;
                    }

                    #ifdef DIGITS_ALWAYS_FIRST
                    case 0x1:        // (digit,alpha) : digit < alpha
                    {
                        return -1;
                        break;
                    }
                    case 0x2:        // (alpha,digit) : alpha > digit
                    {
                        return 1;
                        break;
                    }
                    #else /* DIGITS_ALWAYS_FIRST */
                    case 0x1:        // (digit,alpha)
                    case 0x2:        // (alpha,digit)
                    {
                        // Lexical compare for digits/alpha
                        return (*p1 - *p2);
                        break;
                    }
                    #endif /* DIGITS_ALWAYS_FIRST */

                    default:        // (digit,digit)
                    {
                        state = NUMERIC;

                        // State changed from SCAN to NUMERIC, so skip leading
                        // leading zeros so the numeric comparison is
                        // untainted by them, but capture the first occurrence
                        // of leading zeroes for a final tie-break if needed.

                        if (!zeros1)
                        {
                            // First occurrence of any leading zeroes
                            while (*p1 == '0') { ++p1; ++zeros1; }
                        }
                        else
                        {
                            while (*p1 == '0') { ++p1; }
                        }

                        if (!zeros2)
                        {
                            // First occurrence of any leading zeroes
                            while (*p2 == '0') { ++p2; ++zeros2; }
                        }
                        else
                        {
                            while (*p2 == '0') { ++p2; }
                        }

                        if (zeros1 == zeros2)
                        {
                            // Same number of zeros - so irrelevant
                            zeros1 = zeros2 = 0;
                        }

                        if (!isdigit(*p1)) --p1;
                        if (!isdigit(*p2)) --p2;

                        numbeg1 = numend1 = p1;
                        numbeg2 = numend2 = p2;

                        break;
                    }
                }
                break;
            }

            case ALPHA:
            {
                #ifdef DEBUG_NATSTRCMP
                debugPrint("ALPHA", *p1, *p2);
                #endif

                if (digitMask)
                {
                    state = SCAN;
                }
                else                // (alpha,alpha)
                {
                    if (*p1 == *p2)
                    {
                        ++p1; ++p2;
                    }
                    else
                    {
                        // Lexical compare
                        return (*p1 - *p2);
                    }
                }
                break;
            }

            case NUMERIC:
            {
                while (isdigit(*p1)) numend1 = p1++;
                while (isdigit(*p2)) numend2 = p2++;

                #ifdef DEBUG_NATSTRCMP
                debugPrint("NUMERIC", *p1, *p2);
                debugPrint(numbeg1,numend1, numbeg2,numend2);
                #endif

                // Length (minus 1) of each sequence
                const size_t len1 = (numend1 - numbeg1);
                const size_t len2 = (numend2 - numbeg2);

                if (len1 < len2)
                {
                    return -1;
                }
                else if (len1 > len2)
                {
                    return 1;
                }

                // Same number of digits, leading zeros have been skipped.
                // - so lexical and numerical compares are equivalent.

                const int cmp = strncmp(numbeg1, numbeg2, len1+1);
                if (!cmp)
                {
                    // Identical (digit) sequence - continue
                    state = SCAN;
                }
                else
                {
                    #ifdef MANUAL_STRNCMP
                    return manual_numcompare(numbeg1, numbeg2);
                    #else
                    return cmp;
                    #endif
                }
                break;
            }
        }
    }

    if (zeros1 < zeros2) return -1;
    if (zeros1 > zeros2) return 1;
    if (!*p1 && *p2) return -1;  // s1 shorter than s2
    if (*p1 && !*p2) return 1;   // s1 longer  than s2

    return 0;
}


// ************************************************************************* //
