/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "error.H"
#include "label.H"
#include "Istream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

#if WM_LABEL_SIZE == 32
const char* const Foam::pTraits<int32_t>::typeName = "label";
const char* const Foam::pTraits<int64_t>::typeName = "int64";
#elif WM_LABEL_SIZE == 64
const char* const Foam::pTraits<int32_t>::typeName = "int32";
const char* const Foam::pTraits<int64_t>::typeName = "label";
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::label Foam::readRawLabel(Istream& is)
{
    label val(0);
    readRawLabel(is, &val, 1);
    return val;
}


void Foam::readRawLabel(Istream& is, label* data, size_t nElem)
{
    // No check for binary vs ascii, the caller knows what they are doing

    #if WM_LABEL_SIZE == 32

    // Defined label as int32, non-native type is int64
    // Handle type narrowing limits

    typedef int64_t nonNative;

    if (is.checkLabelSize<nonNative>())
    {
        nonNative parsed;

        for (const label* endData = data + nElem; data != endData; ++data)
        {
            is.readRaw(reinterpret_cast<char*>(&parsed), sizeof(nonNative));

            // Type narrowing
            // Overflow: silently fix, or raise error?
            if (parsed < labelMin)
            {
                *data = labelMin;
            }
            else if (parsed > labelMax)
            {
                *data = labelMax;
            }
            else
            {
                *data = label(parsed);
            }
        }
    }
    else
    {
        // Read with native size
        is.readRaw(reinterpret_cast<char*>(data), nElem*sizeof(label));
    }

    #elif WM_LABEL_SIZE == 64

    // Defined label as int64, non-native type is int32

    typedef int32_t nonNative;

    if (is.checkLabelSize<nonNative>())
    {
        nonNative parsed;

        for (const label* endData = data + nElem; data != endData; ++data)
        {
            is.readRaw(reinterpret_cast<char*>(&parsed), sizeof(nonNative));

            *data = label(parsed);
        }
    }
    else
    {
        // Read with native size
        is.readRaw(reinterpret_cast<char*>(data), nElem*sizeof(label));
    }

    #endif
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::label Foam::pow(label a, label b)
{
    label ans = 1;
    for (label i=0; i<b; i++)
    {
        ans *= a;
    }

    #ifdef FULLDEBUG
    if (b < 0)
    {
        FatalErrorInFunction
            << "negative value for b is not supported"
            << abort(FatalError);
    }
    #endif

    return ans;
}


Foam::label Foam::factorial(label n)
{
    static label factTable[13] =
    {
        1, 1, 2, 6, 24, 120, 720, 5040, 40320,
        362880, 3628800, 39916800, 479001600
    };

    #ifdef FULLDEBUG
    if (n < 0 || n > 12)
    {
        FatalErrorInFunction
            << "n value out of range"
            << abort(FatalError);
    }
    #endif

    return factTable[n];
}


// ************************************************************************* //
