/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2017-2019 OpenCFD Ltd.
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

#include "scalar.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::scalar Foam::readScalar(Istream& is)
{
    scalar val(0);
    is >> val;

    return val;
}


Foam::scalar Foam::readScalarOrDefault(Istream& is, const scalar defaultValue)
{
    if (is.good())
    {
        token tok(is);

        if (tok.isNumber())
        {
            return tok.scalarToken();
        }

        is.putBack(tok);
    }

    return defaultValue;
}


Foam::scalar Foam::readRawScalar(Istream& is)
{
    scalar val(0);
    readRawScalar(is, &val, 1);
    return val;
}


void Foam::readRawScalar(Istream& is, scalar* data, size_t nElem)
{
    // No check for binary vs ascii, the caller knows what they are doing

    #if defined(WM_SP) || defined(WM_SPDP)

    // Defined scalar as a float, non-native type is double
    // Handle type narrowing limits

    typedef double nonNative;

    if (is.checkScalarSize<nonNative>())
    {
        nonNative other;

        for (const scalar* endData = data + nElem; data != endData; ++data)
        {
            is.readRaw(reinterpret_cast<char*>(&other), sizeof(nonNative));

            // Type narrowing
            // Overflow: silently fix, or raise error?

            if (other < -VGREAT)
            {
                *data = -VGREAT;
            }
            else if (other > VGREAT)
            {
                *data = VGREAT;
            }
            else if (other > -VSMALL && other < VSMALL)
            {
                // Underflow: round to zero
                *data = 0;
            }
            else
            {
                *data = scalar(other);
            }
        }
    }
    else
    {
        // Read with native size
        is.readRaw(reinterpret_cast<char*>(data), nElem*sizeof(scalar));
    }

    #elif defined(WM_DP)

    // Defined scalar as a double, non-native type is float

    typedef float nonNative;

    if (is.checkScalarSize<nonNative>())
    {
        nonNative other;

        for (const scalar* endData = data + nElem; data != endData; ++data)
        {
            is.readRaw(reinterpret_cast<char*>(&other), sizeof(nonNative));

            *data = scalar(other);
        }
    }
    else
    {
        // Read with native size
        is.readRaw(reinterpret_cast<char*>(data), nElem*sizeof(scalar));
    }

    #endif
}


// ************************************************************************* //
