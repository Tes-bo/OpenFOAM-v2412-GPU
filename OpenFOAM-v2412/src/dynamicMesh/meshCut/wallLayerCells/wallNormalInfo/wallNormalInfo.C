/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2020-2022 OpenCFD Ltd.
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

#include "wallNormalInfo.H"

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const wallNormalInfo& rhs
)
{
    if (os.format() == IOstreamOption::ASCII)
    {
        os << rhs.normal();
    }
    else
    {
        os.write
        (
            reinterpret_cast<const char*>(&rhs.normal_),
            sizeof(vector)
        );
    }

    os.check(FUNCTION_NAME);
    return os;
}


Foam::Istream& Foam::operator>>
(
    Istream& is,
    wallNormalInfo& rhs
)
{
    if (is.format() == IOstreamOption::ASCII)
    {
        is >> rhs.normal_;
    }
    else if (!is.checkLabelSize<>() || !is.checkScalarSize<>())
    {
        // Non-native label or scalar size
        is.beginRawRead();

        readRawScalar(is, rhs.normal_.data(), vector::nComponents);

        is.endRawRead();
    }
    else
    {
        is.read
        (
            reinterpret_cast<char*>(&rhs.normal_),
            sizeof(vector)
        );
    }

    is.check(FUNCTION_NAME);
    return is;
}


// ************************************************************************* //
