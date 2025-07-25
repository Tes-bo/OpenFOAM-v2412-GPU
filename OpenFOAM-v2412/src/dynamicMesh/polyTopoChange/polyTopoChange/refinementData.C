/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2019-2022,2024 OpenCFD Ltd.
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

#include "refinementData.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::refinementData Foam::transform(const tensor&, const refinementData val)
{
    return val;
}

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const refinementData& rhs
)
{
    if (os.format() == IOstreamOption::ASCII)
    {
        os << rhs.refinementCount_ << token::SPACE << rhs.count_;
    }
    else
    {
        os.write
        (
            reinterpret_cast<const char*>(&rhs.refinementCount_),
            sizeof(refinementData)
        );
    }

    os.check(FUNCTION_NAME);
    return os;
}


Foam::Istream& Foam::operator>>
(
    Istream& is,
    refinementData& rhs
)
{
    if (is.format() == IOstreamOption::ASCII)
    {
        is >> rhs.refinementCount_ >> rhs.count_;
    }
    else
    {
        Detail::readContiguous<refinementData>
        (
            is,
            reinterpret_cast<char*>(&rhs.refinementCount_),
            sizeof(refinementData)
        );
    }

    is.check(FUNCTION_NAME);
    return is;
}


// ************************************************************************* //
