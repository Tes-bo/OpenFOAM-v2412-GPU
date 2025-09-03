/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#include "gradingDescriptors.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::gradingDescriptors::gradingDescriptors()
:
    List<gradingDescriptor>(1, gradingDescriptor())
{}


Foam::gradingDescriptors::gradingDescriptors
(
    const label len
)
:
    List<gradingDescriptor>(len, gradingDescriptor())
{}


Foam::gradingDescriptors::gradingDescriptors(const gradingDescriptor& gd)
:
    List<gradingDescriptor>(1, gd)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::gradingDescriptors::correct()
{
    for (gradingDescriptor& gd : *this)
    {
        gd.correct();
    }
}


void Foam::gradingDescriptors::normalise()
{
    scalar sumBlockFraction = 0;
    scalar sumNDivFraction = 0;

    for (const gradingDescriptor& gd : *this)
    {
        sumBlockFraction += gd.blockFraction_;
        sumNDivFraction += gd.nDivFraction_;
    }

    for (gradingDescriptor& gd : *this)
    {
        gd.blockFraction_ /= sumBlockFraction;
        gd.nDivFraction_  /= sumNDivFraction;
        gd.correct();
    }
}


Foam::gradingDescriptors Foam::gradingDescriptors::inv() const
{
    gradingDescriptors ret(*this);

    forAll(ret, i)
    {
        ret[i] = operator[](ret.size() - i - 1).inv();
    }

    return ret;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, gradingDescriptors& gds)
{
    // Examine next token
    token t(is);

    if (t.isNumber())
    {
        gds = gradingDescriptors(gradingDescriptor(t.number()));
        gds.correct();
    }
    else
    {
        is.putBack(t);

        // Read the list for gradingDescriptors
        is >> static_cast<List<gradingDescriptor>&>(gds);

        gds.normalise();
    }

    is.check(FUNCTION_NAME);
    return is;
}


// ************************************************************************* //
