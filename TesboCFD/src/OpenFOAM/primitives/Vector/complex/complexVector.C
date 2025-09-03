/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "complexVector.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
const char* const Foam::complexVector::vsType::typeName = "complexVector";

template<>
const char* const Foam::complexVector::vsType::componentNames[] =
{
    "x", "y", "z"
};

template<>
const Foam::complexVector Foam::complexVector::vsType::zero
(
    complexVector::uniform(pTraits<complex>::zero)
);

template<>
const Foam::complexVector Foam::complexVector::vsType::one
(
    complexVector::uniform(pTraits<complex>::one)
);

template<>
const Foam::complexVector Foam::complexVector::vsType::max
(
    complexVector::uniform(pTraits<complex>::max)
);

template<>
const Foam::complexVector Foam::complexVector::vsType::min
(
    complexVector::uniform(pTraits<complex>::min)
);

template<>
const Foam::complexVector Foam::complexVector::vsType::rootMax
(
    complexVector::uniform(pTraits<complex>::rootMax)
);

template<>
const Foam::complexVector Foam::complexVector::vsType::rootMin
(
    complexVector::uniform(pTraits<complex>::rootMin)
);


// ************************************************************************* //
