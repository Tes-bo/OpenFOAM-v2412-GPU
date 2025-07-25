/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2015 OpenFOAM Foundation
    Copyright (C) 2019-2022 OpenCFD Ltd.
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

#include "ensightPTraits.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
const char* const
Foam::ensightPTraits<Foam::label>::typeName = "scalar";

template<>
const Foam::direction
Foam::ensightPTraits<Foam::label>::componentOrder[] = {0};


template<>
const char* const
Foam::ensightPTraits<float>::typeName = "scalar";

template<>
const Foam::direction
Foam::ensightPTraits<float>::componentOrder[] = {0};


template<>
const char* const
Foam::ensightPTraits<double>::typeName = "scalar";

template<>
const Foam::direction
Foam::ensightPTraits<double>::componentOrder[] = {0};


template<>
const char* const
Foam::ensightPTraits<Foam::vector>::typeName = "vector";

template<>
const Foam::direction
Foam::ensightPTraits<Foam::vector>::componentOrder[] = {0, 1, 2};


// Use mag(sphericalTensor) instead
template<>
const char* const
Foam::ensightPTraits<Foam::sphericalTensor>::typeName = "scalar";

template<>
const Foam::direction
Foam::ensightPTraits<Foam::sphericalTensor>::componentOrder[] = {0};


template<>
const char* const
Foam::ensightPTraits<Foam::symmTensor>::typeName = "tensor symm";

template<>
const Foam::direction
Foam::ensightPTraits<Foam::symmTensor>::componentOrder[] = {0, 3, 5, 1, 2, 4};


template<>
const char* const
Foam::ensightPTraits<Foam::tensor>::typeName = "tensor asym";

template<>
const Foam::direction
Foam::ensightPTraits<Foam::tensor>::componentOrder[] =
    { 0, 1, 2, 3, 4, 5, 6, 7, 8 };


// ************************************************************************* //
