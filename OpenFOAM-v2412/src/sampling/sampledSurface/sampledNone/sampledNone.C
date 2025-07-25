/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenCFD Ltd.
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

#include "sampledNone.H"
#include "dictionary.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sampledNone, 0);
    addNamedToRunTimeSelectionTable(sampledSurface, sampledNone, word, none);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledNone::sampledNone()
:
    sampledSurface(word::null, nullptr)
{}


Foam::sampledNone::sampledNone(const word& name)
:
    sampledSurface(name, nullptr)
{}


Foam::sampledNone::sampledNone
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    sampledSurface(name, nullptr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::sampledNone::needsUpdate() const
{
    return false;
}


bool Foam::sampledNone::expire()
{
    return false;
}


bool Foam::sampledNone::update()
{
    return false;
}


#undef  makeDummy
#define makeDummy(Func,Type)                                    \
    Foam::tmp<Foam::Field<Foam::Type>>                          \
    Foam::sampledNone::Func(const interpolation<Type>&) const   \
    {                                                           \
        return tmp<Field<Type>>::New();                         \
    }

makeDummy(sample, scalar);
makeDummy(sample, vector);
makeDummy(sample, sphericalTensor);
makeDummy(sample, symmTensor);
makeDummy(sample, tensor);

makeDummy(interpolate, scalar);
makeDummy(interpolate, vector);
makeDummy(interpolate, sphericalTensor);
makeDummy(interpolate, symmTensor);
makeDummy(interpolate, tensor);

#undef makeDummy

// ************************************************************************* //
