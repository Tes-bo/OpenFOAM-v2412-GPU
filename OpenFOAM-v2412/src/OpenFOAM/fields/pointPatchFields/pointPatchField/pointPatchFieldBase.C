/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022-2023 OpenCFD Ltd.
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

#include "pointPatchField.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pointPatchFieldBase, 0);
}

int Foam::pointPatchFieldBase::disallowGenericPatchField
(
    Foam::debug::debugSwitch("disallowGenericPointPatchField", 0)
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointPatchFieldBase::pointPatchFieldBase(const pointPatch& p)
:
    patch_(p),
    updated_(false),
    patchType_()
{}


Foam::pointPatchFieldBase::pointPatchFieldBase
(
    const pointPatch& p,
    const word& patchType
)
:
    pointPatchFieldBase(p)
{
    patchType_ = patchType;
}


Foam::pointPatchFieldBase::pointPatchFieldBase
(
    const pointPatch& p,
    const dictionary& dict
)
:
    pointPatchFieldBase(p)
{
    pointPatchFieldBase::readDict(dict);
}


Foam::pointPatchFieldBase::pointPatchFieldBase
(
    const pointPatchFieldBase& rhs,
    const pointPatch& p
)
:
    patch_(p),
    updated_(false),
    patchType_(rhs.patchType_)
{}


Foam::pointPatchFieldBase::pointPatchFieldBase(const pointPatchFieldBase& rhs)
:
    patch_(rhs.patch_),
    updated_(false),
    patchType_(rhs.patchType_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::pointPatchFieldBase::readDict(const dictionary& dict)
{
    dict.readIfPresent("patchType", patchType_, keyType::LITERAL);
}


const Foam::objectRegistry& Foam::pointPatchFieldBase::db() const
{
    return patch_.boundaryMesh().mesh().thisDb();
}


void Foam::pointPatchFieldBase::checkPatch(const pointPatchFieldBase& rhs) const
{
    if (&patch_ != &(rhs.patch_))
    {
        FatalErrorInFunction
            << "Different patches for pointPatchField"
            << abort(FatalError);
    }
}


// ************************************************************************* //
