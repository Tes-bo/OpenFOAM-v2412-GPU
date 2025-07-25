/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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

#include "faPatchField.H"
#include "faBoundaryMesh.H"
#include "faMesh.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(faPatchFieldBase, 0);
}

int Foam::faPatchFieldBase::disallowGenericPatchField
(
    Foam::debug::debugSwitch("disallowGenericFaPatchField", 0)
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faPatchFieldBase::faPatchFieldBase(const faPatch& p)
:
    patch_(p),
    updated_(false),
    patchType_()
{}


Foam::faPatchFieldBase::faPatchFieldBase
(
    const faPatch& p,
    const word& patchType
)
:
    faPatchFieldBase(p)
{
    patchType_ = patchType;
}


Foam::faPatchFieldBase::faPatchFieldBase
(
    const faPatch& p,
    const dictionary& dict
)
:
    faPatchFieldBase(p)
{
    faPatchFieldBase::readDict(dict);
}


Foam::faPatchFieldBase::faPatchFieldBase
(
    const faPatchFieldBase& rhs,
    const faPatch& p
)
:
    patch_(p),
    updated_(false),
    patchType_(rhs.patchType_)
{}


Foam::faPatchFieldBase::faPatchFieldBase(const faPatchFieldBase& rhs)
:
    patch_(rhs.patch_),
    updated_(false),
    patchType_(rhs.patchType_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::faPatchFieldBase::readDict(const dictionary& dict)
{
    dict.readIfPresent("patchType", patchType_, keyType::LITERAL);
}


const Foam::objectRegistry& Foam::faPatchFieldBase::db() const
{
    return patch_.boundaryMesh().mesh().thisDb();
}


void Foam::faPatchFieldBase::checkPatch(const faPatchFieldBase& rhs) const
{
    if (&patch_ != &(rhs.patch_))
    {
        FatalErrorInFunction
            << "Different patches for faPatchField"
            << abort(FatalError);
    }
}


// ************************************************************************* //
