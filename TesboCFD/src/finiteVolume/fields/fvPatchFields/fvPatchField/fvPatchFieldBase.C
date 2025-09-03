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

#include "fvPatchField.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fvPatchFieldBase, 0);
}

int Foam::fvPatchFieldBase::disallowGenericPatchField
(
    Foam::debug::debugSwitch("disallowGenericFvPatchField", 0)
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvPatchFieldBase::fvPatchFieldBase(const fvPatch& p)
:
    patch_(p),
    updated_(false),
    manipulatedMatrix_(false),
    useImplicit_(false),
    patchType_()
{}


Foam::fvPatchFieldBase::fvPatchFieldBase
(
    const fvPatch& p,
    const word& patchType
)
:
    fvPatchFieldBase(p)
{
    patchType_ = patchType;
}


Foam::fvPatchFieldBase::fvPatchFieldBase
(
    const fvPatch& p,
    const dictionary& dict
)
:
    fvPatchFieldBase(p)
{
    fvPatchFieldBase::readDict(dict);
}


Foam::fvPatchFieldBase::fvPatchFieldBase
(
    const fvPatchFieldBase& rhs,
    const fvPatch& p
)
:
    patch_(p),
    updated_(false),
    manipulatedMatrix_(false),
    useImplicit_(rhs.useImplicit_),
    patchType_(rhs.patchType_)
{}


Foam::fvPatchFieldBase::fvPatchFieldBase(const fvPatchFieldBase& rhs)
:
    patch_(rhs.patch_),
    updated_(false),
    manipulatedMatrix_(false),
    useImplicit_(rhs.useImplicit_),
    patchType_(rhs.patchType_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fvPatchFieldBase::readDict(const dictionary& dict)
{
    dict.readIfPresent("patchType", patchType_, keyType::LITERAL);
    dict.readIfPresent("useImplicit", useImplicit_, keyType::LITERAL);
}


const Foam::objectRegistry& Foam::fvPatchFieldBase::db() const
{
    return patch_.boundaryMesh().mesh().thisDb();
}


void Foam::fvPatchFieldBase::checkPatch(const fvPatchFieldBase& rhs) const
{
    if (&patch_ != &(rhs.patch_))
    {
        FatalErrorInFunction
            << "Different patches for fvPatchField"
            << abort(FatalError);
    }
}


// ************************************************************************* //
