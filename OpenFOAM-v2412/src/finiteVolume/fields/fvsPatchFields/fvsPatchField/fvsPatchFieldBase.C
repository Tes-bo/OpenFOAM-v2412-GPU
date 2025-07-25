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

#include "fvsPatchField.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fvsPatchFieldBase, 0);
}

int Foam::fvsPatchFieldBase::disallowGenericPatchField
(
    Foam::debug::debugSwitch("disallowGenericFvsPatchField", 0)
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvsPatchFieldBase::fvsPatchFieldBase(const fvPatch& p)
:
    patch_(p)
{}


Foam::fvsPatchFieldBase::fvsPatchFieldBase
(
    const fvPatch& p,
    const dictionary& dict
)
:
    patch_(p)
{}


Foam::fvsPatchFieldBase::fvsPatchFieldBase
(
    const fvsPatchFieldBase& rhs,
    const fvPatch& p
)
:
    patch_(p)
{}


Foam::fvsPatchFieldBase::fvsPatchFieldBase(const fvsPatchFieldBase& rhs)
:
    patch_(rhs.patch_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fvsPatchFieldBase::readDict(const dictionary& dict)
{}


const Foam::objectRegistry& Foam::fvsPatchFieldBase::db() const
{
    return patch_.boundaryMesh().mesh().thisDb();
}


void Foam::fvsPatchFieldBase::checkPatch(const fvsPatchFieldBase& rhs) const
{
    if (&patch_ != &(rhs.patch_))
    {
        FatalErrorInFunction
            << "Different patches for fvsPatchField"
            << abort(FatalError);
    }
}


// ************************************************************************* //
