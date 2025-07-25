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

#include "faePatchField.H"
#include "faBoundaryMesh.H"
#include "faMesh.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(faePatchFieldBase, 0);
}

int Foam::faePatchFieldBase::disallowGenericPatchField
(
    Foam::debug::debugSwitch("disallowGenericFaePatchField", 0)
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faePatchFieldBase::faePatchFieldBase(const faPatch& p)
:
    patch_(p)
{}


Foam::faePatchFieldBase::faePatchFieldBase
(
    const faPatch& p,
    const word& patchType
)
:
    faePatchFieldBase(p)
{}


Foam::faePatchFieldBase::faePatchFieldBase
(
    const faPatch& p,
    const dictionary& dict
)
:
    faePatchFieldBase(p)
{
    faePatchFieldBase::readDict(dict);
}


Foam::faePatchFieldBase::faePatchFieldBase
(
    const faePatchFieldBase& rhs,
    const faPatch& p
)
:
    patch_(p)
{}


Foam::faePatchFieldBase::faePatchFieldBase(const faePatchFieldBase& rhs)
:
    patch_(rhs.patch_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::faePatchFieldBase::readDict(const dictionary& dict)
{}


const Foam::objectRegistry& Foam::faePatchFieldBase::db() const
{
    return patch_.boundaryMesh().mesh().thisDb();
}


void Foam::faePatchFieldBase::checkPatch(const faePatchFieldBase& rhs) const
{
    if (&patch_ != &(rhs.patch_))
    {
        FatalErrorInFunction
            << "Different patches for faePatchField"
            << abort(FatalError);
    }
}


// ************************************************************************* //
