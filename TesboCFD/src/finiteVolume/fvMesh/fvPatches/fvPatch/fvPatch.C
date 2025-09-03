/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020-2022 OpenCFD Ltd.
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

#include "fvPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "fvBoundaryMesh.H"
#include "fvMesh.H"
#include "primitiveMesh.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fvPatch, 0);
    defineRunTimeSelectionTable(fvPatch, polyPatch);
    addToRunTimeSelectionTable(fvPatch, fvPatch, polyPatch);
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

const Foam::fvPatch& Foam::fvPatch::lookupPatch(const polyPatch& p)
{
    const fvMesh* meshptr = isA<fvMesh>(p.boundaryMesh().mesh());

    if (!meshptr)
    {
        FatalErrorInFunction
            << "The polyPatch is not attached to a base fvMesh" << nl
            << exit(FatalError);
    }

    return meshptr->boundary()[p.index()];
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvPatch::fvPatch(const polyPatch& p, const fvBoundaryMesh& bm)
:
    polyPatch_(p),
    boundaryMesh_(bm)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvPatch::~fvPatch()
{}  // fvBoundaryMesh was forward declared


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fvPatch::constraintType(const word& patchType)
{
    return
    (
        !patchType.empty()
     && fvPatchField<scalar>::patchConstructorTablePtr_
     && fvPatchField<scalar>::patchConstructorTablePtr_->found(patchType)
    );
}


Foam::wordList Foam::fvPatch::constraintTypes()
{
    const auto& cnstrTable = *polyPatchConstructorTablePtr_;

    wordList cTypes(cnstrTable.size());

    label i = 0;

    forAllConstIters(cnstrTable, iter)
    {
        if (constraintType(iter.key()))
        {
            cTypes[i++] = iter.key();
        }
    }

    cTypes.resize(i);

    return cTypes;
}


const Foam::labelUList& Foam::fvPatch::faceCells() const
{
    return polyPatch_.faceCells();
}


const Foam::vectorField& Foam::fvPatch::Cf() const
{
    return boundaryMesh().mesh().Cf().boundaryField()[index()];
}


Foam::tmp<Foam::vectorField> Foam::fvPatch::Cn() const
{
    return patchInternalField(boundaryMesh().mesh().cellCentres());
}


const Foam::vectorField& Foam::fvPatch::Sf() const
{
    return boundaryMesh().mesh().Sf().boundaryField()[index()];
}


const Foam::scalarField& Foam::fvPatch::magSf() const
{
    return boundaryMesh().mesh().magSf().boundaryField()[index()];
}


Foam::tmp<Foam::vectorField> Foam::fvPatch::unitSf() const
{
    return Sf()/magSf();
}


Foam::tmp<Foam::vectorField> Foam::fvPatch::nf() const
{
    return Sf()/magSf();
}


Foam::tmp<Foam::vectorField> Foam::fvPatch::delta() const
{
    // Use patch-normal delta for all non-coupled BCs
    const vectorField nHat(nf());
    return nHat*(nHat & (Cf() - Cn()));
}


void Foam::fvPatch::makeWeights(scalarField& w) const
{
    w = 1.0;
}


void Foam::fvPatch::makeDeltaCoeffs(scalarField& w) const
{}


void Foam::fvPatch::makeNonOrthoDeltaCoeffs(scalarField& w) const
{}


void Foam::fvPatch::makeNonOrthoCorrVectors(vectorField& w) const
{}


void Foam::fvPatch::initMovePoints()
{}


void Foam::fvPatch::movePoints()
{}


const Foam::scalarField& Foam::fvPatch::deltaCoeffs() const
{
    return boundaryMesh().mesh().deltaCoeffs().boundaryField()[index()];
}


const Foam::scalarField& Foam::fvPatch::weights() const
{
    return boundaryMesh().mesh().weights().boundaryField()[index()];
}


// ************************************************************************* //
