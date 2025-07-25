/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2023 OpenCFD Ltd.
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

#include "oversetFvPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "transform.H"
#include "cellCellStencilObject.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(oversetFvPatch, 0);
    addToRunTimeSelectionTable(fvPatch, oversetFvPatch, polyPatch);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::labelField> Foam::oversetFvPatch::interfaceInternalField
(
    const labelUList& internalData
) const
{
    // Return all internal values
    return patchInternalField(internalData);
}


Foam::tmp<Foam::labelField> Foam::oversetFvPatch::interfaceInternalField
(
    const labelUList& internalData,
    const labelUList& faceCells
) const
{
    auto tpfld = tmp<labelField>::New();
    patchInternalField(internalData, faceCells, tpfld.ref());
    return tpfld;
}


void Foam::oversetFvPatch::initInternalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const labelUList& iF
) const
{}


void Foam::oversetFvPatch::initInternalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const labelUList& iF,
    const labelUList& faceCells
) const
{}


Foam::tmp<Foam::labelField> Foam::oversetFvPatch::internalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const labelUList& iF
) const
{
    return tmp<labelField>::New(iF);
}


Foam::tmp<Foam::labelField> Foam::oversetFvPatch::internalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const labelUList& iF,
    const labelUList& nbrFaceCells
) const
{
    NotImplemented
    return tmp<labelField>::New(iF);
}


// ************************************************************************* //
