/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
    Copyright (C) 2019,2023 OpenCFD Ltd.
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

#include "cyclicGAMGInterfaceField.H"
#include "addToRunTimeSelectionTable.H"
#include "lduMatrix.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cyclicGAMGInterfaceField, 0);
    addToRunTimeSelectionTable
    (
        GAMGInterfaceField,
        cyclicGAMGInterfaceField,
        lduInterface
    );
    addToRunTimeSelectionTable
    (
        GAMGInterfaceField,
        cyclicGAMGInterfaceField,
        lduInterfaceField
    );
    addToRunTimeSelectionTable
    (
        GAMGInterfaceField,
        cyclicGAMGInterfaceField,
        Istream
    );

    // Add under name cyclicSlip
    addNamedToRunTimeSelectionTable
    (
        GAMGInterfaceField,
        cyclicGAMGInterfaceField,
        lduInterface,
        cyclicSlip
    );
    addNamedToRunTimeSelectionTable
    (
        GAMGInterfaceField,
        cyclicGAMGInterfaceField,
        lduInterfaceField,
        cyclicSlip
    );
    addNamedToRunTimeSelectionTable
    (
        GAMGInterfaceField,
        cyclicGAMGInterfaceField,
        Istream,
        cyclicSlip
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cyclicGAMGInterfaceField::cyclicGAMGInterfaceField
(
    const GAMGInterface& GAMGCp,
    const lduInterfaceField& fineInterface
)
:
    GAMGInterfaceField(GAMGCp, fineInterface),
    cyclicInterface_(refCast<const cyclicGAMGInterface>(GAMGCp)),
    doTransform_(false),
    rank_(0)
{
    const cyclicLduInterfaceField& p =
        refCast<const cyclicLduInterfaceField>(fineInterface);

    doTransform_ = p.doTransform();
    rank_ = p.rank();
}


Foam::cyclicGAMGInterfaceField::cyclicGAMGInterfaceField
(
    const GAMGInterface& GAMGCp,
    const bool doTransform,
    const int rank
)
:
    GAMGInterfaceField(GAMGCp, doTransform, rank),
    cyclicInterface_(refCast<const cyclicGAMGInterface>(GAMGCp)),
    doTransform_(doTransform),
    rank_(rank)
{}


Foam::cyclicGAMGInterfaceField::cyclicGAMGInterfaceField
(
    const GAMGInterface& GAMGCp,
    Istream& is
)
:
    GAMGInterfaceField(GAMGCp, is),
    cyclicInterface_(refCast<const cyclicGAMGInterface>(GAMGCp)),
    doTransform_(readBool(is)),
    rank_(readLabel(is))
{}


Foam::cyclicGAMGInterfaceField::cyclicGAMGInterfaceField
(
    const GAMGInterface& GAMGCp,
    const lduInterfaceField& local,
    const UPtrList<lduInterfaceField>& other
)
:
    GAMGInterfaceField(GAMGCp, local),
    cyclicInterface_(refCast<const cyclicGAMGInterface>(GAMGCp)),
    doTransform_(false),
    rank_(0)
{
    const auto& p = refCast<const cyclicLduInterfaceField>(local);

    doTransform_ = p.doTransform();
    rank_ = p.rank();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cyclicGAMGInterfaceField::updateInterfaceMatrix
(
    solveScalarField& result,
    const bool add,
    const lduAddressing& lduAddr,
    const label patchId,
    const solveScalarField& psiInternal,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes
) const
{
    // Get neighbouring field

    const labelList& nbrFaceCells =
        lduAddr.patchAddr
        (
            cyclicInterface_.neighbPatchID()
        );

    solveScalarField pnf(psiInternal, nbrFaceCells);

    transformCoupleField(pnf, cmpt);

    const labelList& faceCells = lduAddr.patchAddr(patchId);

    this->addToInternalField(result, !add, faceCells, coeffs, pnf);
}


void Foam::cyclicGAMGInterfaceField::write(Ostream& os) const
{
    //GAMGInterfaceField::write(os);
    os  << token::SPACE << doTransform()
        << token::SPACE << rank();
}


// ************************************************************************* //
