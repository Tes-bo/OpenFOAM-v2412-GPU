/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2012 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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

#include "mixedEnergyFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "basicThermo.H"

#include "fvMatrix.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mixedEnergyFvPatchScalarField::
mixedEnergyFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF)
{
    valueFraction() = 0.0;
    refValue() = Zero;
    refGrad() = Zero;
    source() = 0.0;
}


Foam::mixedEnergyFvPatchScalarField::
mixedEnergyFvPatchScalarField
(
    const mixedEnergyFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::mixedEnergyFvPatchScalarField::
mixedEnergyFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF, dict)
{}


Foam::mixedEnergyFvPatchScalarField::
mixedEnergyFvPatchScalarField
(
    const mixedEnergyFvPatchScalarField& tppsf
)
:
    mixedFvPatchScalarField(tppsf)
{}


Foam::mixedEnergyFvPatchScalarField::
mixedEnergyFvPatchScalarField
(
    const mixedEnergyFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(tppsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mixedEnergyFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    const basicThermo& thermo = basicThermo::lookupThermo(*this);
    const label patchi = patch().index();

    const scalarField& pw = thermo.p().boundaryField()[patchi];
    mixedFvPatchScalarField& Tw = refCast<mixedFvPatchScalarField>
    (
        const_cast<fvPatchScalarField&>(thermo.T().boundaryField()[patchi])
    );

    Tw.evaluate();

    valueFraction() = Tw.valueFraction();
    refValue() = thermo.he(pw, Tw.refValue(), patchi);
    refGrad() =
        thermo.Cpv(pw, Tw, patchi)*Tw.refGrad()
      + patch().deltaCoeffs()*
        (
            thermo.he(pw, Tw, patchi)
          - thermo.he(pw, Tw, patch().faceCells())
        );

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::mixedEnergyFvPatchScalarField::manipulateMatrix
(
    fvMatrix<scalar>& matrix,
    const label mat,
    const direction cmpt
)
{
    const basicThermo& thermo = basicThermo::lookupThermo(*this);

    label index = this->patch().index();

    const label nbrPatchId =  this->patch().patch().neighbPolyPatchID();

    const label globalPatchID =
        matrix.lduMeshAssembly().patchLocalToGlobalMap()[mat][index];

    const label meshNrbId = matrix.lduMeshAssembly().findNbrMeshId
    (
        this->patch().patch(),
        mat
    );

    const mixedFvPatchField<scalar>& fPatch =
        refCast<const mixedFvPatchField>(thermo.T().boundaryField()[index]);

    const Field<scalar> intCoeffsCmpt
    (
        matrix.internalCoeffs()[globalPatchID].component(cmpt)
    );

    const scalarField sourceCorr(fPatch.source());

    const labelList& faceMap =
        matrix.lduMeshAssembly().faceBoundMap()[mat][index];

    const labelList& myCells =
        matrix.lduMeshAssembly().cellBoundMap()[meshNrbId][nbrPatchId];

    const labelList& nbrCells =
        matrix.lduMeshAssembly().cellBoundMap()[mat][index];

    forAll(faceMap, j)
    {
        label globalFaceI = faceMap[j];

        label myCellI = myCells[j];
        label nbrCellI = nbrCells[j];

        const scalar intCorr = -intCoeffsCmpt[j];
        const scalar srcCorr = -sourceCorr[j];

        if (this->patch().patch().masterImplicit())
        {
            if (myCellI > nbrCellI)
            {
                if (matrix.asymmetric())
                {
                    matrix.lower()[globalFaceI] += intCorr;
                }
            }
            else
            {
                matrix.upper()[globalFaceI] += intCorr;
            }

            matrix.diag()[myCellI] -= intCorr;
            matrix.source()[myCellI] += srcCorr;
        }
        else
        {
            if (myCellI < nbrCellI)
            {
                matrix.upper()[globalFaceI] += intCorr;
            }
            else
            {
                if (matrix.asymmetric())
                {
                    matrix.lower()[globalFaceI] += intCorr;
                }
            }

            matrix.diag()[myCellI] -= intCorr;
            matrix.source()[myCellI] += srcCorr;
        }


//         if (globalFaceI != -1)
//         {
//             const scalar intCorr = -intCoeffsCmpt[j];
//             const scalar srcCorr = -sourceCorr[j];
//
//             if (this->patch().patch().masterImplicit())
//             {
//                 matrix.diag()[u[globalFaceI]] -= intCorr;
//                 if (matrix.asymmetric())
//                 {
//                     matrix.lower()[globalFaceI] += intCorr;
//                 }
//                 matrix.source()[u[globalFaceI]] += srcCorr;
//             }
//             else
//             {
//                 matrix.diag()[l[globalFaceI]] -= intCorr;
//                 matrix.upper()[globalFaceI] += intCorr;
//                 matrix.source()[l[globalFaceI]] += srcCorr;
//             }
//         }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        mixedEnergyFvPatchScalarField
    );
}

// ************************************************************************* //
