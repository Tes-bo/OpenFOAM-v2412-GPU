/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2022 PCOpt/NTUA
    Copyright (C) 2014-2022 FOSS GP
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

#include "waWallFunctionFvPatchScalarField.H"
#include "wallFvPatch.H"
#include "omegaWallFunctionFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void waWallFunctionFvPatchScalarField::checkType()
{
    if (!isA<wallFvPatch>(patch()))
    {
        FatalErrorInFunction
            << "Invalid wall function specification" << nl
            << "    Patch type for patch " << patch().name()
            << " must be wall" << nl
            << "    Current patch type is " << patch().type() << nl << endl
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

waWallFunctionFvPatchScalarField::waWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF),
    adjointScalarBoundaryCondition(p, iF, "wa")
{
    checkType();
}


waWallFunctionFvPatchScalarField::waWallFunctionFvPatchScalarField
(
    const waWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper),
    adjointScalarBoundaryCondition(p, iF, ptf.adjointSolverName_)
{
    checkType();
}


waWallFunctionFvPatchScalarField::waWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p, iF, dict),
    adjointScalarBoundaryCondition(p, iF, dict.get<word>("solverName"))
{
    checkType();
}


waWallFunctionFvPatchScalarField::waWallFunctionFvPatchScalarField
(
    const waWallFunctionFvPatchScalarField& ewfpsf
)
:
    fixedValueFvPatchField<scalar>(ewfpsf),
    adjointScalarBoundaryCondition(ewfpsf)
{
    checkType();
}


waWallFunctionFvPatchScalarField::waWallFunctionFvPatchScalarField
(
    const waWallFunctionFvPatchScalarField& ewfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(ewfpsf, iF),
    adjointScalarBoundaryCondition(ewfpsf)
{
    checkType();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void waWallFunctionFvPatchScalarField::manipulateMatrix
(
    fvMatrix<scalar>& matrix
)
{
    scalarField& Diag = matrix.diag();
    scalarField& lower = matrix.lower();
    scalarField& upper = matrix.upper();
    FieldField<Field, scalar>& internalCoeffs = matrix.internalCoeffs();
    FieldField<Field, scalar>& boundaryCoeffs = matrix.boundaryCoeffs();
    const fvMesh& mesh = patch().boundaryMesh().mesh();
    const labelList& faceCells = patch().faceCells();

    // Add diag term from the omega expression next to the wall
    for (const label celli : faceCells)
    {
        Diag[celli] = 1;
    }

    // We want something similar to setValues, but slightly modified.
    // The solution of the boundary cell should understand contributions from
    // the second cells off the wall but they should see the
    // solution of the boundary cell as zero.
    // Contributions from neighbouring cells with an omegaWallFunction boundary
    // condition should also be zero

    const cellList& cells = mesh.cells();
    const labelUList& own = mesh.owner();

    /*
    const labelUList& nei = mesh.neighbour();
    const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );

    const tmp<volScalarField> tomega = turbModel.omega();
    const volScalarField& omega = tomega();
    typedef omegaWallFunctionFvPatchScalarField omegaWF;
    */

    forAll(faceCells, i)
    {
        const label celli = faceCells[i];
        const cell& c = cells[celli];

        forAll(c, j)
        {
            const label facei = c[j];

            if (mesh.isInternalFace(facei))
            {
                // Neighbouring cells should get no contribution from
                // ourselves in all cases
                //label cellNei(-1);
                if (celli == own[facei])
                {
                    //cellNei = nei[facei];
                    lower[facei] = 0.0;
                }
                else
                {
                    //cellNei = own[facei];
                    upper[facei] = 0.0;
                }
                // Additionally, if the neighbouring cell is also a boundary
                // one with omegaWallFunction in one of its faces,
                // contributions between the two cells should be canceled out
                // as well.
                // Already covered by the above
                /*
                bool neiHasOmegaWFface(false);
                const cell& neiCell = cells[cellNei];
                forAll(neiCell, fNei)
                {
                    const label faceNei = neiCell[fNei];

                    const label patchNei =
                        mesh.boundaryMesh().whichPatch(faceNei);
                    if (patchNei != -1)
                    {
                        const fvPatchField& omegaNei =
                            omega.boundaryField()[patchNei];
                        if (isA<omegaWF>(omegaNei))
                        {
                            neiHasOmegaWFface = true;
                            break;
                        }
                    }
                }
                if (neiHasOmegaWFface)
                {
                    if (celli == own[facei])
                    {
                        upper[facei] = 0.0;
                    }
                    else
                    {
                        lower[facei] = 0.0;
                    }
                }
                */
            }
            // Contributions from boundaries should have already been removed
            // using value*Coeffs and boundary*Coeffs
            // Just to be safe
            else
            {
                const label patchi = mesh.boundaryMesh().whichPatch(facei);
                if (internalCoeffs[patchi].size())
                {
                    label patchFacei =
                        mesh.boundaryMesh()[patchi].whichFace(facei);
                    internalCoeffs[patchi][patchFacei] = Zero;
                    boundaryCoeffs[patchi][patchFacei] = Zero;
                }
            }
        }
    }

    fvPatchField<scalar>::manipulateMatrix(matrix);
}


void waWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    operator == (scalarField(patch().size(), Zero));
    fixedValueFvPatchField<scalar>::updateCoeffs();
}


tmp<Field<scalar>> waWallFunctionFvPatchScalarField::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    return tmp<Field<scalar>>::New(this->size(), Zero);
}


tmp<Field<scalar>> waWallFunctionFvPatchScalarField::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    return tmp<Field<scalar>>::New(this->size(), Zero);
}


tmp<Field<scalar>>
waWallFunctionFvPatchScalarField::gradientBoundaryCoeffs() const
{
    return tmp<Field<scalar>>::New(this->size(), Zero);
}


tmp<Field<scalar>>
waWallFunctionFvPatchScalarField::gradientInternalCoeffs() const
{
    return tmp<Field<scalar>>::New(this->size(), Zero);
}


void waWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fixedValueFvPatchField<scalar>::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    waWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
