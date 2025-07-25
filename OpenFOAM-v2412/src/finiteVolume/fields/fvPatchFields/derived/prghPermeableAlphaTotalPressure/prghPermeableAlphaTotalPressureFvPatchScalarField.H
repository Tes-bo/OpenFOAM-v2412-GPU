/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

Class
    Foam::prghPermeableAlphaTotalPressureFvPatchScalarField

Description
    The \c prghPermeableAlphaTotalPressure is a mixed boundary condition
    for the \c p_rgh variable in multiphase flows.

    It switches between an open boundary using the \c prghTotalPressure
    boundary condition and a wall condition using the \c fixedFluxPressure
    boundary condition.

    In the \c prghPermeableAlphaTotalPressure, an open condition is
    applied when \c alpha is under a user-defined \c alphaMin value and a wall
    condition is applied when \c alpha is larger than the \c alphaMin.

    This boundary condition can be used in conjunction with
    \c pressurePermeableAlphaInletOutletVelocity for the \c U variable.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        // Mandatory entries
        type            prghPermeableAlphaTotalPressure;
        p               uniform 0;

        // Optional entries
        phi             phi;
        rho             rho;
        U               U;
        alphaName       alpha.water;
        alphaMin        0.01;

        // Inherited entries
        value           uniform 0;
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description             | Type | Reqd | Deflt
      p         | Total pressure          | PatchFunction1\<scalar\> | yes | -
      phi       | Name of flux field      | word | no   | phi
      rho       | Name of density field   | word | no   | rho
      U         | Name of velocity field  | word | no   | U
      alpha     | Name of mixture field   | word | no   | none
      alphaMin  | Minimum alpha           | scalar | no | 1
    \endtable

    The inherited entries are elaborated in:
      - \link mixedFvPatchFields.H \endlink
      - \link PatchFunction1.H \endlink

See also
  - Foam::pressurePermeableAlphaInletOutletVelocityFvPatchVectorField
  - Foam::mixedFvPatchVectorField
  - Foam::prghTotalPressureFvPatchScalarField

SourceFiles
    prghPermeableAlphaTotalPressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef prghPermeableAlphaTotalPressureFvPatchScalarField_H
#define prghPermeableAlphaTotalPressureFvPatchScalarField_H

#include "mixedFvPatchField.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "updateableSnGrad.H"
#include "PatchFunction1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
    Class prghPermeableAlphaTotalPressureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class prghPermeableAlphaTotalPressureFvPatchScalarField
:
    public mixedFvPatchField<scalar>,
    public updateablePatchTypes::updateableSnGrad
{
    // Private Data

        //- Total pressure field
        autoPtr<PatchFunction1<scalar>> p0_;

        //- Name of the flux transporting the field
        word phiName_;

        //- Name of the density field used to normalise the mass flux
        word rhoName_;

         //- Name of the velocity field
        word UName_;

        //- Name of the mixture VOF field (if used)
        word alphaName_;

        //- Minimum alpha value to outlet blockage
        scalar alphaMin_;

        //- Current time index (used for updating)
        label curTimeIndex_;


public:

    //- Runtime type information
    TypeName("prghPermeableAlphaTotalPressure");


    // Constructors

        //- Construct from patch and internal field
        prghPermeableAlphaTotalPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        prghPermeableAlphaTotalPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //- prghPermeableAlphaTotalPressureFvPatchScalarField
        //- onto a new patch
        prghPermeableAlphaTotalPressureFvPatchScalarField
        (
            const prghPermeableAlphaTotalPressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        prghPermeableAlphaTotalPressureFvPatchScalarField
        (
            const prghPermeableAlphaTotalPressureFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        prghPermeableAlphaTotalPressureFvPatchScalarField
        (
            const prghPermeableAlphaTotalPressureFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Return a clone
        virtual tmp<fvPatchField<scalar>> clone() const
        {
            return fvPatchField<scalar>::Clone(*this);
        }

        //- Clone with an internal field reference
        virtual tmp<fvPatchField<scalar>> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return fvPatchField<scalar>::Clone(*this, iF);
        }


    // Member Functions

        // Mapping

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchScalarField&,
                const labelList&
            );


        // Evaluation

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //- Update the patch pressure gradient field from the given snGradp
            virtual void updateSnGrad(const scalarField& snGradp);


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
