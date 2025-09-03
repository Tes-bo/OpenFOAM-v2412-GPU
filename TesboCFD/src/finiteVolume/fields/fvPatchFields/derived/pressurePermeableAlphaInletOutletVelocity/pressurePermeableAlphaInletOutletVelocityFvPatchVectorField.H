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
    Foam::pressurePermeableAlphaInletOutletVelocityFvPatchVectorField

Group
    grpInletBoundaryConditions grpOutletBoundaryConditions

Description
    The \c pressurePermeableAlphaInletOutletVelocity is a velocity inlet-outlet
    boundary condition which can be applied to velocity boundaries for
    multiphase flows when the pressure boundary condition is specified.

    In the \c pressurePermeableAlphaInletOutletVelocity, an open condition is
    applied when \c alpha is under a user-defined \c alphaMin value and a wall
    condition is applied when \c alpha is larger than the \c alphaMin.

    This boundary condition can be used in conjunction with
    \c prghPermeableAlphaTotalPressure for the \c p_rgh variable.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        // Mandatory entries
        type            pressurePermeableAlphaInletOutletVelocity;

        // Optional entries
        phi             phi;
        rho             rho;
        alpha           alpha.water;
        alphaMin        0.01;

        // Inherited entries
        value           uniform (0 0 0);
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description             | Type | Reqd | Deflt
      phi          | Name of flux field      | word | no   | phi
      rho          | Name of density field   | word | no   | rho
      alpha        | Name of mixture field   | word | no   | none
      alphaMin     | Minimum alpha           | scalar | no | 1
    \endtable

    The inherited entries are elaborated in:
      - \link mixedFvPatchFields.H \endlink

See also
  - Foam::prghPermeableAlphaTotalPressureFvPatchScalarField
  - Foam::mixedFvPatchVectorField

SourceFiles
    pressurePermeableAlphaInletOutletVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef pressurePermeableAlphaInletOutletVelocityFvPatchVectorField_H
#define pressurePermeableAlphaInletOutletVelocityFvPatchVectorField_H

#include "fvPatchFields.H"
#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
Class pressurePermeableAlphaInletOutletVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class pressurePermeableAlphaInletOutletVelocityFvPatchVectorField
:
    public mixedFvPatchVectorField
{
    // Private Data

        //- Name of flux field
        word phiName_;

        //- Name of density field
        word rhoName_;

        //- Name of the mixture VOF field (if used)
        word alphaName_;

        //- Minimum alpha value to outlet blockage
        scalar alphaMin_;


public:

    //- Runtime type information
    TypeName("permeableAlphaPressureInletOutletVelocity");


    // Constructors

        //- Construct from patch and internal field
        pressurePermeableAlphaInletOutletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        pressurePermeableAlphaInletOutletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //- pressurePermeableAlphaInletOutletVelocityFvPatchVectorField
        //- onto a new patch
        pressurePermeableAlphaInletOutletVelocityFvPatchVectorField
        (
            const pressurePermeableAlphaInletOutletVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        pressurePermeableAlphaInletOutletVelocityFvPatchVectorField
        (
            const pressurePermeableAlphaInletOutletVelocityFvPatchVectorField&
        );

        //- Construct as copy setting internal field reference
        pressurePermeableAlphaInletOutletVelocityFvPatchVectorField
        (
            const pressurePermeableAlphaInletOutletVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Return a clone
        virtual tmp<fvPatchField<vector>> clone() const
        {
            return fvPatchField<vector>::Clone(*this);
        }

        //- Clone with an internal field reference
        virtual tmp<fvPatchField<vector>> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return fvPatchField<vector>::Clone(*this, iF);
        }


    // Member Functions

        //- True: this patch field is altered by assignment
        virtual bool assignable() const { return true; }


        // Access

            //- Return the name of rho
            const word& rhoName() const noexcept
            {
                return rhoName_;
            }

            //- Return reference to the name of rho to allow adjustment
            word& rhoName()
            {
                return rhoName_;
            }

            //- Return the name of phi
            const word& phiName() const noexcept
            {
                return phiName_;
            }

            //- Return reference to the name of phi to allow adjustment
            word& phiName()
            {
                return phiName_;
            }


        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;


    // Member Operators

        //- Copy assignment
        virtual void operator=(const fvPatchField<vector>& pvf);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
