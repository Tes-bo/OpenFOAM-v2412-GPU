/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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

Class
    Foam::fixedShearStressFvPatchVectorField

Group
    grpTurbulenceBoundaryConditions

Description
    This boundary condition sets a user-defined shear stress
    constant and uniform across a given patch by using the expression:

    \f[
        tau_0 = -\nu_{eff} \frac{dU}{dn}
    \f]

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        // Mandatory entries
        type             fixedShearStress;
        tau              <vector>;

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                      | Type | Reqd | Deflt
      type      | Type name: fixedShearStress      | word | yes  | -
      tau       | Shear stress                   | vector | yes  | -
    \endtable

    The inherited entries are elaborated in:
      - \link fixedValueFvPatchFields.H \endlink

SourceFiles
    fixedShearStressFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef fixedShearStressFvPatchVectorField_H
#define fixedShearStressFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class fixedShearStressFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class fixedShearStressFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private Data

        //- Constant shear stress
        const vector tau0_;


public:

    //- Runtime type information
    TypeName("fixedShearStress");


    // Constructors

        //- Construct from patch and internal field
        fixedShearStressFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixedShearStressFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        fixedShearStressFvPatchVectorField
        (
            const fixedShearStressFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fixedShearStressFvPatchVectorField
        (
            const fixedShearStressFvPatchVectorField&
        );

        //- Construct as copy setting internal field reference
        fixedShearStressFvPatchVectorField
        (
            const fixedShearStressFvPatchVectorField&,
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

        // Evaluation

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
