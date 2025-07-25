/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016, 2019 OpenFOAM Foundation
    Copyright (C) 2019-2022 OpenCFD Ltd.
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
    Foam::kLowReWallFunctionFvPatchScalarField

Group
    grpWallFunctions

Description
    This boundary condition provides a wall function for the turbulent kinetic
    energy (i.e. \c k) for low- and high-Reynolds number simulations.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        // Mandatory entries
        type            kLowReWallFunction;

        // Optional entries
        Ceps2           1.9;
        Ck              -0.416;
        Bk              8.366;
        C               11.0;

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                     | Type   | Reqd | Deflt
      type      | Type name: kLowReWallFunction   | word   | yes  | -
      Ceps2     | Model coefficient               | scalar | no   |  1.9
      Ck        | Model coefficient               | scalar | no   | -0.416
      Bk        | Model coefficient               | scalar | no   |  8.366
      C         | Model coefficient               | scalar | no   |  11.0
    \endtable

    The inherited entries are elaborated in:
      - \link fixedValueFvPatchField.H \endlink
      - \link wallFunctionCoefficients.H \endlink

    Viscous and inertial sublayer predictions for \c k are blended in
    a stepwise manner:

    \f[
        k = k_{log} \qquad if \quad y^+ > y^+_{intersection}
    \f]
    \f[
        k = k_{vis} \qquad if \quad y^+ <= y^+_{intersection}
    \f]
    where
    \vartable
      k_{vis}   | k prediction in the viscous sublayer
      k_{log}   | k prediction in the inertial sublayer
      y^+       | estimated wall-normal height of the cell centre in wall units
      y^+_{intersection} | estimated \f$y^+\f$ where sublayers intersect
    \endvartable

SourceFiles
    kLowReWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef kLowReWallFunctionFvPatchScalarField_H
#define kLowReWallFunctionFvPatchScalarField_H

#include "fixedValueFvPatchField.H"
#include "wallFunctionCoefficients.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class kLowReWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class kLowReWallFunctionFvPatchScalarField
:
    public fixedValueFvPatchField<scalar>
{
protected:

    // Protected Data

        //- Ceps2 coefficient
        scalar Ceps2_;

        //- Ck coefficient
        scalar Ck_;

        //- Bk coefficient
        scalar Bk_;

        //- C coefficient
        scalar C_;

        //- Wall-function coefficients
        wallFunctionCoefficients wallCoeffs_;


    // Protected Member Functions

        //- Write local wall function variables
        void writeLocalEntries(Ostream&) const;


public:

    //- Runtime type information
    TypeName("kLowReWallFunction");


    // Constructors

        //- Construct from patch and internal field
        kLowReWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        kLowReWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given kLowReWallFunctionFvPatchScalarField
        //- onto a new patch
        kLowReWallFunctionFvPatchScalarField
        (
            const kLowReWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        kLowReWallFunctionFvPatchScalarField
        (
            const kLowReWallFunctionFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        kLowReWallFunctionFvPatchScalarField
        (
            const kLowReWallFunctionFvPatchScalarField&,
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
