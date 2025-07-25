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


Class
    Foam::waWallFunctionFvPatchScalarField

Description

SourceFiles
    waWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef waWallFunctionFvPatchScalarField_H
#define waWallFunctionFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"
#include "adjointBoundaryConditions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class waWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class waWallFunctionFvPatchScalarField
:
    public fixedValueFvPatchField<scalar>,
    public adjointScalarBoundaryCondition
{

    // Private member functions

        //- Check the type of the patch
        void checkType();


public:

    //- Runtime type information
    TypeName("waWallFunction");


    // Constructors

        //- Construct from patch and internal field
        waWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        waWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  waWallFunctionFvPatchScalarField
        //  onto a new patch
        waWallFunctionFvPatchScalarField
        (
            const waWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        waWallFunctionFvPatchScalarField
        (
            const waWallFunctionFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        waWallFunctionFvPatchScalarField
        (
            const waWallFunctionFvPatchScalarField&,
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


    // Member functions

        // Evaluation functions

            //- Manipulate waEqn at the cells next to the wall
            virtual void manipulateMatrix(fvMatrix<scalar>& matrix);

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //- Return the matrix diagonal coefficients corresponding to the
            //- evaluation of the value of this patchField with given weights
            virtual tmp<Field<scalar>> valueInternalCoeffs
            (
                const tmp<scalarField>&
            ) const;

            //- Return the matrix source coefficients corresponding to the
            //- evaluation of the value of this patchField with given weights
            virtual tmp<Field<scalar>> valueBoundaryCoeffs
            (
                const tmp<scalarField>&

            ) const;

            //- Return the matrix source coefficients corresponding to the
            //- evaluation of the gradient of this patchField
            virtual tmp<Field<scalar>> gradientBoundaryCoeffs() const;

            //- Return the matrix diagonal coefficients corresponding to the
            //- evaluation of the gradient of this patchField
            virtual tmp<Field<scalar>> gradientInternalCoeffs() const;



        // I-O

            //- Write
            void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
