/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2020 PCOpt/NTUA
    Copyright (C) 2013-2020 FOSS GP
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::adjointFarFieldPressureFvPatchScalarField

Description

SourceFiles
    adjointFarFieldPressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef adjointFarFieldPressureFvPatchScalarField_H
#define adjointFarFieldPressureFvPatchScalarField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "adjointBoundaryConditions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class adjointOutletPressureFvPatch Declaration
\*---------------------------------------------------------------------------*/

class adjointFarFieldPressureFvPatchScalarField
:
    public fixedValueFvPatchScalarField,
    public adjointScalarBoundaryCondition
{
public:

    //- Runtime type information
    TypeName("adjointFarFieldPressure");


    // Constructors

        //- Construct from patch and internal field
        adjointFarFieldPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        adjointFarFieldPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given adjointOutletPressureFvPatchScalarField
        // onto a new patch
        adjointFarFieldPressureFvPatchScalarField
        (
            const adjointFarFieldPressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        adjointFarFieldPressureFvPatchScalarField
        (
            const adjointFarFieldPressureFvPatchScalarField&,
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

        //- Return true if this patch field fixes a value.
        //  Needed to check if a level has to be specified while solving
        //  Poissons equations.
        /*
        virtual bool fixesValue() const
        {
            return true;
        }
        */
        //- Return gradient at boundary
        virtual tmp<Field<scalar>> snGrad() const;

        //- Return the matrix diagonal coefficients corresponding to the
        //- evaluation of the value of this patchField with given weights
        virtual tmp<Field<scalar>> valueInternalCoeffs
        (
            const tmp<scalarField>&
        ) const;

        //- Return the matrix source coefficients corresponding to the
        //  evaluation of the value of this patchField with given weights
        virtual tmp<Field<scalar>> valueBoundaryCoeffs
        (
            const tmp<scalarField>&
        ) const;

        //- Return the matrix diagonal coefficients corresponding to the
        //- evaluation of the gradient of this patchField
        virtual tmp<Field<scalar>> gradientInternalCoeffs() const;

        //- Return the matrix source coefficients corresponding to the
        //- evaluation of the gradient of this patchField
        virtual tmp<Field<scalar>> gradientBoundaryCoeffs() const;

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;


    // Member operators

        virtual void operator=(const UList<scalar>&);

        virtual void operator=(const fvPatchField<scalar>&);
        virtual void operator+=(const fvPatchField<scalar>&);
        virtual void operator-=(const fvPatchField<scalar>&);
        virtual void operator*=(const fvPatchField<scalar>&);
        virtual void operator/=(const fvPatchField<scalar>&);

        virtual void operator+=(const Field<scalar>&);
        virtual void operator-=(const Field<scalar>&);

        virtual void operator*=(const Field<scalar>&);
        virtual void operator/=(const Field<scalar>&);

        virtual void operator=(const scalar);
        virtual void operator+=(const scalar);
        virtual void operator-=(const scalar);
        virtual void operator*=(const scalar);
        virtual void operator/=(const scalar);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
