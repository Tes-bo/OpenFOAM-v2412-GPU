/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016 OpenCFD Ltd
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
    Foam::radiation::greyDiffusiveViewFactorFixedValueFvPatchScalarField

Group
    grpThermoBoundaryConditions

Description
    This boundary condition provides a grey-diffuse condition for radiative
    heat flux, \c qr, for use with the view factor model

Usage
    \table
        Property     | Description             | Required    | Default value
        qro          | external radiative heat flux | yes    |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            greyDiffusiveRadiationViewFactor;
        qro             uniform 0;
    }
    \endverbatim

See also
    Foam::radiation::radiationModel
    Foam::radiation::viewFactor
    Foam::radiationCoupledBase
    Foam::fixedValueFvPatchField

SourceFiles
    greyDiffusiveViewFactorFixedValueFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef greyDiffusiveViewFactorFixedValueFvPatchScalarField_H
#define greyDiffusiveViewFactorFixedValueFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
      Class greyDiffusiveViewFactorFixedValueFvPatchScalarField declaration
\*---------------------------------------------------------------------------*/

class greyDiffusiveViewFactorFixedValueFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

        //- External radiative heat flux
        scalarField qro_;


public:

    //- Runtime type information
    TypeName("greyDiffusiveRadiationViewFactor");


    // Constructors

        //- Construct from patch and internal field
        greyDiffusiveViewFactorFixedValueFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        greyDiffusiveViewFactorFixedValueFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given a
        // greyDiffusiveViewFactorFixedValueFvPatchScalarField onto a new patch
        greyDiffusiveViewFactorFixedValueFvPatchScalarField
        (
            const greyDiffusiveViewFactorFixedValueFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        greyDiffusiveViewFactorFixedValueFvPatchScalarField
        (
            const greyDiffusiveViewFactorFixedValueFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        greyDiffusiveViewFactorFixedValueFvPatchScalarField
        (
            const greyDiffusiveViewFactorFixedValueFvPatchScalarField&,
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

        // Access

            //- Return external + solar load radiative heat flux
            tmp<scalarField> qro(label bandI = 0) const;


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        // Mapping functions

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


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
