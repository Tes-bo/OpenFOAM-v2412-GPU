/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 Zeljko Tukovic, FSB Zagreb.
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
    Foam::freeSurfacePressureFvPatchScalarField

Group
    grpGenericBoundaryConditions

Description
    This boundary condition provides static pressure condition for p_rgh,
    calculated as:

        \f[
            p = pa - \vec{g} \cdot \vec{r}
        \f]

    where
    \vartable
        p       | Free surface modified pressure
        pa      | Free surface ambient pressure
        g       | acceleration due to gravity [m/s^2]
    \endtable

Usage
    \table
        Property     | Description             | Required    | Default value
        pa           | static ambient pressure | yes         | 0
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            freeSurfacePressure;
        pa              uniform 0;
        value           uniform 0; // optional initial value
    }
    \endverbatim

See also
    Foam::fixedValueFvPatchScalarField

SourceFiles
    freeSurfacePressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef freeSurfacePressureFvPatchScalarField_H
#define freeSurfacePressureFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
         Class freeSurfacePressureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class freeSurfacePressureFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
protected:

    // Protected data

        //- Ambient pressure
        scalarField pa_;

public:

    //- Runtime type information
    TypeName("freeSurfacePressure");


    // Constructors

        //- Construct from patch and internal field
        freeSurfacePressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        freeSurfacePressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  freeSurfacePressureFvPatchScalarField onto a new patch
        freeSurfacePressureFvPatchScalarField
        (
            const freeSurfacePressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        freeSurfacePressureFvPatchScalarField
        (
            const freeSurfacePressureFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        freeSurfacePressureFvPatchScalarField
        (
            const freeSurfacePressureFvPatchScalarField&,
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

            //- Return the ambient pressure
            const scalarField& pa() const
            {
                return pa_;
            }

            //- Return reference to the ambient pressure to allow adjustment
            scalarField& pa()
            {
                return pa_;
            }


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


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
