/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::acousticWaveTransmissiveFvPatchField

Group
    grpOutletBoundaryConditions

Description
    This boundary condition provides a wave transmissive outflow condition,
    based on solving DDt(W, field) = 0 at the boundary \c W is the wave velocity
    and \c field is the field to which this boundary condition is applied.
    The wave speed is input in the BC.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        // Mandatory entries (unmodifiable)
        type                acousticWaveTransmissive;
        advectiveSpeed      50.0;

        // Mandatory/Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property | Description                          | Type  | Reqd | Dflt
      type     | Type name: acousticWaveTransmissive  | word  | yes  | -
      advectiveSpeed | Advective speed value          | scalar | yes | -
    \endtable

    The inherited entries are elaborated in:
      - \link advectiveFvPatchField.H \endlink

See also
  - Foam::advectiveFvPatchFields

SourceFiles
    acousticWaveTransmissiveFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef acousticWaveTransmissiveFvPatchField_H
#define acousticWaveTransmissiveFvPatchField_H

#include "advectiveFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class acousticWaveTransmissiveFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class acousticWaveTransmissiveFvPatchField
:
    public advectiveFvPatchField<Type>
{
    // Private Data

        //- Advection speed value
        scalar advectiveU_;


public:

    //- Runtime type information
    TypeName("acousticWaveTransmissive");


    // Constructors

        //- Construct from patch and internal field
        acousticWaveTransmissiveFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        acousticWaveTransmissiveFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given acousticWaveTransmissiveFvPatchField
        //- onto a new patch
        acousticWaveTransmissiveFvPatchField
        (
            const acousticWaveTransmissiveFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        acousticWaveTransmissiveFvPatchField
        (
            const acousticWaveTransmissiveFvPatchField&
        );

        //- Construct as copy setting internal field reference
        acousticWaveTransmissiveFvPatchField
        (
            const acousticWaveTransmissiveFvPatchField&,
            const DimensionedField<Type, volMesh>&
        );

        //- Return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return fvPatchField<Type>::Clone(*this);
        }

        //- Clone with an internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return fvPatchField<Type>::Clone(*this, iF);
        }


    // Member Functions

        // Evaluation

            //- Calculate and return the advection speed at the boundary
            virtual tmp<scalarField> advectionSpeed() const;


        // IO

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "acousticWaveTransmissiveFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
