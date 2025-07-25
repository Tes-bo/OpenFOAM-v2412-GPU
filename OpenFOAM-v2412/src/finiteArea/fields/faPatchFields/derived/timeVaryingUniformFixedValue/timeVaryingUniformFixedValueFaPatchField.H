/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
    Copyright (C) 2019-2020 OpenCFD Ltd.
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
    Foam::timeVaryingUniformFixedValueFaPatchField

Description
    A time-varying form of a uniform fixed value finite area
    boundary condition.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        // Mandatory entries (unmodifiable)
        type                timeVaryingUniformFixedValue;
        fileName            "<case>/time-series";
        outOfBounds         clamp;           // (error|warn|clamp|repeat)

        // Mandatory/Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property | Description                             | Type  | Reqd | Dflt
      type     | Type name: timeVaryingUniformFixedValue | word  | yes  | -
      fileName | Name of operand file                    | word  | yes  | -
    \endtable

    The inherited entries are elaborated in:
      - \link faPatchFields.H \endlink

Author
    Zeljko Tukovic, FMENA
    Hrvoje Jasak, Wikki Ltd.

Note
    This class is derived directly from a fixedValue patch rather than from
    a uniformFixedValue patch.

See also
  - Foam::interpolationTable
  - Foam::fixedValueFaPatchField

SourceFiles
    timeVaryingUniformFixedValueFaPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef timeVaryingUniformFixedValueFaPatchField_H
#define timeVaryingUniformFixedValueFaPatchField_H

#include "fixedValueFaPatchField.H"
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class timeVaryingUniformFixedValueFaPatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class timeVaryingUniformFixedValueFaPatchField
:
    public fixedValueFaPatchField<Type>
{
    // Private Data

        //- The time series being used, including the bounding treatment
        interpolationTable<Type> timeSeries_;


public:

    //- Runtime type information
    TypeName("timeVaryingUniformFixedValue");


    // Constructors

        //- Construct from patch and internal field
        timeVaryingUniformFixedValueFaPatchField
        (
            const faPatch&,
            const DimensionedField<Type, areaMesh>&
        );

        //- Construct from patch, internal field and dictionary
        timeVaryingUniformFixedValueFaPatchField
        (
            const faPatch&,
            const DimensionedField<Type, areaMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patch field onto a new patch
        timeVaryingUniformFixedValueFaPatchField
        (
            const timeVaryingUniformFixedValueFaPatchField<Type>&,
            const faPatch&,
            const DimensionedField<Type, areaMesh>&,
            const faPatchFieldMapper&
        );

        //- Construct as copy
        timeVaryingUniformFixedValueFaPatchField
        (
            const timeVaryingUniformFixedValueFaPatchField<Type>&
        );

        //- Construct as copy setting internal field reference
        timeVaryingUniformFixedValueFaPatchField
        (
            const timeVaryingUniformFixedValueFaPatchField<Type>&,
            const DimensionedField<Type, areaMesh>&
        );

        //- Return clone
        virtual tmp<faPatchField<Type>> clone() const
        {
            return faPatchField<Type>::Clone(*this);
        }

        //- Clone with an internal field reference
        virtual tmp<faPatchField<Type>> clone
        (
            const DimensionedField<Type, areaMesh>& iF
        ) const
        {
            return faPatchField<Type>::Clone(*this, iF);
        }


    //- Destructor
    virtual ~timeVaryingUniformFixedValueFaPatchField() = default;


    // Member Functions

        // Access

            //- Return the time series used
            const interpolationTable<Type>& timeSeries() const
            {
                return timeSeries_;
            }


        // Evaluation

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        // IO

            //- Write includes "value" entry
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "timeVaryingUniformFixedValueFaPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
