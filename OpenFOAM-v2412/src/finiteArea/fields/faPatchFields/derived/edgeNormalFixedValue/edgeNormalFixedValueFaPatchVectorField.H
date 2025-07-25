/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
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
    Foam::edgeNormalFixedValueFaPatchVectorField

Description
    Edge normal fixed value vector field finite area boundary condition

    Describes a surface normal vector boundary condition by its magnitude.
    Note: The value is positive for outward-pointing vectors

Author
    Zeljko Tukovic, FMENA
    Hrvoje Jasak, Wikki Ltd.

SourceFiles
    edgeNormalFixedValueFaPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef edgeNormalFixedValueFaPatchVectorField_H
#define edgeNormalFixedValueFaPatchVectorField_H

#include "faPatchFields.H"
#include "fixedValueFaPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class edgeNormalFixedValueFaPatch Declaration
\*---------------------------------------------------------------------------*/

class edgeNormalFixedValueFaPatchVectorField
:
    public fixedValueFaPatchVectorField
{
    // Private Data

        //- Edge-normal velocity value
        scalarField refValue_;


public:

    //- Runtime type information
    TypeName("edgeNormalFixedValue");


    // Constructors

        //- Construct from patch and internal field
        edgeNormalFixedValueFaPatchVectorField
        (
            const faPatch&,
            const DimensionedField<vector, areaMesh>&
        );

        //- Construct from patch, internal field and reference value
        edgeNormalFixedValueFaPatchVectorField
        (
            const faPatch&,
            const DimensionedField<vector, areaMesh>&,
            const scalar refValue
        );

        //- Construct from patch, internal field and dictionary
        edgeNormalFixedValueFaPatchVectorField
        (
            const faPatch&,
            const DimensionedField<vector, areaMesh>&,
            const dictionary&
        );

        //- Construct by mapping given edgeNormalFixedValueFaPatchVectorField
        //- onto a new patch
        edgeNormalFixedValueFaPatchVectorField
        (
            const edgeNormalFixedValueFaPatchVectorField&,
            const faPatch&,
            const DimensionedField<vector, areaMesh>&,
            const faPatchFieldMapper&
        );

        //- Construct as copy
        edgeNormalFixedValueFaPatchVectorField
        (
            const edgeNormalFixedValueFaPatchVectorField&
        );

        //- Construct as copy setting internal field reference
        edgeNormalFixedValueFaPatchVectorField
        (
            const edgeNormalFixedValueFaPatchVectorField&,
            const DimensionedField<vector, areaMesh>&
        );

        //- Return clone
        virtual tmp<faPatchField<vector>> clone() const
        {
            return faPatchField<vector>::Clone(*this);
        }

        //- Clone with an internal field reference
        virtual tmp<faPatchField<vector>> clone
        (
            const DimensionedField<vector, areaMesh>& iF
        ) const
        {
            return faPatchField<vector>::Clone(*this, iF);
        }


    // Member functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const faPatchFieldMapper&
            );

            //- Reverse map the given faPatchField onto this faPatchField
            virtual void rmap
            (
                const faPatchVectorField&,
                const labelList&
            );


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write includes "value" entry
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
