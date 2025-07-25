/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
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
    Foam::solidBodyMotionDisplacementPointPatchVectorField

Description
    Enables the specification of a fixed value boundary condition using the
    solid body motion functions.

SourceFiles
    solidBodyMotionDisplacementPointPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef solidBodyMotionDisplacementPointPatchVectorField_H
#define solidBodyMotionDisplacementPointPatchVectorField_H

#include "solidBodyMotionFunction.H"
#include "fixedValuePointPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class solidBodyMotionDisplacementPointPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class solidBodyMotionDisplacementPointPatchVectorField
:
    public fixedValuePointPatchVectorField
{
    // Private Data

        //- The motion control function
        autoPtr<solidBodyMotionFunction> SBMFPtr_;

        mutable autoPtr<pointField> localPoints0Ptr_;


public:

    //- Runtime type information
    TypeName("solidBodyMotionDisplacement");


    // Constructors

        //- Construct from patch and internal field
        solidBodyMotionDisplacementPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&
        );

        //- Construct from patch, internal field and dictionary
        solidBodyMotionDisplacementPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patch field onto a new patch
        solidBodyMotionDisplacementPointPatchVectorField
        (
            const solidBodyMotionDisplacementPointPatchVectorField&,
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const pointPatchFieldMapper&
        );

        //- Construct as copy
        solidBodyMotionDisplacementPointPatchVectorField
        (
            const solidBodyMotionDisplacementPointPatchVectorField&
        );

        //- Construct as copy setting internal field reference
        solidBodyMotionDisplacementPointPatchVectorField
        (
            const solidBodyMotionDisplacementPointPatchVectorField&,
            const DimensionedField<vector, pointMesh>&
        );


        //- Return a clone
        virtual autoPtr<pointPatchField<vector>> clone() const
        {
            return pointPatchField<vector>::Clone(*this);
        }

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<pointPatchField<vector>> clone
        (
            const DimensionedField<vector, pointMesh>& iF
        ) const
        {
            return pointPatchField<vector>::Clone(*this, iF);
        }


    // Member functions

        // Access

            //- Return the fluctuation scale
            const solidBodyMotionFunction& motion() const
            {
                return *SBMFPtr_;
            }

            const pointField& localPoints0() const;


        // Evaluation Functions

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
