/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 PCOpt/NTUA
    Copyright (C) 2020 FOSS GP
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
    Foam::adjointRotatingWallVelocityFvPatchVectorField

Description
    The same as adjointWallVelocity but additionally computes the sensitivity
    contribution emerging from the change in the positions of the face centres,
    in case rotatingWallVelocity is used for the flow simulation.

SourceFiles
    adjointRotatingWallVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef adjointRotatingWallVelocityFvPatchVectorField_H
#define adjointRotatingWallVelocityFvPatchVectorField_H

#include "adjointWallVelocityFvPatchVectorField.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class adjointRotatingWallVelocity Declaration
\*---------------------------------------------------------------------------*/

class adjointRotatingWallVelocityFvPatchVectorField
:
    public adjointWallVelocityFvPatchVectorField
{
private:

    // Private Data

        //- Origin of the rotation
        vector origin_;

        //- Axis of the rotation
        vector axis_;

        //- Rotational speed
        autoPtr<Function1<scalar>> omega_;


public:

    //- Runtime type information
    TypeName("adjointRotatingWallVelocity");


    // Constructors

        //- Construct from patch and internal field
        adjointRotatingWallVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        adjointRotatingWallVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //- adjointRotatingWallVelocityFvPatchVectorField onto a new patch
        adjointRotatingWallVelocityFvPatchVectorField
        (
            const adjointRotatingWallVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        adjointRotatingWallVelocityFvPatchVectorField
        (
            const adjointRotatingWallVelocityFvPatchVectorField&,
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


    // Member functions

        //- Compute contribution to SDs
        virtual tmp<tensorField> dxdbMult() const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
