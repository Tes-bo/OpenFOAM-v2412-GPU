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
    Foam::adjointOutletVelocityFvPatchVectorField

Description
    Provides the adjoint outlet velocity values (i.e. adjoint velocity in
    case of a zeroGradient U boundary condition). Can have stability issues
    in cases of backflow of the primal velocity.
    The adjointOutletVelocityFlux should preferably be used.


SourceFiles
    adjointOutletVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef adjointOutletVelocityFvPatchVectorField_H
#define adjointOutletVelocityFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "adjointBoundaryConditions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class adjointOutletVelocityFvPatch Declaration
\*---------------------------------------------------------------------------*/

class adjointOutletVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField,
    public adjointVectorBoundaryCondition
{
    // Private Member Functions

        void assignBoundaryValue();


public:

    //- Runtime type information
    TypeName("adjointOutletVelocity");


    // Constructors

        //- Construct from patch and internal field
        adjointOutletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        adjointOutletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given adjointOutletVelocityFvPatchVectorField
        //- onto a new patch
        adjointOutletVelocityFvPatchVectorField
        (
            const adjointOutletVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        adjointOutletVelocityFvPatchVectorField
        (
            const adjointOutletVelocityFvPatchVectorField&,
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

        //- Return true: Allow adjoint solvers to obtain the outlet phia
        //  value through HbyA
        virtual bool assignable() const
        {
            return true;
        }

        //- Update the coefficients associated with the patch field
        //  Apply adjoint BCs through evaluate rather than updateCoeffs
        //  in order to have the correct Ua boundaryField when computing the
        //  adjoint pressure BC
        virtual void evaluate
        (
            const Pstream::commsTypes commsType = Pstream::commsTypes::buffered
        );

        //virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;


    // Member operators

        virtual void operator=(const fvPatchField<vector>& pvf);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
