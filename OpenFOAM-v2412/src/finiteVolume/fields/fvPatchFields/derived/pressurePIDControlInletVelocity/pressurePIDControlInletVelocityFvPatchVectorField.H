/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenCFD Ltd.
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
    Foam::pressurePIDControlInletVelocityFvPatchVectorField

Group
    grpInletBoundaryConditions

Description
    This boundary condition tries to generate an inlet velocity that maintains
    a specified pressure drop between two face zones downstream. The zones
    should fully span a duct through which all the inlet flow passes.

    An incompressible, lossless analysis of the flow between the inlet and the
    two face-zones is performed. An ideal inlet velocity is thereby calculated
    from the user-specified pressure drop. This analysis can include the
    transient effect of the inlet velocity change. In this case, a shape factor
    is included to represent non-linearity in the duct cross section.

    The average pressure drop between the two face zones is measured. The same
    incompressible, lossless analysis is performed using this pressure drop.
    The difference between the two computed velocities is considered as an
    error. The ideal inlet is modified so as to drive this error to zero. This
    is accomplished by means of a PID control algorithm, for which
    non-dimensional gains are specified by the user.

    The shape factor takes a value of 0 for a linear change in cross sectional
    area between the two face zones. A value of 1 represents a step change in
    area at the mid-point between the zones. A smooth cubic or cosine profile
    between two zones with zero divergence is typically represented by a factor
    of between 0.2 and 0.25.

    \heading Patch usage

    \table
        Property    | Description                 | Required | Default value
        upstream    | upstream face zone name     | yes      |
        downstream  | downstream face zone name   | yes      |
        deltaP      | desired pressure drop       | yes      |
        shapeFactor | non-linearity in the nozzle | no       | 0
        p           | pressure field name         | no       | p
        phi         | flux field name             | yes      | phi
        rho         | density field name          | no       | none
        P           | proportional gain           | yes      |
        I           | integral gain               | yes      |
        D           | differential gain           | yes      |
    \endtable

    Example of the boundary condition specification:

    \verbatim
    myPatch
    {
        type            pressurePIDControlInletVelocity;
        upstream        upstream;
        downstream      downstream;
        deltaP          200;
        shapeFactor     0;
        p               p;
        phi             phi;
        rho             none;
        P               0.5;
        I               0.5;
        D               0.1;
        value           uniform (0 0 0);
    }
    \endverbatim

SeeAlso
    Foam::fixedValueFvPatchField

SourceFiles
    pressurePIDControlInletVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef pressurePIDControlInletVelocityFvPatchVectorField_H
#define pressurePIDControlInletVelocityFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
     Class pressurePIDControlInletVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class pressurePIDControlInletVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- Name of the upstream face zone
        const word upstreamName_;

        //- Name of the downstream face zone
        const word downstreamName_;

        //- Desired pressure difference between upstream and downstream
        const scalar deltaP_;

        //- Nozzle shape factor
        const scalar shapeFactor_;

        //- Name of the pressure field
        const word pName_;

        //- Name of the flux field
        const word phiName_;

        //- Name of the density field (if any)
        const word rhoName_;

        //- Proportional gain
        const scalar P_;

        //- Integral gain
        const scalar I_;

        //- Derivative gain
        const scalar D_;

        //- Volumetric flow rate
        scalar Q_;

        //- Error
        scalar error_;

        //- Error integral w.r.t. time
        scalar errorIntegral_;

        //- Old volumetric flow rate
        scalar oldQ_;

        //- Old error
        scalar oldError_;

        //- Old error integral w.r.t. time
        scalar oldErrorIntegral_;

        //- Time index of the last update
        label timeIndex_;


    // Private member functions

        //- Return the pressure interpolated to the faces
        const surfaceScalarField& facePressure() const;

        //- Calculate the average on a face zone
        template <class Type>
        void faceZoneAverage
        (
            const word& name,
            const GeometricField<Type, fvsPatchField, surfaceMesh>& field,
            scalar& area,
            Type& average
        ) const;


public:

   //- Runtime type information
   TypeName("pressurePIDControlInletVelocity");


   // Constructors

        //- Construct from patch and internal field
        pressurePIDControlInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        pressurePIDControlInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  flowRateInletVelocityFvPatchVectorField
        //  onto a new patch
        pressurePIDControlInletVelocityFvPatchVectorField
        (
            const pressurePIDControlInletVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        pressurePIDControlInletVelocityFvPatchVectorField
        (
            const pressurePIDControlInletVelocityFvPatchVectorField&
        );

        //- Construct as copy setting internal field reference
        pressurePIDControlInletVelocityFvPatchVectorField
        (
            const pressurePIDControlInletVelocityFvPatchVectorField&,
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

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "pressurePIDControlInletVelocityFvPatchVectorFieldTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
