/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2019 IH-Cantabria
    Copyright (C) 2018-2019 OpenCFD Ltd.
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
    Foam::waveMakerPointPatchVectorField

Description
    Point motion boundary condition to generate waves based on either piston
    or flap motions.

    Based on the reference
    \verbatim
        Hughes, S.A. (1993).
        Physical Models And Laboratory Techniques In Coastal Engineering.
        Advanced Series On Ocean Engineering, volume 7
    \endverbatim

Usage
    Example patch specification
    \verbatim
    leftwall
    {
        type            waveMaker;
        motionType      flap;
        n               (1 0 0);
        initialDepth    0.25;
        wavePeriod      2.0;
        waveHeight      0.06;
        wavePhase       0;
        rampTime        2.0;
    }
    \endverbatim

    where
    \table
        Property     | Description                  | Required | Default value
        motionType   | See motion types below       | yes      |
        n            | Direction of motion          | yes      |
        initialDepth | Initial depth                | yes      |
        wavePeriod   | wave period                  | yes      |
        waveHeight   | Wave height                  | yes      |
        wavePhase    | wave phase                   | yes      |
        waveAngle    | wave angle                   | no       | 0
        startTime    | Start time                   | no       | case start time
        rampTime     | Time to reach maximum motion | yes      |
        secondOrder  | Second order calculation     | no       | no
        nPaddle      | Number of paddles            | no       | 1
    \endtable

    Available motion types include:
    - piston
    - flap
    - solitary

SourceFiles
    waveMakerPointPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef waveMakerPointPatchVectorField_H
#define waveMakerPointPatchVectorField_H

#include "fixedValuePointPatchField.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class waveMakerPointPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class waveMakerPointPatchVectorField
:
    public fixedValuePointPatchField<vector>
{
    public:

        enum motionTypes
        {
            piston,
            flap,
            solitary
        };

        //- Names for motion types
        static const Enum<motionTypes> motionTypeNames;


    // Private data

        //- Motion type
        motionTypes motionType_;

        //- Patch normal direction
        //  Note: cannot use patch normal of the initial patch unless it is
        //  in its neutral position (flap mode)
        vector n_;

        //- Vertical direction
        vector gHat_;

        //- Initial water depth
        scalar initialDepth_;

        //- Wave period
        scalar wavePeriod_;

        //- Wave height
        scalar waveHeight_;

        //- Wave phase
        scalar wavePhase_;

        //- Wave angle
        scalar waveAngle_;

        //- Wave length
        scalar waveLength_;

        //- Start time
        scalar startTime_;

        //- Ramp time
        scalar rampTime_;

        //- On/off second order calculation switch
        scalar secondOrder_;

        //- Number of wave paddles
        label nPaddle_;

        //- Rotation tensor from global to local system
        tensor Rgl_;

        //- Rotation tensor from local to global system
        tensor Rlg_;

        //- Paddle x co-ordinates / [m]
        scalarField xPaddle_;

        //- Paddle y co-ordinates / [m]
        scalarField yPaddle_;

        //- Addressing from point patch index to paddle index
        labelList pointToPaddle_;

        //- Addressing from patch face index to paddle index
        labelList faceToPaddle_;

        //- Patch face centre x co-ordinates / [m]
        scalarField x_;

        //- Patch face centre y co-ordinates / [m]
        scalarField y_;

        //- Patch face centre z co-ordinates / [m]
        scalarField z_;

        //- Overall (point) span in z-direction / [m]
        scalar zSpan_;

        //- Minimum z (point) height per patch face / [m]
        scalarField zMin_;

        //- Global Minimum z (point) / [m]
        scalar zMinGb_;

        //- Maximum z (point) height per patch face / [m]
        scalarField zMax_;

        //- Calculated water depth at the patch
        scalarField waterDepthRef_;

        //
        scalar firstTime = 0;


    // Protected Member Functions

        //- Return the gravitational acceleration
        const vector& g();

        //- Dispersion equation
        virtual scalar waveLength(const scalar h, const scalar T);

        //- Return the time scaling coefficient
        virtual scalar timeCoeff(const scalar t) const;

        //- Initialise
        virtual void initialiseGeometry();


public:

    //- Runtime type information
    TypeName("waveMaker");


    // Constructors

        //- Construct from patch and internal field
        waveMakerPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&
        );

        //- Construct from patch, internal field and dictionary
        waveMakerPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patch field onto a new patch
        waveMakerPointPatchVectorField
        (
            const waveMakerPointPatchVectorField&,
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const pointPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        waveMakerPointPatchVectorField
        (
            const waveMakerPointPatchVectorField&,
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
