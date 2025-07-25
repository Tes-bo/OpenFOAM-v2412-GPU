/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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
    Foam::NonInertialFrameForce

Group
    grpLagrangianIntermediateForceSubModels

Description
    Calculates particle non-inertial reference frame force.  Variable names as
    from Landau and Lifshitz, Mechanics, 3rd Ed, p126-129.

SourceFiles
    NonInertialFrameForce.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_NonInertialFrameForce_H
#define Foam_NonInertialFrameForce_H

#include "ParticleForce.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class NonInertialFrameForce Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class NonInertialFrameForce
:
    public ParticleForce<CloudType>
{
    // Private Data

        //- Name of the linear acceleration field
        word WName_;

        //- The linear acceleration of the reference frame
        vector W_;

        //- Name of the angular velocity field
        word omegaName_;

        //- The angular velocity of the reference frame
        vector omega_;

        //- Name of the angular acceleration field
        word omegaDotName_;

        //- The angular acceleration of the reference frame
        vector omegaDot_;

        //- Name of the centre of rotation field
        word centreOfRotationName_;

        //- The centre of rotation of the reference frame
        vector centreOfRotation_;


public:

    //- Runtime type information
    TypeName("nonInertialFrame");


    // Constructors

        //- Construct from mesh
        NonInertialFrameForce
        (
            CloudType& owner,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Construct copy
        NonInertialFrameForce(const NonInertialFrameForce& niff);

        //- Construct and return a clone
        virtual autoPtr<ParticleForce<CloudType>> clone() const
        {
            return autoPtr<ParticleForce<CloudType>>
            (
                new ParticleForce<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~NonInertialFrameForce();


    // Member Functions

    // Access

        //- Return the linear acceleration of the reference frame
        const vector& W() const noexcept { return W_; }

        //- Return the angular velocity of the reference frame
        const vector& omega() const noexcept { return omega_; }

        //- Return the angular acceleration of the reference frame
        const vector& omegaDot() const noexcept { return omegaDot_; }

        //- Return the centre of rotation of the reference frame
        const vector& centreOfRotation() const noexcept
        {
            return centreOfRotation_;
        }


        // Evaluation

            //- Cache fields
            virtual void cacheFields(const bool store);

            //- Calculate the non-coupled force
            virtual forceSuSp calcNonCoupled
            (
                const typename CloudType::parcelType& p,
                const typename CloudType::parcelType::trackingData& td,
                const scalar dt,
                const scalar mass,
                const scalar Re,
                const scalar muc
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "NonInertialFrameForce.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
