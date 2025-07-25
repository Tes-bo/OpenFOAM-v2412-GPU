/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2017 OpenFOAM Foundation
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
    Foam::ParticleTrap

Group
    grpLagrangianIntermediateFunctionObjects

Description
    Traps particles within a given phase fraction for multi-phase cases.

    Model is activated using:
    \verbatim
    particleTrap1
    {
        type        particleTrap;
        alphaName   alpha;      // name volume fraction field
        threshold   0.95;       // alpha value below which model is active
    }
    \endverbatim


SourceFiles
    ParticleTrap.C

\*---------------------------------------------------------------------------*/

#ifndef ParticleTrap_H
#define ParticleTrap_H

#include "CloudFunctionObject.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class ParticleTrap Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ParticleTrap
:
    public CloudFunctionObject<CloudType>
{
    // Private Data

        // Typedefs

            //- Convenience typedef for parcel type
            typedef typename CloudType::parcelType parcelType;


        //- Name of vol fraction field
        const word alphaName_;

        //- Pointer to the volume fraction field
        const volScalarField* alphaPtr_;

        //- Gradient of the volume fraction field
        autoPtr<volVectorField> gradAlphaPtr_;

        //- Threshold beyond which model is active
        scalar threshold_;


public:

    //- Runtime type information
    TypeName("particleTrap");


    // Constructors

        //- Construct from dictionary
        ParticleTrap
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );

        //- Construct copy
        ParticleTrap(const ParticleTrap<CloudType>& pe);

        //- Construct and return a clone
        virtual autoPtr<CloudFunctionObject<CloudType>> clone() const
        {
            return autoPtr<CloudFunctionObject<CloudType>>
            (
                new ParticleTrap<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~ParticleTrap() = default;


    // Member Functions

        // Evaluation

            //- Pre-evolve hook
            virtual void preEvolve
            (
                const typename parcelType::trackingData& td
            );

            //- Post-evolve hook
            virtual void postEvolve
            (
                const typename parcelType::trackingData& td
            );

            //- Post-move hook
            virtual bool postMove
            (
                parcelType& p,
                const scalar dt,
                const point& position0,
                const typename parcelType::trackingData& td
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ParticleTrap.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
