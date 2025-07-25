/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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
    Foam::VoidFraction

Group
    grpLagrangianIntermediateFunctionObjects

Description
    Creates particle void fraction field on carrier phase

SourceFiles
    VoidFraction.C

\*---------------------------------------------------------------------------*/

#ifndef VoidFraction_H
#define VoidFraction_H

#include "CloudFunctionObject.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class VoidFraction Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class VoidFraction
:
    public CloudFunctionObject<CloudType>
{
    // Private Data

        // Typedefs

            //- Convenience typedef for parcel type
            typedef typename CloudType::parcelType parcelType;


        //- Void fraction field
        autoPtr<volScalarField> thetaPtr_;


protected:

    // Protected Member Functions

        //- Write post-processing info
        virtual void write();


public:

    //- Runtime type information
    TypeName("voidFraction");


    // Constructors

        //- Construct from dictionary
        VoidFraction
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );

        //- Construct copy
        VoidFraction(const VoidFraction<CloudType>& vf);

        //- Construct and return a clone
        virtual autoPtr<CloudFunctionObject<CloudType>> clone() const
        {
            return autoPtr<CloudFunctionObject<CloudType>>
            (
                new VoidFraction<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~VoidFraction() = default;


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
    #include "VoidFraction.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
