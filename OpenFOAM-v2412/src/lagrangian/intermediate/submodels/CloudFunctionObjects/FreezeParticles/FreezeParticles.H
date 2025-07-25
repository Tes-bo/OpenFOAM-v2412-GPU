/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 OpenCFD Ltd.
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
    Foam::FreezeParticles

Group
    grpLagrangianIntermediateFunctionObjects

Description
    Sets particle active flag to false to freeze particles in space

Usage
    Minimal example by using \c constant/\<CloudProperties\>:
    \verbatim
    cloudFunctions
    {
        FreezeParticles
        {
            // Mandatory entries
            type    FreezeParticles;
        }
    }
    \endverbatim

    where the entries mean:
    \table
      Property      | Description                       | Type   | Reqd | Deflt
      type          | Type name: FreezeParticles        | word   | yes  | -
    \endtable

SourceFiles
    FreezeParticles.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_FreezeParticles_H
#define Foam_FreezeParticles_H

#include "CloudFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class FreezeParticles Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class FreezeParticles
:
    public CloudFunctionObject<CloudType>
{
    // Private Member Functions

        //- Write post-processing info
        virtual void write()
        {}


public:

    typedef typename CloudType::parcelType parcelType;

    //- Runtime type information
    TypeName("FreezeParticles");


    // Constructors

        //- Construct from dictionary
        FreezeParticles
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );

        //- Copy construct
        FreezeParticles(const FreezeParticles<CloudType>& vf);

        //- Construct and return a clone
        virtual autoPtr<CloudFunctionObject<CloudType>> clone() const
        {
            return autoPtr<CloudFunctionObject<CloudType>>
            (
                new FreezeParticles<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~FreezeParticles() = default;


    // Member Functions

        //- Pre-evolve hook
        virtual void preEvolve(const typename parcelType::trackingData& td);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "FreezeParticles.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
