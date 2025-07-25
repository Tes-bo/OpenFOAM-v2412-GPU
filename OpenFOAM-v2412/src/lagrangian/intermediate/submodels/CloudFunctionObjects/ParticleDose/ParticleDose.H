/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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
    Foam::ParticleDose

Group
    grpLagrangianIntermediateFunctionObjects

Description
    Calculate the doses absorbed by a particle as the time integral
    of the particle track along the radiation field G [w/m2].

    Operands:
    \table
      Operand        | Type         | Location
      input          | -            | -
      output file    | -            | -
      output field   | scalarField  | \<time\>/lagrangian/\<cloud\>/D
    \endtable

Usage
    Minimal example by using \c constant/\<CloudProperties\>:
    \verbatim
    cloudFunctions
    {
        ParticleDose1
        {
            // Mandatory entries
            type             ParticleDose;
            GName            G;
        }
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                      | Type   | Reqd | Deflt
      type         | Type name: ParticleDose          | word   | yes  | -
      GName        | Name of the radiation field      | word   | yes  | -
    \endtable

SourceFiles
    ParticleDose.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_ParticleDose_H
#define Foam_ParticleDose_H

#include "CloudFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class ParticleDose Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ParticleDose
:
    public CloudFunctionObject<CloudType>
{
    // Private Data

        // Typedefs

            //- Convenience typedef for parcel type
            typedef typename CloudType::parcelType parcelType;


        //- Incident radiation field name
        word GName_;


public:

    //- Runtime type information
    TypeName("particleDose");


    // Generated Methods

        //- No copy assignment
        void operator=(const ParticleDose<CloudType>&) = delete;


    // Constructors

        //- Construct from dictionary
        ParticleDose
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );

        //- Copy construct
        ParticleDose(const ParticleDose<CloudType>& vf);

        //- Construct and return a clone
        virtual autoPtr<CloudFunctionObject<CloudType>> clone() const
        {
            return autoPtr<CloudFunctionObject<CloudType>>
            (
                new ParticleDose<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~ParticleDose() = default;


    // Member Functions

        //- Post-evolve hook
        virtual void postEvolve(const typename parcelType::trackingData& td);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ParticleDose.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
