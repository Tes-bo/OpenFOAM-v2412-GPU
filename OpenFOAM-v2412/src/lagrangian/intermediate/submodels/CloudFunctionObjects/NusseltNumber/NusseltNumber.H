/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
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
    Foam::NusseltNumber

Group
    grpLagrangianIntermediateFunctionObjects

Description
    Calculates and writes particle Nusselt number field on the cloud.

    Operands:
    \table
      Operand        | Type         | Location
      input          | -            | -
      output file    | -            | -
      output field   | scalarField  | \<time\>/lagrangian/\<cloud\>/Nu
    \endtable

Usage
    Minimal example by using \c constant/<CloudProperties>:
    \verbatim
    cloudFunctionObjects
    {
        NusseltNumber1
        {
            // Mandatory entries
            type             NusseltNumber;
        }
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description               | Type   | Reqd | Deflt
      type         | Type name: NusseltNumber  | word   | yes  | -
    \endtable

SourceFiles
    NusseltNumber.C

\*---------------------------------------------------------------------------*/

#ifndef NusseltNumber_H
#define NusseltNumber_H

#include "CloudFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class NusseltNumber Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class NusseltNumber
:
    public CloudFunctionObject<CloudType>
{
    // Private Data

        // Typedefs

            //- Convenience typedef for parcel type
            typedef typename CloudType::parcelType parcelType;


public:

    //- Runtime type information
    TypeName("NusseltNumber");


    // Generated Methods

        //- No copy assignment
        void operator=(const NusseltNumber<CloudType>&) = delete;


    // Constructors

        //- Construct from dictionary
        NusseltNumber
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );

        //- Copy construct
        NusseltNumber(const NusseltNumber<CloudType>& vf);

        //- Construct and return a clone
        virtual autoPtr<CloudFunctionObject<CloudType>> clone() const
        {
            return autoPtr<CloudFunctionObject<CloudType>>
            (
                new NusseltNumber<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~NusseltNumber() = default;


    // Member Functions

        //- Post-evolve hook
        virtual void postEvolve(const typename parcelType::trackingData& td);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "NusseltNumber.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
