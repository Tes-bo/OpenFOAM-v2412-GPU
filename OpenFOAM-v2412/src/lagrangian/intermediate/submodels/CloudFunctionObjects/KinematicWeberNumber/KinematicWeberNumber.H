/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 OpenCFD Ltd.
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
    Foam::KinematicWeberNumber

Group
    grpLagrangianIntermediateFunctionObjects

Description
    Calculates and writes particle Weber number field on the cloud.

    \f[
        \mathrm{We}_p =
            \frac{\rho_c \, | \mathbf{u}_\mathrm{rel} |^2 \, d_p }{\sigma}
    \f]

    \vartable
        \mathrm{We}_p  | Particle Weber number
        \rho_c         | Density of carrier
        d_p            | Particle diameter
        \mathbf{u}_\mathrm{rel} | Relative velocity between particle and carrier
    \endvartable

    Operands:
    \table
      Operand        | Type         | Location
      input          | -            | -
      output file    | -            | -
      output field   | scalarField  | \<time\>/lagrangian/\<cloud\>/We
    \endtable

Usage
    Minimal example by using \c constant/\<CloudProperties\>:
    \verbatim
    cloudFunctions
    {
        KinematicWeberNumber1
        {
            // Mandatory entries
            type             WeberNumber;
            sigma            <scalar>;
        }
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                 | Type   | Reqd | Deflt
      type         | Type name: WeberNumber      | word   | yes  | -
      sigma        | Surface tension [N/m]       | scalar | yes  | -
    \endtable

See also
  - Foam::ReactingWeberNumber

SourceFiles
    KinematicWeberNumber.C

\*---------------------------------------------------------------------------*/

#ifndef KinematicWeberNumber_H
#define KinematicWeberNumber_H

#include "CloudFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class KinematicWeberNumber Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class KinematicWeberNumber
:
    public CloudFunctionObject<CloudType>
{
    // Private Typedefs

        //- Convenience typedef for parcel type
        typedef typename CloudType::parcelType parcelType;


    // Private Data

        //- Surface tension [N/m]
        scalar sigma_;


public:

    //- Runtime type information
    TypeName("WeberNumber");


    // Generated Methods

        //- No copy assignment
        void operator=(const KinematicWeberNumber<CloudType>&) = delete;


    // Constructors

        //- Construct from dictionary
        KinematicWeberNumber
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );

        //- Copy construct
        KinematicWeberNumber(const KinematicWeberNumber<CloudType>& vf);

        //- Construct and return a clone
        virtual autoPtr<CloudFunctionObject<CloudType>> clone() const
        {
            return autoPtr<CloudFunctionObject<CloudType>>
            (
                new KinematicWeberNumber<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~KinematicWeberNumber() = default;


    // Member Functions

        //- Post-evolve hook
        virtual void postEvolve(const typename parcelType::trackingData& td);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "KinematicWeberNumber.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
