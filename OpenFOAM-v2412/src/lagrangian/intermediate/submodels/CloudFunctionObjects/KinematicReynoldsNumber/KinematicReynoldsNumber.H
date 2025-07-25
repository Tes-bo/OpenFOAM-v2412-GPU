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
    Foam::KinematicReynoldsNumber

Group
    grpLagrangianIntermediateFunctionObjects

Description
    Calculates and writes particle Reynolds number field on the cloud.
    The normalisation factors are calculated without thermal effects.

    \f[
        \mathrm{Re}_p =
            \frac{\rho_c \, | \mathbf{u}_\mathrm{rel} | \, d_p}{\mu_c}
    \f]

    \vartable
        \mathrm{Re}_p  | Particle Reynolds number
        d_p            | Particle diameter
        \rho_c         | Density of carrier
        \mu_c          | Dynamic viscosity of carrier
        \mathbf{u}_\mathrm{rel} | Relative velocity between particle and carrier
    \endvartable

    Operands:
    \table
      Operand        | Type         | Location
      input          | -            | -
      output file    | -            | -
      output field   | scalarField  | \<time\>/lagrangian/\<cloud\>/Re
    \endtable

Usage
    Minimal example by using \c constant/\<CloudProperties\>:
    \verbatim
    cloudFunctions
    {
        KinematicReynoldsNumber1
        {
            // Mandatory entries
            type             ReynoldsNumber;
        }
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                 | Type   | Reqd | Deflt
      type         | Type name: ReynoldsNumber   | word   | yes  | -
    \endtable

See also
  - Foam::ThermoReynoldsNumber

SourceFiles
    KinematicReynoldsNumber.C

\*---------------------------------------------------------------------------*/

#ifndef KinematicReynoldsNumber_H
#define KinematicReynoldsNumber_H

#include "CloudFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class KinematicReynoldsNumber Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class KinematicReynoldsNumber
:
    public CloudFunctionObject<CloudType>
{
    // Private Data

        // Typedefs

            //- Convenience typedef for parcel type
            typedef typename CloudType::parcelType parcelType;


public:

    //- Runtime type information
    TypeName("ReynoldsNumber");


    // Generated Methods

        //- No copy assignment
        void operator=(const KinematicReynoldsNumber<CloudType>&) = delete;


    // Constructors

        //- Construct from dictionary
        KinematicReynoldsNumber
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );

        //- Copy construct
        KinematicReynoldsNumber(const KinematicReynoldsNumber<CloudType>& vf);

        //- Construct and return a clone
        virtual autoPtr<CloudFunctionObject<CloudType>> clone() const
        {
            return autoPtr<CloudFunctionObject<CloudType>>
            (
                new KinematicReynoldsNumber<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~KinematicReynoldsNumber() = default;


    // Member Functions

        //- Post-evolve hook
        virtual void postEvolve(const typename parcelType::trackingData& td);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "KinematicReynoldsNumber.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
