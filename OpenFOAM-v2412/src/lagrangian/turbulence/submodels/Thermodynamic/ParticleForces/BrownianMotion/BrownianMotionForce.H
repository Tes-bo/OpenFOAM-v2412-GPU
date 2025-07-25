/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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
    Foam::BrownianMotionForce

Description
    Calculates particle Brownian motion force.

    Reference:
    \verbatim
        Li, A., & Ahmadi, G. (1992).
        Dispersion and deposition of spherical particles from point sources
        in a turbulent channel flow.
        Aerosol science and technology,
        16(4), 209-226.
    \endverbatim

    Entries
    \table
        Property | Description                            | Type | Reqd | Deflt
        lambda   | Molecular free path length             | scalar | yes  | -
        turbulence | Use turbulence                       | bool | yes | -
        spherical | Spherical instead of cubic distribution | bool | no  | true
    \endtable

Note
    The treatment changed from cubic to spherical around v1606.
    Reinstate support for cubic via spherical (true/false) switch after v2306.

SourceFiles
    BrownianMotionForce.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_BrownianMotionForce_H
#define Foam_BrownianMotionForce_H

#include "ParticleForce.H"
#include "Random.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class BrownianMotionForce Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class BrownianMotionForce
:
    public ParticleForce<CloudType>
{
    // Private Data

        //- Reference to the cloud random number generator
        Random& rndGen_;

        //- Molecular free path length [m]
        const scalar lambda_;

        //- Pointer to the turbulence kinetic energy field
        const volScalarField* kPtr_;

        //- Turbulence flag
        bool turbulence_;

        //- Local ownership of the turbulence k field
        bool ownK_;

        //- Spherical vs cubic (transitional)
        bool useSpherical_;


    // Private Member Functions

        //- Return the k field from the turbulence model
        tmp<volScalarField> kModel() const;


public:

    //- Runtime type information
    TypeName("BrownianMotion");


    // Constructors

        //- Construct from mesh
        BrownianMotionForce
        (
            CloudType& owner,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Construct copy
        BrownianMotionForce(const BrownianMotionForce& bmf);

        //- Construct and return a clone
        virtual autoPtr<ParticleForce<CloudType>> clone() const
        {
            return autoPtr<ParticleForce<CloudType>>
            (
                new BrownianMotionForce<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~BrownianMotionForce();


    // Member Functions

    // Access

        //- The molecular free path length [m]
        scalar lambda() const noexcept { return lambda_; }

        //- The turbulence flag
        bool turbulence() const noexcept { return turbulence_; }


        // Evaluation

            //- Cache fields
            virtual void cacheFields(const bool store);

            //- Calculate the coupled force
            virtual forceSuSp calcCoupled
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
    #include "BrownianMotionForce.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
