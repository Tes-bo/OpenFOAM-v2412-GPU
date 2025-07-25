/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::PilchErdman

Group
    grpLagrangianIntermediateBreakupSubModels

Description
    Particle secondary breakup model based on
    Pilch-Erdman total droplet breakup model.

    Reference:
    \verbatim
        Pilch, M., & Erdman, C. A. (1987).
        Use of breakup time data and velocity history data
        to predict the maximum size of stable fragments for
        acceleration-induced breakup of a liquid drop.
        International journal of multiphase flow, 13(6), 741-757.
        DOI:10.1016/0301-9322(87)90063-2
    \endverbatim

    The droplet fragment velocity is described by the equation:

    \f[
        V_d = V \sqrt(\epsilon)(B_1 T + B_2 T^2)
    \f]

    where
    \vartable
        V_d      | Fragment velocity
        V        | Magnitude of the relative velocity
        \epsilon | Density ratio (rho_carrier/rho_droplet)
        T        | characteristic break-up time
        B_1      | Model coefficient
        B_2      | Model coefficient
    \endvartable
    Where:

    The authors suggest that:
        compressible flow   : B1 = 0.75*1.0; B2 = 3*0.116
        incompressible flow : B1 = 0.75*0.5; B2 = 3*0.0758

SourceFiles
    PilchErdman.C

\*---------------------------------------------------------------------------*/

#ifndef PilchErdman_H
#define PilchErdman_H

#include "BreakupModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class PilchErdman Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class PilchErdman
:
    public BreakupModel<CloudType>
{
    // Private Data

        //- Model coefficient
        scalar B1_;

        //- Model coefficient
        scalar B2_;


public:

    //- Runtime type information
    TypeName("PilchErdman");


    // Constructors

        //- Construct from dictionary
        PilchErdman(const dictionary&, CloudType&);

        //- Construct copy
        PilchErdman(const PilchErdman<CloudType>& bum);

        //- Construct and return a clone
        virtual autoPtr<BreakupModel<CloudType>> clone() const
        {
            return autoPtr<BreakupModel<CloudType>>
            (
                new PilchErdman<CloudType>(*this)
            );
        }

        //- No copy assignment
        void operator=(const PilchErdman<CloudType>&) = delete;


    //- Destructor
    virtual ~PilchErdman() = default;


    // Member Functions

        //- Update the parcel properties
        virtual bool update
        (
            const scalar dt,
            const vector& g,
            scalar& d,
            scalar& tc,
            scalar& ms,
            scalar& nParticle,
            scalar& KHindex,
            scalar& y,
            scalar& yDot,
            const scalar d0,
            const scalar rho,
            const scalar mu,
            const scalar sigma,
            const vector& U,
            const scalar rhoc,
            const scalar muc,
            const vector& Urel,
            const scalar Urmag,
            const scalar tMom,
            scalar& dChild,
            scalar& massChild
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PilchErdman.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
