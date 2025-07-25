/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2019 OpenCFD Ltd.
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
    Foam::NoheterogeneousReacting

Group
    grpLagrangianIntermediateSurfaceReactionSubModels

Description
    Dummy surface reaction model for 'none'

\*---------------------------------------------------------------------------*/

#ifndef NoheterogeneousReacting_H
#define NoheterogeneousReacting_H

#include "HeterogeneousReactingModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class NoheterogeneousReacting Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class NoheterogeneousReacting
:
    public HeterogeneousReactingModel<CloudType>
{
public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from dictionary
        NoheterogeneousReacting(const dictionary& dict, CloudType& owner);

        //- Construct copy
        NoheterogeneousReacting
        (
            const NoheterogeneousReacting<CloudType>& srm
        );

        //- Construct and return a clone
        virtual autoPtr<HeterogeneousReactingModel<CloudType>> clone() const
        {
            return autoPtr<HeterogeneousReactingModel<CloudType>>
            (
                new NoheterogeneousReacting<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~NoheterogeneousReacting() = default;


    // Member Functions

        //- Flag to indicate whether model activates devolatisation model
        virtual bool active() const;

        //- Update surface reactions
        virtual scalar calculate
        (
            const scalar dt,
            const scalar Re,
            const scalar nu,
            const label celli,
            const scalar d,
            const scalar T,
            const scalar Tc,
            const scalar pc,
            const scalar rhoc,
            const scalar mass,
            const scalarField& YSolid,
            scalarField& F,
            const scalar N,
            scalar& NCpW,
            scalarField& dMassSolid,
            scalarField& dMassSRCarrier
        ) const;

        //- Number of reactions in the model
        virtual label nReactions() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "NoheterogeneousReacting.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
