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
    Foam::MUCSheterogeneousRate

Group
    grpLagrangianIntermediateSurfaceReactionSubModels

Description
    Heteregeneous noncatalytic reaction MUCS approach.
    Reference:
        D. Papanastassiou and G. Bitsianes,
        Modelling of Heterogeneous Gas-Solid Reactions,
        Metallurgical Transsactions, 480. Volume 4. 1973

\*---------------------------------------------------------------------------*/

#ifndef MUCSheterogeneousRate_H
#define MUCSheterogeneousRate_H

#include "HeterogeneousReactingModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
template<class CloudType>
class MUCSheterogeneousRate;

/*---------------------------------------------------------------------------*\
                  Class MUCSheterogeneousRate Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class MUCSheterogeneousRate
:
    public HeterogeneousReactingModel<CloudType>
{
    // Private Data

        // Model constants

            //- Binary diffusivity [m2/s]
            const scalar D12_;

            //- Porosity of product layer
            const scalar epsilon_;

            //- Tortuosity
            const scalar gamma_;

            //- Spread factor in pellet
            const scalar sigma_;

            //- Effectiveness factor
            const scalar E_;

            //- Kinetic rate coefficient [m/s]
            const scalar A_;

            //- Effective areas [-]
            const scalar Aeff_;

            //- Kinetic rate activation energy [J/kmol]
            const scalar Ea_;

            //- Stoichiomatric solid reactant
            const scalar nuFuel_;

            //- Stoichiomatric oxydizer reactant
            const scalar nuOx_;

            //- Stoichiomatric solid product
            const scalar nuProd_;


        // Addressing

            //- O2 position in global list
            label O2GlobalId_;

            //- Fuel solid local Id
            label FuelLocalId_;

            //- Product solid local Id
            label ProdLocalId_;


        // Local copies of thermo properties

            //- Molecular weight of O2 [kg/kmol]
            scalar WO2_;



public:

    //- Runtime type information
    TypeName("MUCSheterogeneousRate");


    // Constructors

        //- Construct from dictionary
        MUCSheterogeneousRate
        (
            const dictionary& dict,
            CloudType& owner
        );

        //- Construct copy
        MUCSheterogeneousRate
        (
            const MUCSheterogeneousRate<CloudType>& srm
        );

        //- Construct and return a clone
        virtual autoPtr<HeterogeneousReactingModel<CloudType>> clone() const
        {
            return autoPtr<HeterogeneousReactingModel<CloudType>>
            (
                new MUCSheterogeneousRate<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~MUCSheterogeneousRate() = default;


    // Member Functions

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
    #include "MUCSheterogeneousRate.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
