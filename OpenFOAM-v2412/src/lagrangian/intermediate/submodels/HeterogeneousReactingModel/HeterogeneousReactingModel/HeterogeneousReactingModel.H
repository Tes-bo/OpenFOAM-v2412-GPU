/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenCFD Ltd.
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
    Foam::HeterogeneousReactingModel

Group
    grpLagrangianIntermediateHeterogeneousReactionSubModels

Description
    Base class for heterogeneous reacting models

SourceFiles
    HeterogeneousReactingModel.C
    HeterogeneousReactingModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef HeterogeneousReactingModel_H
#define HeterogeneousReactingModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "CloudSubModelBase.H"
#include "scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class HeterogeneousReactingModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class HeterogeneousReactingModel
:
    public CloudSubModelBase<CloudType>
{
protected:

    // Protected data

        //- Net mass of lagrangian phase consumed
        scalar dMass_;

        //- Number of progress variables
        label nF_;


public:

    //-Runtime type information
    TypeName("heterogeneousReactionModel");


    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        HeterogeneousReactingModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& cloud
        ),
        (dict, cloud)
    );


    // Constructors

        //- Construct null from owner
        HeterogeneousReactingModel(CloudType& owner);

        //- Construct from dictionary
        HeterogeneousReactingModel
        (
            const dictionary& dict,
            CloudType& cloud,
            const word& type
        );

        //- Construct copy
        HeterogeneousReactingModel
        (
            const HeterogeneousReactingModel<CloudType>& srm
        );

        //- Construct and return a clone
        virtual autoPtr
        <
            HeterogeneousReactingModel<CloudType>
        > clone() const = 0;


    //- Destructor
    virtual ~HeterogeneousReactingModel() = default;


    //- Selector
    static autoPtr<HeterogeneousReactingModel<CloudType>> New
    (
        const dictionary& dict,
        CloudType& cloud
    );


    // Member Functions

        //- Update surface reactions
        //  Returns the heat of reaction
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
        ) const = 0;


        //- Add to devolatilisation mass
        void addToSurfaceReactionMass(const scalar dMass);

        //- Number of reactions in the model
        virtual label nReactions() const = 0;

        //- Number of progress variable
        virtual label nF() const;

        //- Write injection info
        virtual void info();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeHeterogeneousReactingModel(CloudType)                              \
                                                                               \
    typedef Foam::CloudType::reactingHeterogeneousCloudType                    \
        reactingHeterogeneousCloudType;                                        \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        Foam::HeterogeneousReactingModel<reactingHeterogeneousCloudType>,      \
        0                                                                      \
    );                                                                         \
    namespace Foam                                                             \
    {                                                                          \
        defineTemplateRunTimeSelectionTable                                    \
        (                                                                      \
            HeterogeneousReactingModel<reactingHeterogeneousCloudType>,        \
            dictionary                                                         \
        );                                                                     \
    }


#define makeHeterogeneousReactingModelType(SS, CloudType)                      \
                                                                               \
    typedef Foam::CloudType::reactingHeterogeneousCloudType                    \
        reactingHeterogeneousCloudType;                                        \
    defineNamedTemplateTypeNameAndDebug                                        \
        (Foam::SS<reactingHeterogeneousCloudType>, 0);                         \
                                                                               \
    Foam::HeterogeneousReactingModel<reactingHeterogeneousCloudType>::         \
        adddictionaryConstructorToTable                                        \
        <Foam::SS<reactingHeterogeneousCloudType>>                             \
        add##SS##CloudType##reactingHeterogeneousCloudType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "HeterogeneousReactingModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
