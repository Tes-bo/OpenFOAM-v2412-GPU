/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::PackingModel

Group
    grpLagrangianIntermediateMPPICPackingSubModels

Description
    Base class for packing models.

SourceFiles
    PackingModel.C
    PackingModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef PackingModel_H
#define PackingModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "CloudSubModelBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

class ParticleStressModel;

class CorrectionLimitingMethod;

template<class Type>
class AveragingMethod;

/*---------------------------------------------------------------------------*\
                         Class PackingModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class PackingModel
:
    public CloudSubModelBase<CloudType>
{
protected:

    //- Protected data

        //- Particle stress model
        autoPtr<ParticleStressModel> particleStressModel_;


public:

    //- Runtime type information
    TypeName("packingModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        PackingModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& owner
        ),
        (dict, owner)
    );


    // Constructors

        //- Construct null from owner
        PackingModel(CloudType& owner);

        //- Construct from components
        PackingModel
        (
            const dictionary& dict,
            CloudType& owner,
            const word& type
        );

        //- Construct copy
        PackingModel(const PackingModel<CloudType>& cm);

        //- Construct and return a clone
        virtual autoPtr<PackingModel<CloudType>> clone() const = 0;


    //- Destructor
    virtual ~PackingModel();


    //- Selector
    static autoPtr<PackingModel<CloudType>> New
    (
        const dictionary& dict,
        CloudType& owner
    );


    // Member Functions

        //- Calculate the velocity correction
        virtual vector velocityCorrection
        (
            typename CloudType::parcelType& p,
            const scalar deltaT
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makePackingModel(CloudType)                                            \
                                                                               \
    typedef Foam::CloudType::kinematicCloudType kinematicCloudType;            \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        Foam::PackingModel<kinematicCloudType>,                                \
        0                                                                      \
    );                                                                         \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        defineTemplateRunTimeSelectionTable                                    \
        (                                                                      \
            PackingModel<kinematicCloudType>,                                  \
            dictionary                                                         \
        );                                                                     \
    }


#define makePackingModelType(SS, CloudType)                                    \
                                                                               \
    typedef Foam::CloudType::kinematicCloudType kinematicCloudType;            \
    defineNamedTemplateTypeNameAndDebug                                        \
        (Foam::PackingModels::SS<kinematicCloudType>, 0);                      \
                                                                               \
    Foam::PackingModel<kinematicCloudType>::                                   \
        adddictionaryConstructorToTable                                        \
        <Foam::PackingModels::SS<kinematicCloudType>>                          \
            add##SS##CloudType##kinematicCloudType##ConstructorToTable_;


// These are used to defined models in MPPIC cloud specifically

#define makePackingModelMPPIC(CloudType)                                   \
                                                                               \
    typedef Foam::CloudType::MPPICCloudType MPPICCloudType;                    \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        Foam::PackingModel<MPPICCloudType>,                                    \
        0                                                                      \
    );                                                                         \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        defineTemplateRunTimeSelectionTable                                    \
        (                                                                      \
            PackingModel<MPPICCloudType>,                                      \
            dictionary                                                         \
        );                                                                     \
    }


#define makePackingModelTypeMPPIC(SS, CloudType)                               \
                                                                               \
    typedef Foam::CloudType::MPPICCloudType MPPICCloudType;                    \
    defineNamedTemplateTypeNameAndDebug                                        \
        (Foam::PackingModels::SS<MPPICCloudType>, 0);                          \
                                                                               \
    Foam::PackingModel<MPPICCloudType>::                                       \
        adddictionaryConstructorToTable                                        \
        <Foam::PackingModels::SS<MPPICCloudType>>                              \
            add##SS##CloudType##MPPICCloudType##ConstructorToTable_;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PackingModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
