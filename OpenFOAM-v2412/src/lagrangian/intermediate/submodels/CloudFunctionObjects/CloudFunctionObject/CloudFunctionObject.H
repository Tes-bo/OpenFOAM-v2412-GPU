/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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
    Foam::CloudFunctionObject

Group
    grpLagrangianIntermediateFunctionObjects

Description
    Templated cloud function object base class

SourceFiles
    CloudFunctionObject.C
    CloudFunctionObjectNew.C

\*---------------------------------------------------------------------------*/

#ifndef CloudFunctionObject_H
#define CloudFunctionObject_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "CloudSubModelBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyPatch;
class tetIndices;

/*---------------------------------------------------------------------------*\
                    Class CloudFunctionObject Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class CloudFunctionObject
:
    public CloudSubModelBase<CloudType>
{
private:

        //- Convenience typedef for parcel type
        typedef typename CloudType::parcelType parcelType;

    // Private data

        //- Output path
        fileName outputDir_;


    // Private Member Functions

        //- Inherited write from CloudSubModelBase
        using CloudSubModelBase<CloudType>::write;

        //- Write post-processing info
        virtual void write();


public:

    //- Runtime type information
    TypeName("cloudFunctionObject");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        CloudFunctionObject,
        dictionary,
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        ),
        (dict, owner, modelName)
    );


    // Constructors

        //- Construct null from owner
        CloudFunctionObject(CloudType& owner);

        //- Construct from dictionary
        CloudFunctionObject
        (
            const dictionary& dict,
            CloudType& owner,
            const word& objectType,
            const word& modelName
        );

        //- Construct copy
        CloudFunctionObject(const CloudFunctionObject<CloudType>& ppm);

        //- Construct and return a clone
        virtual autoPtr<CloudFunctionObject<CloudType>> clone() const
        {
            return autoPtr<CloudFunctionObject<CloudType>>
            (
                new CloudFunctionObject<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~CloudFunctionObject() = default;


    //- Selector
    static autoPtr<CloudFunctionObject<CloudType>> New
    (
        const dictionary& dict,
        CloudType& owner,
        const word& objectType,
        const word& modelName
    );


    // Member Functions

        // Evaluation

            //- Pre-evolve hook
            virtual void preEvolve(const typename parcelType::trackingData& td);

            //- Post-evolve hook
            virtual void postEvolve
            (
                const typename parcelType::trackingData& td
            );

            //- Post-move hook
            virtual bool postMove
            (
                parcelType& p,
                const scalar dt,
                const point& position0,
                const typename parcelType::trackingData& td
            );

            //- Post-patch hook
            virtual bool postPatch
            (
                const parcelType& p,
                const polyPatch& pp,
                const typename parcelType::trackingData& td
            );

            //- Post-face hook
            virtual bool postFace
            (
                const parcelType& p,
                const typename parcelType::trackingData& td
            );


        // Input/output

            //- Return the output path
            const fileName& outputDir() const;

            //- Return the output time path
            fileName writeTimeDir() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeCloudFunctionObject(CloudType)                                     \
                                                                               \
    typedef Foam::CloudType::kinematicCloudType kinematicCloudType;            \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        Foam::CloudFunctionObject<kinematicCloudType>,                         \
        0                                                                      \
    );                                                                         \
    namespace Foam                                                             \
    {                                                                          \
        defineTemplateRunTimeSelectionTable                                    \
        (                                                                      \
            CloudFunctionObject<kinematicCloudType>,                           \
            dictionary                                                         \
        );                                                                     \
    }


#define makeCloudFunctionObjectType(SS, CloudType)                             \
                                                                               \
    typedef Foam::CloudType::kinematicCloudType kinematicCloudType;            \
    defineNamedTemplateTypeNameAndDebug(Foam::SS<kinematicCloudType>, 0);      \
                                                                               \
    Foam::CloudFunctionObject<kinematicCloudType>::                            \
        adddictionaryConstructorToTable<Foam::SS<kinematicCloudType>>          \
            add##SS##CloudType##kinematicCloudType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CloudFunctionObject.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
