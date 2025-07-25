/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::PatchInteractionModel

Group
    grpLagrangianIntermediatePatchInteractionSubModels

Description
    Templated patch interaction model class

SourceFiles
    PatchInteractionModel.C
    PatchInteractionModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef PatchInteractionModel_H
#define PatchInteractionModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "polyPatch.H"
#include "wallPolyPatch.H"
#include "tetIndices.H"
#include "CloudSubModelBase.H"
#include "writeFile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class PatchInteractionModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class PatchInteractionModel
:
    public CloudSubModelBase<CloudType>,
    public functionObjects::writeFile
{
public:

    // Public enumerations

        // Interaction types
        enum interactionType
        {
            itNone,
            itRebound,
            itStick,
            itEscape,
            itOther
        };

        static wordList interactionTypeNames_;


protected:

    // Protected data

        //- Name of velocity field - default = "U"
        const word UName_;


        // Counters

            //- Number of parcels escaped
            label escapedParcels_;

            //- Mass of parcels escaped
            scalar escapedMass_;

            //- Maximum relative U with patch for particle to be removed
            scalar Urmax_;


    // Protected Member Functions

        //- Output file header information
        virtual void writeFileHeader(Ostream& os);


public:

    //- Runtime type information
    TypeName("patchInteractionModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        PatchInteractionModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& owner
        ),
        (dict, owner)
    );


    // Constructors

        //- Construct null from owner
        PatchInteractionModel(CloudType& owner);

        //- Construct from components
        PatchInteractionModel
        (
            const dictionary& dict,
            CloudType& owner,
            const word& type
        );

        //- Construct copy
        PatchInteractionModel(const PatchInteractionModel<CloudType>& pim);

        //- Construct and return a clone
        virtual autoPtr<PatchInteractionModel<CloudType>> clone() const = 0;


    //- Destructor
    virtual ~PatchInteractionModel() = default;


    //- Selector
    static autoPtr<PatchInteractionModel<CloudType>> New
    (
        const dictionary& dict,
        CloudType& owner
    );


    // Access

        //- Return name of velocity field
        const word& UName() const;

        //- Return Urmax
        const scalar& Urmax() const;

    // Member Functions

        //- Convert interaction result to word
        static word interactionTypeToWord(const interactionType& itEnum);

        //- Convert word to interaction result
        static interactionType wordToInteractionType(const word& itWord);

        //- Apply velocity correction
        //  Returns true if particle remains in same cell
        virtual bool correct
        (
            typename CloudType::parcelType& p,
            const polyPatch& pp,
            bool& keepParticle
        ) = 0;

        //- Add to escaped parcels
        virtual void addToEscapedParcels(const scalar mass);

        //- Post-evolve hook
        virtual void postEvolve();

        //- Write patch interaction info
        virtual void info();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makePatchInteractionModel(CloudType)                                   \
                                                                               \
    typedef Foam::CloudType::kinematicCloudType kinematicCloudType;            \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        Foam::PatchInteractionModel<kinematicCloudType>,                       \
        0                                                                      \
    );                                                                         \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        defineTemplateRunTimeSelectionTable                                    \
        (                                                                      \
            PatchInteractionModel<kinematicCloudType>,                         \
            dictionary                                                         \
        );                                                                     \
    }


#define makePatchInteractionModelType(SS, CloudType)                           \
                                                                               \
    typedef Foam::CloudType::kinematicCloudType kinematicCloudType;            \
    defineNamedTemplateTypeNameAndDebug(Foam::SS<kinematicCloudType>, 0);      \
                                                                               \
    Foam::PatchInteractionModel<kinematicCloudType>::                          \
        adddictionaryConstructorToTable<Foam::SS<kinematicCloudType>>          \
            add##SS##CloudType##kinematicCloudType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PatchInteractionModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
