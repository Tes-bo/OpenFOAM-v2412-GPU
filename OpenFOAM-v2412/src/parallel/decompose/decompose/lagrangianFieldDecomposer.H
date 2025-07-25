/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2022 OpenCFD Ltd.
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
    Foam::lagrangianFieldDecomposer

Description
    Lagrangian field decomposer.

SourceFiles
    lagrangianFieldDecomposer.C
    lagrangianFieldDecomposerCache.C
    lagrangianFieldDecomposerReadFields.C
    lagrangianFieldDecomposerTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_lagrangianFieldDecomposer_H
#define Foam_lagrangianFieldDecomposer_H

#include "Cloud.H"
#include "CompactIOField.H"
#include "indexedParticle.H"
#include "passiveParticle.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class IOobjectList;

/*---------------------------------------------------------------------------*\
                    Class lagrangianFieldDecomposer Declaration
\*---------------------------------------------------------------------------*/

class lagrangianFieldDecomposer
{
    // Private Data

        //- Reference to processor mesh
        const polyMesh& procMesh_;

        //- Lagrangian positions for this processor
        Cloud<passiveParticle> positions_;

        //- The indices of the particles on this processor
        labelList particleIndices_;


    // Private Member Functions

        //- No copy construct
        lagrangianFieldDecomposer(const lagrangianFieldDecomposer&) = delete;

        //- No copy assignment
        void operator=(const lagrangianFieldDecomposer&) = delete;


public:

    // Public Classes
    class fieldsCache;


    // Constructors

        //- Construct from components
        lagrangianFieldDecomposer
        (
            const polyMesh& mesh,       //<! unused
            const polyMesh& procMesh,
            const labelList& faceProcAddressing,
            const labelList& cellProcAddressing,
            const word& cloudName,
            const Cloud<indexedParticle>& lagrangianPositions,
            const List<SLList<indexedParticle*>*>& cellParticles
        );


    // Field Reading

        //- Read the fields and store on the pointer list
        template<class Type>
        static void readFields
        (
            const label cloudi,
            const IOobjectList& lagrangianObjects,
            PtrList<PtrList<IOField<Type>>>& cloudFields
        );

        //- Read the field-fields and store on the pointer list
        template<class Type>
        static void readFieldFields
        (
            const label cloudi,
            const IOobjectList& lagrangianObjects,
            PtrList<PtrList<CompactIOField<Field<Type>, Type>>>& cloudFields
        );


    // Member Functions

        //- Decompose volume field
        template<class Type>
        tmp<IOField<Type>> decomposeField
        (
            const word& cloudName,
            const IOField<Type>& field
        ) const;

        template<class Type>
        tmp<CompactIOField<Field<Type>, Type>> decomposeFieldField
        (
            const word& cloudName,
            const CompactIOField<Field<Type>, Type>& field
        ) const;


        template<class GeoField>
        void decomposeFields
        (
            const word& cloudName,
            const PtrList<GeoField>& fields
        ) const;

        template<class GeoField>
        void decomposeFieldFields
        (
            const word& cloudName,
            const PtrList<GeoField>& fields
        ) const;
};


/*---------------------------------------------------------------------------*\
           Class lagrangianFieldDecomposer::fieldsCache Declaration
\*---------------------------------------------------------------------------*/

class lagrangianFieldDecomposer::fieldsCache
{
    // Private Data

        class privateCache;

        //- All field and field-field types for lagrangian
        std::unique_ptr<privateCache> cache_;


    // Private Member Functions

        //- No copy construct
        fieldsCache(const fieldsCache&) = delete;

        //- No copy assignment
        void operator=(const fieldsCache&) = delete;


public:

    // Constructors

        //- Default construct (no clouds)
        fieldsCache();

        //- Construct for given number of clouds
        explicit fieldsCache(const label nClouds);


    //- Destructor
    ~fieldsCache();


    // Member Functions

        //- No clouds
        bool empty() const;

        //- Number of clouds
        label size() const;

        //- Clear out
        void clear();

        //- Resize for the number of clouds
        void resize(const label nClouds);


        //- Read all fields and field-fields for given cloud and objects
        void readAllFields
        (
            const label cloudi,
            const IOobjectList& lagrangianObjects
        );

        //- Decompose and write all fields and field-fields for given cloud
        void decomposeAllFields
        (
            const label cloudi,
            const fileName& cloudDir,
            const lagrangianFieldDecomposer& decomposer,
            bool report = false
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "lagrangianFieldDecomposerReadFields.C"
    #include "lagrangianFieldDecomposerTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
