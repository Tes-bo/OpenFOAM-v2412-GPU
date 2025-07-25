/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
    Copyright (C) 2021-2022 OpenCFD Ltd.
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
    Foam::faFieldDecomposer

Description
    Finite Area area and edge field decomposer.

Author
    Zeljko Tukovic, FSB Zagreb
    Hrvoje Jasak, Wikki Ltd.

SourceFiles
    faFieldDecomposer.C
    fvFieldDecomposerCache.C
    faFieldDecomposerTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_faFieldDecomposer_H
#define Foam_faFieldDecomposer_H

#include "faMesh.H"
#include "faMeshSubset.H"
#include "faPatchFieldMapper.H"
#include "edgeFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class IOobjectList;

/*---------------------------------------------------------------------------*\
                    Class faFieldDecomposer Declaration
\*---------------------------------------------------------------------------*/

class faFieldDecomposer
{
public:

    // Public Classes

        //- Patch field decomposer class
        class patchFieldDecomposer
        :
            public faPatchFieldMapper
        {
            // Private Data

                label sizeBeforeMapping_;
                labelList directAddressing_;

        public:

            // Constructors

                //- Construct given addressing
                patchFieldDecomposer
                (
                    const label sizeBeforeMapping,
                    const labelUList& addressingSlice,
                    const label addressingOffset
                );


            // Member functions

                label size() const
                {
                    return directAddressing_.size();
                }

                virtual label sizeBeforeMapping() const
                {
                    return sizeBeforeMapping_;
                }

                bool direct() const
                {
                    return true;
                }

                virtual bool hasUnmapped() const
                {
                    return false;
                }

                const labelUList& directAddressing() const
                {
                    return directAddressing_;
                }
        };


        //- Processor patch field decomposer class
        class processorAreaPatchFieldDecomposer
        :
            public faPatchFieldMapper
        {
            // Private Data

                label sizeBeforeMapping_;
                labelList directAddressing_;

        public:

            //- Construct addressing from details
            processorAreaPatchFieldDecomposer
            (
                const label nTotalFaces,
                const labelUList& edgeOwner,
                const labelUList& edgeNeigbour,
                const labelUList& addressingSlice,
                const bitSet& flip
            );


            // Member Functions

                label size() const
                {
                    return directAddressing_.size();
                }

                virtual label sizeBeforeMapping() const
                {
                    return sizeBeforeMapping_;
                }

                bool direct() const
                {
                    return true;
                }

                virtual bool hasUnmapped() const
                {
                    return false;
                }

                const labelUList& directAddressing() const
                {
                    return directAddressing_;
                }
        };


        //- Processor patch field decomposer class
        class processorEdgePatchFieldDecomposer
        :
            public faPatchFieldMapper
        {
            label sizeBeforeMapping_;
            labelListList addressing_;
            scalarListList weights_;

        public:

            //- Construct given addressing
            processorEdgePatchFieldDecomposer
            (
                label sizeBeforeMapping,
                const labelUList& addressingSlice
            );


            // Member Functions

                label size() const
                {
                    return addressing_.size();
                }

                virtual label sizeBeforeMapping() const
                {
                    return sizeBeforeMapping_;
                }

                bool direct() const
                {
                    return false;
                }

                virtual bool hasUnmapped() const
                {
                    return false;
                }

                const labelListList& addressing() const
                {
                    return addressing_;
                }

                const scalarListList& weights() const
                {
                    return weights_;
                }
        };


private:

    // Private Data

        //- Reference to processor mesh
        const faMesh& procMesh_;

        //- Reference to edge addressing
        const labelList& edgeAddressing_;

        //- Reference to face addressing
        const labelList& faceAddressing_;

        //- Reference to boundary addressing
        const labelList& boundaryAddressing_;

        //- List of patch field decomposers
        PtrList<patchFieldDecomposer> patchFieldDecomposerPtrs_;

        PtrList<processorAreaPatchFieldDecomposer>
            processorAreaPatchFieldDecomposerPtrs_;

        PtrList<processorEdgePatchFieldDecomposer>
            processorEdgePatchFieldDecomposerPtrs_;


    // Private Member Functions

        //- No copy construct
        faFieldDecomposer(const faFieldDecomposer&) = delete;

        //- No copy assignment
        void operator=(const faFieldDecomposer&) = delete;


public:

    // Public Classes
    class fieldsCache;


    // Constructors

        //- Construct without mappers, added later with reset()
        faFieldDecomposer
        (
            const Foam::zero,
            const faMesh& procMesh,             // Target mesh
            const labelList& edgeAddressing,
            const labelList& faceAddressing,
            const labelList& boundaryAddressing
        );

        //- Construct from components using information from the complete mesh
        faFieldDecomposer
        (
            const faMesh& completeMesh,         // Source mesh
            const faMesh& procMesh,             // Target mesh
            const labelList& edgeAddressing,
            const labelList& faceAddressing,
            const labelList& boundaryAddressing
        );

        //- Construct from components without the complete mesh
        faFieldDecomposer
        (
            // Information about the complete mesh
            const label nTotalFaces,
            const List<labelRange>& boundaryRanges,
            const labelUList& edgeOwner,
            const labelUList& edgeNeigbour,

            // Addressing for processor mesh
            const faMesh& procMesh,             // Target mesh
            const labelList& edgeAddressing,
            const labelList& faceAddressing,
            const labelList& boundaryAddressing
        );


    //- Destructor
    ~faFieldDecomposer() = default;


    // Member Functions

        //- True if no mappers have been allocated
        bool empty() const;

        //- Remove all mappers
        void clear();

        //- Reset mappers using information from the complete mesh
        void reset(const faMesh& completeMesh);

        //- Reset mapper using information about the complete mesh
        void reset
        (
            const label nTotalFaces,
            const List<labelRange>& boundaryRanges,
            const labelUList& edgeOwner,
            const labelUList& edgeNeigbour
        );


    // Mapping

        //- Decompose area field
        template<class Type>
        tmp<GeometricField<Type, faPatchField, areaMesh>>
        decomposeField
        (
            const GeometricField<Type, faPatchField, areaMesh>& field
        ) const;

        //- Decompose surface field
        template<class Type>
        tmp<GeometricField<Type, faePatchField, edgeMesh>>
        decomposeField
        (
            const GeometricField<Type, faePatchField, edgeMesh>& field
        ) const;

        template<class GeoField>
        void decomposeFields(const PtrList<GeoField>& fields) const;


    // Reading helpers

        //- Read the fields and store on the pointer list
        template
        <
            class Type,
            template<class> class PatchField,
            class GeoMesh
        >
        static void readFields
        (
            const typename GeoMesh::Mesh& mesh,
            const IOobjectList& objects,
            PtrList<GeometricField<Type, PatchField, GeoMesh>>& fields,
            const bool readOldTime
        );

        //- Read fields and store on the pointer list
        template<class Mesh, class GeoField>
        static void readFields
        (
            const Mesh& mesh,
            const IOobjectList& objects,
            PtrList<GeoField>& fields
        );
};



/*---------------------------------------------------------------------------*\
               Class faFieldDecomposer::fieldsCache Declaration
\*---------------------------------------------------------------------------*/

class faFieldDecomposer::fieldsCache
{
    // Private Data

        class privateCache;

        //- All field and field-field types for finiteArea
        std::unique_ptr<privateCache> cache_;


    // Private Member Functions

        //- No copy construct
        fieldsCache(const fieldsCache&) = delete;

        //- No copy assignment
        void operator=(const fieldsCache&) = delete;


public:

    // Constructors

        //- Default construct
        fieldsCache();


    //- Destructor
    ~fieldsCache();


    // Member Functions

        //- No fields
        bool empty() const;

        //- Number of fields
        label size() const;

        //- Clear out
        void clear();


        //- Read all fields given mesh and objects
        void readAllFields
        (
            const faMesh& mesh,
            const IOobjectList& objects
        );

        //- Read all fields given mesh and objects.
        //- Supports reading/sending fields
        void readAllFields
        (
            const bitSet& haveMeshOnProc,
            const faMeshSubset* subsetter,
            const faMesh& mesh,
            IOobjectList& objects
        );

        //- Read all fields given mesh and objects.
        //- Supports reading/sending fields
        void readAllFields
        (
            const boolList& haveMeshOnProc,
            const faMeshSubset* subsetter,
            const faMesh& mesh,
            IOobjectList& objects
        );

        //- Decompose and write all fields
        void decomposeAllFields
        (
            const faFieldDecomposer& decomposer,
            bool report = false
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "faFieldDecomposerTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
