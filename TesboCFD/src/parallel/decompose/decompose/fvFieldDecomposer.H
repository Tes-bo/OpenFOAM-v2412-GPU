/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::fvFieldDecomposer

Description
    Finite Volume volume and surface field decomposer.

SourceFiles
    fvFieldDecomposer.C
    fvFieldDecomposerCache.C
    fvFieldDecomposerTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_fvFieldDecomposer_H
#define Foam_fvFieldDecomposer_H

#include "fvMesh.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class IOobjectList;

/*---------------------------------------------------------------------------*\
                    Class fvFieldDecomposer Declaration
\*---------------------------------------------------------------------------*/

class fvFieldDecomposer
{
public:

    // Public Classes

        //- Patch field decomposer class
        class patchFieldDecomposer
        :
            public fvPatchFieldMapper
        {
            // Private Data

                labelList directAddressing_;

        public:

            // Constructors

                //- Construct given addressing
                patchFieldDecomposer
                (
                    const labelUList& addressingSlice,
                    const label addressingOffset
                );


            // Member Functions

                label size() const
                {
                    return directAddressing_.size();
                }

                bool direct() const
                {
                    return true;
                }

                //- Are there unmapped values
                bool hasUnmapped() const
                {
                    return false;
                }

                const labelUList& directAddressing() const
                {
                    return directAddressing_;
                }
        };


        //- Processor patch field decomposer class. Maps either owner or
        //  neighbour data (no interpolate anymore - processorFvPatchField
        //  holds neighbour data)
        class processorVolPatchFieldDecomposer
        :
            public fvPatchFieldMapper
        {
            // Private Data

                labelList directAddressing_;

        public:

            //- Construct addressing from details
            processorVolPatchFieldDecomposer
            (
                const labelUList& faceOwner,
                const labelUList& faceNeigbour,
                const labelUList& addressingSlice
            );

            //- Construct given addressing from complete mesh
            processorVolPatchFieldDecomposer
            (
                const polyMesh& mesh,
                const labelUList& addressingSlice
            );


            // Member Functions

                label size() const
                {
                    return directAddressing_.size();
                }

                bool direct() const
                {
                    return true;
                }

                //- Are there unmapped values
                bool hasUnmapped() const
                {
                    return false;
                }

                const labelUList& directAddressing() const
                {
                    return directAddressing_;
                }
        };


        //- Processor patch field decomposer class. Surface field is assumed
        //  to have direction (so manipulates sign when mapping)
        class processorSurfacePatchFieldDecomposer
        :
            public fvPatchFieldMapper
        {
            labelListList addressing_;
            scalarListList weights_;

        public:

            //- Construct given addressing
            processorSurfacePatchFieldDecomposer
            (
                const labelUList& addressingSlice
            );


            // Member functions

                label size() const
                {
                    return addressing_.size();
                }

                bool direct() const
                {
                    return false;
                }

                //- Are there unmapped values
                bool hasUnmapped() const
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
        const fvMesh& procMesh_;

        //- Reference to face addressing
        const labelList& faceAddressing_;

        //- Reference to cell addressing
        const labelList& cellAddressing_;

        //- Reference to boundary addressing
        const labelList& boundaryAddressing_;

        //- List of patch field decomposers
        PtrList<patchFieldDecomposer> patchFieldDecomposerPtrs_;

        PtrList<processorVolPatchFieldDecomposer>
            processorVolPatchFieldDecomposerPtrs_;

        PtrList<processorSurfacePatchFieldDecomposer>
            processorSurfacePatchFieldDecomposerPtrs_;

        PtrList<scalarField> faceSign_;


    // Private Member Functions

        //- No copy construct
        fvFieldDecomposer(const fvFieldDecomposer&) = delete;

        //- No copy assignment
        void operator=(const fvFieldDecomposer&) = delete;

public:

    // Public Classes
    class fieldsCache;


    // Static Data

        //- Output verbosity when writing
        static int verbose_;


    // Constructors

        //- Construct without mappers, added later with reset()
        fvFieldDecomposer
        (
            const Foam::zero,
            const fvMesh& procMesh,
            const labelList& faceAddressing,
            const labelList& cellAddressing,
            const labelList& boundaryAddressing
        );

        //- Construct from components using information from the complete mesh
        fvFieldDecomposer
        (
            const fvMesh& completeMesh,
            const fvMesh& procMesh,
            const labelList& faceAddressing,
            const labelList& cellAddressing,
            const labelList& boundaryAddressing
        );

        //- Construct from components without the complete mesh
        fvFieldDecomposer
        (
            // Information about the complete mesh
            const List<labelRange>& boundaryRanges,
            const labelUList& faceOwner,
            const labelUList& faceNeigbour,

            // Addressing for processor mesh
            const fvMesh& procMesh,
            const labelList& faceAddressing,
            const labelList& cellAddressing,
            const labelList& boundaryAddressing
        );


    //- Destructor
    ~fvFieldDecomposer() = default;


    // Member Functions

        //- True if no mappers have been allocated
        bool empty() const;

        //- Remove all mappers
        void clear();

        //- Reset mappers using information from the complete mesh
        void reset(const fvMesh& completeMesh);

        //- Reset mapper using information about the complete mesh
        void reset
        (
            const List<labelRange>& boundaryRanges,
            const labelUList& faceOwner,
            const labelUList& faceNeigbour
        );


    // Mapping

        //- Decompose internal field
        template<class Type>
        tmp<DimensionedField<Type, volMesh>>
        decomposeField
        (
            const DimensionedField<Type, volMesh>& field
        ) const;

        //- Decompose volume field
        template<class Type>
        tmp<GeometricField<Type, fvPatchField, volMesh>>
        decomposeField
        (
            const GeometricField<Type, fvPatchField, volMesh>& field,
            const bool allowUnknownPatchFields = false
        ) const;

        //- Decompose surface field
        template<class Type>
        tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        decomposeField
        (
            const GeometricField<Type, fvsPatchField, surfaceMesh>& field
        ) const;

        //- Decompose list of fields
        template<class GeoField>
        void decomposeFields(const PtrList<GeoField>& fields) const;
};


/*---------------------------------------------------------------------------*\
               Class fvFieldDecomposer::fieldsCache Declaration
\*---------------------------------------------------------------------------*/

class fvFieldDecomposer::fieldsCache
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

        //- Default construct
        fieldsCache();


    //- Destructor
    ~fieldsCache();


    // Member Functions

        //- No fields
        bool empty() const;

        //- Total number of fields
        label size() const;

        //- Clear out
        void clear();


        //- Read all fields given mesh and objects
        void readAllFields
        (
            const fvMesh& mesh,
            const IOobjectList& objects
        );

        //- Decompose and write all fields
        void decomposeAllFields
        (
            const fvFieldDecomposer& decomposer,
            bool report = false
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fvFieldDecomposerTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
