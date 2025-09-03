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
    Foam::pointFieldDecomposer

Description
    Point field decomposer.

SourceFiles
    pointFieldDecomposer.C
    pointFieldDecomposerCache.C
    pointFieldDecomposerTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_pointFieldDecomposer_H
#define Foam_pointFieldDecomposer_H

#include "pointMesh.H"
#include "pointPatchFieldMapperPatchRef.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class IOobjectList;

/*---------------------------------------------------------------------------*\
                    Class pointFieldDecomposer Declaration
\*---------------------------------------------------------------------------*/

class pointFieldDecomposer
{
public:

    // Public Classes

        //- Point patch field decomposer class
        class patchFieldDecomposer
        :
            public pointPatchFieldMapperPatchRef
        {
            // Private data

                labelList directAddressing_;

                //- Does map contain any unmapped values
                bool hasUnmapped_;

        public:

            // Constructors

                //- Construct given addressing
                patchFieldDecomposer
                (
                    const pointPatch& completeMeshPatch,
                    const pointPatch& procMeshPatch,
                    const labelList& directAddr
                );


            // Member functions

                label size() const
                {
                    return directAddressing_.size();
                }

                bool direct() const
                {
                    return true;
                }

                bool hasUnmapped() const
                {
                    return hasUnmapped_;
                }

                const labelUList& directAddressing() const
                {
                    return directAddressing_;
                }
        };


private:

    // Private Data

        //- Reference to processor mesh
        const pointMesh& procMesh_;

        //- Reference to point addressing
        const labelList& pointAddressing_;

        //- Reference to boundary addressing
        const labelList& boundaryAddressing_;

        //- List of patch field decomposers
        PtrList<patchFieldDecomposer> patchFieldDecomposerPtrs_;


    // Private Member Functions

        //- No copy construct
        pointFieldDecomposer(const pointFieldDecomposer&) = delete;

        //- No copy assignment
        void operator=(const pointFieldDecomposer&) = delete;

public:

    // Public Classes
    class fieldsCache;


    // Constructors

        //- Construct without mappers, added later with reset()
        pointFieldDecomposer
        (
            const Foam::zero,
            const pointMesh& procMesh,
            const labelList& pointAddressing,
            const labelList& boundaryAddressing
        );

        //- Construct from components
        pointFieldDecomposer
        (
            const pointMesh& completeMesh,
            const pointMesh& procMesh,
            const labelList& pointAddressing,
            const labelList& boundaryAddressing
        );


    //- Destructor
    ~pointFieldDecomposer() = default;


    // Member Functions

        //- True if no mappers have been allocated
        bool empty() const;

        //- Remove all mappers
        void clear();

        //- Reset mappers using information from the complete mesh
        void reset(const pointMesh& completeMesh);


    // Mapping

        //- Decompose point field
        template<class Type>
        tmp<GeometricField<Type, pointPatchField, pointMesh>>
        decomposeField
        (
            const GeometricField<Type, pointPatchField, pointMesh>&
        ) const;

        //- Decompose list of fields
        template<class GeoField>
        void decomposeFields(const PtrList<GeoField>& fields) const;
};


/*---------------------------------------------------------------------------*\
              Class pointFieldDecomposer::fieldsCache Declaration
\*---------------------------------------------------------------------------*/

class pointFieldDecomposer::fieldsCache
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
            const pointMesh& mesh,
            const IOobjectList& objects
        );

        //- Decompose and write all fields
        void decomposeAllFields
        (
            const pointFieldDecomposer& decomposer,
            bool report = false
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "pointFieldDecomposerTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
