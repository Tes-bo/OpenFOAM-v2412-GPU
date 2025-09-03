/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
    Copyright (C) 2024 OpenCFD Ltd.
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
    Foam::faceMapper

Description
    This object provides mapping and fill-in information for face data
    between the two meshes after the topological change.  It is
    constructed from mapPolyMesh.

SourceFiles
    faceMapper.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_faceMapper_H
#define Foam_faceMapper_H

#include "morphFieldMapper.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                           Class faceMapper Declaration
\*---------------------------------------------------------------------------*/

class faceMapper
:
    public morphFieldMapper
{
    // Private Data

        //- Reference to mapPolyMesh
        const mapPolyMesh& mpm_;

        //- The size of the mapper = polyMesh::nFaces()
        const label mapperLen_;

        //- Number of inserted (unmapped) faces
        label nInsertedObjects_;

        //- Is the mapping direct
        bool direct_;


    // Demand-Driven Data

        //- Direct addressing (only one for of addressing is used)
        mutable std::unique_ptr<labelList> directAddrPtr_;

        //- Interpolated addressing (only one for of addressing is used)
        mutable std::unique_ptr<labelListList> interpAddrPtr_;

        //- Interpolation weights
        mutable std::unique_ptr<scalarListList> weightsPtr_;

        //- Inserted faces
        mutable std::unique_ptr<labelList> insertedObjectsPtr_;


    // Private Member Functions

        //- Calculate addressing for mapping with inserted faces
        void calcAddressing() const;

public:

    // Generated Methods

        //- No copy construct
        faceMapper(const faceMapper&) = delete;

        //- No copy assignment
        void operator=(const faceMapper&) = delete;


    // Constructors

        //- Construct from mapPolyMesh
        explicit faceMapper(const mapPolyMesh& mpm);


    //- Destructor
    virtual ~faceMapper();


    // Member Functions

        //- The mapper size
        virtual label size() const;

        //- Return size of field before mapping
        virtual label sizeBeforeMapping() const;

        //- Return number of internal faces before mapping
        virtual label internalSizeBeforeMapping() const;

        //- Is the mapping direct
        virtual bool direct() const
        {
            return direct_;
        }

        virtual bool hasUnmapped() const
        {
            return insertedObjects();
        }

        //- Return direct addressing
        virtual const labelUList& directAddressing() const;

        //- Return interpolated addressing
        virtual const labelListList& addressing() const;

        //- Return interpolaion weights
        virtual const scalarListList& weights() const;

        //- Return flux flip map
        virtual const labelHashSet& flipFaceFlux() const;

        //- Return number of old internalFaces
        virtual label nOldInternalFaces() const;

        //- Return old patch starts
        virtual const labelList& oldPatchStarts() const;

        //- Return old patch sizes
        virtual const labelList& oldPatchSizes() const;

        //- Are there any inserted faces
        virtual bool insertedObjects() const
        {
            return bool(nInsertedObjects_);
        }

        //- Return list of inserted faces
        virtual const labelList& insertedObjectLabels() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
