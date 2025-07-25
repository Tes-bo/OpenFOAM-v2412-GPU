/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Foam::fileFormats::ensightMeshReader

Description

Notes

SourceFiles
    ensightMeshReader.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_ensightMeshReader_H
#define Foam_ensightMeshReader_H

#include "meshReader.H"
//#include "ensightReadFile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class ensightReadFile;

namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                Class fileFormats::ensightMeshReader Declaration
\*---------------------------------------------------------------------------*/

class ensightMeshReader
:
    public meshReader
{
    // Private Data

        //- Merge distance
        const scalar mergeTol_;

        //- Check and correct handedness
        const bool setHandedness_;


protected:

    // Protected Data

        //- mesh point to original node_id
        labelList nodeIds_;

        //- mesh cell to original element_id
        labelList elementIds_;


    // Protected Member Functions

        //- Rotate face so lowest vertex is first
        const face& rotateFace
        (
            const face& f,
            face& rotatedFace
        ) const;

        //- Read set of vertices. Optional mapping
        void readVerts
        (
            ensightReadFile& is,
            const label nVerts,
            const Map<label>& nodeIdToPoints,
            DynamicList<label>& verts
        ) const;

        //- Read set of element/node IDs
        void readIDs
        (
            ensightReadFile& is,
            const bool doRead,
            const label nShapes,
            labelList& foamToElem,
            Map<label>& elemToFoam
        ) const;

        //- Swap handedness of hex if needed
        void setHandedness
        (
            const cellModel& model,
            DynamicList<label>& verts,
            const pointField& points
        ) const;

        //- Read a single part until eof (return true) or until start of next
        //  part (return false)
        bool readGoldPart
        (
            ensightReadFile& is,
            const bool read_node_ids,
            const bool read_elem_ids,

            pointField& points,
            labelList& pointToNodeIds,
            Map<label>& nodeIdToPoints,

            // 3D-elems : cells (cell-to-faces)
            faceListList& cellFaces,
            labelList& cellToElemIds,
            Map<label>& elemIdToCells,

            // 2D-elems : faces
            faceList& faces,
            labelList& faceToElemIDs,
            Map<label>& elemIdToFaces
        ) const;

        //- Read the mesh from the file(s)
        virtual bool readGeometry(const scalar scaleFactor = 1.0);


public:

    //- Runtime type information
    TypeName("ensightMeshReader");


    // Constructors

        //- Construct from case name
        ensightMeshReader
        (
            const fileName& geomFile,
            const objectRegistry& registry,
            const scalar mergeTol = SMALL,
            const scalar scaleFactor = 1.0,
            const bool setHandedness = true
        );


    //- Destructor
    virtual ~ensightMeshReader() = default;


    // Access

        //- Original node id (if supplied) or -1
        const labelList& nodeIds() const noexcept
        {
            return nodeIds_;
        }

        //- Original element id (if supplied) or -1
        const labelList& elementIds() const noexcept
        {
            return elementIds_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileFormats
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
