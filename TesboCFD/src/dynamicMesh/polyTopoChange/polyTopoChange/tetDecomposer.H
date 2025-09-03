/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
    Copyright (C) 2018-2020,2024 OpenCFD Ltd.
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
    Foam::tetDecomposer

Description
    Decomposes polyMesh into tets (or pyramids)

    Cells neighbouring decomposed cells are not decomposed themselves
    so will be polyhedral.

SourceFiles
    tetDecomposer.C

\*---------------------------------------------------------------------------*/

#ifndef tetDecomposer_H
#define tetDecomposer_H

#include "DynamicList.H"
#include "bitSet.H"
#include "boolList.H"
#include "typeInfo.H"
#include "Enum.H"
#include "faceList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyMesh;
class polyTopoChange;
class face;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                           Class tetDecomposer Declaration
\*---------------------------------------------------------------------------*/

class tetDecomposer
{
public:

    // Public data types

        enum decompositionType
        {
            FACE_CENTRE_TRIS, //- Faces decomposed into triangles
                              //  using face-centre
            FACE_DIAG_TRIS,   //- Faces decomposed into triangles diagonally
            PYRAMID,          //- Faces not decomposed (generates pyramids)
            FACE_DIAG_QUADS   //- Like FACE_DIAG_TRIS but does quads mainly
        };
        static const Enum<decompositionType> decompositionTypeNames;


private:

    // Private data

        const polyMesh& mesh_;

        //- From cell to tet point
        labelList cellToPoint_;

        //- From face to tet point
        labelList faceToPoint_;

        // Per face, per point (faceCentre) or triangle (faceDiag)
        // the added tet on the owner side
        labelListList faceOwnerCells_;

        // Per face, per point (faceCentre) or triangle (faceDiag)
        // the added tet on the neighbour side
        labelListList faceNeighbourCells_;


    // Private Member Functions

        //- Modify a face
        void modifyFace
        (
            polyTopoChange& meshMod,
            const face& f,
            const label facei,
            const label own,
            const label nei,
            const label patchi,
            const label zoneI,
            const bool zoneFlip
        ) const;

        //- Add a face
        void addFace
        (
            polyTopoChange& meshMod,
            const face& f,
            const label facei,
            const label own,
            const label nei,
            const label masterPointID,
            const label masterEdgeID,
            const label masterFaceID,
            const label patchi,
            const label zoneI,
            const bool zoneFlip
        ) const;

        //- Work out triangle index given the starting vertex in the face
        label triIndex(const label facei, const label fp) const;

        //- Calculate triangulation of boundary faces
        void splitBoundaryFaces
        (
            List<faceList>& boundaryQuads,
            List<faceList>& boundaryTris
        ) const;

        //- Correct coupled faces to match up
        void relativeIndicesToFace
        (
            const bool doFlip,
            const face& meshFace,
            const faceList& indexLists,
            faceList& faces
        ) const;

        //- Calculate triangulation of any face
        void splitFace
        (
            const List<faceList>& boundaryQuads,
            const List<faceList>& boundaryTris,
            const label facei,
            const label patchi,
            label& quadi,
            faceList& quadFaces,
            label& trii,
            faceList& triFaces
        ) const;

        void splitFace
        (
            const List<faceList>& boundaryQuads,
            const List<faceList>& boundaryTris,
            const label facei,
            const label patchi,
            faceList& quadFaces,
            faceList& triFaces,
            faceList& subFaces
        ) const;

        //- No copy construct
        tetDecomposer(const tetDecomposer&) = delete;

        //- No copy assignment
        void operator=(const tetDecomposer&) = delete;


public:

    //- Runtime type information
    ClassName("tetDecomposer");


    // Constructors

        //- Construct from mesh
        tetDecomposer(const polyMesh&);


    // Member Functions

        // Access

            //- From cell to tet point
            const labelList& cellToPoint() const
            {
                return cellToPoint_;
            }

            //- From face to tet point (only for faceCentre)
            const labelList& faceToPoint() const
            {
                return faceToPoint_;
            }


            //- Per face, per point (faceCentre) or triangle (faceDiag)
            //  the added tet on the owner side. For non-face (pyramid)
            //  size 1.
            const labelListList& faceOwnerCells() const
            {
                return faceOwnerCells_;
            }

            //- Per face, per point (faceCentre) or triangle (faceDiag)
            //  the added tet on the neighbour side. For non-face (pyramid)
            //  size 1.
            const labelListList& faceNeighbourCells() const
            {
                return faceNeighbourCells_;
            }


        // Edit

            //- Insert all changes into meshMod to convert the polyMesh into
            //  subshapes (tets/prisms)
            void setRefinement
            (
                const decompositionType decomposeType,
                const bitSet& decomposeCell,
                polyTopoChange& meshMod
            );

            //- Insert all changes into meshMod to convert the polyMesh into
            //- subshapes (tets/prisms). Explicit control over which faces
            //- get decomposed. Can be used e.g. to not split triangles.
            void setRefinement
            (
                const decompositionType decomposeType,
                const bitSet& decomposeCell,
                const bitSet& decomposeFace,
                polyTopoChange& meshMod
            );

            //- Force recalculation of locally stored data on topological change
            void updateMesh(const mapPolyMesh&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
