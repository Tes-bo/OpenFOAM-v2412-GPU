/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2022 OpenCFD Ltd.
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
    Foam::mapDistributePolyMesh

Description
    Class containing mesh-to-mesh mapping information after a mesh distribution
    where we send parts of meshes (using subsetting) to other processors
    and receive and reconstruct mesh.

    We store mapping from the bits-to-send to the complete starting mesh
    (subXXXMap) and from the received bits to their location in the new
    mesh (constructXXXMap).

SourceFiles
    mapDistributePolyMesh.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_mapDistributePolyMesh_H
#define Foam_mapDistributePolyMesh_H

#include "mapDistribute.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class mapPolyMesh;
class polyMesh;
class mapDistributePolyMesh;

Istream& operator>>(Istream&, mapDistributePolyMesh&);
Ostream& operator<<(Ostream&, const mapDistributePolyMesh&);

template<>
Ostream& operator<<(Ostream&, const InfoProxy<mapDistributePolyMesh>&);


/*---------------------------------------------------------------------------*\
                     Class mapDistributePolyMesh Declaration
\*---------------------------------------------------------------------------*/

class mapDistributePolyMesh
{
    // Private Data

        //- Number of old live points
        label nOldPoints_;

        //- Number of old live faces
        label nOldFaces_;

        //- Number of old live cells
        label nOldCells_;

        //- List of the old patch sizes
        labelList oldPatchSizes_;

        //- List of the old patch start labels
        labelList oldPatchStarts_;

        //- List of numbers of mesh points per old patch
        labelList oldPatchNMeshPoints_;


        //- Point distribute map
        mapDistribute pointMap_;

        //- Face distribute map
        mapDistribute faceMap_;

        //- Cell distribute map
        mapDistribute cellMap_;

        //- Patch distribute map
        mapDistribute patchMap_;


    // Private Member Functions

        void calcPatchSizes();

        void deepCopy(const mapDistributePolyMesh& rhs);

public:

    // Constructors

        //- Default construct - uses worldComm
        mapDistributePolyMesh();

        //- Default construct with specified communicator
        explicit mapDistributePolyMesh(const label comm);

        //- Read construct from dictionary
        explicit mapDistributePolyMesh
        (
            const dictionary& dict,
            const label comm = UPstream::worldComm
        );

        //- Copy construct
        mapDistributePolyMesh(const mapDistributePolyMesh& map);

        //- Move construct
        mapDistributePolyMesh(mapDistributePolyMesh&& map);

        //- Copy/move construct
        mapDistributePolyMesh(mapDistributePolyMesh& map, bool reuse);

        //- Construct from components. Note that mesh has to be changed already
        //  since uses mesh.nPoints etc as the new size.
        mapDistributePolyMesh
        (
            const polyMesh& mesh,

            // mesh before changes
            const label nOldPoints,
            const label nOldFaces,
            const label nOldCells,
            labelList&& oldPatchStarts,
            labelList&& oldPatchNMeshPoints,

            // how to subset pieces of mesh to send across
            labelListList&& subPointMap,
            labelListList&& subFaceMap,
            labelListList&& subCellMap,
            labelListList&& subPatchMap,

            // how to reconstruct received mesh
            labelListList&& constructPointMap,
            labelListList&& constructFaceMap,
            labelListList&& constructCellMap,
            labelListList&& constructPatchMap,

            const bool subFaceHasFlip = false,
            const bool constructFaceHasFlip = false
        );

        //- Move construct from components
        mapDistributePolyMesh
        (
            // mesh before changes
            const label nOldPoints,
            const label nOldFaces,
            const label nOldCells,
            labelList&& oldPatchStarts,
            labelList&& oldPatchNMeshPoints,

            // how to subset pieces of mesh to send across
            mapDistribute&& pointMap,
            mapDistribute&& faceMap,
            mapDistribute&& cellMap,
            mapDistribute&& patchMap
        );

        //- Construct from Istream
        explicit mapDistributePolyMesh(Istream& is);


    // Member Functions

        // Access

            //- Number of points in mesh before distribution
            label nOldPoints() const noexcept
            {
                return nOldPoints_;
            }

            //- Number of faces in mesh before distribution
            label nOldFaces() const noexcept
            {
                return nOldFaces_;
            }

            //- Number of cells in mesh before distribution
            label nOldCells() const noexcept
            {
                return nOldCells_;
            }

            //- List of the old patch sizes
            const labelList& oldPatchSizes() const noexcept
            {
                return oldPatchSizes_;
            }

            //- List of the old patch start labels
            const labelList& oldPatchStarts() const noexcept
            {
                return oldPatchStarts_;
            }

            //- List of numbers of mesh points per old patch
            const labelList& oldPatchNMeshPoints() const noexcept
            {
                return oldPatchNMeshPoints_;
            }

            //- Point distribute map
            const mapDistribute& pointMap() const noexcept
            {
                return pointMap_;
            }

            //- Face distribute map
            const mapDistribute& faceMap() const noexcept
            {
                return faceMap_;
            }

            //- Cell distribute map
            const mapDistribute& cellMap() const noexcept
            {
                return cellMap_;
            }

            //- Patch distribute map
            const mapDistribute& patchMap() const noexcept
            {
                return patchMap_;
            }


        // Edit

            //- Reset to zero size, only retaining communicator(s)
            void clear();


        // Other

            //- Transfer the contents of the argument and annul the argument.
            void transfer(mapDistributePolyMesh& map);

            //- Distribute list of point data
            template<class T>
            void distributePointData(List<T>& values) const
            {
                pointMap_.distribute(values);
            }

            //- Distribute list of face data
            template<class T>
            void distributeFaceData(List<T>& values) const
            {
                faceMap_.distribute(values);
            }

            //- Distribute list of cell data
            template<class T>
            void distributeCellData(List<T>& values) const
            {
                cellMap_.distribute(values);
            }

            //- Distribute list of patch data
            template<class T>
            void distributePatchData(List<T>& values) const
            {
                patchMap_.distribute(values);
            }


            //- Distribute list of point/face/cell/patch indices.
            //  (Converts to boolList, distributes boolList and reconstructs)
            void distributePointIndices(labelList& pointIDs) const;

            void distributeFaceIndices(labelList& faceIDs) const;
            void distributeCellIndices(labelList& cellIDs) const;
            void distributePatchIndices(labelList& patchIDs) const;



    // Member Operators

        //- Copy assignment
        void operator=(const mapDistributePolyMesh& map);

        //- Move assignment
        void operator=(mapDistributePolyMesh&& map);


    // IOstream Operators

        //- Return info proxy,
        //- used to print summary information to a stream
        InfoProxy<mapDistributePolyMesh> info() const noexcept
        {
            return *this;
        }

        //- Read entries from dictionary format
        void readDict(const dictionary& dict);


        //- Write cellMap in dictionary format
        void writeCellMapEntries(Ostream& os) const;

        //- Write faceMap in dictionary format
        void writeFaceMapEntries(Ostream& os) const;

        //- Write pointMap entries in dictionary format
        void writePointMapEntries(Ostream& os) const;

        //- Write patchMap in dictionary format
        void writePatchMapEntries(Ostream& os) const;

        //- Write all map entries in dictionary format
        void writeEntries(Ostream& os) const;


        //- Read plain content (not dictionary) from Istream
        friend Istream& operator>>(Istream&, mapDistributePolyMesh&);

        //- Write plain content (not dictionary) to Ostream
        friend Ostream& operator<<(Ostream&, const mapDistributePolyMesh&);


    // Housekeeping

        //- No correction for topo change
        void updateMesh(const mapPolyMesh&)
        {
            NotImplemented;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
