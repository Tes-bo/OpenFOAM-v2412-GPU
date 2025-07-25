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
    Foam::treeDataCell

Description
    Encapsulation of data needed to search in/for cells. Used to find the
    cell containing a point (e.g. cell-cell mapping).

SourceFiles
    treeDataCell.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_treeDataCell_H
#define Foam_treeDataCell_H

#include "polyMesh.H"
#include "treeBoundBoxList.H"
#include "volumeType.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class Type> class indexedOctree;

/*---------------------------------------------------------------------------*\
                        Class treeDataCell Declaration
\*---------------------------------------------------------------------------*/

class treeDataCell
{
    // Private Data

        //- Reference to the mesh
        const polyMesh& mesh_;

        //- Subset of cells to work on
        const labelList cellLabels_;

        //- Use subset of cells (cellLabels)
        const bool useSubset_;

        //- Whether to precalculate and store cell bounding box
        const bool cacheBb_;

        //- How to decide if point is inside cell
        const polyMesh::cellDecomposition decompMode_;

        //- Cell bounding boxes (valid only if cacheBb_)
        treeBoundBoxList bbs_;


    // Private Member Functions

        //- Get cell bounding box at specified index
        inline treeBoundBox getBounds(const label index) const;

        //- Initialise all member data
        void update();

public:


    class findNearestOp
    {
        const indexedOctree<treeDataCell>& tree_;

    public:

        findNearestOp(const indexedOctree<treeDataCell>& tree);

        void operator()
        (
            const labelUList& indices,
            const point& sample,

            scalar& nearestDistSqr,
            label& minIndex,
            point& nearestPoint
        ) const;

        void operator()
        (
            const labelUList& indices,
            const linePointRef& ln,

            treeBoundBox& tightest,
            label& minIndex,
            point& linePoint,
            point& nearestPoint
        ) const;
    };


    class findIntersectOp
    {
        const indexedOctree<treeDataCell>& tree_;

    public:

        findIntersectOp(const indexedOctree<treeDataCell>& tree);

        bool operator()
        (
            const label index,
            const point& start,
            const point& end,
            point& intersectionPoint
        ) const;
    };


    // Declare name of the class
    ClassNameNoDebug("treeDataCell");


    // Constructors (cachable versions)

        //- Construct from mesh, using all cells in mesh.
        treeDataCell
        (
            const bool cacheBb,
            const polyMesh& mesh,
            polyMesh::cellDecomposition decompMode
        );

        //- Construct from mesh, copying subset of cells.
        treeDataCell
        (
            const bool cacheBb,
            const polyMesh& mesh,
            const labelUList& cellLabels,
            polyMesh::cellDecomposition decompMode
        );

        //- Construct from mesh, moving subset of cells
        treeDataCell
        (
            const bool cacheBb,
            const polyMesh& mesh,
            labelList&& cellLabels,
            polyMesh::cellDecomposition decompMode
        );


    // Constructors (non-caching)

        //- Construct from mesh, using all cells in mesh.
        treeDataCell
        (
            const polyMesh& mesh,
            polyMesh::cellDecomposition decompMode
        )
        :
            treeDataCell(false, mesh, decompMode)
        {}

        //- Construct from mesh, copying subset of cells.
        treeDataCell
        (
            const polyMesh& mesh,
            const labelUList& cellLabels,
            polyMesh::cellDecomposition decompMode
        )
        :
            treeDataCell(false, mesh, cellLabels, decompMode)
        {}

        //- Construct from mesh, moving subset of cells
        treeDataCell
        (
            const polyMesh& mesh,
            labelList&& cellLabels,
            polyMesh::cellDecomposition decompMode
        )
        :
            treeDataCell(false, mesh, std::move(cellLabels), decompMode)
        {}


    // Static Functions

        //- Calculate and return bounding boxes for all mesh cells
        static treeBoundBoxList boxes
        (
            const primitiveMesh& mesh
        );

        //- Calculate and return bounding boxes for specified mesh cells
        static treeBoundBoxList boxes
        (
            const primitiveMesh& mesh,
            const labelUList& cellIds
        );

        //- Return bounding box of specified mesh cells
        static treeBoundBox bounds
        (
            const primitiveMesh& mesh,
            const labelUList& cellIds
        );


    // Member Functions

        //- Object dimension == 3 (volume element)
        int nDim() const noexcept { return 3; }

        //- Return bounding box for the specified cell indices
        treeBoundBox bounds(const labelUList& indices) const;


    // Access

        //- Reference to the supporting mesh
        const polyMesh& mesh() const noexcept { return mesh_; }

        //- The cell decomposition mode used
        polyMesh::cellDecomposition decompMode() const noexcept
        {
            return decompMode_;
        }

        //- The subset of cell ids to use
        const labelList& cellLabels() const noexcept { return cellLabels_; }

        //- Use a subset of cells
        bool useSubset() const noexcept { return useSubset_; }

        //- Is the effective cell selection empty?
        bool empty() const noexcept
        {
            return useSubset_ ? cellLabels_.empty() : !mesh_.nCells();
        }

        //- The size of the cell selection
        label size() const noexcept
        {
            return useSubset_ ? cellLabels_.size() : mesh_.nCells();
        }

        //- Map from shape index to original (non-subset) cell label
        label objectIndex(const label index) const
        {
            return useSubset_ && index >= 0 ? cellLabels_[index] : index;
        }

        //TBD //- Cell at specified shape index
        ///const cell& operator[](const label index) const
        ///{ return mesh_.cells()[objectIndex(index)]; }

        //- Representative point (cell centre) at shape index
        const point& centre(const label index) const
        {
            return mesh_.cellCentres()[objectIndex(index)];
        }

        //- Representative point cloud for contained shapes.
        //- One point per shape, corresponding to the cell centres.
        tmp<pointField> centres() const;


    // Search

        //- Get type (inside,outside,mixed,unknown) of point w.r.t. surface.
        //  Only makes sense for closed surfaces.
        volumeType getVolumeType
        (
            const indexedOctree<treeDataCell>&,
            const point&
        ) const
        {
            NotImplemented;
            return volumeType::UNKNOWN;
        }

        //- Does (bb of) shape at index overlap searchBox
        bool overlaps
        (
            const label index,
            const treeBoundBox& searchBox
        ) const;

        //- Does shape at index contain sample
        bool contains
        (
            const label index,
            const point& sample
        ) const;

        //- Calculates nearest (to sample) point in shape.
        //  Returns actual point and distance (squared)
        void findNearest
        (
            const labelUList& indices,
            const point& sample,

            scalar& nearestDistSqr,
            label& nearestIndex,
            point& nearestPoint
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
