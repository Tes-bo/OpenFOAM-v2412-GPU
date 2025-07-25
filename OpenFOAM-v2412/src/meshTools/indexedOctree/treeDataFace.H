/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2022 OpenCFD Ltd.
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
    Foam::treeDataFace

Description
    Encapsulation of data for searching on faces.

SourceFiles
    treeDataFace.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_treeDataFace_H
#define Foam_treeDataFace_H

#include "face.H"
#include "indexedOctree.H"
#include "treeBoundBoxList.H"
#include "bitSet.H"
#include "primitiveMesh.H"
#include "volumeType.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class polyPatch;

/*---------------------------------------------------------------------------*\
                        Class treeDataFace Declaration
\*---------------------------------------------------------------------------*/

class treeDataFace
{
    // Static Data

        //- Tolerance on linear dimensions
        static scalar tolSqr;


    // Private Data

        //- Reference to the mesh
        const primitiveMesh& mesh_;

        //- Subset of faces to work on
        const labelList faceLabels_;

        //- Inverse of faceLabels. For every mesh whether face is in faceLabels.
        bitSet isTreeFace_;

        //- Use subset of faces (faceLabels)
        const bool useSubset_;

        //- Whether to precalculate and store face bounding box
        const bool cacheBb_;

        //- Face bounding boxes (valid only if cacheBb_)
        treeBoundBoxList bbs_;


    // Private Member Functions

        //- True if specified mesh facei is being used
        inline bool usesFace(const label facei) const;

        //- Get face bounding box at specified index
        inline treeBoundBox getBounds(const label index) const;

        //- Initialise all member data
        void update();

public:


    class findNearestOp
    {
        const indexedOctree<treeDataFace>& tree_;

    public:

        findNearestOp(const indexedOctree<treeDataFace>& tree);

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
        const indexedOctree<treeDataFace>& tree_;

    public:

        findIntersectOp(const indexedOctree<treeDataFace>& tree);

        //- Calculate intersection of triangle with ray.
        //  Sets result accordingly
        bool operator()
        (
            const label index,
            const point& start,
            const point& end,
            point& intersectionPoint
        ) const;
    };


    // Declare name of the class and its debug switch
    ClassName("treeDataFace");


    // Constructors (cachable versions)

        //- Construct from mesh, using all faces
        treeDataFace
        (
            const bool cacheBb,
            const primitiveMesh& mesh
        );

        //- Construct from mesh, using polyPatch faces
        treeDataFace(const bool cacheBb, const polyPatch& patch);

        //- Construct from mesh, using specified range of faces
        treeDataFace
        (
            const bool cacheBb,
            const primitiveMesh& mesh,
            const labelRange& range
        );

        //- Construct from mesh, copying subset of faces
        treeDataFace
        (
            const bool cacheBb,
            const primitiveMesh& mesh,
            const labelUList& faceLabels
        );

        //- Construct from mesh, moving subset of faces
        treeDataFace
        (
            const bool cacheBb,
            const primitiveMesh& mesh,
            labelList&& faceLabels
        );


    // Constructors (non-caching)

        //- Construct from mesh, using all faces
        explicit treeDataFace(const primitiveMesh& mesh)
        :
            treeDataFace(false, mesh)
        {}

        //- Construct from mesh, using polyPatch faces
        explicit treeDataFace(const polyPatch& patch)
        :
            treeDataFace(false, patch)
        {}

        //- Construct from mesh, using specified range of faces
        treeDataFace(const primitiveMesh& mesh, const labelRange& range)
        :
            treeDataFace(false, mesh, range)
        {}

        //- Construct from mesh, copying subset of faces
        treeDataFace(const primitiveMesh& mesh, const labelUList& faceLabels)
        :
            treeDataFace(false, mesh, faceLabels)
        {}

        //- Construct from mesh, moving subset of faces
        treeDataFace(const primitiveMesh& mesh, labelList&& faceLabels)
        :
            treeDataFace(false, mesh, std::move(faceLabels))
        {}


    // Static Functions

        //- Calculate and return bounding boxes for all mesh faces
        static treeBoundBoxList boxes
        (
            const primitiveMesh& mesh
        );

        //- Calculate and return bounding boxes for specified range of faces
        static treeBoundBoxList boxes
        (
            const primitiveMesh& mesh,
            const labelRange& range
        );

        //- Calculate and return bounding boxes for specified mesh faces
        static treeBoundBoxList boxes
        (
            const primitiveMesh& mesh,
            const labelUList& faceIds
        );

        //- Return bounding box of specified range of faces
        static treeBoundBox bounds
        (
            const primitiveMesh& mesh,
            const labelRange& range
        );

        //- Return bounding box of specified mesh faces
        static treeBoundBox bounds
        (
            const primitiveMesh& mesh,
            const labelUList& faceIds
        );


    // Member Functions

        //- Object dimension == 2 (face element)
        int nDim() const noexcept { return 2; }

        //- Return bounding box for the specified face indices
        treeBoundBox bounds(const labelUList& indices) const;


    // Access

        //- Reference to the supporting mesh
        const primitiveMesh& mesh() const noexcept { return mesh_; }

        //- The subset of face ids to use
        const labelList& faceLabels() const noexcept { return faceLabels_; }

        //- Use a subset of faces
        bool useSubset() const noexcept { return useSubset_; }

        //- Is the effective face selection empty?
        bool empty() const noexcept
        {
            return useSubset_ ? faceLabels_.empty() : !mesh_.nFaces();
        }

        //- The size of the face selection
        label size() const noexcept
        {
            return useSubset_ ? faceLabels_.size() : mesh_.nFaces();
        }

        //- Map from shape index to original (non-subset) face label
        label objectIndex(const label index) const
        {
            return useSubset_ && index >= 0 ? faceLabels_[index] : index;
        }

        //- Face at specified shape index
        const face& operator[](const label index) const
        {
            return mesh_.faces()[objectIndex(index)];
        }

        //- Representative point (face centre) at shape index
        const point& centre(const label index) const
        {
            return mesh_.faceCentres()[objectIndex(index)];
        }

        //- Representative point cloud for contained shapes.
        //- One point per shape, corresponding to the face centres.
        tmp<pointField> centres() const;


    // Search

        //- Get type (inside,outside,mixed,unknown) of point w.r.t. surface.
        //  Only makes sense for closed surfaces.
        volumeType getVolumeType
        (
            const indexedOctree<treeDataFace>&,
            const point&
        ) const;

        //- Does (bb of) shape at index overlap searchBox
        bool overlaps
        (
            const label index,
            const treeBoundBox& searchBox
        ) const;

        //- Does shape at index overlap sphere
        bool overlaps
        (
            const label index,
            const point& centre,
            const scalar radiusSqr
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
