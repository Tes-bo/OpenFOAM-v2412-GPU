/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020-2022 OpenCFD Ltd.
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
    Foam::PatchEdgeFaceWave

Description
    Wave propagation of information along patch. Every iteration
    information goes through one layer of faces. Templated on information
    that is transferred.

SourceFiles
    PatchEdgeFaceWave.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_PatchEdgeFaceWave_H
#define Foam_PatchEdgeFaceWave_H

#include "bitSet.H"
#include "scalarField.H"
#include "primitivePatch.H"
#include "vectorTensorTransform.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class polyMesh;

/*---------------------------------------------------------------------------*\
                    Class PatchEdgeFaceWaveBase Declaration
\*---------------------------------------------------------------------------*/

class PatchEdgeFaceWaveBase
{
protected:

    // Protected Static Data

        //- Relative tolerance.
        //  Stop propagation if relative changes less than this tolerance
        //  (responsibility for checking this is up to Type implementation)
        static scalar propagationTol_;


    // Protected Data

        //- Reference to mesh
        const polyMesh& mesh_;

        //- Track if edge has changed
        bitSet changedEdge_;

        //- Track if face has changed
        bitSet changedFace_;

        //- List of changed edges
        DynamicList<label> changedEdges_;

        //- List of changed faces
        DynamicList<label> changedFaces_;

        //- Number of unvisited edges
        label nUnvisitedEdges_;

        //- Number of unvisited faces
        label nUnvisitedFaces_;


public:

    //- Default trackData value (for default template argument)
    static label dummyTrackData_;


    //- Runtime type information
    ClassName("PatchEdgeFaceWave");


    // Constructors

        //- Construct with mesh reference and set initial sizes
        PatchEdgeFaceWaveBase
        (
            const polyMesh& mesh,
            const label nEdges,
            const label nFaces
        );


    // Static Functions

        //- Access to propagation tolerance
        static scalar propagationTol() noexcept
        {
            return propagationTol_;
        }

        //- Change propagation tolerance, return previous value
        static scalar setPropagationTol(const scalar tol) noexcept
        {
            scalar old(propagationTol_);
            propagationTol_ = tol;
            return old;
        }


    // Member Functions

        //- Return access to the mesh
        const polyMesh& mesh() const noexcept { return mesh_; }

        //- Current number of changed edges
        label nChangedEdges() const noexcept { return changedEdges_.size(); }

        //- Current number of changed faces
        label nChangedFaces() const noexcept { return changedFaces_.size(); }

        //- Number of unvisited faces, i.e. faces that were not (yet)
        //- reached from walking across patch.
        //
        //  This can happen from
        //  - not enough iterations done
        //  - a disconnected patch
        //  - a patch without walls in it
        label nUnvisitedFaces() const noexcept { return nUnvisitedFaces_; }

        label nUnvisitedEdges() const noexcept { return nUnvisitedEdges_; }
};


/*---------------------------------------------------------------------------*\
                           Class PatchEdgeFaceWave Declaration
\*---------------------------------------------------------------------------*/

template
<
    class PrimitivePatchType,
    class Type,
    class TrackingData = label
>
class PatchEdgeFaceWave
:
    public PatchEdgeFaceWaveBase
{
    // Private Data

        //- Reference to patch
        const PrimitivePatchType& patch_;

        //- Wall information for all edges
        UList<Type>& allEdgeInfo_;

        //- Information on all patch faces
        UList<Type>& allFaceInfo_;

        //- Additional data to be passed into container
        TrackingData& td_;

        //- Number of evaluations
        label nEvals_;

        // Addressing between edges of patch_ and globalData.coupledPatch()
        labelList patchEdges_;
        labelList coupledEdges_;
        bitSet sameEdgeOrientation_;


    // Private Member Functions

        //- Updates edgeInfo with information from neighbour.
        //  Updates all statistics.
        bool updateEdge
        (
            const label edgeI,
            const label neighbourFacei,
            const Type& neighbourInfo,
            Type& edgeInfo
        );

        //- Updates faceInfo with information from neighbour.
        //  Updates all statistics.
        bool updateFace
        (
            const label facei,
            const label neighbourEdgeI,
            const Type& neighbourInfo,
            Type& faceInfo
        );

        //- Update coupled edges
        void syncEdges();

        //- No copy construct
        PatchEdgeFaceWave(const PatchEdgeFaceWave&) = delete;

        //- No copy assignment
        void operator=(const PatchEdgeFaceWave&) = delete;


public:

    // Constructors

        //- Construct from patch, list of changed edges with the Type
        //- for these edges.
        //  Obtains work arrays to operate on, one of size
        //  number of patch edges, the other number of patch faces.
        //  Iterates until nothing changes or maxIter reached.
        //  (maxIter can be 0)
        PatchEdgeFaceWave
        (
            const polyMesh& mesh,
            const PrimitivePatchType& patch,
            const labelList& initialEdges,
            const List<Type>& initialEdgesInfo,
            UList<Type>& allEdgeInfo,
            UList<Type>& allFaceInfo,
            const label maxIter,
            TrackingData& td = PatchEdgeFaceWaveBase::dummyTrackData_
        );

        //- Construct from patch.
        //  Use setEdgeInfo() and iterate() to do actual calculation
        PatchEdgeFaceWave
        (
            const polyMesh& mesh,
            const PrimitivePatchType& patch,
            UList<Type>& allEdgeInfo,
            UList<Type>& allFaceInfo,
            TrackingData& td = PatchEdgeFaceWaveBase::dummyTrackData_
        );


    // Member Functions

        //- Access allEdgeInfo
        UList<Type>& allEdgeInfo() const noexcept
        {
            return allEdgeInfo_;
        }

        //- Access allFaceInfo
        UList<Type>& allFaceInfo() const noexcept
        {
            return allFaceInfo_;
        }

        //- Additional data to be passed into container
        const TrackingData& data() const noexcept
        {
            return td_;
        }

        //- Copy initial data into allEdgeInfo_
        void setEdgeInfo
        (
            const labelList& changedEdges,
            const List<Type>& changedEdgesInfo
        );

        //- Propagate from edge to face.
        //  \return total number of faces (over all processors) changed.
        label edgeToFace();

        //- Propagate from face to edge.
        //  \return total number of edges (over all processors) changed.
        label faceToEdge();

        //- Iterate until no changes or maxIter reached.
        //  \return actual number of iterations.
        label iterate(const label maxIter);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Update operation
template
<
    class PrimitivePatchType,
    class Type,
    class TrackingData = int
>
class updateOp
{
    //- Additional data to be passed into container
    const polyMesh& mesh_;
    const PrimitivePatchType& patch_;
    const scalar tol_;
    TrackingData& td_;

public:

    updateOp
    (
        const polyMesh& mesh,
        const PrimitivePatchType& patch,
        const scalar tol,
        TrackingData& td
    )
    :
        mesh_(mesh),
        patch_(patch),
        tol_(tol),
        td_(td)
    {}

    void operator()(Type& x, const Type& y) const
    {
        if (y.valid(td_))
        {
            x.updateEdge(mesh_, patch_, y, true, tol_, td_);
        }
    }
};


//- Transform operation
template
<
    class PrimitivePatchType,
    class Type,
    class TrackingData = int
>
class transformOp
{
    //- Additional data to be passed into container
    const polyMesh& mesh_;
    const PrimitivePatchType& patch_;
    const scalar tol_;
    TrackingData& td_;

public:

    transformOp
    (
        const polyMesh& mesh,
        const PrimitivePatchType& patch,
        const scalar tol,
        TrackingData& td
    )
    :
        mesh_(mesh),
        patch_(patch),
        tol_(tol),
        td_(td)
    {}

    void operator()
    (
        const vectorTensorTransform& vt,
        const bool forward,
        List<Type>& fld
    ) const
    {
        if (forward)
        {
            forAll(fld, i)
            {
                fld[i].transform(mesh_, patch_, vt.R(), tol_, td_);
            }
        }
        else
        {
            forAll(fld, i)
            {
                fld[i].transform(mesh_, patch_, vt.R().T(), tol_, td_);
            }
        }
    }
};

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PatchEdgeFaceWave.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
