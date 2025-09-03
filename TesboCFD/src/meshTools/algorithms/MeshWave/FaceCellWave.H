/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2023 OpenCFD Ltd.
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
    Foam::FaceCellWave

Description
    Wave propagation of information through grid. Every iteration
    information goes through one layer of cells. Templated on information
    that is transferred.

    Handles parallel and cyclics and non-parallel cyclics.

    Note: whether to propagate depends on the return value of Type::update
    which returns true (i.e. propagate) if the value changes by more than a
    certain tolerance.
    This tolerance can be very strict for normal face-cell and parallel
    cyclics (we use a value of 0.01 just to limit propagation of small changes)
    but for non-parallel cyclics this tolerance can be critical and if chosen
    too small can lead to non-convergence.

SourceFiles
    FaceCellWave.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_FaceCellWave_H
#define Foam_FaceCellWave_H

#include "bitSet.H"
#include "labelPair.H"
#include "DynamicList.H"
#include "primitiveFieldsFwd.H"
#include "PstreamBuffers.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class polyMesh;
class polyPatch;

/*---------------------------------------------------------------------------*\
                        Class FaceCellWaveBase Declaration
\*---------------------------------------------------------------------------*/

class FaceCellWaveBase
{
protected:

    // Protected Static Data

        static const scalar geomTol_;
        static scalar propagationTol_;


    // Protected Data

        //- Reference to mesh
        const polyMesh& mesh_;

        //- Buffers when updating processor patches
        PstreamBuffers pBufs_;

        //- Track if face has changed
        bitSet changedFace_;

        //- Track if cell has changed
        bitSet changedCell_;

        //- List of changed faces
        DynamicList<label> changedFaces_;

        //- List of changed cells
        DynamicList<label> changedCells_;

        //- Current count of unvisited faces
        label nUnvisitedFaces_;

        //- Current count of unvisited cells
        label nUnvisitedCells_;


public:

    //- Default trackData value (for default template argument)
    static int dummyTrackData_;


    //- Runtime type information
    ClassName("FaceCellWave");


    // Constructors

        //- Construct with mesh reference and set initial sizes
        explicit FaceCellWaveBase(const polyMesh& mesh);


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
        const polyMesh& mesh() const noexcept
        {
            return mesh_;
        }

        //- Current number of changed cells
        label nChangedCells() const noexcept { return changedCells_.size(); }

        //- Current number of changed faces
        label nChangedFaces() const noexcept { return changedFaces_.size(); }

        //- Get number of unvisited cells,
        //- i.e. cells that were not (yet) reached from walking across mesh.
        //  This can happen from
        //  - not enough iterations done
        //  - a disconnected mesh
        //  - a mesh without walls in it
        label nUnvisitedCells() const noexcept { return nUnvisitedCells_; }

        //- Get number of unvisited faces
        label nUnvisitedFaces() const noexcept { return nUnvisitedFaces_; }
};


/*---------------------------------------------------------------------------*\
                           Class FaceCellWave Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class TrackingData = int>
class FaceCellWave
:
    public FaceCellWaveBase
{
protected:

    //- Information tagged with a source or destination id.
    //  With std::pair as lightweight, movable container.
    typedef std::pair<label,Type> taggedInfoType;


    // Protected Data

        //- Optional boundary faces that information should travel through
        const labelPairList explicitConnections_;

        //- Information for all faces
        UList<Type>& allFaceInfo_;

        //- Information for all cells
        UList<Type>& allCellInfo_;

        //- Additional data to be passed into container
        TrackingData& td_;

        // Information exchange for explicit baffle connections
        // Max capacity = 2x number of explicit connections
        DynamicList<taggedInfoType> changedBaffles_;

        //- Contains cyclics
        const bool hasCyclicPatches_;

        //- Contains cyclicAMI
        const bool hasCyclicAMIPatches_;

        //- Number of evaluations
        label nEvals_;


    // Protected Member Functions

        //- Updates cellInfo with information from neighbour.
        //  Updates all statistics.
        bool updateCell
        (
            const label celli,
            const label neighbourFacei,
            const Type& neighbourInfo,
            const scalar tol,
            Type& cellInfo
        );

        //- Updates faceInfo with information from neighbour.
        //  Updates all statistics.
        bool updateFace
        (
            const label facei,
            const label neighbourCelli,
            const Type& neighbourInfo,
            const scalar tol,
            Type& faceInfo
        );

        //- Updates faceInfo with information from same face.
        //  Updates all statistics.
        bool updateFace
        (
            const label facei,
            const Type& neighbourInfo,
            const scalar tol,
            Type& faceInfo
        );


        // Parallel, cyclic

            //- Debugging: check info on both sides of cyclic
            void checkCyclic(const polyPatch& pPatch) const;

            //- Has cyclic patch?
            template<class PatchType>
            bool hasPatch() const;

            //- Merge received patch data into global data
            void mergeFaceInfo
            (
                const polyPatch& patch,
                const label nFaces,
                const labelUList& changedFaces,
                const List<Type>& changedFacesInfo
            );

            //- Extract info for single patch only
            label getChangedPatchFaces
            (
                const polyPatch& patch,
                const label startFacei,
                const label nFaces,
                labelList& changedPatchFaces,
                List<Type>& changedPatchFacesInfo
            ) const;

            //- Handle leaving domain. Implementation referred to Type
            void leaveDomain
            (
                const polyPatch& patch,
                const label nFaces,
                const labelUList& faceLabels,
                List<Type>& faceInfo
            ) const;

            //- Handle leaving domain. Implementation referred to Type
            void enterDomain
            (
                const polyPatch& patch,
                const label nFaces,
                const labelUList& faceLabels,
                List<Type>& faceInfo
            ) const;

            //- Offset face labels by constant value
            static void offset
            (
                const polyPatch& patch,
                const label off,
                const label nFaces,
                labelList& faces
            );

            //- Apply transformation to Type
            void transform
            (
                const tensorField& rotTensor,
                const label nFaces,
                List<Type>& faceInfo
            );

            //- Merge data from across processor boundaries
            //  Transfer changed faces from neighbouring processors.
            void handleProcPatches();

            //- Merge data from across cyclics
            // Transfer changed faces across cyclic halves
            void handleCyclicPatches();

            //- Merge data from across AMI cyclics
            void handleAMICyclicPatches();

            //- Merge data across explicitly provided local connections
            //  These are usually baffles
            void handleExplicitConnections();


        //- No copy construct
        FaceCellWave(const FaceCellWave&) = delete;

        //- No copy assignment
        void operator=(const FaceCellWave&) = delete;


public:

    // Constructors

        //- Construct from mesh.
        //- Use setFaceInfo and iterate() to do actual calculation.
        FaceCellWave
        (
            const polyMesh& mesh,
            UList<Type>& allFaceInfo,
            UList<Type>& allCellInfo,
            TrackingData& td = FaceCellWaveBase::dummyTrackData_
        );

        //- Construct from mesh and list of changed faces with the Type
        //  for these faces. Iterates until nothing changes or maxIter reached.
        //  (maxIter can be 0 or negative). 0 initializes, -1 does not
        FaceCellWave
        (
            const polyMesh& mesh,
            const labelUList& initialChangedFaces,
            const List<Type>& changedFacesInfo,
            UList<Type>& allFaceInfo,
            UList<Type>& allCellInfo,
            const label maxIter,
            TrackingData& td = FaceCellWaveBase::dummyTrackData_
        );

        //- Construct from mesh and explicitly connected boundary faces
        //  and list of changed faces with the Type
        //  for these faces. Iterates until nothing changes or maxIter reached.
        //  (maxIter can be 0 or negative). 0 initializes, -1 does not
        FaceCellWave
        (
            const polyMesh& mesh,
            const labelPairList& explicitConnections,
            const bool handleCyclicAMI,
            const labelUList& initialChangedFaces,
            const List<Type>& changedFacesInfo,
            UList<Type>& allFaceInfo,
            UList<Type>& allCellInfo,
            const label maxIter,
            TrackingData& td = FaceCellWaveBase::dummyTrackData_
        );


    //- Destructor
    virtual ~FaceCellWave() = default;


    // Member Functions

    // Access

        //- Access allFaceInfo
        UList<Type>& allFaceInfo() noexcept
        {
            return allFaceInfo_;
        }

        //- Access allCellInfo
        UList<Type>& allCellInfo() noexcept
        {
            return allCellInfo_;
        }

        //- Additional data to be passed into container
        const TrackingData& data() const noexcept
        {
            return td_;
        }


    // Edit

        //- Set single initial changed face.
        //  This is a noop if the face had already been visited
        void setFaceInfo(const label facei, const Type& faceInfo);

        //- Set initial changed faces
        void setFaceInfo
        (
            const labelUList& changedFaces,
            const List<Type>& changedFacesInfo
        );

        //- Propagate from face to cell.
        //  \return total number of cells (over all processors) changed.
        virtual label faceToCell();

        //- Propagate from cell to face.
        //  \return total number of faces (over all processors) changed.
        //  Note that faces on processor patches are counted twice.
        virtual label cellToFace();

        //- Iterate until no changes or maxIter reached.
        //  \return the number of iterations taken.
        virtual label iterate(const label maxIter);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "FaceCellWave.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
