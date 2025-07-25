/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
    Copyright (C) 2021 OpenCFD Ltd.
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
    Foam::faMeshDecomposition

Description
    Automatic faMesh decomposition class

Author
    Zeljko Tukovic, FSB Zagreb
    Hrvoje Jasak, Wikki Ltd.

SourceFiles
    faMeshDecomposition.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_faMeshDecomposition_H
#define Foam_faMeshDecomposition_H

#include "polyMesh.H"
#include "faMesh.H"
#include "labelList.H"
#include "PtrList.H"
#include "point.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class faMeshDecomposition Declaration
\*---------------------------------------------------------------------------*/

class faMeshDecomposition
:
    public faMesh
{
    // Private Data

        //- Number of processors in decomposition
        label nProcs_;

        //- Is the decomposition data to be distributed for each processor
        bool distributed_;

        //- Are globalFaceZones being used
        bool hasGlobalFaceZones_;

        //- Processor label for each cell
        labelList faceToProc_;

        //- Face labels for each processor mesh
        labelListList procFaceLabels_;

        //-
        List<Map<label>> procMeshEdgesMap_;

        //- Number of internal edges for each processor mesh
        labelList procNInternalEdges_;

        //- Edge labels for patches of processor meshes
        List<labelListList> procPatchEdgeLabels_;

        //- Labels of points for each processor
        labelListList procPatchPointAddressing_;

        //- Labels of edges for each processor
        labelListList procPatchEdgeAddressing_;

        //- Labels of edges for each processor
        labelListList procEdgeAddressing_;

        //- Labels of faces for each processor
        labelListList procFaceAddressing_;

        //- Original patch index for every processor patch
        labelListList procBoundaryAddressing_;

        //- Sizes for processor mesh patches
        // Excludes inter-processor boundaries
        labelListList procPatchSize_;

        //- Start indices for processor patches
        // Excludes inter-processor boundaries
        labelListList procPatchStartIndex_;

        //- Neighbour processor ID for inter-processor boundaries
        labelListList procNeighbourProcessors_;

        //- Sizes for inter-processor patches
        labelListList procProcessorPatchSize_;

        //- Start indices for inter-processor patches
        labelListList procProcessorPatchStartIndex_;

        //- List of globally shared point labels
        labelList globallySharedPoints_;

        //- Are there cyclic-parallel faces
        bool cyclicParallel_;


    // Private Member Functions

        void distributeFaces();


public:

    // Constructors

        //- Construct from components.
        //- Values will come from the volume decomposition
        //  \param mesh the fvMesh
        //  \param nProcessors the number of processors
        //  \param params additional parameters, sent to updateParameters
        faMeshDecomposition
        (
            const polyMesh& mesh,
            const label nProcessors,
            const dictionary& params = dictionary::null
        );


    //- Destructor
    ~faMeshDecomposition() = default;


    // Member Functions

    // Settings

        //- Number of processor in decomposition
        label nProcs() const noexcept
        {
            return nProcs_;
        }

        //- Is decomposition data to be distributed for each processor
        bool distributed() const noexcept
        {
            return distributed_;
        }

        //- Change distributed flag
        bool distributed(const bool on) noexcept
        {
            bool old(distributed_);
            distributed_ = on;
            return old;
        }

        //- Are global face zones used
        bool useGlobalFaceZones() const noexcept
        {
            return distributed_;
        }

        //- Change global face zones flag
        bool useGlobalFaceZones(const bool on) noexcept
        {
            bool old(hasGlobalFaceZones_);
            hasGlobalFaceZones_ = on;
            return old;
        }

        //- Update flags based on the decomposition model settings
        //  Sets "distributed", detects presence of "globalFaceZones"
        void updateParameters(const dictionary& params);


    // Mappings

        //- Face-processor decomposition labels
        const labelList& faceToProc() const noexcept
        {
            return faceToProc_;
        }


    // Decompose

        //- Decompose mesh
        void decomposeMesh();

        //- Write decomposition
        bool writeDecomposition();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
