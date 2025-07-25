/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2024 OpenCFD Ltd.
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
    Foam::faMeshReconstructor

Description
    A bare-bones reconstructor for finiteArea meshes when processor
    meshes are available (in parallel) but an equivalent serial faMesh
    is needed for reconstruction or decomposition.
    In these situations, a serial version of the faMesh is needed,
    but preferably without reconstructing the entire volume mesh.

    It uses the finiteVolume faceProcAddressing in addition to
    the geometric information available from the underlying polyMesh.

    The resulting equivalent faMesh can be used for basic operations,
    but caution should be exercised before attempting large operations.

SourceFiles
    faMeshReconstructor.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_faMeshReconstructor_H
#define Foam_faMeshReconstructor_H

#include "faMesh.H"
#include "primitivePatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class Time;

/*---------------------------------------------------------------------------*\
                     Class faMeshReconstructor Declaration
\*---------------------------------------------------------------------------*/

class faMeshReconstructor
{
    // Private Data

        //- The processor-specific faMesh
        const faMesh& procMesh_;

    // Flags

        //- Problems encountered during construction?
        int errors_;


    // Addressing

        //- Processor face addressing, derived from finite volume information
        labelList faFaceProcAddr_;

        //- Processor boundary addressing
        labelList faBoundaryProcAddr_;

        //- Processor point addressing
        labelList faPointProcAddr_;

        //- Processor edge addressing
        labelList faEdgeProcAddr_;


    // Equivalent surface information

        //- Faces labels for a single patch
        labelList singlePatchFaceLabels_;

        //- Faces for a single patch
        faceList singlePatchFaces_;

        //- Support points for a single patch
        pointField singlePatchPoints_;

        //- Lists of edge-labels (per edge patch) for the single patch
        labelListList singlePatchEdgeLabels_;


    // Demand-driven data

        //- Primitive patch
        mutable autoPtr<primitivePatch> serialPatchPtr_;

        //- Time database for serial meshes
        autoPtr<Time> serialRunTime_;

        //- Dummy volume mesh, used for serial area mesh
        autoPtr<polyMesh> serialVolMesh_;

        //- Equivalent serial area mesh
        autoPtr<faMesh> serialAreaMesh_;


    // Private Member Functions

        //- Calculate all addressing, using finiteVolume faceProcAddressing
        void calcAddressing(const labelUList& fvFaceProcAddr);

        //- Set primitive patch, removing any old one
        void initPatch() const;

        //- Create the serial geometry
        void createMesh();

        //- No copy construct
        faMeshReconstructor(const faMeshReconstructor&) = delete;

        //- No copy assignment
        void operator=(const faMeshReconstructor&) = delete;


public:

    //- Debug flag
    static int debug;


    // Constructors

        //- Construct from components
        explicit faMeshReconstructor
        (
            const faMesh& procMesh,
            IOobjectOption::readOption readVolProcAddr =
                IOobjectOption::MUST_READ
        );

        //- Construct from components
        faMeshReconstructor
        (
            const faMesh& procMesh,
            const labelUList& fvFaceProcAddressing
        );


    //- Destructor
    ~faMeshReconstructor();

        void clearGeom();


    // Member Functions

        //- True if no construct errors encountered
        bool good() const noexcept { return !errors_; }

        //- Processor point addressing
        const labelList& pointProcAddressing() const noexcept
        {
            return faPointProcAddr_;
        }

        //- Processor edge addressing
        const labelList& edgeProcAddressing() const noexcept
        {
            return faEdgeProcAddr_;
        }

        //- Processor face addressing
        const labelList& faceProcAddressing() const noexcept
        {
            return faFaceProcAddr_;
        }

        //- Processor boundary addressing
        const labelList& boundaryProcAddressing() const noexcept
        {
            return faBoundaryProcAddr_;
        }


        //- Serial equivalent patch
        const primitivePatch& patch() const;

        //- Serial equivalent patch
        primitivePatch& patch();

        //- Serial equivalent faMesh
        const faMesh& mesh() const;


    // Write

        //- Write proc addressing
        static void writeAddressing
        (
            const IOobject& io,
            const labelList& faBoundaryProcAddr,
            const labelList& faFaceProcAddr,
            const labelList& faPointProcAddr,
            const labelList& faEdgeProcAddr
        );

        //- Write proc addressing at the polyMesh faceInstances time
        void writeAddressing() const;

        //- Write proc addressing at the given time
        void writeAddressing(const word& timeName) const;

        //- Write mesh information
        static void writeMesh
        (
            const word& timeName,
            const faMesh& fullMesh,
            const labelUList& singlePatchFaceLabels
        );

        //- Write equivalent mesh information at the polyMesh faceInstances time
        void writeMesh() const;

        //- Write equivalent mesh information at the given time
        void writeMesh(const word& timeName) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
