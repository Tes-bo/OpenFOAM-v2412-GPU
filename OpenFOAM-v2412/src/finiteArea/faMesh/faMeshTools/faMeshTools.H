/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022-2024 OpenCFD Ltd.
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
    Foam::faMeshTools

Description
    A collection of tools for operating on an faMesh.

SourceFiles
    faMeshTools.C
    faMeshToolsChecks.C
    faMeshToolsProcAddr.C
    faMeshToolsTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_faMeshTools_H
#define Foam_faMeshTools_H

#include "faMesh.H"
#include "areaFieldsFwd.H"
#include "edgeFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class fileOperation;
class mapDistributePolyMesh;
class polyMesh;
class IOobject;

/*---------------------------------------------------------------------------*\
                         Class faMeshTools Declaration
\*---------------------------------------------------------------------------*/

class faMeshTools
{
    // Private Member Functions

    // Read mesh if available. Otherwise create empty mesh with same non-proc
    // patches as proc0 mesh.
    //  Has two modes of operation.
    //  If the readHandler is non-nullptr, use it to decide on availability.
    //  Requires:
    //  - all processors to have all patches (and in same order).
    //  - io.instance() set to facesInstance
    static autoPtr<faMesh> loadOrCreateMeshImpl
    (
        const IOobject& io,
        refPtr<fileOperation>* readHandlerPtr,  // Can be nullptr
        const polyMesh& pMesh,
        const bool decompose,  // Only used if readHandlerPtr == nullptr
        const bool verbose = false
    );



public:

    //- Unregister the faMesh from its associated polyMesh
    //- to prevent triggering on polyMesh changes etc.
    static void unregisterMesh(const faMesh& mesh);

    //- Force creation of everything that might vaguely be used by patches.
    //  This is fairly horrible.
    static void forceDemandDriven(faMesh& mesh);


    //- Read mesh or create dummy mesh (0 faces, >0 patches).
    //  Works in two modes according to masterOnlyReading:
    //  true : create a dummy mesh for all procs
    //  false: checks locally for mesh directories and only creates dummy mesh
    //         if not present
    static autoPtr<faMesh> newMesh
    (
        //! The IOobject describes the base polyMesh
        const IOobject& io,
        const polyMesh& pMesh,
        const bool masterOnlyReading,
        const bool verbose = false
    );

    // Read mesh if available. Otherwise create empty mesh with same non-proc
    // patches as proc0 mesh. Requires:
    //  - all processors to have all patches (and in same order).
    //  - io.instance() set to facesInstance
    static autoPtr<faMesh> loadOrCreateMesh
    (
        //! The IOobject describes the base polyMesh
        const IOobject& io,
        const polyMesh& pMesh,
        const bool decompose,
        const bool verbose = false
    );

    // Read mesh if available. Otherwise create empty mesh with same non-proc
    // patches as proc0 mesh. Requires:
    //  - all processors to have all patches (and in same order).
    //  - io.instance() set to facesInstance
    static autoPtr<faMesh> loadOrCreateMesh
    (
        //! The IOobject describes the base polyMesh
        const IOobject& io,
        const polyMesh& pMesh,
        //! Non-null reference if a mesh exists on given processor
        refPtr<fileOperation>& readHandler,
        const bool verbose = false
    );


    //- Read decompose/reconstruct addressing
    static mapDistributePolyMesh readProcAddressing
    (
        const faMesh& mesh,
        const autoPtr<faMesh>& baseMeshPtr
    );

    //- Write decompose/reconstruct addressing
    //
    //  \param mesh
    //  \param faDistMap
    //  \param decompose  running in decompose vs reconstruct mode
    //  \param writeHandler file handler
    //  \param procMesh   (optional) processor mesh in reconstruct mode
    //
    //  \note Since the faMesh holds a reference to a polyMesh,
    //  in reconstruct mode it will refer to the base mesh, but
    //  we need a means to proc addressing into the processor locations.
    //  This is the purpose of the additional procMesh reference
    static void writeProcAddressing
    (
        const faMesh& mesh,
        const mapDistributePolyMesh& faDistMap,
        const bool decompose,
        refPtr<fileOperation>& writeHandler,
        const faMesh* procMesh = nullptr
    );


    //- Flatten an edge field into linear addressing
    //  Optionally use primitive patch edge ordering
    template<class Type>
    static tmp<Field<Type>> flattenEdgeField
    (
        const GeometricField<Type, faePatchField, edgeMesh>& fld,
        const bool primitiveOrdering = false
    );

    //- Report mesh information
    static void printMeshChecks
    (
        const faMesh& mesh,
        const int verbose = 1
    );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "faMeshToolsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
