/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
    Copyright (C) 2015-2023 OpenCFD Ltd.
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
    Foam::fvMeshTools

Description
    A collection of tools for operating on an fvMesh.

SourceFiles
    fvMeshTools.C
    fvMeshToolsProcAddr.C
    fvMeshToolsTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_fvMeshTools_H
#define Foam_fvMeshTools_H

#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class fileOperation;
class mapDistributePolyMesh;

/*---------------------------------------------------------------------------*\
                         Class fvMeshTools Declaration
\*---------------------------------------------------------------------------*/

class fvMeshTools
{
    // Private Member Functions

    template<class GeoField>
    static void addPatchFields
    (
        fvMesh&,
        const dictionary& patchFieldDict,
        const word& defaultPatchFieldType,
        const typename GeoField::value_type& defaultPatchValue
    );

    //- Set patchFields according to dictionary
    template<class GeoField>
    static void setPatchFields
    (
        fvMesh& mesh,
        const label patchi,
        const dictionary& patchFieldDict
    );

    //- Set patchFields to value
    template<class GeoField>
    static void setPatchFields
    (
        fvMesh& mesh,
        const label patchi,
        const typename GeoField::value_type& value
    );

    // Remove last patch fields
    template<class GeoField>
    static void trimPatchFields(fvMesh&, const label nPatches);

    template<class GeoField>
    static void reorderPatchFields(fvMesh&, const labelList& oldToNew);

    // Remove trailing patches
    static void trimPatches(fvMesh&, const label nPatches);

    //- Read mesh if available, or create empty mesh with non-proc as per
    //- proc0 mesh.
    //  Has two modes of operation.
    //  If the readHandler is non-nullptr, use it to decide on availability.
    //  Requires:
    //  - all processors to have all patches (and in same order).
    //  - io.instance() set to facesInstance
    static autoPtr<fvMesh> loadOrCreateMeshImpl
    (
        const IOobject& io,
        refPtr<fileOperation>* readHandlerPtr,  // Can be nullptr
        const bool decompose,  // Only used if readHandlerPtr == nullptr
        const bool verbose = false
    );


public:

    // Member Functions

    //- Add patch. Inserts patch before all processor patches.
    //  Supply per field the new patchField per field as a
    //  subdictionary or a default type. If validBoundary call is parallel
    //  synced and all add the same patch with same settings.
    static label addPatch
    (
        fvMesh& mesh,
        const polyPatch& patch,
        const dictionary& patchFieldDict,
        const word& defaultPatchFieldType,
        const bool validBoundary
    );

    //- Change patchField on registered fields according to dictionary
    static void setPatchFields
    (
        fvMesh& mesh,
        const label patchi,
        const dictionary& patchFieldDict
    );

    //- Change patchField to zero on registered fields
    static void zeroPatchFields(fvMesh& mesh, const label patchi);

    //- Reorder and remove trailing patches.
    //  Is parallel synced when validBoundary is true
    static void reorderPatches
    (
        fvMesh&,
        const labelList& oldToNew,
        const label nPatches,
        const bool validBoundary
    );

    //- Remove zero sized patches. All but processor patches are
    //  assumed to be present on all processors (so size will be reduced
    //  if validBoundary). Return map from new
    //  to old patches
    static labelList removeEmptyPatches(fvMesh&, const bool validBoundary);


    //- Set the fvGeometryScheme to basic (to avoid parallel communication)
    static void setBasicGeometry(fvMesh& mesh);

    //- Read mesh or create dummy mesh (0 cells, >0 patches).
    //
    //  Works in two modes according to masterOnlyReading:
    //  true : create dummy mesh on all procs
    //  false: checks locally for mesh directories and only creates dummy mesh
    //         if not present
    static autoPtr<fvMesh> newMesh
    (
        const IOobject& io,
        const bool masterOnlyReading,
        const bool verbose = false
    );

    //- Read mesh if available, or create empty mesh with non-proc as per
    //- proc0 mesh.
    //  Requires:
    //  - all processors to have all patches (and in same order).
    //  - io.instance() set to facesInstance
    static autoPtr<fvMesh> loadOrCreateMesh
    (
        const IOobject& io,
        const bool decompose,
        const bool verbose = false
    );

    //- Read mesh if available, or create empty mesh with non-proc as per
    //- proc0 mesh.
    //  The fileHandler is specified on processors that have a processor mesh.
    //  Generates empty mesh on other processors.
    static autoPtr<fvMesh> loadOrCreateMesh
    (
        const IOobject& io,
        //! Non-null reference if a mesh exists on given processor
        refPtr<fileOperation>& readHandler,
        const bool verbose = true
    );

    //- Create additional fvSchemes/fvSolution files
    static void createDummyFvMeshFiles
    (
        const objectRegistry& parent,
        const word& regionName,
        const bool verbose = false
    );


    //- Read procAddressing components (reconstructing)
    static autoPtr<mapDistributePolyMesh> readProcAddressing
    (
        const fvMesh& procMesh,
        const autoPtr<fvMesh>& baseMeshPtr
    );

    //- Write addressing if decomposing (1 to many) or reconstructing
    //- (many to 1)
    static void writeProcAddressing
    (
        const fvMesh& procMesh,
        const mapDistributePolyMesh& map,
        const bool decompose,
        const fileName& writeHandlerInstance,
        refPtr<fileOperation>& writeHandler
    );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fvMeshToolsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
