/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Foam::faMeshDistributor

Description
    Holds a reference to the original mesh (the baseMesh)
    and optionally to a subset of that mesh (the subMesh)
    with mapping lists for points, faces, and cells.

SourceFiles
    faMeshDistributor.C
    faMeshDistributorNew.C
    faMeshDistributorTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_faMeshDistributor_H
#define Foam_faMeshDistributor_H

#include "faMesh.H"
#include "mapDistributePolyMesh.H"
#include "areaFieldsFwd.H"
#include "edgeFieldsFwd.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class IOobjectList;

/*---------------------------------------------------------------------------*\
                      Class faMeshDistributor Declaration
\*---------------------------------------------------------------------------*/

class faMeshDistributor
{
    // Private Data

        //- The source mesh reference
        const faMesh& srcMesh_;

        //- The destination mesh reference
        const faMesh& tgtMesh_;

        //- Distribution map reference (faMesh)
        const mapDistributePolyMesh& distMap_;

        //- Internal edge mapper
        mutable std::unique_ptr<mapDistributeBase> internalEdgeMapPtr_;

        //- Patch edge mappers
        mutable PtrList<mapDistributeBase> patchEdgeMaps_;

        //- Storage for dummy handler (when using bool control)
        refPtr<fileOperation> dummyHandler_;

        //- Write control via a file handler
        refPtr<fileOperation>& writeHandler_;

        //- Write control as a bool
        Switch isWriteProc_;


    // Private Member Functions

        //- Construct internal edge mapping
        void createInternalEdgeMap() const;

        //- Construct per-patch edge mapping
        void createPatchMaps() const;

        //- Debug: check addressing
        void checkAddressing() const;

        //- Construct reconstruct mapping
        static mapDistributePolyMesh createReconstructMap
        (
            const faMesh& mesh,
            const autoPtr<faMesh>& baseMeshPtr,
            const labelUList& faceProcAddr,
            const labelUList& edgeProcAddr,
            const labelUList& pointProcAddr,
            const labelUList& boundaryProcAddr
        );


public:

    //- Output verbosity when writing
    static int verbose_;


    // Generated Methods

        //- No copy construct
        faMeshDistributor(const faMeshDistributor&) = delete;

        //- No copy assignment
        void operator=(const faMeshDistributor&) = delete;


    // Constructors

        //- Construct from components, using bool to control writing
        faMeshDistributor
        (
            const faMesh& srcMesh,
            const faMesh& tgtMesh,
            const mapDistributePolyMesh& faDistMap,
            const bool isWriteProc
        );

        //- Construct from components, using file handler to control writing
        faMeshDistributor
        (
            const faMesh& srcMesh,
            const faMesh& tgtMesh,
            const mapDistributePolyMesh& faDistMap,
            refPtr<fileOperation>& writeHandler
        );


    // Static Methods

        //- Distribute mesh according to the given (volume) mesh distribution.
        //  Uses 'tgtPolyMesh' for the new mesh
        static mapDistributePolyMesh distribute
        (
            const faMesh& oldMesh,
            const mapDistributePolyMesh& distMap,  //! From polyMesh
            const polyMesh& tgtPolyMesh,
            autoPtr<faMesh>& newMeshPtr
        );

        //- Distribute mesh according to the given (volume) mesh distribution.
        //  Re-uses polyMesh from oldMesh for the new mesh
        static mapDistributePolyMesh distribute
        (
            const faMesh& oldMesh,
            const mapDistributePolyMesh& distMap,  //! From polyMesh
            autoPtr<faMesh>& newMeshPtr
        );


    // Member Functions

    // Field Mapping

        //- Read, distribute and write all/selected point field types
        //- (scalar, vector, ... types)
        label distributeAllFields
        (
            const IOobjectList& objects,
            const wordRes& selectedFields = wordRes()
        ) const;

        //- Distribute area field
        template<class Type>
        tmp<GeometricField<Type, faPatchField, areaMesh>>
        distributeField
        (
            const GeometricField<Type, faPatchField, areaMesh>& fld
        ) const;

        //- Distribute edge field
        template<class Type>
        tmp<GeometricField<Type, faePatchField, edgeMesh>>
        distributeField
        (
            const GeometricField<Type, faePatchField, edgeMesh>& fld
        ) const;

        //- Read and distribute area field
        template<class Type>
        tmp<GeometricField<Type, faPatchField, areaMesh>>
        distributeAreaField
        (
            const IOobject& fieldObject
        ) const;

        //- Read and distribute edge field
        template<class Type>
        tmp<GeometricField<Type, faePatchField, edgeMesh>>
        distributeEdgeField
        (
            const IOobject& fieldObject
        ) const;

        //- Read, distribute and write all/selected area fields
        template<class Type>
        label distributeAreaFields
        (
            const IOobjectList& objects,
            const wordRes& selectedFields = wordRes()
        ) const;

        //- Read, distribute and write all/selected area fields
        template<class Type>
        label distributeEdgeFields
        (
            const IOobjectList& objects,
            const wordRes& selectedFields = wordRes()
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#ifdef NoRepository
#   include "faMeshDistributorTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
