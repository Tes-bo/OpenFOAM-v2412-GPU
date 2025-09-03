/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
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
    Foam::parFvFieldDistributor

Description
    Finite volume reconstructor for volume and surface fields.

    Runs in parallel.
    Reconstructs/redistributes from procMesh to baseMesh.
    baseMesh is non-zero cells on processor0 only.

SourceFiles
    parFvFieldDistributor.C
    parFvFieldDistributorFields.C
    parFvFieldDistributorTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_parFvFieldDistributor_H
#define Foam_parFvFieldDistributor_H

#include "PtrList.H"
#include "fvMesh.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class mapDistributePolyMesh;
class mapDistributeBase;
class IOobjectList;

/*---------------------------------------------------------------------------*\
                    Class parFvFieldDistributor Declaration
\*---------------------------------------------------------------------------*/

class parFvFieldDistributor
{
    // Private Data

        //- Processor mesh reference (source mesh)
        const fvMesh& srcMesh_;

        //- Destination mesh reference (eg, reconstructed mesh)
        fvMesh& tgtMesh_;

        //- Distribution map reference
        const mapDistributePolyMesh& distMap_;

        //- Patch mappers
        PtrList<mapDistributeBase> patchFaceMaps_;

        //- Storage for dummy handler (when using bool control)
        refPtr<fileOperation> dummyHandler_;

        //- Write control via a file handler
        refPtr<fileOperation>& writeHandler_;

        //- Write control as a bool
        Switch isWriteProc_;


    // Private Member Functions

        //- Construct per-patch addressing
        void createPatchFaceMaps();


public:

    //- Output verbosity when writing
    static int verbose_;


    // Generated Methods

        //- No copy construct
        parFvFieldDistributor(const parFvFieldDistributor&) = delete;

        //- No copy assignment
        void operator=(const parFvFieldDistributor&) = delete;


    // Constructors

        //- Construct from components
        //
        //  \param srcMesh  The source mesh (eg, processor)
        //  \param tgtMesh  The target mesh (eg, reconstructed)
        //  \param distMap  The distribution map
        //  \param isWriteProc  Tagged for output writing (on this proc)
        parFvFieldDistributor
        (
            const fvMesh& srcMesh,
            fvMesh& tgtMesh,
            const mapDistributePolyMesh& distMap,
            const bool isWriteProc
        );

        //- Construct from components
        //
        //  \param srcMesh  The source mesh (eg, processor)
        //  \param tgtMesh  The target mesh (eg, reconstructed)
        //  \param distMap  The distribution map
        //  \param writeHandler  Handler for output writing (on this proc)
        parFvFieldDistributor
        (
            const fvMesh& srcMesh,
            fvMesh& tgtMesh,
            const mapDistributePolyMesh& distMap,
            refPtr<fileOperation>& writeHandler
        );


    // Member Functions

        //- Helper: reconstruct and write mesh points
        //  (note: should be moved to something like processorMeshes class)
        void reconstructPoints();

        //- Distribute all fields for known field types
        void distributeAllFields
        (
            const IOobjectList& objects,
            const wordRes& selectedFields
        ) const;

        //- Redistribute volume internal field
        template<class Type>
        tmp<DimensionedField<Type, volMesh>>
        distributeField
        (
            const DimensionedField<Type, volMesh>&
        ) const;

        //- Read and distribute volume internal field
        template<class Type>
        tmp<DimensionedField<Type, volMesh>>
        distributeInternalField(const IOobject& fieldObject) const;


        //- Redistribute volume field
        template<class Type>
        tmp<GeometricField<Type, fvPatchField, volMesh>>
        distributeField
        (
            const GeometricField<Type, fvPatchField, volMesh>& fld
        ) const;

        //- Read and distribute volume field
        template<class Type>
        tmp<GeometricField<Type, fvPatchField, volMesh>>
        distributeVolumeField(const IOobject& fieldObject) const;


        //- Redistribute surface field
        template<class Type>
        tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        distributeField
        (
            const GeometricField<Type, fvsPatchField, surfaceMesh>&
        ) const;

        //- Read and distribute surface field
        template<class Type>
        tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        distributeSurfaceField(const IOobject& fieldObject) const;


        //- Read, redistribute and write all/selected volume internal fields
        template<class Type>
        label distributeInternalFields
        (
            const IOobjectList& objects,
            const wordRes& selectedFields = wordRes()
        ) const;

        //- Read, redistribute and write all/selected volume fields
        template<class Type>
        label distributeVolumeFields
        (
            const IOobjectList& objects,
            const wordRes& selectedFields = wordRes()
        ) const;

        //- Read, reconstruct and write all/selected surface fields
        template<class Type>
        label distributeSurfaceFields
        (
            const IOobjectList& objects,
            const wordRes& selectedFields = wordRes()
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "parFvFieldDistributorTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
