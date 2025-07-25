/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2022 OpenCFD Ltd.
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
    Foam::fvFieldReconstructor

Description
    Finite volume reconstructor for volume and surface fields.

SourceFiles
    fvFieldReconstructor.C
    fvFieldReconstructorTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_fvFieldReconstructor_H
#define Foam_fvFieldReconstructor_H

#include "PtrList.H"
#include "fvMesh.H"
#include "IOobjectList.H"
#include "fvPatchFieldMapper.H"
#include "labelIOList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class fvFieldReconstructor Declaration
\*---------------------------------------------------------------------------*/

class fvFieldReconstructor
{
    // Private Data

        //- Reconstructed mesh reference
        fvMesh& mesh_;

        //- List of processor meshes
        const PtrList<fvMesh>& procMeshes_;

        //- List of processor face addressing lists
        const PtrList<labelIOList>& faceProcAddressing_;

        //- List of processor cell addressing lists
        const PtrList<labelIOList>& cellProcAddressing_;

        //- List of processor boundary addressing lists
        const PtrList<labelIOList>& boundaryProcAddressing_;

        //- Number of fields reconstructed
        label nReconstructed_;


    // Private Member Functions

        //- No copy construct
        fvFieldReconstructor(const fvFieldReconstructor&) = delete;

        //- No copy assignment
        void operator=(const fvFieldReconstructor&) = delete;


public:

    // Public Classes

        //- Mapper for sizing only - does not do any actual mapping.
        class fvPatchFieldReconstructor
        :
            public fvPatchFieldMapper
        {
            label size_;

        public:

            // Constructors

                //- Construct given size
                fvPatchFieldReconstructor(const label size)
                :
                    size_(size)
                {}


            // Member functions

                label size() const
                {
                    return size_;
                }

                bool direct() const
                {
                    return true;
                }

                bool hasUnmapped() const
                {
                    return false;
                }

                const labelUList& directAddressing() const
                {
                    return labelUList::null();
                }
        };


    // Static Data

        //- Output verbosity when writing
        static int verbose_;


    // Constructors

        //- Construct from components
        fvFieldReconstructor
        (
            fvMesh& mesh,
            const PtrList<fvMesh>& procMeshes,
            const PtrList<labelIOList>& faceProcAddressing,
            const PtrList<labelIOList>& cellProcAddressing,
            const PtrList<labelIOList>& boundaryProcAddressing
        );


    // Member Functions

        //- Return number of fields reconstructed
        label nReconstructed() const noexcept
        {
            return nReconstructed_;
        }

        //- Reconstruct volume internal field
        template<class Type>
        tmp<DimensionedField<Type, volMesh>>
        reconstructField
        (
            const IOobject& fieldObject,
            const PtrList<DimensionedField<Type, volMesh>>& procFields
        ) const;

        //- Read and reconstruct volume internal field
        template<class Type>
        tmp<DimensionedField<Type, volMesh>>
        reconstructInternalField(const IOobject& fieldObject) const;


        //- Reconstruct volume field
        template<class Type>
        tmp<GeometricField<Type, fvPatchField, volMesh>>
        reconstructField
        (
            const IOobject& fieldObject,
            const PtrList<GeometricField<Type, fvPatchField, volMesh>>&
        ) const;

        //- Read and reconstruct volume field
        template<class Type>
        tmp<GeometricField<Type, fvPatchField, volMesh>>
        reconstructVolumeField(const IOobject& fieldObject) const;


        //- Reconstruct surface field
        template<class Type>
        tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        reconstructField
        (
            const IOobject& fieldObject,
            const PtrList<GeometricField<Type, fvsPatchField, surfaceMesh>>&
        ) const;

        //- Read and reconstruct surface field
        template<class Type>
        tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        reconstructSurfaceField(const IOobject& fieldObject) const;


        //- Read, reconstruct and write specified volume internal fields
        template<class Type>
        label reconstructInternalFields
        (
            const UPtrList<const IOobject>& fieldObjects
        );

        //- Read, reconstruct and write specified volume fields
        template<class Type>
        label reconstructVolumeFields
        (
            const UPtrList<const IOobject>& fieldObjects
        );

        //- Read, reconstruct and write specified surface fields
        template<class Type>
        label reconstructSurfaceFields
        (
            const UPtrList<const IOobject>& fieldObjects
        );


        //- Read, reconstruct and write all/selected volume internal fields
        //  An empty wordRes corresponds to select ALL.
        template<class Type>
        label reconstructInternalFields
        (
            const IOobjectList& objects,
            const wordRes& selectedFields = wordRes()
        );

        //- Read, reconstruct and write all/selected volume fields
        //  An empty wordRes corresponds to select ALL.
        template<class Type>
        label reconstructVolumeFields
        (
            const IOobjectList& objects,
            const wordRes& selectedFields = wordRes()
        );

        //- Read, reconstruct and write all/selected surface fields
        //  An empty wordRes corresponds to select ALL.
        template<class Type>
        label reconstructSurfaceFields
        (
            const IOobjectList& objects,
            const wordRes& selectedFields = wordRes()
        );

        //- Reconstruct all known field types
        label reconstructAllFields
        (
            const IOobjectList& objects,
            const wordRes& selectedFields = wordRes()
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fvFieldReconstructorTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
