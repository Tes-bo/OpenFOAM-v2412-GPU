/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
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
    Foam::faFieldReconstructor

Description
    Finite Area area and edge fields reconstruction.

Author
    Zeljko Tukovic, FSB Zagreb
    Hrvoje Jasak, Wikki Ltd.

SourceFiles
    faFieldReconstructor.C
    faFieldReconstructorTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_faFieldReconstructor_H
#define Foam_faFieldReconstructor_H

#include "PtrList.H"
#include "faMesh.H"
#include "IOobjectList.H"
#include "faPatchFieldMapper.H"
#include "labelIOList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class faFieldReconstructor Declaration
\*---------------------------------------------------------------------------*/

class faFieldReconstructor
{
    // Private Data

        //- Reconstructed mesh reference
        faMesh& mesh_;

        //- List of processor meshes
        const PtrList<faMesh>& procMeshes_;

        //- List of processor edge addressing lists
        const PtrList<labelIOList>& edgeProcAddressing_;

        //- List of processor face addressing lists
        const PtrList<labelIOList>& faceProcAddressing_;

        //- List of processor boundary addressing lists
        const PtrList<labelIOList>& boundaryProcAddressing_;

        //- Number of fields reconstructed
        label nReconstructed_;


    // Private Member Functions

        //- No copy construct
        faFieldReconstructor(const faFieldReconstructor&) = delete;

        //- No copy assignment
        void operator=(const faFieldReconstructor&) = delete;


public:

    // Public Classes

        class faPatchFieldReconstructor
        :
            public faPatchFieldMapper
        {
            label size_;
            label sizeBeforeMapping_;

        public:

            // Constructors

                //- Construct given size
                faPatchFieldReconstructor
                (
                    const label size,
                    const label sizeBeforeMapping
                )
                :
                    size_(size),
                    sizeBeforeMapping_(sizeBeforeMapping)
                {}


            // Member Functions

                virtual label size() const
                {
                    return size_;
                }

                virtual label sizeBeforeMapping() const
                {
                    return sizeBeforeMapping_;
                }

                virtual bool direct() const
                {
                    return true;
                }

                virtual bool hasUnmapped() const
                {
                    return false;
                }

                virtual const labelUList& directAddressing() const
                {
                    return labelUList::null();
                }
        };


    // Static Data

        //- Output verbosity when writing
        static int verbose_;


    // Constructors

        //- Construct from components
        faFieldReconstructor
        (
            faMesh& mesh,                       // Target mesh
            const PtrList<faMesh>& procMeshes,  // Source meshes
            const PtrList<labelIOList>& edgeProcAddressing,
            const PtrList<labelIOList>& faceProcAddressing,
            const PtrList<labelIOList>& boundaryProcAddressing
        );


    // Member Functions

        //- Return number of fields reconstructed
        label nReconstructed() const noexcept
        {
            return nReconstructed_;
        }

        //- Reconstruct area field
        template<class Type>
        tmp<GeometricField<Type, faPatchField, areaMesh>>
        reconstructField
        (
            const IOobject& fieldObject,
            const PtrList<GeometricField<Type, faPatchField, areaMesh>>&
        ) const;

        //- Read and reconstruct area field
        template<class Type>
        tmp<GeometricField<Type, faPatchField, areaMesh>>
        reconstructAreaField(const IOobject& fieldObject);


        //- Reconstruct edge field
        template<class Type>
        tmp<GeometricField<Type, faePatchField, edgeMesh>>
        reconstructField
        (
            const IOobject& fieldObject,
            const PtrList<GeometricField<Type, faePatchField, edgeMesh>>&
        ) const;

        //- Read and reconstruct edge field
        template<class Type>
        tmp<GeometricField<Type, faePatchField, edgeMesh>>
        reconstructEdgeField(const IOobject& fieldObject);


        //- Read, reconstruct and write specified area fields
        template<class Type>
        label reconstructAreaFields
        (
            const UPtrList<const IOobject>& fieldObjects
        );

        //- Read, reconstruct and write specified edge fields
        template<class Type>
        label reconstructEdgeFields
        (
            const UPtrList<const IOobject>& fieldObjects
        );


        //- Read, reconstruct and write all/selected area fields
        //  An empty wordRes corresponds to select ALL.
        template<class Type>
        label reconstructAreaFields
        (
            const IOobjectList& objects,
            const wordRes& selectedFields = wordRes()
        );

        //- Read, reconstruct and write all/selected edge fields
        //  An empty wordRes corresponds to select ALL.
        template<class Type>
        label reconstructEdgeFields
        (
            const IOobjectList& objects,
            const wordRes& selectedFields = wordRes()
        );

        //- Reconstruct all supported area/edge field types
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
#   include "faFieldReconstructorTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
