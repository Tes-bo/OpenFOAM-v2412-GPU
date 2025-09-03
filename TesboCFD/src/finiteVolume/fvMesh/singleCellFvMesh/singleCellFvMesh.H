/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019,2022 OpenCFD Ltd.
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
    Foam::singleCellFvMesh

Description
    fvMesh as subset of other mesh. Consists of one cell and all original
    boundary faces. Useful when manipulating boundary data. Single internal
    cell only needed to be able to manipulate in a standard way.

SourceFiles
    singleCellFvMesh.C
    singleCellFvMeshInterpolate.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_singleCellFvMesh_H
#define Foam_singleCellFvMesh_H

#include "fvPatchFieldMapper.H"
#include "fvMesh.H"
#include "labelListIOList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class singleCellFvMesh Declaration
\*---------------------------------------------------------------------------*/

class singleCellFvMesh
:
    public fvMesh
{
    // Private Data

        //- Fine patch face to agglomeration index addressing
        const labelListIOList patchFaceAgglomeration_;

        //- From patch faces back to agglomeration or fine mesh
        labelListIOList patchFaceMap_;

        //- From fine mesh faces to coarse mesh
        labelIOList reverseFaceMap_;

        //- From coarse points back to original mesh
        labelIOList pointMap_;

        //- From fine points to coarse mesh
        labelIOList reversePointMap_;


    // Private Member Functions

        //- Calculate agglomerated mesh
        void agglomerateMesh(const fvMesh&, const labelListList&);


        //- No copy construct
        singleCellFvMesh(const singleCellFvMesh&) = delete;

        //- No copy assignment
        void operator=(const singleCellFvMesh&) = delete;


public:

        //- Patch field mapper class for agglomerated meshes
        class agglomPatchFieldMapper
        :
            public fvPatchFieldMapper
        {
            // Private data

                const labelListList& addressing_;
                const scalarListList& weights_;
                bool hasUnmapped_;

        public:

                //- Construct given addressing
                agglomPatchFieldMapper
                (
                    const labelListList& addressing,
                    const scalarListList& weights
                )
                :
                    addressing_(addressing),
                    weights_(weights),
                    hasUnmapped_(false)
                {
                    for (const labelList& addr : addressing)
                    {
                        if (addr.empty())
                        {
                            hasUnmapped_ = true;
                            break;
                        }
                    }
                }

                virtual label size() const
                {
                    return addressing_.size();
                }

                virtual bool direct() const
                {
                    return false;
                }

                bool hasUnmapped() const
                {
                    return hasUnmapped_;
                }

                virtual const labelListList& addressing() const
                {
                    return addressing_;
                }

                virtual const scalarListList& weights() const
                {
                    return weights_;
                }
        };



    // Constructors

        //- Construct from fvMesh and no agglomeration
        singleCellFvMesh
        (
            const IOobject& io,
            const fvMesh&,
            const bool doInit=true
        );

        //- Construct from fvMesh and agglomeration of boundary faces.
        //  Agglomeration is per patch, per patch face index the agglomeration
        //  the face goes into.
        singleCellFvMesh
        (
            const IOobject& io,
            const fvMesh&,
            const labelListList& patchFaceAgglomeration,
            const bool doInit=true
        );

        //- Read from IOobject
        singleCellFvMesh(const IOobject& io, const bool doInit=true);


    // Member Functions

        bool agglomerate() const noexcept
        {
            return !patchFaceAgglomeration_.empty();
        }

        //- Fine patch face to agglomeration index addressing
        const labelListList& patchFaceAgglomeration() const noexcept
        {
            return patchFaceAgglomeration_;
        }

        //- From patchFace on this back to original mesh or agglomeration
        const labelListList& patchFaceMap() const noexcept
        {
            return patchFaceMap_;
        }

        //- From point on this back to original mesh
        const labelList& pointMap() const noexcept
        {
            return pointMap_;
        }

        //- From face on original mesh to face on this
        const labelList& reverseFaceMap() const noexcept
        {
            return reverseFaceMap_;
        }

        //- From point on original mesh to point on this (or -1 for removed
        //- points)
        const labelList& reversePointMap() const noexcept
        {
            return reversePointMap_;
        }

        //- Interpolate for overset (unused)
        using fvMesh::interpolate;

        //- Map volField. Internal field set to average, patch fields straight
        //- copies.
        template<class Type>
        tmp<GeometricField<Type, fvPatchField, volMesh>>
        interpolate
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "singleCellFvMeshInterpolate.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
