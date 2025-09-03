/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2021 OpenCFD Ltd.
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
    Foam::fvMeshSubsetProxy

Description
    Simple proxy for holding a mesh, or mesh-subset.
    The subMeshes are currently limited to cellSet or cellZone definitions.

SourceFiles
    fvMeshSubsetProxy.C
    fvMeshSubsetProxyTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_fvMeshSubsetProxy_H
#define Foam_fvMeshSubsetProxy_H

#include "fvMeshSubset.H"
#include "wordRes.H"
#include "zeroGradientFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class fvMeshSubsetProxy Declaration
\*---------------------------------------------------------------------------*/

class fvMeshSubsetProxy
{
public:

    // Data Types

        //- Internal bookkeeping for subset type
        enum subsetType
        {
            NONE,       //!< No subset
            SET,        //!< Subset with a cellSet
            ZONE,       //!< Subset with a cellZone
            ZONES       //!< Subset with multiple cellZones
        };


private:

    // Private Data

        //- Reference to mesh
        fvMesh& baseMesh_;

        //- Subsetting engine
        fvMeshSubset subsetter_;

        //- Patch ID for exposed internal faces
        label exposedPatchId_;

        //- The subsetting type
        subsetType type_;

        //- Name of the cellSet/cellZone (or empty)
        word name_;

        //- Selection for multiple cell zones
        wordRes names_;

        //- The (cached) cell selection
        bitSet selectedCells_;


    // Private Member Functions

        //- Clear, set to NONE
        void clearOut();

        //- No copy construct
        fvMeshSubsetProxy(const fvMeshSubsetProxy&) = delete;

        //- No copy assignment
        void operator=(const fvMeshSubsetProxy&) = delete;


public:

    // Constructors

        //- Construct a pass-through proxy. No correct() invoked or required.
        explicit fvMeshSubsetProxy(fvMesh& baseMesh);

        //- Construct from components
        fvMeshSubsetProxy
        (
            fvMesh& baseMesh,
            const subsetType type,
            const word& selectionName,
            label exposedPatchId = -1
        );

        //- Construct from components. The subsetType is ZONES.
        fvMeshSubsetProxy
        (
            fvMesh& baseMesh,
            const wordRes& zoneNames,
            label exposedPatchId = -1
        );

        //- Construct from components. The subsetType is ZONES.
        fvMeshSubsetProxy
        (
            fvMesh& baseMesh,
            wordRes&& zoneNames,
            label exposedPatchId = -1
        );


    // Member Functions

    // Access

        //- The entire base mesh
        const fvMesh& baseMesh() const noexcept
        {
            return baseMesh_;
        }

        //- The mesh subsetter
        const fvMeshSubset& subsetter() const noexcept
        {
            return subsetter_;
        }

        //- The mesh subsetter
        fvMeshSubset& subsetter() noexcept
        {
            return subsetter_;
        }

        //- True if sub-mesh should be used
        bool useSubMesh() const noexcept
        {
            return type_ != subsetType::NONE;
        }

        //- Access either base-mesh or sub-mesh
        inline const fvMesh& mesh() const
        {
            if (useSubMesh())
            {
                return subsetter_.subMesh();
            }

            return baseMesh_;
        }

        //- The associated (set or zone) name if any.
        const word& name() const noexcept
        {
            return name_;
        }

        //- The current cell selection, when subsetting is active
        const bitSet& selectedCells() const noexcept
        {
            return selectedCells_;
        }


    // Edit

        //- Define the zones selection, subset the mesh accordingly
        void resetZones(const wordRes& zoneNames);

        //- Update of mesh subset.
        //  Return true if the subset changed from previous call.
        bool correct(bool verbose = false);

        //- Read mesh. Correct on topo-change
        polyMesh::readUpdateState readUpdate();


    // Fields

        //- Construct volField (with zeroGradient) from an internal field
        template<class Type>
        static tmp<GeometricField<Type, fvPatchField, volMesh>>
        zeroGradientField
        (
            const DimensionedField<Type, volMesh>& df
        );

        //- Convert an internal field to a volume field (with zeroGradient)
        template<class Type>
        static tmp<GeometricField<Type, fvPatchField, volMesh>>
        interpolateInternal
        (
            const fvMeshSubset& subsetter,
            const DimensionedField<Type, volMesh>& df
        );

        //- Convert an internal field to a volume field (with zeroGradient)
        //  Currently no proper memory reuse
        template<class Type>
        static tmp<GeometricField<Type, fvPatchField, volMesh>>
        interpolateInternal
        (
            const fvMeshSubset& subsetter,
            const tmp<DimensionedField<Type, volMesh>>& tdf
        );

        //- Wrapper for field or the subsetted field.
        //  Pass through or forward to fvMeshSubset::interpolate()
        template<class GeoField>
        static tmp<GeoField> interpolate
        (
            const fvMeshSubset& subsetter,
            const GeoField& fld
        );

        //- Wrapper for field or the subsetted field.
        //  Pass through or forward to fvMeshSubset::interpolate()
        template<class GeoField>
        static tmp<GeoField> interpolate
        (
            const fvMeshSubset& subsetter,
            const tmp<GeoField>& fld
        );


    // Fields

        //- Convert an internal field to a volume field (with zeroGradient)
        template<class Type>
        tmp<GeometricField<Type, fvPatchField, volMesh>>
        interpolateInternal
        (
            const DimensionedField<Type, volMesh>& df
        ) const;

        //- Convert an internal field to a volume field (with zeroGradient)
        //  Currently no proper memory reuse
        template<class Type>
        tmp<GeometricField<Type, fvPatchField, volMesh>>
        interpolateInternal
        (
            const tmp<DimensionedField<Type, volMesh>>& tdf
        ) const;


        //- Wrapper for field or the subsetted field.
        //  Pass through or forward to fvMeshSubset::interpolate()
        template<class GeoField>
        tmp<GeoField> interpolate(const GeoField& fld) const;

        //- Wrapper for field or the subsetted field.
        //  Pass through or forward to fvMeshSubset::interpolate()
        template<class GeoField>
        tmp<GeoField> interpolate(const tmp<GeoField>& fld) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fvMeshSubsetProxyTemplates.C"
#endif

#endif

// ************************************************************************* //
