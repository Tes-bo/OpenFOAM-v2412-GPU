/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2022 OpenCFD Ltd.
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
    Foam::isoSurfaceBase

Description
    Low-level components common to various iso-surface algorithms.

Note
    The interpolation samplers currently require a volField for the cell
    values. This is largely a restriction imposed by the point algorithm
    and may be revised in the future.

SourceFiles
    isoSurfaceBase.C
    isoSurfaceBaseNew.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_isoSurfaceBase_H
#define Foam_isoSurfaceBase_H

#include "isoSurfaceParams.H"
#include "bitSet.H"
#include "scalarField.H"
#include "volumeType.H"
#include "volFieldsFwd.H"
#include "MeshedSurface.H"
#include "MeshedSurfacesFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class polyMesh;
class tetCell;

/*---------------------------------------------------------------------------*\
                       Class isoSurfaceBase Declaration
\*---------------------------------------------------------------------------*/

class isoSurfaceBase
:
    public meshedSurface,
    public isoSurfaceParams
{
public:

    // Data Types

        //- The type of cell/face cuts
        enum cutType : uint8_t
        {
            NOTCUT = 0,         //!< Not cut
            CUT = 0x1,          //!< Normal cut
            TETCUT = 0x2,       //!< Cell cut is a tet
            SPHERE = 0x4,       //!< All edges to cell centre cut
            ANYCUT = 0xF,       //!< Any cut type (bitmask)

            UNVISITED = 0x10,   //!< Unvisited
            BLOCKED = 0x20,     //!< Blocked (never cut)
            SPECIAL = 0xF0      //!< Bitmask for specials
        };


protected:

    // Protected typedefs for convenience
    typedef meshedSurface Mesh;

    // Typedef for code transition
    typedef cutType cellCutType;


    // Protected Data

        //- Reference to mesh
        const polyMesh& mesh_;

        //- Cell values
        const scalarField& cVals_;

        //- Point values
        const scalarField& pVals_;

        //- Iso value
        const scalar iso_;


    // Controls, restrictions

        //- Optional boundary faces to ignore.
        //  Eg, Used to exclude cyclicACMI (since duplicate faces)
        bitSet ignoreBoundaryFaces_;


    // Sampling information

        //- For every face, the original cell in mesh
        labelList meshCells_;


    // Protected Member Functions

        //- Count the number of cuts matching the mask type
        //  Checks as bitmask or as zero.
        static label countCutType
        (
            const UList<cutType>& cuts,
            const uint8_t maskValue
        );

        //- Dummy templated interpolate method
        template<class Type>
        tmp<Field<Type>> interpolateTemplate
        (
            const VolumeField<Type>& cellValues,
            const Field<Type>& pointValues
        ) const
        {
            return nullptr;
        }

        //- No copy construct
        isoSurfaceBase(const isoSurfaceBase&) = delete;

        //- No copy assignment
        void operator=(const isoSurfaceBase&) = delete;


public:

    // Typedefs for code transition
    using isoSurfaceParams::algorithmType;
    using isoSurfaceParams::filterType;


    // Constructors

        //- Construct with mesh, cell/point values and iso-value
        isoSurfaceBase
        (
            const polyMesh& mesh,
            const scalarField& cellValues,
            const scalarField& pointValues,
            const scalar iso,
            const isoSurfaceParams& params = isoSurfaceParams()
        );


    // Selector

        //- Create for specified algorithm type
        //  Currently uses hard-code lookups based in isoSurfaceParams
        static autoPtr<isoSurfaceBase> New
        (
            const isoSurfaceParams& params,
            const volScalarField& cellValues,
            const scalarField& pointValues,
            const scalar iso,
            const bitSet& ignoreCells = bitSet()
        );


    // Member Functions

    // Access, Edit

        //- The mesh for which the iso-surface is associated
        const polyMesh& mesh() const noexcept
        {
            return mesh_;
        }

        //- The mesh cell values used for creating the iso-surface
        const scalarField& cellValues() const noexcept
        {
            return cVals_;
        }

        //- The mesh point values used for creating the iso-surface
        const scalarField& pointValues() const noexcept
        {
            return pVals_;
        }

        //- The iso-value associated with the surface
        scalar isoValue() const noexcept
        {
            return iso_;
        }

        //- For each face, the original cell in mesh
        const labelList& meshCells() const noexcept
        {
            return meshCells_;
        }

        //- For each face, the original cell in mesh
        labelList& meshCells() noexcept
        {
            return meshCells_;
        }


    // Helpers

        //- Restore non-BLOCKED state to an UNVISITED state
        static void resetCuts(UList<cutType>& cuts);

        //- Mark ignoreCells as BLOCKED
        label blockCells
        (
            UList<cutType>& cuts,
            const bitSet& ignoreCells
        ) const;

        //- Mark cells inside/outside a (valid) bound box as BLOCKED
        //  The volType is INSIDE or OUTSIDE only
        label blockCells
        (
            UList<cutType>& cuts,
            const boundBox& bb,
            const volumeType::type volType
        ) const;


    // Cutting

        //- Set ignoreBoundaryFaces to ignore cyclics (cyclicACMI)
        void ignoreCyclics();

        //- Populate a list of candidate cell cuts using getCellCutType()
        label calcCellCuts(List<cutType>& cuts) const;

        //- Determine face cut for an individual face
        cutType getFaceCutType(const label facei) const;

        //- Cell cut for an individual cell, with special handling
        //- for TETCUT and SPHERE cuts
        cutType getCellCutType(const label celli) const;


    // Sampling

#undef  declareIsoSurfaceInterpolateMethod
#define declareIsoSurfaceInterpolateMethod(Type)                               \
        /*! \brief interpolate Type cellValues, pointValues on iso-surface */  \
        virtual tmp<Field<Type>>                                               \
        interpolate                                                            \
        (                                                                      \
            const VolumeField<Type>& cellValues,                               \
            const Field<Type>& pointValues                                     \
        ) const;

        declareIsoSurfaceInterpolateMethod(scalar);
        declareIsoSurfaceInterpolateMethod(vector);
        declareIsoSurfaceInterpolateMethod(sphericalTensor);
        declareIsoSurfaceInterpolateMethod(symmTensor);
        declareIsoSurfaceInterpolateMethod(tensor);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
