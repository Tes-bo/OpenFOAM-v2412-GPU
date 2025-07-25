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
    Foam::sampledMeshedSurface

Description
    A sampledSurface from a meshed surface.
    It samples on the points/faces of the meshed surface.

    - it either samples cells or (non-coupled) boundary faces

    - 6 different modes:
        - source=cells, interpolate=false:
            finds per surface face centre the \em nearest cell centre
            and uses its value
        - source=cells, interpolate=true:
            finds per surface face centre the \em nearest cell centre.
            Per surface point checks if this nearest cell is the one containing
            point; otherwise projects the point onto the nearest point on
            the boundary of the cell (to make sure interpolateCellPoint
            gets a valid location)

        - source=insideCells, interpolate=false:
            finds per surface face centre the cell containing it
            and uses its value.
            Trims surface faces outside of the mesh.
        - source=insideCells, interpolate=true:
            Per surface point interpolate cell containing it.

        - source=boundaryFaces, interpolate=false:
            finds per surface face centre the \em nearest point on the boundary
            (uncoupled faces only) and uses the value (or 0 if the nearest
            is on an empty boundary)
        - source=boundaryFaces, interpolate=true:
            finds per surface face centre the \em nearest point on the boundary
            (uncoupled faces only).
            Per surface point projects the point onto this boundary face
            (to make sure interpolateCellPoint gets a valid location)

    - since it finds the nearest per surface face, each surface face
      is guaranteed to be on one processor only.
      So after stitching the original surface should be complete.

    This is often embedded as part of a sampled surfaces function object.

Usage
    Example of function object partial specification:
    \verbatim
    surfaces
    {
        surface1
        {
            type    meshedSurface;
            surface something.obj;
            source  cells;

            //- Max sampling distance
            maxDistance 0.005;

            //- Fallback for missed sampling in 'cells' mode
            defaultValue
            {
                "p.*"   1e5;
                T       273.15;
            }
        }
    }
    \endverbatim

    Where the sub-entries comprise:
    \table
        Property | Description                              | Required | Default
        type     | meshedSurface                            | yes |
        surface  | surface name in triSurface/              | yes |
        patches  | Limit to named surface regions (wordRes) | no  |
        source   | cells/insideCells/boundaryFaces          | yes |
        keepIds  | pass through id numbering                | no  | true
        file     | Alternative file name                    | no  |
        fileType | The surface format                       | no  | (extension)
        scale    | Surface scaling factor                   | no  | 0
        maxDistance | Max search distance                   | no  | GREAT
        defaultValue | Value beyond max distance (dictionary) | no  | empty
    \endtable


SourceFiles
    sampledMeshedSurface.C
    sampledMeshedSurfaceTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef sampledMeshedSurface_H
#define sampledMeshedSurface_H

#include "sampledSurface.H"
#include "MeshedSurfaces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class meshSearch;

/*---------------------------------------------------------------------------*\
                       Class sampledMeshedSurface Declaration
\*---------------------------------------------------------------------------*/

class sampledMeshedSurface
:
    public sampledSurface,
    public meshedSurface
{
public:

        //- Types of sampling regions
        enum class samplingSource : unsigned char
        {
            cells,          //!< Use nearest cell value
            insideCells,    //!< Surface face within a cell, or trim
            boundaryFaces   //!< Use nearest boundary face values
        };

private:

    //- The concrete storage type for faces and points
    typedef meshedSurface Mesh;


    // Private Data

        static const Enum<samplingSource> samplingSourceNames_;

        //- The name of the input surface
        word surfaceName_;

        //- Surface to sample on
        meshedSurface surface_;

        //- Whether to sample internal cell values or boundary values
        const samplingSource sampleSource_;

        //- Track if the surface needs an update
        mutable bool needsUpdate_;

        //- Retain element ids/order of original surface
        bool keepIds_;

        //- For compatibility with the meshSurf interface
        labelList zoneIds_;

        //- From local surface triangle to mesh cell/face.
        labelList sampleElements_;

        //- Local points to sample per point
        pointField samplePoints_;

        //- Max search distance squared
        scalar maxDistanceSqr_;

        //- Out-of-range fallback values (when beyond the max search distance)
        dictionary defaultValues_;


    // Private Member Functions

        //- Set zoneIds list based on the surfZone information
        void setZoneMap();

        //- Sample volume field onto surface faces
        template<class Type>
        tmp<Field<Type>> sampleOnFaces
        (
            const interpolation<Type>& sampler
        ) const;

        //- Interpolate volume field onto surface points
        template<class Type>
        tmp<Field<Type>> sampleOnPoints
        (
            const interpolation<Type>& interpolator
        ) const;

        bool update(const meshSearch& meshSearcher);

public:

    //- Declare type-name, virtual type (with debug switch)
    TypeName("sampledMeshedSurface");


    // Constructors

        //- Construct from components
        sampledMeshedSurface
        (
            const word& name,
            const polyMesh& mesh,
            const word& surfaceName,
            const samplingSource sampleSource
        );

        //- Construct from dictionary
        sampledMeshedSurface
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~sampledMeshedSurface() = default;


    // Member Functions

        //- Does the surface need an update?
        virtual bool needsUpdate() const;

        //- Mark the surface as needing an update.
        //  May also free up unneeded data.
        //  Return false if surface was already marked as expired.
        virtual bool expire();

        //- Update the surface as required.
        //  Do nothing (and return false) if no update was needed
        virtual bool update();

        //- Update the surface using a bound box to limit the searching.
        //  For direct use, i.e. not through sample.
        //  Do nothing (and return false) if no update was needed
        bool update(const treeBoundBox& bb);

        //- Points of surface
        virtual const pointField& points() const
        {
            return Mesh::points();
        }

        //- Faces of surface
        virtual const faceList& faces() const
        {
            return Mesh::surfFaces();
        }

        //- Per-face zone/region information
        virtual const labelList& zoneIds() const
        {
            return zoneIds_;
        }

        //- Per-face identifier (eg, element Id)
        virtual const labelList& faceIds() const
        {
            return Mesh::faceIds();
        }

        //- Face area vectors
        virtual const vectorField& Sf() const
        {
            return Mesh::Sf();
        }

        //- Face area magnitudes
        virtual const scalarField& magSf() const
        {
            return Mesh::magSf();
        }

        //- Face centres
        virtual const vectorField& Cf() const
        {
            return Mesh::Cf();
        }

        //- If element ids/order of the original surface are kept
        virtual bool hasFaceIds() const
        {
            return keepIds_ && !Mesh::faceIds().empty();
        }

        //- Sampling boundary values instead of cell values
        bool onBoundary() const
        {
            return sampleSource_ == samplingSource::boundaryFaces;
        }


    // Sample

        //- Sample volume field onto surface faces
        virtual tmp<scalarField> sample
        (
            const interpolation<scalar>& sampler
        ) const;

        //- Sample volume field onto surface faces
        virtual tmp<vectorField> sample
        (
            const interpolation<vector>& sampler
        ) const;

        //- Sample volume field onto surface faces
        virtual tmp<sphericalTensorField> sample
        (
            const interpolation<sphericalTensor>& sampler
        ) const;

        //- Sample volume field onto surface faces
        virtual tmp<symmTensorField> sample
        (
            const interpolation<symmTensor>& sampler
        ) const;

        //- Sample volume field onto surface faces
        virtual tmp<tensorField> sample
        (
            const interpolation<tensor>& sampler
        ) const;


    // Interpolate

        //- Interpolate volume field onto surface points
        virtual tmp<scalarField> interpolate
        (
            const interpolation<scalar>& interpolator
        ) const;

        //- Interpolate volume field onto surface points
        virtual tmp<vectorField> interpolate
        (
            const interpolation<vector>& interpolator
        ) const;

        //- Interpolate volume field onto surface points
        virtual tmp<sphericalTensorField> interpolate
        (
            const interpolation<sphericalTensor>& interpolator
        ) const;

        //- Interpolate volume field onto surface points
        virtual tmp<symmTensorField> interpolate
        (
            const interpolation<symmTensor>& interpolator
        ) const;

        //- Interpolate volume field onto surface points
        virtual tmp<tensorField> interpolate
        (
            const interpolation<tensor>& interpolator
        ) const;


    // Output

        //- Print information
        virtual void print(Ostream& os, int level=0) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "sampledMeshedSurfaceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
