/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2022 OpenCFD Ltd.
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
    Foam::sampledCuttingPlane

Description
    A sampledSurface defined by a plane using an iso-surface algorithm
    to \a cut the mesh.

    This is often embedded as part of a sampled surfaces function object.

Usage
    Example of function object partial specification:
    \verbatim
    surfaces
    {
        surface1
        {
            type    cuttingPlane;
            point   ...;
            normal  ...;
        }
    }
    \endverbatim

    Where the sub-entries comprise:
    \table
        Property | Description                              | Required | Default
        type     | cuttingPlane                             | yes |
        planeType | Plane description (pointAndNormal etc)  | no  |
        offsets   | Offsets of the origin in the normal direction | no | (0)
        isoMethod | Iso-algorithm (cell/topo/point)         | no  | topo
        bounds   | limit with bounding box                  | no  |
        zone     | limit to cell zone (name or regex)       | no  |
        zones    | limit to cell zones (names, regexs)      | no  |
        exposedPatchName | name for zone subset             | optional |
        regularise | Face simplification (enum or bool)     | no  | true
        mergeTol | tolerance for merging points             | no  | 1e-6
        coordinateSystem | define plane within given cartesian system | no |
        transform | define plane within given cartesian system | no |
    \endtable

Note
    The keyword \c zones has priority over \c zone.

SeeAlso
    Foam::plane

SourceFiles
    sampledCuttingPlane.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_sampledCuttingPlane_H
#define Foam_sampledCuttingPlane_H

#include "sampledSurface.H"
#include "plane.H"
#include "fvMeshSubset.H"
#include "isoSurfaceBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class sampledCuttingPlane Declaration
\*---------------------------------------------------------------------------*/

class sampledCuttingPlane
:
    public sampledSurface
{
    // Private Data

        //- Plane
        const plane plane_;

        //- The offsets to the plane - defaults to (0).
        List<scalar> offsets_;

        //- Parameters (filtering etc) for iso-surface
        isoSurfaceParams isoParams_;

        //- Whether to recalculate cell values as average of point values
        bool average_;

        //- Use simple sub-meshing in algorithm itself
        bool simpleSubMesh_;

        //- The zone or zones in which cutting is to occur
        wordRes zoneNames_;

        //- For zones: patch to put exposed faces into
        mutable word exposedPatchName_;

        //- Track if the surface needs an update
        mutable bool needsUpdate_;


    // Sampling geometry. Directly stored or via an iso-surface (ALGO_POINT)

        //- The extracted surface (direct storage)
        mutable meshedSurface surface_;

        //- For every face the original cell in mesh (direct storage)
        mutable labelList meshCells_;

        //- Constructed iso-surface (ALGO_POINT), for interpolators
        autoPtr<isoSurfaceBase> isoSurfacePtr_;


    // Mesh Subsetting

        //- Cached subMesh for (pre-)subset of cell zones
        mutable autoPtr<fvMeshSubset> subMeshPtr_;

        //- Cached ignore cells for (post-)subset of cell zones
        mutable autoPtr<bitSet> ignoreCellsPtr_;


    // Fields

        //- Distance to cell centres
        autoPtr<volScalarField> cellDistancePtr_;

        //- Distance to points
        scalarField pointDistance_;


    // Private Member Functions

        //- Define plane from dictionary entry and apply any coordinate
        //- transformations
        static plane definePlane(const polyMesh& mesh, const dictionary& dict);

        //- Check and warn if bounding box does not intersect mesh or plane
        void checkBoundsIntersection
        (
            const plane& pln,
            const boundBox& meshBb
        ) const;

        //- Fill cellDistance, pointDistance fields for the specified plane
        void setDistanceFields(const plane& pln);

        //- Collect iso-surfaces into a single surface (No point merging)
        void combineSurfaces(PtrList<isoSurfaceBase>& isoSurfPtrs);

        //- Create iso surface
        void createGeometry();

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

        //- Use isoSurfacePtr_ for point interpolation
        template<class Type>
        tmp<Field<Type>> sampleOnIsoSurfacePoints
        (
            const interpolation<Type>& interpolator
        ) const;


protected:

    // Protected Member Functions

        //- Is currently backed by an isoSurfacePtr_
        bool hasIsoSurface() const
        {
            return bool(isoSurfacePtr_);
        }

public:

    //- Runtime type information
    TypeName("sampledCuttingPlane");


    // Constructors

        //- Construct from dictionary
        sampledCuttingPlane
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~sampledCuttingPlane() = default;


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

        //- The current surface geometry
        const meshedSurface& surface() const
        {
            if (isoSurfacePtr_)
            {
                return *isoSurfacePtr_;
            }
            return surface_;
        }

        //- For each face, the original cell in mesh
        const labelList& meshCells() const
        {
            if (isoSurfacePtr_)
            {
                return isoSurfacePtr_->meshCells();
            }
            return meshCells_;
        }

        //- Points of surface
        virtual const pointField& points() const
        {
            return surface().points();
        }

        //- Faces of surface
        virtual const faceList& faces() const
        {
            return surface().surfFaces();
        }

        //- Per-face zone/region information
        virtual const labelList& zoneIds() const
        {
            return labelList::null();
        }

        //- Face area magnitudes
        virtual const vectorField& Sf() const
        {
            return surface().Sf();
        }

        //- Face area magnitudes
        virtual const scalarField& magSf() const
        {
            return surface().magSf();
        }

        //- Face centres
        virtual const vectorField& Cf() const
        {
            return surface().Cf();
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
    #include "sampledCuttingPlaneTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
