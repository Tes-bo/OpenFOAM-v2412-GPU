/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2021 OpenCFD Ltd.
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
    Foam::sampledDistanceSurface

Description
    A sampledSurface defined by a distance to a surface - resolved using
    an iso-surface (algorithm: cell, point, topo)

    This is often embedded as part of a sampled surfaces function object.

Usage
    Example of function object partial specification:
    \verbatim
    surfaces
    {
        surface1
        {
            type        distanceSurface;
            surfaceType triSurfaceMesh;
            surfaceName something.obj;
            topology    proximity;
        }
    }
    \endverbatim

    Where the sub-entries comprise:
    \table
        Property | Description                              | Required | Default
        type     | distanceSurface                          | yes |
        distance | distance from surface                    | no  | 0
        signed   | Use sign when distance is positive       | no  | true
        isoMethod | Iso-algorithm (cell/topo/point)         | no  | default
        regularise | Face simplification (enum or bool)     | no  | true
        average  | Cell values from averaged point values   | no  | false
        bounds   | Limit with bounding box                  | no  |
        surfaceType | Type of surface                       | yes |
        surfaceName | Name of surface in \c triSurface/     | no  | dict name
        topology    | Topology filter name                  | no  | none
        nearestPoints | Points for point-based segmentation | no  |
        maxDistance | Max search distance for nearestPoints | no  | GREAT
        absProximity | Max proximity of face centres        | no  | 1e-5
    \endtable

SourceFiles
    sampledDistanceSurface.C
    sampledDistanceSurfaceTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef sampledDistanceSurface_H
#define sampledDistanceSurface_H

#include "sampledSurface.H"
#include "distanceSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class sampledDistanceSurface Declaration
\*---------------------------------------------------------------------------*/

class sampledDistanceSurface
:
    public sampledSurface,
    public distanceSurface
{
    // Private Data

        //- Whether to recalculate cell values as average of point values
        bool average_;

        //- Track if the surface needs an update
        mutable bool needsUpdate_;


    // Private Member Functions

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
            const interpolation<Type>& interpolation
        ) const;

        //- Use distance surface isoSurfacePtr_ for point interpolation
        template<class Type>
        tmp<Field<Type>> sampleOnIsoSurfacePoints
        (
            const interpolation<Type>& interpolator
        ) const;


public:

    //- Runtime type information
    TypeName("sampledDistanceSurface");


    // Constructors

        //- Construct from dictionary
        sampledDistanceSurface
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );


        sampledDistanceSurface
        (
            const word name,
            const polyMesh& mesh,
            const bool interpolate,
            autoPtr<searchableSurface>&& surface,
            const scalar distance,
            const bool useSignedDistance
            //const isoSurfaceParams& params = isoSurfaceParams()
        );


    //- Destructor
    virtual ~sampledDistanceSurface() = default;


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

        //- Face area vectors
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
    #include "sampledDistanceSurfaceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
