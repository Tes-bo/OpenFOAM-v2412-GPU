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
    Foam::sampledPlane

Description
    A sampledSurface defined by a plane which \a cuts the mesh using the
    cuttingPlane alorithm.  The surface is triangulated by default.

    This is often embedded as part of a sampled surfaces function object.

Usage
    Example of function object partial specification:
    \verbatim
    surfaces
    {
        surface1
        {
            type        plane;
            point       ...;
            normal      ...;
        }
    }
    \endverbatim

    Where the sub-entries comprise:
    \table
        Property | Description                              | Reqd | Default
        type     | plane                                    | yes  |
        planeType | Plane description (pointAndNormal etc)  | no   |
        triangulate | triangulate faces                     | no   | true
        bounds   | limit with bounding box                  | no   |
        zone     | limit to cell zone (name or regex)       | no   |
        zones    | limit to cell zones (names, regexs)      | no   |
        coordinateSystem | define plane within given cartesian system | no |
        transform | define plane within given cartesian system | no |
    \endtable

SeeAlso
    Foam::coordinateSystem
    Foam::cuttingPlane
    Foam::plane

Note
    Does not actually cut until update() called.
    The keyword \c zones has priority over \c zone.

SourceFiles
    sampledPlane.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_sampledPlane_H
#define Foam_sampledPlane_H

#include "sampledSurface.H"
#include "cuttingPlane.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class sampledPlane Declaration
\*---------------------------------------------------------------------------*/

class sampledPlane
:
    public sampledSurface,
    public cuttingPlane
{
    // Private data

        //- The zone or zones in which cutting is to occur
        wordRes zoneNames_;

        //- Optional bounding box to trim against
        const boundBox bounds_;

        //- Triangulate faces or not
        const bool triangulate_;

        //- Track if the surface needs an update
        mutable bool needsUpdate_;


    // Private Member Functions

        //- Define plane from dictionary entry and apply any coordinate
        //- transformations
        static plane definePlane(const polyMesh& mesh, const dictionary& dict);

        //- Define cell selection from zones and bounding box.
        //  Optionally check and warn if the plane does not intersect
        //  with the bounds of the mesh (or submesh) or if the bounding box
        //  does not overlap with the mesh (or submesh)
        bitSet cellSelection(const bool warn=false) const;


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


public:

    //- Runtime type information
    TypeName("sampledPlane");


    // Constructors

        //- Construct from components
        sampledPlane
        (
            const word& name,
            const polyMesh& mesh,
            const plane& planeDesc,
            const wordRes& zones = wordRes(),
            const bool triangulate = true
        );

        //- Construct from dictionary
        sampledPlane
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~sampledPlane() = default;


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
            return cuttingPlane::points();
        }

        //- Faces of surface
        virtual const faceList& faces() const
        {
            return cuttingPlane::surfFaces();
        }

        //- Per-face zone/region information
        //  Could instead return meshCells or cellZoneId of the meshCells.
        virtual const labelList& zoneIds() const
        {
            return labelList::null();
        }

        //- Face area magnitudes
        virtual const vectorField& Sf() const
        {
            return cuttingPlane::Sf();
        }

        //- Face area magnitudes
        virtual const scalarField& magSf() const
        {
            return cuttingPlane::magSf();
        }

        //- Face centres
        virtual const vectorField& Cf() const
        {
            return cuttingPlane::Cf();
        }

        //- For each face, the original cell in mesh
        using cuttingPlane::meshCells;


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
    #include "sampledPlaneTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
