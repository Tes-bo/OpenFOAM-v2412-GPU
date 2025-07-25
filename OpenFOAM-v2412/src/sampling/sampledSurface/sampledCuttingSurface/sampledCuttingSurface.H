/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2019 OpenCFD Ltd.
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
    Foam::sampledCuttingSurface

Description
    A surface define by using an input surface to cut the mesh cells

Usage
    Dictionary controls:
    \table
        Property | Description                             | Required | Default
        type     | surfaceCut                              | yes      |
        surfaceType | type of surface                      | yes      |
        surfaceName | name of surface in \c triSurface/    | no  | dict name
        triangulate | triangulate faces                    | no       | true
        bounds   | limit with bounding box                 | no       |
        zone     | limit to cell zone (name or regex)      | no       |
        zones    | limit to cell zones (names, regexs)     | no       |
    \endtable

Note
    No attempt at resolving degenerate cases.
    Since the cut faces can be quite ugly, they will often be triangulated.

SourceFiles
    sampledCuttingSurface.C
    sampledCuttingSurfaceTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef sampledCuttingSurface_H
#define sampledCuttingSurface_H

#include "sampledSurface.H"
#include "cuttingSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class sampledCuttingSurface Declaration
\*---------------------------------------------------------------------------*/

class sampledCuttingSurface
:
    public sampledSurface,
    public cuttingSurface
{
    // Private Data

        //- The zone or zones in which cutting is to occur
        wordRes zoneNames_;

        //- Optional bounding box to trim against
        const boundBox bounds_;

        //- Triangulate faces or not
        const bool triangulate_;

        //- Track if the surface needs an update
        mutable bool needsUpdate_;


    // Private Member Functions

        //- Define cell selection from zones and bounding box.
        //  Optionally check and warn if the bounding box
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
    TypeName("surfaceCut");


    // Constructors

        //- Construct from components
        sampledCuttingSurface
        (
            const polyMesh& mesh,
            const word& surfaceType,
            const word& surfaceName,
            const bool triangulate = true,
            const boundBox& bounds = boundBox::null()
        );

        //- Construct from dictionary
        sampledCuttingSurface
        (
            const word& defaultSurfaceName,
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~sampledCuttingSurface() = default;


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
            return cuttingSurface::points();
        }

        //- Faces of surface
        virtual const faceList& faces() const
        {
            return cuttingSurface::surfFaces();
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
            return cuttingSurface::Sf();
        }

        //- Face area magnitudes
        virtual const scalarField& magSf() const
        {
            return cuttingSurface::magSf();
        }

        //- Face centres
        virtual const vectorField& Cf() const
        {
            return cuttingSurface::Cf();
        }

        //- For each face, the original cell in mesh
        using cuttingSurface::meshCells;


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
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "sampledCuttingSurfaceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
