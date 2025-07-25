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
    Foam::sampledThresholdCellFaces

Description
    A sampledSurface defined by the cell faces corresponding to a threshold
    value.

    This is often embedded as part of a sampled surfaces function object.

Usage
    Example of function object partial specification:
    \verbatim
    surfaces
    {
        surface1
        {
            type    thresholdCellFaces;
            field   rho;
            lowerLimit 0.1;
        }
    }
    \endverbatim

    Where the sub-entries comprise:
    \table
        Property | Description                             | Required | Default
        type     | thresholdCellFaces                      | yes      |
        field    | field name for threshold                | yes      |
        lowerLimit | lower limit for threshold             | partly   | -Inf
        upperLimit | upper limit for threshold             | partly   | +Inf
        triangulate | triangulate faces                    | no       | false
    \endtable

Note
    Must specify at least one or both of \c lowerLimit or \c upperLimit

SeeAlso
    Foam::thresholdCellFaces

SourceFiles
    sampledThresholdCellFaces.C

\*---------------------------------------------------------------------------*/

#ifndef sampledThresholdCellFaces_H
#define sampledThresholdCellFaces_H

#include "sampledSurface.H"
#include "MeshedSurface.H"
#include "MeshedSurfacesFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class sampledThresholdCellFaces Declaration
\*---------------------------------------------------------------------------*/

class sampledThresholdCellFaces
:
    public sampledSurface,
    public meshedSurface
{
        // Private typedef for convenience
        typedef meshedSurface Mesh;


    // Private Data

        //- Threshold field
        const word fieldName_;

        //- Threshold value
        const scalar lowerThreshold_;

        //- Threshold value
        const scalar upperThreshold_;

        //- Triangulated faces or keep faces as is
        bool triangulate_;

        // Recreated for every time-step

        //- Time at last call, also track it surface needs an update
        mutable label prevTimeIndex_;

        //- For every face the original cell in mesh
        mutable labelList meshCells_;


    // Private Member Functions

        //- Create surface (if time has changed)
        //  Do nothing (and return false) if no update was needed
        bool updateGeometry() const;

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
    TypeName("sampledThresholdCellFaces");


    // Constructors

        //- Construct from dictionary
        sampledThresholdCellFaces
        (
            const word& name,
            const polyMesh&,
            const dictionary&
        );


    //- Destructor
    virtual ~sampledThresholdCellFaces() = default;


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
            return labelList::null();
        }

        //- Face area vectors (normals)
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
    #include "sampledThresholdCellFacesTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
