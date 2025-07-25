/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2020 OpenCFD Ltd.
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
    Foam::sampledMeshedSurfaceNormal

Description
    Variant of sampledMeshedSurface that samples the surface-normal component
    of a vector field.

    Returns a vector field with the value in the first component and sets
    the other two to zero.

    This is often embedded as part of a sampled surfaces function object.

Usage
    Example of function object partial specification:
    \verbatim
    surfaces
    {
        surface1
        {
            type    sampledMeshedSurfaceNormal;
            surface something.obj;
            source  cells;
        }
    }
    \endverbatim

    Where the sub-entries comprise:
    \table
        Property | Description                             | Required | Default
        type     | meshedSurfaceNormal                     | yes      |
        surface  | surface name in triSurface/             | yes      |
        patches  | Limit to named surface regions (wordRes) | no  |
        source   | cells/insideCells/boundaryFaces         | yes      |
        keepIds  | pass through id numbering               | no       | false
    \endtable

SeeAlso
    Foam::sampledMeshedSurface

SourceFiles
    sampledMeshedSurfaceNormal.C
    sampledMeshedSurfaceNormalTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef sampledMeshedSurfaceNormal_H
#define sampledMeshedSurfaceNormal_H

#include "sampledMeshedSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class sampledMeshedSurfaceNormal Declaration
\*---------------------------------------------------------------------------*/

class sampledMeshedSurfaceNormal
:
    public sampledMeshedSurface
{
public:

    //- Runtime type information
    TypeName("sampledMeshedSurfaceNormal");


    // Constructors

        //- Construct from components
        sampledMeshedSurfaceNormal
        (
            const word& name,
            const polyMesh& mesh,
            const word& surfaceName,
            const samplingSource sampleSource
        );

        //- Construct from dictionary
        sampledMeshedSurfaceNormal
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~sampledMeshedSurfaceNormal() = default;


    // Member Functions

    // Sample

        //- Sample volume field onto surface faces
        virtual tmp<scalarField> sample
        (
            const interpolation<scalar>& sampler
        ) const
        {
            NotImplemented;
            return nullptr;
        }

        //- Sample volume field onto surface faces
        virtual tmp<vectorField> sample
        (
            const interpolation<vector>& sampler
        ) const;

        //- Sample volume field onto surface faces
        virtual tmp<sphericalTensorField> sample
        (
            const interpolation<sphericalTensor>& sampler
        ) const
        {
            NotImplemented;
            return nullptr;
        }

        //- Sample volume field onto surface faces
        virtual tmp<symmTensorField> sample
        (
            const interpolation<symmTensor>& sampler
        ) const
        {
            NotImplemented;
            return nullptr;
        }

        //- Sample volume field onto surface faces
        virtual tmp<tensorField> sample
        (
            const interpolation<tensor>& sampler
        ) const
        {
            NotImplemented;
            return nullptr;
        }


    // Interpolate

        //- Interpolate volume field onto surface points
        virtual tmp<scalarField> interpolate
        (
            const interpolation<scalar>& interpolator
        ) const
        {
            NotImplemented;
            return nullptr;
        }

        //- Interpolate volume field onto surface points
        virtual tmp<vectorField> interpolate
        (
            const interpolation<vector>& interpolator
        ) const;

        //- Interpolate volume field onto surface points
        virtual tmp<sphericalTensorField> interpolate
        (
            const interpolation<sphericalTensor>& interpolator
        ) const
        {
            NotImplemented;
            return nullptr;
        }

        //- Interpolate volume field onto surface points
        virtual tmp<symmTensorField> interpolate
        (
            const interpolation<symmTensor>& interpolator
        ) const
        {
            NotImplemented;
            return nullptr;
        }

        //- Interpolate volume field onto surface points
        virtual tmp<tensorField> interpolate
        (
            const interpolation<tensor>& interpolator
        ) const
        {
            NotImplemented;
            return nullptr;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
