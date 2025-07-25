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
    Foam::sampledNone

Description
    A no operation sampledSurface that can be used when a sampler
    is expected but is not desired. For example, to temporarily disable
    a sampling definition, or to provide a boilerplate definition that
    is overwritten at a later stage in a dictionary.

Usage
    \table
        Property | Description                             | Required | Default
        type     | 'none'                                  | yes      |
    \endtable

SourceFiles
    sampledNone.C

\*---------------------------------------------------------------------------*/

#ifndef sampledNone_H
#define sampledNone_H

#include "sampledSurface.H"
#include "MeshedSurfaces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class sampledNone Declaration
\*---------------------------------------------------------------------------*/

class sampledNone
:
    public meshedSurface,
    public sampledSurface
{
    typedef meshedSurface Mesh;

public:

    //- Runtime type information
    TypeName("sampledNone");


    // Constructors

        //- Construct null
        explicit sampledNone();

        //- Construct null, with specified name
        explicit sampledNone(const word& name);

        //- Construct null, with dictionary
        sampledNone(const word& name, const polyMesh&, const dictionary&);


    //- Destructor
    virtual ~sampledNone() = default;


    // Member Functions

        //- Does the surface need an update?
        virtual bool needsUpdate() const;

        //- Mark the surface as needing an update.
        virtual bool expire();

        //- Update the surface as required.
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

        //- Face area magnitudes
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
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
