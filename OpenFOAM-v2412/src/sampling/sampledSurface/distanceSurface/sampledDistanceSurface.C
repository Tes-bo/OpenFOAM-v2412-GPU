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

\*---------------------------------------------------------------------------*/

#include "sampledDistanceSurface.H"
#include "dictionary.H"
#include "volFields.H"
#include "volPointInterpolation.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "volumeType.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sampledDistanceSurface, 0);
    addNamedToRunTimeSelectionTable
    (
        sampledSurface,
        sampledDistanceSurface,
        word,
        distanceSurface
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledDistanceSurface::sampledDistanceSurface
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    sampledSurface(name, mesh, dict),
    distanceSurface(name, mesh, dict),
    average_(dict.getOrDefault("average", false)),
    needsUpdate_(true)
{}


Foam::sampledDistanceSurface::sampledDistanceSurface
(
    const word name,
    const polyMesh& mesh,
    const bool interpolate,
    autoPtr<searchableSurface>&& surface,
    const scalar distance,
    const bool useSignedDistance
)
:
    sampledSurface(name, mesh, interpolate),
    distanceSurface
    (
        mesh,
        interpolate,
        std::move(surface),
        distance,
        useSignedDistance
    ),
    average_(false), // pass as param?
    needsUpdate_(true)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::sampledDistanceSurface::needsUpdate() const
{
    return needsUpdate_;
}


bool Foam::sampledDistanceSurface::expire()
{
    if (debug)
    {
        Pout<< "sampledDistanceSurface::expire :"
            << " needsUpdate:" << needsUpdate_ << endl;
    }

    // Clear derived data
    clearGeom();

    // already marked as expired
    if (needsUpdate_)
    {
        return false;
    }

    needsUpdate_ = true;
    return true;
}


bool Foam::sampledDistanceSurface::update()
{
    if (debug)
    {
        Pout<< "sampledDistanceSurface::update :"
            << " needsUpdate:" << needsUpdate_ << endl;
    }

    if (!needsUpdate_)
    {
        return false;
    }

    distanceSurface::createGeometry();

    needsUpdate_ = false;
    return true;
}


Foam::tmp<Foam::scalarField> Foam::sampledDistanceSurface::sample
(
    const interpolation<scalar>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::vectorField> Foam::sampledDistanceSurface::sample
(
    const interpolation<vector>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::sphericalTensorField> Foam::sampledDistanceSurface::sample
(
    const interpolation<sphericalTensor>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::symmTensorField> Foam::sampledDistanceSurface::sample
(
    const interpolation<symmTensor>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::tensorField> Foam::sampledDistanceSurface::sample
(
    const interpolation<tensor>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::scalarField> Foam::sampledDistanceSurface::interpolate
(
    const interpolation<scalar>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


Foam::tmp<Foam::vectorField> Foam::sampledDistanceSurface::interpolate
(
    const interpolation<vector>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}

Foam::tmp<Foam::sphericalTensorField> Foam::sampledDistanceSurface::interpolate
(
    const interpolation<sphericalTensor>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


Foam::tmp<Foam::symmTensorField> Foam::sampledDistanceSurface::interpolate
(
    const interpolation<symmTensor>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


Foam::tmp<Foam::tensorField> Foam::sampledDistanceSurface::interpolate
(
    const interpolation<tensor>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


void Foam::sampledDistanceSurface::print(Ostream& os, int level) const
{
    os  << "distanceSurface: " << name() << " :";
    distanceSurface::print(os, level);
}


// ************************************************************************* //
