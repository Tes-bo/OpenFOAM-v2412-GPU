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

#include "sampledCuttingSurface.H"
#include "dictionary.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sampledCuttingSurface, 0);
    addToRunTimeSelectionTable
    (
        sampledSurface,
        sampledCuttingSurface,
        word
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::bitSet Foam::sampledCuttingSurface::cellSelection(const bool warn) const
{
    boundBox meshBounds;

    bitSet cellsToSelect =
        cuttingSurface::cellSelection
        (
            mesh(), bounds_, zoneNames_, meshBounds
        );

    if (warn)
    {
        cuttingSurface::checkOverlap(name(), meshBounds, bounds_);
    }

    return cellsToSelect;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledCuttingSurface::sampledCuttingSurface
(
    const polyMesh& mesh,
    const word& surfaceType,
    const word& surfaceName,
    const bool triangulate,
    const boundBox& bounds
)
:
    sampledSurface(surfaceName, mesh),
    cuttingSurface(mesh, surfaceType, surfaceName),
    zoneNames_(),
    bounds_(bounds),
    triangulate_(triangulate),
    needsUpdate_(true)
{}


Foam::sampledCuttingSurface::sampledCuttingSurface
(
    const word& defaultSurfaceName,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    sampledSurface(defaultSurfaceName, mesh, dict),
    cuttingSurface(defaultSurfaceName, mesh, dict),
    zoneNames_(),
    bounds_(dict.getOrDefault("bounds", boundBox::null())),
    triangulate_(dict.getOrDefault("triangulate", true)),
    needsUpdate_(true)
{
    if (!dict.readIfPresent("zones", zoneNames_) && dict.found("zone"))
    {
        zoneNames_.resize(1);
        dict.readEntry("zone", zoneNames_.first());
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::sampledCuttingSurface::needsUpdate() const
{
    return needsUpdate_;
}


bool Foam::sampledCuttingSurface::expire()
{
    // Already marked as expired
    if (needsUpdate_)
    {
        return false;
    }

    sampledSurface::clearGeom();

    needsUpdate_ = true;
    return true;
}


bool Foam::sampledCuttingSurface::update()
{
    if (!needsUpdate_)
    {
        return false;
    }

    sampledSurface::clearGeom();

    performCut(mesh(), triangulate_, cellSelection(true));

    if (debug)
    {
        cuttingSurface::print(Pout, debug);
        Pout<< endl;
    }

    needsUpdate_ = false;
    return true;
}


Foam::tmp<Foam::scalarField> Foam::sampledCuttingSurface::sample
(
    const interpolation<scalar>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::vectorField> Foam::sampledCuttingSurface::sample
(
    const interpolation<vector>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::sphericalTensorField> Foam::sampledCuttingSurface::sample
(
    const interpolation<sphericalTensor>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::symmTensorField> Foam::sampledCuttingSurface::sample
(
    const interpolation<symmTensor>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::tensorField> Foam::sampledCuttingSurface::sample
(
    const interpolation<tensor>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::scalarField> Foam::sampledCuttingSurface::interpolate
(
    const interpolation<scalar>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


Foam::tmp<Foam::vectorField> Foam::sampledCuttingSurface::interpolate
(
    const interpolation<vector>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}

Foam::tmp<Foam::sphericalTensorField> Foam::sampledCuttingSurface::interpolate
(
    const interpolation<sphericalTensor>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


Foam::tmp<Foam::symmTensorField> Foam::sampledCuttingSurface::interpolate
(
    const interpolation<symmTensor>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


Foam::tmp<Foam::tensorField> Foam::sampledCuttingSurface::interpolate
(
    const interpolation<tensor>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


// ************************************************************************* //
