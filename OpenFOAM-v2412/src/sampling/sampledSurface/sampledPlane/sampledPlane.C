/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2022 OpenCFD Ltd.
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

#include "sampledPlane.H"
#include "dictionary.H"
#include "polyMesh.H"
#include "volFields.H"
#include "cartesianCS.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sampledPlane, 0);
    addNamedToRunTimeSelectionTable
    (
        sampledSurface,
        sampledPlane,
        word,
        plane
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::plane Foam::sampledPlane::definePlane
(
    const polyMesh& mesh,
    const dictionary& dict
)
{
    plane pln(dict);

    // Optional (cartesian) coordinate transform.
    // - with registry to allow lookup from globally defined systems

    auto csysPtr = coordinateSystem::NewIfPresent(mesh, dict);

    if (!csysPtr)
    {
        csysPtr = coordinateSystem::NewIfPresent(dict, "transform");
    }

    // Make plane relative to the Cartesian coordinate system
    if (csysPtr)
    {
        coordSystem::cartesian cs(csysPtr());

        const point  orig = cs.globalPosition(pln.origin());
        const vector norm = cs.globalVector(pln.normal());

        DebugInfo
            << "plane "
            << " origin:" << pln.origin()
            << " normal:" << pln.normal()
            << " =>"
            << " origin:" << orig
            << " normal:" << norm
            << endl;

        // Reassign the plane
        pln = plane(orig, norm);
    }

    return pln;
}


Foam::bitSet Foam::sampledPlane::cellSelection(const bool warn) const
{
    return cuttingPlane::cellSelection
    (
        mesh(),
        bounds_,
        zoneNames_,
        name(),
        warn
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledPlane::sampledPlane
(
    const word& name,
    const polyMesh& mesh,
    const plane& planeDesc,
    const wordRes& zones,
    const bool triangulate
)
:
    sampledSurface(name, mesh),
    cuttingPlane(planeDesc),
    zoneNames_(zones),
    bounds_(),
    triangulate_(triangulate),
    needsUpdate_(true)
{
    if (debug)
    {
        if (!zoneNames_.empty())
        {
            Info<< " cellZones " << flatOutput(zoneNames_);

            if (-1 == mesh.cellZones().findIndex(zoneNames_))
            {
                Info<< " not found!";
            }
            Info<< endl;
        }
    }
}


Foam::sampledPlane::sampledPlane
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    sampledSurface(name, mesh, dict),
    cuttingPlane(definePlane(mesh, dict)),
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

    if (debug)
    {
        Info<< "plane " << name << " :"
            << " origin:" << origin()
            << " normal:" << normal();

        if (bounds_.good())
        {
            Info<< " bounds:" << bounds_;
        }

        if (!zoneNames_.empty())
        {
            Info<< " cellZones " << flatOutput(zoneNames_);

            if (-1 == mesh.cellZones().findIndex(zoneNames_))
            {
                Info<< " not found!";
            }
        }
        Info<< endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::sampledPlane::needsUpdate() const
{
    return needsUpdate_;
}


bool Foam::sampledPlane::expire()
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


bool Foam::sampledPlane::update()
{
    if (!needsUpdate_)
    {
        return false;
    }

    sampledSurface::clearGeom();

    performCut(mesh(), triangulate_, cellSelection(true));

    if (debug)
    {
        print(Pout, debug);
        Pout<< endl;
    }

    needsUpdate_ = false;
    return true;
}


Foam::tmp<Foam::scalarField> Foam::sampledPlane::sample
(
    const interpolation<scalar>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::vectorField> Foam::sampledPlane::sample
(
    const interpolation<vector>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::sphericalTensorField> Foam::sampledPlane::sample
(
    const interpolation<sphericalTensor>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::symmTensorField> Foam::sampledPlane::sample
(
    const interpolation<symmTensor>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::tensorField> Foam::sampledPlane::sample
(
    const interpolation<tensor>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::scalarField> Foam::sampledPlane::interpolate
(
    const interpolation<scalar>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


Foam::tmp<Foam::vectorField> Foam::sampledPlane::interpolate
(
    const interpolation<vector>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}

Foam::tmp<Foam::sphericalTensorField> Foam::sampledPlane::interpolate
(
    const interpolation<sphericalTensor>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


Foam::tmp<Foam::symmTensorField> Foam::sampledPlane::interpolate
(
    const interpolation<symmTensor>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


Foam::tmp<Foam::tensorField> Foam::sampledPlane::interpolate
(
    const interpolation<tensor>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


void Foam::sampledPlane::print(Ostream& os, int level) const
{
    os  << "sampledPlane: " << name() << " :"
        << " origin:" << plane::origin()
        << " normal:" << plane::normal()
        << " triangulate:" << triangulate_;

    if (level)
    {
        os << " faces:" << faces().size()
           << " points:" << points().size();
    }
}


// ************************************************************************* //
