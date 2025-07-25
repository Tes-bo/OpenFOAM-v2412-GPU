/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2022 OpenCFD Ltd.
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

#include "sampledIsoSurface.H"
#include "volFieldsFwd.H"
#include "pointFields.H"
#include "volPointInterpolation.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::sampledIsoSurface::sampleOnFaces
(
    const interpolation<Type>& sampler
) const
{
    updateGeometry();  // Recreate geometry if time has changed

    return sampledSurface::sampleOnFaces
    (
        sampler,
        meshCells(),
        surface(),
        points()
    );
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::sampledIsoSurface::sampleOnPoints
(
    const interpolation<Type>& interpolator
) const
{
    updateGeometry();  // Recreate geometry if time has changed

    if (isoSurfacePtr_)
    {
        return this->sampleOnIsoSurfacePoints(interpolator);
    }

    return sampledSurface::sampleOnPoints
    (
        interpolator,
        meshCells(),
        faces(),
        points()
    );
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::sampledIsoSurface::sampleOnIsoSurfacePoints
(
    const interpolation<Type>& interpolator
) const
{
    if (!isoSurfacePtr_)
    {
        FatalErrorInFunction
            << "cannot call without an iso-surface" << nl
            << exit(FatalError);
    }

    // Assume volPointInterpolation for the point field!
    const auto& volFld = interpolator.psi();

    tmp<VolumeField<Type>> tvolFld(volFld);
    tmp<PointField<Type>> tpointFld;

    if (subMeshPtr_)
    {
        // Replace with subset
        tvolFld.reset(subMeshPtr_->interpolate(volFld));
    }

    // Interpolated point field
    tpointFld.reset
    (
        volPointInterpolation::New(tvolFld().mesh()).interpolate(tvolFld())
    );

    if (average_)
    {
        tvolFld.reset(pointAverage(tpointFld()));
    }

    return isoSurfacePtr_->interpolate(tvolFld(), tpointFld());
}


// ************************************************************************* //
