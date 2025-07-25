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

\*---------------------------------------------------------------------------*/

#include "sampledMeshedSurfaceNormal.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sampledMeshedSurfaceNormal, 0);
    // Use shorter name only
    addNamedToRunTimeSelectionTable
    (
        sampledSurface,
        sampledMeshedSurfaceNormal,
        word,
        meshedSurfaceNormal
    );
    // Compatibility name (1912)
    addNamedToRunTimeSelectionTable
    (
        sampledSurface,
        sampledMeshedSurfaceNormal,
        word,
        sampledTriSurfaceMeshNormal
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledMeshedSurfaceNormal::sampledMeshedSurfaceNormal
(
    const word& name,
    const polyMesh& mesh,
    const word& surfaceName,
    const samplingSource sampleSource
)
:
    sampledMeshedSurface(name, mesh, surfaceName, sampleSource)
{}


Foam::sampledMeshedSurfaceNormal::sampledMeshedSurfaceNormal
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    sampledMeshedSurface(name, mesh, dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::Field<Foam::vector>>
Foam::sampledMeshedSurfaceNormal::sample
(
    const interpolation<vector>& sampler
) const
{
    auto tvalues = tmp<Field<vector>>::New(size(), Zero);

    tvalues.ref().replace
    (
        0,
        meshedSurface::faceNormals()
       &sampledMeshedSurface::sample(sampler)
    );

    return tvalues;
}


Foam::tmp<Foam::Field<Foam::vector>>
Foam::sampledMeshedSurfaceNormal::interpolate
(
    const interpolation<vector>& interpolator
) const
{
    auto tvalues = tmp<Field<vector>>::New(points().size(), Zero);

    pointField allNormals(points().size(), Zero);
    UIndirectList<vector>(allNormals, meshPoints()) = pointNormals();

    tvalues.ref().replace
    (
        0,
        allNormals
       &sampledMeshedSurface::interpolate(interpolator)
    );

    return tvalues;
}


// ************************************************************************* //
