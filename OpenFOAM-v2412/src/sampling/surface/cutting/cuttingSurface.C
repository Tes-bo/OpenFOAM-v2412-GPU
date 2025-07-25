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

#include "cuttingSurface.H"
#include "dictionary.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cuttingSurface::cuttingSurface
(
    const polyMesh& mesh,
    const word& surfaceType,
    const word& surfaceName
)
:
    cuttingSurfaceBase(),
    surfPtr_
    (
        searchableSurface::New
        (
            surfaceType,
            IOobject
            (
                surfaceName,            // name
                mesh.time().constant(), // directory
                "triSurface",           // instance
                mesh.time(),            // registry
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            dictionary()
        )
    )
{}


Foam::cuttingSurface::cuttingSurface
(
    const word& defaultSurfaceName,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    cuttingSurfaceBase(),
    surfPtr_
    (
        searchableSurface::New
        (
            dict.get<word>("surfaceType"),
            IOobject
            (
                dict.getOrDefault("surfaceName", defaultSurfaceName),
                mesh.time().constant(), // directory
                "triSurface",           // instance
                mesh.time(),            // registry
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            dict
        )
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cuttingSurface::performCut
(
    const primitiveMesh& mesh,
    const bool triangulate,
    bitSet&& cellIdLabels
)
{
    const fvMesh& fvm = static_cast<const fvMesh&>(mesh);

    Mesh::clear();
    meshCells_.clear();

    // Pre-populate with restriction
    bitSet cellCuts(std::move(cellIdLabels));

    if (cellCuts.size())
    {
        cellCuts.resize(mesh.nCells());
    }

    scalarField pointDist;
    calcCellCuts(fvm, pointDist, cellCuts);


    // Walk cell cuts to create faces

    // Action #1:
    // - Orient edge so it points in the positive gradient direction.
    // - Edge intersection when it spans across point-distance == 0.
    const auto edgeOrientIntersect =
        [=](edge& e) -> bool
        {
            if (pointDist[e.last()] < pointDist[e.first()])
            {
                e.flip();
            }

            const scalar s0 = pointDist[e.first()];
            const scalar s1 = pointDist[e.last()];

            if
            (
                s0 > ROOTVSMALL         // Edge is all positive
             || s1 < ROOTVSMALL         // Edge is all negative
             || Foam::mag(s1-s0) < ROOTVSMALL
            )
            {
                return false;
            }

            return true;
        };

    // Action #2:
    // - The edge intersection alpha for point-distance == 0
    // Return -1 for error.
    // This is like the iso-fraction for an iso-surface.
    const auto edgeAlphaIntersect =
        [=](const edge& e) -> scalar
        {
            const scalar s0 = pointDist[e.first()];
            const scalar s1 = pointDist[e.last()];
            const scalar d = s1-s0;

            return Foam::mag(d) < ROOTVSMALL ? -1 : (-s0/d);
        };

    walkCellCuts
    (
        mesh,
        cellCuts,
        edgeOrientIntersect,
        edgeAlphaIntersect,
        triangulate
    );
}


void Foam::cuttingSurface::print(Ostream& os, int level) const
{
    os  << " surface:" << surfaceName();
    if (level)
    {
        os  << "  faces:" << Mesh::surfFaces().size()
            << "  points:" << Mesh::points().size();
    }
}


// ************************************************************************* //
