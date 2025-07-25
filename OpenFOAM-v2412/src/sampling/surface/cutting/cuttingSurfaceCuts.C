/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "cuttingSurface.H"
#include "fvMesh.H"
#include "volFields.H"
#include "pointFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::cuttingSurface::calcCellCuts
(
    const fvMesh& fvm,
    scalarField& pointDist,
    bitSet& cellCuts
)
{
    const pointField& cc = fvm.C();
    const pointField& pts = fvm.points();

    const label nCells  = fvm.nCells();
    const label nPoints = fvm.nPoints();

    // The point distances

    List<pointIndexHit> nearest;
    surfPtr_().findNearest
    (
        pts,
        scalarField(nPoints, GREAT),
        nearest
    );

    vectorField norms;
    surfPtr_().getNormal(nearest, norms);

    pointDist.resize(nPoints);

    for (label i=0; i < nPoints; ++i)
    {
        const point diff(pts[i] - nearest[i].hitPoint());

        // Normal distance
        pointDist[i] = (diff & norms[i]);
    }


    if (cellCuts.size())
    {
        cellCuts.resize(nCells);    // safety
    }
    else
    {
        cellCuts.resize(nCells, true);
    }


    // Check based on cell distance

    surfPtr_().findNearest
    (
        cc,
        scalarField(cc.size(), GREAT),
        nearest
    );


    for (const label celli : cellCuts)
    {
        if (!fvm.cellBb(celli).contains(nearest[celli].point()))
        {
            cellCuts.unset(celli);
        }
    }

    if (debug)
    {
        volScalarField cCuts
        (
            IOobject
            (
                "cuttingSurface.cellCuts",
                fvm.time().timeName(),
                fvm.time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            fvm,
            dimensionedScalar(dimless, Zero)
        );

        for (const label celli : cellCuts)
        {
            cCuts[celli] = 1;
        }

        Pout<< "Writing cellCuts:" << cCuts.objectPath() << endl;
        cCuts.write();

        pointScalarField pDist
        (
            IOobject
            (
                "cuttingSurface.pointDistance",
                fvm.time().timeName(),
                fvm.time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            pointMesh::New(fvm),
            dimensionedScalar(dimLength, Zero)
        );
        pDist.primitiveFieldRef() = pointDist;

        Pout<< "Writing point distance:" << pDist.objectPath() << endl;
        pDist.write();
    }
}


// ************************************************************************* //
