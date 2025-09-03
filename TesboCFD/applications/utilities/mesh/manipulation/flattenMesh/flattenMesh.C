/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

Application
    flattenMesh

Group
    grpMeshManipulationUtilities

Description
    Flattens the front and back planes of a 2D cartesian mesh.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "emptyPolyPatch.H"
#include "twoDPointCorrector.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Flattens the front and back planes of a 2D cartesian mesh"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createPolyMesh.H"

    pointIOField points
    (
        IOobject
        (
            "points",
            runTime.findInstance(polyMesh::meshSubDir, "points"),
            polyMesh::meshSubDir,
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    );

    boundBox bb(points);

    Info<< "bounding box: min = " << bb.min()
        << " max = " << bb.max() << " metres."
        << endl;


    point midPoint = gAverage(points);

    twoDPointCorrector twoDCorr(mesh);

    direction planeNormalCmpt = twoDCorr.normalDir();

    scalar midCmptVal = midPoint[planeNormalCmpt];
    scalar minCmptVal = bb.min()[planeNormalCmpt];
    scalar maxCmptVal = bb.max()[planeNormalCmpt];

    forAll(points, pointi)
    {
        if (points[pointi][planeNormalCmpt] < midCmptVal)
        {
            points[pointi][planeNormalCmpt] = minCmptVal;
        }
        else
        {
            points[pointi][planeNormalCmpt] = maxCmptVal;
        }
    }

    twoDCorr.correctPoints(points);

    // More precision (for points data)
    IOstream::minPrecision(10);

    Info<< "Writing points into directory " << points.path() << nl << endl;
    points.write();

    Info<< nl << "End" << nl << endl;

    return 0;
}


// ************************************************************************* //
