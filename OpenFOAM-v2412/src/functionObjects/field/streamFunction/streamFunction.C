/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
    Copyright (C) 2019-2024 OpenCFD Ltd.
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

#include "streamFunction.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "emptyPolyPatch.H"
#include "symmetryPlanePolyPatch.H"
#include "symmetryPolyPatch.H"
#include "wedgePolyPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(streamFunction, 0);
    addToRunTimeSelectionTable(functionObject, streamFunction, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::pointScalarField> Foam::functionObjects::streamFunction::calc
(
    const surfaceScalarField& phi
) const
{
    Log << "    functionObjects::" << type() << " " << name()
        << " calculating stream-function" << endl;

    Vector<label> slabNormal((Vector<label>::one - mesh_.geometricD())/2);
    const direction slabDir
    (
        slabNormal
      & Vector<label>(Vector<label>::X, Vector<label>::Y, Vector<label>::Z)
    );

    const pointMesh& pMesh = pointMesh::New(mesh_);

    auto tstreamFunction = tmp<pointScalarField>::New
    (
        IOobject
        (
            "streamFunction",
            time_.timeName(),
            mesh_
        ),
        pMesh,
        dimensionedScalar(phi.dimensions(), Zero)
    );
    auto& streamFunction = tstreamFunction.ref();


    bitSet visitedPoint(mesh_.nPoints());

    label nVisited = 0;
    label nVisitedOld = 0;

    const faceUList& faces = mesh_.faces();
    const pointField& points = mesh_.points();

    const label nInternalFaces = mesh_.nInternalFaces();

    vectorField unitAreas(mesh_.faceAreas());
    unitAreas.normalise();

    const polyPatchList& patches = mesh_.boundaryMesh();

    bool finished = true;

    // Find the boundary face with zero flux. Set the stream function
    // to zero on that face
    bool found = false;

    do
    {
        found = false;

        // Check boundary faces first
        forAll(patches, patchi)
        {
            const auto& pp = patches[patchi];
            const auto& patchPhi = phi.boundaryField()[patchi];

            // Skip empty, symmetry patches etc
            if
            (
                patchPhi.empty()
             || isType<emptyPolyPatch>(pp)
             || isType<symmetryPlanePolyPatch>(pp)
             || isType<symmetryPolyPatch>(pp)
             || isType<wedgePolyPatch>(pp)
            )
            {
                continue;
            }

            forAll(pp, facei)
            {
                const auto& f = pp[facei];

                if (magSqr(patchPhi[facei]) < SMALL)
                {
                    // Zero flux face found
                    found = true;

                    for (const label pointi : f)
                    {
                        if (visitedPoint.test(pointi))
                        {
                            found = false;
                            break;
                        }
                    }

                    if (found)
                    {
                        Log << "        Zero face: patch: " << patchi
                            << " face: " << facei << endl;

                        for (const label pointi : f)
                        {
                            visitedPoint.set(pointi);
                            ++nVisited;

                            streamFunction[pointi] = 0;
                        }

                        break;
                    }
                }
            }

            if (found) break;
        }

        if (!found)
        {
            Log << "        Zero flux boundary face not found. "
                << "Using cell as a reference." << endl;

            for (const cell& c : mesh_.cells())
            {
                labelList zeroPoints = c.labels(mesh_.faces());

                bool found = true;

                for (const label pointi : zeroPoints)
                {
                    if (visitedPoint.test(pointi))
                    {
                        found = false;
                        break;
                    }
                }

                if (found)
                {
                    for (const label pointi : zeroPoints)
                    {
                        visitedPoint.set(pointi);
                        ++nVisited;

                        streamFunction[pointi] = 0;
                    }

                    break;
                }
                else
                {
                    FatalErrorInFunction
                        << "Cannot find initialisation face or a cell."
                        << exit(FatalError);
                }
            }
        }

        // Loop through all faces. If one of the points on
        // the face has the streamfunction value different
        // from -1, all points with -1 ont that face have the
        // streamfunction value equal to the face flux in
        // that point plus the value in the visited point
        do
        {
            finished = true;

            scalar currentStreamValue(0);
            point currentStreamPoint(Zero);

            // Boundary faces first
            forAll(patches, patchi)
            {
                const auto& pp = patches[patchi];
                const auto& patchPhi = phi.boundaryField()[patchi];

                // Skip empty, symmetry patches etc
                if
                (
                    patchPhi.empty()
                 || isType<emptyPolyPatch>(pp)
                 || isType<symmetryPlanePolyPatch>(pp)
                 || isType<symmetryPolyPatch>(pp)
                 || isType<wedgePolyPatch>(pp)
                )
                {
                    continue;
                }

                forAll(pp, facei)
                {
                    const auto& f = pp[facei];

                    // Check if the point has been visited
                    bool pointFound = false;

                    for (const label pointi : f)
                    {
                        if (visitedPoint.test(pointi))
                        {
                            // The point has been visited
                            currentStreamValue = streamFunction[pointi];
                            currentStreamPoint = points[pointi];

                            pointFound = true;
                            break;
                        }
                    }

                    if (!pointFound)
                    {
                        finished = false;
                        continue;
                    }


                    // Sort out other points on the face
                    for (const label pointi : f)
                    {
                        // If the point has not yet been visited
                        if (!visitedPoint.test(pointi))
                        {
                            vector edgeHat =
                            (
                                points[pointi] - currentStreamPoint
                            );
                            edgeHat.replace(slabDir, 0);
                            edgeHat.normalise();

                            const vector& nHat = unitAreas[facei];

                            if (edgeHat.y() > VSMALL)
                            {
                                visitedPoint.set(pointi);
                                ++nVisited;

                                streamFunction[pointi] =
                                (
                                    currentStreamValue
                                  + patchPhi[facei]*sign(nHat.x())
                                );
                            }
                            else if (edgeHat.y() < -VSMALL)
                            {
                                visitedPoint.set(pointi);
                                ++nVisited;

                                streamFunction[pointi] =
                                (
                                    currentStreamValue
                                  - patchPhi[facei]*sign(nHat.x())
                                );
                            }
                            else
                            {
                                if (edgeHat.x() > VSMALL)
                                {
                                    visitedPoint.set(pointi);
                                    ++nVisited;

                                    streamFunction[pointi] =
                                    (
                                        currentStreamValue
                                      + patchPhi[facei]*sign(nHat.y())
                                    );
                                }
                                else if (edgeHat.x() < -VSMALL)
                                {
                                    visitedPoint.set(pointi);
                                    ++nVisited;

                                    streamFunction[pointi] =
                                    (
                                        currentStreamValue
                                      - patchPhi[facei]*sign(nHat.y())
                                    );
                                }
                            }
                        }
                    }
                }
            }

            // Internal faces next
            for (label facei = 0; facei < nInternalFaces; ++facei)
            {
                const auto& f = faces[facei];

                bool pointFound = false;

                for (const label pointi : f)
                {
                    // Check if the point has been visited
                    if (visitedPoint.test(pointi))
                    {
                        currentStreamValue = streamFunction[pointi];
                        currentStreamPoint = points[pointi];
                        pointFound = true;

                        break;
                    }
                }

                if (pointFound)
                {
                    // Sort out other points on the face
                    for (const label pointi : f)
                    {
                        // If the point has not yet been visited
                        if (!visitedPoint.test(pointi))
                        {
                            vector edgeHat =
                            (
                                points[pointi] - currentStreamPoint
                            );

                            edgeHat.replace(slabDir, 0);
                            edgeHat.normalise();

                            const vector& nHat = unitAreas[facei];

                            if (edgeHat.y() > VSMALL)
                            {
                                visitedPoint.set(pointi);
                                ++nVisited;

                                streamFunction[pointi] =
                                (
                                    currentStreamValue
                                  + phi[facei]*sign(nHat.x())
                                );
                            }
                            else if (edgeHat.y() < -VSMALL)
                            {
                                visitedPoint.set(pointi);
                                ++nVisited;

                                streamFunction[pointi] =
                                (
                                    currentStreamValue
                                  - phi[facei]*sign(nHat.x())
                                );
                            }
                        }
                    }
                }
                else
                {
                    finished = false;
                }
            }

            if (nVisited == nVisitedOld)
            {
                // Find new seed.  This must be a
                // multiply connected domain
                Log << "        Exhausted a seed, looking for new seed "
                    << "(this is correct for multiply connected domains).";

                break;
            }
            else
            {
                nVisitedOld = nVisited;
            }
        } while (!finished);
    } while (!finished);

    // Normalise the stream-function by the 2D mesh thickness
    // calculate thickness here to avoid compiler oddness (#2603)
    const scalar thickness = vector(slabNormal) & mesh_.bounds().span();

    streamFunction /= thickness;
    streamFunction.boundaryFieldRef() = Zero;

    return tstreamFunction;
}


bool Foam::functionObjects::streamFunction::calc()
{
    const auto* phiPtr = findObject<surfaceScalarField>(fieldName_);

    if (phiPtr)
    {
        const surfaceScalarField& phi = *phiPtr;

        return store(resultName_, calc(phi));
    }

    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::streamFunction::streamFunction
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression(name, runTime, dict, "phi")
{
    setResultName(typeName, "phi");

    const label nD = mesh_.nGeometricD();

    if (nD != 2)
    {
        FatalErrorInFunction
            << "Case is not 2D, stream-function cannot be computed"
            << exit(FatalError);
    }
}


// ************************************************************************* //
