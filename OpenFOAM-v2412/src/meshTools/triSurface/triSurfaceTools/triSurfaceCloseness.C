/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "triSurfaceTools.H"

#include "triSurface.H"
#include "triSurfaceMesh.H"
#include "triSurfaceFields.H"
#include "MeshedSurface.H"
#include "OFstream.H"
#include "unitConversion.H"
#include "meshTools.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

static void drawHitProblem
(
    label fI,
    const triSurface& surf,
    const pointField& start,
    const pointField& faceCentres,
    const pointField& end,
    const List<pointIndexHit>& hitInfo
)
{
    Info<< nl << "# findLineAll did not hit its own face."
        << nl << "# fI " << fI
        << nl << "# start " << start[fI]
        << nl << "# f centre " << faceCentres[fI]
        << nl << "# end " << end[fI]
        << nl << "# hitInfo " << hitInfo
        << endl;

    meshTools::writeOBJ(Info, start[fI]);
    meshTools::writeOBJ(Info, faceCentres[fI]);
    meshTools::writeOBJ(Info, end[fI]);

    Info<< "l 1 2 3" << endl;

    meshTools::writeOBJ(Info, surf.points()[surf[fI][0]]);
    meshTools::writeOBJ(Info, surf.points()[surf[fI][1]]);
    meshTools::writeOBJ(Info, surf.points()[surf[fI][2]]);

    Info<< "f 4 5 6" << endl;

    forAll(hitInfo, hI)
    {
        label hFI = hitInfo[hI].index();

        meshTools::writeOBJ(Info, surf.points()[surf[hFI][0]]);
        meshTools::writeOBJ(Info, surf.points()[surf[hFI][1]]);
        meshTools::writeOBJ(Info, surf.points()[surf[hFI][2]]);

        Info<< "f "
            << 3*hI + 7 << " "
            << 3*hI + 8 << " "
            << 3*hI + 9
            << endl;
    }
}

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::Pair<Foam::tmp<Foam::scalarField>>
Foam::triSurfaceTools::writeCloseness
(
    const Time& runTime,
    const word& basename,
    const triSurface& surf,
    const scalar internalAngleTolerance,
    const scalar externalAngleTolerance
)
{
    Pair<tmp<scalarField>> tpair
    (
        tmp<scalarField>::New(surf.size(), GREAT),
        tmp<scalarField>::New(surf.size(), GREAT)
    );

    Info<< "Extracting internal and external closeness of surface." << endl;

    triSurfaceMesh searchSurf
    (
        IOobject
        (
            basename + ".closeness",
            runTime.constant(),
            "triSurface",
            runTime,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        surf
    );

    // Prepare start and end points for intersection tests

    const vectorField& normals = searchSurf.faceNormals();
    const scalar span = searchSurf.bounds().mag();

    const scalar externalToleranceCosAngle =
        Foam::cos
        (
            degToRad(180 - externalAngleTolerance)
        );

    const scalar internalToleranceCosAngle =
        Foam::cos
        (
            degToRad(180 - internalAngleTolerance)
        );

    Info<< "externalToleranceCosAngle: " << externalToleranceCosAngle << nl
        << "internalToleranceCosAngle: " << internalToleranceCosAngle << endl;

    // Info<< "span " << span << endl;

    const pointField start(searchSurf.faceCentres() - span*normals);
    const pointField end(searchSurf.faceCentres() + span*normals);
    const pointField& faceCentres = searchSurf.faceCentres();

    List<List<pointIndexHit>> allHitInfo;

    // Find all intersections (in order)
    searchSurf.findLineAll(start, end, allHitInfo);

    scalarField& internalCloseness = tpair[0].ref();
    scalarField& externalCloseness  = tpair[1].ref();

    forAll(allHitInfo, fI)
    {
        const List<pointIndexHit>& hitInfo = allHitInfo[fI];

        if (hitInfo.size() < 1)
        {
            drawHitProblem(fI, surf, start, faceCentres, end, hitInfo);

            // FatalErrorInFunction
            //     << "findLineAll did not hit its own face."
            //     << exit(FatalError);
        }
        else if (hitInfo.size() == 1)
        {
            if (!hitInfo[0].hit())
            {
                // FatalErrorInFunction
                //     << "findLineAll did not hit any face."
                //     << exit(FatalError);
            }
            else if (hitInfo[0].index() != fI)
            {
                drawHitProblem
                (
                    fI,
                    surf,
                    start,
                    faceCentres,
                    end,
                    hitInfo
                );

                // FatalErrorInFunction
                //     << "findLineAll did not hit its own face."
                //     << exit(FatalError);
            }
        }
        else
        {
            label ownHitI = -1;

            forAll(hitInfo, hI)
            {
                // Find the hit on the triangle that launched the ray

                if (hitInfo[hI].index() == fI)
                {
                    ownHitI = hI;

                    break;
                }
            }

            if (ownHitI < 0)
            {
                drawHitProblem
                (
                    fI,
                    surf,
                    start,
                    faceCentres,
                    end,
                    hitInfo
                );

                // FatalErrorInFunction
                //     << "findLineAll did not hit its own face."
                //     << exit(FatalError);
            }
            else if (ownHitI == 0)
            {
                // There are no internal hits, the first hit is the
                // closest external hit

                if
                (
                    (
                        normals[fI]
                      & normals[hitInfo[ownHitI + 1].index()]
                    )
                  < externalToleranceCosAngle
                )
                {
                    externalCloseness[fI] =
                        faceCentres[fI].dist(hitInfo[ownHitI + 1].point());
                }
            }
            else if (ownHitI == hitInfo.size() - 1)
            {
                // There are no external hits, the last but one hit is
                // the closest internal hit

                if
                (
                    (
                        normals[fI]
                      & normals[hitInfo[ownHitI - 1].index()]
                    )
                  < internalToleranceCosAngle
                )
                {
                    internalCloseness[fI] =
                        faceCentres[fI].dist(hitInfo[ownHitI - 1].point());
                }
            }
            else
            {
                if
                (
                    (
                        normals[fI]
                      & normals[hitInfo[ownHitI + 1].index()]
                    )
                  < externalToleranceCosAngle
                )
                {
                    externalCloseness[fI] =
                        faceCentres[fI].dist(hitInfo[ownHitI + 1].point());
                }

                if
                (
                    (
                        normals[fI]
                      & normals[hitInfo[ownHitI - 1].index()]
                    )
                  < internalToleranceCosAngle
                )
                {
                    internalCloseness[fI] =
                        faceCentres[fI].dist(hitInfo[ownHitI - 1].point());
                }
            }
        }
    }

    // write as 'internalCloseness'
    {
        triSurfaceScalarField outputField
        (
            IOobject
            (
                basename + ".internalCloseness",
                runTime.constant(),
                "triSurface",
                runTime,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            surf,
            dimLength,
            scalarField()
        );

        outputField.swap(internalCloseness);
        outputField.write();
        outputField.swap(internalCloseness);
    }

    // write as 'externalCloseness'
    {
        triSurfaceScalarField outputField
        (
            IOobject
            (
                basename + ".externalCloseness",
                runTime.constant(),
                "triSurface",
                runTime,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            surf,
            dimLength,
            scalarField()
        );

        outputField.swap(externalCloseness);
        outputField.write();
        outputField.swap(externalCloseness);
    }

    return tpair;
}


// ************************************************************************* //
