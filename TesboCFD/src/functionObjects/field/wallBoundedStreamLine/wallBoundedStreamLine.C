/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2022 OpenCFD Ltd.
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

#include "wallBoundedStreamLine.H"
#include "wallBoundedStreamLineParticleCloud.H"
#include "sampledSet.H"
#include "faceSet.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(wallBoundedStreamLine, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        wallBoundedStreamLine,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::Tuple2<Foam::tetIndices, Foam::point>
Foam::functionObjects::wallBoundedStreamLine::findNearestTet
(
    const bitSet& isWallPatch,
    const point& seedPt,
    const label celli
) const
{
    const cell& cFaces = mesh_.cells()[celli];

    label minFacei = -1;
    label minTetPti = -1;
    scalar minDistSqr = sqr(GREAT);
    point nearestPt(GREAT, GREAT, GREAT);

    for (label facei : cFaces)
    {
        if (isWallPatch[facei])
        {
            const face& f = mesh_.faces()[facei];
            const label fp0 = mesh_.tetBasePtIs()[facei];
            const point& basePoint = mesh_.points()[f[fp0]];

            label fp = f.fcIndex(fp0);
            for (label i = 2; i < f.size(); i++)
            {
                const point& thisPoint = mesh_.points()[f[fp]];
                const label nextFp = f.fcIndex(fp);
                const point& nextPoint = mesh_.points()[f[nextFp]];

                const triPointRef tri(basePoint, thisPoint, nextPoint);

                //const scalar d2 = magSqr(tri.centre() - seedPt);
                const pointHit nearInfo(tri.nearestPoint(seedPt));
                const scalar d2 = nearInfo.distance();
                if (d2 < minDistSqr)
                {
                    nearestPt = nearInfo.point();
                    minDistSqr = d2;
                    minFacei = facei;
                    minTetPti = i-1;
                }
                fp = nextFp;
            }
        }
    }

    // Return tet and nearest point on wall triangle
    return Tuple2<Foam::tetIndices, Foam::point>
    (
        tetIndices
        (
            celli,
            minFacei,
            minTetPti
        ),
        nearestPt
    );
}


Foam::point Foam::functionObjects::wallBoundedStreamLine::pushIn
(
    const triPointRef& tri,
    const point& pt
) const
{
    //pointHit nearPt(tri.nearestPoint(pt));
    //if (nearPt.distance() > ROOTSMALL)
    //{
    //    FatalErrorInFunction << "tri:" << tri
    //        << " seed:" << pt << exit(FatalError);
    //}

    return (1.0 - ROOTSMALL)*pt + ROOTSMALL*tri.centre();
}


void Foam::functionObjects::wallBoundedStreamLine::track()
{
    // Determine the 'wall' patches
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // These are the faces that need to be followed

    autoPtr<indirectPrimitivePatch> boundaryPatch(wallPatch());
    bitSet isWallPatch(mesh_.nFaces(), boundaryPatch().addressing());


    // Find nearest wall particle for the seedPoints
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    wallBoundedStreamLineParticleCloud particles
    (
        mesh_,
        Foam::zero{},
        cloudName_
    );

    {
        // Get the seed points
        // ~~~~~~~~~~~~~~~~~~~

        const sampledSet& seedPoints = sampledSetPoints();

        forAll(seedPoints, seedi)
        {
            const label celli = seedPoints.cells()[seedi];

            if (celli != -1)
            {
                const point& seedPt = seedPoints[seedi];
                Tuple2<tetIndices, point> nearestId
                (
                    findNearestTet(isWallPatch, seedPt, celli)
                );
                const tetIndices& ids = nearestId.first();

                if (ids.face() != -1 && isWallPatch[ids.face()])
                {
                    //Pout<< "Seeding particle :" << nl
                    //    << "     seedPt:" << seedPt << nl
                    //    << "     face  :" << ids.face() << nl
                    //    << "     at    :" << mesh_.faceCentres()[ids.face()]
                    //    << nl
                    //    << "     cell  :" << mesh_.cellCentres()[ids.cell()]
                    //    << nl << endl;

                    particles.addParticle
                    (
                        new wallBoundedStreamLineParticle
                        (
                            mesh_,
                            // Perturb seed point to be inside triangle
                            pushIn(ids.faceTri(mesh_), nearestId.second()),
                            ids.cell(),
                            ids.face(),     // tetFace
                            ids.tetPt(),
                            -1,             // not on a mesh edge
                            -1,             // not on a diagonal edge
                            (trackDir_ == trackDirType::FORWARD), // forward?
                            lifeTime_       // lifetime
                        )
                    );

                    if (trackDir_ == trackDirType::BIDIRECTIONAL)
                    {
                        // Additional particle for other half of track
                        particles.addParticle
                        (
                            new wallBoundedStreamLineParticle
                            (
                                mesh_,
                                // Perturb seed point to be inside triangle
                                pushIn(ids.faceTri(mesh_), nearestId.second()),
                                ids.cell(),
                                ids.face(),     // tetFace
                                ids.tetPt(),
                                -1,             // not on a mesh edge
                                -1,             // not on a diagonal edge
                                true,           // forward
                                lifeTime_       // lifetime
                            )
                        );
                    }
                }
                else
                {
                    Pout<< type() << " : ignoring seed " << seedPt
                        << " since not in wall cell." << endl;
                }
            }
        }
    }

    const label nSeeds = returnReduce(particles.size(), sumOp<label>());

    Log << type() << " : seeded " << nSeeds << " particles." << endl;


    // Field interpolators
    PtrList<interpolation<scalar>> vsInterp;
    PtrList<interpolation<vector>> vvInterp;

    refPtr<interpolation<vector>> UInterp
    (
        initInterpolations(nSeeds, vsInterp, vvInterp)
    );

    // Additional particle info
    wallBoundedStreamLineParticle::trackingData td
    (
        particles,
        vsInterp,
        vvInterp,
        UInterp.cref(), // velocity interpolator (possibly within vvInterp)
        trackLength_,   // fixed track length
        isWallPatch,    // which faces are to follow

        allTracks_,
        allScalars_,
        allVectors_
    );


    // Set very large dt. Note: cannot use GREAT since 1/GREAT is SMALL
    // which is a trigger value for the tracking...
    const scalar trackTime = Foam::sqrt(GREAT);

    // Track
    particles.move(particles, td, trackTime);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::wallBoundedStreamLine::wallBoundedStreamLine
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    streamLineBase(name, runTime, dict)
{
    read(dict_);
}


Foam::functionObjects::wallBoundedStreamLine::wallBoundedStreamLine
(
    const word& name,
    const Time& runTime,
    const dictionary& dict,
    const wordList& fieldNames
)
:
    streamLineBase(name, runTime, dict, fieldNames)
{
    read(dict_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::wallBoundedStreamLine::read(const dictionary& dict)
{
    if (streamLineBase::read(dict))
    {
        // Make sure that the mesh is trackable
        if (debug)
        {
            // 1. Positive volume decomposition tets
            faceSet faces(mesh_, "lowQualityTetFaces", mesh_.nFaces()/100+1);

            if
            (
                polyMeshTetDecomposition::checkFaceTets
                (
                    mesh_,
                    polyMeshTetDecomposition::minTetQuality,
                    true,
                    &faces
                )
            )
            {
                label nFaces = returnReduce(faces.size(), sumOp<label>());

                WarningInFunction
                    << "Found " << nFaces
                    <<" faces with low quality or negative volume "
                    << "decomposition tets. Writing to faceSet " << faces.name()
                    << endl;
            }

            // 2. All edges on a cell having two faces
            EdgeMap<label> numFacesPerEdge;
            for (const cell& cFaces : mesh_.cells())
            {
                numFacesPerEdge.clear();

                for (const label facei : cFaces)
                {
                    const face& f = mesh_.faces()[facei];
                    forAll(f, fp)
                    {
                        const edge e(f[fp], f.nextLabel(fp));

                        ++(numFacesPerEdge(e, 0));
                    }
                }

                forAllConstIters(numFacesPerEdge, iter)
                {
                    if (iter() != 2)
                    {
                        FatalErrorInFunction
                            << "problem cell:" << cFaces
                            << abort(FatalError);
                    }
                }
            }
        }
    }

    return true;
}


// ************************************************************************* //
