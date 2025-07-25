/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018-2020 OpenCFD Ltd.
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

#include "movingConeTopoFvMesh.H"
#include "Time.H"
#include "mapPolyMesh.H"
#include "layerAdditionRemoval.H"
#include "addToRunTimeSelectionTable.H"
#include "meshTools.H"
#include "OFstream.H"
#include "mathematicalConstants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(movingConeTopoFvMesh, 0);

    addToRunTimeSelectionTable
    (
        topoChangerFvMesh,
        movingConeTopoFvMesh,
        IOobject
    );
    addToRunTimeSelectionTable
    (
        topoChangerFvMesh,
        movingConeTopoFvMesh,
        doInit
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::movingConeTopoFvMesh::vertexMarkup
(
    const pointField& p,
    const scalar curLeft,
    const scalar curRight
) const
{
    Info<< "Updating vertex markup.  curLeft: "
        << curLeft << " curRight: " << curRight << endl;

    auto tvertexMarkup = tmp<scalarField>::New(p.size());
    auto& vertexMarkup = tvertexMarkup.ref();

    forAll(p, pI)
    {
        if (p[pI].x() < curLeft - SMALL)
        {
            vertexMarkup[pI] = -1;
        }
        else if (p[pI].x() > curRight + SMALL)
        {
            vertexMarkup[pI] = 1;
        }
        else
        {
            vertexMarkup[pI] = 0;
        }
    }

    return tvertexMarkup;
}


void Foam::movingConeTopoFvMesh::addZonesAndModifiers()
{
    // Add zones and modifiers for motion action

    if
    (
        pointZones().size()
     || faceZones().size()
     || cellZones().size()
     || topoChanger_.size()
    )
    {
        InfoInFunction
            << "Zones and modifiers already present.  Skipping."
            << endl;

        return;
    }

    Info<< "Time = " << time().timeName() << endl
        << "Adding zones and modifiers to the mesh" << endl;

    const vectorField& fc = faceCentres();
    const vectorField& fa = faceAreas();

    labelList zone1(fc.size());
    boolList flipZone1(fc.size(), false);
    label nZoneFaces1 = 0;

    labelList zone2(fc.size());
    boolList flipZone2(fc.size(), false);
    label nZoneFaces2 = 0;

    forAll(fc, facei)
    {
        if
        (
            fc[facei].x() > -0.003501
         && fc[facei].x() < -0.003499
        )
        {
            if ((fa[facei] & vector(1, 0, 0)) < 0)
            {
                flipZone1[nZoneFaces1] = true;
            }

            zone1[nZoneFaces1] = facei;
            Info<< "face " << facei << " for zone 1.  Flip: "
                << flipZone1[nZoneFaces1] << endl;
            ++nZoneFaces1;
        }
        else if
        (
            fc[facei].x() > -0.00701
         && fc[facei].x() < -0.00699
        )
        {
            zone2[nZoneFaces2] = facei;

            if ((fa[facei] & vector(1, 0, 0)) > 0)
            {
                flipZone2[nZoneFaces2] = true;
            }

            Info<< "face " << facei << " for zone 2.  Flip: "
                << flipZone2[nZoneFaces2] << endl;
            ++nZoneFaces2;
        }
    }

    zone1.setSize(nZoneFaces1);
    flipZone1.setSize(nZoneFaces1);

    zone2.setSize(nZoneFaces2);
    flipZone2.setSize(nZoneFaces2);

    Info<< "zone: " << zone1 << endl;
    Info<< "zone: " << zone2 << endl;

    List<pointZone*> pz(0);
    List<faceZone*> fz(2);
    List<cellZone*> cz(0);

    label nFz = 0;

    fz[nFz] =
        new faceZone
        (
            "rightExtrusionFaces",
            std::move(zone1),
            std::move(flipZone1),
            nFz,
            faceZones()
        );
    ++nFz;

    fz[nFz] =
        new faceZone
        (
            "leftExtrusionFaces",
            std::move(zone2),
            std::move(flipZone2),
            nFz,
            faceZones()
        );
    ++nFz;

    fz.setSize(nFz);

    Info<< "Adding mesh zones." << endl;
    addZones(pz, fz, cz);


    // Add layer addition/removal interfaces

    List<polyMeshModifier*> tm(2);
    label nMods = 0;

    tm[nMods] =
        new layerAdditionRemoval
        (
            "right",
            nMods,
            topoChanger_,
            "rightExtrusionFaces",
            motionDict_.subDict("right").get<scalar>("minThickness"),
            motionDict_.subDict("right").get<scalar>("maxThickness")
        );
    ++nMods;

    tm[nMods] = new layerAdditionRemoval
    (
        "left",
        nMods,
        topoChanger_,
        "leftExtrusionFaces",
        motionDict_.subDict("left").get<scalar>("minThickness"),
        motionDict_.subDict("left").get<scalar>("maxThickness")
    );
    ++nMods;
    tm.setSize(nMods);

    Info<< "Adding " << nMods << " mesh modifiers" << endl;
    topoChanger_.addTopologyModifiers(tm);

    write();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::movingConeTopoFvMesh::movingConeTopoFvMesh
(
    const IOobject& io,
    const bool doInit
)
:
    topoChangerFvMesh(io, doInit),
    motionDict_
    (
        IOdictionary::readContents
        (
            IOobject
            (
                "dynamicMeshDict",
                time().constant(),
                *this,
                IOobject::MUST_READ
            )
        ).optionalSubDict(typeName + "Coeffs")
    )
{
    if (doInit)
    {
        init(false);    // do not initialise lower levels
    }
}


bool Foam::movingConeTopoFvMesh::init(const bool doInit)
{
    if (doInit)
    {
        topoChangerFvMesh::init(doInit);
    }

    motionVelAmplitude_ = motionDict_.get<vector>("motionVelAmplitude");
    motionVelPeriod_ = motionDict_.get<scalar>("motionVelPeriod");
    curMotionVel_ =
        motionVelAmplitude_*sin(time().value()*pi/motionVelPeriod_);
    leftEdge_ = motionDict_.get<scalar>("leftEdge");
    curLeft_ = motionDict_.get<scalar>("leftObstacleEdge");
    curRight_ = motionDict_.get<scalar>("rightObstacleEdge");

    Pout<< "Initial time:" << time().value()
        << " Initial curMotionVel_:" << curMotionVel_
        << endl;

    addZonesAndModifiers();

    curLeft_ = average
    (
        faceZones()
        [
            faceZones().findZoneID("leftExtrusionFaces")
        ]().localPoints()
    ).x() - SMALL;

    curRight_ = average
    (
        faceZones()
        [
            faceZones().findZoneID("rightExtrusionFaces")
        ]().localPoints()
    ).x() + SMALL;

    motionMask_ = vertexMarkup
    (
        points(),
        curLeft_,
        curRight_
    );

    // Assume something changed
    return true;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::movingConeTopoFvMesh::~movingConeTopoFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::movingConeTopoFvMesh::update()
{
    // Do mesh changes (use inflation - put new points in topoChangeMap)
    autoPtr<mapPolyMesh> topoChangeMap = topoChanger_.changeMesh(true);

    // Calculate the new point positions depending on whether the
    // topological change has happened or not
    pointField newPoints;

    vector curMotionVel_ =
        motionVelAmplitude_*sin(time().value()*pi/motionVelPeriod_);

    Pout<< "time:" << time().value() << " curMotionVel_:" << curMotionVel_
        << " curLeft:" << curLeft_ << " curRight:" << curRight_
        << endl;

    if (topoChangeMap)
    {
        Info<< "Topology change. Calculating motion points" << endl;

        if (topoChangeMap().hasMotionPoints())
        {
            Info<< "Topology change. Has premotion points" << endl;

            motionMask_ =
                vertexMarkup
                (
                    topoChangeMap().preMotionPoints(),
                    curLeft_,
                    curRight_
                );

            // Move points inside the motionMask
            newPoints =
                topoChangeMap().preMotionPoints()
              + (
                    pos0(0.5 - mag(motionMask_)) // cells above the body
                )*curMotionVel_*time().deltaTValue();
        }
        else
        {
            Info<< "Topology change. Already set mesh points" << endl;

            motionMask_ =
                vertexMarkup
                (
                    points(),
                    curLeft_,
                    curRight_
                );

            // Move points inside the motionMask
            newPoints =
                points()
              + (
                    pos0(0.5 - mag(motionMask_)) // cells above the body
                )*curMotionVel_*time().deltaTValue();
        }
    }
    else
    {
        Info<< "No topology change" << endl;
        // Set the mesh motion
        newPoints =
            points()
          + (
                pos0(0.5 - mag(motionMask_)) // cells above the body
           )*curMotionVel_*time().deltaTValue();
    }

    // The mesh now contains the cells with zero volume
    Info << "Executing mesh motion" << endl;
    movePoints(newPoints);

    //  The mesh now has got non-zero volume cells

    curLeft_ = average
    (
        faceZones()
        [
            faceZones().findZoneID("leftExtrusionFaces")
        ]().localPoints()
    ).x() - SMALL;

    curRight_ = average
    (
        faceZones()
        [
            faceZones().findZoneID("rightExtrusionFaces")
        ]().localPoints()
    ).x() + SMALL;

    return true;
}


// ************************************************************************* //
