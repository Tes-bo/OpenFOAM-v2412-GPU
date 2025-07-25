/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "mixerFvMesh.H"
#include "Time.H"
#include "regionSplit.H"
#include "slidingInterface.H"
#include "addToRunTimeSelectionTable.H"
#include "mapPolyMesh.H"
#include "unitConversion.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mixerFvMesh, 0);
    addToRunTimeSelectionTable(topoChangerFvMesh, mixerFvMesh, IOobject);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::mixerFvMesh::addZonesAndModifiers()
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

    // Add zones
    List<pointZone*> pz(1);

    // An empty zone for cut points
    pz[0] = new pointZone("cutPointZone", 0, pointZones());


    // Do face zones for slider

    List<faceZone*> fz(3);

    // Inner slider
    const word innerSliderName
    (
        motionDict_.subDict("slider").get<word>("inside")
    );
    const polyPatch& innerSlider = boundaryMesh()[innerSliderName];

    fz[0] = new faceZone
    (
        "insideSliderZone",
        identity(innerSlider.range()),
        false, // none are flipped
        0,
        faceZones()
    );

    // Outer slider
    const word outerSliderName
    (
        motionDict_.subDict("slider").get<word>("outside")
    );
    const polyPatch& outerSlider = boundaryMesh()[outerSliderName];

    fz[1] = new faceZone
    (
        "outsideSliderZone",
        identity(outerSlider.range()),
        false, // none are flipped
        1,
        faceZones()
    );

    // An empty zone for cut faces
    fz[2] = new faceZone("cutFaceZone", 2, faceZones());

    List<cellZone*> cz(1);

    // Mark every cell with its topological region
    regionSplit rs(*this);

    // Get the region of the cell containing the origin.
    const label originRegion = rs[findNearestCell(csys_.origin())];

    labelList movingCells(nCells());
    label nMovingCells = 0;

    forAll(rs, celli)
    {
        if (rs[celli] == originRegion)
        {
            movingCells[nMovingCells] = celli;
            ++nMovingCells;
        }
    }

    movingCells.resize(nMovingCells);
    Info<< "Number of cells in the moving region: " << nMovingCells << endl;

    cz[0] = new cellZone
    (
        "movingCells",
        std::move(movingCells),
        0,
        cellZones()
    );

    Info<< "Adding point, face and cell zones" << endl;
    addZones(pz, fz, cz);

    // Add a topology modifier
    Info<< "Adding topology modifiers" << endl;
    topoChanger_.setSize(1);
    topoChanger_.set
    (
        0,
        new slidingInterface
        (
            "mixerSlider",
            0,
            topoChanger_,
            outerSliderName + "Zone",
            innerSliderName + "Zone",
            "cutPointZone",
            "cutFaceZone",
            outerSliderName,
            innerSliderName,
            slidingInterface::INTEGRAL
        )
    );
    topoChanger_.writeOpt(IOobject::AUTO_WRITE);

    write();
}


void Foam::mixerFvMesh::calcMovingMasks() const
{
    DebugInFunction << "Calculating point and cell masks" << endl;

    if (movingPointsMaskPtr_)
    {
        FatalErrorInFunction
            << "point mask already calculated"
            << abort(FatalError);
    }

    // Set the point mask
    movingPointsMaskPtr_ = new scalarField(points().size(), Zero);
    scalarField& movingPointsMask = *movingPointsMaskPtr_;

    const cellList& c = cells();
    const faceList& f = faces();

    const labelList& cellAddr = cellZones()["movingCells"];

    for (const label celli : cellAddr)
    {
        const cell& curCell = c[celli];

        for (const label facei : curCell)
        {
            // Mark all the points as moving
            const face& curFace = f[facei];

            forAll(curFace, pointi)
            {
                movingPointsMask[curFace[pointi]] = 1;
            }
        }
    }

    const word innerSliderZoneName
    (
        motionDict_.subDict("slider").get<word>("inside") + "Zone"
    );

    const labelList& innerSliderAddr = faceZones()[innerSliderZoneName];

    for (const label facei : innerSliderAddr)
    {
        const face& curFace = f[facei];

        forAll(curFace, pointi)
        {
            movingPointsMask[curFace[pointi]] = 1;
        }
    }

    const word outerSliderZoneName
    (
        motionDict_.subDict("slider").get<word>("outside") + "Zone"
    );

    const labelList& outerSliderAddr = faceZones()[outerSliderZoneName];

    for (const label facei : outerSliderAddr)
    {
        const face& curFace = f[facei];

        forAll(curFace, pointi)
        {
            movingPointsMask[curFace[pointi]] = 0;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mixerFvMesh::mixerFvMesh
(
    const IOobject& io
)
:
    topoChangerFvMesh(io),
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
    ),
    csys_(),
    rpm_(motionDict_.get<scalar>("rpm")),
    movingPointsMaskPtr_(nullptr)
{
    // New() for access to indirect (global) coordSystem.

    auto csysPtr = coordinateSystem::NewIfPresent(*this, dict);

    if (csysPtr)
    {
        static_cast<coordinateSystem&>(csys_) = csysPtr();
    }
    else
    {
        csys_ = coordSystem::cylindrical(motionDict_);
    }

    addZonesAndModifiers();

    Info<< "Mixer mesh:" << nl
        << "    origin: " << csys_.origin() << nl
        << "    axis: " << csys_.e3() << nl
        << "    rpm: " << rpm_ << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mixerFvMesh::~mixerFvMesh()
{
    deleteDemandDrivenData(movingPointsMaskPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Return moving points mask.  Moving points marked with 1
const Foam::scalarField& Foam::mixerFvMesh::movingPointsMask() const
{
    if (!movingPointsMaskPtr_)
    {
        calcMovingMasks();
    }

    return *movingPointsMaskPtr_;
}


bool Foam::mixerFvMesh::update()
{
    // The tangential sweep (radians)
    const vector theta(0, rpmToRads(rpm_)*time().deltaTValue(), 0);

    movePoints
    (
        csys_.globalPosition
        (
            csys_.localPosition(points())
          + theta
            *movingPointsMask()
        )
    );

    // Make changes. Use inflation (so put new points in topoChangeMap)
    autoPtr<mapPolyMesh> topoChangeMap = topoChanger_.changeMesh(true);

    if (topoChangeMap)
    {
        DebugInFunction << "Mesh topology is changing" << nl;

        deleteDemandDrivenData(movingPointsMaskPtr_);
    }

    movePoints
    (
        csys_.globalPosition
        (
            csys_.localPosition(oldPoints())
          + theta
            *movingPointsMask()
        )
    );

    return true;
}


// ************************************************************************* //
