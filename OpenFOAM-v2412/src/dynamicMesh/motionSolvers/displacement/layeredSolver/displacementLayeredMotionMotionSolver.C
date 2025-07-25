/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2023 OpenCFD Ltd.
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

#include "displacementLayeredMotionMotionSolver.H"
#include "addToRunTimeSelectionTable.H"
#include "pointEdgeStructuredWalk.H"
#include "pointFields.H"
#include "PointEdgeWave.H"
#include "syncTools.H"
#include "interpolationTable.H"
#include "mapPolyMesh.H"
#include "pointConstraints.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(displacementLayeredMotionMotionSolver, 0);

    addToRunTimeSelectionTable
    (
        motionSolver,
        displacementLayeredMotionMotionSolver,
        dictionary
    );

    addToRunTimeSelectionTable
    (
        displacementMotionSolver,
        displacementLayeredMotionMotionSolver,
        displacement
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::coordSystem::cylindrical&
Foam::displacementLayeredMotionMotionSolver::getCylindrical
(
    const label cellZoneI,
    const dictionary& zoneDict
)
{
    auto iter = cylSystems_.cfind(cellZoneI);

    if (iter.good())
    {
        return *(iter.val());
    }

    cylSystems_.set(cellZoneI, new coordSystem::cylindrical(zoneDict));

    return *cylSystems_[cellZoneI];
}


void Foam::displacementLayeredMotionMotionSolver::calcZoneMask
(
    const label cellZoneI,
    bitSet& isZonePoint,
    bitSet& isZoneEdge
) const
{
    isZonePoint.resize(mesh().nPoints());
    isZoneEdge.resize(mesh().nEdges());

    if (cellZoneI == -1)
    {
        isZonePoint = true;
        isZoneEdge = true;
        return;
    }


    isZonePoint.reset();
    isZoneEdge.reset();

    const cellZone& cz = mesh().cellZones()[cellZoneI];

    // Mark points, edges inside cellZone
    for (const label celli : cz)
    {
        isZonePoint.set(mesh().cellPoints(celli));
        isZoneEdge.set(mesh().cellEdges(celli));
    }

    syncTools::syncPointList
    (
        mesh(),
        isZonePoint,
        orEqOp<unsigned int>(),
        0
    );

    syncTools::syncEdgeList
    (
        mesh(),
        isZoneEdge,
        orEqOp<unsigned int>(),
        0
    );

    DebugInfo
        << "On cellZone " << cz.name() << " marked "
        << returnReduce(isZonePoint.count(), sumOp<label>()) << " points and "
        << returnReduce(isZoneEdge.count(), sumOp<label>()) << " edges" << nl;
}


// Find distance to starting point
void Foam::displacementLayeredMotionMotionSolver::walkStructured
(
    const label cellZoneI,
    const bitSet& isZonePoint,
    const bitSet& isZoneEdge,
    const labelList& seedPoints,
    const vectorField& seedData,
    scalarField& distance,
    vectorField& data,
    labelField& patchPoints
) const
{
    List<pointEdgeStructuredWalk> seedInfo(seedPoints.size());

    forAll(seedPoints, i)
    {
        const label pointi = seedPoints[i];

        seedInfo[i] = pointEdgeStructuredWalk
        (
            points0()[pointi],  // location of data
            points0()[pointi],  // previous location
            0.0,
            seedData[i],
            pointi              // pass thru seed point id
        );
    }

    // Current info on points
    List<pointEdgeStructuredWalk> allPointInfo(mesh().nPoints());

    // Mark points inside cellZone.
    // Note that we use points0, not mesh.points()
    // so as not to accumulate errors.
    for (const label pointi : isZonePoint)
    {
        allPointInfo[pointi] = pointEdgeStructuredWalk
        (
            points0()[pointi],  // location of data
            point::max,         // not valid
            0.0,
            Zero                // passive data
        );
    }


    // Current info on edges
    List<pointEdgeStructuredWalk> allEdgeInfo(mesh().nEdges());

    // Mark edges inside cellZone
    for (const label edgei : isZoneEdge)
    {
        allEdgeInfo[edgei] = pointEdgeStructuredWalk
        (
            mesh().edges()[edgei].centre(points0()),    // location of data
            point::max,                                 // not valid
            0.0,
            Zero
        );
    }

    // Walk
    PointEdgeWave<pointEdgeStructuredWalk> wallCalc
    (
        mesh(),
        seedPoints,
        seedInfo,

        allPointInfo,
        allEdgeInfo,
        mesh().globalData().nTotalPoints()  // max iterations
    );

    // Extract distance and passive data
    for (const label pointi : isZonePoint)
    {
        const auto& pointInfo = allPointInfo[pointi];

        distance[pointi] = pointInfo.dist();
        data[pointi] = pointInfo.data();

        // Optional information
        if (patchPoints.size())
        {
            patchPoints[pointi] = pointInfo.index();
        }
    }
}


// Evaluate faceZone patch
Foam::tmp<Foam::vectorField>
Foam::displacementLayeredMotionMotionSolver::faceZoneEvaluate
(
    const faceZone& fz,
    const labelList& meshPoints,
    const dictionary& dict,
    const PtrList<pointVectorField>& patchDisp,
    const label patchi
) const
{
    auto tfld = tmp<vectorField>::New(meshPoints.size());
    auto& fld = tfld.ref();

    const word type(dict.get<word>("type"));

    if (type == "fixedValue")
    {
        fld.assign("value", dict, meshPoints.size());
    }
    else if (type == "timeVaryingUniformFixedValue")
    {
        interpolationTable<vector> timeSeries(dict);

        fld = timeSeries(mesh().time().timeOutputValue());
    }
    else if (type == "slip")
    {
        if ((patchi % 2) != 1)
        {
            FatalIOErrorInFunction(*this)
                << "FaceZone:" << fz.name()
                << exit(FatalIOError);
        }
        // Use field set by previous bc
        fld = vectorField(patchDisp[patchi - 1], meshPoints);
    }
    else if (type == "follow")
    {
        // Only on boundary faces - follow boundary conditions
        fld = vectorField(pointDisplacement_, meshPoints);
    }
    else if (type == "uniformFollow")
    {
        // Reads name of name of patch. Then get average point displacement on
        // patch. That becomes the value of fld.
        const word patchName(dict.get<word>("patch"));
        const label patchID = mesh().boundaryMesh().findPatchID(patchName);
        pointField pdf
        (
            pointDisplacement_.boundaryField()[patchID].patchInternalField()
        );
        fld = gAverage(pdf);
    }
    else
    {
        FatalIOErrorInFunction(*this)
            << "Unknown faceZonePatch type " << type
            << " for faceZone " << fz.name() << nl
            << exit(FatalIOError);
    }

    return tfld;
}


void Foam::displacementLayeredMotionMotionSolver::cellZoneSolve
(
    const label cellZoneI,
    const dictionary& zoneDict
)
{
    bitSet isZonePoint, isZoneEdge;
    calcZoneMask(cellZoneI, isZonePoint, isZoneEdge);

    const dictionary& patchesDict = zoneDict.subDict("boundaryField");

    if (patchesDict.size() != 2)
    {
        FatalIOErrorInFunction(*this)
            << "Two faceZones (patches) must be specified per cellZone. "
            << " cellZone:" << cellZoneI
            << " patches:" << patchesDict.toc()
            << exit(FatalIOError);
    }

    PtrList<labelField> patchPoints(patchesDict.size());
    PtrList<scalarField> patchDist(patchesDict.size());
    PtrList<pointVectorField> patchDisp(patchesDict.size());

    // Allocate the fields
    label patchi = 0;
    for (const entry& dEntry : patchesDict)
    {
        const word& faceZoneName = dEntry.keyword();
        const label zoneI = mesh().faceZones().findZoneID(faceZoneName);
        if (zoneI == -1)
        {
            FatalIOErrorInFunction(*this)
                << "Cannot find faceZone " << faceZoneName
                << endl << "Valid zones are " << mesh().faceZones().names()
                << exit(FatalIOError);
        }

        // Determine the points of the faceZone within the cellZone
        const faceZone& fz = mesh().faceZones()[zoneI];

        patchPoints.set(patchi, new labelField(mesh().nPoints(), label(-1)));
        patchDist.set(patchi, new scalarField(mesh().nPoints()));
        patchDisp.set
        (
            patchi,
            new pointVectorField
            (
                mesh().newIOobject
                (
                    mesh().cellZones()[cellZoneI].name() + "_" + fz.name()
                ),
                pointDisplacement_  // to inherit the boundary conditions
            )
        );

        ++patchi;
    }



    // 'correctBoundaryConditions'
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Loops over all the faceZones and walks their boundary values

    // Make sure we can pick up bc values from field
    pointDisplacement_.correctBoundaryConditions();

    patchi = 0;
    for (const entry& dEntry : patchesDict)
    {
        const word& faceZoneName = dEntry.keyword();
        const dictionary& faceZoneDict = dEntry.dict();

        // Determine the points of the faceZone within the cellZone
        const faceZone& fz = mesh().faceZones()[faceZoneName];
        const labelList& fzMeshPoints = fz().meshPoints();
        DynamicList<label> meshPoints(fzMeshPoints.size());
        forAll(fzMeshPoints, i)
        {
            if (isZonePoint[fzMeshPoints[i]])
            {
                meshPoints.append(fzMeshPoints[i]);
            }
        }

        // Get initial value for all the faceZone points
        tmp<vectorField> tseed = faceZoneEvaluate
        (
            fz,
            meshPoints,
            faceZoneDict,
            patchDisp,
            patchi
        );

        DebugInfo
            << "For cellZone:" << cellZoneI
            << " for faceZone:" << fz.name()
            << " nPoints:" << tseed().size()
            << " have patchField:"
            << " max:" << gMax(tseed())
            << " min:" << gMin(tseed())
            << " avg:" << gAverage(tseed())
            << endl;


        // Set distance and transported value
        walkStructured
        (
            cellZoneI,
            isZonePoint,
            isZoneEdge,

            meshPoints,
            tseed,
            patchDist[patchi],
            patchDisp[patchi],
            patchPoints[patchi]
        );

        // Implement real bc.
        patchDisp[patchi].correctBoundaryConditions();

        ++patchi;
    }


    // Solve
    // ~~~~~

    if (debug)
    {
        // Normalised distance
        auto tdistance = pointScalarField::New
        (
            IOobject::scopedName
            (
                mesh().cellZones()[cellZoneI].name(),
                "distance"
            ),
            IOobject::NO_REGISTER,
            pointMesh::New(mesh()),
            dimensionedScalar(dimLength, Zero)
        );
        auto& distance = tdistance.ref();

        for (const label pointi : isZonePoint)
        {
            const scalar d1 = patchDist[0][pointi];
            const scalar d2 = patchDist[1][pointi];
            if (d1 + d2 > SMALL)
            {
                distance[pointi] = d1/(d1 + d2);
            }
        }

        Info<< "Writing " << pointScalarField::typeName
            << distance.name() << " to "
            << mesh().time().timeName() << endl;
        distance.write();
    }


    const word interpolationScheme(zoneDict.get<word>("interpolationScheme"));

    if (interpolationScheme == "oneSided")
    {
        for (const label pointi : isZonePoint)
        {
            pointDisplacement_[pointi] = patchDisp[0][pointi];
        }
    }
    else if (interpolationScheme == "linear")
    {
        for (const label pointi : isZonePoint)
        {
            const scalar d1 = patchDist[0][pointi];
            const scalar d2 = patchDist[1][pointi];
            const scalar s = d1/(d1 + d2 + VSMALL);

            const vector& pd1 = patchDisp[0][pointi];
            const vector& pd2 = patchDisp[1][pointi];

            pointDisplacement_[pointi] = (1 - s)*pd1 + s*pd2;
        }
    }
    else if (interpolationScheme == "cylindrical")
    {
        const coordSystem::cylindrical& cs =
            this->getCylindrical(cellZoneI, zoneDict);

        // May wish to implement alternative distance calculation here


        FatalErrorInFunction
            << "cylindrical interpolation not yet available" << nl
            << "coordinate system " << cs << nl
            << exit(FatalError);
    }
    else
    {
        FatalErrorInFunction
            << "Invalid interpolationScheme: " << interpolationScheme
            << ". Valid schemes: (oneSided linear cylindrical)" << nl
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::displacementLayeredMotionMotionSolver::
displacementLayeredMotionMotionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict
)
:
    displacementMotionSolver(mesh, dict, typeName)
{}


Foam::displacementLayeredMotionMotionSolver::
displacementLayeredMotionMotionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict,
    const pointVectorField& pointDisplacement,
    const pointIOField& points0
)
:
    displacementMotionSolver(mesh, dict, pointDisplacement, points0, typeName)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::displacementLayeredMotionMotionSolver::curPoints() const
{
    tmp<pointField> tcurPoints
    (
        points0() + pointDisplacement_.primitiveField()
    );

    return tcurPoints;
}


void Foam::displacementLayeredMotionMotionSolver::solve()
{
    // The points have moved so before interpolation update the motionSolver
    movePoints(mesh().points());

    // Apply boundary conditions
    pointDisplacement_.boundaryFieldRef().updateCoeffs();

    // Solve motion on all regions (=cellZones)
    for (const entry& dEntry : coeffDict().subDict("regions"))
    {
        const word& cellZoneName = dEntry.keyword();
        const dictionary& regionDict = dEntry.dict();

        const label zoneI = mesh().cellZones().findZoneID(cellZoneName);

        Info<< "solving for zone: " << cellZoneName << endl;

        if (zoneI == -1)
        {
            FatalIOErrorInFunction(*this)
                << "Cannot find cellZone " << cellZoneName
                << endl << "Valid zones are " << mesh().cellZones().names()
                << exit(FatalIOError);
        }

        cellZoneSolve(zoneI, regionDict);
    }

    // Update pointDisplacement for solved values
    const pointConstraints& pcs =
        pointConstraints::New(pointDisplacement_.mesh());
    pcs.constrainDisplacement(pointDisplacement_, false);
}


void Foam::displacementLayeredMotionMotionSolver::updateMesh
(
    const mapPolyMesh& mpm
)
{
    displacementMotionSolver::updateMesh(mpm);

    const vectorField displacement(this->newPoints() - points0_);

    forAll(points0_, pointi)
    {
        const label oldPointi = mpm.pointMap()[pointi];

        if (oldPointi >= 0)
        {
            const label masterPointi = mpm.reversePointMap()[oldPointi];

            if ((masterPointi != pointi))
            {
                // newly inserted point in this cellZone

                // need to set point0 so that it represents the position that
                // it would have had if it had existed for all time
                points0_[pointi] -= displacement[pointi];
            }
        }
    }
}


// ************************************************************************* //
