/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2022 OpenCFD Ltd.
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

#include "dynamicMultiMotionSolverFvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "bitSet.H"
#include "syncTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dynamicMultiMotionSolverFvMesh, 0);
    addToRunTimeSelectionTable
    (
        dynamicFvMesh,
        dynamicMultiMotionSolverFvMesh,
        IOobject
    );
    addToRunTimeSelectionTable
    (
        dynamicFvMesh,
        dynamicMultiMotionSolverFvMesh,
        doInit
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicMultiMotionSolverFvMesh::dynamicMultiMotionSolverFvMesh
(
    const IOobject& io,
    const bool doInit
)
:
    dynamicFvMesh(io, doInit)
{
    if (doInit)
    {
        init(false);    // do not initialise lower levels
    }
}


bool Foam::dynamicMultiMotionSolverFvMesh::init(const bool doInit)
{
    if (doInit)
    {
        dynamicFvMesh::init(doInit);
    }

    IOdictionary dynDict
    (
        IOobject
        (
            "dynamicMeshDict",
            time().constant(),
            *this,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    );
    const dictionary& dynamicMeshCoeffs = dynDict.subDict(typeName + "Coeffs");

    motionPtr_.resize(dynamicMeshCoeffs.size());
    pointIDs_.resize(dynamicMeshCoeffs.size());

    label zonei = 0;

    bitSet movePts;

    for (const entry& dEntry : dynamicMeshCoeffs)
    {
        if (dEntry.isDict())
        {
            const dictionary& subDict = dEntry.dict();

            wordRe cellZoneName;
            subDict.readEntry("cellZone", cellZoneName);

            // Also handles groups, multiple zones (as wordRe match) ...
            labelList zoneIDs = cellZones().indices(cellZoneName);

            if (zoneIDs.empty())
            {
                FatalIOErrorInFunction(dynamicMeshCoeffs)
                    << "No matching cellZones: " << cellZoneName << nl
                    << "    Valid zones : "
                    << flatOutput(cellZones().names()) << nl
                    << "    Valid groups: "
                    << flatOutput(cellZones().groupNames())
                    << nl
                    << exit(FatalIOError);
            }

            IOobject io(dynDict, IOobject::NO_READ, IOobject::NO_WRITE);

            motionPtr_.set
            (
                zonei,
                motionSolver::New
                (
                    *this,
                    IOdictionary(io, subDict)
                )
            );


            // Markup points associated with cell zone(s)

            movePts.reset();
            movePts.resize(nPoints());

            for (const label zoneID : zoneIDs)
            {
                for (const label celli : cellZones()[zoneID])
                {
                    for (const label facei : cells()[celli])
                    {
                        movePts.set(faces()[facei]);
                    }
                }
            }

            syncTools::syncPointList
            (
                *this, movePts, orEqOp<unsigned int>(), 0u
            );

            pointIDs_[zonei] = movePts.sortedToc();

            Info<< "Applying motionSolver " << motionPtr_[zonei].type()
                << " to "
                << returnReduce(pointIDs_[zonei].size(), sumOp<label>())
                << " points of cellZone " << cellZoneName << endl;

            ++zonei;
        }
    }

    motionPtr_.resize(zonei);
    pointIDs_.resize(zonei);

    // Assume changed ...
    return true;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::dynamicMultiMotionSolverFvMesh::update()
{
    pointField transformedPts(points());

    forAll(motionPtr_, zonei)
    {
        const labelList& zonePoints = pointIDs_[zonei];

        const pointField newPoints(motionPtr_[zonei].newPoints());

        for (const label pointi : zonePoints)
        {
            transformedPts[pointi] = newPoints[pointi];
        }
    }

    fvMesh::movePoints(transformedPts);

    static bool hasWarned = false;

    volVectorField* Uptr = getObjectPtr<volVectorField>("U");

    if (Uptr)
    {
        Uptr->correctBoundaryConditions();
    }
    else if (!hasWarned)
    {
        hasWarned = true;

        WarningInFunction
            << "Did not find volVectorField U."
            << " Not updating U boundary conditions." << endl;
    }

    return true;
}


// ************************************************************************* //
