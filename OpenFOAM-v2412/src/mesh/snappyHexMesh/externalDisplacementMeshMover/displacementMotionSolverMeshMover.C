/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2015 OpenCFD Ltd.
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

#include "displacementMotionSolverMeshMover.H"
#include "addToRunTimeSelectionTable.H"
#include "pointConstraints.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(displacementMotionSolverMeshMover, 1);

    addToRunTimeSelectionTable
    (
        externalDisplacementMeshMover,
        displacementMotionSolverMeshMover,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::displacementMotionSolverMeshMover::moveMesh
(
    const dictionary& moveDict,
    const label nAllowableErrors,
    labelList& checkFaces
)
{
    const label nRelaxIter = moveDict.get<label>("nRelaxIter");

    meshMover_.setDisplacementPatchFields();

    Info<< typeName << " : Moving mesh ..." << endl;

    scalar oldErrorReduction = -1;

    bool meshOk = false;

    for (label iter = 0; iter < 2*nRelaxIter; ++ iter)
    {
        Info<< typeName << " : Iteration " << iter << endl;

        if (iter == nRelaxIter)
        {
            Info<< typeName
                << " : Displacement scaling for error reduction set to 0."
                << endl;
            oldErrorReduction = meshMover_.setErrorReduction(0.0);
        }

        if
        (
            meshMover_.scaleMesh
            (
                checkFaces,
                baffles_,
                meshMover_.paramDict(),
                moveDict,
                true,
                nAllowableErrors
            )
        )
        {
            Info<< typeName << " : Successfully moved mesh" << endl;
            meshOk = true;
            break;
        }
    }

    if (oldErrorReduction >= 0)
    {
        meshMover_.setErrorReduction(oldErrorReduction);
    }

    Info<< typeName << " : Finished moving mesh ..." << endl;

    return meshOk;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::displacementMotionSolverMeshMover::displacementMotionSolverMeshMover
(
    const dictionary& dict,
    const List<labelPair>& baffles,
    pointVectorField& pointDisplacement,
    const bool dryRun
)
:
    externalDisplacementMeshMover(dict, baffles, pointDisplacement, dryRun),

    solverPtr_
    (
        displacementMotionSolver::New
        (
            dict.get<word>("solver"),
            pointDisplacement.mesh()(),
            IOdictionary
            (
                IOobject
                (
                    "motionSolverDict",
                    pointDisplacement.mesh().time().constant(),
                    pointDisplacement.db(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                ),
                dict
            ),
            pointDisplacement,
            pointIOField
            (
                IOobject
                (
                    "points0",
                    pointDisplacement.mesh().time().constant(),
                    pointDisplacement.db(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                ),
                pointDisplacement.mesh()().points()
            )
        )
    ),

    adaptPatchIDs_(getFixedValueBCs(pointDisplacement)),
    adaptPatchPtr_(getPatch(mesh(), adaptPatchIDs_)),

    scale_
    (
        IOobject
        (
            "scale",
            pointDisplacement.time().timeName(),
            pointDisplacement.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        pMesh(),
        dimensionedScalar("scale", dimless, 1.0)
    ),

    oldPoints_(mesh().points()),

    meshMover_
    (
        const_cast<polyMesh&>(mesh()),
        const_cast<pointMesh&>(pMesh()),
        adaptPatchPtr_(),
        pointDisplacement,
        scale_,
        oldPoints_,
        adaptPatchIDs_,
        dict,
        dryRun
    ),

    fieldSmoother_(mesh())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::displacementMotionSolverMeshMover::~displacementMotionSolverMeshMover()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::displacementMotionSolverMeshMover::move
(
    const dictionary& moveDict,
    const label nAllowableErrors,
    labelList& checkFaces
)
{
    // Correct and smooth the patch displacements so points next to
    // points where the extrusion was disabled also use less extrusion.
    // Note that this has to update the pointDisplacement boundary conditions
    // as well, not just the internal field.
    {
        const label nSmoothPatchThickness = meshRefinement::get<label>
        (
            moveDict, "nSmoothThickness", dryRun_, keyType::REGEX
        );

        const word minThicknessName = meshRefinement::get<word>
        (
            moveDict, "minThicknessName", dryRun_, keyType::REGEX, word::null
        );

        scalarField zeroMinThickness;

        if (minThicknessName == "none")
        {
            zeroMinThickness = scalarField(adaptPatchPtr_().nPoints(), Zero);
        }

        const scalarField& minThickness =
        (
            (minThicknessName == "none")
          ? zeroMinThickness
          : mesh().lookupObject<scalarField>(minThicknessName)
        );

        const bitSet isPatchMasterPoint
        (
            meshRefinement::getMasterPoints
            (
                mesh(),
                adaptPatchPtr_().meshPoints()
            )
        );

        const bitSet isPatchMasterEdge
        (
            meshRefinement::getMasterEdges
            (
                mesh(),
                adaptPatchPtr_().meshEdges
                (
                    mesh().edges(),
                    mesh().pointEdges()
                )
            )
        );

        // Smooth patch displacement

        vectorField displacement
        (
            pointDisplacement().internalField(),
            adaptPatchPtr_().meshPoints()
        );

        fieldSmoother_.minSmoothField
        (
            nSmoothPatchThickness,
            isPatchMasterPoint,
            isPatchMasterEdge,
            adaptPatchPtr_(),
            minThickness,
            displacement
        );


        // scalar resid = 0;

        forAll(displacement, patchPointI)
        {
            const label pointI(adaptPatchPtr_().meshPoints()[patchPointI]);

            // resid += mag(pointDisplacement()[pointI]-displacement[patchPointI]);

            pointDisplacement()[pointI] = displacement[patchPointI];
        }

        // Take over smoothed displacements on bcs
        meshMover_.setDisplacementPatchFields();
    }

    // Use motionSolver to calculate internal displacement
    {
        solverPtr_->pointDisplacement() == pointDisplacement();
        // Force solving and constraining - just so its pointDisplacement gets
        // the correct value
        (void)solverPtr_->newPoints();
        pointDisplacement() == solverPtr_->pointDisplacement();
    }

    return moveMesh(moveDict, nAllowableErrors, checkFaces);
}


void Foam::displacementMotionSolverMeshMover::movePoints(const pointField& p)
{
    externalDisplacementMeshMover::movePoints(p);

    // Update motion solver for new geometry
    solverPtr_->movePoints(p);

    // Update motionSmoother for new geometry (moves adaptPatchPtr_)
    meshMover_.movePoints();

    // Assume current mesh location is correct (reset oldPoints, scale)
    meshMover_.correct();
}


// ************************************************************************* //
