/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2019 PCOpt/NTUA
    Copyright (C) 2013-2019 FOSS GP
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#include "displacementMethoddisplacementLaplacian.H"
#include "displacementLaplacianFvMotionSolver.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(displacementMethoddisplacementLaplacian, 1);
addToRunTimeSelectionTable
(
    displacementMethod,
    displacementMethoddisplacementLaplacian,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

displacementMethoddisplacementLaplacian::displacementMethoddisplacementLaplacian
(
    fvMesh& mesh,
    const labelList& patchIDs
)
:
    displacementMethod(mesh, patchIDs),
    pointDisplacement_
    (
        refCast<displacementLaplacianFvMotionSolver>
        (
            motionPtr_()
        ).pointDisplacement()
    ),
    cellDisplacement_
    (
        refCast<displacementLaplacianFvMotionSolver>
        (
            motionPtr_()
        ).cellDisplacement()
    ),
    /*
    // Getting a reference from the database is dangerous since multiple
    // fields with the same name might be registered
    pointDisplacement_
    (
        const_cast<pointVectorField&>
        (
            mesh_.lookupObject<pointVectorField>("pointDisplacement")
        )
    ),
    cellDisplacement_
    (
        const_cast<volVectorField&>
        (
            mesh_.lookupObject<volVectorField>("cellDisplacement")
        )
    ),
    */
    resetFields_
    (
        IOdictionary::readContents
        (
            IOobject
            (
                "dynamicMeshDict",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ
            )
        ).subDict("displacementLaplacianCoeffs").getOrDefault<bool>
        (
            "resetFields",
            true
        )
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void displacementMethoddisplacementLaplacian::setMotionField
(
    const pointVectorField& pointMovement
)
{
    Info<< "Resetting mesh motion fields to zero " << endl;

    if (resetFields_)
    {
        pointDisplacement_.primitiveFieldRef() = vector::zero;
        cellDisplacement_.primitiveFieldRef() = vector::zero;
        cellDisplacement_.correctBoundaryConditions();
    }

    // Set boundary mesh movement and calculate
    // max current boundary displacement
    for (label patchI : patchIDs_)
    {
        // Set boundary field. Needed for the motionSolver class
        pointDisplacement_.boundaryFieldRef()[patchI] ==
            pointMovement.boundaryField()[patchI].patchInternalField()();

        // Set boundary field values as seen from the internalField!
        // Needed for determining the max displacement
        pointDisplacement_.boundaryFieldRef()[patchI].setInInternalField
        (
            pointDisplacement_.primitiveFieldRef(),
            pointMovement.boundaryField()[patchI].patchInternalField()()
        );

        // Find max value
        maxDisplacement_ =
            max
            (
                maxDisplacement_,
                gMax
                (
                    mag
                    (
                        pointDisplacement_.boundaryField()[patchI].
                            patchInternalField()
                    )
                )
            );
    }
}


void displacementMethoddisplacementLaplacian::setMotionField
(
    const volVectorField& cellMovement
)
{
    NotImplemented;
}


void displacementMethoddisplacementLaplacian::setControlField
(
    const vectorField& controlField
)
{
    NotImplemented;
}


void displacementMethoddisplacementLaplacian::setControlField
(
    const scalarField& controlField
)
{
    NotImplemented;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
