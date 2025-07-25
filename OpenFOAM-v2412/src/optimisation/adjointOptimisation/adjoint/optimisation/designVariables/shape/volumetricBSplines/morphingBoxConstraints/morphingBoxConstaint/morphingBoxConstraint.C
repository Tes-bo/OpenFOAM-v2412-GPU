/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2023 PCOpt/NTUA
    Copyright (C) 2021-2023 FOSS GP
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

#include "morphingBoxConstraint.H"
#include "volumetricBSplinesDesignVariables.H"
#include "createZeroField.H"
#include "IOmanip.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * * //

namespace Foam
{
    defineRunTimeSelectionTable(morphingBoxConstraint, dictionary);
    defineTypeNameAndDebug(morphingBoxConstraint, 0);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::morphingBoxConstraint::writeDVSensitivities
(
    const scalarField& sens,
    const word& solverName
)
{
    if (Pstream::master())
    {
        OFstream derivFile
            (derivativesFolder_/solverName + mesh_.time().timeName());

        unsigned int width = IOstream::defaultPrecision() + 7;
        derivFile
            << setw(width) << "#varID" << " "
            << setw(width) << "adjointSensitivity"
            << endl;

        const labelList& activeVars = designVariables_.activeDesignVariables();
        forAll(activeVars, varI)
        {
            const label activeVarI = activeVars[varI];
            derivFile
                << setw(width) << activeVarI << " "
                << setw(width) << sens[activeVarI] << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::morphingBoxConstraint::morphingBoxConstraint
(
    const fvMesh& mesh,
    const dictionary& dict,
    volumetricBSplinesDesignVariables& designVariables
)
:
    mesh_(mesh),
    dict_(dict),
    designVariables_(designVariables),
    volBSplinesBase_(designVariables.getVolBSplinesBase()),
    initialCPs_(3*volBSplinesBase_.getTotalControlPointsNumber()),
    initialiseVars_(true),
    derivativesFolder_("optimisation"/type() + "Derivatives")
{
    // Store initial control points
    const PtrList<NURBS3DVolume>& boxes = volBSplinesBase_.boxes();
    label varID(0);
    for (const NURBS3DVolume& boxI : boxes)
    {
        const vectorField& cps = boxI.getControlPoints();
        for (const vector& cpI : cps)
        {
            initialCPs_[varID++] = cpI.x();
            initialCPs_[varID++] = cpI.y();
            initialCPs_[varID++] = cpI.z();
        }
    }

    // Create sensitivities folder
    mkDir(derivativesFolder_);
}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::morphingBoxConstraint> Foam::morphingBoxConstraint::New
(
    const fvMesh& mesh,
    const dictionary& dict,
    volumetricBSplinesDesignVariables& designVariables
)
{
    const word modelType(dict.getOrDefault<word>("constraintType", "none"));

    Info<< "morphingBoxConstraint type : " << modelType << endl;

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalErrorInLookup
        (
            "constraintType",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalError);
    }

    return autoPtr<morphingBoxConstraint>
    (
        ctorPtr(mesh, dict, designVariables)
    );
}


// * * * * * * * * * * * * * * *  Member Functions   * * * * * * * * * * * * //

void Foam::morphingBoxConstraint::computeBounds
(
    autoPtr<scalarField>& lowerBounds,
    autoPtr<scalarField>& upperBounds
)
{
    if (lowerBounds || upperBounds)
    {
        NotImplemented;
    }
}


void Foam::morphingBoxConstraint::updateBounds
(
    autoPtr<scalarField>& lowerBounds,
    autoPtr<scalarField>& upperBounds
)
{
    if (designVariables_.updateBounds() && (lowerBounds || upperBounds))
    {
        NotImplemented;
    }
}


Foam::tmp<Foam::scalarField> Foam::morphingBoxConstraint::postProcessSens
(
    const scalarField& controlPointSens,
    const word& adjointSolverName
)
{
    // Sensitivities w.r.t. the design variables
    auto tdvSens
        (tmp<scalarField>::New(designVariables_.scalarField::size(), Zero));
    scalarField& dvSens = tdvSens.ref();
    computeDVsSensitivities(dvSens, controlPointSens);
    writeDVSensitivities(dvSens, adjointSolverName);

    return tdvSens;
}


Foam::scalar Foam::morphingBoxConstraint::computeEta
(
    scalarField& correction,
    const scalar maxInitChange
)
{
    vectorField cpMovement(designVariables_.controlPointMovement(correction));
    const scalar maxDisplacement
    (
        volBSplinesBase_.computeMaxBoundaryDisplacement
        (
            cpMovement,
            designVariables_.getPatches().toc()
        )
    );

    Info<< "maxAllowedDisplacement/maxDisplacement of boundary\t"
        << maxInitChange << "/" << maxDisplacement << endl;
    const scalar eta(maxInitChange/ maxDisplacement);

    Info<< "Setting eta value to " << eta << endl;
    correction *= eta;

    return eta;
}


bool Foam::morphingBoxConstraint::writeData(Ostream& os) const
{
    return true;
}


// ************************************************************************* //
