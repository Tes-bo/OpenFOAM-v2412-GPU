/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2023 PCOpt/NTUA
    Copyright (C) 2013-2023 FOSS GP
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

#include "localIOdictionary.H"
#include "topODesignVariables.H"
#include "wallDist.H"
#include "wallFvPatch.H"
#include "cutFaceIso.H"
#include "cutCellIso.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(topODesignVariables, 0);
    addToRunTimeSelectionTable
    (
        designVariables,
        topODesignVariables,
        designVariables
    );
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::topODesignVariables::updateField
(
    const scalarField& correction,
    const label fluidID
)
{
    DebugInfo
        << "Updating design variables for field " << fluidID << endl;
    const label n = mesh_.nCells();
    SubField<scalar> localCorrection(correction, n, fluidID*n);
    SubField<scalar> field(*this, n, fluidID*n);

    // Update porosity in adjoint porous cells
    if (zones_.adjointPorousZoneIDs().empty())
    {
        forAll(field, cellI)
        {
            field[cellI] +=
                min
                (
                    max
                    (
                        field[cellI] + localCorrection[cellI],
                        scalar(0)
                    ),
                    1.
                )
              - field[cellI];
        }
    }
    else
    {
        for (label cellZoneID : zones_.adjointPorousZoneIDs())
        {
            const labelList& zoneCells = mesh_.cellZones()[cellZoneID];
            for (label cellI : zoneCells)
            {
                field[cellI] +=
                    min
                    (
                        max
                        (
                            field[cellI] + localCorrection[cellI],
                            scalar(0)
                        ),
                        1.
                    )
                  - field[cellI];
            }
        }
    }
}


void Foam::topODesignVariables::applyFixedValues()
{
    SubField<scalar> alpha(*this, mesh_.nCells());
    // Zero porosity in the cells next to IO
    for (label cellI : zones_.IOCells())
    {
        alpha[cellI] = 0.;
    }

    // Apply fixed porosity
    forAll(zones_.fixedPorousZoneIDs(), zI)
    {
        const label cellZoneID = zones_.fixedPorousZoneIDs()[zI];
        const labelList& zoneCells = mesh_.cellZones()[cellZoneID];
        const scalar alphaValue(zones_.fixedPorousValues()[zI]);
        for (label cellI : zoneCells)
        {
            alpha[cellI] = alphaValue;
        }
    }

    // Apply fixed zero porosity
    for (label cellZoneID : zones_.fixedZeroPorousZoneIDs())
    {
        const labelList& zoneCells = mesh_.cellZones()[cellZoneID];
        for (label cellI : zoneCells)
        {
            alpha[cellI] = 0.;
        }
    }
}


Foam::scalar Foam::topODesignVariables::computeEta(scalarField& correction)
{
    const scalar maxChange(gMax(mag(correction)));
    Info<< "maxInitChange/maxChange \t"
        << maxInitChange_() << "/" << maxChange << endl;
    const scalar eta(maxInitChange_() / maxChange);
    Info<< "Setting eta value to " << eta << endl;
    correction *= eta;

    return eta;
}


void Foam::topODesignVariables::setActiveDesignVariables
(
    const label fluidID,
    const bool activeIO
)
{
    const label offset(fluidID*mesh_.nCells());
    label varI(activeDesignVariables_.size());
    activeDesignVariables_.setSize(offset + mesh_.nCells(), -1);
    // Set active design variables
    // If specific porosity zones are prescribed, use them directly
    if (!zones_.adjointPorousZoneIDs().empty())
    {
        for (label cellZoneID : zones_.adjointPorousZoneIDs())
        {
            for (const label var : mesh_.cellZones()[cellZoneID])
            {
                activeDesignVariables_[varI++] = var + offset;
            }
        }
    }
    // Else, pick up all cells in non-constant porosity zones
    else
    {
        boolList isActiveDV(mesh_.nCells(), true);
        // Exclude cells with fixed porosity
        for (label cellZoneID : zones_.fixedPorousZoneIDs())
        {
            for (label cellI : mesh_.cellZones()[cellZoneID])
            {
                isActiveDV[cellI] = false;
            }
        }
        for (label cellZoneID : zones_.fixedZeroPorousZoneIDs())
        {
            for (label cellI : mesh_.cellZones()[cellZoneID])
            {
                isActiveDV[cellI] = false;
            }
        }
        if (!activeIO)
        {
            for (label cellI : zones_.IOCells())
            {
                isActiveDV[cellI] = false;
            }
        }

        // Set active design variables
        forAll(isActiveDV, cellI)
        {
            if (isActiveDV[cellI])
            {
                activeDesignVariables_[varI++] = offset + cellI;
            }
        }
    }
    activeDesignVariables_.setSize(varI);
}


void Foam::topODesignVariables::readField
(
    const word& name,
    const label fluidID,
    const bool setIOValues
)
{
    const label offset(fluidID*mesh_.nCells());
    if (localIOdictionary::found(name))
    {
        SubField<scalar>(*this, mesh_.nCells(), offset) =
            scalarField(name, *this, mesh_.nCells());
    }
    else
    {
        IOobject header
        (
            name,
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        );
        if (header.typeHeaderOk<volScalarField>())
        {
            Info<< "Setting design variables based on the alpha field "
                << nl << endl;
            volScalarField volField
            (
                header,
                mesh_
            );
            const scalarField& field = volField.primitiveField();
            forAll(field, cI)
            {
                scalarField::operator[](offset + cI) = field[cI];
            }
        }
    }
}


void Foam::topODesignVariables::initialize()
{
    // Set active design variables
    setActiveDesignVariables();

    // Read in values from file, if present
    readField("alpha", 0, true);

    if (regularisation_.growFromWalls())
    {
        scalarField& alpha = *this;
        for (const fvPatch& patch : mesh_.boundary())
        {
            if (isA<wallFvPatch>(patch))
            {
                const labelList& faceCells = patch.faceCells();
                forAll(faceCells, cI)
                {
                    alpha[faceCells[cI]] = 1.;
                }
            }
        }
    }

    // Make sure alpha has fixed values where it should
    scalarField dummyCorrection(mesh_.nCells(), Zero);
    update(dummyCorrection);

    // Read bounds for design variables, if present
    readBounds(autoPtr<scalar>::New(0), autoPtr<scalar>::New(1));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::topODesignVariables::topODesignVariables
(
    fvMesh& mesh,
    const dictionary& dict
)
:
    topODesignVariables(mesh, dict, mesh.nCells())
{}


Foam::topODesignVariables::topODesignVariables
(
    fvMesh& mesh,
    const dictionary& dict,
    const label size
)
:
    topOVariablesBase(mesh, dict),
    designVariables(mesh, dict, size),
    alpha_(SubField<scalar>(*this, mesh.nCells(), 0)),
    regularisation_
    (
        mesh,
        alpha_,
        zones_,
        dict_.subDict("regularisation")
    ),
    writeAllFields_
    (
        dict.getOrDefaultCompat<bool>
        (
            "writeAllFields", {{"writeAllAlphaFields", 2306}}, false
        )
    ),
    addFvOptions_(dict.getOrDefault<bool>("addFvOptions", false))
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::topODesignVariables> Foam::topODesignVariables::New
(
    fvMesh& mesh,
    const dictionary& dict
)
{
    return autoPtr<topODesignVariables>::New(mesh, dict);
}


// * * * * * * * * * * * * * * *  Member Functions   * * * * * * * * * * * * //

const Foam::volScalarField& Foam::topODesignVariables::beta() const
{
    return regularisation_.beta();
}


const Foam::scalarField& Foam::topODesignVariables::interpolationField
(
    const word& interpolationFieldName
) const
{
    return beta().primitiveField();
}


void Foam::topODesignVariables::interpolate
(
    volScalarField& field,
    const topOInterpolationFunction& interpolationFunc,
    const FieldField<Field, scalar>& fluidValues,
    const scalarField& solidValues,
    const label fieldi,
    const word& interpolationFieldName
) const
{
    const scalarField& indicator = interpolationField(interpolationFieldName);
    scalarField interpolant(indicator.size(), Zero);
    interpolationFunc.interpolate(indicator, interpolant);

    // Interpolate field values
    const scalar diff(solidValues[fieldi] - fluidValues[0][fieldi]);
    field.primitiveFieldRef() = fluidValues[0][fieldi] + interpolant*diff;
    field.correctBoundaryConditions();
}


void Foam::topODesignVariables::interpolationSensitivities
(
    scalarField& sens,
    const topOInterpolationFunction& interpolationFunc,
    const FieldField<Field, scalar>& fluidValues,
    const scalarField& solidValues,
    const label fieldi,
    const word& designVariablesName,
    const word& interpolationFieldName
) const
{
    const scalarField& indicator = interpolationField(interpolationFieldName);
    sens *=
        (solidValues[fieldi] - fluidValues[0][fieldi])
       *interpolationFunc.derivative(indicator);
}


void Foam::topODesignVariables::nullifyInSolid
(
    scalarField& field,
    const topOInterpolationFunction& interpolationFunc
) const
{
    const scalarField& beta = this->beta().primitiveField();
    scalarField interpolant(beta.size(), Zero);
    interpolationFunc.interpolate(beta, interpolant);
    field *= scalar(1) - interpolant;
}


void Foam::topODesignVariables::nullifyInSolidSensitivities
(
    scalarField& sens,
    const topOInterpolationFunction& interpolationFunc,
    const word& designVariablesName
) const
{
    const scalarField& beta = this->beta().primitiveField();
    sens *= - interpolationFunc.derivative(beta);
}


Foam::tmp<Foam::scalarField> Foam::topODesignVariables::penalty
(
    const word& interpolationFieldName,
    const topOInterpolationFunction& interpolationFunc
) const
{
    const scalarField& beta = this->beta().primitiveField();
    tmp<scalarField> tinterpolant(tmp<scalarField>::New(beta.size(), Zero));
    interpolationFunc.interpolate(beta, tinterpolant.ref());
    return tinterpolant;
}


Foam::tmp<Foam::scalarField> Foam::topODesignVariables::penaltySensitivities
(
    const word& interpolationFieldName,
    const topOInterpolationFunction& interpolationFunc
) const
{
    const scalarField& beta = this->beta().primitiveField();
    return interpolationFunc.derivative(beta);
}


void Foam::topODesignVariables::update(scalarField& correction)
{
    // Update alpha values
    updateField(correction);

    // Fix alpha in zones
    applyFixedValues();

    // Update the beta field
    regularisation_.updateBeta();

    // Though the mesh is kept constant, the distance from wall may change
    // if the method computing it includes fvOptions that depend on the
    // indicator field.
    // Trick wallDist into updating it

    wallDist::try_movePoints(mesh_);


    // Write the 0.5 beta iso-line to files, as an indication of the
    // fluid-solid interface
    labelList changedFaces(mesh_.nFaces(), -1);
    List<wallPointData<label>> changedFacesInfo(mesh_.nFaces());
    writeFluidSolidInterface(-beta(), -0.5, changedFaces, changedFacesInfo);
}


bool Foam::topODesignVariables::globalSum() const
{
    return true;
}


Foam::tmp<Foam::scalarField> Foam::topODesignVariables::assembleSensitivities
(
    adjointSensitivity& adjointSens
)
{
    // Raw sensitivities field.
    // Does not include the adjoint to the regularisation and projection steps
    const scalarField& fieldSens = adjointSens.fieldSensPtr()->primitiveField();

    // Return field (complete sensitivities)
    auto tobjectiveSens(tmp<scalarField>::New(fieldSens));
    scalarField& objectiveSens = tobjectiveSens.ref();

    // Add part due to regularisation and projection
    regularisation_.postProcessSens(objectiveSens);

    // Write final sensitivities field
    if (writeAllFields_ && mesh_.time().writeTime())
    {
        volScalarField sens
        (
            IOobject
            (
                "topOSens" + adjointSens.getAdjointSolver().solverName(),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar(dimless, Zero)
        );
        sens.primitiveFieldRef() = objectiveSens;
        sens.write();
    }

    return tobjectiveSens;
}


void Foam::topODesignVariables::setInitialValues()
{
    // Rest of the contrsuctor initialisation
    initialize();
}


void Foam::topODesignVariables::addFvOptions
(
    const PtrList<primalSolver>& primalSolvers,
    const PtrList<adjointSolverManager>& adjointSolverManagers
)
{
    // WIP
    if (addFvOptions_)
    {
        for (const primalSolver& solver : primalSolvers)
        {
            solver.addTopOFvOptions();
        }
        for (const adjointSolverManager& manager : adjointSolverManagers)
        {
            const PtrList<adjointSolver>& adjointSolvers =
                manager.adjointSolvers();
            for (const adjointSolver& solver : adjointSolvers)
            {
                solver.addTopOFvOptions();
            }
        }
    }
}


void Foam::topODesignVariables::writeDesignVars()
{
    if (writeAllFields_ && mesh_.time().writeTime())
    {
        volScalarField alpha
        (
            IOobject
            (
                "alpha",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar(dimless/dimTime, Zero)
        );
        alpha.primitiveFieldRef() = alpha_;

        alpha.write();
    }
}


bool Foam::topODesignVariables::writeData(Ostream& os) const
{
    const scalarField& alpha = alpha_;
    alpha.writeEntry("alpha", os);

    return true;
}


// ************************************************************************* //
