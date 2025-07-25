/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022-2023 PCOpt/NTUA
    Copyright (C) 2022-2023 FOSS GP
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

#include "levelSetDesignVariables.H"
#include "wallDist.H"
#include "zeroGradientFvPatchField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * * //

defineTypeNameAndDebug(levelSetDesignVariables, 1);
addToRunTimeSelectionTable
(
    designVariables,
    levelSetDesignVariables,
    designVariables
);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void levelSetDesignVariables::readField()
{
    scalarField& vars = *this;
    if (localIOdictionary::found("alpha"))
    {
        vars = (scalarField("alpha", *this, vars.size()));
    }
    else
    {
        // Initialise as the distance from the wall patches of the initial
        // domain
        const labelHashSet wallPatchIDs =
            mesh_.boundaryMesh().findPatchIDs<wallPolyPatch>();
        volScalarField y
        (
            IOobject
            (
                "yLevelSet",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar(dimLength, Zero),
            patchDistMethod::patchTypes<scalar>(mesh_, wallPatchIDs)
        );
        patchDistMethod::New
        (
            dict_.subDict("initialisation"),
            mesh_,
            wallPatchIDs
        )->correct(y);
        vars = y.primitiveField();

        if (debug)
        {
            writeDesignVars();
        }
    }
}


void levelSetDesignVariables::applyFixedPorosityValues()
{
    scalarField& betaIf = beta_.primitiveFieldRef();

    // Safety - set beta equal to zero next to the IO cells
    for (const label IOcell : zones_.IOCells())
    {
        betaIf[IOcell] = Zero;
    }

    const labelList& fixedZeroPorousZones = zones_.fixedZeroPorousZoneIDs();
    const labelList& fixedPorousZones = zones_.fixedPorousZoneIDs();
    const scalarList& fixedPorousValues = zones_.fixedPorousValues();

    // Apply fixed porosity
    for (const label zoneID : fixedZeroPorousZones)
    {
        const labelList& zone = mesh_.cellZones()[zoneID];
        for (const label cellI : zone)
        {
            betaIf[cellI] = Zero;
        }
    }

    // Apply fixed porosity
    forAll(fixedPorousZones, zI)
    {
        const label zoneID = fixedPorousZones[zI];
        const scalar value = fixedPorousValues[zI];
        const labelList& zone = mesh_.cellZones()[zoneID];
        for (const label cellI : zone)
        {
            betaIf[cellI] = value >= 0 ? 0 : 1;
        }
    }

    beta_.correctBoundaryConditions();
}


void levelSetDesignVariables::setActiveDesignVariables(bool activeIO)
{
    activeDesignVariables_.setSize(mesh_.nCells(), -1);
    if (zones_.adjointPorousZoneIDs().empty())
    {
        boolList isActiveVar(mesh_.nCells(), true);

        const labelList& fixedZeroPorousZones =
            zones_.fixedZeroPorousZoneIDs();
        for (const label zoneID : fixedZeroPorousZones)
        {
            const labelList& zone = mesh_.cellZones()[zoneID];
            for (const label cellI : zone)
            {
                isActiveVar[cellI] = false;
            }
        }

        const labelList& fixedPorousZones = zones_.fixedPorousZoneIDs();
        for (const label zoneID : fixedPorousZones)
        {
            const labelList& zone = mesh_.cellZones()[zoneID];
            for (const label cellI : zone)
            {
                isActiveVar[cellI] = false;
            }
        }

        if (!activeIO)
        {
            for (label cellI : zones_.IOCells())
            {
                isActiveVar[cellI] = false;
            }
        }

        label iVar(0);
        forAll(isActiveVar, cI)
        {
            if (isActiveVar[cI])
            {
                activeDesignVariables_[iVar++] = cI;
            }
        }
        activeDesignVariables_.setSize(iVar);
    }
    else
    {
        const labelList& adjointPorousZoneIDs = zones_.adjointPorousZoneIDs();

        label iVar(0);
        for (const label cellZoneID : adjointPorousZoneIDs)
        {
            const labelList& zone = mesh_.cellZones()[cellZoneID];

            for (const label cellI : zone)
            {
                activeDesignVariables_[iVar] = cellI;
                iVar++;
            }
        }
        activeDesignVariables_.setSize(iVar);
    }
}


void levelSetDesignVariables::updateBeta()
{
    // Compute the beta field by passing the distance field through
    // a Heaviside function
    scalarField& beta = beta_.primitiveFieldRef();
    interpolation_->interpolate(signedDistances_.primitiveField(), beta);
    beta = 1 - beta;
    // Apply fixed values if necessary
    applyFixedPorosityValues();

    beta_.correctBoundaryConditions();
}


void Foam::levelSetDesignVariables::updateSignedDistances()
{
    Info<< "Re-initilising the level-set distance field" << nl << endl;

    volScalarField y
    (
        IOobject
        (
            "yLevelSet",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimLength, Zero),
        fvPatchFieldBase::zeroGradientType()
    );
    y.primitiveFieldRef() = aTilda_.primitiveFieldRef();
    y.correctBoundaryConditions();

    changedFaces_.clear();
    changedFaces_.setSize(mesh_.nFaces(), -1);

    changedFacesInfo_.clear();
    changedFacesInfo_.setSize(mesh_.nFaces());

    writeFluidSolidInterface(aTilda_, 0, changedFaces_, changedFacesInfo_);

    List<wallPointData<label>> allFaceInfo(mesh_.nFaces());
    allCellInfo_.clear();
    allCellInfo_.setSize(mesh_.nCells());

    FaceCellWave<wallPointData<label>> wave
    (
        mesh_,
        changedFaces_,
        changedFacesInfo_,
        allFaceInfo,
        allCellInfo_,
        mesh_.globalData().nTotalCells() + 1
    );

    // Transfer the distance from cellInfo to the alphaTilda field
    forAll(allCellInfo_, celli)
    {
        if (allCellInfo_[celli].valid(wave.data()))
        {
            signedDistances_[celli] =
                sign(aTilda_[celli])*Foam::sqrt(allCellInfo_[celli].distSqr());
        }
    }
    signedDistances_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

levelSetDesignVariables::levelSetDesignVariables
(
    fvMesh& mesh,
    const dictionary& dict
)
:
    topOVariablesBase(mesh, dict),
    designVariables(mesh, dict, mesh.nCells()),
    radius_
        (regularisationRadius::New(mesh, dict.subDict("regularisation"), false)),
    regularisation_
        (regularisationPDE::New(mesh, dict.subDict("regularisation"), zones_)),
    aTilda_
    (
        IOobject
        (
            "aTilda",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero),
        fvPatchFieldBase::zeroGradientType()
    ),
    signedDistances_
    (
        IOobject
        (
            "signedDistances",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero),
        fvPatchFieldBase::zeroGradientType()
    ),
    interpolation_
    (
        topOInterpolationFunction::New(mesh_, dict_.subDict("interpolation"))
    ),
    beta_
    (
        IOobject
        (
            "beta",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    ),
    fixATildaValues_(dict.getOrDefault<bool>("fixATildaValues", true)),
    writeAllDistanceFields_
        (dict.getOrDefault<bool>("writeAllDistanceFields", false)),
    changedFaces_(),
    changedFacesInfo_(),
    allCellInfo_()
{
    // Read the alpha field if present, or set it based on the distance field
    readField();

    // Read bounds of design variables if present.
    // If not, use the maximum distance field in the fluid to set them.
    scalar maxDist = gMax(*this);
    scalar lowerBound =
        localIOdictionary::getOrDefault<scalar>("lowerBound", -maxDist - SMALL);
    scalar upperBound =
        localIOdictionary::getOrDefault<scalar>("upperBound", maxDist + SMALL);
    readBounds
    (
        autoPtr<scalar>(new scalar(lowerBound)),
        autoPtr<scalar>(new scalar(upperBound))
    );
    DebugInfo
        << "Using lower/upper bounds "
        << lowerBounds_()[0] << "/" << upperBounds_()[0]
        << endl;

    // Update the beta field based on the initial design vars
    scalarField zeroUpdate(scalarField::size(), Zero);
    update(zeroUpdate);

    // Determine which design variables are active
    setActiveDesignVariables();
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<levelSetDesignVariables> levelSetDesignVariables::New
(
    fvMesh& mesh,
    const dictionary& dict
)
{
    return autoPtr<levelSetDesignVariables>
    (
        new levelSetDesignVariables(mesh, dict)
    );
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const volScalarField& levelSetDesignVariables::beta() const
{
    return beta_;
}


void levelSetDesignVariables::update(scalarField& correction)
{
    scalarField::operator+=(correction);

    // Compute the regularised design variables
    regularisation_->regularise
    (
        aTilda_, *this, aTilda_.primitiveFieldRef(),
        true, radius_(), upperBounds_()[0], fixATildaValues_
    );
    aTilda_.correctBoundaryConditions();

    if (writeAllDistanceFields_)
    {
        writeDesignVars();
        aTilda_.write();
    }

    // Compute signed distances based on aTilda
    updateSignedDistances();

    // Set beta based on the signed distances
    updateBeta();

    if (writeAllDistanceFields_)
    {
        signedDistances_.write();
        beta_.write();
    }

    // Though the mesh is kept constant, the distance from wall may change
    // due to fvOptions depending on beta. Trick wallDist into updating it

    wallDist::try_movePoints(mesh_);
}


scalar levelSetDesignVariables::computeEta(scalarField& correction)
{
    const scalar maxChange(gMax(mag(correction)));
    Info<< "maxInitChange/maxChange \t"
        << maxInitChange_() << "/" << maxChange << endl;
    const scalar eta(maxInitChange_() / maxChange);
    Info<< "Setting eta value to " << eta << endl;
    correction *= eta;

    return eta;
}


bool levelSetDesignVariables::globalSum() const
{
    return true;
}


tmp<scalarField> levelSetDesignVariables::assembleSensitivities
(
    adjointSensitivity& adjointSens
)
{
    // Raw sensitivities field
    const scalarField& fieldSens = adjointSens.fieldSensPtr()->primitiveField();

    // Return field
    auto tobjectiveSens(tmp<scalarField>::New(fieldSens));
    scalarField& objectiveSens = tobjectiveSens.ref();

    // Multiply with dBetadAtilda
    objectiveSens *=
        -interpolation_->derivative(signedDistances_.primitiveField());

    // Solve the adjoint to the regularisation equation
    regularisation_->
        regularise(aTilda_, objectiveSens, objectiveSens, false, radius_());

    objectiveSens *= mesh_.V();

    if (writeAllDistanceFields_)
    {
        volScalarField sens
        (
            IOobject
            (
                "sens" + adjointSens.getAdjointSolver().solverName(),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar(dimless, Zero)
        );
        sens.primitiveFieldRef() = objectiveSens;
        sens.write();
    }

    return tobjectiveSens;
}


void levelSetDesignVariables::writeDesignVars()
{
    if (writeAllDistanceFields_ || mesh_.time().writeTime())
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
            dimensionedScalar(dimLength, Zero)
        );
        alpha.primitiveFieldRef() = *this;
        alpha.correctBoundaryConditions();

        alpha.write();
    }
}


bool levelSetDesignVariables::writeData(Ostream& os) const
{
    // Lower and upper bound values might be computed based on the initial
    // distance field. Write to file to enable continuation
    os.writeEntry("lowerBound", lowerBounds_()[0]);
    os.writeEntry("upperBound", upperBounds_()[0]);

    scalarField::writeEntry("alpha", os);

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
