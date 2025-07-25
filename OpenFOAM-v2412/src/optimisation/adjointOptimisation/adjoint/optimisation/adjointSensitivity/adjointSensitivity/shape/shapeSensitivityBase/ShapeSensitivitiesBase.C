/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2023 PCOpt/NTUA
    Copyright (C) 2013-2023 FOSS GP
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

#include "HashSet.H"
#include "ShapeSensitivitiesBase.H"
#include "adjointSensitivity.H"
#include "adjointSolver.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ShapeSensitivitiesBase, 0);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::ShapeSensitivitiesBase::allocateEikonalSolver()
{
    // Allocate distance solver if needed
    if (includeDistance_ && !eikonalSolver_)
    {
        eikonalSolver_.reset
        (
            new adjointEikonalSolver
            (
                mesh_,
                this->dict(),
                adjointSolver_,
                geometryVariationIntegrationPatches()
            )
        );
    }
}


bool Foam::ShapeSensitivitiesBase::hasMultiplier
(
    bool (objective::*hasFunction)() const
)
{
    bool hasMult(false);
    const PtrList<objective>& objectives =
        adjointSolver_.getObjectiveManager().getObjectiveFunctions();
    for (const objective& func : objectives)
    {
        hasMult = hasMult || (func.*hasFunction)();
    }
    return hasMult;
}


void Foam::ShapeSensitivitiesBase::writeFaceBasedSens() const
{
    // Wall face sensitivity projected to normal
    if (wallFaceSensNormalPtr_)
    {
        constructAndWriteSensitivityField<scalar>
        (
            wallFaceSensNormalPtr_,
            "faceSensNormal" + suffix_
        );
    }

    if (writeAllSurfaceFiles_)
    {
        // Wall face sensitivity vectors
        if (wallFaceSensVecPtr_)
        {
            constructAndWriteSensitivityField<vector>
            (
                wallFaceSensVecPtr_,
                "faceSensVec" + suffix_
            );
        }

        // Normal sens as vectors
        if (wallFaceSensNormalVecPtr_)
        {
            constructAndWriteSensitivityField<vector>
            (
                wallFaceSensNormalVecPtr_,
                "faceSensNormalVec" + suffix_
            );
        }
    }
}


void Foam::ShapeSensitivitiesBase::writePointBasedSens() const
{
    // Wall point sensitivity projected to normal
    if (wallPointSensNormalPtr_)
    {
        constructAndWriteSensitivtyPointField<scalar>
        (
            wallPointSensNormalPtr_,
            "pointSensNormal" + suffix_
        );
    }

    // Write point-based sensitivities, if present
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (writeAllSurfaceFiles_)
    {
        // Wall point sensitivity vectors
        if (wallPointSensVecPtr_)
        {
            constructAndWriteSensitivtyPointField<vector>
            (
                wallPointSensVecPtr_,
                "pointSensVec" + suffix_
            );
        }

        // Normal point as vectors
        if (wallPointSensNormalVecPtr_)
        {
            constructAndWriteSensitivtyPointField<vector>
            (
                wallPointSensNormalVecPtr_,
                "pointSensNormalVec" + suffix_
            );
        }
    }
}


void Foam::ShapeSensitivitiesBase::clearSurfaceFields()
{
    // Face-based boundary sens
    if (wallFaceSensVecPtr_)
    {
        wallFaceSensVecPtr_() = vector::zero;
    }
    if (wallFaceSensNormalVecPtr_)
    {
        wallFaceSensNormalVecPtr_() = vector::zero;
    }
    if (wallFaceSensNormalPtr_)
    {
        wallFaceSensNormalPtr_() = scalar(0);
    }

    // Point-based boundary sens
    if (wallPointSensVecPtr_)
    {
        for (vectorField& patchSens : wallPointSensVecPtr_())
        {
            patchSens = vector::zero;
        }
    }
    if (wallPointSensNormalVecPtr_)
    {
        for (vectorField& patchSens : wallPointSensNormalVecPtr_())
        {
            patchSens = vector::zero;
        }
    }
    if (wallPointSensNormalPtr_)
    {
        for (scalarField& patchSens : wallPointSensNormalPtr_())
        {
            patchSens = scalar(0);
        }
    }
}


void Foam::ShapeSensitivitiesBase::allocateMultipliers()
{
    gradDxDbMult_.reset
    (
        new volTensorField
        (
            IOobject
            (
                "gradDxDbMult",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedTensor(sqr(dimLength)/pow3(dimTime), Zero)
        )
    );
    if (hasMultiplier(&objective::hasDivDxDbMult))
    {
        divDxDbMult_.reset(new scalarField(mesh_.nCells(), Zero));
    }
    if (hasMultiplier(&objective::hasdSdbMult))
    {
        dSfdbMult_.reset(createZeroBoundaryPtr<vector>(mesh_));
    }
    if (hasMultiplier(&objective::hasdndbMult))
    {
        dnfdbMult_.reset(createZeroBoundaryPtr<vector>(mesh_));
    }
    if (hasMultiplier(&objective::hasdxdbDirectMult))
    {
        dxdbDirectMult_.reset(createZeroBoundaryPtr<vector>(mesh_));
    }
    bcDxDbMult_.reset(createZeroBoundaryPtr<vector>(mesh_));
    optionsDxDbMult_.reset(new vectorField(mesh_.nCells(), Zero));
}


void Foam::ShapeSensitivitiesBase::clearMultipliers()
{
    gradDxDbMult_() = dimensionedTensor(gradDxDbMult_().dimensions(), Zero);
    if (divDxDbMult_)
    {
        divDxDbMult_() = Zero;
    }
    if (eikonalSolver_)
    {
        eikonalSolver_->reset();
    }
    if (dxdbMult_)
    {
        dxdbMult_() = Zero;
    }
    if (dSfdbMult_)
    {
        dSfdbMult_() = Zero;
    }
    if (dnfdbMult_)
    {
        dnfdbMult_() = Zero;
    }
    if (dxdbDirectMult_)
    {
        dxdbDirectMult_() = Zero;
    }
    if (pointDxDbDirectMult_)
    {
        for (vectorField& field : pointDxDbDirectMult_())
        {
            field = Zero;
        }
    }
    bcDxDbMult_() = Zero;
    optionsDxDbMult_() = Zero;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ShapeSensitivitiesBase::ShapeSensitivitiesBase
(
    const fvMesh& mesh,
    const dictionary& dict,
    class adjointSolver& adjointSolver
)
:
    adjointSensitivity(mesh, dict, adjointSolver),
    sensitivityPatchIDs_
    (
        mesh.boundaryMesh().patchSet
        (
            dict.optionalSubDict(mesh.name()).
                get<wordRes>("patches", keyType::REGEX_RECURSIVE)
        )
    ),
    writeAllSurfaceFiles_
    (
        dict.getOrDefault<bool>("writeAllSurfaceFiles", false)
    ),
    wallFaceSensVecPtr_(nullptr),
    wallFaceSensNormalPtr_(nullptr),
    wallFaceSensNormalVecPtr_(nullptr),

    wallPointSensVecPtr_(nullptr),
    wallPointSensNormalPtr_(nullptr),
    wallPointSensNormalVecPtr_(nullptr)
{
    allocateEikonalSolver();
    allocateMultipliers();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::ShapeSensitivitiesBase::readDict(const dictionary& dict)
{
    if (adjointSensitivity::readDict(dict))
    {
        sensitivityPatchIDs_ =
            mesh_.boundaryMesh().patchSet
            (
                dict_.optionalSubDict(mesh_.name()).
                    get<wordRes>("patches", keyType::REGEX_RECURSIVE)
            );
        writeAllSurfaceFiles_ =
            dict_.getOrDefault<bool>("writeAllSurfaceFiles", false);

        if (includeDistance_)
        {
            if (eikonalSolver_)
            {
                eikonalSolver_().readDict(dict);
            }
            else
            {
                allocateEikonalSolver();
            }
        }

        return true;
    }

    return false;
}


const Foam::labelHashSet&
Foam::ShapeSensitivitiesBase::geometryVariationIntegrationPatches() const
{
    return sensitivityPatchIDs_;
}


void Foam::ShapeSensitivitiesBase::accumulateIntegrand(const scalar dt)
{
    // Accumulate multiplier of grad(dxdb)
    adjointSolver_.accumulateGradDxDbMultiplier(gradDxDbMult_(), dt);

    // Accumulate multiplier of div(dxdb)
    adjointSolver_.accumulateDivDxDbMultiplier(divDxDbMult_, dt);

    // Terms from fvOptions - missing contributions from turbulence models
    adjointSolver_.accumulateOptionsDxDbMultiplier(optionsDxDbMult_(), dt);

    // Accumulate source for the adjoint to the eikonal equation
    if (eikonalSolver_)
    {
        eikonalSolver_->accumulateIntegrand(dt);
    }

    // Accumulate direct sensitivities
    adjointSolver_.accumulateGeometryVariationsMultipliers
    (
        dSfdbMult_,
        dnfdbMult_,
        dxdbDirectMult_,
        pointDxDbDirectMult_,
        geometryVariationIntegrationPatches(),
        dt
    );

    // Accumulate sensitivities due to boundary conditions
    adjointSolver_.accumulateBCSensitivityIntegrand
        (bcDxDbMult_, geometryVariationIntegrationPatches(), dt);
}


void Foam::ShapeSensitivitiesBase::clearSensitivities()
{
    adjointSensitivity::clearSensitivities();
    clearSurfaceFields();
    clearMultipliers();
}


void Foam::ShapeSensitivitiesBase::write(const word& baseName)
{
    adjointSensitivity::write(baseName);
    writeFaceBasedSens();
    writePointBasedSens();
}


Foam::tmp<Foam::volVectorField>
Foam::ShapeSensitivitiesBase::getWallFaceSensVec()
{
    if (wallFaceSensVecPtr_)
    {
        return
            constructVolSensitivtyField<vector>
            (
                wallFaceSensVecPtr_,
                "faceSensVec" + suffix_
            );
    }
    else
    {
        WarningInFunction
            << " no faceSensVec boundary field. Returning zero" << endl;

        return
            tmp<volVectorField>
            (
                createZeroFieldPtr<vector>
                (
                    mesh_,
                    "faceSensVec" + suffix_,
                    dimless
                ).ptr()
            );
    }
}


Foam::tmp<Foam::volScalarField>
Foam::ShapeSensitivitiesBase::getWallFaceSensNormal()
{
    if (wallFaceSensNormalPtr_)
    {
        return
            constructVolSensitivtyField<scalar>
            (
                wallFaceSensNormalPtr_,
                "faceSensNormal" + suffix_
            );
    }
    else
    {
        WarningInFunction
            << " no wallFaceSensNormal boundary field. Returning zero" << endl;

        return
            tmp<volScalarField>
            (
                createZeroFieldPtr<scalar>
                (
                    mesh_,
                    "faceSensNormal" + suffix_, dimless
                ).ptr()
            );
    }
}


Foam::tmp<Foam::volVectorField>
Foam::ShapeSensitivitiesBase::getWallFaceSensNormalVec()
{
    if (wallFaceSensNormalVecPtr_)
    {
        return
            constructVolSensitivtyField<vector>
            (
                wallFaceSensNormalVecPtr_,
                "faceSensNormalVec" + suffix_
            );
    }
    else
    {
        WarningInFunction
            << " no wallFaceSensNormalVec boundary field. Returning zero"
            << endl;

        return
            tmp<volVectorField>
            (
                createZeroFieldPtr<vector>
                (
                    mesh_,
                    "faceSensNormalVec" + suffix_,
                    dimless
                ).ptr()
            );
    }
}


Foam::tmp<Foam::pointVectorField>
Foam::ShapeSensitivitiesBase::getWallPointSensVec()
{
    tmp<volVectorField> tWallFaceSensVec = getWallFaceSensVec();
    volPointInterpolation volPointInter(mesh_);

    return (volPointInter.interpolate(tWallFaceSensVec));
}


Foam::tmp<Foam::pointScalarField>
Foam::ShapeSensitivitiesBase::getWallPointSensNormal()
{
    tmp<volScalarField> tWallFaceSensNormal = getWallFaceSensNormal();
    volPointInterpolation volPointInter(mesh_);

    return (volPointInter.interpolate(tWallFaceSensNormal));
}


Foam::tmp<Foam::pointVectorField>
Foam::ShapeSensitivitiesBase::getWallPointSensNormalVec()
{
    tmp<volVectorField> tWallFaceSensNormalVec = getWallFaceSensNormalVec();
    volPointInterpolation volPointInter(mesh_);

    return (volPointInter.interpolate(tWallFaceSensNormalVec));
}


const Foam::boundaryVectorField&
Foam::ShapeSensitivitiesBase::getWallFaceSensVecBoundary() const
{
    return wallFaceSensVecPtr_();
}


const Foam::boundaryScalarField&
Foam::ShapeSensitivitiesBase::getWallFaceSensNormalBoundary() const
{
    return wallFaceSensNormalPtr_();
}


const Foam::boundaryVectorField&
Foam::ShapeSensitivitiesBase::getWallFaceSensNormalVecBoundary() const
{
    return wallFaceSensNormalVecPtr_();
}


// ************************************************************************* //
