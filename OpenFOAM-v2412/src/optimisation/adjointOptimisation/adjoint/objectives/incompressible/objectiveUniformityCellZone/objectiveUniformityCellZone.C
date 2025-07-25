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

#include "objectiveUniformityCellZone.H"
#include "createZeroField.H"
#include "IOmanip.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace objectives
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(objectiveUniformityCellZone, 0);
addToRunTimeSelectionTable
(
    objectiveIncompressible,
    objectiveUniformityCellZone,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

objectiveUniformityCellZone::objectiveUniformityCellZone
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& adjointSolverName,
    const word& primalSolverName
)
:
    objectiveIncompressible(mesh, dict, adjointSolverName, primalSolverName),
    zones_(mesh_.cellZones().indices(dict.get<wordRes>("zones"))),
    UMean_(zones_.size(), Zero),
    UVar_(zones_.size(), Zero),
    volZone_(zones_.size(), Zero)
{
    // Append Ua name to fieldNames
    fieldNames_.setSize
    (
        1,
        mesh_.lookupObject<solver>(adjointSolverName_).
            extendedVariableName("Ua")
    );

    // Check if cellZones provided include at least one cell
    checkCellZonesSize(zones_);

    // Allocate source term to the adjoint momentum equations
    dJdvPtr_.reset
    (
        createZeroFieldPtr<vector>
        (
            mesh_,
            "dJdv" + type(),
            dimLength/sqr(dimTime)
        )
    );
    // Allocate term to be added to volume-based sensitivity derivatives
    divDxDbMultPtr_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "divDxDbMult" + objectiveName_,
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar(sqr(dimLength)/pow3(dimTime), Zero),
            fvPatchFieldBase::zeroGradientType()
        )
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar objectiveUniformityCellZone::J()
{
    J_ = Zero;

    // References
    const volVectorField& U = vars_.UInst();
    const scalarField& V = mesh_.V().field();

    for (const label zI : zones_)
    {
        const cellZone& zoneI = mesh_.cellZones()[zI];
        scalarField VZone(V, zoneI);
        vectorField UZone(U.primitiveField(), zoneI);
        volZone_[zI] = gSum(VZone);
        UMean_[zI] = gSum(UZone*VZone)/volZone_[zI];
        UVar_[zI] = gSum(magSqr(UZone - UMean_[zI])*VZone)/volZone_[zI];

        J_ += 0.5*UVar_[zI];
    }

    return J_;
}


void objectiveUniformityCellZone::update_dJdv()
{
    const volVectorField& U = vars_.U();

    for (const label zI : zones_)
    {
        const cellZone& zoneI = mesh_.cellZones()[zI];
        for (const label cellI : zoneI)
        {
            dJdvPtr_()[cellI] = (U[cellI] - UMean_[zI])/volZone_[zI];
        }
    }
}


void objectiveUniformityCellZone::update_divDxDbMultiplier()
{
    volScalarField& divDxDbMult = divDxDbMultPtr_();
    const volVectorField& U = vars_.U();

    for (const label zI : zones_)
    {
        const cellZone& zoneI = mesh_.cellZones()[zI];
        for (const label cellI : zoneI)
        {
            divDxDbMult[cellI] =
                0.5*(magSqr(U[cellI] - UMean_[zI]) - UVar_[zI])/volZone_[zI];
        }
    }
    divDxDbMult.correctBoundaryConditions();
}


void objectiveUniformityCellZone::addHeaderColumns() const
{
    OFstream& file = objFunctionFilePtr_();
    for (const label zI : zones_)
    {
        const word zoneName(mesh_.cellZones()[zI].name());
        file<< setw(width_) << word(zoneName + "-" + "UMean") <<  " ";
        file<< setw(width_) << word(zoneName + "-" + "UVar") <<  " ";
        file<< setw(width_) << word(zoneName + "-" + "UStd") <<  " ";
        file<< setw(width_) << word(zoneName + "-" + "Vol") <<  " ";
    }
}


void objectiveUniformityCellZone::addColumnValues() const
{
    OFstream& file = objFunctionFilePtr_();
    forAll(UMean_, zI)
    {
        file<< setw(width_) << mag(UMean_[zI]) << " ";
        file<< setw(width_) << UVar_[zI] << " ";
        file<< setw(width_) << sqrt(UVar_[zI]) << " ";
        file<< setw(width_) << volZone_[zI] << " ";
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace objectives
} // End namespace Foam

// ************************************************************************* //
