/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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

#include "referenceTemperature.H"
#include "fvMesh.H"
#include "basicThermo.H"
#include "surfaceInterpolate.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace heatExchangerModels
{
    defineTypeNameAndDebug(referenceTemperature, 0);
    addToRunTimeSelectionTable
    (
        heatExchangerModel,
        referenceTemperature,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar
Foam::heatExchangerModels::referenceTemperature::primaryNetMassFlux() const
{
    const auto& phi = mesh_.lookupObject<surfaceScalarField>(phiName_);

    scalar sumPhi = 0;

    forAll(faceId_, i)
    {
        const label facei = faceId_[i];
        if (facePatchId_[i] != -1)
        {
            const label patchi = facePatchId_[i];
            sumPhi += phi.boundaryField()[patchi][facei]*faceSign_[i];
        }
        else
        {
            sumPhi += phi[facei]*faceSign_[i];
        }
    }
    reduce(sumPhi, sumOp<scalar>());

    return sumPhi;
}


Foam::scalar
Foam::heatExchangerModels::referenceTemperature::primaryInletTemperature() const
{
    const auto& phi = mesh_.lookupObject<surfaceScalarField>(phiName_);
    const auto& T = mesh_.lookupObject<volScalarField>(TName_);

    const surfaceScalarField Tf(fvc::interpolate(T));

    scalar sumMagPhi = 0;
    scalar primaryInletTfMean = 0;

    forAll(faceId_, i)
    {
        const label facei = faceId_[i];
        if (facePatchId_[i] != -1)
        {
            const label patchi = facePatchId_[i];
            const scalar phii = phi.boundaryField()[patchi][facei]*faceSign_[i];
            const scalar magPhii = mag(phii);

            sumMagPhi += magPhii;

            const scalar Tfi = Tf.boundaryField()[patchi][facei];

            primaryInletTfMean += Tfi*magPhii;
        }
        else
        {
            const scalar phii = phi[facei]*faceSign_[i];
            const scalar magPhii = mag(phii);

            sumMagPhi += magPhii;

            primaryInletTfMean += Tf[facei]*magPhii;
        }
    }
    reduce(sumMagPhi, sumOp<scalar>());
    reduce(primaryInletTfMean, sumOp<scalar>());

    return primaryInletTfMean/(sumMagPhi + ROOTVSMALL);
}


Foam::scalarField
Foam::heatExchangerModels::referenceTemperature::temperatureDifferences
(
    const labelList& cells
) const
{
    const auto& T = mesh_.lookupObject<volScalarField>(TName_);
    const scalarField TCells(T, cells);
    scalarField deltaTCells(cells.size(), Zero);

    if (Qt_ > 0)
    {
        forAll(deltaTCells, i)
        {
            deltaTCells[i] = max(Tref_ - TCells[i], scalar(0));
        }
    }
    else
    {
        forAll(deltaTCells, i)
        {
            deltaTCells[i] = max(TCells[i] - Tref_, scalar(0));
        }
    }

    return deltaTCells;
}


Foam::scalar
Foam::heatExchangerModels::referenceTemperature::weight
(
    const labelList& cells,
    const scalarField& deltaTCells
) const
{
    scalar sumWeight = 0;

    const auto& U = mesh_.lookupObject<volVectorField>(UName_);
    const scalarField& V = mesh_.V();

    forAll(cells, i)
    {
        const label celli = cells[i];
        sumWeight += V[celli]*mag(U[celli])*deltaTCells[i];
    }
    reduce(sumWeight, sumOp<scalar>());

    return sumWeight;
}


void Foam::heatExchangerModels::referenceTemperature::writeFileHeader
(
    Ostream& os
) const
{
    writeFile::writeHeader(os, "Effectiveness heat exchanger source");
    writeFile::writeCommented(os, "Time");
    writeFile::writeTabbed(os, "Net mass flux [kg/s]");
    writeFile::writeTabbed(os, "Total heat exchange [W]");
    writeFile::writeTabbed(os, "Reference T [K]");
    os  << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatExchangerModels::referenceTemperature::referenceTemperature
(
    const fvMesh& mesh,
    const word& name,
    const dictionary& coeffs
)
:
    heatExchangerModel(mesh, name, coeffs),
    targetQdotPtr_
    (
        Function1<scalar>::New
        (
            "targetQdot",
            coeffs,
            &mesh_
        )
    ),
    TrefTablePtr_(nullptr),
    sumPhi_(0),
    Qt_(0),
    Tref_(0)
{
    writeFileHeader(file());
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::heatExchangerModels::referenceTemperature::initialise()
{
    heatExchangerModel::initialise();
}


Foam::tmp<Foam::scalarField>
Foam::heatExchangerModels::referenceTemperature::energyDensity
(
    const labelList& cells
)
{
    sumPhi_ = primaryNetMassFlux();

    Qt_ = targetQdotPtr_->value(mesh_.time().value());

    if (TrefTablePtr_)
    {
        const scalar primaryInletT = primaryInletTemperature();
        Tref_ = TrefTablePtr_()(mag(sumPhi_), primaryInletT);
    }

    const scalarField deltaTCells(temperatureDifferences(cells));

    const scalar sumWeight = weight(cells, deltaTCells);

    return Qt_*deltaTCells/(sumWeight + ROOTVSMALL);
}


bool Foam::heatExchangerModels::referenceTemperature::read
(
    const dictionary& dict
)
{
    if (!writeFile::read(dict))
    {
        return false;
    }

    Info<< incrIndent << indent << "- using model: " << type() << endl;

    if (coeffs_.readIfPresent("Tref", Tref_))
    {
        Info<< indent << "- using constant reference temperature: " << Tref_
            << endl;
    }
    else
    {
        TrefTablePtr_.reset(new interpolation2DTable<scalar>(coeffs_));

        Info<< indent << "- using reference temperature table"
            << endl;
    }

    UName_ = coeffs_.getOrDefault<word>("U", "U");
    TName_ = coeffs_.getOrDefault<word>("T", "T");
    phiName_ = coeffs_.getOrDefault<word>("phi", "phi");
    coeffs_.readEntry("faceZone", faceZoneName_);

    Info<< decrIndent;

    return true;
}


void Foam::heatExchangerModels::referenceTemperature::write(const bool log)
{
    if (log)
    {
        Info<< nl
            << type() << ": " << name_ << nl << incrIndent
            << indent << "Net mass flux [kg/s]      : " << sumPhi_ << nl
            << indent << "Total heat exchange [W]   : " << Qt_ << nl
            << indent << "Reference T [K]           : " << Tref_
            << decrIndent;
    }

    if (Pstream::master())
    {
        Ostream& os = file();
        writeCurrentTime(os);

        os  << tab << sumPhi_
            << tab << Qt_
            << tab << Tref_
            << endl;
    }

    if (log) Info<< nl << endl;
}


// ************************************************************************* //
