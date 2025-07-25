/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022-2023 OpenCFD Ltd.
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

#include "effectivenessTable.H"
#include "basicThermo.H"
#include "surfaceInterpolate.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace heatExchangerModels
{
    defineTypeNameAndDebug(effectivenessTable, 0);
    addToRunTimeSelectionTable
    (
        heatExchangerModel,
        effectivenessTable,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::heatExchangerModels::effectivenessTable::writeFileHeader
(
    Ostream& os
) const
{
    writeFile::writeHeader(os, "Heat exchanger source");
    writeFile::writeCommented(os, "Time");
    writeFile::writeTabbed(os, "Net mass flux [kg/s]");
    writeFile::writeTabbed(os, "Total heat exchange [W]");
    writeFile::writeTabbed(os, "Secondary inlet T [K]");
    writeFile::writeTabbed(os, "Reference T [K]");
    writeFile::writeTabbed(os, "Effectiveness [-]");

    if (secondaryCpPtr_)
    {
        writeFile::writeTabbed(os, "Secondary outlet T [K]");
    }

    os  << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatExchangerModels::effectivenessTable::effectivenessTable
(
    const fvMesh& mesh,
    const word& name,
    const dictionary& coeffs
)
:
    heatExchangerModel(mesh, name, coeffs),
    userPrimaryInletT_(false),
    targetQdotActive_(false),
    secondaryCpPtr_
    (
        Function1<scalar>::NewIfPresent
        (
            "secondaryCp",
            coeffs,
            word::null,
            &mesh
        )
    ),
    eTable_(),
    targetQdotCalcInterval_(5),
    secondaryMassFlowRate_(0),
    secondaryInletT_(0),
    primaryInletT_(0),
    targetQdot_(0),
    targetQdotRelax_(0.5),
    sumPhi_(0),
    Qt_(0),
    Tref_(0),
    effectiveness_(0)
{
    writeFileHeader(file());
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::heatExchangerModels::effectivenessTable::initialise()
{
    eTable_.reset(new interpolation2DTable<scalar>(coeffs_));

    heatExchangerModel::initialise();
}


Foam::tmp<Foam::scalarField>
Foam::heatExchangerModels::effectivenessTable::energyDensity
(
    const labelList& cells
)
{
    const auto& thermo = mesh_.lookupObject<basicThermo>(basicThermo::dictName);
    const auto& phi = mesh_.lookupObject<surfaceScalarField>(phiName_);
    const auto& T = mesh_.lookupObject<volScalarField>(TName_);

    const surfaceScalarField Cpf(fvc::interpolate(thermo.Cp()));
    const surfaceScalarField Tf(fvc::interpolate(T));

    sumPhi_ = 0;
    scalar sumMagPhi = 0;
    scalar CpfMean = 0;
    scalar primaryInletTfMean = 0;
    forAll(faceId_, i)
    {
        const label facei = faceId_[i];
        if (facePatchId_[i] != -1)
        {
            const label patchi = facePatchId_[i];
            const scalar phii = phi.boundaryField()[patchi][facei]*faceSign_[i];
            const scalar magPhii = mag(phii);

            sumPhi_ += phii;
            sumMagPhi += magPhii;

            const scalar Cpfi = Cpf.boundaryField()[patchi][facei];
            const scalar Tfi = Tf.boundaryField()[patchi][facei];

            CpfMean += Cpfi*magPhii;
            primaryInletTfMean += Tfi*magPhii;
        }
        else
        {
            const scalar phii = phi[facei]*faceSign_[i];
            const scalar magPhii = mag(phii);

            sumPhi_ += phii;
            sumMagPhi += magPhii;

            CpfMean += Cpf[facei]*magPhii;
            primaryInletTfMean += Tf[facei]*magPhii;
        }
    }
    reduce(CpfMean, sumOp<scalar>());
    reduce(sumPhi_, sumOp<scalar>());
    reduce(sumMagPhi, sumOp<scalar>());

    CpfMean /= sumMagPhi + ROOTVSMALL;

    scalar primaryInletT = primaryInletT_;
    if (!userPrimaryInletT_)
    {
        reduce(primaryInletTfMean, sumOp<scalar>());
        primaryInletT = primaryInletTfMean/(sumMagPhi + ROOTVSMALL);
    }

    effectiveness_ = eTable_()(mag(sumPhi_), secondaryMassFlowRate_);

    const scalar alpha = effectiveness_*CpfMean*mag(sumPhi_);

    Qt_ = alpha*(secondaryInletT_ - primaryInletT);

    if
    (
        targetQdotActive_
     && (mesh_.time().timeIndex() % targetQdotCalcInterval_ == 0)
    )
    {
        const scalar dT = (targetQdot_ - Qt_)/(alpha + ROOTVSMALL);
        secondaryInletT_ += targetQdotRelax_*dT;
    }

    const scalarField TCells(T, cells);
    scalarField deltaTCells(cells.size(), Zero);
    Tref_ = 0;
    if (Qt_ > 0)
    {
        Tref_ = gMax(TCells);
        forAll(deltaTCells, i)
        {
            deltaTCells[i] = max(Tref_ - TCells[i], scalar(0));
        }
    }
    else
    {
        Tref_ = gMin(TCells);
        forAll(deltaTCells, i)
        {
            deltaTCells[i] = max(TCells[i] - Tref_, scalar(0));
        }
    }

    const auto& U = mesh_.lookupObject<volVectorField>(UName_);
    const scalarField& V = mesh_.V();
    scalar sumWeight = 0;
    forAll(cells, i)
    {
        const label celli = cells[i];
        sumWeight += V[celli]*mag(U[celli])*deltaTCells[i];
    }
    reduce(sumWeight, sumOp<scalar>());

    return Qt_*deltaTCells/(sumWeight + ROOTVSMALL);
}


bool Foam::heatExchangerModels::effectivenessTable::read(const dictionary& dict)
{
    if (!writeFile::read(dict))
    {
        return false;
    }

    Info<< incrIndent << indent << "- using model: " << type() << endl;

    coeffs_.readEntry("secondaryMassFlowRate", secondaryMassFlowRate_);
    coeffs_.readEntry("secondaryInletT", secondaryInletT_);

    if (coeffs_.readIfPresent("primaryInletT", primaryInletT_))
    {
        userPrimaryInletT_ = true;

        Info<< indent
            << "- using user-specified primary flow inlet temperature: "
            << primaryInletT_ << endl;
    }
    else
    {
        Info<< indent
            << "- using flux-weighted primary flow inlet temperature"
            << endl;
    }

    if (coeffs_.readIfPresent("targetQdot", targetQdot_))
    {
        targetQdotActive_ = true;

        Info<< indent << "- using target heat rejection: " << targetQdot_ << nl;

        coeffs_.readIfPresent
        (
            "targetQdotCalcInterval",
            targetQdotCalcInterval_
        );

        Info<< indent << "- updating secondary inlet temperature every "
            << targetQdotCalcInterval_ << " iterations" << nl;

        coeffs_.readIfPresent("targetQdotRelax", targetQdotRelax_);

        Info<< indent << "- temperature relaxation: "
            << targetQdotRelax_ << endl;
    }

    UName_ = coeffs_.getOrDefault<word>("U", "U");
    TName_ = coeffs_.getOrDefault<word>("T", "T");
    phiName_ = coeffs_.getOrDefault<word>("phi", "phi");
    coeffs_.readEntry("faceZone", faceZoneName_);

    Info<< decrIndent;

    return true;
}


void Foam::heatExchangerModels::effectivenessTable::write
(
    const bool log
)
{
    if (log)
    {
        Info<< nl
            << type() << ": " << name_ << nl << incrIndent
            << indent << "Net mass flux [kg/s]      : " << sumPhi_ << nl
            << indent << "Total heat exchange [W]   : " << Qt_ << nl
            << indent << "Secondary inlet T [K]     : " << secondaryInletT_<< nl
            << indent << "Reference T [K]           : " << Tref_ << nl
            << indent << "Effectiveness [-]         : " << effectiveness_
            << decrIndent;
    }

    if (Pstream::master())
    {
        Ostream& os = file();
        writeCurrentTime(os);

        os  << tab << sumPhi_
            << tab << Qt_
            << tab << secondaryInletT_
            << tab << Tref_
            << tab << effectiveness_;

        if (secondaryCpPtr_)
        {
            // Secondary Cp as a function of the starting secondary temperature
            const scalar secondaryCp = secondaryCpPtr_->value(secondaryInletT_);
            const scalar secondaryOutletT =
                secondaryInletT_ - Qt_/(secondaryMassFlowRate_*secondaryCp);

            if (log)
            {
                Info << nl << incrIndent << indent
                    << "Secondary outlet T [K]    : " << secondaryOutletT
                    << decrIndent;
            }

            os  << tab << secondaryOutletT;
        }
        os  << endl;
    }

    Info<< nl << endl;
}


// ************************************************************************* //
