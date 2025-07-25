/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2022 OpenCFD Ltd.
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

#include "contactHeatFluxSource.H"
#include "faMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "famSup.H"
#include "zeroGradientFaPatchFields.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fa
{
    defineTypeNameAndDebug(contactHeatFluxSource, 0);
    addToRunTimeSelectionTable(option, contactHeatFluxSource, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fa::contactHeatFluxSource::contactHeatFluxSource
(
    const word& sourceName,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fa::faceSetOption(sourceName, modelType, dict, mesh),
    TName_(dict.getOrDefault<word>("T", "T")),
    TprimaryName_(dict.get<word>("Tprimary")),
    Tprimary_(mesh_.lookupObject<volScalarField>(TprimaryName_)),
    thicknessLayers_(),
    kappaLayers_(),
    contactRes_(0),
    curTimeIndex_(-1),
    coupling_()
{
    fieldNames_.resize(1, TName_);

    fa::option::resetApplied();

    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::areaMesh>>
Foam::fa::contactHeatFluxSource::htc() const
{
    auto thtc = DimensionedField<scalar, areaMesh>::New
    (
        "htc_" + option::name(),
        regionMesh(),
        dimensionedScalar(dimPower/dimArea/dimTemperature, Zero)
    );
    auto& htc = thtc.ref();

    // Set up values per coupled patch

    PtrList<scalarField> patchValues(coupling_.size());

    forAll(coupling_, patchi)
    {
        const auto* tempCoupled = coupling_.get(patchi);

        if (tempCoupled)
        {
            const fvPatch& p = tempCoupled->patch();

            patchValues.set
            (
                patchi,
                (
                    tempCoupled->kappa(p.patchInternalField(Tprimary_))
                  * p.deltaCoeffs()
                )
            );
        }
    }

    vsm().mapToSurface<scalar>(patchValues, htc.field());

    if (contactRes_ != 0)
    {
        htc.field() += contactRes_;
    }

    // Zero htc for non-mapped faces
    faceSetOption::subsetFilter(htc.field());

    return thtc;
}


void Foam::fa::contactHeatFluxSource::addSup
(
    const areaScalarField& h,
    const areaScalarField& rhoCph,
    faMatrix<scalar>& eqn,
    const label fieldi
)
{
    if (isActive())
    {
        DebugInfo
            << name() << ": applying source to "
            << eqn.psi().name() << endl;

        if (curTimeIndex_ != mesh().time().timeIndex())
        {
            tmp<DimensionedField<scalar, areaMesh>> htcw(htc());

            // Wall temperature - mapped from primary field to finite-area
            auto Twall = DimensionedField<scalar, areaMesh>::New
            (
                "Tw_" + option::name(),
                regionMesh(),
                dimensionedScalar(dimTemperature, Zero)
            );

            vsm().mapInternalToSurface<scalar>(Tprimary_, Twall.ref().field());

            eqn += -fam::Sp(htcw(), eqn.psi()) + htcw()*Twall;

            curTimeIndex_ = mesh().time().timeIndex();
        }
    }
}


bool Foam::fa::contactHeatFluxSource::read(const dictionary& dict)
{
    if (fa::option::read(dict))
    {
        coeffs_.readIfPresent("T", TName_);

        contactRes_ = 0;

        if (dict.readIfPresent("thicknessLayers", thicknessLayers_))
        {
            dict.readEntry("kappaLayers", kappaLayers_);

            // Calculate effective thermal resistance by harmonic averaging

            forAll(thicknessLayers_, iLayer)
            {
                contactRes_ += thicknessLayers_[iLayer]/kappaLayers_[iLayer];
            }

            if (thicknessLayers_.size())
            {
                contactRes_ = scalar(1)/contactRes_;
            }
        }


        // Set up coupling (per-patch) for referenced polyPatches (sorted order)
        // - size is maxPolyPatch+1

        const labelList& patches = regionMesh().whichPolyPatches();

        coupling_.clear();
        coupling_.resize(patches.empty() ? 0 : (patches.last()+1));

        for (const label patchi : patches)
        {
            const fvPatch& p = mesh_.boundary()[patchi];

            coupling_.set
            (
                patchi,
                new temperatureCoupling(p, dict)
            );
        }

        return true;
    }

    return false;
}


// ************************************************************************* //
