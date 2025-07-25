/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2018 OpenFOAM Foundation
    Copyright (C) 2016-2022 OpenCFD Ltd.
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

#include "greyDiffusiveRadiationMixedFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "boundaryRadiationProperties.H"

#include "fvDOM.H"
#include "constants.H"
#include "unitConversion.H"

using namespace Foam::constant;
using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::greyDiffusiveRadiationMixedFvPatchScalarField::
greyDiffusiveRadiationMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    TName_("T"),
    qRadExt_(0),
    qRadExtDir_(Zero)
{
    refValue() = Zero;
    refGrad() = Zero;
    valueFraction() = 1.0;
}


Foam::radiation::greyDiffusiveRadiationMixedFvPatchScalarField::
greyDiffusiveRadiationMixedFvPatchScalarField
(
    const greyDiffusiveRadiationMixedFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    TName_(ptf.TName_),
    qRadExt_(ptf.qRadExt_),
    qRadExtDir_(ptf.qRadExtDir_)
{}


Foam::radiation::greyDiffusiveRadiationMixedFvPatchScalarField::
greyDiffusiveRadiationMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    TName_(dict.getOrDefault<word>("T", "T")),
    qRadExt_(dict.getOrDefault<scalar>("qRadExt", 0)),
    qRadExtDir_(dict.getOrDefault<vector>("qRadExtDir", Zero))
{
    if (this->readMixedEntries(dict))
    {
        // Full restart
        this->readValueEntry(dict, IOobjectOption::MUST_READ);
    }
    else
    {
        refValue() = Zero;
        refGrad() = Zero;
        valueFraction() = 1.0;

        fvPatchScalarField::operator=(refValue());
    }
}


Foam::radiation::greyDiffusiveRadiationMixedFvPatchScalarField::
greyDiffusiveRadiationMixedFvPatchScalarField
(
    const greyDiffusiveRadiationMixedFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf),
    TName_(ptf.TName_),
    qRadExt_(ptf.qRadExt_),
    qRadExtDir_(ptf.qRadExtDir_)
{}


Foam::radiation::greyDiffusiveRadiationMixedFvPatchScalarField::
greyDiffusiveRadiationMixedFvPatchScalarField
(
    const greyDiffusiveRadiationMixedFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    TName_(ptf.TName_),
    qRadExt_(ptf.qRadExt_),
    qRadExtDir_(ptf.qRadExtDir_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::radiation::greyDiffusiveRadiationMixedFvPatchScalarField::
updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    const int oldTag = UPstream::incrMsgType();

    const auto& Tp = patch().lookupPatchField<volScalarField>(TName_);

    const fvDOM& dom = db().lookupObject<fvDOM>("radiationProperties");

    label rayId = -1;
    label lambdaId = -1;
    dom.setRayIdLambdaId(internalField().name(), rayId, lambdaId);

    const label patchi = patch().index();

    if (dom.nLambda() != 1)
    {
        FatalErrorInFunction
            << " a grey boundary condition is used with a non-grey "
            << "absorption model" << nl << exit(FatalError);
    }

    scalarField& Iw = *this;

    const vectorField n(patch().nf());

    radiativeIntensityRay& ray =
        const_cast<radiativeIntensityRay&>(dom.IRay(rayId));

    const scalarField nAve(n & ray.dAve());

    ray.qr().boundaryFieldRef()[patchi] += Iw*nAve;

    const boundaryRadiationProperties& boundaryRadiation =
        boundaryRadiationProperties::New(internalField().mesh());

    const tmp<scalarField> temissivity
    (
        boundaryRadiation.emissivity(patch().index(), 0, nullptr, &Tp)
    );

    const scalarField& emissivity = temissivity();

    const tmp<scalarField> ttransmissivity
    (
        boundaryRadiation.transmissivity(patch().index(), 0, nullptr, &Tp)
    );

    const scalarField& transmissivity = ttransmissivity();

    scalarField& qem = ray.qem().boundaryFieldRef()[patchi];
    scalarField& qin = ray.qin().boundaryFieldRef()[patchi];

    const vector& myRayId = dom.IRay(rayId).d();

    scalarField Ir(patch().size(), Zero);
    forAll(Iw, facei)
    {
        for (label rayi=0; rayi < dom.nRay(); rayi++)
        {
            const vector& d = dom.IRay(rayi).d();

            if ((-n[facei] & d) < 0.0)
            {
                // q into the wall
                const scalarField& IFace =
                    dom.IRay(rayi).ILambda(lambdaId).boundaryField()[patchi];

                const vector& rayDave = dom.IRay(rayi).dAve();
                Ir[facei] += IFace[facei]*(n[facei] & rayDave);
            }
        }
    }

    if (dom.useSolarLoad())
    {
        // Looking for primary heat flux single band
        Ir += patch().lookupPatchField<volScalarField>
        (
            dom.primaryFluxName_ + "_0"
        );

        word qSecName = dom.relfectedFluxName_ + "_0";

        if (this->db().foundObject<volScalarField>(qSecName))
        {
             const volScalarField& qSec =
                this->db().lookupObject<volScalarField>(qSecName);

            Ir += qSec.boundaryField()[patch().index()];
        }
    }

    scalarField Iexternal(this->size(), 0.0);

    if (dom.useExternalBeam())
    {
        const vector sunDir = dom.solarCalc().direction();
        const scalar directSolarRad = dom.solarCalc().directSolarRad();

        //label nRaysBeam = dom.nRaysBeam();
        label SunRayId(-1);
        scalar maxSunRay = -GREAT;

        // Looking for the ray closest to the Sun direction
        for (label rayI=0; rayI < dom.nRay(); rayI++)
        {
            const vector& iD = dom.IRay(rayI).d();
            scalar dir = sunDir & iD;
            if (dir > maxSunRay)
            {
                maxSunRay = dir;
                SunRayId = rayI;
            }
        }

        if (rayId == SunRayId)
        {
            const scalarField nAve(n & dom.IRay(rayId).dAve());
            forAll(Iexternal, faceI)
            {
                Iexternal[faceI] = directSolarRad/mag(dom.IRay(rayId).dAve());
            }
        }
    }

    scalarField Isource(this->size(), 0.0);

    if (qRadExt_ > 0)
    {
        if (mag(qRadExtDir_) > 0)
        {
            label rayqoId = -1;
            scalar maxRay = -GREAT;

            // Looking for the ray closest to the Sun direction
            for (label rayI = 0; rayI < dom.nRay(); ++rayI)
            {
                const vector& iD = dom.IRay(rayI).d();
                const scalar dir = qRadExtDir_ & iD;

                if (dir > maxRay)
                {
                    maxRay = dir;
                    rayqoId = rayI;
                }
            }

            if (rayId == rayqoId)
            {
                forAll(Isource, faceI)
                {
                    Isource[faceI] += qRadExt_/mag(dom.IRay(rayId).dAve());
                }
            }
        }
        else
        {
            forAll(Iw, faceI)
            {
                label rayqoId = -1;
                scalar maxRay = -GREAT;

                // Looking for the ray closest to the Sun direction
                for (label rayI = 0; rayI < dom.nRay(); ++rayI)
                {
                    const vector& iD = dom.IRay(rayI).d();
                    const scalar dir = -n[faceI] & iD;

                    if (dir > maxRay)
                    {
                        maxRay = dir;
                        rayqoId = rayI;
                    }
                }

                if (rayId == rayqoId)
                {
                    Isource[faceI] += qRadExt_/mag(dom.IRay(rayId).dAve());
                }
            }
        }
    }

    forAll(Iw, faceI)
    {
        if ((-n[faceI] & myRayId) > 0.0)
        {
            // direction out of the wall
            refGrad()[faceI] = 0.0;
            valueFraction()[faceI] = 1.0;
            refValue()[faceI] =
                Isource[faceI]
              + Iexternal[faceI]*transmissivity[faceI]
              + (
                    Ir[faceI]*(scalar(1) - emissivity[faceI])
                  + emissivity[faceI]*physicoChemical::sigma.value()
                  * pow4(Tp[faceI])
                )/pi;

            // Emitted heat flux from this ray direction
            qem[faceI] = refValue()[faceI]*nAve[faceI];
        }
        else
        {
            // direction into the wall
            valueFraction()[faceI] = 0.0;
            refGrad()[faceI] = 0.0;
            refValue()[faceI] = 0.0; //not used

            // Incident heat flux on this ray direction
            qin[faceI] = Iw[faceI]*nAve[faceI];
        }
    }

    UPstream::msgType(oldTag);  // Restore tag

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::radiation::greyDiffusiveRadiationMixedFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchField<scalar>::write(os);
    os.writeEntryIfDifferent<word>("T", "T", TName_);
    os.writeEntryIfDifferent<scalar>("qRadExt", Zero, qRadExt_);
    os.writeEntryIfDifferent<vector>("qRadExtDir", Zero, qRadExtDir_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{
    makePatchTypeField
    (
        fvPatchScalarField,
        greyDiffusiveRadiationMixedFvPatchScalarField
    );
}
}


// ************************************************************************* //
