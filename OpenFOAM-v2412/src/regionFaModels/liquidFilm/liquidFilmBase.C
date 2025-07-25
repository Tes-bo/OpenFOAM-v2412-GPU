/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2023 OpenCFD Ltd.
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

#include "liquidFilmBase.H"
#include "gravityMeshObject.H"
#include "movingWallVelocityFvPatchVectorField.H"
#include "turbulentFluidThermoModel.H"
#include "turbulentTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(liquidFilmBase, 0);

defineRunTimeSelectionTable(liquidFilmBase, dictionary);

const Foam::word liquidFilmName("liquidFilm");

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

liquidFilmBase::liquidFilmBase
(
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    regionFaModel(mesh, liquidFilmName, modelType, dict, true),
    momentumPredictor_
    (
        this->solution().subDict("PIMPLE").get<bool>("momentumPredictor")
    ),
    nOuterCorr_
    (
        this->solution().subDict("PIMPLE").get<label>("nOuterCorr")
    ),
    nCorr_(this->solution().subDict("PIMPLE").get<label>("nCorr")),
    nFilmCorr_
    (
        this->solution().subDict("PIMPLE").get<label>("nFilmCorr")
    ),
    h0_("h0", dimLength, 1e-7, dict),
    deltaWet_("deltaWet", dimLength, 1e-4, dict),
    UName_(dict.get<word>("U")),
    pName_(dict.getOrDefault<word>("p",  word::null)),
    pRef_(dict.get<scalar>("pRef")),
    h_
    (
        IOobject
        (
            "hf_" + regionName_,
            regionMesh().time().timeName(),
            regionMesh().thisDb(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh()
    ),
    Uf_
    (
        IOobject
        (
            "Uf_" + regionName_,
            regionMesh().time().timeName(),
            regionMesh().thisDb(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh()
    ),
    pf_
    (
        IOobject
        (
            "pf_" + regionName_,
            regionMesh().time().timeName(),
            regionMesh().thisDb(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimPressure, Zero)
    ),
    ppf_
    (
        IOobject
        (
            "ppf_" + regionName_,
            regionMesh().time().timeName(),
            regionMesh().thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimPressure, Zero)
    ),
    phif_
    (
        IOobject
        (
            "phif_" + regionName_,
            regionMesh().time().timeName(),
            regionMesh().thisDb(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fac::interpolate(Uf_) & regionMesh().Le()
    ),
    phi2s_
    (
        IOobject
        (
            "phi2s_" + regionName_,
            regionMesh().time().timeName(),
            regionMesh().thisDb(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fac::interpolate(h_*Uf_) & regionMesh().Le()
    ),
    gn_
    (
        IOobject
        (
            "gn",
            regionMesh().time().timeName(),
            regionMesh().thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimAcceleration, Zero)
    ),
    g_(meshObjects::gravity::New(primaryMesh().time())),
    massSource_
    (
        IOobject
        (
            "massSource",
            primaryMesh().time().timeName(),
            primaryMesh().thisDb()
        ),
        primaryMesh(),
        dimensionedScalar(dimMass, Zero)
    ),
    momentumSource_
    (
        IOobject
        (
            "momentumSource",
            primaryMesh().time().timeName(),
            primaryMesh().thisDb()
        ),
        primaryMesh(),
        dimensionedVector(dimPressure, Zero)
    ),
    pnSource_
    (
        IOobject
        (
            "pnSource",
            primaryMesh().time().timeName(),
            primaryMesh().thisDb()
        ),
        primaryMesh(),
        dimensionedScalar(dimPressure, Zero)
    ),
    addedMassTotal_(0),
    faOptions_(Foam::fa::options::New(primaryMesh()))
{
    const areaVectorField& ns = regionMesh().faceAreaNormals();

    gn_ = g_ & ns;

    if (!faOptions_.optionList::size())
    {
        Info << "No finite area options present" << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

liquidFilmBase::~liquidFilmBase()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

scalar liquidFilmBase::CourantNumber() const
{
    scalar CoNum = 0.0;
    scalar velMag = 0.0;

    edgeScalarField SfUfbyDelta
    (
        regionMesh().edgeInterpolation::deltaCoeffs()*mag(phif_)
    );

    CoNum =
        max(SfUfbyDelta/regionMesh().magLe()).value()*time().deltaTValue();

    velMag = max(mag(phif_)/regionMesh().magLe()).value();

    reduce(CoNum, maxOp<scalar>());
    reduce(velMag, maxOp<scalar>());

    Info<< "Max film Courant Number: " << CoNum
        << " Film velocity magnitude: " << velMag << endl;

    return CoNum;
}


Foam::tmp<Foam::areaVectorField> liquidFilmBase::Uw() const
{
    auto tUw = areaVectorField::New
    (
        "tUw",
        IOobjectOption::NO_REGISTER,
        regionMesh(),
        dimensionedVector(dimVelocity, Zero)
    );
    auto& Uw = tUw.ref();

    if (primaryMesh().moving())
    {
        const labelList& patches = regionMesh().whichPolyPatches();

        // Set up mapping values per patch

        PtrMap<vectorField> patchValues(2*patches.size());

        for (const label patchi : patches)
        {
            const auto* wpp = isA<movingWallVelocityFvPatchVectorField>
            (
                primaryMesh().boundaryMesh()[patchi]
            );

            if (wpp)
            {
                patchValues.set(patchi, wpp->Uwall());
            }
        }

        if (patchValues.size())
        {
            // Map Up to area
            tmp<vectorField> UsWall = vsmPtr_->mapToSurface(patchValues);

            const vectorField& nHat =
                regionMesh().faceAreaNormals().internalField();

            Uw.primitiveFieldRef() = UsWall() - nHat*(UsWall() & nHat);
        }
    }

    return tUw;
}


Foam::tmp<Foam::areaVectorField> liquidFilmBase::Us() const
{
    auto tUs = areaVectorField::New
    (
        "tUs",
        IOobjectOption::NO_REGISTER,
        regionMesh(),
        dimensionedVector(dimVelocity, Zero)
    );

    // Uf quadratic profile
    tUs.ref() = Foam::sqrt(2.0)*Uf_;

    return tUs;
}


Foam::tmp<Foam::areaVectorField> liquidFilmBase::Up() const
{
    const volVectorField& Uprimary =
        primaryMesh().lookupObject<volVectorField>(UName_);

    auto tUp = areaVectorField::New
    (
        "tUp",
        IOobjectOption::NO_REGISTER,
        regionMesh(),
        dimensionedVector(dimVelocity, Zero)
    );
    auto& Up = tUp.ref();


    // Set up mapping values per patch

    const labelList& patches = regionMesh().whichPolyPatches();

    PtrMap<vectorField> patchValues(2*patches.size());

    // U tangential on primary
    for (const label patchi : patches)
    {
        const fvPatchVectorField& Uw = Uprimary.boundaryField()[patchi];

        patchValues.set(patchi, -Uw.snGrad());
    }


    // Map U tang to surface
    vsmPtr_->mapToSurface(patchValues, Up.primitiveFieldRef());

    Up.primitiveFieldRef() *= h_.primitiveField();

    // U tangent on surface
    const vectorField& nHat = regionMesh().faceAreaNormals().internalField();

    Up.primitiveFieldRef() -= nHat*(Up.primitiveField() & nHat);

    return tUp;
}


tmp<areaScalarField> liquidFilmBase::pg() const
{
    auto tpg = areaScalarField::New
    (
        "tpg",
        IOobjectOption::NO_REGISTER,
        regionMesh(),
        dimensionedScalar(dimPressure, Zero)
    );
    auto& pfg = tpg.ref();

    if (!pName_.empty())
    {
        vsmPtr_->mapInternalToSurface
        (
            primaryMesh().lookupObject<volScalarField>(pName_),
            pfg.primitiveFieldRef()
        );
    }

    return tpg;
}


tmp<areaScalarField> liquidFilmBase::alpha() const
{
    auto talpha = areaScalarField::New
    (
        "talpha",
        IOobjectOption::NO_REGISTER,
        regionMesh(),
        dimensionedScalar(dimless, Zero)
    );
    auto& alpha = talpha.ref();

    alpha = pos0(h_ - deltaWet_);

    return talpha;
}


void liquidFilmBase::addSources
(
    const label patchi,
    const label facei,
    const scalar massSource,
    const vector& momentumSource,
    const scalar pressureSource,
    const scalar energySource
)
{
    massSource_.boundaryFieldRef()[patchi][facei] += massSource;
    pnSource_.boundaryFieldRef()[patchi][facei] += pressureSource;
    momentumSource_.boundaryFieldRef()[patchi][facei] += momentumSource;
}


void liquidFilmBase::preEvolveRegion()
{
    regionFaModel::preEvolveRegion();
}


void liquidFilmBase::postEvolveRegion()
{
    if (debug && primaryMesh().time().writeTime())
    {
        massSource_.write();
        pnSource_.write();
        momentumSource_.write();
    }

    massSource_.boundaryFieldRef() = Zero;
    pnSource_.boundaryFieldRef() = Zero;
    momentumSource_.boundaryFieldRef() = Zero;

    regionFaModel::postEvolveRegion();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
