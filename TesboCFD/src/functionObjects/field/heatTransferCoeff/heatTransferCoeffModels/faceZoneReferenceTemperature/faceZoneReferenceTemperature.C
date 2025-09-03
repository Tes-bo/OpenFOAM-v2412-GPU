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

#include "faceZoneReferenceTemperature.H"
#include "surfaceInterpolate.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace heatTransferCoeffModels
{
    defineTypeNameAndDebug(faceZoneReferenceTemperature, 0);
    addToRunTimeSelectionTable
    (
        heatTransferCoeffModel,
        faceZoneReferenceTemperature,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::heatTransferCoeffModels::faceZoneReferenceTemperature::
setFaceZoneFaces(const dictionary& dict)
{
    const auto& mesh =
        mesh_.objectRegistry::db().lookupObject<fvMesh>(refRegionName_);

    const word faceZoneName(dict.get<word>("referenceFaceZone"));

    faceZonei_ = mesh.faceZones().findZoneID(faceZoneName);

    if (faceZonei_ < 0)
    {
        FatalIOErrorInFunction(dict)
            << "referenceFaceZone: " << faceZoneName
            << " does not exist in referenceRegion: " << refRegionName_
            << exit(FatalIOError);
    }

    const faceZone& fZone = mesh.faceZones()[faceZonei_];

    label numFaces = fZone.size();

    if (!returnReduceOr(numFaces))
    {
        FatalIOErrorInFunction(dict)
            << "referenceFaceZone: " << faceZoneName
            << " contains no faces."
            << exit(FatalIOError);
    }

    faceId_.resize_nocopy(numFaces);
    facePatchId_.resize_nocopy(numFaces);

    numFaces = 0;

    // TDB: handle multiple zones
    {
        forAll(fZone, i)
        {
            const label meshFacei = fZone[i];

            // Internal faces
            label faceId = meshFacei;
            label facePatchId = -1;

            // Boundary faces
            if (!mesh.isInternalFace(meshFacei))
            {
                facePatchId = mesh.boundaryMesh().whichPatch(meshFacei);
                const polyPatch& pp = mesh.boundaryMesh()[facePatchId];

                if (isA<emptyPolyPatch>(pp))
                {
                    continue;  // Ignore empty patch
                }

                const auto* cpp = isA<coupledPolyPatch>(pp);

                if (cpp && !cpp->owner())
                {
                    continue;  // Ignore neighbour side
                }

                faceId = pp.whichFace(meshFacei);
            }

            if (faceId >= 0)
            {
                faceId_[numFaces] = faceId;
                facePatchId_[numFaces] = facePatchId;

                ++numFaces;
            }
        }
    }

    // Shrink to size used
    faceId_.resize(numFaces);
    facePatchId_.resize(numFaces);
}


Foam::scalar Foam::heatTransferCoeffModels::faceZoneReferenceTemperature::
faceZoneAverageTemperature()
{
    const auto& mesh =
        mesh_.objectRegistry::db().lookupObject<fvMesh>(refRegionName_);

    const auto& T = mesh.lookupObject<volScalarField>(TName_);
    const surfaceScalarField Tf(fvc::interpolate(T));

    const surfaceScalarField& magSf = mesh.magSf();

    scalar Tmean = 0;
    scalar sumMagSf = 0;

    forAll(faceId_, i)
    {
        const label facei = faceId_[i];
        if (facePatchId_[i] != -1)
        {
            const label patchi = facePatchId_[i];
            const scalar sf = magSf.boundaryField()[patchi][facei];

            Tmean += Tf.boundaryField()[patchi][facei]*sf;
            sumMagSf += sf;
        }
        else
        {
            const scalar sf = magSf[facei];
            Tmean += Tf[facei]*sf;
            sumMagSf += sf;
        }
    }
    reduce(Tmean, sumOp<scalar>());
    reduce(sumMagSf, sumOp<scalar>());

    Tmean /= sumMagSf;

    return Tmean;
}


void Foam::heatTransferCoeffModels::faceZoneReferenceTemperature::htc
(
    volScalarField& htc,
    const FieldField<Field, scalar>& q
)
{
    // Retrieve temperature boundary fields for current region
    const auto& T = mesh_.lookupObject<volScalarField>(TName_);
    const volScalarField::Boundary& Tbf = T.boundaryField();

    // Retrieve heat-transfer coefficient boundary fields for current region
    volScalarField::Boundary& htcBf = htc.boundaryFieldRef();

    // Calculate area-averaged temperature field
    // for the reference face zone and region
    // (reference region can be different from current region)
    const scalar Tref = faceZoneAverageTemperature();

    // Calculate heat-transfer coefficient boundary fields for current region
    for (const label patchi : patchIDs_)
    {
        htcBf[patchi] = q[patchi]/(Tref - Tbf[patchi] + ROOTVSMALL);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatTransferCoeffModels::faceZoneReferenceTemperature::
faceZoneReferenceTemperature
(
    const dictionary& dict,
    const fvMesh& mesh,
    const word& TName
)
:
    heatTransferCoeffModel(dict, mesh, TName),
    faceZonei_(-1),
    refRegionName_(polyMesh::defaultRegion),
    faceId_(),
    facePatchId_()
{
    read(dict);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::heatTransferCoeffModels::faceZoneReferenceTemperature::read
(
    const dictionary& dict
)
{
    if (!heatTransferCoeffModel::read(dict))
    {
        return false;
    }

    dict.readIfPresent("referenceRegion", refRegionName_);

    setFaceZoneFaces(dict);

    return true;
}


// ************************************************************************* //
