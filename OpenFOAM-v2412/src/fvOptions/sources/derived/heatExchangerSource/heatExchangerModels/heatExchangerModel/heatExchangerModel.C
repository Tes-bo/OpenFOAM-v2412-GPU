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

#include "heatExchangerModel.H"
#include "coupledPolyPatch.H"
#include "surfaceInterpolate.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(heatExchangerModel, 0);
    defineRunTimeSelectionTable(heatExchangerModel, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatExchangerModel::heatExchangerModel
(
    const fvMesh& mesh,
    const word& name,
    const dictionary& coeffs
)
:
    writeFile(mesh, name, typeName, coeffs),
    mesh_(mesh),
    coeffs_(coeffs),
    name_(name),
    UName_("U"),
    TName_("T"),
    phiName_("phi"),
    faceZoneName_("unknown-faceZone"),
    faceId_(),
    facePatchId_(),
    faceSign_()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::heatExchangerModel::initialise()
{
    const label zoneID = mesh_.faceZones().findZoneID(faceZoneName_);

    if (zoneID < 0)
    {
        FatalErrorInFunction
            << type() << " " << name_ << ": "
            << "    Unknown face zone name: " << faceZoneName_
            << ". Valid face zones are: " << mesh_.faceZones().names()
            << exit(FatalError);
    }

    const faceZone& fZone = mesh_.faceZones()[zoneID];

    // Total number of faces selected
    label numFaces = fZone.size();

    faceId_.resize_nocopy(numFaces);
    facePatchId_.resize_nocopy(numFaces);
    faceSign_.resize_nocopy(numFaces);

    numFaces = 0;

    // TDB: handle multiple zones
    {
        forAll(fZone, i)
        {
            const label meshFacei = fZone[i];
            const label flipSign = (fZone.flipMap()[i] ? -1 : 1);

            // Internal faces
            label faceId = meshFacei;
            label facePatchId = -1;

            // Boundary faces
            if (!mesh_.isInternalFace(meshFacei))
            {
                facePatchId = mesh_.boundaryMesh().whichPatch(meshFacei);
                const polyPatch& pp = mesh_.boundaryMesh()[facePatchId];

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
                faceSign_[numFaces] = flipSign;

                ++numFaces;
            }
        }
    }

    // Shrink to size used
    faceId_.resize(numFaces);
    facePatchId_.resize(numFaces);
    faceSign_.resize(numFaces);
}


// ************************************************************************* //
