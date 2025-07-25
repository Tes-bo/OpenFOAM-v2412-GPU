/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 Louis Vittoz, SimScale GmbH
    Copyright (C) 2023 OpenCFD Ltd.
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

#include "fanMomentumSource.H"
#include "IFstream.H"
#include "addToRunTimeSelectionTable.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(fanMomentumSource, 0);
    addToRunTimeSelectionTable(option, fanMomentumSource, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions   * * * * * * * * * * * //

void Foam::fv::fanMomentumSource::writeProps
(
    const scalar gradP,
    const scalar flowRate
) const
{
    // Only write on output time
    if (mesh_.time().writeTime())
    {
        IOdictionary propsDict
        (
            IOobject
            (
                name_ + "Properties",
                mesh_.time().timeName(),
                "uniform",
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            )
        );
        propsDict.add("gradient", gradP);
        propsDict.add("flow_rate", flowRate);
        propsDict.regIOobject::write();
    }
}


void Foam::fv::fanMomentumSource::initUpstreamFaces()
{
    // First calculate the centre of gravity of the cell zone
    const vectorField& cellCentres = mesh_.cellCentres();
    const scalarField& cellVolumes = mesh_.cellVolumes();

    vector centreGravityCellZone(Zero);
    scalar cellZoneVolume = 0;
    for (const label celli : cells_)
    {
        const scalar cellVolume = cellVolumes[celli];
        centreGravityCellZone += cellCentres[celli]*cellVolume;
        cellZoneVolume += cellVolume;
    }

    reduce(centreGravityCellZone, sumOp<vector>());
    reduce(cellZoneVolume, sumOp<scalar>());

    if (cellZoneVolume < SMALL)
    {
        FatalErrorInFunction
            << "Cannot initialize upstream faces because the total cell zone"
            << " volume is zero or negative: " << cellZoneVolume << "." << nl
            << "Check whether there are cells in fan momentum source cell zone."
            << abort(FatalError);
    }

    centreGravityCellZone /= cellZoneVolume;

    // Collect faces upstream of the centre of gavity
    const faceZone& fZone = mesh_.faceZones()[surroundingFaceZoneID_];
    const vectorField& faceCentreGravity = mesh_.faceCentres();

    upstreamPatchFaceInfo_.resize_nocopy(fZone.size());

    label count = 0;
    for (const label facei : fZone)
    {
        if
        (
            (flowDir_ & (faceCentreGravity[facei] - centreGravityCellZone)) < 0
        )
        {
            labelPair patchFaceInfo(-1, -1);

            if (mesh_.isInternalFace(facei))
            {
                // Patch ID already set to -1, set only the face ID
                patchFaceInfo.second() = facei;
            }
            else
            {
                patchFaceInfo.first() = mesh_.boundaryMesh().whichPatch(facei);
                const polyPatch& pp = mesh_.boundaryMesh()[patchFaceInfo.first()];
                const auto* cpp = isA<coupledPolyPatch>(pp);

                if (cpp)
                {
                    patchFaceInfo.second() =
                        cpp->owner()
                      ? pp.whichFace(facei)
                      : -1;
                }
                else if (!isA<emptyPolyPatch>(pp))
                {
                    patchFaceInfo.second() = pp.whichFace(facei);
                }
                // else both face ID and patch ID remain at -1
            }

            // If this is an upstream face, set it in the list
            if (patchFaceInfo.second() >= 0)
            {
                upstreamPatchFaceInfo_[count] = patchFaceInfo;
                count++;
            }
        }
    }

    upstreamPatchFaceInfo_.setSize(count);

    // Fill cellsInZones_ with all cell IDs
    for (const label celli : cells_)
    {
        cellsInZones_.insert(celli);
    }

    // Sanity check
    const labelUList& owners = mesh_.owner();
    const labelUList& neighbours = mesh_.neighbour();
    for (const labelPair& patchFaceInfo : upstreamPatchFaceInfo_)
    {
        if (patchFaceInfo.first() == -1)
        {
            const label facei = patchFaceInfo.second();
            const label own = owners[facei];
            const label nei = neighbours[facei];

            // To be valid: one cell has to be inside the cellZone and the other
            // one, outside
            if (cellsInZones_.found(own) == cellsInZones_.found(nei))
            {
                FatalErrorInFunction
                    << "It seems that the faceZone is not part of the cellZone "
                    << "boundaries."
                    << abort(FatalError);
            }
        }
    }
}


Foam::scalar
Foam::fv::fanMomentumSource::calcFlowRate
(
    const surfaceScalarField& phi,
    const surfaceScalarField& rhof
) const
{
    // Sanity check to make sure we're not left in an inconsistent state because
    // here we're operating on primitive (non-dimensional) fields
    if (isNull(rhof) != (phi.dimensions() == dimVolume/dimTime))
    {
        FatalErrorInFunction
            << "Incorrect usage of the function. " << nl
            << "For incompressible flow, pass surfaceScalarField::null for rhof"
            << " and volumetric flux for phi." << nl
            << "For compressible flow, pass face-interpolated density as rhof"
            << " and mass flux for phi."
            << exit(FatalError);
    }

    // Calculate the flow rate through the upstream faces
    scalarList phif(upstreamPatchFaceInfo_.size());

    const labelUList& owners = mesh_.owner();

    forAll(upstreamPatchFaceInfo_, i)
    {
        const labelPair& patchFaceInfo = upstreamPatchFaceInfo_[i];

        if (isNull(rhof))
        {
            // Incompressible variant (phi = volumetric flux)
            phif[i] =
                patchFaceInfo.first() < 0
              ? phi[patchFaceInfo.second()]
              : phi.boundaryField()[patchFaceInfo];
        }
        else
        {
            // Compressible variant (phi = mass flux)
            phif[i] =
               patchFaceInfo.first() < 0
             ? phi[patchFaceInfo.second()]/
               rhof.internalField()[patchFaceInfo.second()]
             : phi.boundaryField()[patchFaceInfo]/
               rhof.boundaryField()[patchFaceInfo];
        }

        // Sign of the flux needs to be flipped if this is an internal face
        // whose owner is found in the cell zone
        if
        (
            patchFaceInfo.first() < 0
         && cellsInZones_.found(owners[patchFaceInfo.second()])
        )
        {
            phif[i] *= -1;
        }
    }

    return gSum(phif);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::fanMomentumSource::fanMomentumSource
(
    const word& sourceName,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fv::cellSetOption(sourceName, modelType, dict, mesh),
    upstreamPatchFaceInfo_(),
    cellsInZones_(),
    fanCurvePtr_(Function1<scalar>::New("fanCurve", coeffs_, &mesh)),
    UName_(coeffs_.getOrDefault<word>("U", "U")),
    flowDir_(coeffs_.get<vector>("flowDir")),
    thickness_(coeffs_.get<scalar>("thickness")),
    gradPFan_(0.0),
    surroundingFaceZoneID_(-1),
    rho_(coeffs_.getOrDefault<scalar>("rho", SMALL)),
    useRefRho_(coeffs_.found("rho"))
{
    // Skip all the checks if the source term has been deactivated
    // because there are no selected cells.
    // Such a situation typically occurs for multiple fluid regions
    if (fv::option::isActive())
    {
        const word faceZoneName = coeffs_.getWord("faceZone");

        surroundingFaceZoneID_ = mesh_.faceZones().findZoneID(faceZoneName);

        if (surroundingFaceZoneID_ < 0)
        {
            FatalErrorInFunction
                << type() << " " << this->name() << ": "
                << "    Unknown face zone name: " << faceZoneName
                << ". Valid face zones are: " << mesh_.faceZones().names()
                << exit(FatalError);
        }

        flowDir_.normalise();
        if (mag(flowDir_) < SMALL)
        {
            FatalIOErrorInFunction(coeffs_)
                << "'flowDir' cannot be a zero-valued vector."
                << exit(FatalIOError);
        }

        if (thickness_ < VSMALL)
        {
            FatalIOErrorInFunction(coeffs_)
                << "'thickness' cannot be non-positive."
                << exit(FatalIOError);
        }

        if (coeffs_.found("fields"))
        {
            IOWarningInFunction(coeffs_)
                << "Found 'fields' entry, which will be ignored. This fvOption"
                << " will operate only on field 'U' = " << UName_
                << endl;
        }
        fieldNames_.resize(1, UName_);

        fv::option::resetApplied();

        // Read the initial pressure gradient from file if it exists
        IFstream propsFile
        (
            mesh_.time().timePath()/"uniform"/(name_ + "Properties")
        );

        if (propsFile.good())
        {
            Info<< "    Reading pressure gradient from file" << endl;
            dictionary propsDict(propsFile);
            propsDict.readEntry("gradient", gradPFan_);
        }

        Info<< "    Initial pressure gradient = " << gradPFan_ << nl << endl;

        if (rho_ < SMALL)
        {
            FatalIOErrorInFunction(coeffs_)
                << "'rho' cannot be zero or negative."
                << exit(FatalIOError);
        }

        initUpstreamFaces();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::fanMomentumSource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    volVectorField::Internal Su
    (
        IOobject
        (
            name_ + fieldNames_[fieldi] + "Sup",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh_,
        dimensionedVector(eqn.dimensions()/dimVolume, Zero)
    );

    const auto& phi = mesh().lookupObject<surfaceScalarField>("phi");

    if (phi.dimensions() != dimVelocity*dimArea)
    {
        FatalErrorInFunction
            << "You called incompressible variant of addSup for case with "
            << "a mass flux and not volumetric flux. This is not allowed."
            << abort(FatalError);
    }

    if (!useRefRho_)
    {
        FatalErrorInFunction
            << "You called incompressible addSup without providing a "
            << "reference density. Add 'rho' entry to the dictionary."
            << abort(FatalError);
    }

    const scalar flowRate = calcFlowRate(phi, surfaceScalarField::null());

    // Pressure drop for this flow rate
    // if flow rate is negative, pressure is clipped at the static pressure
    gradPFan_ =
        fanCurvePtr_->value(max(flowRate, scalar(0)))/thickness_/rho_;

    // Create the source term
    UIndirectList<vector>(Su, cells_) = flowDir_*gradPFan_;

    eqn += Su;

    writeProps(gradPFan_, flowRate);
}


void Foam::fv::fanMomentumSource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    volVectorField::Internal Su
    (
        IOobject
        (
            name_ + fieldNames_[fieldi] + "Sup",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh_,
        dimensionedVector(eqn.dimensions()/dimVolume, Zero)
    );

    const auto& phi = mesh().lookupObject<surfaceScalarField>("phi");

    if (phi.dimensions() != dimMass/dimTime)
    {
        FatalErrorInFunction
            << "You called compressible variant of addSup for case with "
            << "a volumetric flux and not mass flux. This is not allowed."
            << abort(FatalError);
    }

    const surfaceScalarField rhof(fvc::interpolate(rho));
    const scalar flowRate = calcFlowRate(phi, rhof);

    // Pressure drop for this flow rate
    // if flow rate is negative, pressure is clipped at the static pressure
    gradPFan_ = fanCurvePtr_->value(max(flowRate, scalar(0)))/thickness_;

    // Create the source term
    UIndirectList<vector>(Su, cells_) = flowDir_*gradPFan_;

    eqn += Su;

    writeProps(gradPFan_, flowRate);
}


bool Foam::fv::fanMomentumSource::read(const dictionary& dict)
{
    NotImplemented;

    return false;
}


// ************************************************************************* //
