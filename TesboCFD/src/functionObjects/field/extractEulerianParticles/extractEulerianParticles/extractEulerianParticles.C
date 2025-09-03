/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2023 OpenCFD Ltd.
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

#include "extractEulerianParticles.H"
#include "regionSplit2D.H"
#include "mathematicalConstants.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "surfaceInterpolate.H"
#include "pairPatchAgglomeration.H"
#include "emptyPolyPatch.H"
#include "coupledPolyPatch.H"
#include "binned.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(extractEulerianParticles, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        extractEulerianParticles,
        dictionary
    );
}
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::extractEulerianParticles::checkFaceZone()
{
    DebugInFunction << endl;

    zoneID_ = mesh_.faceZones().findZoneID(faceZoneName_);
    if (zoneID_ == -1)
    {
        FatalErrorInFunction
            << "Unable to find faceZone " << faceZoneName_
            << ".  Available faceZones are: " << mesh_.faceZones().names()
            << exit(FatalError);
    }

    const faceZone& fz = mesh_.faceZones()[zoneID_];
    const label nFaces = fz.size();
    const label allFaces = returnReduce(nFaces, sumOp<label>());

    if (allFaces < nInjectorLocations_)
    {
        FatalErrorInFunction
            << "faceZone " << faceZoneName_
            << ": Number of faceZone faces (" << allFaces
            << ") is less than the number of requested locations ("
            << nInjectorLocations_ << ")."
            << exit(FatalError);
    }

    Info<< type() << " " << name() << " output:" << nl
        << "    faceZone : " << faceZoneName_ << nl
        << "    faces    : " << allFaces << nl
        << endl;

    // Initialise old iteration blocked faces
    // Note: for restart, this info needs to be written/read
    regions0_.setSize(fz.size(), -1);
}


void Foam::functionObjects::extractEulerianParticles::initialiseBins()
{
    DebugInFunction << endl;

    if (!nInjectorLocations_)
    {
        return;
    }

    const faceZone& fz = mesh_.faceZones()[zoneID_];

    // Agglomerate faceZone faces into nInjectorLocations_ global locations
    const indirectPrimitivePatch patch
    (
        IndirectList<face>(mesh_.faces(), fz),
        mesh_.points()
    );

    const label nFaces = fz.size();
    label nLocations = nInjectorLocations_;

    if (Pstream::parRun())
    {
        label nGlobalFaces = returnReduce(nFaces, sumOp<label>());
        scalar fraction = scalar(nFaces)/scalar(nGlobalFaces);
        nLocations = ceil(fraction*nInjectorLocations_);
        if (debug)
        {
            Pout<< "nFaces:" << nFaces
                << ", nGlobalFaces:" << nGlobalFaces
                << ", fraction:" << fraction
                << ", nLocations:" << nLocations
                << endl;
        }
    }

    pairPatchAgglomeration ppa
    (
        patch.localFaces(),
        patch.localPoints(),
        10,
        50,
        nLocations,
        labelMax,
        180
    );

    ppa.agglomerate();

    label nCoarseFaces = 0;
    if (nFaces != 0)
    {
        fineToCoarseAddr_ = ppa.restrictTopBottomAddressing();
        nCoarseFaces = max(fineToCoarseAddr_) + 1;
    }

    globalCoarseFaces_ = globalIndex(nCoarseFaces);

    Info<< "Created " << returnReduce(nCoarseFaces, sumOp<label>())
        << " coarse faces" << endl;
}


Foam::tmp<Foam::surfaceScalarField>
Foam::functionObjects::extractEulerianParticles::phiU() const
{
    DebugInFunction << endl;

    const auto& phi = mesh_.lookupObject<surfaceScalarField>(phiName_);

    if (phi.dimensions() == dimMass/dimTime)
    {
        const auto& rho = mesh_.lookupObject<volScalarField>(rhoName_);

        return phi/fvc::interpolate(rho);
    }

    return phi;
}


void Foam::functionObjects::extractEulerianParticles::setBlockedFaces
(
    const surfaceScalarField& alphaf,
    const faceZone& fz,
    boolList& blockedFaces
)
{
    DebugInFunction << endl;

    // Initialise storage for patch and patch-face indices where faceZone
    // intersects mesh patch(es)
    patchIDs_.setSize(fz.size(), -1);
    patchFaceIDs_.setSize(fz.size(), -1);

    label nBlockedFaces = 0;
    forAll(fz, localFacei)
    {
        const label facei = fz[localFacei];

        if (mesh_.isInternalFace(facei))
        {
            if (alphaf[facei] > alphaThreshold_)
            {
                blockedFaces[localFacei] = true;
            }
        }
        else
        {
            label patchi = mesh_.boundaryMesh().whichPatch(facei);
            label patchFacei = -1;

            const polyPatch& pp = mesh_.boundaryMesh()[patchi];
            const scalarField& alphafp = alphaf.boundaryField()[patchi];
            const auto* cpp = isA<coupledPolyPatch>(pp);

            if (cpp)
            {
                patchFacei = (cpp->owner() ? pp.whichFace(facei) : -1);
            }
            else if (!isA<emptyPolyPatch>(pp))
            {
                patchFacei = pp.whichFace(facei);
            }

            if (patchFacei == -1)
            {
                patchi = -1;
            }
            else if (alphafp[patchFacei] > alphaThreshold_)
            {
                blockedFaces[localFacei] = true;
            }

            patchIDs_[localFacei] = patchi;
            patchFaceIDs_[localFacei] = patchFacei;
        }
    }

    DebugInFunction << "Number of blocked faces: " << nBlockedFaces << endl;
}


void Foam::functionObjects::extractEulerianParticles::collectParticle
(
    const scalar time,
    const label regioni
)
{
    DebugInFunction << "collectParticle: " << regioni << endl;

    const label particlei = regionToParticleMap_[regioni];
    eulerianParticle p = particles_[particlei];

    if (p.faceIHit != -1 && nInjectorLocations_)
    {
        // Use coarse face index for tag output
        label coarseFacei = fineToCoarseAddr_[p.faceIHit];
        p.faceIHit = globalCoarseFaces_.toGlobal(coarseFacei);
    }

    reduce(p, sumParticleOp<eulerianParticle>());

    const scalar pDiameter = cbrt(6.0*p.V/constant::mathematical::pi);

    if ((pDiameter > minDiameter_) && (pDiameter < maxDiameter_))
    {
        if (Pstream::master())
        {
            const scalar d = cbrt(6.0*p.V/constant::mathematical::pi);
            const point position = p.VC/(p.V + ROOTVSMALL);
            const vector U = p.VU/(p.V + ROOTVSMALL);
            label tag = -1;
            if (nInjectorLocations_)
            {
                tag = p.faceIHit;
            }

            injectedParticle* ip = new injectedParticle
            (
                mesh_,
                position,
                tag,
                time,
                d,
                U,
                false // not looking to set cell owner etc.
            );

            cloud_.addParticle(ip);

            collectedVolume_ += p.V;
        }

        ++nCollectedParticles_;
    }
    else
    {
        // Discard particles over/under diameter threshold
        ++nDiscardedParticles_;
        discardedVolume_ += p.V;
    }
}


void Foam::functionObjects::extractEulerianParticles::calculateAddressing
(
    const label nNewRegions,
    const scalar time,
    labelList& regionFaceIDs
)
{
    DebugInFunction << endl;

    // Determine mapping between old and new regions so that we can
    // accumulate particle info
    labelList oldToNewRegion(particles_.size(), -1);
    labelList newToNewRegion(identity(nNewRegions));

    forAll(regionFaceIDs, facei)
    {
        label newRegioni = regionFaceIDs[facei];
        label oldRegioni = regions0_[facei];

        if (newRegioni != -1 && oldRegioni != -1)
        {
            // If old region has split into multiple regions we need to
            // renumber new regions to maintain connectivity with old regions
            newToNewRegion[newRegioni] =
                max(newRegioni, oldToNewRegion[oldRegioni]);
            oldToNewRegion[oldRegioni] = newRegioni;
        }
    }

    // Create map from new regions to slots in particles list
    // - filter through new-to-new addressing to identify new particles
    Pstream::listCombineReduce(newToNewRegion, maxEqOp<label>());

    label nParticle = -1;
    labelHashSet newRegions;
    Map<label> newRegionToParticleMap;
    forAll(newToNewRegion, newRegioni0)
    {
        label newRegioni = newToNewRegion[newRegioni0];
        if (newRegions.insert(newRegioni))
        {
            ++nParticle;
        }

        // New particle slot
        newRegionToParticleMap.insert(newRegioni0, nParticle);
    }

    // Accumulate old region data or create a new particle if there is no
    // mapping from the old-to-new region
    Pstream::listCombineReduce(oldToNewRegion, maxEqOp<label>());

    List<eulerianParticle> newParticles(newRegionToParticleMap.size());
    forAll(oldToNewRegion, oldRegioni)
    {
        label newRegioni = oldToNewRegion[oldRegioni];
        if (newRegioni == -1)
        {
            // No mapping from old-to-new - collect new particle
            DebugInfo
                << "Collecting particle from oldRegion:" << oldRegioni
                << endl;

            collectParticle(time, oldRegioni);
        }
        else
        {
            // Combine existing particle into new particle
            label newParticlei = newRegionToParticleMap[newRegioni];
            label oldParticlei = regionToParticleMap_[oldRegioni];

            DebugInfo
                << "Combining newRegioni: " << newRegioni
                << "(p:" << newParticlei << ") and "
                << "oldRegioni: " << oldRegioni
                << "(p:" << oldParticlei << ")"
                << endl;

            newParticles[newParticlei] =
                sumParticleOp<eulerianParticle>()
                (
                    newParticles[newParticlei],
                    particles_[oldParticlei]
                );
        }
    }

    // Reset the particles list and addressing for latest available info
    particles_.transfer(newParticles);
    regionToParticleMap_ = newRegionToParticleMap;

    // Reset the region IDs for the next integration step
    // - these become the oldRegioni's
    regions0_ = regionFaceIDs;
}


void Foam::functionObjects::extractEulerianParticles::accumulateParticleInfo
(
    const surfaceScalarField& alphaf,
    const surfaceScalarField& phi,
    const labelList& regionFaceIDs,
    const faceZone& fz
)
{
    DebugInFunction << endl;

    const volVectorField& U = mesh_.lookupObject<volVectorField>(UName_);
    const surfaceVectorField Uf(fvc::interpolate(U));

    const scalar deltaT = mesh_.time().deltaTValue();
    const pointField& faceCentres = mesh_.faceCentres();

    forAll(regionFaceIDs, localFacei)
    {
        const label newRegioni = regionFaceIDs[localFacei];
        if (newRegioni != -1)
        {
            const label particlei = regionToParticleMap_[newRegioni];
            const label meshFacei = fz[localFacei];
            eulerianParticle& p = particles_[particlei];

            if (p.faceIHit < 0)
            {
                // New particle - does not exist in particles_ list
                p.faceIHit = localFacei;
                p.V = 0;
                p.VC = vector::zero;
                p.VU = vector::zero;
            }

            // Accumulate particle properties
            scalar magPhii = mag(faceValue(phi, localFacei, meshFacei));
            vector Ufi = faceValue(Uf, localFacei, meshFacei);
            scalar dV = magPhii*deltaT;
            p.V += dV;
            p.VC += dV*faceCentres[meshFacei];
            p.VU += dV*Ufi;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

Foam::functionObjects::extractEulerianParticles::extractEulerianParticles
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(runTime, name),
    cloud_(mesh_, "eulerianParticleCloud"),
    faceZoneName_(),
    zoneID_(-1),
    patchIDs_(),
    patchFaceIDs_(),
    alphaName_("alpha"),
    alphaThreshold_(0.1),
    UName_("U"),
    rhoName_("rho"),
    phiName_("phi"),
    nInjectorLocations_(0),
    fineToCoarseAddr_(),
    globalCoarseFaces_(),
    regions0_(),
    particles_(),
    regionToParticleMap_(),
    minDiameter_(ROOTVSMALL),
    maxDiameter_(GREAT),
    nCollectedParticles_(getProperty<label>("nCollectedParticles", 0)),
    collectedVolume_(getProperty<scalar>("collectedVolume", 0)),
    nDiscardedParticles_(getProperty<label>("nDiscardedParticles", 0)),
    discardedVolume_(getProperty<scalar>("discardedVolume", 0))
{
    if (mesh_.nSolutionD() != 3)
    {
        FatalErrorInFunction
            << name << " function object only applicable to 3-D cases"
            << exit(FatalError);
    }

    read(dict);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::functionObjects::extractEulerianParticles::read
(
    const dictionary& dict
)
{
    DebugInFunction << endl;

    if (fvMeshFunctionObject::read(dict) && writeFile::read(dict))
    {
        dict.readEntry("faceZone", faceZoneName_);
        dict.readEntry("alpha", alphaName_);

        dict.readIfPresent("alphaThreshold", alphaThreshold_);
        dict.readIfPresent("U", UName_);
        dict.readIfPresent("rho", rhoName_);
        dict.readIfPresent("phi", phiName_);
        dict.readIfPresent("nLocations", nInjectorLocations_);
        dict.readIfPresent("minDiameter", minDiameter_);
        dict.readIfPresent("maxDiameter", maxDiameter_);

        checkFaceZone();

        if (nInjectorLocations_)
        {
            initialiseBins();
        }

        return true;
    }

    return false;
}


bool Foam::functionObjects::extractEulerianParticles::execute()
{
    DebugInFunction << endl;

    Log << type() << " " << name() << " output:" << nl;

    const volScalarField& alpha =
        mesh_.lookupObject<volScalarField>(alphaName_);

    auto talphaf = surfaceScalarField::New
    (
        IOobject::scopedName(typeName, "alphaf"),
        IOobject::NO_REGISTER,
        fvc::interpolate(alpha)
    );

    const faceZone& fz = mesh_.faceZones()[zoneID_];
    const indirectPrimitivePatch patch
    (
        IndirectList<face>(mesh_.faces(), fz),
        mesh_.points()
    );

    // Set the blocked faces, i.e. where alpha > alpha threshold value
    boolList blockedFaces(fz.size(), false);
    setBlockedFaces(talphaf(), fz, blockedFaces);

    // Split the  faceZone according to the blockedFaces
    // - Returns a list of (disconnected) region index per face zone face
    regionSplit2D regionFaceIDs(mesh_, patch, blockedFaces);

    // Global number of regions
    const label nRegionsNew = regionFaceIDs.nRegions();

    // Calculate the addressing between the old and new region information
    // Also collects particles that have traversed the faceZone
    // - Note: may also update regionFaceIDs
    calculateAddressing
    (
        nRegionsNew,
        mesh_.time().value(),
        regionFaceIDs
    );

    // Process latest region information
    tmp<surfaceScalarField> tphi = phiU();
    accumulateParticleInfo(talphaf(), tphi(), regionFaceIDs, fz);

    Log << "    Collected particles   : " << nCollectedParticles_ << nl
        << "    Collected volume      : " << collectedVolume_ << nl
        << "    Discarded particles   : " << nDiscardedParticles_ << nl
        << "    Discarded volume      : " << discardedVolume_ << nl
        << "    Particles in progress : " << particles_.size() << nl
        << endl;

    return true;
}


bool Foam::functionObjects::extractEulerianParticles::write()
{
    DebugInFunction << endl;

    cloud_.write();

    setProperty("nCollectedParticles", nCollectedParticles_);
    setProperty("collectedVolume", collectedVolume_);
    setProperty("nDiscardedParticles", nDiscardedParticles_);
    setProperty("discardedVolume", discardedVolume_);

    return true;
}


// ************************************************************************* //
