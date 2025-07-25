/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2022 OpenCFD Ltd.
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

#include "sampledFaceZone.H"
#include "dictionary.H"
#include "polyMesh.H"
#include "polyPatch.H"
#include "processorPolyPatch.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "volPointInterpolation.H"
#include "indirectPrimitivePatch.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sampledFaceZone, 0);
    addNamedToRunTimeSelectionTable
    (
        sampledSurface,
        sampledFaceZone,
        word,
        faceZone
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledFaceZone::sampledFaceZone
(
    const word& name,
    const polyMesh& mesh,
    const UList<wordRe>& zoneNames,
    const bool triangulate
)
:
    sampledSurface(name, mesh),
    selectionNames_(zoneNames),
    triangulate_(triangulate),
    needsUpdate_(true)
{}


Foam::sampledFaceZone::sampledFaceZone
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    sampledSurface(name, mesh, dict),
    selectionNames_(dict.get<wordRes>("zones")),
    triangulate_(dict.getOrDefault("triangulate", false)),
    needsUpdate_(true)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::labelList& Foam::sampledFaceZone::zoneIDs() const
{
    if (zoneIds_.empty())
    {
        // Zone indices for all matches, already sorted
        zoneIds_ = mesh().faceZones().indices(selectionNames_);
    }

    return zoneIds_;
}


bool Foam::sampledFaceZone::needsUpdate() const
{
    return needsUpdate_;
}


bool Foam::sampledFaceZone::expire()
{
    // already marked as expired
    if (needsUpdate_)
    {
        return false;
    }

    sampledSurface::clearGeom();
    Mesh::clear();

    zoneIds_.clear();

    faceId_.clear();
    facePatchId_.clear();

    needsUpdate_ = true;
    return true;
}


bool Foam::sampledFaceZone::update()
{
    if (!needsUpdate_)
    {
        return false;
    }

    // Total number of faces selected
    label numFaces = 0;
    for (const label zonei : zoneIDs())
    {
        numFaces += mesh().faceZones()[zonei].size();
    }

    if (zoneIDs().empty())
    {
        WarningInFunction
            << type() << ' ' << name() << ": "
            << "    No matching face zone(s): "
            << flatOutput(selectionNames_)  << nl
            << "    Known face zones: "
            << flatOutput(mesh().faceZones().names()) << nl;
    }

    // Could also check numFaces

    // The mesh face or local patch face and the patch id
    faceId_.resize_nocopy(numFaces);
    facePatchId_.resize_nocopy(numFaces);

    IndirectList<face> selectedFaces(mesh().faces(), labelList());
    labelList& meshFaceIds = selectedFaces.addressing();
    meshFaceIds.resize_nocopy(numFaces);

    numFaces = 0;

    for (const label zoneId : zoneIDs())
    {
        const faceZone& fZone = mesh().faceZones()[zoneId];

        for (const label meshFacei : fZone)
        {
            // Internal faces
            label faceId = meshFacei;
            label facePatchId = -1;

            // Boundary faces
            if (!mesh().isInternalFace(meshFacei))
            {
                facePatchId = mesh().boundaryMesh().whichPatch(meshFacei);
                const polyPatch& pp = mesh().boundaryMesh()[facePatchId];

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
                meshFaceIds[numFaces] = meshFacei;

                ++numFaces;
            }
        }
    }

    // Shrink to size used
    faceId_.resize(numFaces);
    facePatchId_.resize(numFaces);
    meshFaceIds.resize(numFaces);

    uindirectPrimitivePatch zoneFaces(selectedFaces, mesh().points());

    this->storedPoints() = zoneFaces.localPoints();
    this->storedFaces()  = zoneFaces.localFaces();

    // triangulate - uses remapFaces()
    if (triangulate_)
    {
        Mesh::triangulate();
    }

    needsUpdate_ = false;
    return true;
}


// remap action on triangulation
void Foam::sampledFaceZone::remapFaces(const labelUList& faceMap)
{
    if (!faceMap.empty())
    {
        Mesh::remapFaces(faceMap);
        faceId_ = labelList
        (
            labelUIndList(faceId_, faceMap)
        );
        facePatchId_ = labelList
        (
            labelUIndList(facePatchId_, faceMap)
        );
    }
}


Foam::tmp<Foam::scalarField> Foam::sampledFaceZone::sample
(
    const interpolation<scalar>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::vectorField> Foam::sampledFaceZone::sample
(
    const interpolation<vector>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::sphericalTensorField> Foam::sampledFaceZone::sample
(
    const interpolation<sphericalTensor>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::symmTensorField> Foam::sampledFaceZone::sample
(
    const interpolation<symmTensor>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::tensorField> Foam::sampledFaceZone::sample
(
    const interpolation<tensor>& sampler
) const
{
    return sampleOnFaces(sampler);
}


bool Foam::sampledFaceZone::withSurfaceFields() const
{
    return true;
}


Foam::tmp<Foam::scalarField> Foam::sampledFaceZone::sample
(
    const surfaceScalarField& sField
) const
{
    return sampleOnFaces(sField);
}


Foam::tmp<Foam::vectorField> Foam::sampledFaceZone::sample
(
    const surfaceVectorField& sField
) const
{
    return sampleOnFaces(sField);
}


Foam::tmp<Foam::sphericalTensorField> Foam::sampledFaceZone::sample
(
    const surfaceSphericalTensorField& sField
) const
{
    return sampleOnFaces(sField);
}


Foam::tmp<Foam::symmTensorField> Foam::sampledFaceZone::sample
(
    const surfaceSymmTensorField& sField
) const
{
    return sampleOnFaces(sField);
}


Foam::tmp<Foam::tensorField> Foam::sampledFaceZone::sample
(
    const surfaceTensorField& sField
) const
{
    return sampleOnFaces(sField);
}


Foam::tmp<Foam::scalarField> Foam::sampledFaceZone::interpolate
(
    const interpolation<scalar>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


Foam::tmp<Foam::vectorField> Foam::sampledFaceZone::interpolate
(
    const interpolation<vector>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}

Foam::tmp<Foam::sphericalTensorField>
Foam::sampledFaceZone::interpolate
(
    const interpolation<sphericalTensor>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


Foam::tmp<Foam::symmTensorField> Foam::sampledFaceZone::interpolate
(
    const interpolation<symmTensor>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


Foam::tmp<Foam::tensorField> Foam::sampledFaceZone::interpolate
(
    const interpolation<tensor>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


void Foam::sampledFaceZone::print(Ostream& os, int level) const
{
    os  << "faceZone: " << name() << " :"
        << " zones:" << flatOutput(selectionNames_);

    if (level)
    {
        os  << "  faces:" << faces().size()
            << "  points:" << points().size();
    }
}


// ************************************************************************* //
