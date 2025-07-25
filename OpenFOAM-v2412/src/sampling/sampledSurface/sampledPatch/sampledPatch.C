/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2021 OpenCFD Ltd.
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

#include "sampledPatch.H"
#include "dictionary.H"
#include "polyMesh.H"
#include "polyPatch.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "indirectPrimitivePatch.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sampledPatch, 0);
    addNamedToRunTimeSelectionTable(sampledSurface, sampledPatch, word, patch);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledPatch::sampledPatch
(
    const word& name,
    const polyMesh& mesh,
    const UList<wordRe>& patchNames,
    const bool triangulate
)
:
    sampledSurface(name, mesh),
    selectionNames_(patchNames),
    triangulate_(triangulate),
    needsUpdate_(true)
{}


Foam::sampledPatch::sampledPatch
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    sampledSurface(name, mesh, dict),
    selectionNames_(dict.get<wordRes>("patches")),
    triangulate_(dict.getOrDefault("triangulate", false)),
    needsUpdate_(true)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::labelList& Foam::sampledPatch::patchIDs() const
{
    if (patchIDs_.empty())
    {
        labelList selected
        (
            mesh().boundaryMesh().patchSet(selectionNames_).sortedToc()
        );

        DynamicList<label> bad;
        for (const label patchi : selected)
        {
            const polyPatch& pp = mesh().boundaryMesh()[patchi];

            if (isA<emptyPolyPatch>(pp))
            {
                bad.append(patchi);
            }
        }

        if (bad.size())
        {
            label nGood = (selected.size() - bad.size());

            auto& os = nGood > 0 ? WarningInFunction : FatalErrorInFunction;

            os  << "Cannot sample an empty patch" << nl;

            for (const label patchi : bad)
            {
                os  << "    "
                    << mesh().boundaryMesh()[patchi].name() << nl;
            }

            if (nGood)
            {
                os  << "No non-empty patches selected" << endl
                    << exit(FatalError);
            }
            else
            {
                os  << "Selected " << nGood << " non-empty patches" << nl;
            }

            patchIDs_.resize(nGood);
            nGood = 0;
            for (const label patchi : selected)
            {
                if (!bad.found(patchi))
                {
                    patchIDs_[nGood] = patchi;
                    ++nGood;
                }
            }
        }
        else
        {
            patchIDs_ = std::move(selected);
        }
    }

    return patchIDs_;
}


bool Foam::sampledPatch::needsUpdate() const
{
    return needsUpdate_;
}


bool Foam::sampledPatch::expire()
{
    // already marked as expired
    if (needsUpdate_)
    {
        return false;
    }

    sampledSurface::clearGeom();
    Mesh::clear();

    patchIDs_.clear();
    patchStart_.clear();

    patchIndex_.clear();
    patchFaceLabels_.clear();

    needsUpdate_ = true;
    return true;
}


bool Foam::sampledPatch::update()
{
    if (!needsUpdate_)
    {
        return false;
    }

    // Total number of faces selected
    label numFaces = 0;
    for (const label patchi : patchIDs())
    {
        const polyPatch& pp = mesh().boundaryMesh()[patchi];
        numFaces += pp.size();
    }

    patchStart_.resize(patchIDs().size());

    // The originating patch and local face in the patch.
    patchIndex_.resize(numFaces);
    patchFaceLabels_.resize(numFaces);

    IndirectList<face> selectedFaces(mesh().faces(), labelList());
    labelList& meshFaceIds = selectedFaces.addressing();
    meshFaceIds.resize(numFaces);

    numFaces = 0;

    forAll(patchIDs(), idx)
    {
        const label patchi = patchIDs()[idx];
        const polyPatch& pp = mesh().boundaryMesh()[patchi];
        const label len = pp.size();

        patchStart_[idx] = numFaces;

        SubList<label>(patchIndex_, len, numFaces) = idx;

        SubList<label>(patchFaceLabels_, len, numFaces) = identity(len);

        SubList<label>(meshFaceIds, len, numFaces) = identity(len, pp.start());

        numFaces += len;
    }


    uindirectPrimitivePatch allPatches(selectedFaces, mesh().points());

    this->storedPoints() = allPatches.localPoints();
    this->storedFaces()  = allPatches.localFaces();


    // triangulate uses remapFaces()
    // - this is somewhat less efficient since it recopies the faces
    // that we just created, but we probably don't want to do this
    // too often anyhow.
    if (triangulate_)
    {
        Mesh::triangulate();
    }

    if (debug)
    {
        print(Pout, debug);
        Pout<< endl;
    }

    needsUpdate_ = false;
    return true;
}


// remap action on triangulation
void Foam::sampledPatch::remapFaces(const labelUList& faceMap)
{
    if (!faceMap.empty())
    {
        Mesh::remapFaces(faceMap);
        patchFaceLabels_ = labelList
        (
            labelUIndList(patchFaceLabels_, faceMap)
        );
        patchIndex_ = labelList
        (
            labelUIndList(patchIndex_, faceMap)
        );

        // Update patchStart
        if (patchIndex_.size())
        {
            patchStart_[patchIndex_[0]] = 0;
            for (label i = 1; i < patchIndex_.size(); ++i)
            {
                if (patchIndex_[i] != patchIndex_[i-1])
                {
                    patchStart_[patchIndex_[i]] = i;
                }
            }
        }
    }
}


Foam::tmp<Foam::scalarField> Foam::sampledPatch::sample
(
    const interpolation<scalar>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::vectorField> Foam::sampledPatch::sample
(
    const interpolation<vector>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::sphericalTensorField> Foam::sampledPatch::sample
(
    const interpolation<sphericalTensor>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::symmTensorField> Foam::sampledPatch::sample
(
    const interpolation<symmTensor>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::tensorField> Foam::sampledPatch::sample
(
    const interpolation<tensor>& sampler
) const
{
    return sampleOnFaces(sampler);
}


bool Foam::sampledPatch::withSurfaceFields() const
{
    return true;
}


Foam::tmp<Foam::scalarField> Foam::sampledPatch::sample
(
    const surfaceScalarField& sField
) const
{
    return sampleOnFaces(sField);
}


Foam::tmp<Foam::vectorField> Foam::sampledPatch::sample
(
    const surfaceVectorField& sField
) const
{
    return sampleOnFaces(sField);
}


Foam::tmp<Foam::sphericalTensorField> Foam::sampledPatch::sample
(
    const surfaceSphericalTensorField& sField
) const
{
    return sampleOnFaces(sField);
}


Foam::tmp<Foam::symmTensorField> Foam::sampledPatch::sample
(
    const surfaceSymmTensorField& sField
) const
{
    return sampleOnFaces(sField);
}


Foam::tmp<Foam::tensorField> Foam::sampledPatch::sample
(
    const surfaceTensorField& sField
) const
{
    return sampleOnFaces(sField);
}


Foam::tmp<Foam::scalarField> Foam::sampledPatch::interpolate
(
    const interpolation<scalar>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


Foam::tmp<Foam::vectorField> Foam::sampledPatch::interpolate
(
    const interpolation<vector>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


Foam::tmp<Foam::sphericalTensorField> Foam::sampledPatch::interpolate
(
    const interpolation<sphericalTensor>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


Foam::tmp<Foam::symmTensorField> Foam::sampledPatch::interpolate
(
    const interpolation<symmTensor>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


Foam::tmp<Foam::tensorField> Foam::sampledPatch::interpolate
(
    const interpolation<tensor>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


void Foam::sampledPatch::print(Ostream& os, int level) const
{
    os  << "sampledPatch: " << name() << " :"
        << " patches:" << flatOutput(selectionNames_);

    if (level)
    {
        os  << "  faces:" << faces().size()
            << "  points:" << points().size();
    }
}


// ************************************************************************* //
