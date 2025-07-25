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

#include "faceSetOption.H"
#include "faceSet.H"
#include "areaFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace fa
    {
        defineTypeNameAndDebug(faceSetOption, 0);
    }
}


const Foam::Enum
<
    Foam::fa::faceSetOption::selectionModeType
>
Foam::fa::faceSetOption::selectionModeTypeNames_
({
    { selectionModeType::smAll, "all" },
    { selectionModeType::smFaceSet, "faceSet" },
    { selectionModeType::smFaceZone, "faceZone" },
    { selectionModeType::smPatch, "patch" },
    { selectionModeType::smFaceZone, "volFaceZone" }  // Compat?
});


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fa::faceSetOption::setSelection(const dictionary& dict)
{
    selectionNames_.clear();

    switch (selectionMode_)
    {
        case smAll:
        {
            break;
        }
        case smFaceSet:
        {
            selectionNames_.resize(1);
            dict.readEntry("faceSet", selectionNames_.first());
            break;
        }
        case smFaceZone:
        {
            if
            (
                !dict.readIfPresent("faceZones", selectionNames_)
             || selectionNames_.empty()
            )
            {
                selectionNames_.resize(1);
                dict.readEntry("faceZone", selectionNames_.first());
            }
            break;
        }
        case smPatch:
        {
            if
            (
                !dict.readIfPresent("patches", selectionNames_)
             || selectionNames_.empty()
            )
            {
                selectionNames_.resize(1);
                dict.readEntry("patch", selectionNames_.first());
            }
            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unknown selectionMode "
                << selectionModeTypeNames_[selectionMode_]
                << ". Valid selectionMode types : "
                << selectionModeTypeNames_
                << exit(FatalError);
        }
    }
}


void Foam::fa::faceSetOption::setArea()
{
    // Set area information

    scalar sumArea = 0;
    for (const label facei : faces_)
    {
        sumArea += regionMesh().S()[facei];
    }
    reduce(sumArea, sumOp<scalar>());

    const scalar old(A_);
    A_ = sumArea;

    // Compare area values, stringified using current write precision
    if
    (
        Time::timeName(old, IOstream::defaultPrecision())
     != Time::timeName(A_, IOstream::defaultPrecision())
    )
    {
        Info<< indent
            << "- selected " << returnReduce(faces_.size(), sumOp<label>())
            << " face(s) with area " << A_ << endl;
    }
}


void Foam::fa::faceSetOption::setFaceSelection()
{
    switch (selectionMode_)
    {
        case smAll:
        {
            Info<< indent << "- selecting all faces" << endl;
            faces_ = identity(regionMesh().nFaces());

            break;
        }

        case smFaceSet:
        {
            Info<< indent
                << "- selecting face subset using volume-mesh faceSet "
                << zoneName() << nl;

            const faceSet subset(mesh_, zoneName());

            const labelUList& faceLabels = regionMesh().faceLabels();

            faces_.resize_nocopy(faceLabels.size());

            label nUsed = 0;
            forAll(faceLabels, facei)
            {
                const label meshFacei = faceLabels[facei];

                if (subset.test(meshFacei))
                {
                    faces_[nUsed] = facei;
                    ++nUsed;
                }
            }
            faces_.resize(nUsed);
            break;
        }

        case smFaceZone:
        {
            Info<< indent
                << "- selecting face subset using volume-mesh faceZones "
                << flatOutput(selectionNames_) << nl;

            const auto& zones = mesh_.faceZones();

            // Also handles groups, multiple zones etc ...
            labelList zoneIDs = zones.indices(selectionNames_);

            if (zoneIDs.empty())
            {
                FatalErrorInFunction
                    << "No matching faceZones: "
                    << flatOutput(selectionNames_) << nl
                    << "Valid zones : "
                    << flatOutput(zones.names()) << nl
                    << "Valid groups: "
                    << flatOutput(zones.groupNames())
                    << nl
                    << exit(FatalError);
            }

            const bitSet subset(mesh_.faceZones().selection(zoneIDs));

            const labelUList& faceLabels = regionMesh().faceLabels();

            faces_.resize_nocopy(faceLabels.size());

            label nUsed = 0;
            forAll(faceLabels, facei)
            {
                const label meshFacei = faceLabels[facei];

                if (subset.test(meshFacei))
                {
                    faces_[nUsed] = facei;
                    ++nUsed;
                }
            }
            faces_.resize(nUsed);
            break;
        }

        case smPatch:
        {
            Info<< indent
                << "- selecting face subset using volume-mesh patches "
                << flatOutput(selectionNames_) << nl;

            const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

            // Also handles groups, multiple patches etc ...
            labelList patchIDs = pbm.indices(selectionNames_);

            if (patchIDs.empty())
            {
                FatalErrorInFunction
                    << "No matching patches: "
                    << flatOutput(selectionNames_) << nl
                    << "Valid patches : "
                    << flatOutput(pbm.names()) << nl
                    << "Valid groups: "
                    << flatOutput(pbm.groupNames()) << nl
                    << exit(FatalError);
            }

            const List<labelPair>& patchFaces = regionMesh().whichPatchFaces();

            faces_.resize_nocopy(patchFaces.size());

            label nUsed = 0;
            forAll(patchFaces, facei)
            {
                const label patchi = patchFaces[facei].first();

                if (patchIDs.found(patchi))
                {
                    faces_[nUsed] = facei;
                    ++nUsed;
                }
            }
            faces_.resize(nUsed);
            break;
        }

        default:
        {
            FatalErrorInFunction
                << "Unknown selectionMode "
                << selectionModeTypeNames_[selectionMode_]
                << ". Valid selectionMode types are "
                << selectionModeTypeNames_
                << exit(FatalError);
        }
    }

    if (smAll != selectionMode_ && returnReduceAnd(faces_.empty()))
    {
        WarningInFunction
            << "No faces selected!" << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fa::faceSetOption::faceSetOption
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fa::option(name, modelType, dict, mesh),
    timeStart_(-1),
    duration_(0),
    selectionMode_(selectionModeTypeNames_.get("selectionMode", coeffs_)),
    selectionNames_(),
    A_(0)
{
    if (isActive())
    {
        Info<< incrIndent;
        read(dict);
        setSelection(coeffs_);
        setFaceSelection();
        setArea();
        Info<< decrIndent;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fa::faceSetOption::isActive()
{
    if (fa::option::isActive() && inTimeLimits(mesh_.time().value()))
    {
        // Update the face set if the mesh is changing
        if (mesh_.changing())
        {
            if (mesh_.topoChanging())
            {
                setArea();
                // Force printing of new set area
                A_ = -GREAT;
            }

            // Report new area (if changed)
            setArea();
        }

        return true;
    }

    return false;
}


bool Foam::fa::faceSetOption::read(const dictionary& dict)
{
    if (fa::option::read(dict))
    {
        timeStart_ = -1;

        if (coeffs_.readIfPresent("timeStart", timeStart_))
        {
            coeffs_.readEntry("duration", duration_);
        }

        return true;
    }

    return false;
}


// ************************************************************************* //
