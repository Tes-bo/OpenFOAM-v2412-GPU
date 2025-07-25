/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2021 OpenCFD Ltd.
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

#include "sampledIsoSurface.H"
#include "dictionary.H"
#include "fvMesh.H"
#include "volFields.H"
#include "volPointInterpolation.H"
#include "addToRunTimeSelectionTable.H"
#include "PtrList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sampledIsoSurface, 0);
    addNamedToRunTimeSelectionTable
    (
        sampledSurface,
        sampledIsoSurface,
        word,
        isoSurface
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sampledIsoSurface::getIsoFields() const
{
    const fvMesh& fvm = static_cast<const fvMesh&>(mesh());

    // Get volField
    // ~~~~~~~~~~~~

    volFieldPtr_ = fvm.getObjectPtr<volScalarField>(isoField_);

    if (volFieldPtr_)
    {
        DebugInFunction
            << "Lookup volField " << isoField_ << endl;

        storedVolFieldPtr_.clear();
    }
    else
    {
        // Bit of a hack. Read field and store.

        DebugInFunction
            << "Checking " << isoField_
            << " for same time " << fvm.time().timeName() << endl;

        if
        (
            !storedVolFieldPtr_
         || (fvm.time().timeName() != storedVolFieldPtr_().instance())
        )
        {
            DebugInFunction
                << "Reading volField " << isoField_
                << " from time " << fvm.time().timeName() << endl;

            IOobject vfHeader
            (
                isoField_,
                fvm.time().timeName(),
                fvm,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            );

            if (vfHeader.typeHeaderOk<volScalarField>(true))
            {
                storedVolFieldPtr_.reset
                (
                    new volScalarField
                    (
                        vfHeader,
                        fvm
                    )
                );
                volFieldPtr_ = storedVolFieldPtr_.get(); // get(), not release()
            }
            else
            {
                FatalErrorInFunction
                    << "Cannot find isosurface field " << isoField_
                    << " in database or directory " << vfHeader.path()
                    << exit(FatalError);
            }
        }
    }


    // Get pointField
    // ~~~~~~~~~~~~~~

    // In case of multiple iso values we don't want to calculate multiple e.g.
    // "volPointInterpolate(p)" so register it and re-use it. This is the
    // same as the 'cache' functionality from volPointInterpolate but
    // unfortunately that one does not guarantee that the field pointer
    // remain: e.g. some other functionObject might delete the cached version.
    // (volPointInterpolation::interpolate with cache=false deletes any
    //  registered one or if mesh.changing())

    if (!subMeshPtr_)
    {
        const word pointFldName =
            "volPointInterpolate_"
          + type()
          + "("
          + isoField_
          + ')';


        pointFieldPtr_ = fvm.getObjectPtr<pointScalarField>(pointFldName);

        if (pointFieldPtr_)
        {
            DebugInFunction
                << "lookup pointField " << pointFldName << endl;

            if (!pointFieldPtr_->upToDate(*volFieldPtr_))
            {
                DebugInFunction
                    << "updating pointField " << pointFldName << endl;

                // Update the interpolated value
                volPointInterpolation::New(fvm).interpolate
                (
                    *volFieldPtr_,
                    const_cast<pointScalarField&>(*pointFieldPtr_)
                );
            }
        }
        else
        {
            // Not in registry. Interpolate.

            DebugInFunction
                << "creating pointField " << pointFldName << endl;

            // Interpolate without cache. Note that we're registering it
            // below so next time round it goes into the condition
            // above.
            pointFieldPtr_ =
                volPointInterpolation::New(fvm).interpolate
                (
                    *volFieldPtr_,
                    pointFldName,
                    false
                ).ptr();

            const_cast<pointScalarField*>(pointFieldPtr_)->store();
        }


        // If averaging redo the volField.
        // Can only be done now since needs the point field.
        if (average_)
        {
            storedVolFieldPtr_.reset
            (
                pointAverage(*pointFieldPtr_).ptr()
            );
            volFieldPtr_ = storedVolFieldPtr_.get(); // get(), not release()
        }


        DebugInFunction
            << "volField " << volFieldPtr_->name()
            << " min:" << min(*volFieldPtr_).value()
            << " max:" << max(*volFieldPtr_).value() << nl
            << "pointField " << pointFieldPtr_->name()
            << " min:" << gMin(pointFieldPtr_->primitiveField())
            << " max:" << gMax(pointFieldPtr_->primitiveField()) << endl;
    }
    else
    {
        // Get subMesh variants
        const fvMesh& subFvm = subMeshPtr_().subMesh();

        // Either lookup on the submesh or subset the whole-mesh volField

        volSubFieldPtr_ = subFvm.getObjectPtr<volScalarField>(isoField_);

        if (volSubFieldPtr_)
        {
            DebugInFunction
                << "Sub-mesh lookup volField " << isoField_ << endl;

            storedVolSubFieldPtr_.clear();
        }
        else
        {
            DebugInFunction
                << "Sub-setting volField " << isoField_ << endl;

            storedVolSubFieldPtr_.reset
            (
                subMeshPtr_().interpolate
                (
                    *volFieldPtr_
                ).ptr()
            );
            storedVolSubFieldPtr_->checkOut();
            volSubFieldPtr_ = storedVolSubFieldPtr_.get(); // not release()
        }


        // The point field on subMesh

        const word pointFldName =
            "volPointInterpolate_"
          + type()
          + "("
          + volSubFieldPtr_->name()
          + ')';


        pointFieldPtr_ = subFvm.getObjectPtr<pointScalarField>(pointFldName);

        if (pointFieldPtr_)
        {
            DebugInFunction
                << "Sub-mesh lookup pointField " << pointFldName << endl;

            if (!pointFieldPtr_->upToDate(*volSubFieldPtr_))
            {
                DebugInFunction
                    << "Updating submesh pointField " << pointFldName << endl;

                // Update the interpolated value
                volPointInterpolation::New(subFvm).interpolate
                (
                    *volSubFieldPtr_,
                    const_cast<pointScalarField&>(*pointFieldPtr_)
                );
            }
        }
        else
        {
            DebugInFunction
                << "Interpolating submesh volField "
                << volSubFieldPtr_->name()
                << " to get submesh pointField " << pointFldName << endl;

            pointSubFieldPtr_ =
                volPointInterpolation::New
                (
                    subFvm
                ).interpolate(*volSubFieldPtr_).ptr();

            const_cast<pointScalarField*>(pointSubFieldPtr_)->store();
        }


        // If averaging redo the volField. Can only be done now since needs the
        // point field.
        if (average_)
        {
            storedVolSubFieldPtr_.reset
            (
                pointAverage(*pointSubFieldPtr_).ptr()
            );
            volSubFieldPtr_ = storedVolSubFieldPtr_.get(); // not release()
        }


        DebugInFunction
            << "volSubField "
            << volSubFieldPtr_->name()
            << " min:" << min(*volSubFieldPtr_).value()
            << " max:" << max(*volSubFieldPtr_).value() << nl
            << "pointSubField "
            << pointSubFieldPtr_->name()
            << " min:" << gMin(pointSubFieldPtr_->primitiveField())
            << " max:" << gMax(pointSubFieldPtr_->primitiveField()) << endl;
    }
}


void Foam::sampledIsoSurface::combineSurfaces
(
    PtrList<isoSurfaceBase>& isoSurfPtrs
)
{
    isoSurfacePtr_.reset(nullptr);

    // Already checked previously for ALGO_POINT, but do it again
    // - ALGO_POINT still needs fields (for interpolate)
    // The others can do straight transfer
    if
    (
        isoParams_.algorithm() == isoSurfaceParams::ALGO_POINT
     && isoSurfPtrs.size() == 1
    )
    {
        // Shift from list to autoPtr
        isoSurfacePtr_.reset(isoSurfPtrs.release(0));
    }
    else if (isoSurfPtrs.size() == 1)
    {
        autoPtr<isoSurfaceBase> ptr(isoSurfPtrs.release(0));
        auto& surf = *ptr;

        surface_.transfer(static_cast<meshedSurface&>(surf));
        meshCells_.transfer(surf.meshCells());
    }
    else
    {
        // Combine faces with point offsets
        //
        // Note: use points().size() from surface, not nPoints()
        // since there may be uncompacted dangling nodes

        label nFaces = 0, nPoints = 0;

        for (const auto& surf : isoSurfPtrs)
        {
            nFaces += surf.size();
            nPoints += surf.points().size();
        }

        faceList newFaces(nFaces);
        pointField newPoints(nPoints);
        meshCells_.resize(nFaces);

        surfZoneList newZones(isoSurfPtrs.size());

        nFaces = 0;
        nPoints = 0;
        forAll(isoSurfPtrs, surfi)
        {
            autoPtr<isoSurfaceBase> ptr(isoSurfPtrs.release(surfi));
            auto& surf = *ptr;

            SubList<face> subFaces(newFaces, surf.size(), nFaces);
            SubList<point> subPoints(newPoints, surf.points().size(), nPoints);
            SubList<label> subCells(meshCells_, surf.size(), nFaces);

            newZones[surfi] = surfZone
            (
                surfZoneIdentifier::defaultName(surfi),
                subFaces.size(),    // size
                nFaces,             // start
                surfi               // index
            );

            subFaces = surf.surfFaces();
            subPoints = surf.points();
            subCells = surf.meshCells();

            if (nPoints)
            {
                for (face& f : subFaces)
                {
                    for (label& pointi : f)
                    {
                        pointi += nPoints;
                    }
                }
            }

            nFaces += subFaces.size();
            nPoints += subPoints.size();
        }

        meshedSurface combined
        (
            std::move(newPoints),
            std::move(newFaces),
            newZones
        );

        surface_.transfer(combined);
    }

    // Addressing into the full mesh
    if (subMeshPtr_ && meshCells_.size())
    {
        meshCells_ =
            UIndirectList<label>(subMeshPtr_->cellMap(), meshCells_);
    }
}


bool Foam::sampledIsoSurface::updateGeometry() const
{
    const fvMesh& fvm = static_cast<const fvMesh&>(mesh());

    // No update needed
    if (fvm.time().timeIndex() == prevTimeIndex_)
    {
        return false;
    }

    prevTimeIndex_ = fvm.time().timeIndex();

    // Clear any previously stored topologies
    surface_.clear();
    meshCells_.clear();
    isoSurfacePtr_.reset(nullptr);

    // Clear derived data
    sampledSurface::clearGeom();

    const bool hasCellZones =
        (-1 != mesh().cellZones().findIndex(zoneNames_));

    // Geometry
    if
    (
        simpleSubMesh_
     && isoParams_.algorithm() != isoSurfaceParams::ALGO_POINT
    )
    {
        subMeshPtr_.reset(nullptr);

        // Handle cell zones as inverse (blocked) selection
        if (!ignoreCellsPtr_)
        {
            ignoreCellsPtr_.reset(new bitSet);

            if (hasCellZones)
            {
                bitSet select(mesh().cellZones().selection(zoneNames_));

                if (select.any() && !select.all())
                {
                    // From selection to blocking
                    select.flip();

                    *ignoreCellsPtr_ = std::move(select);
                }
            }
        }
    }
    else
    {
        // A standard subMesh treatment

        if (ignoreCellsPtr_)
        {
            ignoreCellsPtr_->clearStorage();
        }
        else
        {
            ignoreCellsPtr_.reset(new bitSet);
        }

        // Get sub-mesh if any
        if (!subMeshPtr_ && hasCellZones)
        {
            const label exposedPatchi =
                mesh().boundaryMesh().findPatchID(exposedPatchName_);

            DebugInfo
                << "Allocating subset of size "
                << mesh().cellZones().selection(zoneNames_).count()
                << " with exposed faces into patch "
                << exposedPatchi << endl;

            subMeshPtr_.reset
            (
                new fvMeshSubset
                (
                    fvm,
                    mesh().cellZones().selection(zoneNames_),
                    exposedPatchi
                )
            );
        }
    }


    // The fields
    getIsoFields();

    refPtr<volScalarField> tvolFld(*volFieldPtr_);
    refPtr<pointScalarField> tpointFld(*pointFieldPtr_);

    if (subMeshPtr_)
    {
        tvolFld.cref(*volSubFieldPtr_);
        tpointFld.cref(*pointSubFieldPtr_);
    }


    // Create surfaces for each iso level

    PtrList<isoSurfaceBase> isoSurfPtrs(isoValues_.size());

    forAll(isoValues_, surfi)
    {
        isoSurfPtrs.set
        (
            surfi,
            isoSurfaceBase::New
            (
                isoParams_,
                tvolFld(),
                tpointFld().primitiveField(),
                isoValues_[surfi],
                *ignoreCellsPtr_
            )
        );
    }

    // And flatten
    const_cast<sampledIsoSurface&>(*this)
        .combineSurfaces(isoSurfPtrs);


    // triangulate uses remapFaces()
    // - this is somewhat less efficient since it recopies the faces
    // that we just created, but we probably don't want to do this
    // too often anyhow.
    if
    (
        triangulate_
     && surface_.size()
     && (isoParams_.algorithm() == isoSurfaceParams::ALGO_TOPO)
    )
    {
        labelList faceMap;
        surface_.triangulate(faceMap);
        meshCells_ = UIndirectList<label>(meshCells_, faceMap)();
    }

    if (debug)
    {
        Pout<< "isoSurface::updateGeometry() : constructed iso:" << nl
            << "    field          : " << isoField_ << nl
            << "    value          : " << flatOutput(isoValues_) << nl
            << "    average        : " << Switch(average_) << nl
            << "    filter         : "
            << Switch(bool(isoParams_.filter())) << nl
            << "    bounds         : " << isoParams_.getClipBounds() << nl;
        if (subMeshPtr_)
        {
            Pout<< "    zone size      : "
                << subMeshPtr_->subMesh().nCells() << nl;
        }
        Pout<< "    points         : " << points().size() << nl
            << "    faces          : " << surface().size() << nl
            << "    cut cells      : " << meshCells().size()
            << endl;
    }

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledIsoSurface::sampledIsoSurface
(
    const isoSurfaceParams& params,
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    sampledSurface(name, mesh, dict),
    isoField_(dict.get<word>("isoField")),
    isoValues_(),
    isoParams_(dict, params),
    average_(dict.getOrDefault("average", false)),
    triangulate_(dict.getOrDefault("triangulate", false)),
    simpleSubMesh_(dict.getOrDefault("simpleSubMesh", false)),
    zoneNames_(),
    exposedPatchName_(),
    prevTimeIndex_(-1),

    surface_(),
    meshCells_(),
    isoSurfacePtr_(nullptr),

    subMeshPtr_(nullptr),
    ignoreCellsPtr_(nullptr),

    storedVolFieldPtr_(nullptr),
    volFieldPtr_(nullptr),
    pointFieldPtr_(nullptr),

    storedVolSubFieldPtr_(nullptr),
    volSubFieldPtr_(nullptr),
    pointSubFieldPtr_(nullptr)
{
    if (params.algorithm() != isoSurfaceParams::ALGO_DEFAULT)
    {
        // Forced use of specified algorithm (ignore dictionary entry)
        isoParams_.algorithm(params.algorithm());
    }

    // The isoValues or isoValue

    if (!dict.readIfPresent("isoValues", isoValues_))
    {
        isoValues_.resize(1);
        dict.readEntry("isoValue", isoValues_.first());
    }

    if (isoValues_.empty())
    {
        FatalIOErrorInFunction(dict)
            << "No isoValue or isoValues specified." << nl
            << exit(FatalIOError);
    }

    if (isoValues_.size() > 1)
    {
        const label nOrig = isoValues_.size();

        inplaceUniqueSort(isoValues_);

        if (nOrig != isoValues_.size())
        {
            IOWarningInFunction(dict)
                << "Removed non-unique isoValues" << nl;
        }
    }

    if (isoParams_.algorithm() == isoSurfaceParams::ALGO_POINT)
    {
        // Not possible for ALGO_POINT
        simpleSubMesh_ = false;

        // Previously emitted an error about using ALGO_POINT with
        // non-interpolated, but that was before we had "sampleScheme"
        // at the top level

        if (isoValues_.size() > 1)
        {
            FatalIOErrorInFunction(dict)
                << "Multiple values on iso-surface (point) not supported"
                << " since needs original interpolators." << nl
                << exit(FatalIOError);
        }
    }

    if (isoParams_.algorithm() == isoSurfaceParams::ALGO_TOPO)
    {
        if
        (
            triangulate_
         && (isoParams_.filter() == isoSurfaceParams::filterType::NONE)
        )
        {
            FatalIOErrorInFunction(dict)
                << "Cannot triangulate without a regularise filter" << nl
                << exit(FatalIOError);
        }
    }


    // Zones

    if (!dict.readIfPresent("zones", zoneNames_) && dict.found("zone"))
    {
        zoneNames_.resize(1);
        dict.readEntry("zone", zoneNames_.first());
    }

    if (-1 != mesh.cellZones().findIndex(zoneNames_))
    {
        dict.readIfPresent("exposedPatchName", exposedPatchName_);

        DebugInfo
            << "Restricting to cellZone(s) " << flatOutput(zoneNames_)
            << " with exposed internal faces into patch "
            << mesh.boundaryMesh().findPatchID(exposedPatchName_) << endl;
    }
}


Foam::sampledIsoSurface::sampledIsoSurface
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    sampledIsoSurface(isoSurfaceParams(), name, mesh, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sampledIsoSurface::~sampledIsoSurface()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::sampledIsoSurface::needsUpdate() const
{
    const fvMesh& fvm = static_cast<const fvMesh&>(mesh());

    return fvm.time().timeIndex() != prevTimeIndex_;
}


bool Foam::sampledIsoSurface::expire()
{
    surface_.clear();
    meshCells_.clear();
    isoSurfacePtr_.reset(nullptr);
    subMeshPtr_.reset(nullptr);

    // Clear derived data
    sampledSurface::clearGeom();

    // Already marked as expired
    if (prevTimeIndex_ == -1)
    {
        return false;
    }

    // Force update
    prevTimeIndex_ = -1;
    return true;
}


bool Foam::sampledIsoSurface::update()
{
    return updateGeometry();
}


Foam::tmp<Foam::scalarField>
Foam::sampledIsoSurface::sample
(
    const interpolation<scalar>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::vectorField>
Foam::sampledIsoSurface::sample
(
    const interpolation<vector>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::sphericalTensorField>
Foam::sampledIsoSurface::sample
(
    const interpolation<sphericalTensor>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::symmTensorField>
Foam::sampledIsoSurface::sample
(
    const interpolation<symmTensor>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::tensorField>
Foam::sampledIsoSurface::sample
(
    const interpolation<tensor>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::scalarField>
Foam::sampledIsoSurface::interpolate
(
    const interpolation<scalar>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


Foam::tmp<Foam::vectorField>
Foam::sampledIsoSurface::interpolate
(
    const interpolation<vector>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}

Foam::tmp<Foam::sphericalTensorField>
Foam::sampledIsoSurface::interpolate
(
    const interpolation<sphericalTensor>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


Foam::tmp<Foam::symmTensorField>
Foam::sampledIsoSurface::interpolate
(
    const interpolation<symmTensor>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


Foam::tmp<Foam::tensorField>
Foam::sampledIsoSurface::interpolate
(
    const interpolation<tensor>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


void Foam::sampledIsoSurface::print(Ostream& os, int level) const
{
    os  << "isoSurface: " << name() << " :";
    isoParams_.print(os);
    os  << " field:" << isoField_
        << " value:" << flatOutput(isoValues_);
}


// ************************************************************************* //
