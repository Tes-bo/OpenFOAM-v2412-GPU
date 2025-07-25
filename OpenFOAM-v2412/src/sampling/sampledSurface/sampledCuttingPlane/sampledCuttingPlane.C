/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "sampledCuttingPlane.H"
#include "dictionary.H"
#include "fvMesh.H"
#include "volFields.H"
#include "volPointInterpolation.H"
#include "cartesianCS.H"
#include "addToRunTimeSelectionTable.H"
#include "PtrList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sampledCuttingPlane, 0);
    addNamedToRunTimeSelectionTable
    (
        sampledSurface,
        sampledCuttingPlane,
        word,
        cuttingPlane
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::plane Foam::sampledCuttingPlane::definePlane
(
    const polyMesh& mesh,
    const dictionary& dict
)
{
    plane pln(dict);

    // Optional (cartesian) coordinate transform.
    // - with registry to allow lookup from globally defined systems

    auto csysPtr = coordinateSystem::NewIfPresent(mesh, dict);

    if (!csysPtr)
    {
        csysPtr = coordinateSystem::NewIfPresent(dict, "transform");
    }

    // Make plane relative to the Cartesian coordinate system
    if (csysPtr)
    {
        coordSystem::cartesian cs(csysPtr());

        const point  orig = cs.globalPosition(pln.origin());
        const vector norm = cs.globalVector(pln.normal());

        DebugInfo
            << "plane "
            << " origin:" << pln.origin()
            << " normal:" << pln.normal()
            << " =>"
            << " origin:" << orig << " normal:" << norm
            << endl;

        // Reassign the plane
        pln = plane(orig, norm);
    }

    return pln;
}


void Foam::sampledCuttingPlane::checkBoundsIntersection
(
    const plane& pln,
    const boundBox& meshBb
) const
{
    // Verify specified bounding box
    const boundBox& clipBb = isoParams_.getClipBounds();

    if (clipBb.good())
    {
        // Bounding box does not overlap with (global) mesh!
        if (!clipBb.overlaps(meshBb))
        {
            WarningInFunction
                << nl
                << name() << " : "
                << "Bounds " << clipBb
                << " do not overlap the mesh bounding box " << meshBb
                << nl << endl;
        }

        // Plane does not intersect the bounding box
        if (!clipBb.intersects(pln))
        {
            WarningInFunction
                << nl
                << name() << " : "
                << "Plane "<< pln << " does not intersect the bounds "
                << clipBb
                << nl << endl;
        }
    }

    // Plane does not intersect the (global) mesh!
    if (!meshBb.intersects(pln))
    {
        WarningInFunction
            << nl
            << name() << " : "
            << "Plane "<< pln << " does not intersect the mesh bounds "
            << meshBb
            << nl << endl;
    }
}


void Foam::sampledCuttingPlane::setDistanceFields(const plane& pln)
{
    volScalarField& cellDistance = cellDistancePtr_();

    // Get mesh from volField,
    // so automatically submesh or baseMesh

    const fvMesh& mesh = cellDistance.mesh();

    // Distance to cell centres
    // ~~~~~~~~~~~~~~~~~~~~~~~~

    // Internal field
    {
        const auto& cc = mesh.cellCentres();
        scalarField& fld = cellDistance.primitiveFieldRef();

        forAll(cc, i)
        {
            fld[i] = pln.signedDistance(cc[i]);
        }
    }

    // Patch fields
    {
        volScalarField::Boundary& cellDistanceBf =
            cellDistance.boundaryFieldRef();

        forAll(cellDistanceBf, patchi)
        {
            if
            (
                isA<emptyFvPatchScalarField>
                (
                    cellDistanceBf[patchi]
                )
            )
            {
                cellDistanceBf.set
                (
                    patchi,
                    new calculatedFvPatchScalarField
                    (
                        mesh.boundary()[patchi],
                        cellDistance
                    )
                );

                const polyPatch& pp = mesh.boundary()[patchi].patch();
                pointField::subField cc = pp.patchSlice(mesh.faceCentres());

                fvPatchScalarField& fld = cellDistanceBf[patchi];
                fld.setSize(pp.size());
                forAll(fld, i)
                {
                    fld[i] = pln.signedDistance(cc[i]);
                }
            }
            else
            {
                // Other side cell centres?
                const pointField& cc = mesh.C().boundaryField()[patchi];
                fvPatchScalarField& fld = cellDistanceBf[patchi];

                forAll(fld, i)
                {
                    fld[i] = pln.signedDistance(cc[i]);
                }
            }
        }
    }


    // On processor patches the mesh.C() will already be the cell centre
    // on the opposite side so no need to swap cellDistance.

    // Distance to points
    pointDistance_.resize(mesh.nPoints());
    {
        const pointField& pts = mesh.points();

        forAll(pointDistance_, i)
        {
            pointDistance_[i] = pln.signedDistance(pts[i]);
        }
    }
}


void Foam::sampledCuttingPlane::combineSurfaces
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
        // Shift ownership from list to autoPtr
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


void Foam::sampledCuttingPlane::createGeometry()
{
    if (debug)
    {
        Pout<< "sampledCuttingPlane::createGeometry :updating geometry."
            << endl;
    }

    // Clear any previously stored topologies
    surface_.clear();
    meshCells_.clear();
    isoSurfacePtr_.reset(nullptr);

    // Clear derived data
    sampledSurface::clearGeom();

    // Clear any stored fields
    pointDistance_.clear();
    cellDistancePtr_.clear();

    const bool hasCellZones =
        (-1 != mesh().cellZones().findIndex(zoneNames_));

    const fvMesh& fvm = static_cast<const fvMesh&>(this->mesh());

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

            bitSet cellsToSelect(mesh().cellZones().selection(zoneNames_));

            DebugInfo
                << "Allocating subset of size "
                << cellsToSelect.count() << " with exposed faces into patch "
                << exposedPatchi << endl;


            // If we will use a fvMeshSubset so can apply bounds as well to make
            // the initial selection smaller.

            const boundBox& clipBb = isoParams_.getClipBounds();
            if (clipBb.good() && cellsToSelect.any())
            {
                const auto& cellCentres = fvm.C();

                for (const label celli : cellsToSelect)
                {
                    const point& cc = cellCentres[celli];

                    if (!clipBb.contains(cc))
                    {
                        cellsToSelect.unset(celli);
                    }
                }

                DebugInfo
                    << "Bounded subset of size "
                    << cellsToSelect.count() << endl;
            }

            subMeshPtr_.reset
            (
                new fvMeshSubset(fvm, cellsToSelect, exposedPatchi)
            );
        }
    }


    // Select either the submesh or the underlying mesh
    const fvMesh& mesh =
    (
        subMeshPtr_
      ? subMeshPtr_->subMesh()
      : fvm
    );

    checkBoundsIntersection(plane_, mesh.bounds());


    // Distance to cell centres
    // ~~~~~~~~~~~~~~~~~~~~~~~~

    cellDistancePtr_.reset
    (
        new volScalarField
        (
            mesh.newIOobject("cellDistance"),
            mesh,
            dimensionedScalar(dimLength, Zero)
        )
    );
    const auto& cellDistance = *cellDistancePtr_;

    setDistanceFields(plane_);

    if (debug)
    {
        Pout<< "Writing cell distance:" << cellDistance.objectPath() << endl;
        cellDistance.write();
        auto tpointDist = pointScalarField::New
        (
            "pointDistance",
            IOobject::NO_REGISTER,
            pointMesh::New(mesh),
            dimensionedScalar(dimLength, Zero)
        );
        auto& pointDist = tpointDist.ref();

        pointDist.primitiveFieldRef() = pointDistance_;

        Pout<< "Writing point distance:" << pointDist.objectPath() << endl;
        pointDist.write();
    }


    // Create surfaces for each offset

    PtrList<isoSurfaceBase> isoSurfPtrs(offsets_.size());

    forAll(offsets_, surfi)
    {
        isoSurfPtrs.set
        (
            surfi,
            isoSurfaceBase::New
            (
                isoParams_,
                cellDistance,
                pointDistance_,
                offsets_[surfi],
                *ignoreCellsPtr_
            )
        );
    }

    // And flatten
    combineSurfaces(isoSurfPtrs);


    // Discard fields if not required by an iso-surface
    if (!isoSurfacePtr_)
    {
        cellDistancePtr_.reset(nullptr);
        pointDistance_.clear();
    }

    if (debug)
    {
        print(Pout, debug);
        Pout<< endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledCuttingPlane::sampledCuttingPlane
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    sampledSurface(name, mesh, dict),
    plane_(definePlane(mesh, dict)),
    offsets_(),
    isoParams_
    (
        dict,
        isoSurfaceParams::ALGO_TOPO,
        isoSurfaceParams::filterType::DIAGCELL
    ),
    average_(dict.getOrDefault("average", false)),
    simpleSubMesh_(dict.getOrDefault("simpleSubMesh", false)),
    zoneNames_(),
    exposedPatchName_(),
    needsUpdate_(true),

    surface_(),
    meshCells_(),
    isoSurfacePtr_(nullptr),

    subMeshPtr_(nullptr),
    ignoreCellsPtr_(nullptr),
    cellDistancePtr_(nullptr),
    pointDistance_()
{
    dict.readIfPresent("offsets", offsets_);

    if (offsets_.empty())
    {
        offsets_.resize(1);
        offsets_.first() = Zero;
    }

    if (offsets_.size() > 1)
    {
        const label nOrig = offsets_.size();

        inplaceUniqueSort(offsets_);

        if (nOrig != offsets_.size())
        {
            IOWarningInFunction(dict)
                << "Removed non-unique offsets" << nl;
        }
    }

    if (isoParams_.algorithm() == isoSurfaceParams::ALGO_POINT)
    {
        // Not possible for ALGO_POINT
        simpleSubMesh_ = false;

        // Not possible for ALGO_POINT
        if (offsets_.size() > 1)
        {
            FatalIOErrorInFunction(dict)
                << "Multiple offsets with iso-surface (point) not supported"
                << " since needs original interpolators." << nl
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


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::sampledCuttingPlane::needsUpdate() const
{
    return needsUpdate_;
}


bool Foam::sampledCuttingPlane::expire()
{
    if (debug)
    {
        Pout<< "sampledCuttingPlane::expire :"
            << " needsUpdate:" << needsUpdate_ << endl;
    }

    surface_.clear();
    meshCells_.clear();
    isoSurfacePtr_.reset(nullptr);

    // Clear derived data
    sampledSurface::clearGeom();

    // Already marked as expired
    if (needsUpdate_)
    {
        return false;
    }

    needsUpdate_ = true;
    return true;
}


bool Foam::sampledCuttingPlane::update()
{
    if (debug)
    {
        Pout<< "sampledCuttingPlane::update :"
            << " needsUpdate:" << needsUpdate_ << endl;
    }

    if (!needsUpdate_)
    {
        return false;
    }

    createGeometry();

    needsUpdate_ = false;
    return true;
}


Foam::tmp<Foam::scalarField>
Foam::sampledCuttingPlane::sample
(
    const interpolation<scalar>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::vectorField>
Foam::sampledCuttingPlane::sample
(
    const interpolation<vector>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::sphericalTensorField>
Foam::sampledCuttingPlane::sample
(
    const interpolation<sphericalTensor>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::symmTensorField>
Foam::sampledCuttingPlane::sample
(
    const interpolation<symmTensor>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::tensorField>
Foam::sampledCuttingPlane::sample
(
    const interpolation<tensor>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::scalarField>
Foam::sampledCuttingPlane::interpolate
(
    const interpolation<scalar>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


Foam::tmp<Foam::vectorField>
Foam::sampledCuttingPlane::interpolate
(
    const interpolation<vector>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


Foam::tmp<Foam::sphericalTensorField>
Foam::sampledCuttingPlane::interpolate
(
    const interpolation<sphericalTensor>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


Foam::tmp<Foam::symmTensorField>
Foam::sampledCuttingPlane::interpolate
(
    const interpolation<symmTensor>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


Foam::tmp<Foam::tensorField>
Foam::sampledCuttingPlane::interpolate
(
    const interpolation<tensor>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


void Foam::sampledCuttingPlane::print(Ostream& os, int level) const
{
    os  << "sampledCuttingPlane: " << name() << " :"
        << " plane:" << plane_
        << " offsets:" << flatOutput(offsets_);

    if (level)
    {
        os  << "  faces:" << faces().size()
            << "  points:" << points().size();
    }
}


// ************************************************************************* //
