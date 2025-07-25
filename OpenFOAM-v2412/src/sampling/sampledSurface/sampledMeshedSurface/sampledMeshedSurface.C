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

#include "sampledMeshedSurface.H"
#include "meshSearch.H"
#include "Tuple2.H"
#include "globalIndex.H"
#include "treeDataCell.H"
#include "treeDataFace.H"
#include "meshTools.H"
#include "addToRunTimeSelectionTable.H"
#include "stringListOps.H"   // For stringListOps::findMatching()

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::sampledMeshedSurface::samplingSource
>
Foam::sampledMeshedSurface::samplingSourceNames_
({
    { samplingSource::cells, "cells" },
    { samplingSource::insideCells, "insideCells" },
    { samplingSource::boundaryFaces, "boundaryFaces" },
});


namespace Foam
{
    defineTypeNameAndDebug(sampledMeshedSurface, 0);
    // Use shorter name only
    addNamedToRunTimeSelectionTable
    (
        sampledSurface,
        sampledMeshedSurface,
        word,
        meshedSurface
    );
    // Compatibility name (1912)
    addNamedToRunTimeSelectionTable
    (
        sampledSurface,
        sampledMeshedSurface,
        word,
        sampledTriSurfaceMesh
    );

} // End namespace Foam


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// The IOobject for reading
inline static IOobject selectReadIO(const word& name, const Time& runTime)
{
    return IOobject
    (
        name,
        runTime.constant(),     // instance
        "triSurface",           // local
        runTime,                // registry
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        IOobject::NO_REGISTER
    );
}

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sampledMeshedSurface::setZoneMap()
{
    // Populate zoneIds_ based on surfZone information

    const meshedSurface& s = static_cast<const meshedSurface&>(*this);

    const auto& zones = s.surfZones();

    zoneIds_.resize(s.size());

    // Trivial case
    if (zoneIds_.empty() || zones.size() <= 1)
    {
        zoneIds_ = 0;
        return;
    }


    label beg = 0;

    forAll(zones, zonei)
    {
        const label len = min(zones[zonei].size(), zoneIds_.size() - beg);

        // Assign sub-zone Ids
        SubList<label>(zoneIds_, len, beg) = zonei;

        beg += len;
    }

    // Anything remaining? Should not happen
    {
        const label len = (zoneIds_.size() - beg);

        if (len > 0)
        {
            SubList<label>(zoneIds_, len, beg) = max(0, zones.size()-1);
        }
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::sampledMeshedSurface::update(const meshSearch& meshSearcher)
{
    // Global numbering for cells/faces
    // - only used to uniquely identify local elements
    globalIndex globalCells(onBoundary() ? mesh().nFaces() : mesh().nCells());

    // Find the cells the triangles of the surface are in.
    // Does approximation by looking at the face centres only
    const pointField& fc = surface_.faceCentres();

    // sqr(distance), global index
    typedef Tuple2<scalar, label> nearInfo;
    List<nearInfo> nearest(fc.size(), nearInfo(Foam::sqr(GREAT), labelMax));

    if (sampleSource_ == samplingSource::cells)
    {
        // Search for nearest cell

        const auto& cellTree = meshSearcher.cellTree();
        const auto& treeData = cellTree.shapes();

        forAll(fc, facei)
        {
            const point& pt = fc[facei];
            auto& near = nearest[facei];

            pointIndexHit info = cellTree.findNearest(pt, Foam::sqr(GREAT));

            if (info.hit())
            {
                const label objectIndex = treeData.objectIndex(info.index());

                near.first()  = info.point().distSqr(pt);
                near.second() = globalCells.toGlobal(objectIndex);
            }
        }
    }
    else if (sampleSource_ == samplingSource::insideCells)
    {
        // Search for cell containing point

        const auto& cellTree = meshSearcher.cellTree();
        const auto& treeData = cellTree.shapes();

        forAll(fc, facei)
        {
            const point& pt = fc[facei];
            auto& near = nearest[facei];

            if (cellTree.bb().contains(pt))
            {
                const label index = cellTree.findInside(pt);

                if (index != -1)
                {
                    const label objectIndex = treeData.objectIndex(index);

                    near.first()  = 0;
                    near.second() = globalCells.toGlobal(objectIndex);
                }
            }
        }
    }
    else  // samplingSource::boundaryFaces
    {
        // Search for nearest boundary face
        // on all non-coupled boundary faces

        const auto& bndTree = meshSearcher.nonCoupledBoundaryTree();
        const auto& treeData = bndTree.shapes();

        forAll(fc, facei)
        {
            const point& pt = fc[facei];
            auto& near = nearest[facei];

            pointIndexHit info = bndTree.findNearest(pt, sqr(GREAT));

            if (info.hit())
            {
                const label objectIndex = treeData.objectIndex(info.index());

                near.first()  = info.point().distSqr(pt);
                near.second() = globalCells.toGlobal(objectIndex);
            }
        }
    }


    // See which processor has the nearest. Mark and subset
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Pstream::listCombineReduce(nearest, minFirstEqOp<scalar>{});

    labelList cellOrFaceLabels(fc.size(), -1);

    bitSet facesToSubset(fc.size());

    forAll(nearest, facei)
    {
        const auto& near = nearest[facei];

        const label index = near.second();

        if (index == labelMax)
        {
            // Not found on any processor, or simply too far.
            // How to map?
        }
        else if (globalCells.isLocal(index))
        {
            facesToSubset.set(facei);

            // Mark as special if too far away
            cellOrFaceLabels[facei] =
            (
                (near.first() < maxDistanceSqr_)
              ? globalCells.toLocal(index)
              : -1
            );
        }
    }


    if (debug)
    {
        Pout<< "Local out of faces:" << cellOrFaceLabels.size()
            << " keeping:" << facesToSubset.count() << endl;
    }


    // Subset the surface
    meshedSurface& s = static_cast<meshedSurface&>(*this);

    labelList pointMap;
    labelList faceMap;

    s = surface_.subsetMesh(facesToSubset, pointMap, faceMap);


    // Ensure zoneIds_ are indeed correct
    setZoneMap();


    // Subset cellOrFaceLabels (for compact faces)
    cellOrFaceLabels = labelUIndList(cellOrFaceLabels, faceMap)();


    // Collect the samplePoints and sampleElements
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (sampledSurface::isPointData())
    {
        // With point interpolation

        samplePoints_ = points();
        sampleElements_.resize(pointMap.size(), -1);

        // Store any face per point (without using pointFaces())
        labelList pointToFace(std::move(pointMap));

        forAll(s, facei)
        {
            const face& f = s[facei];

            for (const label labi : f)
            {
                pointToFace[labi] = facei;
            }
        }


        if (sampleSource_ == samplingSource::cells)
        {
            // sampleElements_ : per surface point, the cell
            // samplePoints_   : per surface point, a location inside the cell

            forAll(samplePoints_, pointi)
            {
                // Use _copy_ of point
                const point pt = samplePoints_[pointi];

                const label celli = cellOrFaceLabels[pointToFace[pointi]];

                sampleElements_[pointi] = celli;

                if
                (
                    celli >= 0
                 && !mesh().pointInCell(pt, celli, meshSearcher.decompMode())
                )
                {
                    // Point not inside cell
                    // Find nearest point on faces of cell

                    scalar minDistSqr = VGREAT;

                    for (const label facei : mesh().cells()[celli])
                    {
                        const face& f = mesh().faces()[facei];

                        pointHit info =
                            f.nearestPoint
                            (
                                pt,
                                mesh().points()
                            );

                        if (info.distance() < minDistSqr)
                        {
                            minDistSqr = info.distance();
                            samplePoints_[pointi] = info.point();
                        }
                    }
                }
            }
        }
        else if (sampleSource_ == samplingSource::insideCells)
        {
            // sampleElements_ : per surface point the cell
            // samplePoints_   : per surface point a location inside the cell

            forAll(samplePoints_, pointi)
            {
                const label celli = cellOrFaceLabels[pointToFace[pointi]];

                sampleElements_[pointi] = celli;
            }
        }
        else  // samplingSource::boundaryFaces
        {
            // sampleElements_ : per surface point, the boundary face containing
            //                   the location
            // samplePoints_   : per surface point, a location on the boundary

            forAll(samplePoints_, pointi)
            {
                const point& pt = samplePoints_[pointi];

                const label facei = cellOrFaceLabels[pointToFace[pointi]];

                sampleElements_[pointi] = facei;

                if (facei >= 0)
                {
                    samplePoints_[pointi] =
                        mesh().faces()[facei].nearestPoint
                        (
                            pt,
                            mesh().points()
                        ).point();
                }
            }
        }
    }
    else
    {
        // if sampleSource_ == cells:
        //      sampleElements_ : per surface face, the cell
        //      samplePoints_   : n/a
        // if sampleSource_ == insideCells:
        //      sampleElements_ : -1 or per surface face, the cell
        //      samplePoints_   : n/a
        // if sampleSource_ == boundaryFaces:
        //      sampleElements_ : per surface face, the boundary face
        //      samplePoints_   : n/a

        sampleElements_.transfer(cellOrFaceLabels);
        samplePoints_.clear();
    }


    if (debug)
    {
        this->clearOut();
        OFstream str(mesh().time().path()/"surfaceToMesh.obj");
        Info<< "Dumping correspondence from local surface (points or faces)"
            << " to mesh (cells or faces) to " << str.name() << endl;

        const vectorField& centres =
        (
            onBoundary()
          ? mesh().faceCentres()
          : mesh().cellCentres()
        );

        if (sampledSurface::isPointData())
        {
            label vertI = 0;
            forAll(samplePoints_, pointi)
            {
                meshTools::writeOBJ(str, points()[pointi]);
                ++vertI;

                meshTools::writeOBJ(str, samplePoints_[pointi]);
                ++vertI;

                label elemi = sampleElements_[pointi];
                meshTools::writeOBJ(str, centres[elemi]);
                ++vertI;
                str << "l " << vertI-2 << ' ' << vertI-1 << ' ' << vertI << nl;
            }
        }
        else
        {
            label vertI = 0;
            forAll(sampleElements_, triI)
            {
                meshTools::writeOBJ(str, faceCentres()[triI]);
                ++vertI;

                label elemi = sampleElements_[triI];
                meshTools::writeOBJ(str, centres[elemi]);
                ++vertI;
                str << "l " << vertI-1 << ' ' << vertI << nl;
            }
        }
    }

    needsUpdate_ = false;
    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledMeshedSurface::sampledMeshedSurface
(
    const word& name,
    const polyMesh& mesh,
    const word& surfaceName,
    const samplingSource sampleSource
)
:
    sampledSurface(name, mesh),
    meshedSurface(),
    surfaceName_(surfaceName),
    surface_
    (
        selectReadIO(surfaceName, mesh.time()),
        dictionary::null
    ),
    sampleSource_(sampleSource),
    needsUpdate_(true),
    keepIds_(true),
    zoneIds_(),
    sampleElements_(),
    samplePoints_(),
    maxDistanceSqr_(Foam::sqr(GREAT)),
    defaultValues_()
{}


Foam::sampledMeshedSurface::sampledMeshedSurface
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    sampledSurface(name, mesh, dict),
    meshedSurface(),
    surfaceName_
    (
        meshedSurface::findFile
        (
            selectReadIO(dict.get<word>("surface"), mesh.time()),
            dict
        ).name()
    ),
    surface_
    (
        selectReadIO(dict.get<word>("surface"), mesh.time()),
        dict
    ),
    sampleSource_(samplingSourceNames_.get("source", dict)),
    needsUpdate_(true),
    keepIds_(dict.getOrDefault("keepIds", true)),
    zoneIds_(),
    sampleElements_(),
    samplePoints_(),
    maxDistanceSqr_(Foam::sqr(GREAT)),
    defaultValues_(dict.subOrEmptyDict("defaultValue"))
{
    if (dict.readIfPresent("maxDistance", maxDistanceSqr_))
    {
        // Info<< "Limit to maxDistance " << maxDistanceSqr_ << nl;
        // if (defaultValues_.empty())
        // {
        //     Info<< "defaultValues = zero" << nl;
        // }
        // else
        // {
        //     defaultValues_.writeEntry(Info);
        // }

        maxDistanceSqr_ = Foam::sqr(maxDistanceSqr_);
    }

    wordRes includePatches;
    dict.readIfPresent("patches", includePatches);
    includePatches.uniq();

    // Could also shift this to the reader itself,
    // but not yet necessary.

    if (!includePatches.empty())
    {
        Info<< "Subsetting surface " << surfaceName_
            << " to patches: " << flatOutput(includePatches) << nl;

        const surfZoneList& zones = surface_.surfZones();

        const labelList zoneIndices
        (
            stringListOps::findMatching
            (
                zones,
                includePatches,
                wordRes(),
                nameOp<surfZone>()
            )
        );

        // Faces to subset
        bitSet includeMap(surface_.size());

        for (const label zonei : zoneIndices)
        {
            const surfZone& zn = zones[zonei];
            includeMap.set(zn.range());
        }

        if (includeMap.none())
        {
            WarningInFunction
                << "Patch selection results in an empty surface"
                << " - ignoring" << nl;
        }
        else if (!includeMap.all())
        {
            meshedSurface subSurf(surface_.subsetMesh(includeMap));

            if (subSurf.empty())
            {
                WarningInFunction
                    << "Bad surface subset (empty)"
                    << " - skip and hope for the best" << nl;
            }
            else
            {
                // Replace
                surface_.transfer(subSurf);
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::sampledMeshedSurface::needsUpdate() const
{
    return needsUpdate_;
}


bool Foam::sampledMeshedSurface::expire()
{
    // already marked as expired
    if (needsUpdate_)
    {
        return false;
    }

    sampledSurface::clearGeom();
    meshedSurface::clear();
    zoneIds_.clear();

    sampleElements_.clear();
    samplePoints_.clear();

    needsUpdate_ = true;
    return true;
}


bool Foam::sampledMeshedSurface::update()
{
    if (!needsUpdate_)
    {
        return false;
    }

    // Calculate surface and mesh overlap bounding box
    treeBoundBox bb(surface_.points(), surface_.meshPoints());

    // Restrict surface to (global!) mesh bound box
    bb &= mesh().bounds();

    if (bb.good())
    {
        // Extend a bit
        bb.grow(0.5*bb.span());
        bb.inflate(1e-6);
    }
    else
    {
        // Surface and mesh do not overlap at all. Guarantee a valid
        // bounding box so we don't get any 'invalid bounding box' errors.

        WarningInFunction
            << "Surface " << surfaceName_
            << " does not overlap bounding box of mesh " << mesh().bounds()
            << endl;

        bb.reset(mesh().bounds().centre());
        bb.grow(1e-6*mesh().bounds().span());
    }

    // Mesh search engine, no triangulation of faces.
    meshSearch meshSearcher(mesh(), bb, polyMesh::FACE_PLANES);

    return update(meshSearcher);
}


bool Foam::sampledMeshedSurface::update(const treeBoundBox& bb)
{
    if (!needsUpdate_)
    {
        return false;
    }

    // Mesh search engine on subset, no triangulation of faces.
    meshSearch meshSearcher(mesh(), bb, polyMesh::FACE_PLANES);

    return update(meshSearcher);
}


Foam::tmp<Foam::scalarField> Foam::sampledMeshedSurface::sample
(
    const interpolation<scalar>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::vectorField> Foam::sampledMeshedSurface::sample
(
    const interpolation<vector>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::sphericalTensorField> Foam::sampledMeshedSurface::sample
(
    const interpolation<sphericalTensor>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::symmTensorField> Foam::sampledMeshedSurface::sample
(
    const interpolation<symmTensor>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::tensorField> Foam::sampledMeshedSurface::sample
(
    const interpolation<tensor>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::scalarField> Foam::sampledMeshedSurface::interpolate
(
    const interpolation<scalar>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


Foam::tmp<Foam::vectorField> Foam::sampledMeshedSurface::interpolate
(
    const interpolation<vector>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}

Foam::tmp<Foam::sphericalTensorField> Foam::sampledMeshedSurface::interpolate
(
    const interpolation<sphericalTensor>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


Foam::tmp<Foam::symmTensorField> Foam::sampledMeshedSurface::interpolate
(
    const interpolation<symmTensor>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


Foam::tmp<Foam::tensorField> Foam::sampledMeshedSurface::interpolate
(
    const interpolation<tensor>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


void Foam::sampledMeshedSurface::print(Ostream& os, int level) const
{
    os  << "meshedSurface: " << name() << " :"
        << " surface:" << surfaceName_;

    if (level)
    {
        os  << " faces:"   << faces().size()
            << " points:"  << points().size()
            << " zoneids:" << zoneIds().size();
    }
}


// ************************************************************************* //
