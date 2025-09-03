/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2024 OpenCFD Ltd.
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

#include "triSurface.H"
#include "Time.H"
#include "surfZoneList.H"
#include "MeshedSurface.H"
#include "ListOps.H"
#include "stringListOps.H"   // For stringListOps::findMatching()

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(triSurface, 0);
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Helper function to print triangle info
static void printTriangle
(
    Ostream& os,
    const string& pre,
    const labelledTri& f,
    const pointField& points
)
{
    os
        << pre.c_str() << "vertex numbers:"
        << f[0] << ' ' << f[1] << ' ' << f[2] << nl

        << pre.c_str() << "vertex coords :"
        << points[f[0]] << ' ' << points[f[1]] << ' ' << points[f[2]]

        << pre.c_str() << "region        :" << f.region() << nl
        << endl;
}

} // End namespace Foam


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::fileName Foam::triSurface::triSurfInstance(const Time& d)
{
    fileName foamName(d.caseName() + ".ftr");

    // Search back through the time directories list to find the time
    // closest to and lower than current time

    instantList ts = d.times();
    label i;

    for (i=ts.size()-1; i>=0; i--)
    {
        if (ts[i].value() <= d.timeOutputValue())
        {
            break;
        }
    }

    // Noting that the current directory has already been searched
    // for mesh data, start searching from the previously stored time directory

    if (i>=0)
    {
        for (label j=i; j>=0; j--)
        {
            if (isFile(d.path()/ts[j].name()/typeName/foamName))
            {
                if (debug)
                {
                    Pout<< " triSurface::triSurfInstance(const Time& d)"
                        << "reading " << foamName
                        << " from " << ts[j].name()/typeName
                        << endl;
                }

                return ts[j].name();
            }
        }
    }

    if (debug)
    {
        Pout<< " triSurface::triSurfInstance(const Time& d)"
            << "reading " << foamName
            << " from constant/" << endl;
    }
    return d.constant();
}


Foam::List<Foam::labelledTri> Foam::triSurface::convertToTri
(
    const faceList& faces,
    const label defaultRegion
)
{
    List<labelledTri> triFaces(faces.size());

    forAll(triFaces, facei)
    {
        const face& f = faces[facei];

        if (f.size() != 3)
        {
            FatalErrorInFunction
                << "Face at position " << facei
                << " does not have three vertices:" << f
                << abort(FatalError);
        }

        labelledTri& tri = triFaces[facei];

        tri[0] = f[0];
        tri[1] = f[1];
        tri[2] = f[2];
        tri.region() = defaultRegion;
    }

    return triFaces;
}


Foam::List<Foam::labelledTri> Foam::triSurface::convertToTri
(
    const triFaceList& faces,
    const label defaultRegion
)
{
    List<labelledTri> triFaces(faces.size());

    forAll(triFaces, facei)
    {
        const triFace& f = faces[facei];

        labelledTri& tri = triFaces[facei];

        tri[0] = f[0];
        tri[1] = f[1];
        tri[2] = f[2];
        tri.region() = defaultRegion;
    }

    return triFaces;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Remove non-triangles, double triangles.
void Foam::triSurface::checkTriangles(const bool verbose)
{
    // Simple check on indices ok.
    const label maxPointi = points().size() - 1;

    for (const auto& f : *this)
    {
        for (const label verti : f)
        {
            if (verti < 0 || verti > maxPointi)
            {
                FatalErrorInFunction
                    << "triangle " << f
                    << " uses point indices outside point range 0.."
                    << maxPointi
                    << exit(FatalError);
            }
        }
    }

    // Two phase process
    //   1. mark invalid faces
    //   2. pack
    // Done to keep numbering constant in phase 1

    // List of valid triangles
    bitSet valid(size(), true);

    forAll(*this, facei)
    {
        const labelledTri& f = (*this)[facei];

        if (!f.good())
        {
            // 'degenerate' triangle check
            valid.unset(facei);

            if (verbose)
            {
                WarningInFunction
                    << "triangle " << facei
                    << " does not have three unique vertices:\n";
                printTriangle(Warning, "    ", f, points());
            }
        }
        else
        {
            // duplicate triangle check
            const labelList& fEdges = faceEdges()[facei];

            // Check if faceNeighbours use same points as this face.
            // Note: discards normal information - sides of baffle are merged.

            for (const label edgei : fEdges)
            {
                const labelList& eFaces = edgeFaces()[edgei];

                for (const label neighbour : eFaces)
                {
                    if (neighbour > facei)
                    {
                        // lower numbered faces already checked
                        const labelledTri& n = (*this)[neighbour];

                        if
                        (
                            ((f[0] == n[0]) || (f[0] == n[1]) || (f[0] == n[2]))
                         && ((f[1] == n[0]) || (f[1] == n[1]) || (f[1] == n[2]))
                         && ((f[2] == n[0]) || (f[2] == n[1]) || (f[2] == n[2]))
                        )
                        {
                            valid.unset(facei);

                            if (verbose)
                            {
                                WarningInFunction
                                    << "triangles share the same vertices:\n"
                                    << "    face 1 :" << facei << endl;
                                printTriangle(Warning, "    ", f, points());

                                Warning
                                    << endl
                                    << "    face 2 :"
                                    << neighbour << endl;
                                printTriangle(Warning, "    ", n, points());
                            }

                            break;
                        }
                    }
                }
            }
        }
    }

    if (!valid.all())
    {
        // Compact
        label newFacei = 0;
        for (const label facei : valid)
        {
            (*this)[newFacei++] = (*this)[facei];
        }

        if (verbose)
        {
            WarningInFunction
                << "Removing " << size() - newFacei
                << " illegal faces." << endl;
        }
        (*this).setSize(newFacei);

        // Topology can change because of renumbering
        clearOut();
    }
}


// Check/fix edges with more than two triangles
void Foam::triSurface::checkEdges(const bool verbose)
{
    const labelListList& eFaces = edgeFaces();

    forAll(eFaces, edgei)
    {
        const labelList& myFaces = eFaces[edgei];

        if (myFaces.empty())
        {
            FatalErrorInFunction
                << "Edge " << edgei << " with vertices " << edges()[edgei]
                << " has no edgeFaces"
                << exit(FatalError);
        }
        else if (myFaces.size() > 2 && verbose)
        {
            WarningInFunction
                << "Edge " << edgei << " with vertices " << edges()[edgei]
                << " has more than 2 faces connected to it : " << myFaces
                << endl;
        }
    }
}


// Returns patch info. Sets faceMap to the indexing according to patch
// numbers. Patch numbers start at 0.
Foam::surfacePatchList
Foam::triSurface::calcPatches(labelList& faceMap) const
{
    // Determine the sorted order:
    // use sortedOrder directly (the intermediate list is discarded anyhow)

    labelList regions(size());
    forAll(regions, facei)
    {
        regions[facei] = operator[](facei).region();
    }
    sortedOrder(regions, faceMap);
    regions.clear();

    // Extend regions
    label maxRegion = patches_.size()-1;    // for non-compacted regions

    if (!faceMap.empty())
    {
        maxRegion = max
        (
            maxRegion,
            operator[](faceMap.back()).region()
        );
    }

    // Get new region list
    surfacePatchList newPatches(maxRegion + 1);

    // Fill patch sizes
    forAll(*this, facei)
    {
        label region = operator[](facei).region();

        newPatches[region].size()++;
    }


    // Fill rest of patch info

    label startFacei = 0;
    forAll(newPatches, newPatchi)
    {
        surfacePatch& newPatch = newPatches[newPatchi];

        newPatch.index() = newPatchi;
        newPatch.start() = startFacei;

        // Take over any information from existing patches
        if
        (
            newPatchi < patches_.size()
         && !patches_[newPatchi].name().empty()
        )
        {
            newPatch.name() = patches_[newPatchi].name();
        }
        else
        {
            newPatch.name() = surfacePatch::defaultName(newPatchi);
        }

        if
        (
            newPatchi < patches_.size()
         && !patches_[newPatchi].geometricType().empty()
        )
        {
            newPatch.geometricType() = patches_[newPatchi].geometricType();
        }
        else
        {
            newPatch.geometricType() = surfacePatch::emptyType;
        }

        startFacei += newPatch.size();
    }

    return newPatches;
}


void Foam::triSurface::setDefaultPatches()
{
    labelList faceMap;

    // Get names, types and sizes
    surfacePatchList newPatches(calcPatches(faceMap));

    // Take over names and types (but not size)
    patches_.setSize(newPatches.size());

    forAll(newPatches, patchi)
    {
        patches_[patchi].index() = patchi;
        patches_[patchi].name() = newPatches[patchi].name();
        patches_[patchi].geometricType() = newPatches[patchi].geometricType();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::triSurface::triSurface()
:
    MeshReference(List<labelledTri>(), pointField()),
    patches_(),
    sortedEdgeFacesPtr_(nullptr),
    edgeOwnerPtr_(nullptr)
{}


Foam::triSurface::triSurface(const triSurface& surf)
:
    MeshReference(surf, surf.points()),
    patches_(surf.patches()),
    sortedEdgeFacesPtr_(nullptr),
    edgeOwnerPtr_(nullptr)
{}


Foam::triSurface::triSurface(triSurface&& surf)
:
    triSurface()
{
    transfer(surf);
}


Foam::triSurface::triSurface
(
    const List<labelledTri>& triangles,
    const geometricSurfacePatchList& patches,
    const pointField& pts
)
:
    MeshReference(triangles, pts),
    patches_(patches),
    sortedEdgeFacesPtr_(nullptr),
    edgeOwnerPtr_(nullptr)
{}


Foam::triSurface::triSurface
(
    List<labelledTri>& triangles,
    const geometricSurfacePatchList& patches,
    pointField& pts,
    const bool reuse
)
:
    MeshReference(triangles, pts, reuse),
    patches_(patches),
    sortedEdgeFacesPtr_(nullptr),
    edgeOwnerPtr_(nullptr)
{}


Foam::triSurface::triSurface
(
    const List<labelledTri>& triangles,
    const pointField& pts
)
:
    MeshReference(triangles, pts),
    patches_(),
    sortedEdgeFacesPtr_(nullptr),
    edgeOwnerPtr_(nullptr)
{
    setDefaultPatches();
}


Foam::triSurface::triSurface
(
    const triFaceList& triangles,
    const pointField& pts
)
:
    MeshReference(convertToTri(triangles, 0), pts),
    patches_(),
    sortedEdgeFacesPtr_(nullptr),
    edgeOwnerPtr_(nullptr)
{
    setDefaultPatches();
}


Foam::triSurface::triSurface
(
    const fileName& name,
    const scalar scaleFactor
)
:
    triSurface(name, name.ext(), scaleFactor)
{}


Foam::triSurface::triSurface
(
    const fileName& name,
    const word& fileType,
    const scalar scaleFactor
)
:
    triSurface()
{
    read(name, fileType);
    scalePoints(scaleFactor);
    setDefaultPatches();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::triSurface::~triSurface()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::triSurface::clearTopology()
{
    MeshReference::clearTopology();
    sortedEdgeFacesPtr_.reset(nullptr);
    edgeOwnerPtr_.reset(nullptr);
}


void Foam::triSurface::clearPatchMeshAddr()
{
    MeshReference::clearPatchMeshAddr();
}


void Foam::triSurface::clearOut()
{
    MeshReference::clearOut();
    clearTopology();
    clearPatchMeshAddr();
}


void Foam::triSurface::swap(triSurface& surf)
{
    if (this == &surf)
    {
        return;  // Self-swap is a no-op
    }

    clearOut();
    surf.clearOut();

    storedFaces().swap(surf.storedFaces());
    storedPoints().swap(surf.storedPoints());
    patches_.swap(surf.patches());
}


const Foam::labelListList& Foam::triSurface::sortedEdgeFaces() const
{
    if (!sortedEdgeFacesPtr_)
    {
        calcSortedEdgeFaces();
    }

    return *sortedEdgeFacesPtr_;
}


const Foam::labelList& Foam::triSurface::edgeOwner() const
{
    if (!edgeOwnerPtr_)
    {
        calcEdgeOwner();
    }

    return *edgeOwnerPtr_;
}


void Foam::triSurface::movePoints(const pointField& pts)
{
    // Remove all geometry dependent data
    sortedEdgeFacesPtr_.reset(nullptr);

    // Adapt for new point positions
    MeshReference::movePoints(pts);

    // Copy new points
    storedPoints() = pts;
}


void Foam::triSurface::swapPoints(pointField& pts)
{
    // Remove all geometry dependent data
    sortedEdgeFacesPtr_.reset(nullptr);

    // Adapt for new point positions
    MeshReference::movePoints(pts);

    // Move/swap new points
    storedPoints().swap(pts);
}


void Foam::triSurface::scalePoints(const scalar scaleFactor)
{
    // Avoid bad or no scaling
    if (scaleFactor > SMALL && !equal(scaleFactor, 1))
    {
        // Remove all geometry dependent data
        this->clearTopology();

        // Adapt for new point positions
        MeshReference::movePoints(pointField());

        this->storedPoints() *= scaleFactor;
    }
}


// Remove non-triangles, double triangles.
void Foam::triSurface::cleanup(const bool verbose)
{
    // Merge points (already done for STL, TRI)
    stitchTriangles(SMALL, verbose);

    // Merging points might have changed geometric factors
    clearOut();

    checkTriangles(verbose);

    checkEdges(verbose);
}


void Foam::triSurface::compactPoints(labelList& pointMap)
{
    this->clearOut();   // Topology changes

    // Remove unused points while walking and renumbering faces
    // in visit order - walk order as per localFaces()

    labelList oldToCompact(this->points().size(), -1);
    DynamicList<label> compactPointMap(oldToCompact.size());

    for (auto& f : this->storedFaces())
    {
        for (label& pointi : f)
        {
            label compacti = oldToCompact[pointi];
            if (compacti == -1)
            {
                compacti = compactPointMap.size();
                oldToCompact[pointi] = compacti;
                compactPointMap.append(pointi);
            }
            pointi = compacti;
        }
    }

    pointField newPoints
    (
        UIndirectList<point>(this->points(), compactPointMap)
    );

    this->swapPoints(newPoints);

    if (notNull(pointMap))
    {
        pointMap.transfer(compactPointMap);
    }
}


Foam::List<Foam::surfZone>
Foam::triSurface::sortedZones(labelList& faceMap) const
{
    surfacePatchList patches(calcPatches(faceMap));

    surfZoneList zones(patches.size());
    forAll(patches, patchi)
    {
        zones[patchi] = surfZone(patches[patchi]);
    }

    return zones;
}


void Foam::triSurface::triFaceFaces(List<face>& plainFaces) const
{
    plainFaces.setSize(size());

    forAll(*this, facei)
    {
        plainFaces[facei] = this->operator[](facei);
    }
}


// Finds area, starting at facei, delimited by borderEdge. Marks all visited
// faces (from face-edge-face walk) with currentZone.
void Foam::triSurface::markZone
(
    const boolList& borderEdge,
    const label facei,
    const label currentZone,
    labelList& faceZone
) const
{
    // List of faces whose faceZone has been set.
    labelList changedFaces(1, facei);

    while (true)
    {
        // Pick up neighbours of changedFaces
        DynamicList<label> newChangedFaces(2*changedFaces.size());

        for (const label facei : changedFaces)
        {
            const labelList& fEdges = faceEdges()[facei];

            for (const label edgei : fEdges)
            {
                if (!borderEdge[edgei])
                {
                    const labelList& eFaces = edgeFaces()[edgei];

                    for (const label nbrFacei : eFaces)
                    {
                        if (faceZone[nbrFacei] == -1)
                        {
                            faceZone[nbrFacei] = currentZone;
                            newChangedFaces.append(nbrFacei);
                        }
                        else if (faceZone[nbrFacei] != currentZone)
                        {
                            FatalErrorInFunction
                                << "Zones " << faceZone[nbrFacei]
                                << " at face " << nbrFacei
                                << " connects to zone " << currentZone
                                << " at face " << facei
                                << abort(FatalError);
                        }
                    }
                }
            }
        }

        if (newChangedFaces.empty())
        {
            break;
        }

        changedFaces.transfer(newChangedFaces);
    }
}


// Finds areas delimited by borderEdge (or 'real' edges).
// Fills faceZone accordingly
Foam::label Foam::triSurface::markZones
(
    const boolList& borderEdge,
    labelList& faceZone
) const
{
    faceZone.setSize(size());
    faceZone = -1;

    if (borderEdge.size() != nEdges())
    {
        FatalErrorInFunction
            << "borderEdge boolList not same size as number of edges" << endl
            << "borderEdge:" << borderEdge.size() << endl
            << "nEdges    :" << nEdges()
            << exit(FatalError);
    }

    label zoneI = 0;

    for (label startFacei = 0;; zoneI++)
    {
        // Find first non-coloured face
        for (; startFacei < size(); startFacei++)
        {
            if (faceZone[startFacei] == -1)
            {
                break;
            }
        }

        if (startFacei >= size())
        {
            break;
        }

        faceZone[startFacei] = zoneI;

        markZone(borderEdge, startFacei, zoneI, faceZone);
    }

    return zoneI;
}


Foam::triSurface Foam::triSurface::subsetMeshImpl
(
    const labelList& pointMap,
    const labelList& faceMap
) const
{
    const pointField& locPoints = localPoints();
    const List<labelledTri>& locFaces = localFaces();

    // Subset of points (compact)
    pointField newPoints(UIndirectList<point>(locPoints, pointMap));

    // Inverse point mapping - same as ListOps invert() without checks
    labelList oldToNew(locPoints.size(), -1);
    forAll(pointMap, pointi)
    {
        oldToNew[pointMap[pointi]] = pointi;
    }

    // Subset of faces
    List<labelledTri> newFaces(UIndirectList<labelledTri>(locFaces, faceMap));

    // Renumber face node labels
    for (auto& f : newFaces)
    {
        for (label& vert : f)
        {
            vert = oldToNew[vert];
        }
    }
    oldToNew.clear();

    // Construct sub-surface
    return triSurface(newFaces, patches(), newPoints, true);
}


Foam::triSurface
Foam::triSurface::subsetMesh
(
    const UList<bool>& include,
    labelList& pointMap,
    labelList& faceMap
) const
{
    this->subsetMeshMap(include, pointMap, faceMap);
    return this->subsetMeshImpl(pointMap, faceMap);
}


Foam::triSurface
Foam::triSurface::subsetMesh
(
    const bitSet& include,
    labelList& pointMap,
    labelList& faceMap
) const
{
    this->subsetMeshMap(include, pointMap, faceMap);
    return this->subsetMeshImpl(pointMap, faceMap);
}


Foam::triSurface
Foam::triSurface::subsetMesh(const UList<bool>& include) const
{
    labelList pointMap, faceMap;
    return this->subsetMesh(include, pointMap, faceMap);
}


Foam::triSurface
Foam::triSurface::subsetMesh(const bitSet& include) const
{
    labelList pointMap, faceMap;
    return this->subsetMesh(include, pointMap, faceMap);
}


Foam::triSurface
Foam::triSurface::subsetMesh
(
    const wordRes& includeNames,
    const wordRes& excludeNames
) const
{
    const bitSet selectPatches
    (
        stringListOps::findMatching
        (
            patches_,
            includeNames,
            excludeNames,
            nameOp<geometricSurfacePatch>()
        )
    );

    bitSet include(this->size());

    forAll(*this, facei)
    {
        const label patchi = (*this)[facei].region();

        if (selectPatches.test(patchi))
        {
            include.set(facei);
        }
    }

    return this->subsetMesh(include);
}


void Foam::triSurface::swapFaces(List<labelledTri>& faceLst)
{
    clearOut();   // Topology changes

    this->storedFaces().swap(faceLst);
}


void Foam::triSurface::transfer(triSurface& surf)
{
    clearOut();

    storedFaces().transfer(surf.storedFaces());
    storedPoints().transfer(surf.storedPoints());
    patches_.transfer(surf.patches());

    surf.clearOut();
}


void Foam::triSurface::transfer(MeshedSurface<labelledTri>& surf)
{
    // Transcribe zone -> patch info
    auto patches = ListOps::create<geometricSurfacePatch>
    (
        surf.surfZones(),
        geometricSurfacePatch::fromIdentifier()
    );

    // Fairly ugly, but the only way to get at the data safely
    List<labelledTri> fcs;
    pointField pts;

    surf.swapFaces(fcs);
    surf.swapPoints(pts);

    clearOut();
    surf.clear();

    triSurface s(fcs, patches, pts, true);
    swap(s);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::triSurface::operator=(const triSurface& surf)
{
    clearOut();

    storedFaces() = surf;
    storedPoints() = surf.points();
    patches_ = surf.patches();
}


void Foam::triSurface::operator=(triSurface&& surf)
{
    transfer(surf);
}


void Foam::triSurface::operator=(MeshedSurface<labelledTri>&& surf)
{
    transfer(surf);
}


// ************************************************************************* //
