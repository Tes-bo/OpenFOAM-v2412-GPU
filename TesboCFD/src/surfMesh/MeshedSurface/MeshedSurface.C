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

#include "MeshedSurface.H"
#include "UnsortedMeshedSurface.H"
#include "MeshedSurfaceProxy.H"
#include "mergePoints.H"
#include "Time.H"
#include "ListOps.H"
#include "polyBoundaryMesh.H"
#include "polyMesh.H"
#include "surfMesh.H"
#include "primitivePatch.H"
#include "faceTraits.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class Face>
Foam::wordHashSet Foam::MeshedSurface<Face>::readTypes()
{
    return wordHashSet(*fileExtensionConstructorTablePtr_);
}


template<class Face>
Foam::wordHashSet Foam::MeshedSurface<Face>::writeTypes()
{
    return wordHashSet(*writefileExtensionMemberFunctionTablePtr_);
}


template<class Face>
bool Foam::MeshedSurface<Face>::canReadType
(
    const word& fileType,
    bool verbose
)
{
    return fileFormats::surfaceFormatsCore::checkSupport
    (
        readTypes() | FriendType::readTypes(),
        fileType,
        verbose,
        "reading"
    );
}


template<class Face>
bool Foam::MeshedSurface<Face>::canWriteType
(
    const word& fileType,
    bool verbose
)
{
    return fileFormats::surfaceFormatsCore::checkSupport
    (
        writeTypes() | ProxyType::writeTypes(),
        fileType,
        verbose,
        "writing"
    );
}


template<class Face>
bool Foam::MeshedSurface<Face>::canRead
(
    const fileName& name,
    bool verbose
)
{
    const word ext =
    (
        name.has_ext("gz")
      ? name.stem().ext()
      : name.ext()
    );

    return canReadType(ext, verbose);
}


template<class Face>
void Foam::MeshedSurface<Face>::write
(
    const fileName& name,
    const MeshedSurface<Face>& surf,
    IOstreamOption streamOpt,
    const dictionary& options
)
{
    write(name, name.ext(), surf, streamOpt, options);
}


template<class Face>
void Foam::MeshedSurface<Face>::write
(
    const fileName& name,
    const word& fileType,
    const MeshedSurface<Face>& surf,
    IOstreamOption streamOpt,
    const dictionary& options
)
{
    if (fileType.empty())
    {
        // Handle empty/missing type

        const word ext(name.ext());

        if (ext.empty())
        {
            FatalErrorInFunction
                << "Cannot determine format from filename" << nl
                << "    " << name << nl
                << exit(FatalError);
        }

        write(name, ext, surf, streamOpt, options);
        return;
    }


    DebugInFunction << "Writing to " << name << nl;

    auto* mfuncPtr = writefileExtensionMemberFunctionTable(fileType);

    if (!mfuncPtr)
    {
        // Delegate to proxy if possible
        const wordHashSet delegate(ProxyType::writeTypes());

        if (!delegate.found(fileType))
        {
            FatalErrorInFunction
                << "Unknown write format " << fileType << nl << nl
                << "Valid types:" << nl
                << flatOutput((delegate | writeTypes()).sortedToc()) << nl
                << exit(FatalError);
        }

        MeshedSurfaceProxy<Face>(surf).write
        (
            name, fileType, streamOpt, options
        );
    }
    else
    {
        mfuncPtr(name, surf, streamOpt, options);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface()
:
    MeshReference(List<Face>(), pointField()),
    faceIds_(),
    zones_()
{}


template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface
(
    const MeshedSurface<Face>& surf
)
:
    MeshReference(surf.surfFaces(), surf.points()),
    faceIds_(surf.faceIds_),
    zones_(surf.zones_)
{}


template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface
(
    const UnsortedMeshedSurface<Face>& surf
)
:
    MeshReference(List<Face>(), surf.points()),  // Copy points only
    faceIds_(),
    zones_()
{
    labelList faceMap;
    this->storedZones() = surf.sortedZones(faceMap);

    // Faces, in the sorted order
    const List<Face>& origFaces = surf;
    {
        List<Face> newFaces(origFaces.size());

        forAll(newFaces, facei)
        {
            newFaces[faceMap[facei]] = origFaces[facei];
        }

        this->storedFaces().transfer(newFaces);
    }

    // FaceIds, in the sorted order
    const labelList& origIds = surf.faceIds();

    if (origIds.size() == origFaces.size())
    {
        labelList newFaceIds(origIds.size());

        forAll(newFaceIds, facei)
        {
            newFaceIds[faceMap[facei]] = origIds[facei];
        }

        this->storedFaceIds().transfer(newFaceIds);
    }
}


template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface
(
    MeshedSurface<Face>&& surf
)
:
    MeshedSurface<Face>()
{
    transfer(surf);
}


template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface
(
    UnsortedMeshedSurface<Face>&& surf
)
:
    MeshedSurface<Face>()
{
    transfer(surf);
}


template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface
(
    const pointField& pointLst,
    const UList<Face>& faceLst,
    const UList<surfZone>& zoneLst
)
:
    MeshReference(faceLst, pointLst), // Copy construct
    faceIds_(),
    zones_(zoneLst)
{
    this->checkZones(false);  // Non-verbose fix zones
}


template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface
(
    pointField&& pointLst,
    List<Face>&& faceLst,
    const UList<surfZone>& zoneLst
)
:
    MeshReference(faceLst, pointLst, true), // Move construct
    faceIds_(),
    zones_(zoneLst)
{
    this->checkZones(false);  // Non-verbose fix zones
}


template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface
(
    const pointField& pointLst,
    const UList<Face>& faceLst,
    const labelUList& zoneSizes,
    const UList<word>& zoneNames
)
:
    MeshReference(faceLst, pointLst), // Copy construct
    faceIds_(),
    zones_()
{
    if (zoneSizes.size())
    {
        if (zoneNames.size())
        {
            addZones(zoneSizes, zoneNames);
        }
        else
        {
            addZones(zoneSizes);
        }
    }
}


template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface
(
    pointField&& pointLst,
    List<Face>&& faceLst,
    const labelUList& zoneSizes,
    const UList<word>& zoneNames
)
:
    MeshReference(faceLst, pointLst, true), // Move construct
    faceIds_(),
    zones_()
{
    if (zoneSizes.size())
    {
        if (zoneNames.size())
        {
            addZones(zoneSizes, zoneNames);
        }
        else
        {
            addZones(zoneSizes);
        }
    }
}


template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface(const surfMesh& mesh)
:
    MeshedSurface<Face>()
{
    // Need same face type as surfMesh
    MeshedSurface<face> surf
    (
        mesh.points(),
        mesh.faces(),
        mesh.surfZones()
    );

    this->transcribe(surf);
}


template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface
(
    const polyBoundaryMesh& bMesh,
    const bool useGlobalPoints
)
:
    MeshedSurface<Face>()
{
    const polyMesh& mesh = bMesh.mesh();
    const polyPatchList& bPatches = bMesh;

    // Get a single patch for all boundaries
    primitivePatch allBoundary
    (
        SubList<face>
        (
            mesh.faces(),
            mesh.nBoundaryFaces(),
            mesh.nInternalFaces()
        ),
        mesh.points()
    );

    // use global/local points:
    const pointField& bPoints =
    (
        useGlobalPoints ? mesh.points() : allBoundary.localPoints()
    );

    // global/local face addressing:
    const List<Face>& bFaces =
    (
        useGlobalPoints ? allBoundary : allBoundary.localFaces()
    );


    // create zone list
    surfZoneList newZones(bPatches.size());

    label startFacei = 0;
    label nZone = 0;
    for (const polyPatch& p : bPatches)
    {
        if (p.size())
        {
            newZones[nZone] = surfZone
            (
                p.name(),
                p.size(),
                startFacei,
                nZone
            );

            ++nZone;
            startFacei += p.size();
        }
    }

    newZones.setSize(nZone);

    // Face type as per polyBoundaryMesh
    MeshedSurface<face> surf(bPoints, bFaces, newZones);

    this->transcribe(surf);
}


template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface
(
    const fileName& name,
    const word& fileType
)
:
    MeshedSurface<Face>()
{
    read(name, fileType);
}


template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface(const fileName& name)
:
    MeshedSurface<Face>()
{
    read(name);
}


template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface(Istream& is)
:
    MeshedSurface<Face>()
{
    read(is);
}


template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface
(
    const Time& runTime
)
:
    MeshedSurface<Face>(runTime, word::null)
{}


template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface
(
    const Time& runTime,
    const word& surfName
)
:
    MeshedSurface<Face>()
{
    surfMesh mesh
    (
        IOobject
        (
            "dummyName",
            runTime.timeName(),
            runTime,
            IOobjectOption::MUST_READ,
            IOobjectOption::NO_WRITE,
            IOobjectOption::NO_REGISTER
        ),
        surfName
    );

    // The geometry components, returned via autoPtr
    MeshedSurface<face> surf
    (
        std::move(*(mesh.releaseGeom()))
    );

    this->transcribe(surf);
}


template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface
(
    const IOobject& io,
    const dictionary& dict,
    const bool isGlobal
)
:
    MeshedSurface<Face>()
{
    fileName fName
    (
        fileFormats::surfaceFormatsCore::checkFile(io, dict, isGlobal)
    );

    this->read(fName, dict.getOrDefault<word>("fileType", word::null));

    this->scalePoints(dict.getOrDefault<scalar>("scale", 0));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Face>
Foam::MeshedSurface<Face>::~MeshedSurface()
{
    clear();
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Face>
void Foam::MeshedSurface<Face>::remapFaces
(
    const labelUList& faceMapNewToOld
)
{
    if (faceMapNewToOld.empty())
    {
        return;
    }

    surfZoneList& zones = storedZones();

    if (zones.size() == 1)
    {
        // Single zone case is trivial
        zones[0].size() = faceMapNewToOld.size();
        return;
    }

    // Recalculate the zone start/size
    label newFacei = 0;
    label origEndi = 0;

    for (surfZone& zone : zones)
    {
        // Adjust zone start
        zone.start() = newFacei;
        origEndi += zone.size();

        for (label facei = newFacei; facei < faceMapNewToOld.size(); ++facei)
        {
            if (faceMapNewToOld[facei] < origEndi)
            {
                ++newFacei;
            }
            else
            {
                break;
            }
        }

        // Adjust zone size
        zone.size() = newFacei - zone.start();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Face>
void Foam::MeshedSurface<Face>::clear()
{
    MeshReference::clearOut();  // Topology changes

    storedPoints().clear();
    storedFaces().clear();
    storedFaceIds().clear();
    storedZones().clear();
}


template<class Face>
void Foam::MeshedSurface<Face>::movePoints(const pointField& newPoints)
{
    MeshReference::clearGeom();  // Changes areas, normals etc.

    // Adapt for new point positions
    MeshReference::movePoints(newPoints);

    // Copy new points
    storedPoints() = newPoints;
}


template<class Face>
void Foam::MeshedSurface<Face>::scalePoints(const scalar scaleFactor)
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


// Remove badly degenerate faces, double faces.
template<class Face>
void Foam::MeshedSurface<Face>::cleanup(const bool verbose)
{
    // Merge points (already done for STL, TRI)
    stitchFaces(SMALL, verbose);

    checkFaces(verbose);
    this->checkTopology(verbose);
}


template<class Face>
void Foam::MeshedSurface<Face>::compactPoints(labelList& pointMap)
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


template<class Face>
bool Foam::MeshedSurface<Face>::stitchFaces
(
    const scalar tol,
    const bool verbose
)
{
    pointField& ps = this->storedPoints();

    // Merge points (inplace)
    labelList pointMap;
    label nChanged = Foam::inplaceMergePoints(ps, tol, verbose, pointMap);

    if (!nChanged)
    {
        return false;
    }

    if (verbose)
    {
        InfoInFunction<< "Renumbering all faces" << endl;
    }

    List<Face>& faceLst = this->storedFaces();

    labelList faceMap(faceLst.size(), -1);

    // Reset the point labels to the unique points array
    label newFacei = 0;
    forAll(faceLst, facei)
    {
        Face& f = faceLst[facei];
        for (label& vert : f)
        {
            vert = pointMap[vert];
        }

        // For extra safety: collapse face as well
        if (f.collapse() >= 3)
        {
            if (newFacei != facei)
            {
                faceLst[newFacei] = f;
            }
            faceMap[newFacei] = facei;
            ++newFacei;
        }
        else if (verbose)
        {
            Pout<< "MeshedSurface::stitchFaces : "
                << "Removing collapsed face " << facei << endl
                << "    vertices   :" << f << endl;
        }
    }
    pointMap.clear();

    if (newFacei != faceLst.size())
    {
        if (verbose)
        {
            Pout<< "MeshedSurface::stitchFaces : "
                << "Removed " << faceLst.size() - newFacei
                << " faces" << endl;
        }
        faceMap.resize(newFacei);
        faceLst.resize(newFacei);

        // The faceMap is a newToOld mapping and only removes elements
        if (faceIds_.size())
        {
            forAll(faceMap, facei)
            {
                faceIds_[facei] = faceIds_[faceMap[facei]];
            }

            faceIds_.resize(newFacei);
        }

        remapFaces(faceMap);
    }
    faceMap.clear();

    // Topology can change when points are merged, etc
    MeshReference::clearOut();

    return true;
}


// Remove badly degenerate faces and double faces.
template<class Face>
bool Foam::MeshedSurface<Face>::checkFaces
(
    const bool verbose
)
{
    bool changed = false;
    List<Face>& faceLst = this->storedFaces();

    labelList faceMap(faceLst.size());

    label newFacei = 0;
    const label maxPointi = this->points().size();

    // Detect badly labelled faces and mark degenerate faces
    forAll(faceLst, facei)
    {
        Face& f = faceLst[facei];

        // Avoid degenerate faces
        if (f.collapse() >= 3)
        {
            for (const label vert : f)
            {
                if (vert < 0 || vert >= maxPointi)
                {
                    FatalErrorInFunction
                        << "face " << f
                        << " uses point indices outside point range 0.."
                        << (maxPointi-1)
                        << exit(FatalError);
                }
            }

            faceMap[facei] = facei;
            ++newFacei;
        }
        else
        {
            // Mark as bad face
            faceMap[facei] = -1;

            changed = true;
            if (verbose)
            {
                WarningInFunction
                    << "face[" << facei << "] = " << f
                    << " does not have three unique vertices" << endl;
            }
        }
    }

    // Detect doubled faces
    // do not touch the faces
    const labelListList& fFaces = this->faceFaces();
    newFacei = 0;
    forAll(faceLst, facei)
    {
        // Skip already collapsed faces
        if (faceMap[facei] < 0)
        {
            continue;
        }

        const Face& f = faceLst[facei];

        // Duplicate face check
        bool okay = true;
        const labelList& neighbours = fFaces[facei];

        // Check if faceNeighbours use same points as this face.
        // Note: discards normal information - sides of baffle are merged.
        for (const label neiFacei : neighbours)
        {
            if (neiFacei <= facei || faceMap[neiFacei] < 0)
            {
                // lower numbered faces already checked
                // skip neighbours that are themselves collapsed
                continue;
            }

            const Face& nei = faceLst[neiFacei];

            if (f == nei)
            {
                okay = false;

                if (verbose)
                {
                    WarningInFunction
                        << "faces share the same vertices:" << nl
                        << "    face[" << facei << "] : " << f << nl
                        << "    face[" << neiFacei << "] : " << nei << endl;
                    // printFace(Warning, "    ", f, points());
                    // printFace(Warning, "    ", nei, points());
                }

                break;
            }
        }

        if (okay)
        {
            faceMap[facei] = facei;
            ++newFacei;
        }
        else
        {
            faceMap[facei] = -1;
        }
    }


    // Until now, faceMap is an identity for good faces and -1 for bad faces

    // Phase 1: pack
    // Done to keep numbering constant in phase 1

    if (changed || newFacei < faceLst.size())
    {
        changed = true;

        if (verbose)
        {
            WarningInFunction
                << "Removed " << faceLst.size() - newFacei
                << " illegal faces." << endl;
        }

        // Compress the face list
        newFacei = 0;
        forAll(faceLst, facei)
        {
            if (faceMap[facei] >= 0)
            {
                if (newFacei != facei)
                {
                    faceLst[newFacei] = std::move(faceLst[facei]);
                }
                faceMap[newFacei] = facei;
                ++newFacei;
            }
        }

        faceMap.resize(newFacei);
        faceLst.resize(newFacei);

        // The faceMap is a newToOld mapping and only removes elements
        if (faceIds_.size())
        {
            forAll(faceMap, facei)
            {
                faceIds_[facei] = faceIds_[faceMap[facei]];
            }

            faceIds_.resize(newFacei);
        }

        remapFaces(faceMap);
    }
    faceMap.clear();

    // Topology can change because of renumbering
    MeshReference::clearOut();
    return changed;
}


template<class Face>
Foam::label Foam::MeshedSurface<Face>::nTriangles() const
{
    if (faceTraits<Face>::isTri())
    {
        return MeshReference::size();
    }

    return nTriangles
    (
        const_cast<labelList&>(labelList::null())
    );
}


template<class Face>
Foam::label Foam::MeshedSurface<Face>::nTriangles
(
    labelList& faceMap
) const
{
    label nTri = 0;
    const List<Face>& faceLst = surfFaces();

    // Count triangles needed
    for (const auto& f : faceLst)
    {
        nTri += f.nTriangles();
    }

    // Nothing to do
    if (nTri <= faceLst.size())
    {
        if (notNull(faceMap))
        {
            faceMap.clear();
        }
    }
    else if (notNull(faceMap))
    {
        // Face map requested
        faceMap.resize(nTri);

        nTri = 0;
        forAll(faceLst, facei)
        {
            label n = faceLst[facei].nTriangles();
            while (n-- > 0)
            {
                faceMap[nTri++] = facei;
            }
        }

        faceMap.resize(nTri);
    }

    return nTri;
}


template<class Face>
Foam::label Foam::MeshedSurface<Face>::triangulate()
{
    if (faceTraits<Face>::isTri())
    {
        // Inplace triangulation of triFace/labelledTri surface = no-op
        return 0;
    }
    else
    {
        return triangulate
        (
            const_cast<labelList&>(labelList::null())
        );
    }
}


template<class Face>
Foam::label Foam::MeshedSurface<Face>::triangulate
(
    labelList& faceMapOut
)
{
    labelList dummyFaceMap;

    labelList& faceMap =
    (
        notNull(faceMapOut)
      ? faceMapOut
      : dummyFaceMap
    );

    if (faceTraits<Face>::isTri())
    {
        // Inplace triangulation of triFace/labelledTri surface = no-op
        faceMap.clear();
        return 0;
    }

    label nTri = 0;
    label maxTri = 0;  // the maximum number of triangles for any single face
    List<Face>& faceLst = this->storedFaces();

    // How many triangles will be needed
    for (const auto& f : faceLst)
    {
        const label n = f.nTriangles();
        if (maxTri < n)
        {
            maxTri = n;
        }
        nTri += n;
    }

    // Nothing to do
    if (nTri <= faceLst.size())
    {
        faceMap.clear();
        return 0;
    }

    this->storedFaceIds().clear();  // Invalid or misleading

    List<Face> newFaces(nTri);
    faceMap.resize(nTri);

    if (this->points().empty())
    {
        // triangulate without points
        // simple face triangulation around f[0]
        nTri = 0;
        forAll(faceLst, facei)
        {
            const Face& f = faceLst[facei];

            for (label fp = 1; fp < f.size() - 1; ++fp)
            {
                const label fp1 = f.fcIndex(fp);

                newFaces[nTri] = Face{f[0], f[fp], f[fp1]};
                faceMap[nTri] = facei;
                ++nTri;
            }
        }
    }
    else
    {
        // triangulate with points
        List<face> tmpTri(maxTri);

        nTri = 0;
        forAll(faceLst, facei)
        {
            // 'face' not '<Face>'
            const face& f = faceLst[facei];

            label nTmp = 0;
            f.triangles(this->points(), nTmp, tmpTri);
            for (label triI = 0; triI < nTmp; triI++)
            {
                newFaces[nTri] = Face
                (
                    static_cast<labelUList&>(tmpTri[triI])
                );
                faceMap[nTri] = facei;
                ++nTri;
            }
        }
    }

    // The number of *additional* faces
    nTri -= faceLst.size();

    faceLst.transfer(newFaces);
    remapFaces(faceMap);

    // Topology can change because of renumbering
    MeshReference::clearOut();

    return nTri;
}


template<class Face>
Foam::MeshedSurface<Face>
Foam::MeshedSurface<Face>::subsetMeshImpl
(
    const labelList& pointMap,
    const labelList& faceMap
) const
{
    const pointField& locPoints = this->localPoints();
    const List<Face>& locFaces  = this->localFaces();

    // Subset of points (compact)
    pointField newPoints(UIndirectList<point>(locPoints, pointMap));

    // Inverse point mapping - same as ListOps invert() without checks
    labelList oldToNew(locPoints.size(), -1);
    forAll(pointMap, pointi)
    {
        oldToNew[pointMap[pointi]] = pointi;
    }

    // Subset of faces
    List<Face> newFaces(UIndirectList<Face>(locFaces, faceMap));

    // Renumber face node labels
    for (auto& f : newFaces)
    {
        for (label& vert : f)
        {
            vert = oldToNew[vert];
        }
    }
    oldToNew.clear();

    // Deep copy of zones, leave start/size intact!!
    surfZoneList newZones(zones_);

    // Recalculate the zone start/size
    label newFacei = 0;
    label origEndi = 0;

    for (surfZone& zone : newZones)
    {
        // The old zone ending
        origEndi += zone.size();

        // The new zone start
        zone.start() = newFacei;

        for (label facei = newFacei; facei < faceMap.size(); ++facei)
        {
            if (faceMap[facei] < origEndi)
            {
                ++newFacei;
            }
            else
            {
                break;
            }
        }

        // The new zone size
        zone.size() = newFacei - zone.start();
    }


    // Subset of faceIds. Can be empty.
    labelList newFaceIds;
    if (faceIds_.size())
    {
        newFaceIds = labelUIndList(faceIds_, faceMap);
    }

    // Construct the sub-surface
    MeshedSurface<Face> newSurf;
    newSurf.storedFaces().transfer(newFaces);
    newSurf.storedPoints().transfer(newPoints);
    newSurf.storedZones().transfer(newZones);
    newSurf.storedFaceIds().transfer(newFaceIds);

    return newSurf;
}


template<class Face>
Foam::MeshedSurface<Face>
Foam::MeshedSurface<Face>::subsetMesh
(
    const UList<bool>& include,
    labelList& pointMap,
    labelList& faceMap
) const
{
    this->subsetMeshMap(include, pointMap, faceMap);
    return this->subsetMeshImpl(pointMap, faceMap);
}


template<class Face>
Foam::MeshedSurface<Face>
Foam::MeshedSurface<Face>::subsetMesh
(
    const bitSet& include,
    labelList& pointMap,
    labelList& faceMap
) const
{
    this->subsetMeshMap(include, pointMap, faceMap);
    return this->subsetMeshImpl(pointMap, faceMap);
}


template<class Face>
Foam::MeshedSurface<Face>
Foam::MeshedSurface<Face>::subsetMesh
(
    const UList<bool>& include
) const
{
    labelList pointMap, faceMap;
    return this->subsetMesh(include, pointMap, faceMap);
}


template<class Face>
Foam::MeshedSurface<Face>
Foam::MeshedSurface<Face>::subsetMesh
(
    const bitSet& include
) const
{
    labelList pointMap, faceMap;
    return this->subsetMesh(include, pointMap, faceMap);
}


template<class Face>
Foam::MeshedSurface<Face> Foam::MeshedSurface<Face>::subsetMesh
(
    const wordRes& includeNames,
    const wordRes& excludeNames
) const
{
    bitSet include(this->size());

    for
    (
        const label zonei
      : fileFormats::surfaceFormatsCore::getSelectedPatches
        (
            zones_,
            includeNames,
            excludeNames
        )
    )
    {
        include.set(zones_[zonei].range());
    }

    return this->subsetMesh(include);
}


template<class Face>
void Foam::MeshedSurface<Face>::swap
(
    MeshedSurface<Face>& surf
)
{
    if (this == &surf)
    {
        return;  // Self-swap is a no-op
    }

    MeshReference::clearOut(); // Topology changes
    surf.clearOut();        // Topology changes

    this->storedPoints().swap(surf.storedPoints());
    this->storedFaces().swap(surf.storedFaces());
    this->storedZones().swap(surf.storedZones());
    this->storedFaceIds().swap(surf.storedFaceIds());
}


template<class Face>
void Foam::MeshedSurface<Face>::transfer
(
    pointField& pointLst,
    List<Face>& faceLst
)
{
    MeshReference::clearOut();  // Topology changes

    this->storedPoints().transfer(pointLst);
    this->storedFaces().transfer(faceLst);
    this->storedZones().clear();
    this->storedFaceIds().clear();  // Likely to be invalid
}


template<class Face>
void Foam::MeshedSurface<Face>::transfer
(
    MeshedSurface<Face>& surf
)
{
    if (this == &surf)
    {
        return;  // Self-assigment is a no-op
    }

    MeshReference::clearOut();  // Topology changes

    this->storedPoints().transfer(surf.storedPoints());
    this->storedFaces().transfer(surf.storedFaces());
    this->storedZones().transfer(surf.storedZones());
    this->storedFaceIds().transfer(surf.storedFaceIds());

    surf.clear();
}


template<class Face>
void Foam::MeshedSurface<Face>::transfer
(
    UnsortedMeshedSurface<Face>& surf
)
{
    // Clear everything
    this->clear();

    labelList faceMap;
    surfZoneList zoneLst = surf.sortedZones(faceMap);

    List<Face>& faceLst = surf.storedFaces();

    if (zoneLst.size() > 1)
    {
        // Unknown if we really need to sort the faces
        List<Face> sortedFaces(faceMap.size());

        forAll(faceMap, facei)
        {
            sortedFaces[faceMap[facei]].transfer(faceLst[facei]);
        }

        faceLst.swap(sortedFaces);  // Replace with sorted faces
    }

    MeshedSurface<Face> newSurf
    (
        std::move(surf.storedPoints()),
        std::move(faceLst),
        zoneLst
    );

    surf.clear();

    this->swap(newSurf);
}


template<class Face>
Foam::autoPtr<Foam::MeshedSurface<Face>>
Foam::MeshedSurface<Face>::releaseGeom()
{
    return autoPtr<MeshedSurface<Face>>::New(std::move(*this));
}


template<class Face>
void Foam::MeshedSurface<Face>::swapFaces(List<Face>& faces)
{
    MeshReference::clearOut();  // Topology changes

    this->storedFaceIds().clear();  // Likely to be invalid

    this->storedFaces().swap(faces);

    this->checkZones(false);  // Non-verbose fix zones
}


template<class Face>
void Foam::MeshedSurface<Face>::swapPoints(pointField& points)
{
    // Adapt for new point positions
    MeshReference::movePoints(points);

    this->storedPoints().swap(points);
}


template<class Face>
bool Foam::MeshedSurface<Face>::read(const fileName& name)
{
    this->clear();
    transfer(*New(name));
    return true;
}


template<class Face>
bool Foam::MeshedSurface<Face>::read
(
    const fileName& name,
    const word& fileType
)
{
    this->clear();
    transfer(*New(name, fileType));
    return true;
}


template<class Face>
void Foam::MeshedSurface<Face>::write
(
    const Time& t,
    const word& surfName
) const
{
    MeshedSurfaceProxy<Face>(*this).write(t, surfName);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Face>
void Foam::MeshedSurface<Face>::operator=(const MeshedSurface<Face>& surf)
{
    if (this == &surf)
    {
        return;  // Self-assignment is a no-op
    }

    // Clear everything
    this->clear();

    this->storedPoints() = surf.points();
    this->storedFaces()  = surf.surfFaces();
    this->storedFaceIds() = surf.faceIds();
    this->storedZones()  = surf.surfZones();
}


template<class Face>
void Foam::MeshedSurface<Face>::operator=(MeshedSurface<Face>&& surf)
{
    transfer(surf);
}


template<class Face>
Foam::MeshedSurface<Face>::operator Foam::MeshedSurfaceProxy<Face>() const
{
    return MeshedSurfaceProxy<Face>
    (
        this->points(),
        this->surfFaces(),
        this->surfZones(),
        labelUList::null(), // faceMap = none
        this->faceIds()
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "MeshedSurfaceZones.C"
#include "MeshedSurfaceIO.C"
#include "MeshedSurfaceNew.C"

// ************************************************************************* //
