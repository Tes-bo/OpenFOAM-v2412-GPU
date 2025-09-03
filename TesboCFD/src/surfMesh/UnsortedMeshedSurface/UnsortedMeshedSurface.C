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
#include "Fstream.H"
#include "Time.H"
#include "ListOps.H"
#include "polyBoundaryMesh.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class Face>
Foam::wordHashSet Foam::UnsortedMeshedSurface<Face>::readTypes()
{
    return wordHashSet(*fileExtensionConstructorTablePtr_);
}


template<class Face>
Foam::wordHashSet Foam::UnsortedMeshedSurface<Face>::writeTypes()
{
    return wordHashSet(*writefileExtensionMemberFunctionTablePtr_);
}


template<class Face>
bool Foam::UnsortedMeshedSurface<Face>::canReadType
(
    const word& fileType,
    bool verbose
)
{
   return fileFormats::surfaceFormatsCore::checkSupport
   (
       readTypes() | MeshReference::readTypes(),
       fileType,
       verbose,
       "reading"
   );
}


template<class Face>
bool Foam::UnsortedMeshedSurface<Face>::canWriteType
(
    const word& fileType,
    bool verbose
)
{
    return fileFormats::surfaceFormatsCore::checkSupport
    (
        writeTypes(),
        fileType,
        verbose,
        "writing"
    );
}


template<class Face>
bool Foam::UnsortedMeshedSurface<Face>::canRead
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
void Foam::UnsortedMeshedSurface<Face>::write
(
    const fileName& name,
    const UnsortedMeshedSurface<Face>& surf,
    IOstreamOption streamOpt,
    const dictionary& options
)
{
    write(name, name.ext(), surf, streamOpt, options);
}


template<class Face>
void Foam::UnsortedMeshedSurface<Face>::write
(
    const fileName& name,
    const word& fileType,
    const UnsortedMeshedSurface<Face>& surf,
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
Foam::UnsortedMeshedSurface<Face>::UnsortedMeshedSurface()
:
    MeshReference()
{}


template<class Face>
Foam::UnsortedMeshedSurface<Face>::UnsortedMeshedSurface
(
    const UnsortedMeshedSurface<Face>& surf
)
:
    MeshReference(surf.points(), surf.surfFaces()), // Copy construct (no zones)
    zoneIds_(surf.zoneIds()),
    zoneToc_(surf.zoneToc())
{}


template<class Face>
Foam::UnsortedMeshedSurface<Face>::UnsortedMeshedSurface
(
    const MeshedSurface<Face>& surf
)
:
    MeshReference(surf.points(), surf.surfFaces()), // Copy construct (no zones)
    zoneIds_(),
    zoneToc_()
{
    setZones(surf.surfZones());
}


template<class Face>
Foam::UnsortedMeshedSurface<Face>::UnsortedMeshedSurface
(
    UnsortedMeshedSurface<Face>&& surf
)
:
    UnsortedMeshedSurface<Face>()
{
    transfer(surf);
}


template<class Face>
Foam::UnsortedMeshedSurface<Face>::UnsortedMeshedSurface
(
    MeshedSurface<Face>&& surf
)
:
    UnsortedMeshedSurface<Face>()
{
    transfer(surf);
}


template<class Face>
Foam::UnsortedMeshedSurface<Face>::UnsortedMeshedSurface
(
    pointField&& pointLst,
    List<Face>&& faceLst,
    List<label>&& zoneIds,
    UList<surfZoneIdentifier>& tocInfo
)
:
    MeshReference(std::move(pointLst), std::move(faceLst)),
    zoneIds_(std::move(zoneIds)),
    zoneToc_(tocInfo)
{}


template<class Face>
Foam::UnsortedMeshedSurface<Face>::UnsortedMeshedSurface
(
    const fileName& name,
    const word& ext
)
:
    UnsortedMeshedSurface<Face>()
{
    read(name, ext);
}


template<class Face>
Foam::UnsortedMeshedSurface<Face>::UnsortedMeshedSurface
(
    const fileName& name
)
:
    UnsortedMeshedSurface<Face>()
{
    read(name);
}


template<class Face>
Foam::UnsortedMeshedSurface<Face>::UnsortedMeshedSurface
(
    Istream& is
)
:
    UnsortedMeshedSurface<Face>()
{
    readIstream(is);
}


template<class Face>
Foam::UnsortedMeshedSurface<Face>::UnsortedMeshedSurface
(
    const Time& runTime
)
:
    UnsortedMeshedSurface<Face>()
{
    MeshedSurface<Face> surf(runTime);
    transfer(surf);
}


template<class Face>
Foam::UnsortedMeshedSurface<Face>::UnsortedMeshedSurface
(
    const Time& runTime,
    const word& surfName
)
:
    UnsortedMeshedSurface<Face>()
{
    MeshedSurface<Face> surf(runTime, surfName);
    transfer(surf);
}


template<class Face>
Foam::UnsortedMeshedSurface<Face>::UnsortedMeshedSurface
(
    const IOobject& io,
    const dictionary& dict,
    const bool isGlobal
)
:
    UnsortedMeshedSurface<Face>()
{
    fileName fName
    (
        fileFormats::surfaceFormatsCore::checkFile(io, dict, isGlobal)
    );

    this->read(fName, dict.getOrDefault<word>("fileType", word::null));

    this->scalePoints(dict.getOrDefault<scalar>("scale", 0));
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Face>
void Foam::UnsortedMeshedSurface<Face>::setOneZone()
{
    this->removeZones(); // Parent information is unreliable

    zoneIds_.resize(size());
    zoneIds_ = 0;

    // Assign single default zone
    zoneToc_.resize(1);

    zoneToc_[0].index() = 0;

    if (zoneToc_[0].name().empty())
    {
        zoneToc_[0].name() = "zone0";
    }
}


template<class Face>
void Foam::UnsortedMeshedSurface<Face>::setZones
(
    const surfZoneList& zoneLst
)
{
    this->removeZones(); // Parent information is unreliable

    zoneIds_.resize(size());
    zoneToc_.resize(zoneLst.size());

    forAll(zoneToc_, zonei)
    {
        const surfZone& zone = zoneLst[zonei];
        zoneToc_[zonei] = zone;

        // Assign sub-zone Ids
        SubList<label>(zoneIds_, zone.range()) = zonei;
    }
}


template<class Face>
void Foam::UnsortedMeshedSurface<Face>::setZones
(
    const labelUList& sizes,
    const UList<word>& names
)
{
    this->removeZones(); // Parent information is unreliable

    zoneIds_.resize(size());
    zoneToc_.resize(sizes.size());

    label start = 0;
    forAll(zoneToc_, zonei)
    {
        zoneToc_[zonei] = surfZoneIdentifier(names[zonei], zonei);

        // Assign sub-zone Ids
        SubList<label>(zoneIds_, sizes[zonei], start) = zonei;

        start += sizes[zonei];
    }
}


template<class Face>
void Foam::UnsortedMeshedSurface<Face>::setZones
(
    const labelUList& sizes
)
{
    this->removeZones(); // Parent information is unreliable

    zoneIds_.resize(size());
    zoneToc_.resize(sizes.size());

    label start = 0;
    forAll(zoneToc_, zonei)
    {
        zoneToc_[zonei] = surfZoneIdentifier
        (
            surfZoneIdentifier::defaultName(zonei),
            zonei
        );

        // Assign sub-zone Ids
        SubList<label>(zoneIds_, sizes[zonei], start) = zonei;

        start += sizes[zonei];
    }
}


template<class Face>
void Foam::UnsortedMeshedSurface<Face>::remapFaces
(
    const labelUList& faceMapNewToOld
)
{
    // Re-assign the zone Ids
    if (faceMapNewToOld.empty())
    {
        return;
    }

    if (zoneToc_.empty())
    {
        setOneZone();
    }
    else if (zoneToc_.size() == 1)
    {
        zoneIds_ = 0;  // Optimized for single-zone case
    }
    else
    {
        List<label> newZonesIds(faceMapNewToOld.size());

        forAll(faceMapNewToOld, facei)
        {
            newZonesIds[facei] = zoneIds_[faceMapNewToOld[facei]];
        }
        zoneIds_.transfer(newZonesIds);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Face>
bool Foam::UnsortedMeshedSurface<Face>::readIstream(Istream& is)
{
    is  >> this->storedZoneIds()
        >> this->storedPoints()
        >> this->storedFaces();

    is.check(FUNCTION_NAME);
    return true;
}


template<class Face>
void Foam::UnsortedMeshedSurface<Face>::writeOstream(Ostream& os) const
{
    os  << this->zoneIds()
        << this->points()
        << this->surfFaces();

    os.check(FUNCTION_NAME);
}


template<class Face>
void Foam::UnsortedMeshedSurface<Face>::setSize(const label s)
{
    this->storedFaces().resize(s);
    // if zones extend: set with last zoneId
    zoneIds_.resize(s, zoneToc_.size() - 1);
}


template<class Face>
void Foam::UnsortedMeshedSurface<Face>::clear()
{
    MeshReference::clear();
    zoneIds_.clear();
    zoneToc_.clear();
}


template<class Face>
Foam::surfZoneList Foam::UnsortedMeshedSurface<Face>::sortedZones
(
    labelList& faceMap
) const
{
    // supply some zone names
    Map<word> zoneNames;
    forAll(zoneToc_, zonei)
    {
        zoneNames.insert(zonei, zoneToc_[zonei].name());
    }

    // std::sort() really seems to mix up the order.
    // and std::stable_sort() might take too long / too much memory

    // Assuming that we have relatively fewer zones compared to the
    // number of items, just do it ourselves

    // Step 1: get zone sizes and store (origId => zoneI)
    Map<label> lookup;
    for (const label origId : zoneIds_)
    {
        ++(lookup(origId, 0));
    }

    // Step 2: assign start/size (and name) to the newZones
    // re-use the lookup to map (zoneId => zoneI)
    surfZoneList zoneLst(lookup.size());
    label start = 0;
    label zonei = 0;
    forAllIters(lookup, iter)
    {
        const label origId = iter.key();

        const word zoneName =
            zoneNames.lookup
            (
                origId,
                surfZoneIdentifier::defaultName(zonei)
            );

        zoneLst[zonei] = surfZone
        (
            zoneName,
            0,           // initialize with zero size
            start,
            zonei
        );

        // increment the start for the next zone
        // and save the (zoneId => zoneI) mapping
        start += iter();
        iter() = zonei++;
    }


    // Step 3: build the re-ordering
    faceMap.resize(zoneIds_.size());

    forAll(zoneIds_, facei)
    {
        const label zonei = lookup[zoneIds_[facei]];
        faceMap[facei] = zoneLst[zonei].start() + zoneLst[zonei].size()++;
    }

    // With reordered faces registered in faceMap
    return zoneLst;
}


template<class Face>
Foam::UnsortedMeshedSurface<Face>
Foam::UnsortedMeshedSurface<Face>::subsetMeshImpl
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

    // Subset of zones
    List<label> newZones(UIndirectList<label>(zoneIds_, faceMap));

    // Retain the same zone toc information
    List<surfZoneIdentifier> subToc(zoneToc_);

    // Construct the sub-surface
    return UnsortedMeshedSurface<Face>
    (
        std::move(newPoints),
        std::move(newFaces),
        std::move(newZones),
        std::move(subToc)
    );
}


template<class Face>
Foam::UnsortedMeshedSurface<Face>
Foam::UnsortedMeshedSurface<Face>::subsetMesh
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
Foam::UnsortedMeshedSurface<Face>
Foam::UnsortedMeshedSurface<Face>::subsetMesh
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
Foam::UnsortedMeshedSurface<Face>
Foam::UnsortedMeshedSurface<Face>::subsetMesh
(
    const UList<bool>& include
) const
{
    labelList pointMap, faceMap;
    return this->subsetMesh(include, pointMap, faceMap);
}


template<class Face>
Foam::UnsortedMeshedSurface<Face>
Foam::UnsortedMeshedSurface<Face>::subsetMesh
(
    const bitSet& include
) const
{
    labelList pointMap, faceMap;
    return this->subsetMesh(include, pointMap, faceMap);
}


template<class Face>
void Foam::UnsortedMeshedSurface<Face>::swap
(
    UnsortedMeshedSurface<Face>& surf
)
{
    if (this == &surf)
    {
        return;  // Self-swap is a no-op
    }

    this->clearOut();  // Topology changes
    surf.clearOut();   // Topology changes

    this->storedPoints().swap(surf.storedPoints());
    this->storedFaces().swap(surf.storedFaces());
    zoneIds_.swap(surf.zoneIds_);
    zoneToc_.swap(surf.zoneToc_);

    this->storedZones().clear(); // Should not be there anyhow
    surf.storedZones().clear();
}


template<class Face>
void Foam::UnsortedMeshedSurface<Face>::transfer
(
    UnsortedMeshedSurface<Face>& surf
)
{
    if (this == &surf)
    {
        return;  // Self-assignment is a no-op
    }

    this->clear();

    this->storedPoints().transfer(surf.storedPoints());
    this->storedFaces().transfer(surf.storedFaces());
    zoneIds_.transfer(surf.zoneIds_);
    zoneToc_.transfer(surf.zoneToc_);

    surf.clear();
}


template<class Face>
void Foam::UnsortedMeshedSurface<Face>::transfer
(
    MeshedSurface<Face>& surf
)
{
    surfZoneList zoneInfo(surf.surfZones());

    this->clear();

    MeshReference::transfer(surf);

    setZones(zoneInfo);
}


template<class Face>
Foam::autoPtr<Foam::labelList>
Foam::UnsortedMeshedSurface<Face>::releaseZoneIds()
{
    return autoPtr<labelList>::New(this->storedZoneIds());
}


template<class Face>
bool Foam::UnsortedMeshedSurface<Face>::read(const fileName& name)
{
    this->clear();
    transfer(*New(name));
    return true;
}


template<class Face>
bool Foam::UnsortedMeshedSurface<Face>::read
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
void Foam::UnsortedMeshedSurface<Face>::write
(
    const Time& t,
    const word& surfName
) const
{
    MeshedSurfaceProxy<Face>(*this).write(t, surfName);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Face>
void Foam::UnsortedMeshedSurface<Face>::operator=
(
    const UnsortedMeshedSurface<Face>& surf
)
{
    if (&surf == this)
    {
        return;  // Self-assignment is a no-op
    }

    clear();

    this->storedPoints() = surf.points();
    this->storedFaces()  = surf.surfFaces();
    zoneIds_ = surf.zoneIds_;
    zoneToc_ = surf.zoneToc_;
}


template<class Face>
void Foam::UnsortedMeshedSurface<Face>::operator=
(
    UnsortedMeshedSurface<Face>&& surf
)
{
    transfer();
}


template<class Face>
Foam::UnsortedMeshedSurface<Face>::operator
Foam::MeshedSurfaceProxy<Face>() const
{
    labelList faceMap;
    List<surfZone> zoneLst = this->sortedZones(faceMap);

    return MeshedSurfaceProxy<Face>
    (
        this->points(),
        this->surfFaces(),
        zoneLst,
        faceMap
    );
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Face>
Foam::Istream& Foam::operator>>
(
    Istream& is,
    UnsortedMeshedSurface<Face>& surf
)
{
    surf.readIstream(is);
    return is;
}


template<class Face>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const UnsortedMeshedSurface<Face>& surf
)
{
    surf.writeOstream(os);
    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "UnsortedMeshedSurfaceNew.C"

// ************************************************************************* //
