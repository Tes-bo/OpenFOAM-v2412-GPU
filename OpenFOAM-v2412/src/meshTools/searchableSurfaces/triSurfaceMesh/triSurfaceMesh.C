/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2020,2022 OpenCFD Ltd.
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

#include "triSurfaceMesh.H"
#include "Random.H"
#include "addToRunTimeSelectionTable.H"
#include "edgeHashes.H"
#include "triSurfaceFields.H"
#include "Time.H"
#include "PatchTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(triSurfaceMesh, 0);
    addToRunTimeSelectionTable(searchableSurface, triSurfaceMesh, dict);
}

Foam::word Foam::triSurfaceMesh::meshSubDir = "triSurface";


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::triSurfaceMesh::addFaceToEdge
(
    const edge& e,
    EdgeMap<label>& facesPerEdge
)
{
    //label& count = facesPerEdge(e, 0);  // lookup or new entry
    //if (count == 2)
    //{
    //    return false;
    //}
    //++count;
    //return true;

    auto fnd = facesPerEdge.find(e);
    if (fnd)
    {
        label& count = fnd.val();
        const int dir = edge::compare(e, fnd.key());
        if (dir == 1)
        {
            // Incorrect order. Mark with special value
            count = -1;
        }
        else if (dir == 0)
        {
            FatalErrorInFunction
                << "Incorrect matched edge " << fnd.key()
                << " to edge " << e << exit(FatalError);
            return false;
        }
        else if (count != -1)
        {
            if (count == 2)
            {
                return false;
            }
            ++count;
        }
    }
    else
    {
        facesPerEdge.insert(e, 1);
    }

    return true;
}


bool Foam::triSurfaceMesh::isSurfaceClosed() const
{
    if (debug)
    {
        Pout<< "triSurfaceMesh::isSurfaceClosed:"
            << " determining closedness for surface with "
            << triSurface::size() << " triangles" << endl;
    }

    const pointField& pts = triSurface::points();

/*
    // Construct pointFaces. Let's hope surface has compact point
    // numbering ...
    labelListList pointFaces;
    invertManyToMany(pts.size(), *this, pointFaces);

    // Loop over all faces surrounding point. Count edges emanating from point.
    // Every edge should be used by two faces exactly.
    // To prevent doing work twice per edge only look at edges to higher
    // point
    EdgeMap<label> facesPerEdge;
    forAll(pointFaces, pointi)
    {
        const labelList& pFaces = pointFaces[pointi];

        facesPerEdge.clear();
        for (const label facei : pFaces)
        {
            const triSurface::face_type& f = triSurface::operator[](facei);
            const label fp = f.find(pointi);

            // Something weird: if I expand the code of addFaceToEdge in both
            // below instances it gives a segmentation violation on some
            // surfaces. Compiler (4.3.2) problem?


            // Forward edge
            const label nextPointi = f[f.fcIndex(fp)];

            if (nextPointi > pointi)
            {
                bool okFace = addFaceToEdge
                (
                    edge(pointi, nextPointi),
                    facesPerEdge
                );

                if (!okFace)
                {
                    if (debug)
                    {
                        Pout<< "triSurfaceMesh::isSurfaceClosed :"
                            << " surface is open" << endl;
                    }
                    return false;
                }
            }

            // Reverse edge
            const label prevPointi = f[f.rcIndex(fp)];

            if (prevPointi > pointi)
            {
                bool okFace = addFaceToEdge
                (
                    edge(pointi, prevPointi),
                    facesPerEdge
                );

                if (!okFace)
                {
                    if (debug)
                    {
                        Pout<< "triSurfaceMesh::isSurfaceClosed :"
                            << " surface is open" << endl;
                    }
                    return false;
                }
            }
        }

        // Check for any edges used only once.
        forAllConstIters(facesPerEdge, iter)
        {
            if (iter.val() == -1)
            {
                // Special value for incorrect orientation
                if (debug)
                {
                    Pout<< "triSurfaceMesh::isSurfaceClosed :"
                        << " surface is closed but has inconsistent"
                        << " face orientation" << endl;
                }
                WarningInFunction << "Surface " << searchableSurface::name()
                    << " is closed but has inconsistent face orientation"
                    << " at edge " << pts[iter.key().first()]
                    << pts[iter.key().second()] << endl;
                return false;
            }
            else if (iter.val() != 2)
            {
                if (debug)
                {
                    Pout<< "triSurfaceMesh::isSurfaceClosed :"
                        << " surface is open" << endl;
                }
                return false;
            }
        }
    }
*/

    const triSurface& ts = *this;
    EdgeMap<label> facesPerEdge(2*ts.size());
    for (const auto& f : ts)
    {
        forAll(f, fp)
        {
            // Count number of occurences of edge between fp and fp+1
            const bool okFace = addFaceToEdge(f.edge(fp), facesPerEdge);

            if (!okFace)
            {
                if (debug)
                {
                    Pout<< "triSurfaceMesh::isSurfaceClosed :"
                        << " surface is non-manifold" << endl;
                }
                return false;
            }
        }
    }


    // Check for any edges used only once.
    bool haveWarned = false;
    forAllConstIters(facesPerEdge, iter)
    {
        if (iter.val() != 2 && iter.val() != -1)
        {
            if (debug)
            {
                Pout<< "triSurfaceMesh::isSurfaceClosed :"
                    << " surface is open" << endl;
            }
            return false;
        }
    }

    // Check for any edges with inconsistent normal orientation.
    forAllConstIters(facesPerEdge, iter)
    {
        if (iter.val() == -1)
        {
            // Special value for incorrect orientation
            if (!haveWarned)
            {
                WarningInFunction
                    << "Surface " << searchableSurface::name()
                    << " is closed but has inconsistent face orientation"
                    << nl
                    << "    at edge " << pts[iter.key().first()]
                    << pts[iter.key().second()]
                    << nl
                    << "    This means it probably cannot be used"
                    << " for inside/outside queries."
                    << " Suppressing further messages." << endl;
                haveWarned = true;
            }
            //- Return open?
            //return false;
        }
    }

    if (debug)
    {
        Pout<< "triSurfaceMesh::isSurfaceClosed :"
            << " surface is closed" << endl;
    }
    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::triSurfaceMesh::triSurfaceMesh(const IOobject& io, const triSurface& s)
:
    searchableSurface(io),
    objectRegistry
    (
        IOobject
        (
            io.name(),
            io.instance(),
            io.local(),
            io.db(),
            io.readOpt(),
            io.writeOpt(),
            false       // searchableSurface already registered under name
        )
    ),
    triSurface(s),
    triSurfaceRegionSearch(static_cast<const triSurface&>(*this)),
    minQuality_(-1),
    surfaceClosed_(-1),
    outsideVolType_(volumeType::UNKNOWN)
{
    bounds() = boundBox(triSurface::points(), false);
}


Foam::triSurfaceMesh::triSurfaceMesh(const IOobject& io)
:
    // Find instance for triSurfaceMesh
    searchableSurface(io),
    // Reused found instance in objectRegistry
    objectRegistry
    (
        IOobject
        (
            io.name(),
            searchableSurface::instance(),
            io.local(),
            io.db(),
            io.readOpt(),
            io.writeOpt(),
            false       // searchableSurface already registered under name
        )
    ),
    triSurface(static_cast<const searchableSurface&>(*this), dictionary::null),
    triSurfaceRegionSearch(static_cast<const triSurface&>(*this)),
    minQuality_(-1),
    surfaceClosed_(-1),
    outsideVolType_(volumeType::UNKNOWN)
{
    bounds() = boundBox(triSurface::points(), false);
}


Foam::triSurfaceMesh::triSurfaceMesh
(
    const IOobject& io,
    const dictionary& dict
)
:
    searchableSurface(io),
    // Reused found instance in objectRegistry
    objectRegistry
    (
        IOobject
        (
            io.name(),
            searchableSurface::instance(),
            io.local(),
            io.db(),
            io.readOpt(),
            io.writeOpt(),
            false       // searchableSurface already registered under name
        )
    ),
    triSurface(static_cast<const searchableSurface&>(*this), dict),
    triSurfaceRegionSearch(static_cast<const triSurface&>(*this), dict),
    minQuality_(-1),
    surfaceClosed_(-1),
    outsideVolType_(volumeType::UNKNOWN)
{
    // Read with supplied file name instead of objectPath/filePath
    if (dict.readIfPresent("file", fName_, keyType::LITERAL))
    {
        fName_ = triSurface::relativeFilePath
        (
            static_cast<const searchableSurface&>(*this),
            fName_,
            true
        );
    }

    // Report optional scale factor (eg, mm to m),
    // which was already applied during triSurface construction
    scalar scaleFactor{0};
    if (dict.getOrDefault("scale", scaleFactor) && scaleFactor > 0)
    {
        Info<< searchableSurface::name()
            << " : using scale " << scaleFactor << endl;
    }

    bounds() = boundBox(triSurface::points(), false);

    // Have optional minimum quality for normal calculation
    if (dict.readIfPresent("minQuality", minQuality_) && minQuality_ > 0)
    {
        Info<< searchableSurface::name()
            << " : ignoring triangles with quality < "
            << minQuality_ << " for normals calculation." << endl;
    }
}


Foam::triSurfaceMesh::triSurfaceMesh(const IOobject& io, const readAction r)
:
    // Find instance for triSurfaceMesh
    searchableSurface(io),
    // Reused found instance in objectRegistry
    objectRegistry
    (
        IOobject
        (
            io.name(),
            searchableSurface::instance(),
            io.local(),
            io.db(),
            io.readOpt(),
            io.writeOpt(),
            false       // searchableSurface already registered under name
        )
    ),
    triSurface(),
    triSurfaceRegionSearch(static_cast<const triSurface&>(*this)),
    minQuality_(-1),
    surfaceClosed_(-1),
    outsideVolType_(volumeType::UNKNOWN)
{
    // Check IO flags
    if (io.isAnyRead())
    {
        const bool searchGlobal(r == localOrGlobal || r == masterOnly);

        const fileName actualFile
        (
            searchGlobal
          ? io.globalFilePath(typeName)
          : io.localFilePath(typeName)
        );

        // const fileName actualFile
        // (
        //     io.filePath(searchGlobal, typeName)
        // );

        if (debug)
        {
            Pout<< "triSurfaceMesh(const IOobject& io) :"
                << " loading surface " << io.objectPath()
                << " local filePath:" << io.localFilePath(typeName)
                << " from:" << actualFile << endl;
        }

        if (searchGlobal && Pstream::parRun())
        {
            // Check where surface was found
            const fileName localFile(io.localFilePath(typeName));

            if
            (
                r == masterOnly
             && ((actualFile.empty() || actualFile != localFile))
            )
            {
                // Found undecomposed surface. Load on master only
                if (Pstream::master())
                {
                    triSurface s2(actualFile);
                    triSurface::transfer(s2);
                }
                Pstream::broadcast(triSurface::patches());
                if (debug)
                {
                    Pout<< "triSurfaceMesh(const IOobject& io) :"
                        << " loaded triangles:" << triSurface::size() << endl;
                }
            }
            else
            {
                // Read on all processors
                triSurface s2(actualFile);
                triSurface::transfer(s2);
                if (debug)
                {
                    Pout<< "triSurfaceMesh(const IOobject& io) :"
                        << " loaded triangles:" << triSurface::size() << endl;
                }
            }
        }
        else
        {
            // Read on all processors
            triSurface s2(actualFile);
            triSurface::transfer(s2);
            if (debug)
            {
                Pout<< "triSurfaceMesh(const IOobject& io) :"
                    << " loaded triangles:" << triSurface::size() << endl;
            }
        }
    }

    bounds() = boundBox(triSurface::points(), false);
}


Foam::triSurfaceMesh::triSurfaceMesh
(
    const IOobject& io,
    const dictionary& dict,
    const readAction r
)
:
    searchableSurface(io),
    // Reused found instance in objectRegistry
    objectRegistry
    (
        IOobject
        (
            io.name(),
            searchableSurface::instance(),
            io.local(),
            io.db(),
            io.readOpt(),
            io.writeOpt(),
            false       // searchableSurface already registered under name
        )
    ),
    triSurface(),
    triSurfaceRegionSearch(static_cast<const triSurface&>(*this), dict),
    minQuality_(-1),
    surfaceClosed_(-1),
    outsideVolType_(volumeType::UNKNOWN)
{
    if (io.isAnyRead())
    {
        // Surface type (optional)
        const word surfType(dict.getOrDefault<word>("fileType", word::null));

        // Scale factor (optional)
        const scalar scaleFactor(dict.getOrDefault<scalar>("scale", 0));

        const bool searchGlobal(r == localOrGlobal || r == masterOnly);

        fileName actualFile
        (
            searchGlobal
          ? io.globalFilePath(typeName)
          : io.localFilePath(typeName)
        );

        // const fileName actualFile
        // (
        //     io.filePath(searchGlobal, typeName)
        // );

        // Reading from supplied file name instead of objectPath/filePath
        if (dict.readIfPresent("file", fName_, keyType::LITERAL))
        {
            fName_ = relativeFilePath
            (
                static_cast<const searchableSurface&>(*this),
                fName_,
                searchGlobal
            );
            actualFile = fName_;
        }

        if (debug)
        {
            Pout<< "triSurfaceMesh(const IOobject& io, const dictionary&) :"
                << " loading surface " << io.objectPath()
                << " local filePath:" << io.localFilePath(typeName)
                << " from:" << actualFile << endl;
        }


        if (searchGlobal && Pstream::parRun())
        {
            // Check where surface was found. Bit tricky:
            // - master will have actualFile (in parent directory)
            //   different from localFilePath (in processor0/)
            // - slave might have actualFile empty and localFile empty

            const fileName localFile(io.localFilePath(typeName));

            if
            (
                r == masterOnly
             && ((actualFile.empty() || actualFile != localFile))
            )
            {
                // Surface not loaded from processor directories -> undecomposed
                // surface. Load on master only
                if (Pstream::master())
                {
                    triSurface s2(actualFile, surfType, scaleFactor);
                    triSurface::transfer(s2);
                }
                Pstream::broadcast(triSurface::patches());
                if (debug)
                {
                    Pout<< "triSurfaceMesh(const IOobject& io) :"
                        << " loaded triangles:" << triSurface::size() << endl;
                }
            }
            else
            {
                // Read on all processors
                triSurface s2(actualFile, surfType, scaleFactor);
                triSurface::transfer(s2);
                if (debug)
                {
                    Pout<< "triSurfaceMesh(const IOobject& io) :"
                        << " loaded triangles:" << triSurface::size() << endl;
                }
            }
        }
        else
        {
            // Read on all processors
            triSurface s2(actualFile, surfType, scaleFactor);
            triSurface::transfer(s2);
            if (debug)
            {
                Pout<< "triSurfaceMesh(const IOobject& io) :"
                    << " loaded triangles:" << triSurface::size() << endl;
            }
        }


        // Report optional scale factor (eg, mm to m),
        // which was already applied during triSurface reading
        if (scaleFactor > 0)
        {
            Info<< searchableSurface::name()
                << " : using scale " << scaleFactor << endl;
        }
    }


    bounds() = boundBox(triSurface::points(), false);

    // Have optional minimum quality for normal calculation
    if (dict.readIfPresent("minQuality", minQuality_) && minQuality_ > 0)
    {
        Info<< searchableSurface::name()
            << " : ignoring triangles with quality < "
            << minQuality_ << " for normals calculation." << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::triSurfaceMesh::~triSurfaceMesh()
{
    clearOut();
}


void Foam::triSurfaceMesh::clearOut()
{
    // Do not clear closedness status
    triSurfaceRegionSearch::clearOut();
    edgeTree_.clear();
    triSurface::clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField> Foam::triSurfaceMesh::coordinates() const
{
    auto tpts = tmp<pointField>::New();
    auto& pts = tpts.ref();

    if (triSurface::hasFaceCentres())
    {
        // Can reuse existing values
        pts = triSurface::faceCentres();
    }
    else
    {
        typedef SubList<labelledTri> FaceListType;

        // Calculate face centres from a copy to avoid incurring
        // additional storage
        pts = PrimitivePatch<FaceListType, const pointField&>
        (
            FaceListType(*this, triSurface::size()),
            triSurface::points()
        ).faceCentres();
    }

    return tpts;
}


void Foam::triSurfaceMesh::boundingSpheres
(
    pointField& centres,
    scalarField& radiusSqr
) const
{
    centres = coordinates();
    radiusSqr.setSize(size());
    radiusSqr = 0.0;

    const pointField& pts = triSurface::points();

    forAll(*this, facei)
    {
        const labelledTri& f = triSurface::operator[](facei);
        const point& fc = centres[facei];
        for (const label pointi : f)
        {
            radiusSqr[facei] = max(radiusSqr[facei], fc.distSqr(pts[pointi]));
        }
    }

    // Add a bit to make sure all points are tested inside
    radiusSqr += Foam::sqr(SMALL);
}


Foam::tmp<Foam::pointField> Foam::triSurfaceMesh::points() const
{
    return triSurface::points();
}


bool Foam::triSurfaceMesh::overlaps(const boundBox& bb) const
{
     const indexedOctree<treeDataTriSurface>& octree = tree();

     labelList indices = octree.findBox(treeBoundBox(bb));

     return !indices.empty();
}


void Foam::triSurfaceMesh::movePoints(const pointField& newPoints)
{
    if (debug)
    {
        Pout<< "triSurfaceMesh::movePoints :"
            << " moving at time " << objectRegistry::time().timeName()
            << endl;
    }

    // Preserve topological point status (surfaceClosed_, outsideVolType_)

    // Update local information (instance, event number)
    searchableSurface::instance() = objectRegistry::time().timeName();
    objectRegistry::instance() = searchableSurface::instance();

    const label event = getEvent();
    searchableSurface::eventNo() = event;
    objectRegistry::eventNo() = searchableSurface::eventNo();

    // Clear additional addressing
    triSurfaceRegionSearch::clearOut();
    edgeTree_.clear();
    triSurface::movePoints(newPoints);

    bounds() = boundBox(triSurface::points(), false);
    if (debug)
    {
        Pout<< "triSurfaceMesh::movePoints: finished moving points" << endl;
    }
}


const Foam::indexedOctree<Foam::treeDataEdge>&
Foam::triSurfaceMesh::edgeTree() const
{
    if (!edgeTree_)
    {
        if (debug)
        {
            Pout<< "triSurfaceMesh::edgeTree :"
                << " constructing tree for " << nEdges() - nInternalEdges()
                << " boundary edges" << endl;
        }

        // Boundary edges
        const labelRange bEdges(nInternalEdges(), nBoundaryEdges());

        treeBoundBox bb(point::zero);

        if (bEdges.size())
        {
            label nPoints;
            PatchTools::calcBounds
            (
                *this,
                bb,
                nPoints
            );

            // Random number generator. Bit dodgy since not exactly random ;-)
            Random rndGen(65431);

            // Slightly extended bb. Slightly off-centred just so on symmetric
            // geometry there are less face/edge aligned items.

            bb.inflate(rndGen, 1e-4, ROOTVSMALL);
        }


        if (debug)
        {
            Pout<< "triSurfaceMesh::edgeTree : "
                << "calculating edge tree for bb:" << bb << endl;
        }

        const scalar oldTol =
            indexedOctree<treeDataEdge>::perturbTol(tolerance());

        edgeTree_.reset
        (
            new indexedOctree<treeDataEdge>
            (
                // Boundary edges
                treeDataEdge(edges(), localPoints(), bEdges),

                bb,                 // bb
                maxTreeDepth(),     // maxLevel
                10,                 // leafsize
                3.0                 // duplicity
            )
        );

        indexedOctree<treeDataEdge>::perturbTol(oldTol);

        if (debug)
        {
            Pout<< "triSurfaceMesh::edgeTree :"
                << " finished constructing tree for "
                << nEdges() - nInternalEdges()
                << " boundary edges" << endl;
        }
    }

    return *edgeTree_;
}


const Foam::wordList& Foam::triSurfaceMesh::regions() const
{
    if (regions_.empty())
    {
        regions_.setSize(patches().size());
        forAll(regions_, regioni)
        {
            regions_[regioni] = patches()[regioni].name();
        }
    }
    return regions_;
}


bool Foam::triSurfaceMesh::hasVolumeType() const
{
    if (surfaceClosed_ == -1)
    {
        if (isSurfaceClosed())
        {
            surfaceClosed_ = 1;
        }
        else
        {
            surfaceClosed_ = 0;
        }
    }

    return surfaceClosed_ == 1;
}


Foam::volumeType Foam::triSurfaceMesh::outsideVolumeType() const
{
    if (outsideVolType_ == volumeType::UNKNOWN)
    {
        // Get point outside bounds()
        const point outsidePt(bounds().max() + 0.5*bounds().span());

        if (debug)
        {
            Pout<< "triSurfaceMesh::outsideVolumeType :"
                << " triggering outsidePoint:" << outsidePt
                << " orientation" << endl;
        }

        //outsideVolType_ = tree().shapes().getVolumeType(tree(), outsidePt);
        // Note: do not use tree directly so e.g. distributedTriSurfaceMesh
        //       has opportunity to intercept
        List<volumeType> outsideVolTypes;
        getVolumeType(pointField(1, outsidePt), outsideVolTypes);
        outsideVolType_ = outsideVolTypes[0];

        if (debug)
        {
            Pout<< "triSurfaceMesh::outsideVolumeType :"
                << " finished outsidePoint:" << outsidePt
                << " orientation:" << volumeType::names[outsideVolType_]
                << endl;
        }
    }

    return outsideVolType_;
}


void Foam::triSurfaceMesh::findNearest
(
    const pointField& samples,
    const scalarField& nearestDistSqr,
    List<pointIndexHit>& info
) const
{
    if (debug)
    {
        Pout<< "triSurfaceMesh::findNearest :"
            << " trying to find nearest for " << samples.size()
            << " samples with max sphere "
            << (samples.size() ? Foam::sqrt(max(nearestDistSqr)) : Zero)
            << endl;
    }
    triSurfaceSearch::findNearest(samples, nearestDistSqr, info);
    if (debug)
    {
        Pout<< "triSurfaceMesh::findNearest :"
            << " finished trying to find nearest for " << samples.size()
            << " samples" << endl;
    }
}


void Foam::triSurfaceMesh::findNearest
(
    const pointField& samples,
    const scalarField& nearestDistSqr,
    const labelList& regionIndices,
    List<pointIndexHit>& info
) const
{
    if (debug)
    {
        Pout<< "triSurfaceMesh::findNearest :"
            << " trying to find nearest and region for " << samples.size()
            << " samples with max sphere "
            << (samples.size() ? Foam::sqrt(max(nearestDistSqr)) : Zero)
            << endl;
    }
    triSurfaceRegionSearch::findNearest
    (
        samples,
        nearestDistSqr,
        regionIndices,
        info
    );
    if (debug)
    {
        Pout<< "triSurfaceMesh::findNearest :"
            << " finished trying to find nearest and region for "
            << samples.size() << " samples" << endl;
    }
}


void Foam::triSurfaceMesh::findLine
(
    const pointField& start,
    const pointField& end,
    List<pointIndexHit>& info
) const
{
    if (debug)
    {
        Pout<< "triSurfaceMesh::findLine :"
            << " intersecting with "
            << start.size() << " rays" << endl;
    }
    triSurfaceSearch::findLine(start, end, info);
    if (debug)
    {
        Pout<< "triSurfaceMesh::findLine :"
            << " finished intersecting with "
            << start.size() << " rays" << endl;
    }
}


void Foam::triSurfaceMesh::findLineAny
(
    const pointField& start,
    const pointField& end,
    List<pointIndexHit>& info
) const
{
    if (debug)
    {
        Pout<< "triSurfaceMesh::findLineAny :"
            << " intersecting with "
            << start.size() << " rays" << endl;
    }
    triSurfaceSearch::findLineAny(start, end, info);
    if (debug)
    {
        Pout<< "triSurfaceMesh::findLineAny :"
            << " finished intersecting with "
            << start.size() << " rays" << endl;
    }
}


void Foam::triSurfaceMesh::findLineAll
(
    const pointField& start,
    const pointField& end,
    List<List<pointIndexHit>>& info
) const
{
    if (debug)
    {
        Pout<< "triSurfaceMesh::findLineAll :"
            << " intersecting with "
            << start.size() << " rays" << endl;
    }
    triSurfaceSearch::findLineAll(start, end, info);
    if (debug)
    {
        Pout<< "triSurfaceMesh::findLineAll :"
            << " finished intersecting with "
            << start.size() << " rays" << endl;
    }
}


void Foam::triSurfaceMesh::getRegion
(
    const List<pointIndexHit>& info,
    labelList& region
) const
{
    if (debug)
    {
        Pout<< "triSurfaceMesh::getRegion :"
            << " getting region for "
            << info.size() << " triangles" << endl;
    }
    region.setSize(info.size());
    forAll(info, i)
    {
        if (info[i].hit())
        {
            region[i] = triSurface::operator[](info[i].index()).region();
        }
        else
        {
            region[i] = -1;
        }
    }
    if (debug)
    {
        Pout<< "triSurfaceMesh::getRegion :"
            << " finished getting region for "
            << info.size() << " triangles" << endl;
    }
}


void Foam::triSurfaceMesh::getNormal
(
    const List<pointIndexHit>& info,
    vectorField& normal
) const
{
    if (debug)
    {
        Pout<< "triSurfaceMesh::getNormal :"
            << " getting normal for "
            << info.size() << " triangles" << endl;
    }

    const triSurface& s = *this;
    const pointField& pts = s.points();

    normal.setSize(info.size());

    if (minQuality_ >= 0)
    {
        // Make sure we don't use triangles with low quality since
        // normal is not reliable.

        const labelListList& faceFaces = s.faceFaces();

        forAll(info, i)
        {
            if (info[i].hit())
            {
                const label facei = info[i].index();
                normal[i] = s[facei].unitNormal(pts);

                scalar qual = s[facei].tri(pts).quality();

                if (qual < minQuality_)
                {
                    // Search neighbouring triangles
                    const labelList& fFaces = faceFaces[facei];

                    for (const label nbri : fFaces)
                    {
                        scalar nbrQual = s[nbri].tri(pts).quality();
                        if (nbrQual > qual)
                        {
                            qual = nbrQual;
                            normal[i] = s[nbri].unitNormal(pts);
                        }
                    }
                }
            }
            else
            {
                // Set to what?
                normal[i] = Zero;
            }
        }
    }
    else
    {
        forAll(info, i)
        {
            if (info[i].hit())
            {
                const label facei = info[i].index();

                // Uncached
                normal[i] = s[facei].unitNormal(pts);
            }
            else
            {
                // Set to what?
                normal[i] = Zero;
            }
        }
    }
    if (debug)
    {
        Pout<< "triSurfaceMesh::getNormal :"
            << " finished getting normal for "
            << info.size() << " triangles" << endl;
    }
}


void Foam::triSurfaceMesh::setField(const labelList& values)
{
    auto* fldPtr = getObjectPtr<triSurfaceLabelField>("values");

    if (fldPtr)
    {
        (*fldPtr).field() = values;
    }
    else
    {
        fldPtr = new triSurfaceLabelField
        (
            IOobject
            (
                "values",
                objectRegistry::time().timeName(),  // instance
                meshSubDir,                         // local
                *this,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                IOobject::REGISTER
            ),
            *this,
            dimless,
            labelField(values)
        );

        // Store field on triMesh
        fldPtr->store();
    }
    if (debug)
    {
        Pout<< "triSurfaceMesh::setField :"
            << " finished setting field for "
            << values.size() << " triangles" << endl;
    }
}


void Foam::triSurfaceMesh::getField
(
    const List<pointIndexHit>& info,
    labelList& values
) const
{
    const auto* fldPtr = getObjectPtr<triSurfaceLabelField>("values");

    if (fldPtr)
    {
        const auto& fld = *fldPtr;

        values.setSize(info.size());

        forAll(info, i)
        {
            if (info[i].hit())
            {
                values[i] = fld[info[i].index()];
            }
        }
    }
    if (debug)
    {
        Pout<< "triSurfaceMesh::setField :"
            << " finished getting field for "
            << info.size() << " triangles" << endl;
    }
}


void Foam::triSurfaceMesh::getVolumeType
(
    const pointField& points,
    List<volumeType>& volType
) const
{
    const scalar oldTol =
        indexedOctree<treeDataTriSurface>::perturbTol(tolerance());

    if (debug)
    {
        Pout<< "triSurfaceMesh::getVolumeType :"
            << " finding orientation for " << points.size()
            << " samples" << endl;
    }

    volType.setSize(points.size());

    forAll(points, pointi)
    {
        const point& pt = points[pointi];

        if (tree().bb().contains(pt))
        {
            // Use cached volume type per each tree node
            volType[pointi] = tree().getVolumeType(pt);
        }
        else if (hasVolumeType())
        {
            // Precalculate and cache value for this outside point
            if (outsideVolType_ == volumeType::UNKNOWN)
            {
                outsideVolType_ = tree().shapes().getVolumeType(tree(), pt);
            }
            volType[pointi] = outsideVolType_;
        }
        else
        {
            // Have to calculate directly as outside the octree
            volType[pointi] = tree().shapes().getVolumeType(tree(), pt);
        }
    }

    indexedOctree<treeDataTriSurface>::perturbTol(oldTol);

    if (debug)
    {
        Pout<< "triSurfaceMesh::getVolumeType :"
            << " finished finding orientation for " << points.size()
            << " samples" << endl;
    }
}


bool Foam::triSurfaceMesh::writeObject
(
    IOstreamOption,
    const bool writeOnProc
) const
{
    const Time& runTime = searchableSurface::time();
    const fileName& instance = searchableSurface::instance();

    if
    (
        instance != runTime.timeName()
     && instance != runTime.system()
     && instance != runTime.caseSystem()
     && instance != runTime.constant()
     && instance != runTime.caseConstant()
    )
    {
        const_cast<triSurfaceMesh&>(*this).searchableSurface::instance() =
            runTime.timeName();
        const_cast<triSurfaceMesh&>(*this).objectRegistry::instance() =
            runTime.timeName();
    }

    fileName fullPath;
    if (fName_.size())
    {
        // Override file name

        fullPath = fName_;

        fullPath.expand();
        if (!fullPath.isAbsolute())
        {
            // Add directory from regIOobject
            fullPath = searchableSurface::objectPath().path()/fullPath;
        }
    }
    else
    {
        fullPath = searchableSurface::objectPath();
    }

    if (!mkDir(fullPath.path()))
    {
        return false;
    }

    triSurface::write(fullPath);

    if (!isFile(fullPath))
    {
        return false;
    }

    return true;
}


// ************************************************************************* //
