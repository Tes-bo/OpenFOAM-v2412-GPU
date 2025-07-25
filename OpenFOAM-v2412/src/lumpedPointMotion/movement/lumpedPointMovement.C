/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "lumpedPointMovement.H"
#include "lumpedPointIOMovement.H"
#include "Fstream.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "PtrMap.H"
#include "triFace.H"
#include "labelPair.H"
#include "indexedOctree.H"
#include "treeDataPoint.H"
#include "pointIndexHit.H"
#include "pointPatch.H"
#include "PstreamReduceOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

int Foam::lumpedPointMovement::debug
(
    ::Foam::debug::debugSwitch("lumpedPointMovement", 0)
);


const Foam::word
Foam::lumpedPointMovement::canonicalName("lumpedPointMovement");


const Foam::Enum
<
    Foam::lumpedPointMovement::outputFormatType
>
Foam::lumpedPointMovement::formatNames
({
    { outputFormatType::PLAIN, "plain" },
    { outputFormatType::DICTIONARY, "dictionary" },
});


const Foam::Enum
<
    Foam::lumpedPointMovement::scalingType
>
Foam::lumpedPointMovement::scalingNames
({
    { scalingType::LENGTH, "length" },
    { scalingType::FORCE, "force" },
    { scalingType::MOMENT, "moment" },
});



// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

//! \cond fileScope
//- Space-separated vector value (ASCII)
static inline Ostream& putPlain(Ostream& os, const vector& v)
{
    return os << v.x() << ' ' << v.y() << ' ' << v.z();
}


//! \cond fileScope
//- Write list content with size, bracket, content, bracket one-per-line.
//  This makes for consistent for parsing, regardless of the list length.
template <class T>
static void writeList(Ostream& os, const string& header, const UList<T>& list)
{
    const label len = list.size();

    // Header string
    os  << header.c_str() << nl;

    // Write size and start delimiter
    os  << len << nl << token::BEGIN_LIST << nl;

    // Write contents
    for (label i=0; i < len; ++i)
    {
        os << list[i] << nl;
    }

    // Write end delimiter
    os << token::END_LIST;
    os.endEntry() << nl;
}
//! \endcond

} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::lumpedPointMovement::lumpedPointMovement()
:
    origin_(Zero),
    state0_(),
    state_(),
    originalIds_(),
    controllers_(),
    patchControls_(),
    relax_(1),
    ownerId_(-1),
    forcesDict_(),
    coupler_(),
    inputName_("positions.in"),
    outputName_("forces.out"),
    logName_("movement.log"),
    inputFormat_(lumpedPointState::inputFormatType::DICTIONARY),
    outputFormat_(outputFormatType::DICTIONARY),
    scaleInput_(-1),
    scaleOutput_(-1),
    calcFrequency_(1),
    lastTrigger_(-1)
{}


Foam::lumpedPointMovement::lumpedPointMovement
(
    const dictionary& dict,
    label ownerId
)
:
    origin_(Zero),
    state0_(),
    state_(),
    originalIds_(),
    controllers_(),
    patchControls_(),
    relax_(1),
    ownerId_(ownerId),
    forcesDict_(),
    coupler_(),
    inputName_("positions.in"),
    outputName_("forces.out"),
    logName_("movement.log"),
    scaleInput_(-1),
    scaleOutput_(-1),
    calcFrequency_(1),
    lastTrigger_(-1)
{
    readDict(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::lumpedPointMovement::couplingPending(const label timeIndex) const
{
    return (timeIndex >= lastTrigger_ + calcFrequency_);
}


void Foam::lumpedPointMovement::couplingCompleted(const label timeIndex) const
{
    lastTrigger_ = timeIndex;
}


void Foam::lumpedPointMovement::readDict(const dictionary& dict)
{
    origin_ = dict.getOrDefault<point>("origin", Zero);

    // Initial point locations (zero rotation)
    auto tpoints0 = tmp<pointField>::New();
    auto& points0 = tpoints0.ref();

    dict.readEntry("points", points0);
    points0 += origin_;

    originalIds_.clear();
    controllers_.clear();
    patchControls_.clear();


    // The FEA ids for the points (optional)
    Map<label> pointIdMap;

    if (dict.readIfPresent("pointLabels", originalIds_) && originalIds_.size())
    {
        if (originalIds_.size() != points0.size())
        {
            FatalIOErrorInFunction(dict)
                << "Incorrect number of pointLabels. Had "
                << originalIds_.size() << " for " << points0.size() << " points"
                << nl << endl
                << exit(FatalIOError);
        }

        labelHashSet duplicates;

        label pointi = 0;

        for (const label id : originalIds_)
        {
            if (!pointIdMap.insert(id, pointi))
            {
                duplicates.set(id);
            }

            ++pointi;
        }

        if (!duplicates.empty())
        {
            FatalIOErrorInFunction(dict)
                << "Found duplicate point ids "
                << flatOutput(duplicates.sortedToc()) << nl << endl
                << exit(FatalIOError);
        }
    }

    const dictionary* dictptr = dict.findDict("controllers");
    if (dictptr)
    {
        controllers_ = HashPtrTable<lumpedPointController>(*dictptr);

        // Verify the input
        forAllIters(controllers_, iter)
        {
            (*iter)->remapPointLabels(points0.size(), pointIdMap);
        }
    }
    else
    {
        // Add single global controller
        // Warning?

        controllers_.clear();
    }

    relax_ = dict.getOrDefault<scalar>("relax", 1);
    relax_ = clamp(relax_, zero_one{});

    forcesDict_.merge(dict.subOrEmptyDict("forces"));

    const dictionary& commDict = dict.subDict("communication");
    coupler_.readDict(commDict);

    calcFrequency_ = commDict.getOrDefault<label>("calcFrequency", 1);

    // Leave trigger intact

    commDict.readEntry("inputName", inputName_);
    commDict.readEntry("outputName", outputName_);
    commDict.readIfPresent("logName", logName_);

    inputFormat_ = lumpedPointState::formatNames.get("inputFormat", commDict);
    outputFormat_ = formatNames.get("outputFormat", commDict);

    scaleInput_  = -1;
    scaleOutput_ = -1;

    const dictionary* scaleDict = nullptr;

    if ((scaleDict = commDict.findDict("scaleInput")) != nullptr)
    {
        for (int i=0; i < scaleInput_.size(); ++i)
        {
            const word& key = scalingNames[scalingType(i)];

            if
            (
                scaleDict->readIfPresent(key, scaleInput_[i])
             && scaleInput_[i] > 0
            )
            {
                Info<<"Using input " << key << " multiplier: "
                    << scaleInput_[i] << nl;
            }
        }
    }

    if ((scaleDict = commDict.findDict("scaleOutput")) != nullptr)
    {
        for (int i=0; i < scaleOutput_.size(); ++i)
        {
            const word& key = scalingNames[scalingType(i)];

            if
            (
                scaleDict->readIfPresent(key, scaleOutput_[i])
             && scaleOutput_[i] > 0
            )
            {
                Info<<"Using output " << key << " multiplier: "
                    << scaleOutput_[i] << nl;
            }
        }
    }

    state0_ = lumpedPointState
    (
        tpoints0,
        quaternion::eulerOrderNames.getOrDefault
        (
            "rotationOrder",
            dict,
            quaternion::eulerOrder::ZXZ
        ),
        dict.getOrDefault("degrees", false)
    );

    state0_.scalePoints(scaleInput_[scalingType::LENGTH]);

    state_ = state0_;
}


bool Foam::lumpedPointMovement::hasPatchControl
(
    const polyPatch& pp
) const
{
    return hasPatchControl(pp.index());
}


bool Foam::lumpedPointMovement::hasInterpolator
(
    const pointPatch& fpatch
) const
{
    return hasInterpolator(fpatch.index());
}


void Foam::lumpedPointMovement::checkPatchControl
(
    const polyPatch& pp
) const
{
    const auto ctrlIter = patchControls_.cfind(pp.index());

    if (!ctrlIter.good())
    {
        FatalErrorInFunction
            << "No controllers for patch " << pp.name()
            << exit(FatalError);
    }

    const patchControl& ctrl = *ctrlIter;

    for (const word& ctrlName : ctrl.names_)
    {
        const auto iter = controllers_.cfind(ctrlName);

        if (!iter.good())
        {
            FatalErrorInFunction
                << "No controller: " << ctrlName << nl
                << " For patch " << pp.name()
                << exit(FatalError);
        }
    }
}


void Foam::lumpedPointMovement::setPatchControl
(
    const polyPatch& pp,
    const wordList& ctrlNames,
    const pointField& points0
)
{
    // Reference mass centres
    const pointField& lumpedCentres0 = state0().points();

    const label patchIndex = pp.index();

    // Info<<"Add patch control for patch " << patchIndex << " "
    //     << pp.name() << nl;

    patchControl& ctrl = patchControls_(patchIndex);
    ctrl.names_ = ctrlNames;

    labelList& faceToPoint = ctrl.faceToPoint_;
    faceToPoint.resize(pp.size(), -1);

    checkPatchControl(pp);

    const faceList& faces = pp.boundaryMesh().mesh().faces();

    // Subset of points to search (if specified)
    labelHashSet subsetPointIds;

    for (const word& ctrlName : ctrl.names_)
    {
        const auto iter = controllers_.cfind(ctrlName);

        if (!iter.good())
        {
            FatalErrorInFunction
                << "No controller: " << ctrlName << nl
                << exit(FatalError);
        }

        const labelList& pointLabels = (*iter)->pointLabels();

        subsetPointIds.insert(pointLabels);
    }

    if (ctrl.names_.size() && subsetPointIds.empty())
    {
        FatalErrorInFunction
            << "Controllers specified, but without any points" << nl
            << exit(FatalError);
    }


    treeBoundBox bb(lumpedCentres0);
    bb.inflate(0.01);

    indexedOctree<treeDataPoint> ppTree
    (
        treeDataPoint
        (
            lumpedCentres0,
            subsetPointIds.sortedToc(),
            subsetPointIds.size()   // subset is optional/required
        ),
        bb,     // overall search domain
        8,      // maxLevel
        10,     // leafsize
        3.0     // duplicity
    );
    const auto& treePoints = ppTree.shapes();

    const scalar searchDistSqr(sqr(GREAT));

    const label patchStart = pp.start();
    forAll(pp, patchFacei)
    {
        const point fc(faces[patchStart + patchFacei].centre(points0));

        // Store the original pointId, not subset id
        faceToPoint[patchFacei] =
            treePoints.objectIndex
            (
                ppTree.findNearest(fc, searchDistSqr).index()
            );
    }

    if (debug)
    {
        Pout<<"Added face mapping for patch: " << patchIndex << endl;
    }
}


void Foam::lumpedPointMovement::setInterpolator
(
    const pointPatch& fpatch,
    const pointField& points0
)
{
    // Reference mass centres
    const pointField& lumpedCentres0 = state0().points();

    const label patchIndex = fpatch.index();

    patchControl& ctrl = patchControls_(patchIndex);

    List<lumpedPointInterpolator>& interpList = ctrl.interp_;
    interpList.clear();

    // The connectivity, adjacency list
    Map<labelHashSet> adjacency;

    // Subset of points to search (if specified)
    labelHashSet subsetPointIds;

    for (const word& ctrlName : ctrl.names_)
    {
        const auto iter = controllers_.cfind(ctrlName);

        if (!iter.good())
        {
            FatalErrorInFunction
                << "No controller: " << ctrlName << nl
                << exit(FatalError);
        }

        const labelList& pointLabels = (*iter)->pointLabels();

        subsetPointIds.insert(pointLabels);

        // Adjacency lists
        forAll(pointLabels, i)
        {
            const label curr = pointLabels[i];

            if (i)
            {
                const label prev = pointLabels[i-1];
                adjacency(prev).insert(curr);
                adjacency(curr).insert(prev);
            }
            else if (!adjacency.found(curr))
            {
                adjacency(curr).clear();
            }
        }
    }

    if (ctrl.names_.empty())
    {
        // Adjacency lists
        const label len = state0().size();

        for (label i=0; i < len; ++i)
        {
            const label curr = i;

            if (i)
            {
                const label prev = i-1;
                adjacency(prev).insert(curr);
                adjacency(curr).insert(prev);
            }
            else if (!adjacency.found(curr))
            {
                adjacency(curr).clear();
            }
        }
    }

    if (ctrl.names_.size() && subsetPointIds.empty())
    {
        FatalErrorInFunction
            << "Controllers specified, but without any points" << nl
            << exit(FatalError);
    }


    // Pairs defining connecting points as triangles
    Map<labelPairList> adjacencyPairs;

    barycentric2D bary;

    {
        // Pairs for the ends
        DynamicList<labelPair> pairs;

        // Mag sin(angle) around the triangle point 0,
        // used to sort the generated triangles according to the acute angle
        DynamicList<scalar> acuteAngles;

        forAllConstIters(adjacency, iter)
        {
            const label nearest = iter.key();

            labelList neighbours(iter.val().sortedToc());

            const label len = neighbours.size();

            pairs.clear();
            acuteAngles.clear();

            const point& nearPt = lumpedCentres0[nearest];

            for (label j=1; j < len; ++j)
            {
                for (label i=0; i < j; ++i)
                {
                    labelPair neiPair(neighbours[i], neighbours[j]);

                    triPointRef tri
                    (
                        nearPt,
                        lumpedCentres0[neiPair.first()],
                        lumpedCentres0[neiPair.second()]
                    );

                    // Require non-degenerate triangles
                    if (tri.pointToBarycentric(tri.a(), bary) > SMALL)
                    {
                        // Triangle OK
                        pairs.append(neiPair);

                        vector ab(normalised(tri.b() - tri.a()));
                        vector ac(normalised(tri.c() - tri.a()));

                        // Angle between neighbouring edges
                        // Use negative cosine to map [0-180] -> [-1 .. +1]

                        acuteAngles.append(-(ab & ac));
                    }
                }
            }

            if (pairs.size() > 1)
            {
                // Sort by acute angle
                labelList order(sortedOrder(acuteAngles));
                inplaceReorder(order, pairs);
            }

            adjacencyPairs.insert(nearest, pairs);
        }
    }

    if (debug & 2)
    {
        Info<< "Adjacency table for patch: " << fpatch.name() << nl;

        for (const label own : adjacency.sortedToc())
        {
            Info<< own << " =>";
            for (const label nei : adjacency[own].sortedToc())
            {
                Info<< ' ' << nei;
            }

            if (originalIds_.size())
            {
                Info<< "  # " << originalIds_[own] << " =>";
                for (const label nei : adjacency[own].sortedToc())
                {
                    Info<< ' ' << originalIds_[nei];
                }
            }

            Info<< "  # tri " << flatOutput(adjacencyPairs[own]);
            Info<< nl;
        }
    }

    treeBoundBox bb(lumpedCentres0);
    bb.inflate(0.01);

    indexedOctree<treeDataPoint> ppTree
    (
        treeDataPoint
        (
            lumpedCentres0,
            subsetPointIds.sortedToc(),
            subsetPointIds.size()   // subset is optional/required
        ),
        bb,     // overall search domain
        8,      // maxLevel
        10,     // leafsize
        3.0     // duplicity
    );
    const auto& treePoints = ppTree.shapes();


    // Searching

    const scalar searchDistSqr(sqr(GREAT));

    const labelList& meshPoints = fpatch.meshPoints();

    interpList.resize(meshPoints.size());

    DynamicList<scalar> unsortedNeiWeightDist;
    DynamicList<label>  unsortedNeighbours;

    forAll(meshPoints, pointi)
    {
        const point& ptOnMesh = points0[meshPoints[pointi]];

        // Nearest (original) point id
        const label nearest =
            treePoints.objectIndex
            (
                ppTree.findNearest(ptOnMesh, searchDistSqr).index()
            );

        interpList[pointi].nearest(nearest);

        if (nearest == -1)
        {
            // Should not really happen
            continue;
        }

        // Have the nearest lumped point, now find the next nearest
        // but check that the direction is also correct.

        // OK: within the 0-1 bounds
        // 1+----+0
        //      .
        //     .
        //    +pt

        const point& nearPt = lumpedCentres0[nearest];

        const linePointRef toMeshPt(nearPt, ptOnMesh);

        const labelPairList& adjacentPairs = adjacencyPairs[nearest];

        unsortedNeighbours = adjacency[nearest].toc();
        unsortedNeiWeightDist.resize(unsortedNeighbours.size());

        forAll(unsortedNeighbours, nbri)
        {
            unsortedNeiWeightDist[nbri] =
                magSqr(ptOnMesh - lumpedCentres0[unsortedNeighbours[nbri]]);
        }

        // Sort by distance
        labelList distOrder(sortedOrder(unsortedNeiWeightDist));

        label ngood = 0;

        // Recalculate distance as weighting
        for (const label nbri : distOrder)
        {
            const label nextPointi = unsortedNeighbours[nbri];

            const point& nextPt = lumpedCentres0[nextPointi];

            const linePointRef toNextPt(nearPt, nextPt);

            const scalar weight =
                (toMeshPt.vec() & toNextPt.unitVec()) / toNextPt.mag();

            if (weight < 0)
            {
                // Reject: wrong direction or other bad value
                continue;
            }
            else
            {
                // Store weight
                unsortedNeiWeightDist[nbri] = weight;

                // Retain good weight
                distOrder[ngood] = nbri;
                ++ngood;
            }
        }

        distOrder.resize(ngood);

        if (distOrder.size() < 1)
        {
            continue;
        }

        UIndirectList<label>  neighbours(unsortedNeighbours, distOrder);
        UIndirectList<scalar> neiWeight(unsortedNeiWeightDist, distOrder);

        bool useFirst = true;

        if (neighbours.size() > 1 && adjacentPairs.size())
        {
            // Check for best two neighbours

            bitSet neiPointid;
            neiPointid.set(neighbours);

            for (const labelPair& ends : adjacentPairs)
            {
                label nei1 = ends.first();
                label nei2 = ends.second();

                if (!neiPointid.test(nei1) || !neiPointid.test(nei2))
                {
                    // Reject, invalid combination for this point location
                    continue;
                }
                else if (neighbours.find(nei2) < neighbours.find(nei1))
                {
                    // Order by distance, which is not really needed,
                    // but helps with diagnostics
                    std::swap(nei1, nei2);
                }


                triFace triF(nearest, nei1, nei2);

                if
                (
                    triF.tri(lumpedCentres0).pointToBarycentric(ptOnMesh, bary)
                  > SMALL
                 && !bary.outside()
                )
                {
                    // Use barycentric weights
                    interpList[pointi].set(triF, bary);

                    useFirst = false;
                    break;
                }
            }
        }

        if (useFirst)
        {
            // Use geometrically closest neighbour
            interpList[pointi].next(neighbours.first(), neiWeight.first());
        }
    }
}


Foam::List<Foam::scalar> Foam::lumpedPointMovement::areas
(
    const polyMesh& pmesh
) const
{
    const label nLumpedPoints = state0().size();

    List<scalar> zoneAreas(nLumpedPoints, Zero);

    if (patchControls_.empty())
    {
        WarningInFunction
            << "Attempted to calculate areas without setMapping()"
            << nl << endl;
        return zoneAreas;
    }

    const polyBoundaryMesh& patches = pmesh.boundaryMesh();

    // fvMesh and has pressure field
    if (isA<fvMesh>(pmesh))
    {
        const auto& mesh = refCast<const fvMesh>(pmesh);

        // Face areas (on patches)
        const surfaceVectorField::Boundary& patchSf =
            mesh.Sf().boundaryField();

        forAllConstIters(patchControls_, iter)
        {
            const label patchIndex = iter.key();
            const patchControl& ctrl = iter.val();

            const labelList& faceToPoint = ctrl.faceToPoint_;

            const polyPatch& pp = patches[patchIndex];

            forAll(pp, patchFacei)
            {
                // Force from the patch-face into sum
                const label pointIndex = faceToPoint[patchFacei];

                if (pointIndex < 0)
                {
                    // Unmapped, for whatever reason?
                    continue;
                }

                // Force from the patch-face into sum
                zoneAreas[pointIndex] += mag(patchSf[patchIndex][patchFacei]);
            }
        }
    }

    Pstream::listCombineReduce(zoneAreas, plusEqOp<scalar>());

    return zoneAreas;
}


bool Foam::lumpedPointMovement::forcesAndMoments
(
    const polyMesh& pmesh,
    List<vector>& forces,
    List<vector>& moments
) const
{
    const label nLumpedPoints = state0().size();

    forces.resize(nLumpedPoints);
    moments.resize(nLumpedPoints);

    if (patchControls_.empty())
    {
        WarningInFunction
            << "Attempted to calculate forces without setMapping()"
            << nl << endl;

        forces.resize(nLumpedPoints, Zero);
        moments.resize(nLumpedPoints, Zero);
        return false;
    }

    // Initialize with zero
    forces = Zero;
    moments = Zero;

    // Current mass centres
    const pointField& lumpedCentres = state().points();

    const polyBoundaryMesh& patches = pmesh.boundaryMesh();

    const word pName(forcesDict_.getOrDefault<word>("p", "p"));
    scalar pRef   = forcesDict_.getOrDefault<scalar>("pRef",   0);
    scalar rhoRef = forcesDict_.getOrDefault<scalar>("rhoRef", 1);


    // Calculated force per patch - cache
    PtrMap<vectorField> forceOnPatches;

    const volScalarField* pPtr = pmesh.findObject<volScalarField>(pName);

    // fvMesh and has pressure field
    if (isA<fvMesh>(pmesh) && pPtr)
    {
        const auto& mesh = refCast<const fvMesh>(pmesh);
        const volScalarField& p = *pPtr;

        // Face centres (on patches)
        const surfaceVectorField::Boundary& patchCf = mesh.Cf().boundaryField();

        // Face areas (on patches)
        const surfaceVectorField::Boundary& patchSf = mesh.Sf().boundaryField();

        // Pressure (on patches)
        const volScalarField::Boundary& patchPress = p.boundaryField();

        // rhoRef if the pressure field is dynamic, i.e. p/rho otherwise 1
        rhoRef = (p.dimensions() == dimPressure ? 1.0 : rhoRef);

        // Scale pRef by density for incompressible simulations
        pRef /= rhoRef;

        forAllConstIters(patchControls_, iter)
        {
            const label patchIndex = iter.key();
            const patchControl& ctrl = iter.val();

            const labelList& faceToPoint = ctrl.faceToPoint_;

            if (!forceOnPatches.found(patchIndex))
            {
                // Patch faces are +ve outwards,
                // so the forces (exerted by fluid on solid)
                // already have the correct sign
                forceOnPatches.set
                (
                    patchIndex,
                    (
                        rhoRef
                      * patchSf[patchIndex] * (patchPress[patchIndex] - pRef)
                    ).ptr()
                );
            }

            const vectorField& forceOnPatch = *forceOnPatches[patchIndex];

            const polyPatch& pp = patches[patchIndex];

            forAll(pp, patchFacei)
            {
                // Force from the patch-face into sum
                const label pointIndex = faceToPoint[patchFacei];

                if (pointIndex < 0)
                {
                    // Unmapped, for whatever reason?
                    continue;
                }

                // Force from the patch-face into sum
                forces[pointIndex] += forceOnPatch[patchFacei];

                // Effective torque arm:
                // - translated into the lumped-points coordinate system
                //   prior to determining the distance
                const vector lever
                (
                    patchCf[patchIndex][patchFacei]
                  - lumpedCentres[pointIndex]
                );

                // Moment from the patch-face into sum
                moments[pointIndex] += lever ^ forceOnPatch[patchFacei];
            }
        }
    }
    else
    {
        Info<<"No pressure field" << endl;
    }

    Pstream::listCombineReduce(forces, plusEqOp<vector>());
    Pstream::listCombineReduce(moments, plusEqOp<vector>());

    return true;
}


Foam::tmp<Foam::pointField>
Foam::lumpedPointMovement::pointsDisplacement
(
    const pointPatch& fpatch,
    const pointField& points0
) const
{
    return pointsDisplacement(state(), fpatch, points0);
}


Foam::tmp<Foam::pointField>
Foam::lumpedPointMovement::pointsDisplacement
(
    const lumpedPointState& state,
    const pointPatch& fpatch,
    const pointField& points0
) const
{
    const label patchIndex = fpatch.index();

    // Reference mass centres
    const pointField& lumpedCentres0 = state0().points();

    // Current mass centres
    const pointField& lumpedCentres = state.points();

    // The local-to-global transformation tensor
    const tensorField& localToGlobal = state.rotations();

    const labelList& meshPoints = fpatch.meshPoints();

    // Could also verify the sizes (state vs original)

    auto tdisp = tmp<pointField>::New(fpatch.size());
    auto& disp = tdisp.ref();

    // The interpolator
    const List<lumpedPointInterpolator>& interpList
        = patchControls_[patchIndex].interp_;

    forAll(meshPoints, pointi)
    {
        const lumpedPointInterpolator& interp = interpList[pointi];

        const point& p0 = points0[meshPoints[pointi]];

        const vector origin0 = interp.interpolate(lumpedCentres0);
        const vector origin = interp.interpolate(lumpedCentres);
        const tensor rotTensor = interp.interpolate(localToGlobal);

        disp[pointi] = (rotTensor & (p0 - origin0)) + origin - p0;
    }

    return tdisp;
}


Foam::tmp<Foam::pointField>
Foam::lumpedPointMovement::pointsPosition
(
    const lumpedPointState& state,
    const pointPatch& fpatch,
    const pointField& points0
) const
{
    const label patchIndex = fpatch.index();

    // Reference mass centres
    const pointField& lumpedCentres0 = state0().points();

    // Current mass centres
    const pointField& lumpedCentres = state.points();

    // The local-to-global transformation tensor
    const tensorField& localToGlobal = state.rotations();

    const labelList& meshPoints = fpatch.meshPoints();

    // Could also verify the sizes (state vs original)

    auto tdisp = tmp<pointField>::New(fpatch.size());
    auto& disp = tdisp.ref();

    // The interpolator
    const List<lumpedPointInterpolator>& interpList =
        patchControls_[patchIndex].interp_;

    forAll(meshPoints, pointi)
    {
        const lumpedPointInterpolator& interp = interpList[pointi];

        const point& p0 = points0[meshPoints[pointi]];

        const vector origin0 = interp.interpolate(lumpedCentres0);
        const vector origin = interp.interpolate(lumpedCentres);
        const tensor rotTensor = interp.interpolate(localToGlobal);

        disp[pointi] = (rotTensor & (p0 - origin0)) + origin;
    }

    return tdisp;
}


void Foam::lumpedPointMovement::writeDict(Ostream& os) const
{
    // os.writeEntry("axis", axis_);
    // os.writeEntry("locations", locations_);
}


bool Foam::lumpedPointMovement::readState()
{
    lumpedPointState prev = state_;

    const bool status = state_.readData
    (
        inputFormat_,
        coupler().resolveFile(inputName_),
        state0().rotationOrder(),
        state0().degrees()
    );

    scalePoints(state_);

    state_.relax(relax_, prev);

    return status;
}


bool Foam::lumpedPointMovement::writeData
(
    Ostream& os,
    const UList<vector>& forces,
    const UList<vector>& moments,
    const outputFormatType fmt,
    const Tuple2<scalar, scalar>* timesWritten
) const
{
    const bool writeMoments = (moments.size() == forces.size());

    if (fmt == outputFormatType::PLAIN)
    {
        os  <<"########" << nl;
        if (timesWritten)
        {
            os  << "# Time value=" << timesWritten->first() << nl
                << "# Time prev=" << timesWritten->second() << nl;
        }
        os  << "# size=" << this->size() << nl
            << "# columns (points) (forces)";

        if (writeMoments)
        {
            os << " (moments)";
        }

        os << nl;

        bool report = false;
        scalar scaleLength = scaleOutput_[scalingType::LENGTH];
        scalar scaleForce  = scaleOutput_[scalingType::FORCE];
        scalar scaleMoment = scaleOutput_[scalingType::MOMENT];

        if (scaleLength > 0)
        {
            report = true;
        }
        else
        {
            scaleLength = 1.0;
        }

        if (scaleForce > 0)
        {
            report = true;
        }
        else
        {
            scaleForce = 1.0;
        }

        if (writeMoments)
        {
            if (scaleMoment > 0)
            {
                report = true;
            }
            else
            {
                scaleMoment = 1.0;
            }
        }

        if (report)
        {
            os  <<"# scaling points=" << scaleLength
                <<" forces=" << scaleForce;

            if (writeMoments)
            {
                os  <<" moments=" << scaleMoment;
            }

            os << nl;
        }

        os <<"########" << nl;

        forAll(state0().points(), i)
        {
            const point& pt = state0().points()[i];

            putPlain(os, scaleLength * pt) << ' ';

            if (i < forces.size())
            {
                const vector val(scaleForce * forces[i]);
                putPlain(os, val);
            }
            else
            {
                putPlain(os, vector::zero);
            }

            if (writeMoments)
            {
                os << ' ';
                if (i < moments.size())
                {
                    const vector val(scaleMoment * moments[i]);
                    putPlain(os, val);
                }
                else
                {
                    putPlain(os, vector::zero);
                }
            }

            os << nl;
        }
    }
    else
    {
        // Make it easier for external programs to parse
        // - exclude the usual OpenFOAM 'FoamFile' header
        // - ensure lists have consistent format

        os  <<"////////" << nl;
        if (timesWritten)
        {
            os.writeEntry("time", timesWritten->first());
            os.writeEntry("prevTime", timesWritten->second());
        }
        os  << nl;

        writeList(os, "points", state0().points());
        writeList(os, "forces", forces);

        if (writeMoments)
        {
            writeList(os, "moments", moments);
        }
    }

    return true;
}


bool Foam::lumpedPointMovement::writeData
(
    const UList<vector>& forces,
    const UList<vector>& moments,
    const Tuple2<scalar, scalar>* timesWritten
) const
{
    if (!UPstream::master())
    {
        return false;
    }

    // Regular output
    {
        OFstream os
        (
            coupler().resolveFile(outputName_)
        );

        writeData(os, forces, moments, outputFormat_, timesWritten);
    }

    // Log output - simple append to existing
    {
        OFstream os
        (
            coupler().resolveFile(logName_),
            IOstreamOption(),
            IOstreamOption::APPEND_ATE
        );

        writeData(os, forces, moments, outputFormatType::PLAIN, timesWritten);
    }

    return true;
}


// ************************************************************************* //
