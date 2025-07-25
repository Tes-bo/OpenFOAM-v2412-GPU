/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2015-2022 OpenCFD Ltd.
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

#include "shellSurfaces.H"
#include "searchableSurface.H"
#include "boundBox.H"
#include "triSurfaceMesh.H"
#include "refinementSurfaces.H"
#include "searchableSurfaces.H"
#include "orientedSurface.H"
#include "pointIndexHit.H"
#include "volumeType.H"
#include "distributedTriSurfaceMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::shellSurfaces::refineMode
>
Foam::shellSurfaces::refineModeNames_
({
    { refineMode::INSIDE, "inside" },
    { refineMode::OUTSIDE, "outside" },
    { refineMode::DISTANCE, "distance" },
});


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::shellSurfaces::setAndCheckLevels
(
    const label shellI,
    const List<Tuple2<scalar, label>>& distLevels
)
{
    const searchableSurface& shell = allGeometry_[shells_[shellI]];

    if (modes_[shellI] != DISTANCE && distLevels.size() != 1)
    {
        FatalErrorInFunction
            << "For refinement mode "
            << refineModeNames_[modes_[shellI]]
            << " specify only one distance+level."
            << " (its distance gets discarded)"
            << exit(FatalError);
    }
    // Extract information into separate distance and level
    distances_[shellI].setSize(distLevels.size());
    levels_[shellI].setSize(distLevels.size());

    forAll(distLevels, j)
    {
        distances_[shellI][j] = distLevels[j].first();
        levels_[shellI][j] = distLevels[j].second();

        if (levels_[shellI][j] < -1)
        {
            FatalErrorInFunction
                << "Shell " << shell.name()
                << " has illegal refinement level "
                << levels_[shellI][j]
                << exit(FatalError);
        }


        // Check in incremental order
        if (j > 0)
        {
            if
            (
                (distances_[shellI][j] <= distances_[shellI][j-1])
             || (levels_[shellI][j] > levels_[shellI][j-1])
            )
            {
                FatalErrorInFunction
                    << "For refinement mode "
                    << refineModeNames_[modes_[shellI]]
                    << " : Refinement should be specified in order"
                    << " of increasing distance"
                    << " (and decreasing refinement level)." << endl
                    << "Distance:" << distances_[shellI][j]
                    << " refinementLevel:" << levels_[shellI][j]
                    << exit(FatalError);
            }
        }
    }

    if (modes_[shellI] == DISTANCE)
    {
        if (!dryRun_)
        {
            Info<< "Refinement level according to distance to "
                << shell.name() << endl;
            forAll(levels_[shellI], j)
            {
                Info<< "    level " << levels_[shellI][j]
                    << " for all cells within " << distances_[shellI][j]
                    << " metre." << endl;
            }
        }
    }
    else
    {
        if (!shell.hasVolumeType())
        {
            FatalErrorInFunction
                << "Shell " << shell.name()
                << " does not support testing for "
                << refineModeNames_[modes_[shellI]] << endl
                << "Probably it is not closed."
                << exit(FatalError);
        }

        if (!dryRun_)
        {
            if (modes_[shellI] == INSIDE)
            {
                Info<< "Refinement level " << levels_[shellI][0]
                    << " for all cells inside " << shell.name() << endl;
            }
            else
            {
                Info<< "Refinement level " << levels_[shellI][0]
                    << " for all cells outside " << shell.name() << endl;
            }
        }
    }
}


// Specifically orient triSurfaces using a calculated point outside.
// Done since quite often triSurfaces not of consistent orientation which
// is (currently) necessary for sideness calculation
void Foam::shellSurfaces::orient()
{
    // Determine outside point.
    boundBox overallBb;

    bool hasSurface = false;

    forAll(shells_, shellI)
    {
        const searchableSurface& s = allGeometry_[shells_[shellI]];

        if
        (
            modes_[shellI] != DISTANCE
        &&  isA<triSurfaceMesh>(s)
        && !isA<distributedTriSurfaceMesh>(s)
        )
        {
            hasSurface = true;

            const triSurfaceMesh& shell = refCast<const triSurfaceMesh>(s);
            if (shell.triSurface::size())
            {
                // Assume surface is compact!
                overallBb.add(s.bounds());
            }
        }
    }

    if (returnReduceOr(hasSurface))
    {
        const point outsidePt = overallBb.max() + overallBb.span();

        //Info<< "Using point " << outsidePt << " to orient shells" << endl;

        forAll(shells_, shellI)
        {
            const searchableSurface& s = allGeometry_[shells_[shellI]];

            if
            (
                modes_[shellI] != DISTANCE
            &&  isA<triSurfaceMesh>(s)
            && !isA<distributedTriSurfaceMesh>(s)
            )
            {
                triSurfaceMesh& shell = const_cast<triSurfaceMesh&>
                (
                    refCast<const triSurfaceMesh>(s)
                );

                // Flip surface so outsidePt is outside.
                bool anyFlipped = orientedSurface::orient
                (
                    shell,
                    outsidePt,
                    true
                );

                if (anyFlipped && !dryRun_)
                {
                    // orientedSurface will have done a clearOut of the surface.
                    // we could do a clearout of the triSurfaceMeshes::trees()
                    // but these aren't affected by orientation
                    // (except for cached
                    // sideness which should not be set at this point.
                    // !!Should check!)

                    Info<< "shellSurfaces : Flipped orientation of surface "
                        << s.name()
                        << " so point " << outsidePt << " is outside." << endl;
                }
            }
        }
    }
}


// Find maximum level of a shell.
void Foam::shellSurfaces::findHigherLevel
(
    const pointField& pt,
    const label shellI,
    labelList& maxLevel
) const
{
    const labelList& levels = levels_[shellI];

    if (modes_[shellI] == DISTANCE)
    {
        // Distance mode.

        const scalarField& distances = distances_[shellI];

        // Collect all those points that have a current maxLevel less than
        // (any of) the shell. Also collect the furthest distance allowable
        // to any shell with a higher level.

        labelList candidateMap(pt.size());
        scalarField candidateDistSqr(pt.size());
        label candidateI = 0;

        forAll(maxLevel, pointi)
        {
            forAllReverse(levels, levelI)
            {
                if (levels[levelI] > maxLevel[pointi])
                {
                    candidateMap[candidateI] = pointi;
                    candidateDistSqr[candidateI] = sqr(distances[levelI]);
                    candidateI++;
                    break;
                }
            }
        }
        candidateMap.setSize(candidateI);
        candidateDistSqr.setSize(candidateI);

        // Do the expensive nearest test only for the candidate points.
        List<pointIndexHit> nearInfo;
        allGeometry_[shells_[shellI]].findNearest
        (
            pointField(pt, candidateMap),
            candidateDistSqr,
            nearInfo
        );

        // Update maxLevel
        forAll(nearInfo, i)
        {
            if (nearInfo[i].hit())
            {
                const label pointi = candidateMap[i];

                // Check which level it actually is in.
                const label minDistI = findLower
                (
                    distances,
                    nearInfo[i].point().dist(pt[pointi])
                );

                // pt is inbetween shell[minDistI] and shell[minDistI+1]
                maxLevel[pointi] = levels[minDistI+1];
            }
        }
    }
    else
    {
        // Inside/outside mode

        // Collect all those points that have a current maxLevel less than the
        // shell.

        pointField candidates(pt.size());
        labelList candidateMap(pt.size());
        label candidateI = 0;

        forAll(maxLevel, pointi)
        {
            if (levels[0] > maxLevel[pointi])
            {
                candidates[candidateI] = pt[pointi];
                candidateMap[candidateI] = pointi;
                candidateI++;
            }
        }
        candidates.setSize(candidateI);
        candidateMap.setSize(candidateI);

        // Do the expensive nearest test only for the candidate points.
        List<volumeType> volType;
        allGeometry_[shells_[shellI]].getVolumeType(candidates, volType);

        forAll(volType, i)
        {
            label pointi = candidateMap[i];

            if
            (
                (
                    modes_[shellI] == INSIDE
                 && volType[i] == volumeType::INSIDE
                )
             || (
                    modes_[shellI] == OUTSIDE
                 && volType[i] == volumeType::OUTSIDE
                )
            )
            {
                maxLevel[pointi] = levels[0];
            }
        }
    }
}


void Foam::shellSurfaces::findHigherGapLevel
(
    const pointField& pt,
    const labelList& ptLevel,
    const label shellI,
    labelList& gapShell,
    List<FixedList<label, 3>>& gapInfo,
    List<volumeType>& gapMode
) const
{
    //TBD: hardcoded for region 0 information
    const FixedList<label, 3>& info = extendedGapLevel_[shellI][0];
    const volumeType mode = extendedGapMode_[shellI][0];

    if (info[2] == 0)
    {
        return;
    }

    if (modes_[shellI] == DISTANCE)
    {
        // Distance mode.
        const scalar distance = distances_[shellI][0];

        labelList candidateMap(pt.size());
        scalarField candidateDistSqr(pt.size());
        label candidateI = 0;

        forAll(ptLevel, pointI)
        {
            if (ptLevel[pointI] >= info[1] && ptLevel[pointI] < info[2])
            {
                candidateMap[candidateI] = pointI;
                candidateDistSqr[candidateI] = sqr(distance);
                candidateI++;
            }
        }
        candidateMap.setSize(candidateI);
        candidateDistSqr.setSize(candidateI);

        // Do the expensive nearest test only for the candidate points.
        List<pointIndexHit> nearInfo;
        allGeometry_[shells_[shellI]].findNearest
        (
            pointField(pt, candidateMap),
            candidateDistSqr,
            nearInfo
        );

        // Update info
        forAll(nearInfo, i)
        {
            if (nearInfo[i].hit())
            {
                const label pointI = candidateMap[i];
                gapShell[pointI] = shellI;
                gapInfo[pointI] = info;
                gapMode[pointI] = mode;
            }
        }
    }
    else
    {
        // Collect all those points that have a current maxLevel less than the
        // shell.

        labelList candidateMap(pt.size());
        label candidateI = 0;

        forAll(ptLevel, pointI)
        {
            if (ptLevel[pointI] >= info[1] && ptLevel[pointI] < info[2])
            {
                candidateMap[candidateI++] = pointI;
            }
        }
        candidateMap.setSize(candidateI);

        // Do the expensive nearest test only for the candidate points.
        List<volumeType> volType;
        allGeometry_[shells_[shellI]].getVolumeType
        (
            pointField(pt, candidateMap),
            volType
        );

        forAll(volType, i)
        {
            const label pointI = candidateMap[i];
            const bool isInside = (volType[i] == volumeType::INSIDE);

            if
            (
                (
                    (modes_[shellI] == INSIDE && isInside)
                 || (modes_[shellI] == OUTSIDE && !isInside)
                )
             && info[2] > gapInfo[pointI][2]
            )
            {
                gapShell[pointI] = shellI;
                gapInfo[pointI] = info;
                gapMode[pointI] = mode;
            }
        }
    }
}


void Foam::shellSurfaces::findLevel
(
    const pointField& pt,
    const label shellI,
    labelList& minLevel,
    labelList& shell
) const
{
    const labelList& levels = levels_[shellI];

    if (modes_[shellI] == DISTANCE)
    {
        // Distance mode.

        const scalarField& distances = distances_[shellI];

        // Collect all those points that have a current level equal/greater
        // (any of) the shell. Also collect the furthest distance allowable
        // to any shell with a higher level.

        pointField candidates(pt.size());
        labelList candidateMap(pt.size());
        scalarField candidateDistSqr(pt.size());
        label candidateI = 0;

        forAll(shell, pointI)
        {
            if (shell[pointI] == -1)
            {
                forAllReverse(levels, levelI)
                {
                    if (levels[levelI] <= minLevel[pointI])
                    {
                        candidates[candidateI] = pt[pointI];
                        candidateMap[candidateI] = pointI;
                        candidateDistSqr[candidateI] = sqr(distances[levelI]);
                        candidateI++;
                        break;
                    }
                }
            }
        }
        candidates.setSize(candidateI);
        candidateMap.setSize(candidateI);
        candidateDistSqr.setSize(candidateI);

        // Do the expensive nearest test only for the candidate points.
        List<pointIndexHit> nearInfo;
        allGeometry_[shells_[shellI]].findNearest
        (
            candidates,
            candidateDistSqr,
            nearInfo
        );

        // Update maxLevel
        forAll(nearInfo, i)
        {
            if (nearInfo[i].hit())
            {
                // Check which level it actually is in.
                label minDistI = findLower
                (
                    distances,
                    nearInfo[i].point().dist(candidates[i])
                );

                label pointI = candidateMap[i];

                // pt is inbetween shell[minDistI] and shell[minDistI+1]
                shell[pointI] = shellI;
                minLevel[pointI] = levels[minDistI+1];
            }
        }
    }
    else
    {
        // Inside/outside mode

        // Collect all those points that have a current maxLevel less than the
        // shell.

        pointField candidates(pt.size());
        labelList candidateMap(pt.size());
        label candidateI = 0;

        forAll(shell, pointI)
        {
            if (shell[pointI] == -1 && levels[0] <= minLevel[pointI])
            {
                candidates[candidateI] = pt[pointI];
                candidateMap[candidateI] = pointI;
                candidateI++;
            }
        }
        candidates.setSize(candidateI);
        candidateMap.setSize(candidateI);

        // Do the expensive nearest test only for the candidate points.
        List<volumeType> volType;
        allGeometry_[shells_[shellI]].getVolumeType(candidates, volType);

        forAll(volType, i)
        {
            if
            (
                (
                    modes_[shellI] == INSIDE
                 && volType[i] == volumeType::INSIDE
                )
             || (
                    modes_[shellI] == OUTSIDE
                 && volType[i] == volumeType::OUTSIDE
                )
            )
            {
                label pointI = candidateMap[i];
                shell[pointI] = shellI;
                minLevel[pointI] = levels[0];
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::shellSurfaces::shellSurfaces
(
    const searchableSurfaces& allGeometry,
    const dictionary& shellsDict,
    const bool dryRun
)
:
    allGeometry_(allGeometry),
    dryRun_(dryRun)
{
    // Wildcard specification : loop over all surfaces and try to find a match.

    // Count number of shells.
    label shellI = 0;
    for (const word& geomName : allGeometry_.names())
    {
        if (shellsDict.found(geomName))
        {
            ++shellI;
        }
    }


    // Size lists
    shells_.setSize(shellI);
    modes_.setSize(shellI);
    distances_.setSize(shellI);
    levels_.setSize(shellI);
    dirLevels_.setSize(shellI);
    smoothDirection_.setSize(shellI);
    nSmoothExpansion_.setSize(shellI);
    nSmoothPosition_.setSize(shellI);

    extendedGapLevel_.setSize(shellI);
    extendedGapMode_.setSize(shellI);
    selfProximity_.setSize(shellI);

    const FixedList<label, 3> nullGapLevel({0, 0, 0});

    wordHashSet unmatchedKeys(shellsDict.toc());
    shellI = 0;

    forAll(allGeometry_.names(), geomI)
    {
        const word& geomName = allGeometry_.names()[geomI];

        const entry* eptr = shellsDict.findEntry(geomName, keyType::REGEX);

        if (eptr)
        {
            const dictionary& dict = eptr->dict();
            unmatchedKeys.erase(eptr->keyword());

            shells_[shellI] = geomI;
            modes_[shellI] = refineModeNames_.get("mode", dict);

            // Read pairs of distance+level
            setAndCheckLevels(shellI, dict.lookup("levels"));


            // Directional refinement
            // ~~~~~~~~~~~~~~~~~~~~~~

            dirLevels_[shellI] = Tuple2<labelPair,labelVector>
            (
                labelPair(labelMax, labelMin),
                labelVector::zero
            );
            const entry* levelPtr =
                dict.findEntry("levelIncrement", keyType::REGEX);

            if (levelPtr)
            {
                // Do reading ourselves since using labelPair would require
                // additional bracket pair
                ITstream& is = levelPtr->stream();

                is.readBegin("levelIncrement");
                is  >> dirLevels_[shellI].first().first()
                    >> dirLevels_[shellI].first().second()
                    >> dirLevels_[shellI].second();
                is.readEnd("levelIncrement");

                if (modes_[shellI] == INSIDE)
                {
                    if (!dryRun_)
                    {
                        Info<< "Additional directional refinement level"
                            << " for all cells inside " << geomName << endl;
                    }
                }
                else if (modes_[shellI] == OUTSIDE)
                {
                    if (!dryRun_)
                    {
                        Info<< "Additional directional refinement level"
                            << " for all cells outside " << geomName << endl;
                    }
                }
                else
                {
                    FatalIOErrorInFunction(shellsDict)
                        << "Unsupported mode "
                        << refineModeNames_[modes_[shellI]]
                        << exit(FatalIOError);
                }
            }

            // Directional smoothing
            // ~~~~~~~~~~~~~~~~~~~~~

            nSmoothExpansion_[shellI] = 0;
            nSmoothPosition_[shellI] = 0;
            smoothDirection_[shellI] =
                dict.getOrDefault("smoothDirection", vector::zero);

            if (smoothDirection_[shellI] != vector::zero)
            {
                dict.readEntry("nSmoothExpansion", nSmoothExpansion_[shellI]);
                dict.readEntry("nSmoothPosition", nSmoothPosition_[shellI]);
            }


            // Gap specification
            // ~~~~~~~~~~~~~~~~~


            // Shell-wide gap level specification
            const searchableSurface& surface = allGeometry_[geomI];
            const wordList& regionNames = surface.regions();

            FixedList<label, 3> gapSpec
            (
                dict.getOrDefault
                (
                    "gapLevel",
                    nullGapLevel
                )
            );
            extendedGapLevel_[shellI].setSize(regionNames.size());
            extendedGapLevel_[shellI] = gapSpec;

            extendedGapMode_[shellI].setSize(regionNames.size());
            extendedGapMode_[shellI] =
                volumeType("gapMode", dict, volumeType::MIXED);

            // Detect self-intersections
            selfProximity_[shellI].setSize
            (
                regionNames.size(),
                dict.getOrDefault<bool>("gapSelf", true)
            );


            // Override on a per-region basis?

            if (dict.found("regions"))
            {
                const dictionary& regionsDict = dict.subDict("regions");
                forAll(regionNames, regionI)
                {
                    if (regionsDict.found(regionNames[regionI]))
                    {
                        // Get the dictionary for region
                        const dictionary& regionDict = regionsDict.subDict
                        (
                            regionNames[regionI]
                        );
                        FixedList<label, 3> gapSpec
                        (
                            regionDict.getOrDefault
                            (
                                "gapLevel",
                                nullGapLevel
                            )
                        );
                        extendedGapLevel_[shellI][regionI] = gapSpec;

                        extendedGapMode_[shellI][regionI] =
                            volumeType
                            (
                                "gapMode",
                                regionDict,
                                volumeType::MIXED
                            );

                        selfProximity_[shellI][regionI] =
                            regionDict.getOrDefault<bool>
                            (
                                "gapSelf",
                                true
                            );
                    }
                }
            }

            if (extendedGapLevel_[shellI][0][0] > 0)
            {
                Info<< "Refinement level up to "
                    << extendedGapLevel_[shellI][0][2]
                    << " for all cells in gaps for shell "
                    << geomName << endl;

                if (distances_[shellI].size() > 1)
                {
                    WarningInFunction << "Using first distance only out of "
                        << distances_[shellI] << " to detect gaps for shell "
                        << geomName << endl;
                }
            }

            shellI++;
        }
    }

    if (unmatchedKeys.size() > 0)
    {
        IOWarningInFunction(shellsDict)
            << "Not all entries in refinementRegions dictionary were used."
            << " The following entries were not used : "
            << unmatchedKeys.sortedToc()
            << endl;
    }

    // Orient shell surfaces before any searching is done. Note that this
    // only needs to be done for inside or outside. Orienting surfaces
    // constructs lots of addressing which we want to avoid.
    orient();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Highest shell level
Foam::label Foam::shellSurfaces::maxLevel() const
{
    label overallMax = 0;
    forAll(levels_, shellI)
    {
        overallMax = max(overallMax, max(levels_[shellI]));
    }
    return overallMax;
}


Foam::labelList Foam::shellSurfaces::maxGapLevel() const
{
    labelList surfaceMax(extendedGapLevel_.size(), Zero);

    forAll(extendedGapLevel_, shelli)
    {
        const List<FixedList<label, 3>>& levels = extendedGapLevel_[shelli];
        forAll(levels, i)
        {
            surfaceMax[shelli] = max(surfaceMax[shelli], levels[i][2]);
        }
    }
    return surfaceMax;
}


Foam::labelPairList Foam::shellSurfaces::directionalSelectLevel() const
{
    labelPairList levels(dirLevels_.size());
    forAll(dirLevels_, shelli)
    {
        levels[shelli] = dirLevels_[shelli].first();
    }
    return levels;
}


const Foam::labelList& Foam::shellSurfaces::nSmoothExpansion() const
{
    return nSmoothExpansion_;
}


const Foam::vectorField& Foam::shellSurfaces::smoothDirection() const
{
    return smoothDirection_;
}


const Foam::labelList& Foam::shellSurfaces::nSmoothPosition() const
{
    return nSmoothPosition_;
}


void Foam::shellSurfaces::findHigherLevel
(
    const pointField& pt,
    const labelList& ptLevel,
    labelList& maxLevel
) const
{
    // Maximum level of any shell. Start off with level of point.
    maxLevel = ptLevel;

    forAll(shells_, shelli)
    {
        findHigherLevel(pt, shelli, maxLevel);
    }
}


void Foam::shellSurfaces::findHigherGapLevel
(
    const pointField& pt,
    const labelList& ptLevel,
    labelList& gapShell,
    List<FixedList<label, 3>>& gapInfo,
    List<volumeType>& gapMode
) const
{
    gapShell.setSize(pt.size());
    gapShell = -1;

    gapInfo.setSize(pt.size());
    gapInfo = FixedList<label, 3>({0, 0, 0});

    gapMode.setSize(pt.size());
    gapMode = volumeType::MIXED;

    forAll(shells_, shelli)
    {
        findHigherGapLevel(pt, ptLevel, shelli, gapShell, gapInfo, gapMode);
    }
}


void Foam::shellSurfaces::findHigherGapLevel
(
    const pointField& pt,
    const labelList& ptLevel,
    List<FixedList<label, 3>>& gapInfo,
    List<volumeType>& gapMode
) const
{
    labelList gapShell;
    findHigherGapLevel(pt, ptLevel, gapShell, gapInfo, gapMode);
}


void Foam::shellSurfaces::findLevel
(
    const pointField& pt,
    const labelList& ptLevel,
    labelList& shell
) const
{
    shell.setSize(pt.size());
    shell = -1;

    labelList minLevel(ptLevel);

    forAll(shells_, shelli)
    {
        findLevel(pt, shelli, minLevel, shell);
    }
}


void Foam::shellSurfaces::findDirectionalLevel
(
    const pointField& pt,
    const labelList& ptLevel,
    const labelList& dirLevel,  // directional level
    const direction dir,
    labelList& shell
) const
{
    shell.setSize(pt.size());
    shell = -1;

    List<volumeType> volType;

    // Current back to original
    DynamicList<label> candidateMap(pt.size());

    forAll(shells_, shelli)
    {
        if (modes_[shelli] == INSIDE || modes_[shelli] == OUTSIDE)
        {
            const labelPair& selectLevels = dirLevels_[shelli].first();
            const label addLevel = dirLevels_[shelli].second()[dir];

            // Collect the cells that are of the right original level
            candidateMap.clear();
            forAll(ptLevel, celli)
            {
                label level = ptLevel[celli];

                if
                (
                    level >= selectLevels.first()
                 && level <= selectLevels.second()
                 && dirLevel[celli] < level+addLevel
                )
                {
                    candidateMap.append(celli);
                }
            }

            // Do geometric test
            pointField candidatePt(pt, candidateMap);
            allGeometry_[shells_[shelli]].getVolumeType(candidatePt, volType);

            // Extract selected cells
            forAll(candidateMap, i)
            {
                if
                (
                    (
                        modes_[shelli] == INSIDE
                     && volType[i] == volumeType::INSIDE
                    )
                 || (
                        modes_[shelli] == OUTSIDE
                     && volType[i] == volumeType::OUTSIDE
                    )
                )
                {
                    shell[candidateMap[i]] = shelli;
                }
            }
        }
    }
}


// ************************************************************************* //
