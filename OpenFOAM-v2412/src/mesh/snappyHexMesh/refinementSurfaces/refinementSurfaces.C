/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2015-2022,2024 OpenCFD Ltd.
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

#include "refinementSurfaces.H"
#include "Time.H"
#include "searchableSurfaces.H"
#include "shellSurfaces.H"
#include "triSurfaceMesh.H"
#include "labelPair.H"
#include "searchableSurfacesQueries.H"
#include "UPtrList.H"
#include "volumeType.H"
// For dictionary::get wrapper
#include "meshRefinement.H"

#include "OBJstream.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::labelList Foam::refinementSurfaces::findHigherLevel
(
    const searchableSurface& geom,
    const shellSurfaces& shells,
    const List<pointIndexHit>& intersectionInfo,
    const labelList& surfaceLevel       // current level
) const
{
    // See if a cached level field available
    labelList minLevelField;
    geom.getField(intersectionInfo, minLevelField);


    // Detect any uncached values and do proper search
    labelList localLevel(surfaceLevel);
    {
        // Check hits:
        // 1. cached value == -1 : store for re-testing
        // 2. cached value != -1 : use
        // 3. uncached : use region 0 value

        DynamicList<label> retestSet;
        // label nHits = 0;

        forAll(intersectionInfo, i)
        {
            if (intersectionInfo[i].hit())
            {
                // ++nHits;

                // Check if minLevelField for this surface.
                if (minLevelField.size())
                {
                    if (minLevelField[i] == -1)
                    {
                        retestSet.append(i);
                    }
                    else
                    {
                        localLevel[i] = max(localLevel[i], minLevelField[i]);
                    }
                }
                else
                {
                    retestSet.append(i);
                }
            }
        }

        if (returnReduceOr(retestSet.size()))
        {
            //Info<< "Retesting "
            //    << returnReduce(retestSet.size(), sumOp<label>())
            //    << " out of " << returnReduce(nHits, sumOp<label>())
            //    << " intersections on uncached elements on geometry "
            //    << geom.name() << endl;

            pointField samples(retestSet.size());
            forAll(retestSet, i)
            {
                samples[i] = intersectionInfo[retestSet[i]].hitPoint();
            }
            labelList shellLevel;
            shells.findHigherLevel
            (
                samples,
                labelUIndList(surfaceLevel, retestSet)(),
                shellLevel
            );
            forAll(retestSet, i)
            {
                label sampleI = retestSet[i];
                localLevel[sampleI] = max(localLevel[sampleI], shellLevel[i]);
            }
        }
    }

    return localLevel;
}


Foam::labelList Foam::refinementSurfaces::calcSurfaceIndex
(
    const searchableSurfaces& allGeometry,
    const labelList& surfaces
)
{
    // Determine overall number of global regions
    label globalI = 0;
    forAll(surfaces, surfI)
    {
        globalI += allGeometry[surfaces[surfI]].regions().size();
    }

    labelList regionToSurface(globalI);
    globalI = 0;
    forAll(surfaces, surfI)
    {
        const label nLocal = allGeometry[surfaces[surfI]].regions().size();
        for (label i = 0; i < nLocal; i++)
        {
            regionToSurface[globalI++] = surfI;
        }
    }

    return regionToSurface;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::refinementSurfaces::refinementSurfaces
(
    const searchableSurfaces& allGeometry,
    const dictionary& surfacesDict,
    const label gapLevelIncrement,
    const bool dryRun
)
:
    allGeometry_(allGeometry),
    surfaces_(surfacesDict.size()),
    names_(surfacesDict.size()),
    surfZones_(surfacesDict.size()),
    regionOffset_(surfacesDict.size()),
    dryRun_(dryRun)
{
    // Wildcard specification : loop over all surface, all regions
    // and try to find a match.

    // Count number of surfaces.
    label surfI = 0;
    forAll(allGeometry_.names(), geomI)
    {
        const word& geomName = allGeometry_.names()[geomI];

        if (surfacesDict.found(geomName))
        {
            surfI++;
        }
    }

    // Size lists
    surfaces_.setSize(surfI);
    names_.setSize(surfI);
    surfZones_.setSize(surfI);
    regionOffset_.setSize(surfI);

    // Per surface data
    labelList globalMinLevel(surfI, Zero);
    labelList globalMaxLevel(surfI, Zero);
    labelList globalLevelIncr(surfI, Zero);

    const FixedList<label, 3> nullGapLevel({0, 0, 0});
    List<FixedList<label, 3>> globalGapLevel(surfI);
    List<volumeType> globalGapMode(surfI);
    boolList globalGapSelf(surfI);

    const FixedList<label, 4> nullCurvLevel({0, 0, 0, -1});
    List<FixedList<label, 4>> globalCurvLevel(surfI);

    scalarField globalAngle(surfI, -GREAT);
    PtrList<dictionary> globalPatchInfo(surfI);

    labelList globalBlockLevel(surfI, labelMax);
    labelList globalLeakLevel(surfI, labelMax);

    // Supported in buffer-layer mode only
    boolList globalBufferLayers(surfI, true);

    // Per surface, per region data
    List<Map<label>> regionMinLevel(surfI);
    List<Map<label>> regionMaxLevel(surfI);
    List<Map<label>> regionLevelIncr(surfI);
    List<Map<FixedList<label, 3>>> regionGapLevel(surfI);
    List<Map<volumeType>> regionGapMode(surfI);
    List<Map<bool>> regionGapSelf(surfI);
    List<Map<FixedList<label, 4>>> regionCurvLevel(surfI);
    List<Map<scalar>> regionAngle(surfI);
    List<Map<autoPtr<dictionary>>> regionPatchInfo(surfI);
    List<Map<label>> regionBlockLevel(surfI);
    List<Map<label>> regionLeakLevel(surfI);
    List<Map<label>> regionBufferLayers(surfI);

    wordHashSet unmatchedKeys(surfacesDict.toc());

    surfI = 0;
    forAll(allGeometry_.names(), geomI)
    {
        const word& geomName = allGeometry_.names()[geomI];

        const entry* ePtr =
            surfacesDict.findEntry(geomName, keyType::REGEX);

        if (ePtr)
        {
            const dictionary& dict = ePtr->dict();
            unmatchedKeys.erase(ePtr->keyword());

            names_[surfI] = geomName;
            surfaces_[surfI] = geomI;

            const labelPair refLevel
            (
                meshRefinement::get<labelPair>
                (
                    dict,
                    "level",
                    dryRun_,
                    keyType::REGEX,
                    labelPair(0, 0)
                )
            );

            globalMinLevel[surfI] = refLevel[0];
            globalMaxLevel[surfI] = refLevel[1];
            globalLevelIncr[surfI] = dict.getOrDefault
            (
                "gapLevelIncrement",
                gapLevelIncrement
            );

            if
            (
                globalMinLevel[surfI] < 0
             || globalMaxLevel[surfI] < globalMinLevel[surfI]
             || globalMaxLevel[surfI] < 0
             || globalLevelIncr[surfI] < 0
            )
            {
                FatalIOErrorInFunction(dict)
                    << "Illegal level specification for surface "
                    << names_[surfI]
                    << " : minLevel:" << globalMinLevel[surfI]
                    << " maxLevel:" << globalMaxLevel[surfI]
                    << " levelIncrement:" << globalLevelIncr[surfI]
                    << exit(FatalIOError);
            }


            // Optional gapLevel specification

            globalGapLevel[surfI] = dict.getOrDefault
            (
                "gapLevel",
                nullGapLevel
            );
            globalGapMode[surfI] =
                volumeType("gapMode", dict, volumeType::MIXED);

            if
            (
                globalGapMode[surfI] == volumeType::UNKNOWN
             || globalGapLevel[surfI][0] < 0
             || globalGapLevel[surfI][1] < 0
             || globalGapLevel[surfI][2] < 0
             || globalGapLevel[surfI][1] > globalGapLevel[surfI][2]
            )
            {
                FatalIOErrorInFunction(dict)
                    << "Illegal gapLevel specification for surface "
                    << names_[surfI]
                    << " : gapLevel:" << globalGapLevel[surfI]
                    << " gapMode:" << globalGapMode[surfI].str()
                    << exit(FatalIOError);
            }

            globalGapSelf[surfI] =
                dict.getOrDefault<bool>("gapSelf", true);

            globalCurvLevel[surfI] = nullCurvLevel;
            if (dict.readIfPresent("curvatureLevel", globalCurvLevel[surfI]))
            {
                if (globalCurvLevel[surfI][0] <= 0)
                {
                    FatalIOErrorInFunction(dict)
                        << "Illegal curvatureLevel specification for surface "
                        << names_[surfI]
                        << " : curvatureLevel:" << globalCurvLevel[surfI]
                        << exit(FatalIOError);
                }
            }

            const searchableSurface& surface = allGeometry_[surfaces_[surfI]];

            // Surface zones
            surfZones_.set
            (
                surfI,
                new surfaceZonesInfo
                (
                    surface,
                    dict,
                    allGeometry_.regionNames()[surfaces_[surfI]]
                )
            );

            // Global perpendicular angle
            if (dict.found("patchInfo"))
            {
                globalPatchInfo.set
                (
                    surfI,
                    dict.subDict("patchInfo").clone()
                );
            }
            dict.readIfPresent("perpendicularAngle", globalAngle[surfI]);
            dict.readIfPresent("blockLevel", globalBlockLevel[surfI]);
            dict.readIfPresent("leakLevel", globalLeakLevel[surfI]);
            dict.readIfPresent("addBufferLayers", globalBufferLayers[surfI]);

            if (dict.found("regions"))
            {
                const dictionary& regionsDict = dict.subDict("regions");
                const wordList& regionNames = surface.regions();

                forAll(regionNames, regionI)
                {
                    if (regionsDict.found(regionNames[regionI]))
                    {
                        // Get the dictionary for region
                        const dictionary& regionDict = regionsDict.subDict
                        (
                            regionNames[regionI]
                        );

                        const labelPair refLevel
                        (
                            meshRefinement::get<labelPair>
                            (
                                regionDict,
                                "level",
                                dryRun_,
                                keyType::REGEX,
                                //labelPair(0, 0)
                                labelPair
                                (
                                    globalMinLevel[surfI],
                                    globalMaxLevel[surfI]
                                )
                            )
                        );


                        regionMinLevel[surfI].insert(regionI, refLevel[0]);
                        regionMaxLevel[surfI].insert(regionI, refLevel[1]);
                        label levelIncr = regionDict.getOrDefault
                        (
                            "gapLevelIncrement",
                            gapLevelIncrement
                        );
                        regionLevelIncr[surfI].insert(regionI, levelIncr);

                        if
                        (
                            refLevel[0] < 0
                         || refLevel[1] < refLevel[0]
                         || levelIncr < 0
                        )
                        {
                            FatalIOErrorInFunction(dict)
                                << "Illegal level specification for surface "
                                << names_[surfI] << " region "
                                << regionNames[regionI]
                                << " : minLevel:" << refLevel[0]
                                << " maxLevel:" << refLevel[1]
                                << " levelIncrement:" << levelIncr
                                << exit(FatalIOError);
                        }



                        // Optional gapLevel specification

                        FixedList<label, 3> gapSpec
                        (
                            regionDict.getOrDefault
                            (
                                "gapLevel",
                                globalGapLevel[surfI]   //nullGapLevel
                            )
                        );
                        regionGapLevel[surfI].insert(regionI, gapSpec);
                        volumeType gapModeSpec
                        (
                            "gapMode",
                            regionDict,
                            globalGapMode[surfI]   //volumeType::MIXED
                        );
                        regionGapMode[surfI].insert(regionI, gapModeSpec);
                        if
                        (
                            gapModeSpec == volumeType::UNKNOWN
                         || gapSpec[0] < 0
                         || gapSpec[1] < 0
                         || gapSpec[2] < 0
                         || gapSpec[1] > gapSpec[2]
                        )
                        {
                            FatalIOErrorInFunction(dict)
                                << "Illegal gapLevel specification for surface "
                                << names_[surfI]
                                << " : gapLevel:" << gapSpec
                                << " gapMode:" << gapModeSpec.str()
                                << exit(FatalIOError);
                        }
                        regionGapSelf[surfI].insert
                        (
                            regionI,
                            regionDict.getOrDefault<bool>
                            (
                                "gapSelf",
                                globalGapSelf[surfI]    //true
                            )
                        );

                        // Start off with the per-surface level
                        FixedList<label, 4> curvSpec(nullCurvLevel);
                        if
                        (
                            regionDict.readIfPresent("curvatureLevel", curvSpec)
                        )
                        {
                            if (curvSpec[0] <= 0)
                            {
                                FatalIOErrorInFunction(dict)
                                    << "Illegal curvatureLevel"
                                    << " specification for surface "
                                    << names_[surfI]
                                    << " : curvatureLevel:" << curvSpec
                                    << exit(FatalIOError);
                            }
                            regionCurvLevel[surfI].insert(regionI, curvSpec);
                        }

                        if (regionDict.found("perpendicularAngle"))
                        {
                            regionAngle[surfI].insert
                            (
                                regionI,
                                regionDict.get<scalar>("perpendicularAngle")
                            );
                        }

                        if (regionDict.found("patchInfo"))
                        {
                            regionPatchInfo[surfI].insert
                            (
                                regionI,
                                regionDict.subDict("patchInfo").clone()
                            );
                        }

                        label l;
                        if (regionDict.readIfPresent<label>("blockLevel", l))
                        {
                            regionBlockLevel[surfI].insert(regionI, l);
                        }
                        if (regionDict.readIfPresent<label>("leakLevel", l))
                        {
                            regionLeakLevel[surfI].insert(regionI, l);
                        }
                        bool s;
                        if
                        (
                            regionDict.readIfPresent<bool>("addBufferLayers", s)
                        )
                        {
                            regionBufferLayers[surfI].insert(regionI, s);
                        }
                    }
                }
            }
            surfI++;
        }
    }

    if (unmatchedKeys.size() > 0)
    {
        IOWarningInFunction(surfacesDict)
            << "Not all entries in refinementSurfaces dictionary were used."
            << " The following entries were not used : "
            << unmatchedKeys.sortedToc()
            << endl;
    }


    // Calculate local to global region offset
    label nRegions = 0;

    forAll(surfaces_, surfI)
    {
        regionOffset_[surfI] = nRegions;
        nRegions += allGeometry_[surfaces_[surfI]].regions().size();
    }

    // Rework surface specific information into information per global region

    regionToSurface_ = calcSurfaceIndex(allGeometry_, surfaces_);


    minLevel_.setSize(nRegions);
    minLevel_ = 0;
    maxLevel_.setSize(nRegions);
    maxLevel_ = 0;
    gapLevel_.setSize(nRegions);
    gapLevel_ = -1;
    extendedGapLevel_.setSize(nRegions);
    extendedGapLevel_ = nullGapLevel;
    extendedGapMode_.setSize(nRegions);
    extendedGapMode_ = volumeType::UNKNOWN;
    selfProximity_.setSize(nRegions);
    selfProximity_ = true;
    extendedCurvatureLevel_.setSize(nRegions);
    extendedCurvatureLevel_ = nullCurvLevel;
    perpendicularAngle_.setSize(nRegions);
    perpendicularAngle_ = -GREAT;
    patchInfo_.setSize(nRegions);
    blockLevel_.setSize(nRegions);
    blockLevel_ = labelMax;
    leakLevel_.setSize(nRegions);
    leakLevel_ = labelMax;
    addBufferLayers_.setSize(nRegions);
    addBufferLayers_ = false;


    forAll(globalMinLevel, surfI)
    {
        label nRegions = allGeometry_[surfaces_[surfI]].regions().size();

        // Initialise to global (i.e. per surface)
        for (label i = 0; i < nRegions; i++)
        {
            const label globalRegionI = regionOffset_[surfI] + i;

            minLevel_[globalRegionI] = globalMinLevel[surfI];
            maxLevel_[globalRegionI] = globalMaxLevel[surfI];
            gapLevel_[globalRegionI] =
                maxLevel_[globalRegionI]
              + globalLevelIncr[surfI];
            extendedGapLevel_[globalRegionI] = globalGapLevel[surfI];
            extendedGapMode_[globalRegionI] = globalGapMode[surfI];
            selfProximity_[globalRegionI] = globalGapSelf[surfI];
            extendedCurvatureLevel_[globalRegionI] = globalCurvLevel[surfI];
            perpendicularAngle_[globalRegionI] = globalAngle[surfI];
            if (globalPatchInfo.set(surfI))
            {
                patchInfo_.set
                (
                    globalRegionI,
                    globalPatchInfo[surfI].clone()
                );
            }
            blockLevel_[globalRegionI] = globalBlockLevel[surfI];
            leakLevel_[globalRegionI] = globalLeakLevel[surfI];
            addBufferLayers_[globalRegionI] = globalBufferLayers[surfI];
        }

        // Overwrite with region specific information
        forAllConstIters(regionMinLevel[surfI], iter)
        {
            const label globalRegionI = regionOffset_[surfI] + iter.key();

            minLevel_[globalRegionI] = iter.val();
            maxLevel_[globalRegionI] = regionMaxLevel[surfI][iter.key()];
            gapLevel_[globalRegionI] =
                maxLevel_[globalRegionI]
              + regionLevelIncr[surfI][iter.key()];
        }

        forAllConstIters(regionGapLevel[surfI], iter)
        {
            const label globalRegionI = regionOffset_[surfI] + iter.key();

            extendedGapLevel_[globalRegionI] =
                regionGapLevel[surfI][iter.key()];
            extendedGapMode_[globalRegionI] =
                regionGapMode[surfI][iter.key()];
            selfProximity_[globalRegionI] =
                regionGapSelf[surfI][iter.key()];
        }

        forAllConstIters(regionCurvLevel[surfI], iter)
        {
            const label globalRegionI = regionOffset_[surfI] + iter.key();
            extendedCurvatureLevel_[globalRegionI] = iter.val();
        }

        forAllConstIters(regionAngle[surfI], iter)
        {
            const label globalRegionI = regionOffset_[surfI] + iter.key();
            perpendicularAngle_[globalRegionI] = iter.val();
        }

        const Map<autoPtr<dictionary>>& localInfo = regionPatchInfo[surfI];
        forAllConstIters(localInfo, iter)
        {
            const label globalRegionI = regionOffset_[surfI] + iter.key();
            const dictionary& dict = *(iter.val());

            patchInfo_.set(globalRegionI, dict.clone());
        }

        forAllConstIters(regionBlockLevel[surfI], iter)
        {
            const label globalRegionI = regionOffset_[surfI] + iter.key();

            blockLevel_[globalRegionI] = iter.val();
            leakLevel_[globalRegionI] = iter.val();
        }
        forAllConstIters(regionBufferLayers[surfI], iter)
        {
            const label globalRegionI = regionOffset_[surfI] + iter.key();
            addBufferLayers_[globalRegionI] = iter.val();
        }
    }
}


Foam::refinementSurfaces::refinementSurfaces
(
    const searchableSurfaces& allGeometry,
    const labelList& surfaces,
    const wordList& names,
    const PtrList<surfaceZonesInfo>& surfZones,
    const labelList& regionOffset,
    const labelList& minLevel,
    const labelList& maxLevel,
    const labelList& gapLevel,
    const scalarField& perpendicularAngle,
    PtrList<dictionary>& patchInfo,
    const bool dryRun
)
:
    allGeometry_(allGeometry),
    surfaces_(surfaces),
    names_(names),
    surfZones_(surfZones),
    regionOffset_(regionOffset),
    regionToSurface_(calcSurfaceIndex(allGeometry, surfaces)),
    minLevel_(minLevel),
    maxLevel_(maxLevel),
    gapLevel_(gapLevel),
    perpendicularAngle_(perpendicularAngle),
    patchInfo_(patchInfo.size()),
    dryRun_(dryRun)
{
    forAll(patchInfo_, pI)
    {
        if (patchInfo.set(pI))
        {
            patchInfo_.set(pI, patchInfo.set(pI, nullptr));
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelPair Foam::refinementSurfaces::whichSurface
(
    const label globalRegionI
) const
{
    const label surfI = regionToSurface_[globalRegionI];
    const label localI = globalRegionI-regionOffset_[surfI];
    return labelPair(surfI, localI);
}


// // Count number of triangles per surface region
// Foam::labelList Foam::refinementSurfaces::countRegions(const triSurface& s)
// {
//     const geometricSurfacePatchList& regions = s.patches();
//
//     labelList nTris(regions.size(), Zero);
//
//     forAll(s, triI)
//     {
//         nTris[s[triI].region()]++;
//     }
//     return nTris;
// }


// // Pre-calculate the refinement level for every element
// void Foam::refinementSurfaces::wantedRefinementLevel
// (
//     const shellSurfaces& shells,
//     const label surfI,
//     const List<pointIndexHit>& info,    // Indices
//     const pointField& ctrs,             // Representative coordinate
//     labelList& minLevelField
// ) const
// {
//     const searchableSurface& geom = allGeometry_[surfaces_[surfI]];
//
//     // Get per element the region
//     labelList region;
//     geom.getRegion(info, region);
//
//     // Initialise fields to region wise minLevel
//     minLevelField.setSize(ctrs.size());
//     minLevelField = -1;
//
//     forAll(minLevelField, i)
//     {
//         if (info[i].hit())
//         {
//             minLevelField[i] = minLevel(surfI, region[i]);
//         }
//     }
//
//     // Find out if triangle inside shell with higher level
//     // What level does shell want to refine fc to?
//     labelList shellLevel;
//     shells.findHigherLevel(ctrs, minLevelField, shellLevel);
//
//     forAll(minLevelField, i)
//     {
//         minLevelField[i] = max(minLevelField[i], shellLevel[i]);
//     }
// }


Foam::labelList Foam::refinementSurfaces::maxGapLevel() const
{
    labelList surfaceMax(surfaces_.size(), Zero);

    forAll(surfaces_, surfI)
    {
        const wordList& regionNames = allGeometry_[surfaces_[surfI]].regions();

        forAll(regionNames, regionI)
        {
            label globalI = globalRegion(surfI, regionI);
            const FixedList<label, 3>& gapInfo = extendedGapLevel_[globalI];
            surfaceMax[surfI] = max(surfaceMax[surfI], gapInfo[2]);
        }
    }
    return surfaceMax;
}


Foam::labelList Foam::refinementSurfaces::maxCurvatureLevel() const
{
    labelList surfaceMax(surfaces_.size(), Zero);

    forAll(surfaces_, surfI)
    {
        const wordList& regionNames = allGeometry_[surfaces_[surfI]].regions();

        forAll(regionNames, regionI)
        {
            label globalI = globalRegion(surfI, regionI);
            const FixedList<label, 4>& gapInfo =
                extendedCurvatureLevel_[globalI];
            surfaceMax[surfI] = max(surfaceMax[surfI], gapInfo[2]);
        }
    }
    return surfaceMax;
}


// Precalculate the refinement level for every element of the searchable
// surface.
void Foam::refinementSurfaces::setMinLevelFields(const shellSurfaces& shells)
{
    forAll(surfaces_, surfI)
    {
        const searchableSurface& geom = allGeometry_[surfaces_[surfI]];

        // Cache the refinement level (max of surface level and shell level)
        // on a per-element basis. Only makes sense if there are lots of
        // elements. Possibly should have 'enough' elements to have fine
        // enough resolution but for now just make sure we don't catch e.g.
        // searchableBox (size=6)
        if (geom.globalSize() > 10)
        {
            // Representative local coordinates and bounding sphere
            pointField ctrs;
            scalarField radiusSqr;
            geom.boundingSpheres(ctrs, radiusSqr);

            labelList minLevelField(ctrs.size(), Zero);
            {
                // Get the element index in a roundabout way. Problem is e.g.
                // distributed surface where local indices differ from global
                // ones (needed for getRegion call)
                List<pointIndexHit> info;
                geom.findNearest(ctrs, radiusSqr, info);

                // Get per element the region
                labelList region;
                geom.getRegion(info, region);

                // From the region get the surface-wise refinement level
                forAll(minLevelField, i)
                {
                    if (info[i].hit()) //Note: should not be necessary
                    {
                        minLevelField[i] = minLevel(surfI, region[i]);
                    }
                }
            }

            // Find out if triangle inside shell with higher level
            // What level does shell want to refine fc to?
            labelList shellLevel;
            shells.findHigherLevel(ctrs, minLevelField, shellLevel);


            // In case of triangulated surfaces only cache value if triangle
            // centre and vertices are in same shell
            const auto* tsPtr = isA<triSurface>(geom);

            if (tsPtr)
            {
                label nUncached = 0;

                // Check if points differing from ctr level

                const triSurface& ts = *tsPtr;
                const pointField& points = ts.points();

                // Determine minimum expected level to avoid having to
                // test lots of points
                labelList minPointLevel(points.size(), labelMax);
                forAll(shellLevel, triI)
                {
                    const labelledTri& t = ts[triI];
                    label level = shellLevel[triI];
                    forAll(t, tI)
                    {
                        minPointLevel[t[tI]] = min(minPointLevel[t[tI]], level);
                    }
                }


                // See if inside any shells with higher refinement level
                labelList pointLevel;
                shells.findHigherLevel(points, minPointLevel, pointLevel);


                // See if triangle centre values differ from triangle points
                forAll(shellLevel, triI)
                {
                    const labelledTri& t = ts[triI];
                    label fLevel = shellLevel[triI];
                    if
                    (
                        (pointLevel[t[0]] != fLevel)
                     || (pointLevel[t[1]] != fLevel)
                     || (pointLevel[t[2]] != fLevel)
                    )
                    {
                        //Pout<< "Detected triangle " << t.tri(ts.points())
                        //    << " partially inside/partially outside" << endl;

                        // Mark as uncached
                        shellLevel[triI] = -1;
                        nUncached++;
                    }
                }

                if (!dryRun_)
                {
                    Info<< "For geometry " << geom.name()
                        << " detected "
                        << returnReduce(nUncached, sumOp<label>())
                        << " uncached triangles out of " << geom.globalSize()
                        << endl;
                }
            }


            // Combine overall level field with current shell level. Make sure
            // to preserve -1 (from triSurfaceMeshes with triangles partly
            // inside/outside
            forAll(minLevelField, i)
            {
                if (min(minLevelField[i], shellLevel[i]) < 0)
                {
                    minLevelField[i] = -1;
                }
                else
                {
                    minLevelField[i] = max(minLevelField[i], shellLevel[i]);
                }
            }

            // Store minLevelField on surface
            const_cast<searchableSurface&>(geom).setField(minLevelField);
        }
    }
}


// Precalculate the refinement level for every element of the searchable
// surface.
void Foam::refinementSurfaces::setCurvatureMinLevelFields
(
    const scalar cosAngle,
    const scalar level0EdgeLength
)
{
    const labelList maxCurvLevel(maxCurvatureLevel());


    forAll(surfaces_, surfI)
    {
        // Check if there is a specification of the extended curvature
        if (maxCurvLevel[surfI] <= 0)
        {
            continue;
        }

        const searchableSurface& geom = allGeometry_[surfaces_[surfI]];

        const auto* tsPtr = isA<triSurface>(geom);

        // Cache the refinement level (max of surface level and shell level)
        // on a per-element basis. Only makes sense if there are lots of
        // elements. Possibly should have 'enough' elements to have fine
        // enough resolution but for now just make sure we don't catch e.g.
        // searchableBox (size=6)
        if (tsPtr)
        {
            // Representative local coordinates and bounding sphere
            pointField ctrs;
            scalarField radiusSqr;
            geom.boundingSpheres(ctrs, radiusSqr);

            labelList minLevelField(ctrs.size(), Zero);
            //labelList surfMin(ctrs.size(), Zero);
            labelList surfMax(ctrs.size(), Zero);
            labelList nCurvCells(ctrs.size(), Zero);
            labelList curvIgnore(ctrs.size(), -1);
            {
                // Get the element index in a roundabout way. Problem is e.g.
                // distributed surface where local indices differ from global
                // ones (needed for getRegion call)
                List<pointIndexHit> info;
                geom.findNearest(ctrs, radiusSqr, info);

                // Get per element the region
                labelList region;
                geom.getRegion(info, region);

                // See if a cached level field available (from e.g. shells)
                labelList cachedField;
                geom.getField(info, cachedField);

                // From the region get the surface-wise refinement level
                forAll(minLevelField, i)
                {
                    if (info[i].hit()) //Note: should not be necessary
                    {
                        const label globali = globalRegion(surfI, region[i]);
                        curvIgnore[i] = extendedCurvatureLevel_[globali][3];
                        nCurvCells[i] = extendedCurvatureLevel_[globali][0];
                        //surfMin[i] = extendedCurvatureLevel_[globali][1];
                        surfMax[i] = extendedCurvatureLevel_[globali][2];

                        minLevelField[i] = minLevel(surfI, region[i]);
                        if (cachedField.size() > i)
                        {
                            minLevelField[i] =
                                max(minLevelField[i], cachedField[i]);
                        }
                    }
                }
            }

            // Calculate per-triangle curvature. This is the max of the
            // measured point-based curvature + some constraints.
            scalarField cellCurv(ctrs.size(), Zero);
            {
                // Walk surface and detect sharp features. Returns maximum
                // curvature (per surface point. Note: returns per point, not
                // per localpoint)
                const auto& ts = *tsPtr;
                auto tcurv(triSurfaceTools::curvatures(ts));
                auto& curv = tcurv.ref();

                // Reset curvature on sharp edges (and neighbours since
                // curvature uses extended stencil)
                {
                    const auto& edgeFaces = ts.edgeFaces();
                    const auto& edges = ts.edges();
                    const auto& points = ts.points();
                    const auto& mp = ts.meshPoints();

                    bitSet isOnSharpEdge(points.size());
                    forAll(edgeFaces, edgei)
                    {
                        const auto& eFaces = edgeFaces[edgei];
                        const edge meshE(mp, edges[edgei]);

                        if (eFaces.size() == 2)
                        {
                            const auto& f0 = ts[eFaces[0]];
                            const auto& f1 = ts[eFaces[1]];

                            const vector n0 = f0.unitNormal(points);
                            vector n1 = f1.unitNormal(points);

                            const auto dir0 = f0.edgeDirection(meshE);
                            const auto dir1 = f1.edgeDirection(meshE);

                            if (dir0 && ((dir0 > 0) == (dir1 > 0)))
                            {
                                // Flip since use edge in same direction
                                // (should not be the case for 'proper'
                                //  surfaces)
                                n1 = -n1;
                            }

                            if ((n0&n1) < cosAngle)
                            {
                                isOnSharpEdge.set(meshE[0]);
                                isOnSharpEdge.set(meshE[1]);
                            }
                        }
                    }

                    // Extend by one layer
                    {
                        bitSet oldOnSharpEdge(isOnSharpEdge);
                        isOnSharpEdge = false;
                        for (const auto& f : ts)
                        {
                            for (const label pointi : f)
                            {
                                if (oldOnSharpEdge[pointi])
                                {
                                    // Mark all points on triangle
                                    isOnSharpEdge.set(f);
                                    break;
                                }
                            }
                        }
                    }


                    // Unmark curvature
                    autoPtr<OBJstream> str;
                    //if (debug)
                    //{
                    //    str.reset
                    //    (
                    //        new OBJstream
                    //        (
                    //            "sharpEdgePoints_"
                    //          + geom.name()
                    //          + ".obj"
                    //        )
                    //    );
                    //    Info<< "Writing sharp edge points to "
                    //        << str().name() << endl;
                    //}

                    for (const label pointi : isOnSharpEdge)
                    {
                        // Reset measured point-based curvature
                        curv[pointi] = 0.0;
                        if (str)
                        {
                            str().write(points[pointi]);
                        }
                    }
                }

                // Reset curvature on -almost- sharp edges.
                // This resets the point-based curvature if the edge
                // is considered to be a sharp edge based on its actual
                // curvature. This is only used if the 'ignore' level is
                // given.
                {
                    // Pass 1: accumulate constraints on the points - get
                    //         the minimum of curvature constraints on the
                    //         connected triangles. Looks at user-specified
                    //         min curvature - does not use actual measured
                    //         curvature
                    scalarField pointMinCurv(ts.nPoints(), VGREAT);

                    forAll(ts, i)
                    {
                        // Is ignore level given for surface
                        const label level = curvIgnore[i];
                        if (level >= 0)
                        {
                            // Convert to (inv) size
                            const scalar length = level0EdgeLength/(2<<level);
                            const scalar invLength = 1.0/length;
                            for (const label pointi : ts[i])
                            {
                                if
                                (
                                    invLength < pointMinCurv[pointi]
                                 && curv[pointi] > SMALL
                                )
                                {
                                    //Pout<< "** at location:"
                                    //    << ts.points()[pointi]
                                    //    << " measured curv:" << curv[pointi]
                                    //    << " radius:" << 1.0/curv[pointi]
                                    //    << " ignore level:" << level
                                    //    << " ignore radius:" << length
                                    //    << " resetting minCurv to "
                                    //    << invLength
                                    //    << endl;
                                }

                                pointMinCurv[pointi] =
                                    min(pointMinCurv[pointi], invLength);
                            }
                        }
                    }

                    // Clip curvature (should do nothing for most points unless
                    // ignore-level is triggered)
                    forAll(pointMinCurv, pointi)
                    {
                        if (pointMinCurv[pointi] < curv[pointi])
                        {
                            //Pout<< "** at location:" << ts.points()[pointi]
                            //    << " measured curv:" << curv[pointi]
                            //    << " radius:" << 1.0/curv[pointi]
                            //    << " cellLimit:" << pointMinCurv[pointi]
                            //    << endl;

                            // Set up to ignore point
                            //curv[pointi] = pointMinCurv[pointi];
                            curv[pointi] = 0.0;
                        }
                    }
                }


                forAll(ts, i)
                {
                    const auto& f = ts[i];
                    // Take max curvature (= min radius of curvature)
                    cellCurv[i] = max(curv[f[0]], max(curv[f[1]], curv[f[2]]));
                }
            }


            //if(debug)
            //{
            //    const scalar maxCurv = gMax(cellCurv);
            //    if (maxCurv > SMALL)
            //    {
            //        const scalar r = scalar(1.0)/maxCurv;
            //
            //        Pout<< "For geometry " << geom.name()
            //            << " have curvature max " << maxCurv
            //            << " which equates to radius:" << r
            //            << " which equates to refinement level "
            //            << log2(level0EdgeLength/r)
            //            << endl;
            //    }
            //}

            forAll(cellCurv, i)
            {
                if (cellCurv[i] > SMALL && nCurvCells[i] > 0)
                {
                    //- ?If locally have a cached field override the
                    //  surface-based level ignore any curvature?
                    //if (minLevelField[i] > surfMin[i])
                    //{
                    //    // Ignore curvature
                    //}
                    //else
                    //if (surfMin[i] == surfMax[i])
                    //{
                    //    // Ignore curvature. Bypass calculation below.
                    //}
                    //else
                    {
                        // Re-work the curvature into a radius and into a
                        // number of cells
                        const scalar r = scalar(1.0)/cellCurv[i];
                        const scalar level =
                            log2(nCurvCells[i]*level0EdgeLength/r);
                        const label l = round(level);

                        if (l > minLevelField[i] && l <= surfMax[i])
                        {
                            minLevelField[i] = l;
                        }
                    }
                }
            }


            // Store minLevelField on surface
            const_cast<searchableSurface&>(geom).setField(minLevelField);
        }
    }
}


// Find intersections of edge. Return -1 or first surface with higher minLevel
// number.
void Foam::refinementSurfaces::findHigherIntersection
(
    const shellSurfaces& shells,

    const pointField& start,
    const pointField& end,
    const labelList& currentLevel,   // current cell refinement level

    labelList& surfaces,
    labelList& surfaceLevel
) const
{
    surfaces.setSize(start.size());
    surfaces = -1;
    surfaceLevel.setSize(start.size());
    surfaceLevel = -1;

    if (surfaces_.empty())
    {
        return;
    }

    if (surfaces_.size() == 1)
    {
        // Optimisation: single segmented surface. No need to duplicate
        // point storage.

        label surfI = 0;

        const searchableSurface& geom = allGeometry_[surfaces_[surfI]];

        // Do intersection test
        List<pointIndexHit> intersectionInfo(start.size());
        geom.findLineAny(start, end, intersectionInfo);


        // Surface-based refinement level
        labelList surfaceOnlyLevel(start.size(), -1);
        {
            // Get per intersection the region
            labelList region;
            geom.getRegion(intersectionInfo, region);

            forAll(intersectionInfo, i)
            {
                if (intersectionInfo[i].hit())
                {
                    surfaceOnlyLevel[i] = minLevel(surfI, region[i]);
                }
            }
        }


        // Get shell refinement level if higher
        const labelList localLevel
        (
            findHigherLevel
            (
                geom,
                shells,
                intersectionInfo,
                surfaceOnlyLevel // starting level
            )
        );


        // Combine localLevel with current level
        forAll(localLevel, i)
        {
            if (localLevel[i] > currentLevel[i])
            {
                surfaces[i] = surfI;    // index of surface
                surfaceLevel[i] = localLevel[i];
            }
        }

        return;
    }



    // Work arrays
    pointField p0(start);
    pointField p1(end);
    labelList intersectionToPoint(identity(start.size()));
    List<pointIndexHit> intersectionInfo(start.size());

    forAll(surfaces_, surfI)
    {
        const searchableSurface& geom = allGeometry_[surfaces_[surfI]];

        // Do intersection test
        geom.findLineAny(p0, p1, intersectionInfo);


        // Surface-based refinement level
        labelList surfaceOnlyLevel(intersectionInfo.size(), -1);
        {
            // Get per intersection the region
            labelList region;
            geom.getRegion(intersectionInfo, region);

            forAll(intersectionInfo, i)
            {
                if (intersectionInfo[i].hit())
                {
                    surfaceOnlyLevel[i] = minLevel(surfI, region[i]);
                }
            }
        }


        // Get shell refinement level if higher
        const labelList localLevel
        (
            findHigherLevel
            (
                geom,
                shells,
                intersectionInfo,
                surfaceOnlyLevel
            )
        );


        // Combine localLevel with current level
        label missI = 0;
        forAll(localLevel, i)
        {
            label pointI = intersectionToPoint[i];

            if (localLevel[i] > currentLevel[pointI])
            {
                // Mark point for refinement
                surfaces[pointI] = surfI;
                surfaceLevel[pointI] = localLevel[i];
            }
            else
            {
                p0[missI] = start[pointI];
                p1[missI] = end[pointI];
                intersectionToPoint[missI] = pointI;
                missI++;
            }
        }


        // All done? Note that this decision should be synchronised
        if (returnReduceAnd(missI == 0))
        {
            break;
        }

        // Trim misses
        p0.setSize(missI);
        p1.setSize(missI);
        intersectionToPoint.setSize(missI);
        intersectionInfo.setSize(missI);
    }
}


void Foam::refinementSurfaces::findAllIntersections
(
    const pointField& start,
    const pointField& end,
    const labelList& currentLevel,   // current cell refinement level

    const labelList& globalMinLevel,
    const labelList& globalMaxLevel,

    List<vectorList>& surfaceNormal,
    labelListList& surfaceLevel
) const
{
    surfaceLevel.setSize(start.size());
    surfaceNormal.setSize(start.size());

    if (surfaces_.empty())
    {
        return;
    }

    // Work arrays
    List<List<pointIndexHit>> hitInfo;

    forAll(surfaces_, surfI)
    {
        const searchableSurface& surface = allGeometry_[surfaces_[surfI]];

        surface.findLineAll(start, end, hitInfo);

        // Repack hits for surface into flat list
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // To avoid overhead of calling getRegion for every point

        label n = 0;
        forAll(hitInfo, pointI)
        {
            n += hitInfo[pointI].size();
        }

        List<pointIndexHit> surfInfo(n);
        labelList pointMap(n);
        n = 0;

        forAll(hitInfo, pointI)
        {
            const List<pointIndexHit>& pHits = hitInfo[pointI];

            forAll(pHits, i)
            {
                surfInfo[n] = pHits[i];
                pointMap[n] = pointI;
                n++;
            }
        }

        labelList surfRegion(n);
        vectorField surfNormal(n);
        surface.getRegion(surfInfo, surfRegion);
        surface.getNormal(surfInfo, surfNormal);

        surfInfo.clear();


        // Extract back into pointwise
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~

        forAll(surfRegion, i)
        {
            label region = globalRegion(surfI, surfRegion[i]);
            label pointI = pointMap[i];

            if
            (
                currentLevel[pointI] >= globalMinLevel[region]
             && currentLevel[pointI] < globalMaxLevel[region]
            )
            {
                // Append to pointI info
                label sz = surfaceNormal[pointI].size();
                surfaceNormal[pointI].setSize(sz+1);
                surfaceNormal[pointI][sz] = surfNormal[i];

                surfaceLevel[pointI].setSize(sz+1);
                surfaceLevel[pointI][sz] = globalMaxLevel[region];
            }
        }
    }
}


void Foam::refinementSurfaces::findAllIntersections
(
    const pointField& start,
    const pointField& end,
    const labelList& currentLevel,   // current cell refinement level

    const labelList& globalMinLevel,
    const labelList& globalMaxLevel,

    List<pointList>& surfaceLocation,
    List<vectorList>& surfaceNormal,
    labelListList& surfaceLevel
) const
{
    surfaceLevel.setSize(start.size());
    surfaceNormal.setSize(start.size());
    surfaceLocation.setSize(start.size());

    if (surfaces_.empty())
    {
        return;
    }

    // Work arrays
    List<List<pointIndexHit>> hitInfo;

    forAll(surfaces_, surfI)
    {
        const searchableSurface& surface = allGeometry_[surfaces_[surfI]];

        surface.findLineAll(start, end, hitInfo);

        // Repack hits for surface into flat list
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // To avoid overhead of calling getRegion for every point

        label n = 0;
        forAll(hitInfo, pointI)
        {
            n += hitInfo[pointI].size();
        }

        List<pointIndexHit> surfInfo(n);
        labelList pointMap(n);
        n = 0;

        forAll(hitInfo, pointI)
        {
            const List<pointIndexHit>& pHits = hitInfo[pointI];

            forAll(pHits, i)
            {
                surfInfo[n] = pHits[i];
                pointMap[n] = pointI;
                n++;
            }
        }

        labelList surfRegion(n);
        vectorField surfNormal(n);
        surface.getRegion(surfInfo, surfRegion);
        surface.getNormal(surfInfo, surfNormal);

        // Extract back into pointwise
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~

        forAll(surfRegion, i)
        {
            label region = globalRegion(surfI, surfRegion[i]);
            label pointI = pointMap[i];

            if
            (
                currentLevel[pointI] >= globalMinLevel[region]
             && currentLevel[pointI] < globalMaxLevel[region]
            )
            {
                // Append to pointI info
                label sz = surfaceNormal[pointI].size();
                surfaceLocation[pointI].setSize(sz+1);
                surfaceLocation[pointI][sz] = surfInfo[i].hitPoint();

                surfaceNormal[pointI].setSize(sz+1);
                surfaceNormal[pointI][sz] = surfNormal[i];

                surfaceLevel[pointI].setSize(sz+1);
                // Level should just be higher than provided point level.
                // Actual value is not important.
                surfaceLevel[pointI][sz] = globalMaxLevel[region];
            }
        }
    }
}


//void Foam::refinementSurfaces::findAllIntersections
//(
//    const labelList& surfacesToTest,
//    const pointField& start,
//    const pointField& end,
//
//    List<labelList>& hitSurface,            // surface index
//    List<pointList>& hitLocation,           // surface location
//    List<labelList>& hitRegion,
//    List<vectorField>& hitNormal
//) const
//{
//    hitSurface.setSize(start.size());
//    hitLocation.setSize(start.size());
//    hitRegion.setSize(start.size());
//    hitNormal.setSize(start.size());
//
//    if (surfaces_.empty())
//    {
//        return;
//    }
//
//    // Work arrays
//    List<List<pointIndexHit>> hitInfo;
//
//    for (const label surfI : surfacesToTest)
//    {
//        const searchableSurface& surface = allGeometry_[surfaces_[surfI]];
//
//        surface.findLineAll(start, end, hitInfo);
//
//        // Repack hits for surface into flat list
//        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//        // To avoid overhead of calling getRegion for every point
//
//        label n = 0;
//        forAll(hitInfo, pointI)
//        {
//            n += hitInfo[pointI].size();
//        }
//
//        List<pointIndexHit> surfInfo(n);
//        labelList pointMap(n);
//        n = 0;
//
//        forAll(hitInfo, pointI)
//        {
//            const List<pointIndexHit>& pHits = hitInfo[pointI];
//
//            forAll(pHits, i)
//            {
//                surfInfo[n] = pHits[i];
//                pointMap[n] = pointI;
//                n++;
//            }
//        }
//
//        labelList surfRegion(n);
//        vectorField surfNormal(n);
//        surface.getRegion(surfInfo, surfRegion);
//        surface.getNormal(surfInfo, surfNormal);
//
//        surfInfo.clear();
//
//
//        // Extract back into pointwise
//        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//        forAll(surfRegion, i)
//        {
//            const label pointI = pointMap[i];
//
//            // Append to pointI info
//            hitSurface[pointI].append(surfI);
//            hitRegion[pointI].append(globalRegion(surfI, surfRegion[i]));
//            hitLocation[pointI].append(surfInfo[i].hitPoint());
//            hitNormal[pointI].append(surfNormal[i]);
//        }
//    }
//}


void Foam::refinementSurfaces::findNearestIntersection
(
    const labelList& surfacesToTest,
    const pointField& start,
    const pointField& end,

    labelList& surface1,
    List<pointIndexHit>& hit1,
    labelList& region1,
    labelList& surface2,
    List<pointIndexHit>& hit2,
    labelList& region2
) const
{
    // 1. intersection from start to end
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Initialize arguments
    surface1.setSize(start.size());
    surface1 = -1;
    hit1.setSize(start.size());
    region1.setSize(start.size());

    // Current end of segment to test.
    pointField nearest(end);
    // Work array
    List<pointIndexHit> nearestInfo(start.size());
    labelList region;

    forAll(surfacesToTest, testI)
    {
        label surfI = surfacesToTest[testI];

        const searchableSurface& surface = allGeometry_[surfaces_[surfI]];

        // See if any intersection between start and current nearest
        surface.findLine
        (
            start,
            nearest,
            nearestInfo
        );
        surface.getRegion
        (
            nearestInfo,
            region
        );

        forAll(nearestInfo, pointI)
        {
            if (nearestInfo[pointI].hit())
            {
                hit1[pointI] = nearestInfo[pointI];
                surface1[pointI] = surfI;
                region1[pointI] = region[pointI];
                nearest[pointI] = hit1[pointI].hitPoint();
            }
        }
    }


    // 2. intersection from end to last intersection
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Find the nearest intersection from end to start. Note that we initialize
    // to the first intersection (if any).
    surface2 = surface1;
    hit2 = hit1;
    region2 = region1;

    // Set current end of segment to test.
    forAll(nearest, pointI)
    {
        if (hit1[pointI].hit())
        {
            nearest[pointI] = hit1[pointI].point();
        }
        else
        {
            // Disable testing by setting to end.
            nearest[pointI] = end[pointI];
        }
    }

    forAll(surfacesToTest, testI)
    {
        label surfI = surfacesToTest[testI];

        const searchableSurface& surface = allGeometry_[surfaces_[surfI]];

        // See if any intersection between end and current nearest
        surface.findLine
        (
            end,
            nearest,
            nearestInfo
        );
        surface.getRegion
        (
            nearestInfo,
            region
        );

        forAll(nearestInfo, pointI)
        {
            if (nearestInfo[pointI].hit())
            {
                hit2[pointI] = nearestInfo[pointI];
                surface2[pointI] = surfI;
                region2[pointI] = region[pointI];
                nearest[pointI] = hit2[pointI].hitPoint();
            }
        }
    }


    // Make sure that if hit1 has hit something, hit2 will have at least the
    // same point (due to tolerances it might miss its end point)
    forAll(hit1, pointI)
    {
        if (hit1[pointI].hit() && !hit2[pointI].hit())
        {
            hit2[pointI] = hit1[pointI];
            surface2[pointI] = surface1[pointI];
            region2[pointI] = region1[pointI];
        }
    }
}


void Foam::refinementSurfaces::findNearestIntersection
(
    const labelList& surfacesToTest,
    const pointField& start,
    const pointField& end,

    labelList& surface1,
    List<pointIndexHit>& hit1,
    labelList& region1,
    vectorField& normal1,

    labelList& surface2,
    List<pointIndexHit>& hit2,
    labelList& region2,
    vectorField& normal2
) const
{
    // 1. intersection from start to end
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Initialize arguments
    surface1.setSize(start.size());
    surface1 = -1;
    hit1.setSize(start.size());
    region1.setSize(start.size());
    region1 = -1;
    normal1.setSize(start.size());
    normal1 = Zero;

    // Current end of segment to test.
    pointField nearest(end);
    // Work array
    List<pointIndexHit> nearestInfo(start.size());
    labelList region;
    vectorField normal;

    forAll(surfacesToTest, testI)
    {
        label surfI = surfacesToTest[testI];
        const searchableSurface& geom = allGeometry_[surfaces_[surfI]];

        // See if any intersection between start and current nearest
        geom.findLine(start, nearest, nearestInfo);
        geom.getRegion(nearestInfo, region);
        geom.getNormal(nearestInfo, normal);

        forAll(nearestInfo, pointI)
        {
            if (nearestInfo[pointI].hit())
            {
                hit1[pointI] = nearestInfo[pointI];
                surface1[pointI] = surfI;
                region1[pointI] = region[pointI];
                normal1[pointI] = normal[pointI];
                nearest[pointI] = hit1[pointI].hitPoint();
            }
        }
    }


    // 2. intersection from end to last intersection
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Find the nearest intersection from end to start. Note that we initialize
    // to the first intersection (if any).
    surface2 = surface1;
    hit2 = hit1;
    region2 = region1;
    normal2 = normal1;

    // Set current end of segment to test.
    forAll(nearest, pointI)
    {
        if (hit1[pointI].hit())
        {
            nearest[pointI] = hit1[pointI].point();
        }
        else
        {
            // Disable testing by setting to end.
            nearest[pointI] = end[pointI];
        }
    }

    forAll(surfacesToTest, testI)
    {
        label surfI = surfacesToTest[testI];
        const searchableSurface& geom = allGeometry_[surfaces_[surfI]];

        // See if any intersection between end and current nearest
        geom.findLine(end, nearest, nearestInfo);
        geom.getRegion(nearestInfo, region);
        geom.getNormal(nearestInfo, normal);

        forAll(nearestInfo, pointI)
        {
            if (nearestInfo[pointI].hit())
            {
                hit2[pointI] = nearestInfo[pointI];
                surface2[pointI] = surfI;
                region2[pointI] = region[pointI];
                normal2[pointI] = normal[pointI];
                nearest[pointI] = hit2[pointI].hitPoint();
            }
        }
    }


    // Make sure that if hit1 has hit something, hit2 will have at least the
    // same point (due to tolerances it might miss its end point)
    forAll(hit1, pointI)
    {
        if (hit1[pointI].hit() && !hit2[pointI].hit())
        {
            hit2[pointI] = hit1[pointI];
            surface2[pointI] = surface1[pointI];
            region2[pointI] = region1[pointI];
            normal2[pointI] = normal1[pointI];
        }
    }
}


void Foam::refinementSurfaces::findNearestIntersection
(
    const pointField& start,
    const pointField& end,

    labelList& surface1,
    vectorField& normal1
) const
{
    // Initialize arguments
    surface1.setSize(start.size());
    surface1 = -1;
    normal1.setSize(start.size());
    normal1 = Zero;

    // Current end of segment to test.
    pointField nearest(end);
    // Work array
    List<pointIndexHit> nearestInfo(start.size());
    labelList region;
    vectorField normal;

    forAll(surfaces_, surfI)
    {
        const searchableSurface& geom = allGeometry_[surfaces_[surfI]];

        // See if any intersection between start and current nearest
        geom.findLine(start, nearest, nearestInfo);
        geom.getNormal(nearestInfo, normal);

        forAll(nearestInfo, pointI)
        {
            if (nearestInfo[pointI].hit())
            {
                surface1[pointI] = surfI;
                normal1[pointI] = normal[pointI];
                nearest[pointI] = nearestInfo[pointI].point();
            }
        }
    }
}


void Foam::refinementSurfaces::findNearestIntersection
(
//    const labelList& surfacesToTest,
    const pointField& start,
    const pointField& end,

    labelList& surface1,
    labelList& region1,
    List<pointIndexHit>& hitInfo1,
    vectorField& normal1
) const
{
    // Initialize arguments
    surface1.setSize(start.size());
    surface1 = -1;
    region1.setSize(start.size());
    region1 = -1;
    hitInfo1.setSize(start.size());
    hitInfo1 = pointIndexHit();
    normal1.setSize(start.size());
    normal1 = Zero;

    // Current end of segment to test.
    pointField nearest(end);
    // Work array
    List<pointIndexHit> nearestInfo(start.size());
    labelList region;
    vectorField normal;

    forAll(surfaces_, surfI)
    {
        const searchableSurface& geom = allGeometry_[surfaces_[surfI]];

        // See if any intersection between start and current nearest
        geom.findLine(start, nearest, nearestInfo);
        geom.getRegion(nearestInfo, region);
        geom.getNormal(nearestInfo, normal);

        forAll(nearestInfo, pointI)
        {
            if (nearestInfo[pointI].hit())
            {
                surface1[pointI] = surfI;
                region1[pointI] = region[pointI];
                hitInfo1[pointI] = nearestInfo[pointI];
                normal1[pointI] = normal[pointI];
                nearest[pointI] = nearestInfo[pointI].point();
            }
        }
    }
}


void Foam::refinementSurfaces::findNearestIntersection
(
    const pointField& start,
    const pointField& end,

    labelList& surface1,
    List<pointIndexHit>& hitInfo1,
    vectorField& normal1
) const
{
    // 1. intersection from start to end
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Initialize arguments
    surface1.setSize(start.size());
    surface1 = -1;
    hitInfo1.setSize(start.size());
    hitInfo1 = pointIndexHit();
    normal1.setSize(start.size());
    normal1 = Zero;

    // Current end of segment to test.
    pointField nearest(end);
    // Work array
    List<pointIndexHit> nearestInfo(start.size());
    labelList region;
    vectorField normal;

    forAll(surfaces_, surfI)
    {
        const searchableSurface& geom = allGeometry_[surfaces_[surfI]];

        // See if any intersection between start and current nearest
        geom.findLine(start, nearest, nearestInfo);
        geom.getNormal(nearestInfo, normal);

        forAll(nearestInfo, pointI)
        {
            if (nearestInfo[pointI].hit())
            {
                surface1[pointI] = surfI;
                hitInfo1[pointI] = nearestInfo[pointI];
                normal1[pointI] = normal[pointI];
                nearest[pointI] = nearestInfo[pointI].point();
            }
        }
    }
}


void Foam::refinementSurfaces::findAnyIntersection
(
    const pointField& start,
    const pointField& end,

    labelList& hitSurface,
    List<pointIndexHit>& hitInfo
) const
{
    searchableSurfacesQueries::findAnyIntersection
    (
        allGeometry_,
        surfaces_,
        start,
        end,
        hitSurface,
        hitInfo
    );
}


void Foam::refinementSurfaces::findNearest
(
    const labelList& surfacesToTest,
    const pointField& samples,
    const scalarField& nearestDistSqr,
    labelList& hitSurface,
    List<pointIndexHit>& hitInfo
) const
{
    labelList geometries(labelUIndList(surfaces_, surfacesToTest));

    // Do the tests. Note that findNearest returns index in geometries.
    searchableSurfacesQueries::findNearest
    (
        allGeometry_,
        geometries,
        samples,
        nearestDistSqr,
        hitSurface,
        hitInfo
    );

    // Rework the hitSurface to be surface (i.e. index into surfaces_)
    forAll(hitSurface, pointI)
    {
        if (hitSurface[pointI] != -1)
        {
            hitSurface[pointI] = surfacesToTest[hitSurface[pointI]];
        }
    }
}


void Foam::refinementSurfaces::findNearestRegion
(
    const labelList& surfacesToTest,
    const pointField& samples,
    const scalarField& nearestDistSqr,
    labelList& hitSurface,
    labelList& hitRegion
) const
{
    labelList geometries(labelUIndList(surfaces_, surfacesToTest));

    // Do the tests. Note that findNearest returns index in geometries.
    List<pointIndexHit> hitInfo;
    searchableSurfacesQueries::findNearest
    (
        allGeometry_,
        geometries,
        samples,
        nearestDistSqr,
        hitSurface,
        hitInfo
    );

    // Rework the hitSurface to be surface (i.e. index into surfaces_)
    forAll(hitSurface, pointI)
    {
        if (hitSurface[pointI] != -1)
        {
            hitSurface[pointI] = surfacesToTest[hitSurface[pointI]];
        }
    }

    // Collect the region
    hitRegion.setSize(hitSurface.size());
    hitRegion = -1;

    forAll(surfacesToTest, i)
    {
        label surfI = surfacesToTest[i];

        // Collect hits for surfI
        const labelList localIndices(findIndices(hitSurface, surfI));

        List<pointIndexHit> localHits
        (
            UIndirectList<pointIndexHit>
            (
                hitInfo,
                localIndices
            )
        );

        labelList localRegion;
        allGeometry_[surfaces_[surfI]].getRegion(localHits, localRegion);

        forAll(localIndices, i)
        {
            hitRegion[localIndices[i]] = localRegion[i];
        }
    }
}


void Foam::refinementSurfaces::findNearestRegion
(
    const labelList& surfacesToTest,
    const pointField& samples,
    const scalarField& nearestDistSqr,
    labelList& hitSurface,
    List<pointIndexHit>& hitInfo,
    labelList& hitRegion,
    vectorField& hitNormal
) const
{
    labelList geometries(labelUIndList(surfaces_, surfacesToTest));

    // Do the tests. Note that findNearest returns index in geometries.
    searchableSurfacesQueries::findNearest
    (
        allGeometry_,
        geometries,
        samples,
        nearestDistSqr,
        hitSurface,
        hitInfo
    );

    // Rework the hitSurface to be surface (i.e. index into surfaces_)
    forAll(hitSurface, pointI)
    {
        if (hitSurface[pointI] != -1)
        {
            hitSurface[pointI] = surfacesToTest[hitSurface[pointI]];
        }
    }

    // Collect the region
    hitRegion.setSize(hitSurface.size());
    hitRegion = -1;
    hitNormal.setSize(hitSurface.size());
    hitNormal = Zero;

    forAll(surfacesToTest, i)
    {
        label surfI = surfacesToTest[i];

        // Collect hits for surfI
        const labelList localIndices(findIndices(hitSurface, surfI));

        List<pointIndexHit> localHits
        (
            UIndirectList<pointIndexHit>
            (
                hitInfo,
                localIndices
            )
        );

        // Region
        labelList localRegion;
        allGeometry_[surfaces_[surfI]].getRegion(localHits, localRegion);

        forAll(localIndices, i)
        {
            hitRegion[localIndices[i]] = localRegion[i];
        }

        // Normal
        vectorField localNormal;
        allGeometry_[surfaces_[surfI]].getNormal(localHits, localNormal);

        forAll(localIndices, i)
        {
            hitNormal[localIndices[i]] = localNormal[i];
        }
    }
}


//// Find intersection with max of edge. Return -1 or the surface
//// with the highest maxLevel above currentLevel
//Foam::label Foam::refinementSurfaces::findHighestIntersection
//(
//    const point& start,
//    const point& end,
//    const label currentLevel,   // current cell refinement level
//
//    pointIndexHit& maxHit
//) const
//{
//    // surface with highest maxlevel
//    label maxSurface = -1;
//    // maxLevel of maxSurface
//    label maxLevel = currentLevel;
//
//    forAll(*this, surfI)
//    {
//        pointIndexHit hit = operator[](surfI).findLineAny(start, end);
//
//        if (hit.hit())
//        {
//            const triSurface& s = operator[](surfI);
//
//            label region = globalRegion(surfI, s[hit.index()].region());
//
//            if (maxLevel_[region] > maxLevel)
//            {
//                maxSurface = surfI;
//                maxLevel = maxLevel_[region];
//                maxHit = hit;
//            }
//        }
//    }
//
//    if (maxSurface == -1)
//    {
//        // maxLevel unchanged. No interesting surface hit.
//        maxHit.setMiss();
//    }
//
//    return maxSurface;
//}


void Foam::refinementSurfaces::findInside
(
    const labelList& testSurfaces,
    const pointField& pt,
    labelList& insideSurfaces
) const
{
    insideSurfaces.setSize(pt.size());
    insideSurfaces = -1;

    forAll(testSurfaces, i)
    {
        label surfI = testSurfaces[i];

        const searchableSurface& surface = allGeometry_[surfaces_[surfI]];

        const surfaceZonesInfo::areaSelectionAlgo selectionMethod =
            surfZones_[surfI].zoneInside();

        if
        (
            selectionMethod != surfaceZonesInfo::INSIDE
         && selectionMethod != surfaceZonesInfo::OUTSIDE
        )
        {
            FatalErrorInFunction
                << "Trying to use surface "
                << surface.name()
                << " which has non-geometric inside selection method "
                << surfaceZonesInfo::areaSelectionAlgoNames[selectionMethod]
                << exit(FatalError);
        }

        if (surface.hasVolumeType())
        {
            List<volumeType> volType;
            surface.getVolumeType(pt, volType);

            forAll(volType, pointI)
            {
                if (insideSurfaces[pointI] == -1)
                {
                    if
                    (
                        (
                            volType[pointI] == volumeType::INSIDE
                         && selectionMethod == surfaceZonesInfo::INSIDE
                        )
                     || (
                            volType[pointI] == volumeType::OUTSIDE
                         && selectionMethod == surfaceZonesInfo::OUTSIDE
                        )
                    )
                    {
                        insideSurfaces[pointI] = surfI;
                    }
                }
            }
        }
    }
}


void Foam::refinementSurfaces::findNearest
(
    const labelList& surfacesToTest,
    const labelListList& regions,

    const pointField& samples,
    const scalarField& nearestDistSqr,

    labelList& hitSurface,
    List<pointIndexHit>& hitInfo
) const
{
    labelList geometries(labelUIndList(surfaces_, surfacesToTest));

    // Do the tests. Note that findNearest returns index in geometries.
    searchableSurfacesQueries::findNearest
    (
        allGeometry_,
        geometries,
        regions,
        samples,
        nearestDistSqr,
        hitSurface,
        hitInfo
    );

    // Rework the hitSurface to be surface (i.e. index into surfaces_)
    forAll(hitSurface, pointI)
    {
        if (hitSurface[pointI] != -1)
        {
            hitSurface[pointI] = surfacesToTest[hitSurface[pointI]];
        }
    }
}


void Foam::refinementSurfaces::findNearestRegion
(
    const labelList& surfacesToTest,
    const labelListList& regions,

    const pointField& samples,
    const scalarField& nearestDistSqr,

    labelList& hitSurface,
    List<pointIndexHit>& hitInfo,
    labelList& hitRegion,
    vectorField& hitNormal
) const
{
    labelList geometries(labelUIndList(surfaces_, surfacesToTest));

    // Do the tests. Note that findNearest returns index in geometries.
    searchableSurfacesQueries::findNearest
    (
        allGeometry_,
        geometries,
        regions,
        samples,
        nearestDistSqr,
        hitSurface,
        hitInfo
    );

    // Rework the hitSurface to be surface (i.e. index into surfaces_)
    forAll(hitSurface, pointI)
    {
        if (hitSurface[pointI] != -1)
        {
            hitSurface[pointI] = surfacesToTest[hitSurface[pointI]];
        }
    }

    // Collect the region
    hitRegion.setSize(hitSurface.size());
    hitRegion = -1;
    hitNormal.setSize(hitSurface.size());
    hitNormal = Zero;

    forAll(surfacesToTest, i)
    {
        label surfI = surfacesToTest[i];

        // Collect hits for surfI
        const labelList localIndices(findIndices(hitSurface, surfI));

        List<pointIndexHit> localHits
        (
            UIndirectList<pointIndexHit>
            (
                hitInfo,
                localIndices
            )
        );

        // Region
        labelList localRegion;
        allGeometry_[surfaces_[surfI]].getRegion(localHits, localRegion);

        forAll(localIndices, i)
        {
            hitRegion[localIndices[i]] = localRegion[i];
        }

        // Normal
        vectorField localNormal;
        allGeometry_[surfaces_[surfI]].getNormal(localHits, localNormal);

        forAll(localIndices, i)
        {
            hitNormal[localIndices[i]] = localNormal[i];
        }
    }
}


// ************************************************************************* //
