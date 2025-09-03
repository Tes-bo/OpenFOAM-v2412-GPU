/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2015 OpenFOAM Foundation
    Copyright (C) 2015-2023 OpenCFD Ltd.
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

#include "surfaceZonesInfo.H"
#include "searchableSurface.H"
#include "searchableSurfaces.H"
#include "polyMesh.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::surfaceZonesInfo::areaSelectionAlgo
>
Foam::surfaceZonesInfo::areaSelectionAlgoNames
({
    { areaSelectionAlgo::INSIDE, "inside" },
    { areaSelectionAlgo::OUTSIDE, "outside" },
    { areaSelectionAlgo::INSIDEPOINT, "insidePoint" },
    { areaSelectionAlgo::NONE, "none" },
});


const Foam::Enum
<
    Foam::surfaceZonesInfo::faceZoneNaming
>
Foam::surfaceZonesInfo::faceZoneNamingNames
({
    { faceZoneNaming::NOZONE, "none" },
    { faceZoneNaming::SINGLE, "single" },
    { faceZoneNaming::REGION, "region" }
});


const Foam::Enum
<
    Foam::surfaceZonesInfo::faceZoneType
>
Foam::surfaceZonesInfo::faceZoneTypeNames
({
    { faceZoneType::INTERNAL, "internal" },
    { faceZoneType::BAFFLE, "baffle" },
    { faceZoneType::BOUNDARY, "boundary" },
});


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceZonesInfo::surfaceZonesInfo
(
    const searchableSurface& surface,
    const dictionary& surfacesDict,
    const wordList& regionNames
)
:
    faceZoneNames_(),
    cellZoneName_(),
    zoneInside_(NONE),
    zoneInsidePoints_(),
    faceType_(INTERNAL)
{
    const label nRegions = surface.regions().size();

    // Old syntax
    surfaceZonesInfo::faceZoneNaming namingType = faceZoneNaming::NOZONE;

    word namingMethod;
    word faceZoneName;
    if (surfacesDict.readIfPresent("faceZone", faceZoneName))
    {
        // Single zone name per surface
        if (surfacesDict.found("faceZoneNaming"))
        {
            FatalIOErrorInFunction(surfacesDict)
                << "Cannot provide both \"faceZone\" and \"faceZoneNaming\""
                << exit(FatalIOError);
        }

        namingType = faceZoneNaming::SINGLE;
        faceZoneNames_.setSize(nRegions, faceZoneName);
    }
    else if (surfacesDict.readIfPresent("faceZoneNaming", namingMethod))
    {
        //namingType = faceZoneNamingNames.get("faceZoneNaming", surfacesDict);
        namingType = faceZoneNamingNames[namingMethod];

        // Generate faceZone names. Maybe make runtime-selection table?
        switch (namingType)
        {
            case faceZoneNaming::NOZONE:
            break;

            case faceZoneNaming::SINGLE:
            {
                // Should already be handled above
                faceZoneNames_.setSize
                (
                    nRegions,
                    surfacesDict.get<word>("faceZone")
                );
            }
            break;

            case faceZoneNaming::REGION:
            {
                faceZoneNames_ = regionNames;
            }
            break;
        }
    }

    if (faceZoneNames_.size())
    {
        if (faceZoneNames_.size() != nRegions)
        {
            FatalIOErrorInFunction(surfacesDict)
                << "Number of faceZones (through 'faceZones' keyword)"
                << " does not correspond to the number of regions "
                << nRegions << " in surface " << surface.name()
                << exit(FatalIOError);
        }

        // Read optional entry to determine inside of faceZone

        word method;
        bool hasSide = surfacesDict.readIfPresent("cellZoneInside", method);
        if (hasSide)
        {
            zoneInside_ = areaSelectionAlgoNames[method];
            if (zoneInside_ == INSIDEPOINT)
            {
                const bool foundPoints = surfacesDict.readIfPresent
                (
                    "insidePoints",
                    zoneInsidePoints_,
                    keyType::LITERAL
                );

                if (foundPoints)
                {
                    if (surfacesDict.found("insidePoint", keyType::LITERAL))
                    {
                        FatalIOErrorInFunction(surfacesDict)
                            << "Cannot supply both 'insidePoint'"
                            << " and 'insidePoints'" << exit(FatalIOError);
                    }
                }
                else
                {
                    zoneInsidePoints_ = pointField
                    (
                        1,
                        surfacesDict.get<point>("insidePoint", keyType::LITERAL)
                    );
                }
            }
        }
        else
        {
            // Check old syntax
            bool inside;
            if (surfacesDict.readIfPresent("zoneInside", inside))
            {
                hasSide = true;
                zoneInside_ = (inside ? INSIDE : OUTSIDE);
            }
        }

        // Read optional cellZone name

        if (surfacesDict.readIfPresent("cellZone", cellZoneName_))
        {
            if
            (
                (
                    zoneInside_ == INSIDE
                 || zoneInside_ == OUTSIDE
                )
            && !surface.hasVolumeType()
            )
            {
                IOWarningInFunction(surfacesDict)
                    << "Illegal entry zoneInside "
                    << areaSelectionAlgoNames[zoneInside_]
                    << " for faceZones "
                    << faceZoneNames_
                    << " since surface is not closed." << endl;
            }
        }
        else if (hasSide)
        {
            IOWarningInFunction(surfacesDict)
                << "Unused entry zoneInside for faceZone "
                << faceZoneNames_
                << " since no cellZone specified."
                << endl;
        }

        // How to handle faces on faceZone
        word faceTypeMethod;
        if (surfacesDict.readIfPresent("faceType", faceTypeMethod))
        {
            faceType_ = faceZoneTypeNames[faceTypeMethod];
        }
    }
}


Foam::surfaceZonesInfo::surfaceZonesInfo
(
    const wordList& faceZoneNames,
    const word& cellZoneName,
    const areaSelectionAlgo& zoneInside,
    const pointField& zoneInsidePoints,
    const faceZoneType& faceType
)
:
    faceZoneNames_(faceZoneNames),
    cellZoneName_(cellZoneName),
    zoneInside_(zoneInside),
    zoneInsidePoints_(zoneInsidePoints),
    faceType_(faceType)
{}


Foam::surfaceZonesInfo::surfaceZonesInfo(const surfaceZonesInfo& surfZone)
:
    faceZoneNames_(surfZone.faceZoneNames()),
    cellZoneName_(surfZone.cellZoneName()),
    zoneInside_(surfZone.zoneInside()),
    zoneInsidePoints_(surfZone.zoneInsidePoints()),
    faceType_(surfZone.faceType())
{}


Foam::labelList Foam::surfaceZonesInfo::getUnnamedSurfaces
(
    const PtrList<surfaceZonesInfo>& surfList
)
{
    labelList anonymousSurfaces(surfList.size());

    label i = 0;
    forAll(surfList, surfI)
    {
        if (surfList[surfI].faceZoneNames().empty())
        {
            anonymousSurfaces[i++] = surfI;
        }
    }
    anonymousSurfaces.setSize(i);

    return anonymousSurfaces;
}


Foam::labelList Foam::surfaceZonesInfo::getNamedSurfaces
(
    const PtrList<surfaceZonesInfo>& surfList
)
{
   labelList namedSurfaces(surfList.size());

    label namedI = 0;
    forAll(surfList, surfI)
    {
        if
        (
            surfList.set(surfI)
         && surfList[surfI].faceZoneNames().size()
        )
        {
            namedSurfaces[namedI++] = surfI;
        }
    }
    namedSurfaces.setSize(namedI);

    return namedSurfaces;
}


Foam::labelList Foam::surfaceZonesInfo::getStandaloneNamedSurfaces
(
    const PtrList<surfaceZonesInfo>& surfList
)
{
   labelList namedSurfaces(surfList.size());

    label namedI = 0;
    forAll(surfList, surfI)
    {
        if
        (
            surfList.set(surfI)
        &&  surfList[surfI].faceZoneNames().size()
        && !surfList[surfI].cellZoneName().size()
        )
        {
            namedSurfaces[namedI++] = surfI;
        }
    }
    namedSurfaces.setSize(namedI);

    return namedSurfaces;
}


Foam::labelList Foam::surfaceZonesInfo::getClosedNamedSurfaces
(
    const PtrList<surfaceZonesInfo>& surfList,
    const searchableSurfaces& allGeometry,
    const labelList& surfaces
)
{
    labelList closed(surfList.size());

    label closedI = 0;
    forAll(surfList, surfI)
    {
        if
        (
            surfList.set(surfI)
         && surfList[surfI].cellZoneName().size()
         && (
                surfList[surfI].zoneInside() == surfaceZonesInfo::INSIDE
             || surfList[surfI].zoneInside() == surfaceZonesInfo::OUTSIDE
            )
         && allGeometry[surfaces[surfI]].hasVolumeType()
        )
        {
            closed[closedI++] = surfI;
        }
    }
    closed.setSize(closedI);

    return closed;
}


Foam::labelList Foam::surfaceZonesInfo::getUnclosedNamedSurfaces
(
    const PtrList<surfaceZonesInfo>& surfList,
    const searchableSurfaces& allGeometry,
    const labelList& surfaces
)
{
    labelList unclosed(surfList.size());

    label unclosedI = 0;
    forAll(surfList, surfI)
    {
        if
        (
            surfList.set(surfI)
         && !allGeometry[surfaces[surfI]].hasVolumeType()
        )
        {
            unclosed[unclosedI++] = surfI;
        }
    }
    unclosed.setSize(unclosedI);

    return unclosed;
}


Foam::labelList Foam::surfaceZonesInfo::getAllClosedNamedSurfaces
(
    const PtrList<surfaceZonesInfo>& surfList,
    const searchableSurfaces& allGeometry,
    const labelList& surfaces
)
{
    labelList closed(surfList.size());

    label closedI = 0;
    forAll(surfList, surfI)
    {
        if
        (
            surfList.set(surfI)
         && surfList[surfI].cellZoneName().size()
         && allGeometry[surfaces[surfI]].hasVolumeType()
        )
        {
            closed[closedI++] = surfI;
        }
    }
    closed.setSize(closedI);

    return closed;
}


Foam::labelList Foam::surfaceZonesInfo::getInsidePointNamedSurfaces
(
    const PtrList<surfaceZonesInfo>& surfList
)
{
    labelList closed(surfList.size());

    label closedI = 0;
    forAll(surfList, surfI)
    {
        if
        (
            surfList.set(surfI)
         && surfList[surfI].cellZoneName().size()
         && surfList[surfI].zoneInside() == surfaceZonesInfo::INSIDEPOINT
        )
        {
            closed[closedI++] = surfI;
        }
    }
    closed.setSize(closedI);

    return closed;
}


Foam::label Foam::surfaceZonesInfo::addCellZone
(
    const word& name,
    const labelList& addressing,
    polyMesh& mesh
)
{
    cellZoneMesh& cellZones = mesh.cellZones();

    label zoneID = cellZones.findZoneID(name);

    if (zoneID == -1)
    {
        zoneID = cellZones.size();

        cellZones.emplace_back
        (
            name,
            addressing,
            zoneID,
            cellZones
        );
    }

    return zoneID;
}


Foam::labelList Foam::surfaceZonesInfo::addCellZonesToMesh
(
    const PtrList<surfaceZonesInfo>& surfList,
    const labelList& namedSurfaces,
    polyMesh& mesh
)
{
    labelList surfaceToCellZone(surfList.size(), -1);

    forAll(namedSurfaces, i)
    {
        label surfI = namedSurfaces[i];

        const word& cellZoneName = surfList[surfI].cellZoneName();

        if (!cellZoneName.empty())
        {
            label zoneI = addCellZone
            (
                cellZoneName,
                labelList(),   // addressing
                mesh
            );

            surfaceToCellZone[surfI] = zoneI;
        }
    }

    // Check they are synced
    List<wordList> allCellZones(Pstream::nProcs());
    allCellZones[Pstream::myProcNo()] = mesh.cellZones().names();
    Pstream::allGatherList(allCellZones);

    for (label proci = 1; proci < allCellZones.size(); proci++)
    {
        if (allCellZones[proci] != allCellZones[0])
        {
            FatalErrorInFunction
                << "Zones not synchronised among processors." << nl
                << " Processor0 has cellZones:" << allCellZones[0]
                << " , processor" << proci
                << " has cellZones:" << allCellZones[proci]
                << exit(FatalError);
        }
    }

    return surfaceToCellZone;
}



Foam::label Foam::surfaceZonesInfo::addFaceZone
(
    const word& name,
    const labelList& addressing,
    const boolList& flipMap,
    polyMesh& mesh
)
{
    faceZoneMesh& faceZones = mesh.faceZones();

    label zoneID = faceZones.findZoneID(name);

    if (zoneID == -1)
    {
        zoneID = faceZones.size();

        faceZones.emplace_back
        (
            name,
            addressing,
            flipMap,
            zoneID,
            faceZones
        );
    }
//    else
//    {
//        // Already have faceZone. Add to addressing (if necessary)
//
//        faceZone& fz = faceZones[zoneI];
//
//        DebugVar(fz.size());
//        DebugVar(fz.addressing().size());
//
//        if (fz.size() != addressing.size())
//        {
//            faceZones.clearAddressing();
//            fz.resetAddressing(addressing, flipMap);
//        }
//        else
//        {
//            const labelList& oldAddressing = fz;
//            const boolList& oldFlipMap = fz.flipMap();
//
//            bitSet isZoneFace(mesh.nFaces(), oldAddressing);
//            bitSet isZoneFlip(mesh.nFaces());
//            forAll(oldAddressing, i)
//            {
//                const label facei = oldAddressing[i];
//                isZoneFlip[facei] = oldFlipMap[i];
//            }
//
//            const bitSet newZoneFace(mesh.nFaces(), addressing);
//            bitSet newZoneFlip(mesh.nFaces());
//            forAll(addressing, i)
//            {
//                if (flipMap[i])
//                {
//                    newZoneFlip.set(addressing[i]);
//                }
//            }
//
//            bool isChanged = false;
//            forAll(isZoneFace, facei)
//            {
//                if
//                (
//                    isZoneFace[facei] != newZoneFace[facei]
//                 || isZoneFlip[facei] != newZoneFlip[facei]
//                )
//                {
//                    isZoneFace[facei] = newZoneFace[facei];
//                    isZoneFlip[facei] = newZoneFlip[facei];
//                    isChanged = true;
//                }
//            }
//
//            if (isChanged)
//            {
//                labelList newAddressing(isZoneFace.sortedToc());
//                boolList newFlip(newAddressing.size(), false);
//                forAll(newAddressing, i)
//                {
//                    newFlip[i] = isZoneFlip[newAddressing[i]];
//                }
//                faceZones.clearAddressing();
//                fz.resetAddressing
//                (
//                    std::move(newAddressing),
//                    std::move(newFlip)
//                );
//            }
//        }
//    }

    return zoneID;
}


Foam::labelListList Foam::surfaceZonesInfo::addFaceZonesToMesh
(
    const PtrList<surfaceZonesInfo>& surfList,
    const labelList& namedSurfaces,
    polyMesh& mesh
)
{
    labelListList surfaceToFaceZones(surfList.size());

    faceZoneMesh& faceZones = mesh.faceZones();

    forAll(namedSurfaces, i)
    {
        label surfI = namedSurfaces[i];

        const wordList& faceZoneNames = surfList[surfI].faceZoneNames();

        surfaceToFaceZones[surfI].setSize(faceZoneNames.size(), -1);
        forAll(faceZoneNames, j)
        {
            const word& faceZoneName = faceZoneNames[j];

            label zoneI = addFaceZone
            (
                faceZoneName,   //name
                labelList(),    //addressing
                boolList(),     //flipmap
                mesh
            );

            surfaceToFaceZones[surfI][j] = zoneI;
        }
    }

    // Check they are synced
    List<wordList> allFaceZones(Pstream::nProcs());
    allFaceZones[Pstream::myProcNo()] = faceZones.names();
    Pstream::allGatherList(allFaceZones);

    for (label proci = 1; proci < allFaceZones.size(); proci++)
    {
        if (allFaceZones[proci] != allFaceZones[0])
        {
            FatalErrorInFunction
                << "Zones not synchronised among processors." << nl
                << " Processor0 has faceZones:" << allFaceZones[0]
                << " , processor" << proci
                << " has faceZones:" << allFaceZones[proci]
                << exit(FatalError);
        }
    }

    return surfaceToFaceZones;
}


// ************************************************************************* //
