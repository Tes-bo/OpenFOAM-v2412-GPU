/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 Shell Research Ltd.
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#include "PDRsetFields.H"
#include "PDRobstacle.H"
#include "volumeType.H"

using namespace Foam;
using namespace Foam::constant;

#undef USE_ZERO_INSTANCE_GROUPS
// #define USE_ZERO_INSTANCE_GROUPS


// Counting
Foam::labelPair Foam::PDRlegacy::readObstacleFiles
(
    const fileName& obsFileDir,
    const wordList& obsFileNames,
    Map<obstacleGrouping>& groups
)
{
    // Default group with single instance and position (0,0,0)
    groups(0).clear();
    groups(0).append(point::zero);

    string buffer;

    if (!obsFileNames.empty())
    {
        Info<< "Counting groups in obstacle files" << nl;
    }
    for (const word& inputFile : obsFileNames)
    {
        Info<< "    file: " << inputFile << nl;

        fileName path = (obsFileDir / inputFile);

        IFstream is(path);

        while (is.good())
        {
            // Process each line of obstacle files
            is.getLine(buffer);

            const auto firstCh = buffer.find_first_not_of(" \t\n\v\f\r");

            if
            (
                firstCh == std::string::npos
             || buffer[firstCh] == '#'
            )
            {
                // Empty line or comment line
                continue;
            }

            int typeId;
            double x, y, z; // Double (not scalar) to match sscanf spec

            if
            (
                sscanf(buffer.c_str(), "%d %lf %lf %lf", &typeId, &x, &y, &z)<4
             || typeId == 0
             || typeId == PDRobstacle::MESH_PLANE
            )
            {
                continue;
            }

            x *= pars.scale;
            y *= pars.scale;
            z *= pars.scale;

            const label groupId = typeId / 100;
            typeId %= 100;

            if (typeId == PDRobstacle::OLD_INLET)
            {
                Info<< "Ignored old-inlet type" << nl;
                continue;
            }
            else if (typeId == PDRobstacle::GRATING && pars.ignoreGratings)
            {
                Info<< "Ignored grating" << nl;
                continue;
            }

            if (typeId == 0)
            {
                // Defining a group location
                groups(groupId).append(x, y, z);
            }
            else if (PDRobstacle::isCylinder(typeId))
            {
                // Increment cylinder count for the group
                groups(groupId).addCylinder();
            }
            else
            {
                // Increment obstacle count for the group
                groups(groupId).addObstacle();
            }
        }
    }


    label nTotalObs = 0;
    label nTotalCyl = 0;

    label nMissedObs = 0;
    label nMissedCyl = 0;

    forAllConstIters(groups, iter)
    {
        const auto& group = iter.val();

        nTotalObs += group.nTotalObstacle();
        nTotalCyl += group.nTotalCylinder();

        if (group.empty())
        {
            nMissedObs += group.nObstacle();
            nMissedCyl += group.nCylinder();
        }
    }

    for (const label groupId : groups.sortedToc())
    {
        const auto& group = groups[groupId];

        if (groupId)
        {
            if (group.size())
            {
                Info<< "Found " << group.size()
                    << " instances of group " << groupId << " ("
                    << group.nObstacle() << " obstacles "
                    << group.nCylinder() << " cylinders)"
                    << nl;
            }
        }
        else
        {
            // The group 0 is for ungrouped obstacles
            Info<< "Found "
                << group.nObstacle() << " obstacles "
                << group.nCylinder() << " cylinders not in groups" << nl;
        }
    }

    Info<< "Number of obstacles: "
        << (nTotalObs + nTotalCyl) << " ("
        << nTotalCyl << " cylinders)" << nl;

    if (nMissedObs + nMissedCyl)
    {
        #ifdef USE_ZERO_INSTANCE_GROUPS

        nTotalObs += nMissedObs;
        nTotalCyl += nMissedCyl;
        Info<< "Adding " << (nMissedObs + nMissedCyl)
            << " obstacles in groups without instances to default group" << nl;

        #else

        Warning
            << nl << "Found " << (nMissedObs + nMissedCyl)
            << " obstacles in groups without instances" << nl << nl;

        if (pars.debugLevel > 1)
        {
            for (const label groupId : groups.sortedToc())
            {
                const auto& group = groups[groupId];

                if
                (
                    groupId && group.empty()
                 && (group.nObstacle() || group.nCylinder())
                )
                {
                    Info<< "    Group " << groupId << " ("
                        << group.nObstacle() << " obstacles "
                        << group.nCylinder() << " cylinders)"
                        << nl;
                }
            }
        }
        #endif
    }

    return labelPair(nTotalObs, nTotalCyl);
}


Foam::scalar Foam::PDRlegacy::readObstacleFiles
(
    const fileName& obsFileDir,
    const wordList& obsFileNames,
    const Map<obstacleGrouping>& groups,
    const boundBox& meshBb,

    DynamicList<PDRobstacle>& blocks,
    DynamicList<PDRobstacle>& cylinders
)
{
    // Catch programming errors
    if (!groups.found(0))
    {
        FatalErrorInFunction
            << "No default group 0 defined!" << nl
            << exit(FatalError);
    }

    scalar totVolume = 0;
    label nOutside = 0;
    label nProtruding = 0;

    scalar shift = pars.obs_expand;

    string buffer;

    if (!obsFileNames.empty())
    {
        Info<< "Reading obstacle files" << nl;
    }

    for (const word& inputFile : obsFileNames)
    {
        Info<< "    file: " << inputFile << nl;

        fileName path = (obsFileDir / inputFile);

        IFstream is(path);

        label lineNo = 0;
        while (is.good())
        {
            // Process each line of obstacle files
            ++lineNo;
            is.getLine(buffer);

            const auto firstCh = buffer.find_first_not_of(" \t\n\v\f\r");

            if
            (
                firstCh == std::string::npos
             || buffer[firstCh] == '#'
            )
            {
                // Empty line or comment line
                continue;
            }

            // Quick reject

            int typeId;     // Int (not label) to match sscanf spec
            double x, y, z; // Double (not scalar) to match sscanf spec

            if
            (
                sscanf(buffer.c_str(), "%d %lf %lf %lf", &typeId, &x, &y, &z) < 4
             || typeId == 0
             || typeId == PDRobstacle::MESH_PLANE
            )
            {
                continue;
            }

            int groupId = typeId / 100;
            typeId %= 100;

            if
            (
                typeId == PDRobstacle::OLD_INLET
             || (typeId == PDRobstacle::GRATING && pars.ignoreGratings)
            )
            {
                // Silent - already warned during counting
                continue;
            }

            if (typeId == 0)
            {
                // Group location - not an obstacle
                continue;
            }

            if (!groups.found(groupId))
            {
                // Catch programming errors.
                // - group should be there after the previous read
                Warning
                    << "Encountered undefined group: " << groupId << nl;
                continue;
            }

            #ifdef USE_ZERO_INSTANCE_GROUPS
            const obstacleGrouping& group =
            (
                groups[groups[groupId].size() ? groupId : 0]
            );
            #else
            const obstacleGrouping& group = groups[groupId];
            #endif

            // Add the obstacle to the list with different position
            // offsets according to its group. Non-group obstacles
            // are treated as group 0, which has a single instance
            // with position (0,0,0) and are added only once.

            PDRobstacle scanObs;

            if
            (
                !scanObs.setFromLegacy
                (
                    (groupId * 100) + typeId,
                    buffer,
                    lineNo,
                    inputFile
                )
            )
            {
                continue;
            }

            scanObs.scale(pars.scale);

            // Ignore anything below minimum width
            if (scanObs.tooSmall(pars.min_width))
            {
                continue;
            }


            for (const point& origin : group)
            {
                // A different (very small) shift for each obstacle
                // so that faces cannot be coincident
                shift += floatSMALL;
                const scalar shift2 = shift * 2.0;

                switch (typeId)
                {
                    case PDRobstacle::CYLINDER:
                    {
                        // Make a copy
                        PDRobstacle obs(scanObs);

                        // Offset for the position of the group
                        obs.pt += origin;

                        // Shift the end outwards so, if exactly on
                        // cell boundary, now overlap cell.
                        // So included in Aw.
                        obs.pt -= point::uniform(shift);
                        obs.len() += shift2;
                        obs.dia() -= floatSMALL;


                        // Trim against the mesh bounds.
                        // Ignore if it doesn't overlap, or bounds are invalid
                        const volumeType vt = obs.trim(meshBb);

                        switch (vt)
                        {
                            case volumeType::OUTSIDE:
                                ++nOutside;
                                continue; // Can ignore the rest
                                break;

                            case volumeType::MIXED:
                                ++nProtruding;
                                break;

                            default:
                                break;
                        }

                        // Later for position sorting
                        switch (obs.orient)
                        {
                            case vector::X:
                                obs.sortBias = obs.len();
                                break;
                            case vector::Y:
                                obs.sortBias = 0.5*obs.dia();
                                break;
                            case vector::Z:
                                obs.sortBias = 0.5*obs.dia();
                                break;
                        }

                        totVolume += obs.volume();
                        cylinders.append(obs);

                        break;
                    }

                    case PDRobstacle::DIAG_BEAM:
                    {
                        // Make a copy
                        PDRobstacle obs(scanObs);

                        // Offset for the position of the group
                        obs.pt += origin;

                        // Shift the end outwards so, if exactly on
                        // cell boundary, now overlap cell.
                        // So included in Aw.
                        obs.pt -= point::uniform(shift);
                        obs.len() += shift2;
                        obs.wa += shift2;
                        obs.wb += shift2;

                        totVolume += obs.volume();
                        cylinders.append(obs);

                        break;
                    }

                    case PDRobstacle::CUBOID_1:       // Cuboid "Type 1"
                    case PDRobstacle::LOUVRE_BLOWOFF: // Louvred wall or blow-off panel
                    case PDRobstacle::CUBOID:         // Cuboid
                    case PDRobstacle::WALL_BEAM:      // Beam against wall (treated here as normal cuboid)
                    case PDRobstacle::GRATING:        // Grating
                    case PDRobstacle::RECT_PATCH:     // Inlet, outlet or ather b.c. (rectangular)
                    {
                        // Make a copy
                        PDRobstacle obs(scanObs);

                        // Offset for the position of the group
                        obs.pt += origin;

                        if (typeId == PDRobstacle::GRATING)
                        {
                            if (obs.slat_width <= 0)
                            {
                                obs.slat_width = pars.def_grating_slat_w;
                            }
                        }

                        // Shift the end outwards so, if exactly on
                        // cell boundary, now overlap cell.
                        // So included in Aw.
                        obs.pt -= point::uniform(shift);
                        obs.span += point::uniform(shift2);


                        // Trim against the mesh bounds.
                        // Ignore if it doesn't overlap, or bounds are invalid
                        const volumeType vt = obs.trim(meshBb);

                        switch (vt)
                        {
                            case volumeType::OUTSIDE:
                                ++nOutside;
                                continue; // Can ignore the rest
                                break;

                            case volumeType::MIXED:
                                ++nProtruding;
                                break;

                            default:
                                break;
                        }

                        totVolume += obs.volume();

                        blocks.append(obs);

                        break;
                    }
                }
            }
        }

        if (nOutside || nProtruding)
        {
            Info<< "Warning: " << nOutside << " obstacles outside domain, "
                << nProtruding << " obstacles partly outside domain" << nl;
        }
    }


    // #ifdef FULLDEBUG
    // Info<< blocks << nl << cylinders << nl;
    // #endif

    return totVolume;
}


Foam::scalar Foam::PDRobstacle::legacyReadFiles
(
    const fileName& obsFileDir,
    const wordList& obsFileNames,
    const boundBox& meshBb,
    DynamicList<PDRobstacle>& blocks,
    DynamicList<PDRobstacle>& cylinders
)
{
    // Still just with legacy reading

    // Count the obstacles and get the group locations
    Map<PDRlegacy::obstacleGrouping> groups;

    const labelPair obsCounts =
        PDRlegacy::readObstacleFiles(obsFileDir, obsFileNames, groups);

    const label nObstacle = obsCounts.first();
    const label nCylinder = obsCounts.second();

    // Info<< "grouping: " << groups << endl;

    if (!nObstacle && !nCylinder)
    {
        FatalErrorInFunction
            << "No obstacles in domain" << nl
            << exit(FatalError);
    }

    blocks.clear();
    blocks.reserve(4 * max(4, nObstacle));

    cylinders.clear();
    cylinders.reserve(4 * max(4, nCylinder));

    return PDRlegacy::readObstacleFiles
    (
        obsFileDir, obsFileNames, groups,
        meshBb,
        blocks,
        cylinders
    );
}


// ************************************************************************* //
