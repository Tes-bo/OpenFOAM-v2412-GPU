/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::PDRobstacle::read(Istream& is)
{
    this->clear();

    const word obsType(is);
    const dictionary dict(is);

    auto* mfuncPtr = readdictionaryMemberFunctionTable(obsType);

    if (!mfuncPtr)
    {
        FatalIOErrorInLookup
        (
            is,
            "obstacle",
            obsType,
            *readdictionaryMemberFunctionTablePtr_
        ) << exit(FatalIOError);
    }

    mfuncPtr(*this, dict);

    return true;
}


Foam::scalar Foam::PDRobstacle::readFiles
(
    const fileName& obsFileDir,
    const wordList& obsFileNames,
    const boundBox& meshBb,

    DynamicList<PDRobstacle>& blocks,
    DynamicList<PDRobstacle>& cylinders
)
{
    blocks.clear();
    cylinders.clear();

    scalar totVolume = 0;
    label nOutside = 0;
    label nProtruding = 0;

    scalar shift = pars.obs_expand;

    if (!obsFileNames.empty())
    {
        Info<< "Reading obstacle files" << nl;
    }

    label maxGroup = -1;

    for (const word& inputFile : obsFileNames)
    {
        Info<< "    file: " << inputFile << nl;

        fileName path = (obsFileDir / inputFile);

        IFstream is(path);
        dictionary inputDict(is);

        const scalar scaleFactor = inputDict.getOrDefault<scalar>("scale", 0);

        const label verbose = inputDict.getOrDefault<label>("verbose", 0);

        for (const entry& dEntry : inputDict)
        {
            if (!dEntry.isDict())
            {
                // ignore non-dictionary entry
                continue;
            }

            const dictionary& dict = dEntry.dict();

            if (!dict.getOrDefault("enabled", true))
            {
                continue;
            }

            label obsGroupId = 0;
            if (dict.readIfPresent("groupId", obsGroupId))
            {
                maxGroup = max(maxGroup, obsGroupId);
            }
            else
            {
                obsGroupId = ++maxGroup;
            }


            pointField pts;
            dict.readIfPresent("locations", pts);
            if (pts.empty())
            {
                pts.resize(1, Zero);
            }

            List<PDRobstacle> obsInput;
            dict.readEntry("obstacles", obsInput);

            label nCyl = 0; // The number of cylinders vs blocks

            for (PDRobstacle& obs : obsInput)
            {
                obs.groupId = obsGroupId;
                obs.scale(scaleFactor);

                if (obs.isCylinder())
                {
                    ++nCyl;
                }
            }

            const label nBlock = (obsInput.size() - nCyl);

            blocks.reserve(blocks.size() + nBlock*pts.size());
            cylinders.reserve(cylinders.size() + nCyl*pts.size());

            if (verbose)
            {
                Info<< "Read " << obsInput.size() << " obstacles ("
                    << nCyl << " cylinders) with "
                    << pts.size() << " locations" << nl;

                if (verbose > 1)
                {
                    Info<< "locations " << pts << nl
                        << "obstacles " << obsInput << nl;
                }
            }

            for (const PDRobstacle& scanObs : obsInput)
            {
                // Reject anything below minimum width
                if (scanObs.tooSmall(pars.min_width))
                {
                    continue;
                }

                for (const point& origin : pts)
                {
                    // A different (very small) shift for each obstacle
                    // so that faces cannot be coincident

                    shift += floatSMALL;
                    const scalar shift2 = shift * 2.0;


                switch (scanObs.typeId)
                {
                    case PDRobstacle::CYLINDER:
                    {
                        // Make a copy
                        PDRobstacle obs(scanObs);

                        // Offset for the group position
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

                        // Offset for the group position
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

                    case PDRobstacle::CUBOID_1:
                    case PDRobstacle::LOUVRE_BLOWOFF:
                    case PDRobstacle::CUBOID:
                    case PDRobstacle::WALL_BEAM:
                    case PDRobstacle::GRATING:
                    case PDRobstacle::RECT_PATCH:
                    {
                        // Make a copy
                        PDRobstacle obs(scanObs);

                        // Offset for the position of the group
                        obs.pt += origin;

                        if (obs.typeId == PDRobstacle::GRATING)
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

            // Info<< "Cylinders: " << cylinders << nl;
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


    Info<< "Number of obstacles: "
        << (blocks.size() + cylinders.size()) << " ("
        << cylinders.size() << " cylinders)" << nl;

    return totVolume;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, PDRobstacle& obs)
{
    obs.read(is);

    return is;
}


// ************************************************************************* //
