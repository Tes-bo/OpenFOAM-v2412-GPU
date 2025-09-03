/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2023 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM, distributed under GPL-3.0-or-later.

Description
    Gets single or multiple region names based on command-line options:
    (-allRegions | -regions | -regions)

Priority
    1. -allRegions
    2. -regions = specify multiple regions to select, or a single region
    3. -region = specify a single region

Note
    There is no semantical difference between "-regions name"
    and "-region name"

Required Classes
    - Foam::polyMesh
    - Foam::regionProperties
    - Foam::IOobjectOption

Required Variables
    - args [argList]
    - runTime [Time]

Provides Variables
    - regionNames [wordList]

See Also
    addAllRegionOptions.H

\*---------------------------------------------------------------------------*/

wordList regionNames;
{
    wordRes selectByName;

    if (args.found("allRegions"))
    {
        regionNames =
        (
            Foam::regionProperties(runTime, IOobjectOption::READ_IF_PRESENT)
            .names()
        );

        if (regionNames.empty())
        {
            InfoErr
                << "Warning: No regionProperties, assume default region"
                << nl << endl;
        }
        else
        {
            Info<< "Using all regions: " << flatOutput(regionNames) << nl;
        }
    }
    else if (args.readListIfPresent<wordRe>("regions", selectByName))
    {
        if
        (
            selectByName.size() == 1
         && selectByName.front().isLiteral()
        )
        {
            // Identical to -region NAME
            regionNames.resize(1);
            regionNames.front() = selectByName.front();
        }
        else if (selectByName.size())
        {
            regionNames =
            (
                Foam::regionProperties(runTime, IOobjectOption::READ_IF_PRESENT)
                .names()
            );

            if (regionNames.empty())
            {
                Info<< "Warning: No regionProperties, assume default region"
                    << nl << endl;
            }
            else
            {
                const labelList matching = selectByName.matching(regionNames);

                if (matching.empty())
                {
                    InfoErr
                        << "No match in regions: "
                        << flatOutput(regionNames) << nl
                        << "... stopping" << nl << endl;
                    return 1;
                }

                regionNames = wordList(regionNames, matching);

                Info<< "Using regions: " << flatOutput(regionNames) << nl;
            }
        }
    }
    else if (args.found("region"))
    {
        regionNames.resize(1);
        regionNames.front() = args.get<word>("region");
    }

    // Fallback to defaultRegion
    if (regionNames.empty())
    {
        regionNames.resize(1);
        regionNames.front() = Foam::polyMesh::defaultRegion;
    }
}


// ************************************************************************* //
