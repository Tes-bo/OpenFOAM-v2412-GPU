/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM, distributed under GPL-3.0-or-later.

Description
    Add multi-region command-line options: -allRegions, -regions, -region

Required Classes
    - Foam::argList

See Also
    getAllRegionOptions.H

\*---------------------------------------------------------------------------*/

{
    Foam::argList::addBoolOption
    (
        "allRegions",
        "Use all regions in regionProperties"
    );

    Foam::argList::addOption
    (
        "regions",
        "wordRes",
        "Use specified mesh region. Eg, -regions gas\n"
        "Or from regionProperties.  Eg, -regions '(gas \"solid.*\")'"
    );

    Foam::argList::addOption
    (
        "region",
        "name",
        "Use specified mesh region. Eg, -region gas"
    );
}


// ************************************************************************* //
