/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2023 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM, distributed under GPL-3.0-or-later.

Description
    Add various time-related options (see Foam::timeSelector)

Required Classes
    - Foam::argList

\*---------------------------------------------------------------------------*/

Foam::argList::addBoolOption
(
    "constant",
    "Include 'constant/' dir in the times"
);

Foam::argList::addBoolOption
(
    "latestTime",
    "Select the latest time"
);

Foam::argList::addBoolOption
(
    "noZero",
    "Exclude '0/' dir from the times"
);

Foam::argList::addOption
(
    "time",
    "value",
    "Select the nearest time to the specified value"
);


// ************************************************************************* //
