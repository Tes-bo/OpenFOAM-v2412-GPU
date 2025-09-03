/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2021 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM, distributed under GPL-3.0-or-later.

Description
    Add -dry-run and -dry-run-write options

Required Classes
    - Foam::argList

\*---------------------------------------------------------------------------*/

Foam::argList::addDryRunOption
(
    "Check case set-up only using a single time step"
);
Foam::argList::addBoolOption
(
    "dry-run-write",
    "Check case set-up and write only using a single time step"
);

// ************************************************************************* //
