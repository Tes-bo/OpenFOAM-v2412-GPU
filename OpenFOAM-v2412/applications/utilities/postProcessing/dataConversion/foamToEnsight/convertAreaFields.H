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
    Code chunk for converting area fields
    included by foamToEnsight.

    Possible objects to convert are given in faObjects

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Area fields (finiteArea)
if (doFiniteArea && ensFaCasePtr && ensFaMeshPtr)
{
    Info<< "    area   field (";

    writeAllAreaFields(*ensFaCasePtr, *ensFaMeshPtr, faObjects);

    Info<< " )" << nl;
}


// ************************************************************************* //
