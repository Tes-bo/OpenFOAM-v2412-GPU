/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM, distributed under GPL-3.0-or-later.

Typedef
    Foam::meshSubsetHelper

Description
    Compatibility name. Renamed (JUL-2018) as Foam::fvMeshSubsetProxy

\*---------------------------------------------------------------------------*/

#ifndef meshSubsetHelper_H
#define meshSubsetHelper_H

#include "fvMeshSubsetProxy.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
     typedef fvMeshSubsetProxy meshSubsetHelper;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
