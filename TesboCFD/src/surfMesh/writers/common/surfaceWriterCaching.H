/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM, distributed under GPL-3.0-or-later.

Typedef
    Foam::surfaceWriters::writerCaching

Description
    Compatibility name. Replaced (JAN-2022) by
    Foam::ensightOutput::writerCaching

\*---------------------------------------------------------------------------*/

#ifndef Foam_surfaceWriters_writerCaching_H
#define Foam_surfaceWriters_writerCaching_H

#include "ensightWriterCaching.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceWriters
{

typedef Foam::ensightOutput::writerCaching writerCaching;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceWriters
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
