/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
    Copyright (C) 2022 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM, distributed under GPL-3.0-or-later.

Description
    Compatibility include (AUG-2022).
    Foam::combineReduce wrapper for Pstream::combineReduce

\*---------------------------------------------------------------------------*/

#ifndef FoamCompat_PstreamCombineReduceOps_H
#define FoamCompat_PstreamCombineReduceOps_H

#include "Pstream.H"
#include "ops.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Compatibility wrapper for Pstream::combineReduce
template<class T, class CombineOp>
void combineReduce
(
    T& value,
    const CombineOp& cop,
    const int tag = UPstream::msgType(),
    const label comm = UPstream::worldComm
)
{
    Pstream::combineReduce(value, cop, tag, comm);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
