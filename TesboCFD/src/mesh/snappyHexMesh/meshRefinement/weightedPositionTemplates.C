/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "vectorTensorTransform.H"
#include "coupledPolyPatch.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<template<class> class Container>
void Foam::weightedPosition::operator()
(
    const coupledPolyPatch& cpp,
    Container<weightedPosition>& map
)
const
{
    Field<point> fld(map.size());
    label i = 0;
    forAllConstIters(map, iter)
    {
        point pt(iter->second());
        if (mag(iter->first()) > VSMALL)
        {
            pt /= iter->first();
        }
        fld[i++] = pt;
    }
    cpp.transformPosition(fld);
    i = 0;
    forAllIters(map, iter)
    {
        iter->second() = fld[i++]*iter->first();
    }
}


// ************************************************************************* //
