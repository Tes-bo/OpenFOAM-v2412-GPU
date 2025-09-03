/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "vectorFieldField.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<template<class> class Field, class Cmpt>
void Foam::zip
(
    FieldField<Field, Vector<Cmpt>>& result,
    const FieldField<Field, Cmpt>& x,
    const FieldField<Field, Cmpt>& y,
    const FieldField<Field, Cmpt>& z
)
{
    forAll(result, i)
    {
        Foam::zip(result[i], x[i], y[i], z[i]);
    }
}


template<template<class> class Field, class Cmpt>
void Foam::unzip
(
    const FieldField<Field, Vector<Cmpt>>& input,
    FieldField<Field, Cmpt>& x,
    FieldField<Field, Cmpt>& y,
    FieldField<Field, Cmpt>& z
)
{
    forAll(input, i)
    {
        Foam::unzip(input[i], x[i], y[i], z[i]);
    }
}


// ************************************************************************* //
