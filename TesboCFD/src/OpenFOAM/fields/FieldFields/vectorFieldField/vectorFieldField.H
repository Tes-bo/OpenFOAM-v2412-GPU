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

InClass
    Foam::vectorFieldField

Description
    Specialisation of FieldField\<T\> for vector.

SourceFiles
    vectorFieldField.C

\*---------------------------------------------------------------------------*/

#ifndef vectorFieldField_H
#define vectorFieldField_H

#include "scalarFieldField.H"
#include "vectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global Functions * * * * * * * * * * * * * * //

//- Compose a vector field field from components
template<template<class> class Field, class Cmpt>
void zip
(
    FieldField<Field, Vector<Cmpt>>& result,
    const FieldField<Field, Cmpt>& x,
    const FieldField<Field, Cmpt>& y,
    const FieldField<Field, Cmpt>& z
);

//- Unzip vector field field into components
template<template<class> class Field, class Cmpt>
void unzip
(
    const FieldField<Field, Vector<Cmpt>>& input,
    FieldField<Field, Cmpt>& x,
    FieldField<Field, Cmpt>& y,
    FieldField<Field, Cmpt>& z
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "vectorFieldField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
