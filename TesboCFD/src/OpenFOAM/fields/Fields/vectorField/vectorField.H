/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2022 OpenCFD Ltd.
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

Typedef
    Foam::vectorField

Description
    Specialisation of Field\<T\> for vector.

SourceFiles
    vectorFieldTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_vectorField_H
#define Foam_vectorField_H

#include "scalarField.H"
#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef Field<vector> vectorField;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Inplace normalise (float) vector field
template<> void Field<Vector<float>>::normalise();

//- Inplace normalise (double) vector field
template<> void Field<Vector<double>>::normalise();


//- Zip together vector field from components
template<class Cmpt>
void zip
(
    Field<Vector<Cmpt>>& result,
    const UList<Cmpt>& x,
    const UList<Cmpt>& y,
    const UList<Cmpt>& z
);

//- Zip together vector field from components
template<class Cmpt>
tmp<Field<Vector<Cmpt>>> zip
(
    const Field<Cmpt>& x,
    const Field<Cmpt>& y,
    const Field<Cmpt>& z
);

//- Unzip vector field into components
template<class Cmpt>
void unzip
(
    const UList<Vector<Cmpt>>& input,
    Field<Cmpt>& x,
    Field<Cmpt>& y,
    Field<Cmpt>& z
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "vectorFieldTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
