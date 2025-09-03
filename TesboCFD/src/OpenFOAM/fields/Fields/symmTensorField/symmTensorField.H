/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2019-2023 OpenCFD Ltd.
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
    Foam::symmTensorField

Description
    Specialisation of Field\<T\> for symmTensor.

SourceFiles
    symmTensorField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_symmTensorField_H
#define Foam_symmTensorField_H

#include "scalarField.H"
#include "vectorField.H"
#include "sphericalTensor.H"
#include "symmTensor.H"
#include "tensor.H"

#define TEMPLATE
#include "FieldFunctionsM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

typedef Field<symmTensor> symmTensorField;

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

//- Zip together symmTensor field from components
template<class Cmpt>
void zip
(
    Field<SymmTensor<Cmpt>>& result,
    const UList<Cmpt>& xx, const UList<Cmpt>& xy, const UList<Cmpt>& xz,
    const UList<Cmpt>& yy, const UList<Cmpt>& yz,
    const UList<Cmpt>& zz
);

//- Unzip symmTensor field into components
template<class Cmpt>
void unzip
(
    const UList<SymmTensor<Cmpt>>& input,
    Field<Cmpt>& xx, Field<Cmpt>& xy, Field<Cmpt>& xz,
    Field<Cmpt>& yy, Field<Cmpt>& yz,
    Field<Cmpt>& zz
);

//- Zip together symmTensor field from components
template<class Cmpt>
tmp<Field<SymmTensor<Cmpt>>> zip
(
    const Field<Cmpt>& xx, const Field<Cmpt>& xy, const Field<Cmpt>& xz,
    const Field<Cmpt>& yy, const Field<Cmpt>& yz,
    const Field<Cmpt>& zz
);


//- Zip together symmTensor field from row components
template<class Cmpt>
void zipRows
(
    Field<SymmTensor<Cmpt>>& result,
    const UList<Vector<Cmpt>>& x,
    const UList<Vector<Cmpt>>& y,
    const UList<Vector<Cmpt>>& z
);

//- Zip together symmTensor field from column components
template<class Cmpt>
void zipCols
(
    Field<SymmTensor<Cmpt>>& result,
    const UList<Vector<Cmpt>>& x,
    const UList<Vector<Cmpt>>& y,
    const UList<Vector<Cmpt>>& z
)
{
    zipRows(result, x, y, z);
}


//- Extract symmTensor field rows
template<class Cmpt>
void unzipRows
(
    const UList<SymmTensor<Cmpt>>& input,
    Field<Vector<Cmpt>>& x,
    Field<Vector<Cmpt>>& y,
    Field<Vector<Cmpt>>& z
);

//- Extract symmTensor field columns
template<class Cmpt>
void unzipCols
(
    const UList<SymmTensor<Cmpt>>& input,
    Field<Vector<Cmpt>>& x,
    Field<Vector<Cmpt>>& y,
    Field<Vector<Cmpt>>& z
)
{
    unzipRows(input, x, y, z);
}


//- Extract a symmTensor field row (x,y,z) == (0,1,2)
template<class Cmpt>
void unzipRow
(
    const UList<SymmTensor<Cmpt>>& input,
    const direction idx,  //!< vector::components
    Field<Vector<Cmpt>>& result
);

//- Extract a symmTensor field column (x,y,z) == (0,1,2)
template<class Cmpt>
void unzipCol
(
    const UList<SymmTensor<Cmpt>>& input,
    const direction idx,  //!< vector::components
    Field<Vector<Cmpt>>& result
)
{
    unzipRow(input, idx, result);
}


//- Extract a symmTensor field diagonal
template<class Cmpt>
void unzipDiag
(
    const UList<SymmTensor<Cmpt>>& input,
    Field<Vector<Cmpt>>& result
);

//- Extract a symmTensor field row (x,y,z) == (0,1,2)
template<class Cmpt>
tmp<Field<Vector<Cmpt>>> unzipRow
(
    const Field<SymmTensor<Cmpt>>& input,
    const direction idx  //!< vector::components
);

//- Extract a symmTensor field column (x,y,z) == (0,1,2)
template<class Cmpt>
tmp<Field<Vector<Cmpt>>> unzipCol
(
    const Field<SymmTensor<Cmpt>>& input,
    const direction idx  //!< vector::components
)
{
    return unzipRow(input, idx);
}

//- Extract a symmTensor field diagonal
template<class Cmpt>
tmp<Field<Vector<Cmpt>>> unzipDiag
(
    const Field<SymmTensor<Cmpt>>& input
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

UNARY_FUNCTION(symmTensor, vector, sqr)
UNARY_FUNCTION(symmTensor, symmTensor, innerSqr)

UNARY_FUNCTION(scalar, symmTensor, tr)
UNARY_FUNCTION(sphericalTensor, symmTensor, sph)
UNARY_FUNCTION(symmTensor, symmTensor, symm)
UNARY_FUNCTION(symmTensor, symmTensor, twoSymm)
UNARY_FUNCTION(symmTensor, symmTensor, devSymm)
UNARY_FUNCTION(symmTensor, symmTensor, devTwoSymm)
UNARY_FUNCTION(symmTensor, symmTensor, dev)
UNARY_FUNCTION(symmTensor, symmTensor, dev2)
UNARY_FUNCTION(scalar, symmTensor, det)
UNARY_FUNCTION(symmTensor, symmTensor, cof)
UNARY_FUNCTION(symmTensor, symmTensor, inv)
UNARY_FUNCTION(symmTensor, symmTensor, pinv)


// * * * * * * * * * * * * * * * global operators  * * * * * * * * * * * * * //

UNARY_OPERATOR(vector, symmTensor, *, hdual)

BINARY_OPERATOR(tensor, symmTensor, symmTensor, &, dot)
BINARY_TYPE_OPERATOR(tensor, symmTensor, symmTensor, &, dot)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "undefFieldFunctionsM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "symmTensorFieldTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
