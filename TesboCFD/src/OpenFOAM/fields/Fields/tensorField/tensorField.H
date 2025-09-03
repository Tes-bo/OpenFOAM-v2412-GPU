/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::tensorField

Description
    Specialisation of Field\<T\> for tensor.

SourceFiles
    tensorField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_tensorField_H
#define Foam_tensorField_H

#include "scalarField.H"
#include "vectorField.H"
#include "sphericalTensorField.H"
#include "symmTensorField.H"
#include "tensor.H"

#define TEMPLATE
#include "FieldFunctionsM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

typedef Field<tensor> tensorField;

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

//- Zip together tensor field from components
template<class Cmpt>
void zip
(
    Field<Tensor<Cmpt>>& result,
    const UList<Cmpt>& xx, const UList<Cmpt>& xy, const UList<Cmpt>& xz,
    const UList<Cmpt>& yx, const UList<Cmpt>& yy, const UList<Cmpt>& yz,
    const UList<Cmpt>& zx, const UList<Cmpt>& zy, const UList<Cmpt>& zz
);

//- Unzip tensor field into components
template<class Cmpt>
void unzip
(
    const UList<Tensor<Cmpt>>& input,
    Field<Cmpt>& xx, Field<Cmpt>& xy, Field<Cmpt>& xz,
    Field<Cmpt>& yx, Field<Cmpt>& yy, Field<Cmpt>& yz,
    Field<Cmpt>& zx, Field<Cmpt>& zy, Field<Cmpt>& zz
);

//- Zip together tensor field from components
template<class Cmpt>
tmp<Field<Tensor<Cmpt>>> zip
(
    const Field<Cmpt>& xx, const Field<Cmpt>& xy, const Field<Cmpt>& xz,
    const Field<Cmpt>& yx, const Field<Cmpt>& yy, const Field<Cmpt>& yz,
    const Field<Cmpt>& zx, const Field<Cmpt>& zy, const Field<Cmpt>& zz
);


//- Zip together tensor field from row components
template<class Cmpt>
void zipRows
(
    Field<Tensor<Cmpt>>& result,
    const UList<Vector<Cmpt>>& x,
    const UList<Vector<Cmpt>>& y,
    const UList<Vector<Cmpt>>& z
);

//- Zip together tensor field from column components
template<class Cmpt>
void zipCols
(
    Field<Tensor<Cmpt>>& result,
    const UList<Vector<Cmpt>>& x,
    const UList<Vector<Cmpt>>& y,
    const UList<Vector<Cmpt>>& z
);


//- Extract tensor field rows
template<class Cmpt>
void unzipRows
(
    const UList<Tensor<Cmpt>>& input,
    Field<Vector<Cmpt>>& x,
    Field<Vector<Cmpt>>& y,
    Field<Vector<Cmpt>>& z
);

//- Extract tensor field columns
template<class Cmpt>
void unzipCols
(
    const UList<Tensor<Cmpt>>& input,
    Field<Vector<Cmpt>>& x,
    Field<Vector<Cmpt>>& y,
    Field<Vector<Cmpt>>& z
);


//- Extract a tensor field row (x,y,z) == (0,1,2)
template<class Cmpt>
void unzipRow
(
    const UList<Tensor<Cmpt>>& input,
    const direction idx,  //!< vector::components
    Field<Vector<Cmpt>>& result
);

//- Extract a tensor field column (x,y,z) == (0,1,2)
template<class Cmpt>
void unzipCol
(
    const UList<Tensor<Cmpt>>& input,
    const direction idx,  //!< vector::components
    Field<Vector<Cmpt>>& result
);

//- Extract a tensor field diagonal
template<class Cmpt>
void unzipDiag
(
    const UList<Tensor<Cmpt>>& input,
    Field<Vector<Cmpt>>& result
);

//- Extract a tensor field row (x,y,z) == (0,1,2)
template<class Cmpt>
tmp<Field<Vector<Cmpt>>> unzipRow
(
    const Field<Tensor<Cmpt>>& input,
    const direction idx  //!< vector::components
);

//- Extract a tensor field column (x,y,z) == (0,1,2)
template<class Cmpt>
tmp<Field<Vector<Cmpt>>> unzipCol
(
    const Field<Tensor<Cmpt>>& input,
    const direction idx  //!< vector::components
);

//- Extract a tensor field diagonal
template<class Cmpt>
tmp<Field<Vector<Cmpt>>> unzipDiag
(
    const Field<Tensor<Cmpt>>& input
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

UNARY_FUNCTION(scalar, tensor, tr)
UNARY_FUNCTION(sphericalTensor, tensor, sph)
UNARY_FUNCTION(symmTensor, tensor, symm)
UNARY_FUNCTION(symmTensor, tensor, twoSymm)
UNARY_FUNCTION(symmTensor, tensor, devSymm)
UNARY_FUNCTION(symmTensor, tensor, devTwoSymm)
UNARY_FUNCTION(tensor, tensor, skew)
UNARY_FUNCTION(tensor, tensor, dev)
UNARY_FUNCTION(tensor, tensor, dev2)
UNARY_FUNCTION(scalar, tensor, det)
UNARY_FUNCTION(tensor, tensor, cof)
UNARY_FUNCTION(tensor, tensor, inv)
UNARY_FUNCTION(tensor, tensor, pinv)

UNARY_FUNCTION(vector, symmTensor, eigenValues)
UNARY_FUNCTION(tensor, symmTensor, eigenVectors)


// * * * * * * * * * * * * * * * global operators  * * * * * * * * * * * * * //

UNARY_OPERATOR(vector, tensor, *, hdual)
UNARY_OPERATOR(tensor, vector, *, hdual)

BINARY_OPERATOR(vector, vector, tensor, /, divide)
BINARY_TYPE_OPERATOR(vector, vector, tensor, /, divide)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "undefFieldFunctionsM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "tensorFieldTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
