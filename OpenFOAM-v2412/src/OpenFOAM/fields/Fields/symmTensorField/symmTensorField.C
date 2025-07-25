/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "symmTensorField.H"
#include "transformField.H"

#define TEMPLATE
#include "FieldFunctionsM.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

UNARY_FUNCTION(symmTensor, vector, sqr)
UNARY_FUNCTION(symmTensor, symmTensor, innerSqr)

UNARY_FUNCTION(scalar, symmTensor, tr)
UNARY_FUNCTION(sphericalTensor, symmTensor, sph)
UNARY_FUNCTION(symmTensor, symmTensor, symm)
UNARY_FUNCTION(symmTensor, symmTensor, twoSymm)
UNARY_FUNCTION(symmTensor, symmTensor, dev)
UNARY_FUNCTION(symmTensor, symmTensor, dev2)
UNARY_FUNCTION(scalar, symmTensor, det)
UNARY_FUNCTION(symmTensor, symmTensor, cof)

void inv(Field<symmTensor>& result, const UList<symmTensor>& f1)
{
    // With 'failsafe' invert
    // std::transform
    TFOR_ALL_F_OP_F_FUNC(symmTensor, result, =, symmTensor, f1, safeInv)
}

tmp<symmTensorField> inv(const UList<symmTensor>& tf)
{
    auto tresult = tmp<symmTensorField>::New(tf.size());
    inv(tresult.ref(), tf);
    return tresult;
}

tmp<symmTensorField> inv(const tmp<symmTensorField>& tf)
{
    tmp<symmTensorField> tresult = New(tf);
    inv(tresult.ref(), tf());
    tf.clear();
    return tresult;
}

UNARY_FUNCTION(symmTensor, symmTensor, pinv)


template<>
tmp<Field<symmTensor>> transformFieldMask<symmTensor>
(
    const tensorField& tf
)
{
    return symm(tf);
}

template<>
tmp<Field<symmTensor>> transformFieldMask<symmTensor>
(
    const tmp<tensorField>& ttf
)
{
    tmp<Field<symmTensor>> ret = transformFieldMask<symmTensor>(ttf());
    ttf.clear();
    return ret;
}


template<>
tmp<Field<symmTensor>> transformFieldMask<symmTensor>
(
    const symmTensorField& stf
)
{
    return stf;
}

template<>
tmp<Field<symmTensor>> transformFieldMask<symmTensor>
(
    const tmp<symmTensorField>& tstf
)
{
    return tstf;
}


// * * * * * * * * * * * * * * * global operators  * * * * * * * * * * * * * //

UNARY_OPERATOR(vector, symmTensor, *, hdual)

BINARY_OPERATOR(tensor, symmTensor, symmTensor, &, dot)
BINARY_TYPE_OPERATOR(tensor, symmTensor, symmTensor, &, dot)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "undefFieldFunctionsM.H"

// ************************************************************************* //
