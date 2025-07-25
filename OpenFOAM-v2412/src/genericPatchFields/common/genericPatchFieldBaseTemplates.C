/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2023 OpenCFD Ltd.
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

#include "genericPatchFieldBase.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class MapperType>
void Foam::genericPatchFieldBase::mapGeneric
(
    const genericPatchFieldBase& rhs,
    const MapperType& mapper
)
{
    #undef  doLocalCode
    #define doLocalCode(ValueType, Member)                                    \
    forAllIters(rhs.Member, iter)                                             \
    {                                                                         \
        this->Member.insert                                                   \
        (                                                                     \
            iter.key(),                                                       \
            autoPtr<Field<ValueType>>::New(*iter.val(), mapper)               \
        );                                                                    \
    }

    //doLocalCode(label, labelFields_);
    doLocalCode(scalar, scalarFields_);
    doLocalCode(vector, vectorFields_);
    doLocalCode(sphericalTensor, sphTensorFields_);
    doLocalCode(symmTensor, symmTensorFields_);
    doLocalCode(tensor, tensorFields_);
    #undef doLocalCode
}


template<class MapperType>
void Foam::genericPatchFieldBase::autoMapGeneric
(
    const MapperType& mapper
)
{
    #undef  doLocalCode
    #define doLocalCode(ValueType, Member)                                    \
    forAllIters(this->Member, iter)                                           \
    {                                                                         \
        iter.val()->autoMap(mapper);                                          \
    }

    //doLocalCode(label, labelFields_);
    doLocalCode(scalar, scalarFields_);
    doLocalCode(vector, vectorFields_);
    doLocalCode(sphericalTensor, sphTensorFields_);
    doLocalCode(symmTensor, symmTensorFields_);
    doLocalCode(tensor, tensorFields_);
    #undef doLocalCode
}


// ************************************************************************* //
