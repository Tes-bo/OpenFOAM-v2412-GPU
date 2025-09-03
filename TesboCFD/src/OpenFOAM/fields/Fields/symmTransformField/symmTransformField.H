/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2023 OpenCFD Ltd.
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
    Foam::symmTransformField

Description
    Spatial transformation functions for primitive fields.

SourceFiles
    symmTransformField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_symmTransformField_H
#define Foam_symmTransformField_H

#include "symmTransform.H"
#include "symmTensorField.H"
#include "sphericalTensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// transform()

template<class Type>
void transform
(
    Field<Type>& result,
    const symmTensor& rot,
    const Field<Type>& fld
);

template<class Type>
void transform
(
    Field<Type>& result,
    const symmTensorField& rot,
    const Field<Type>& fld
);


template<class Type>
tmp<Field<Type>> transform
(
    const symmTensorField& rot,
    const Field<Type>& fld
);

template<class Type>
tmp<Field<Type>> transform
(
    const symmTensorField& rot,
    const tmp<Field<Type>>& tfld
);

template<class Type>
tmp<Field<Type>> transform
(
    const tmp<symmTensorField>& trot,
    const Field<Type>& tfld
);

template<class Type>
tmp<Field<Type>> transform
(
    const tmp<symmTensorField>& trot,
    const tmp<Field<Type>>& tfld
);


template<class Type>
tmp<Field<Type>> transform
(
    const symmTensor& rot,
    const Field<Type>& fld
);

template<class Type>
tmp<Field<Type>> transform
(
    const symmTensor& rot,
    const tmp<Field<Type>>& tfld
);


// Specializations

template<>
tmp<Field<sphericalTensor>>
transformFieldMask<sphericalTensor>(const symmTensorField&);

template<>
tmp<Field<sphericalTensor>>
transformFieldMask<sphericalTensor>(const tmp<symmTensorField>&);


template<>
tmp<Field<symmTensor>>
transformFieldMask<symmTensor>(const symmTensorField&);

template<>
tmp<Field<symmTensor>>
transformFieldMask<symmTensor>(const tmp<symmTensorField>&);


template<>
tmp<Field<tensor>>
transformFieldMask<tensor>(const symmTensorField&);

template<>
tmp<Field<tensor>>
transformFieldMask<tensor>(const tmp<symmTensorField>&);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "symmTransformField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
