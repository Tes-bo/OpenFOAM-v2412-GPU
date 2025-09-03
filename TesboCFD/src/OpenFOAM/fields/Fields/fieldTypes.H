/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2024 OpenCFD Ltd.
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

Namespace
    Foam::fieldTypes

Description
    Collection of common field types

Note
    The fieldTypes header contains macros for primitive types
    that Fields are instantiated for.

\*---------------------------------------------------------------------------*/

#ifndef Foam_fieldTypes_H
#define Foam_fieldTypes_H

#include "label.H"
#include "scalar.H"
#include "vector.H"
#include "sphericalTensor.H"
#include "symmTensor.H"
#include "tensor.H"
#include "triad.H"
#include "macros.H"
#include "wordList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define CAPITALIZE_bool Bool
#define CAPITALIZE_label Label
#define CAPITALIZE_scalar Scalar
#define CAPITALIZE_complex Complex
#define CAPITALIZE_vector Vector
#define CAPITALIZE_sphericalTensor SphericalTensor
#define CAPITALIZE_symmTensor SymmTensor
#define CAPITALIZE_tensor Tensor

#define FOR_ALL_FIELD_TYPES(Macro, ...)                                        \
    Macro(scalar, __VA_ARGS__)                                                 \
    Macro(vector, __VA_ARGS__)                                                 \
    Macro(sphericalTensor, __VA_ARGS__)                                        \
    Macro(symmTensor, __VA_ARGS__)                                             \
    Macro(tensor, __VA_ARGS__)


/*---------------------------------------------------------------------------*\
                       Namespace fieldTypes Declaration
\*---------------------------------------------------------------------------*/

namespace Foam
{
namespace fieldTypes
{

//- Standard basic field types (label, scalar, vector, tensor, etc)
//  These also correspond to cloud output fields.
extern const wordList basic;


// Commonly used patch field types

//- An \c empty patch field type
inline const char* emptyTypeName_() noexcept { return "empty"; }

//- An \c empty patch field type
extern const word emptyType;

//- A \c calculated patch field type
inline const char* calculatedTypeName_() noexcept { return "calculated"; }

//- A \c calculated patch field type
extern const word calculatedType;

//- A combined \c zero-gradient and \c calculated patch field type
inline const char* extrapolatedCalculatedTypeName_() noexcept
{
    return "extrapolatedCalculated";
}

//- A combined \c zero-gradient and \c calculated patch field type
extern const word extrapolatedCalculatedType;

//- A \c processor patch field type
inline const char* processorTypeName_() noexcept { return "processor"; }

//- A \c processor patch field type
extern const word processorType;

//- A \c zeroGradient patch field type
inline const char* zeroGradientTypeName_() noexcept { return "zeroGradient"; }

//- A \c zeroGradient patch field type
extern const word zeroGradientType;

} // End namespace fieldTypes
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
