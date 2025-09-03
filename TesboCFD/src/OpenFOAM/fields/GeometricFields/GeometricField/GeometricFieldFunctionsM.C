/*---------------------------------------------------------------------------* \
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

#include "GeometricFieldReuseFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define UNARY_FUNCTION(ReturnType, Type1, Func, Dfunc)                         \
                                                                               \
TEMPLATE                                                                       \
void Func                                                                      \
(                                                                              \
    GeometricField<ReturnType, PatchField, GeoMesh>& result,                   \
    const GeometricField<Type1, PatchField, GeoMesh>& f1                       \
)                                                                              \
{                                                                              \
    Foam::Func(result.primitiveFieldRef(), f1.primitiveField());               \
    Foam::Func(result.boundaryFieldRef(), f1.boundaryField());                 \
    result.oriented() = f1.oriented();                                         \
    result.correctLocalBoundaryConditions();                                   \
    if (GeometricBoundaryField<ReturnType, PatchField, GeoMesh>::debug)        \
    {                                                                          \
        result.boundaryField().check();                                        \
    }                                                                          \
}                                                                              \
                                                                               \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const GeometricField<Type1, PatchField, GeoMesh>& f1                       \
)                                                                              \
{                                                                              \
    auto tres =                                                                \
        reuseTmpGeometricField<ReturnType, Type1, PatchField, GeoMesh>::New    \
        (                                                                      \
            f1,                                                                \
            #Func "(" + f1.name() + ')',                                       \
            Dfunc(f1.dimensions())                                             \
        );                                                                     \
                                                                               \
    Foam::Func(tres.ref(), f1);                                                \
    return tres;                                                               \
}                                                                              \
                                                                               \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const tmp<GeometricField<Type1, PatchField, GeoMesh>>& tf1                 \
)                                                                              \
{                                                                              \
    const auto& f1 = tf1();                                                    \
                                                                               \
    auto tres =                                                                \
        reuseTmpGeometricField<ReturnType, Type1, PatchField, GeoMesh>::New    \
        (                                                                      \
            tf1,                                                               \
            #Func "(" + f1.name() + ')',                                       \
            Dfunc(f1.dimensions())                                             \
        );                                                                     \
                                                                               \
    Foam::Func(tres.ref(), f1);                                                \
    tf1.clear();                                                               \
    return tres;                                                               \
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define UNARY_OPERATOR(ReturnType, Type1, Op, OpFunc, Dfunc)                   \
                                                                               \
TEMPLATE                                                                       \
void OpFunc                                                                    \
(                                                                              \
    GeometricField<ReturnType, PatchField, GeoMesh>& result,                   \
    const GeometricField<Type1, PatchField, GeoMesh>& f1                       \
)                                                                              \
{                                                                              \
    Foam::OpFunc(result.primitiveFieldRef(), f1.primitiveField());             \
    Foam::OpFunc(result.boundaryFieldRef(), f1.boundaryField());               \
    result.oriented() = f1.oriented();                                         \
    result.correctLocalBoundaryConditions();                                   \
    if (GeometricBoundaryField<ReturnType, PatchField, GeoMesh>::debug)        \
    {                                                                          \
        result.boundaryField().check();                                        \
    }                                                                          \
}                                                                              \
                                                                               \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> operator Op               \
(                                                                              \
    const GeometricField<Type1, PatchField, GeoMesh>& f1                       \
)                                                                              \
{                                                                              \
    auto tres =                                                                \
        reuseTmpGeometricField<ReturnType, Type1, PatchField, GeoMesh>::New    \
        (                                                                      \
            f1,                                                                \
            #Op + f1.name(),                                                   \
            Dfunc(f1.dimensions())                                             \
        );                                                                     \
                                                                               \
    Foam::OpFunc(tres.ref(), f1);                                              \
    return tres;                                                               \
}                                                                              \
                                                                               \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> operator Op               \
(                                                                              \
    const tmp<GeometricField<Type1, PatchField, GeoMesh>>& tf1                 \
)                                                                              \
{                                                                              \
    const auto& f1 = tf1();                                                    \
                                                                               \
    auto tres =                                                                \
        reuseTmpGeometricField<ReturnType, Type1, PatchField, GeoMesh>::New    \
        (                                                                      \
            tf1,                                                               \
            #Op + f1.name(),                                                   \
            Dfunc(f1.dimensions())                                             \
        );                                                                     \
                                                                               \
    Foam::OpFunc(tres.ref(), f1);                                              \
    tf1.clear();                                                               \
    return tres;                                                               \
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define BINARY_FUNCTION(ReturnType, Type1, Type2, Func)                        \
                                                                               \
TEMPLATE                                                                       \
void Func                                                                      \
(                                                                              \
    GeometricField<ReturnType, PatchField, GeoMesh>& result,                   \
    const GeometricField<Type1, PatchField, GeoMesh>& f1,                      \
    const GeometricField<Type2, PatchField, GeoMesh>& f2                       \
)                                                                              \
{                                                                              \
    Foam::Func                                                                 \
    (                                                                          \
        result.primitiveFieldRef(),                                            \
        f1.primitiveField(),                                                   \
        f2.primitiveField()                                                    \
    );                                                                         \
    Foam::Func                                                                 \
    (                                                                          \
        result.boundaryFieldRef(),                                             \
        f1.boundaryField(),                                                    \
        f2.boundaryField()                                                     \
    );                                                                         \
    result.oriented() = Func(f1.oriented(), f2.oriented());                    \
    result.correctLocalBoundaryConditions();                                   \
    if (GeometricBoundaryField<ReturnType, PatchField, GeoMesh>::debug)        \
    {                                                                          \
        result.boundaryField().check();                                        \
    }                                                                          \
}                                                                              \
                                                                               \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const GeometricField<Type1, PatchField, GeoMesh>& f1,                      \
    const GeometricField<Type2, PatchField, GeoMesh>& f2                       \
)                                                                              \
{                                                                              \
    auto tres =                                                                \
        reuseTmpGeometricField<ReturnType, Type1, PatchField, GeoMesh>::New    \
        (                                                                      \
            f1,                                                                \
            #Func "(" + f1.name() + ',' + f2.name() + ')',                     \
            Func(f1.dimensions(), f2.dimensions())                             \
        );                                                                     \
                                                                               \
    Foam::Func(tres.ref(), f1, f2);                                            \
    return tres;                                                               \
}                                                                              \
                                                                               \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const GeometricField<Type1, PatchField, GeoMesh>& f1,                      \
    const tmp<GeometricField<Type2, PatchField, GeoMesh>>& tf2                 \
)                                                                              \
{                                                                              \
    const auto& f2 = tf2();                                                    \
                                                                               \
    auto tres =                                                                \
        reuseTmpGeometricField<ReturnType, Type2, PatchField, GeoMesh>::New    \
        (                                                                      \
            tf2,                                                               \
            #Func "(" + f1.name() + ',' + f2.name() + ')',                     \
            Func(f1.dimensions(), f2.dimensions())                             \
        );                                                                     \
                                                                               \
    Foam::Func(tres.ref(), f1, f2);                                            \
    tf2.clear();                                                               \
    return tres;                                                               \
}                                                                              \
                                                                               \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const tmp<GeometricField<Type1, PatchField, GeoMesh>>& tf1,                \
    const GeometricField<Type2, PatchField, GeoMesh>& f2                       \
)                                                                              \
{                                                                              \
    const auto& f1 = tf1();                                                    \
                                                                               \
    auto tres =                                                                \
        reuseTmpGeometricField<ReturnType, Type1, PatchField, GeoMesh>::New    \
        (                                                                      \
            tf1,                                                               \
            #Func "(" + f1.name() + ',' + f2.name() + ')',                     \
            Func(f1.dimensions(), f2.dimensions())                             \
        );                                                                     \
                                                                               \
    Foam::Func(tres.ref(), f1, f2);                                            \
    tf1.clear();                                                               \
    return tres;                                                               \
}                                                                              \
                                                                               \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const tmp<GeometricField<Type1, PatchField, GeoMesh>>& tf1,                \
    const tmp<GeometricField<Type2, PatchField, GeoMesh>>& tf2                 \
)                                                                              \
{                                                                              \
    const auto& f1 = tf1();                                                    \
    const auto& f2 = tf2();                                                    \
                                                                               \
    auto tres =                                                                \
        reuseTmpTmpGeometricField                                              \
            <ReturnType, Type1, Type1, Type2, PatchField, GeoMesh>             \
        ::New                                                                  \
        (                                                                      \
            tf1,                                                               \
            tf2,                                                               \
            #Func "(" + f1.name() + ',' + f2.name() + ')',                     \
            Func(f1.dimensions(), f2.dimensions())                             \
        );                                                                     \
                                                                               \
    Foam::Func(tres.ref(), f1, f2);                                            \
    tf1.clear();                                                               \
    tf2.clear();                                                               \
    return tres;                                                               \
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define BINARY_TYPE_FUNCTION_SF(ReturnType, Type1, Type2, Func)                \
                                                                               \
TEMPLATE                                                                       \
void Func                                                                      \
(                                                                              \
    GeometricField<ReturnType, PatchField, GeoMesh>& result,                   \
    const dimensioned<Type1>& dt1,                                             \
    const GeometricField<Type2, PatchField, GeoMesh>& f2                       \
)                                                                              \
{                                                                              \
    Foam::Func(result.primitiveFieldRef(), dt1.value(), f2.primitiveField());  \
    Foam::Func(result.boundaryFieldRef(), dt1.value(), f2.boundaryField());    \
    result.oriented() = f2.oriented();                                         \
    result.correctLocalBoundaryConditions();                                   \
    if (GeometricBoundaryField<ReturnType, PatchField, GeoMesh>::debug)        \
    {                                                                          \
        result.boundaryField().check();                                        \
    }                                                                          \
}                                                                              \
                                                                               \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const dimensioned<Type1>& dt1,                                             \
    const GeometricField<Type2, PatchField, GeoMesh>& f2                       \
)                                                                              \
{                                                                              \
    auto tres =                                                                \
        reuseTmpGeometricField<ReturnType, Type2, PatchField, GeoMesh>::New    \
        (                                                                      \
            f2,                                                                \
            #Func "(" + dt1.name() + ',' + f2.name() + ')',                    \
            Func(dt1.dimensions(), f2.dimensions())                            \
        );                                                                     \
                                                                               \
    Foam::Func(tres.ref(), dt1, f2);                                           \
    return tres;                                                               \
}                                                                              \
                                                                               \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const Type1& s1,                                                           \
    const GeometricField<Type2, PatchField, GeoMesh>& f2                       \
)                                                                              \
{                                                                              \
    return Func(dimensioned<Type1>(s1), f2);                                   \
}                                                                              \
                                                                               \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const dimensioned<Type1>& dt1,                                             \
    const tmp<GeometricField<Type2, PatchField, GeoMesh>>& tf2                 \
)                                                                              \
{                                                                              \
    const auto& f2 = tf2();                                                    \
                                                                               \
    auto tres =                                                                \
        reuseTmpGeometricField<ReturnType, Type2, PatchField, GeoMesh>::New    \
        (                                                                      \
            tf2,                                                               \
            #Func "(" + dt1.name() + ',' + f2.name() + ')',                    \
            Func(dt1.dimensions(), f2.dimensions())                            \
        );                                                                     \
                                                                               \
    Foam::Func(tres.ref(), dt1, f2);                                           \
    tf2.clear();                                                               \
    return tres;                                                               \
}                                                                              \
                                                                               \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const Type1& s1,                                                           \
    const tmp<GeometricField<Type2, PatchField, GeoMesh>>& tf2                 \
)                                                                              \
{                                                                              \
    return Func(dimensioned<Type1>(s1), tf2);                                  \
}


#define BINARY_TYPE_FUNCTION_FS(ReturnType, Type1, Type2, Func)                \
                                                                               \
TEMPLATE                                                                       \
void Func                                                                      \
(                                                                              \
    GeometricField<ReturnType, PatchField, GeoMesh>& result,                   \
    const GeometricField<Type1, PatchField, GeoMesh>& f1,                      \
    const dimensioned<Type2>& dt2                                              \
)                                                                              \
{                                                                              \
    Foam::Func(result.primitiveFieldRef(), f1.primitiveField(), dt2.value());  \
    Foam::Func(result.boundaryFieldRef(), f1.boundaryField(), dt2.value());    \
    result.oriented() = f1.oriented();                                         \
    result.correctLocalBoundaryConditions();                                   \
    if (GeometricBoundaryField<ReturnType, PatchField, GeoMesh>::debug)        \
    {                                                                          \
        result.boundaryField().check();                                        \
    }                                                                          \
}                                                                              \
                                                                               \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const GeometricField<Type1, PatchField, GeoMesh>& f1,                      \
    const dimensioned<Type2>& dt2                                              \
)                                                                              \
{                                                                              \
    auto tres =                                                                \
        reuseTmpGeometricField<ReturnType, Type1, PatchField, GeoMesh>::New    \
        (                                                                      \
            f1,                                                                \
            #Func "(" + f1.name() + ',' + dt2.name() + ')',                    \
            Func(f1.dimensions(), dt2.dimensions())                            \
        );                                                                     \
                                                                               \
    Foam::Func(tres.ref(), f1, dt2);                                           \
    return tres;                                                               \
}                                                                              \
                                                                               \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const GeometricField<Type1, PatchField, GeoMesh>& f1,                      \
    const Type2& s2                                                            \
)                                                                              \
{                                                                              \
    return Func(f1, dimensioned<Type2>(s2));                                   \
}                                                                              \
                                                                               \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const tmp<GeometricField<Type1, PatchField, GeoMesh>>& tf1,                \
    const dimensioned<Type2>& dt2                                              \
)                                                                              \
{                                                                              \
    const auto& f1 = tf1();                                                    \
                                                                               \
    auto tres =                                                                \
        reuseTmpGeometricField<ReturnType, Type1, PatchField, GeoMesh>::New    \
        (                                                                      \
            tf1,                                                               \
            #Func "(" + f1.name() + ',' + dt2.name() + ')',                    \
            Func(f1.dimensions(), dt2.dimensions())                            \
        );                                                                     \
                                                                               \
    Foam::Func(tres.ref(), f1, dt2);                                           \
    tf1.clear();                                                               \
    return tres;                                                               \
}                                                                              \
                                                                               \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const tmp<GeometricField<Type1, PatchField, GeoMesh>>& tf1,                \
    const Type2& s2                                                            \
)                                                                              \
{                                                                              \
    return Func(tf1, dimensioned<Type2>(s2));                                  \
}


#define BINARY_TYPE_FUNCTION(ReturnType, Type1, Type2, Func)                   \
    BINARY_TYPE_FUNCTION_SF(ReturnType, Type1, Type2, Func)                    \
    BINARY_TYPE_FUNCTION_FS(ReturnType, Type1, Type2, Func)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define BINARY_OPERATOR(ReturnType, Type1, Type2, Op, OpName, OpFunc)          \
                                                                               \
TEMPLATE                                                                       \
void OpFunc                                                                    \
(                                                                              \
    GeometricField<ReturnType, PatchField, GeoMesh>& result,                   \
    const GeometricField<Type1, PatchField, GeoMesh>& f1,                      \
    const GeometricField<Type2, PatchField, GeoMesh>& f2                       \
)                                                                              \
{                                                                              \
    Foam::OpFunc                                                               \
    (                                                                          \
        result.primitiveFieldRef(),                                            \
        f1.primitiveField(),                                                   \
        f2.primitiveField()                                                    \
    );                                                                         \
    Foam::OpFunc                                                               \
    (                                                                          \
        result.boundaryFieldRef(),                                             \
        f1.boundaryField(),                                                    \
        f2.boundaryField()                                                     \
    );                                                                         \
    result.oriented() = (f1.oriented() Op f2.oriented());                      \
    result.correctLocalBoundaryConditions();                                   \
    if (GeometricBoundaryField<ReturnType, PatchField, GeoMesh>::debug)        \
    {                                                                          \
        result.boundaryField().check();                                        \
    }                                                                          \
}                                                                              \
                                                                               \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> operator Op               \
(                                                                              \
    const GeometricField<Type1, PatchField, GeoMesh>& f1,                      \
    const GeometricField<Type2, PatchField, GeoMesh>& f2                       \
)                                                                              \
{                                                                              \
    auto tres =                                                                \
        reuseTmpGeometricField<ReturnType, Type1, PatchField, GeoMesh>::New    \
        (                                                                      \
            f1,                                                                \
            '(' + f1.name() + OpName + f2.name() + ')',                        \
            (f1.dimensions() Op f2.dimensions())                               \
        );                                                                     \
                                                                               \
    Foam::OpFunc(tres.ref(), f1, f2);                                          \
    return tres;                                                               \
}                                                                              \
                                                                               \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> operator Op               \
(                                                                              \
    const GeometricField<Type1, PatchField, GeoMesh>& f1,                      \
    const tmp<GeometricField<Type2, PatchField, GeoMesh>>& tf2                 \
)                                                                              \
{                                                                              \
    const auto& f2 = tf2();                                                    \
                                                                               \
    auto tres =                                                                \
        reuseTmpGeometricField<ReturnType, Type2, PatchField, GeoMesh>::New    \
        (                                                                      \
            tf2,                                                               \
            '(' + f1.name() + OpName + f2.name() + ')',                        \
            (f1.dimensions() Op f2.dimensions())                               \
        );                                                                     \
                                                                               \
    Foam::OpFunc(tres.ref(), f1, f2);                                          \
    tf2.clear();                                                               \
    return tres;                                                               \
}                                                                              \
                                                                               \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> operator Op               \
(                                                                              \
    const tmp<GeometricField<Type1, PatchField, GeoMesh>>& tf1,                \
    const GeometricField<Type2, PatchField, GeoMesh>& f2                       \
)                                                                              \
{                                                                              \
    const auto& f1 = tf1();                                                    \
                                                                               \
    auto tres =                                                                \
        reuseTmpGeometricField<ReturnType, Type1, PatchField, GeoMesh>::New    \
        (                                                                      \
            tf1,                                                               \
            '(' + f1.name() + OpName + f2.name() + ')',                        \
            (f1.dimensions() Op f2.dimensions())                               \
        );                                                                     \
                                                                               \
    Foam::OpFunc(tres.ref(), f1, f2);                                          \
    tf1.clear();                                                               \
    return tres;                                                               \
}                                                                              \
                                                                               \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> operator Op               \
(                                                                              \
    const tmp<GeometricField<Type1, PatchField, GeoMesh>>& tf1,                \
    const tmp<GeometricField<Type2, PatchField, GeoMesh>>& tf2                 \
)                                                                              \
{                                                                              \
    const auto& f1 = tf1();                                                    \
    const auto& f2 = tf2();                                                    \
                                                                               \
    auto tres =                                                                \
        reuseTmpTmpGeometricField                                              \
            <ReturnType, Type1, Type1, Type2, PatchField, GeoMesh>::New        \
        (                                                                      \
            tf1,                                                               \
            tf2,                                                               \
            '(' + f1.name() + OpName + f2.name() + ')',                        \
            (f1.dimensions() Op f2.dimensions())                               \
        );                                                                     \
                                                                               \
    Foam::OpFunc(tres.ref(), f1, f2);                                          \
    tf1.clear();                                                               \
    tf2.clear();                                                               \
    return tres;                                                               \
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define BINARY_TYPE_OPERATOR_SF(ReturnType, Type1, Type2, Op, OpName, OpFunc)  \
                                                                               \
TEMPLATE                                                                       \
void OpFunc                                                                    \
(                                                                              \
    GeometricField<ReturnType, PatchField, GeoMesh>& result,                   \
    const dimensioned<Type1>& dt1,                                             \
    const GeometricField<Type2, PatchField, GeoMesh>& f2                       \
)                                                                              \
{                                                                              \
    Foam::OpFunc(result.primitiveFieldRef(), dt1.value(), f2.primitiveField());\
    Foam::OpFunc(result.boundaryFieldRef(), dt1.value(), f2.boundaryField());  \
    result.oriented() = f2.oriented();                                         \
    result.correctLocalBoundaryConditions();                                   \
    if (GeometricBoundaryField<ReturnType, PatchField, GeoMesh>::debug)        \
    {                                                                          \
        result.boundaryField().check();                                        \
    }                                                                          \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> operator Op               \
(                                                                              \
    const dimensioned<Type1>& dt1,                                             \
    const GeometricField<Type2, PatchField, GeoMesh>& f2                       \
)                                                                              \
{                                                                              \
    auto tres =                                                                \
        reuseTmpGeometricField<ReturnType, Type2, PatchField, GeoMesh>::New    \
        (                                                                      \
            f2,                                                                \
            '(' + dt1.name() + OpName + f2.name() + ')',                       \
            (dt1.dimensions() Op f2.dimensions())                              \
        );                                                                     \
                                                                               \
    Foam::OpFunc(tres.ref(), dt1, f2);                                         \
    return tres;                                                               \
}                                                                              \
                                                                               \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> operator Op               \
(                                                                              \
    const Type1& t1,                                                           \
    const GeometricField<Type2, PatchField, GeoMesh>& f2                       \
)                                                                              \
{                                                                              \
    return dimensioned<Type1>(t1) Op f2;                                       \
}                                                                              \
                                                                               \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> operator Op               \
(                                                                              \
    const dimensioned<Type1>& dt1,                                             \
    const tmp<GeometricField<Type2, PatchField, GeoMesh>>& tf2                 \
)                                                                              \
{                                                                              \
    const auto& f2 = tf2();                                                    \
                                                                               \
    auto tres =                                                                \
        reuseTmpGeometricField<ReturnType, Type2, PatchField, GeoMesh>::New    \
        (                                                                      \
            tf2,                                                               \
            '(' + dt1.name() + OpName + f2.name() + ')',                       \
            (dt1.dimensions() Op f2.dimensions())                              \
        );                                                                     \
                                                                               \
    Foam::OpFunc(tres.ref(), dt1, f2);                                         \
    tf2.clear();                                                               \
    return tres;                                                               \
}                                                                              \
                                                                               \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> operator Op               \
(                                                                              \
    const Type1& s1,                                                           \
    const tmp<GeometricField<Type2, PatchField, GeoMesh>>& tf2                 \
)                                                                              \
{                                                                              \
    return dimensioned<Type1>(s1) Op tf2;                                      \
}


#define BINARY_TYPE_OPERATOR_FS(ReturnType, Type1, Type2, Op, OpName, OpFunc)  \
                                                                               \
TEMPLATE                                                                       \
void OpFunc                                                                    \
(                                                                              \
    GeometricField<ReturnType, PatchField, GeoMesh>& result,                   \
    const GeometricField<Type1, PatchField, GeoMesh>& f1,                      \
    const dimensioned<Type2>& dt2                                              \
)                                                                              \
{                                                                              \
    Foam::OpFunc(result.primitiveFieldRef(), f1.primitiveField(), dt2.value());\
    Foam::OpFunc(result.boundaryFieldRef(), f1.boundaryField(), dt2.value());  \
    result.oriented() = f1.oriented();                                         \
    result.correctLocalBoundaryConditions();                                   \
    if (GeometricBoundaryField<ReturnType, PatchField, GeoMesh>::debug)        \
    {                                                                          \
        result.boundaryField().check();                                        \
    }                                                                          \
}                                                                              \
                                                                               \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> operator Op               \
(                                                                              \
    const GeometricField<Type1, PatchField, GeoMesh>& f1,                      \
    const dimensioned<Type2>& dt2                                              \
)                                                                              \
{                                                                              \
    auto tres =                                                                \
        reuseTmpGeometricField<ReturnType, Type1, PatchField, GeoMesh>::New    \
        (                                                                      \
            f1,                                                                \
            '(' + f1.name() + OpName + dt2.name() + ')',                       \
            (f1.dimensions() Op dt2.dimensions())                              \
        );                                                                     \
                                                                               \
    Foam::OpFunc(tres.ref(), f1, dt2);                                         \
    return tres;                                                               \
}                                                                              \
                                                                               \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> operator Op               \
(                                                                              \
    const GeometricField<Type1, PatchField, GeoMesh>& f1,                      \
    const Type2& s2                                                            \
)                                                                              \
{                                                                              \
    return f1 Op dimensioned<Type2>(s2);                                       \
}                                                                              \
                                                                               \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> operator Op               \
(                                                                              \
    const tmp<GeometricField<Type1, PatchField, GeoMesh>>& tf1,                \
    const dimensioned<Type2>& dt2                                              \
)                                                                              \
{                                                                              \
    const auto& f1 = tf1();                                                    \
                                                                               \
    auto tres =                                                                \
        reuseTmpGeometricField<ReturnType, Type1, PatchField, GeoMesh>::New    \
        (                                                                      \
            tf1,                                                               \
            '(' + f1.name() + OpName + dt2.name() + ')',                       \
            (f1.dimensions() Op dt2.dimensions())                              \
        );                                                                     \
                                                                               \
    Foam::OpFunc(tres.ref(), f1, dt2);                                         \
    tf1.clear();                                                               \
    return tres;                                                               \
}                                                                              \
                                                                               \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> operator Op               \
(                                                                              \
    const tmp<GeometricField<Type1, PatchField, GeoMesh>>& tf1,                \
    const Type2& s2                                                            \
)                                                                              \
{                                                                              \
    return tf1 Op dimensioned<Type2>(s2);                                      \
}


#define BINARY_TYPE_OPERATOR(ReturnType, Type1, Type2, Op, OpName, OpFunc)     \
    BINARY_TYPE_OPERATOR_SF(ReturnType, Type1, Type2, Op, OpName, OpFunc)      \
    BINARY_TYPE_OPERATOR_FS(ReturnType, Type1, Type2, Op, OpName, OpFunc)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define TERNARY_FUNCTION(ReturnType, Type1, Type2, Type3, Func)                \
                                                                               \
TEMPLATE                                                                       \
void Func                                                                      \
(                                                                              \
    GeometricField<ReturnType, PatchField, GeoMesh>& result,                   \
    const GeometricField<Type1, PatchField, GeoMesh>& f1,                      \
    const GeometricField<Type2, PatchField, GeoMesh>& f2,                      \
    const GeometricField<Type3, PatchField, GeoMesh>& f3                       \
)                                                                              \
{                                                                              \
    Foam::Func                                                                 \
    (                                                                          \
        result.primitiveFieldRef(),                                            \
        f1.primitiveField(),                                                   \
        f2.primitiveField(),                                                   \
        f3.primitiveField()                                                    \
    );                                                                         \
    Foam::Func                                                                 \
    (                                                                          \
        result.boundaryFieldRef(),                                             \
        f1.boundaryField(),                                                    \
        f2.boundaryField(),                                                    \
        f3.boundaryField()                                                     \
    );                                                                         \
    result.oriented() = Func(f1.oriented(), f2.oriented());                    \
    result.correctLocalBoundaryConditions();                                   \
    if (GeometricBoundaryField<ReturnType, PatchField, GeoMesh>::debug)        \
    {                                                                          \
        result.boundaryField().check();                                        \
    }                                                                          \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const GeometricField<Type1, PatchField, GeoMesh>& f1,                      \
    const GeometricField<Type2, PatchField, GeoMesh>& f2,                      \
    const GeometricField<Type3, PatchField, GeoMesh>& f3                       \
)                                                                              \
{                                                                              \
    auto tres =                                                                \
        reuseTmpGeometricField<ReturnType, Type1, PatchField, GeoMesh>::New    \
        (                                                                      \
            f1,                                                                \
            #Func "(" + f1.name() + ',' + f2.name() + ',' + f3.name() + ')',   \
            Func(f1.dimensions(), f2.dimensions())                             \
        );                                                                     \
                                                                               \
    Foam::Func(tres.ref(), f1, f2, f3);                                        \
    return tres;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const tmp<GeometricField<Type1, PatchField, GeoMesh>>& tf1,                \
    const GeometricField<Type2, PatchField, GeoMesh>& f2,                      \
    const GeometricField<Type3, PatchField, GeoMesh>& f3                       \
)                                                                              \
{                                                                              \
    const auto& f1 = tf1();                                                    \
                                                                               \
    auto tres =                                                                \
        reuseTmpGeometricField<ReturnType, Type1, PatchField, GeoMesh>::New    \
        (                                                                      \
            f1,                                                                \
            #Func "(" + f1.name() + ',' + f2.name() + ',' + f3.name() + ')',   \
            Func(f1.dimensions(), f2.dimensions())                             \
        );                                                                     \
                                                                               \
    Foam::Func(tres.ref(), f1, f2, f3);                                        \
    tf1.clear();                                                               \
    return tres;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const GeometricField<Type1, PatchField, GeoMesh>& f1,                      \
    const tmp<GeometricField<Type2, PatchField, GeoMesh>>& tf2,                \
    const GeometricField<Type3, PatchField, GeoMesh>& f3                       \
)                                                                              \
{                                                                              \
    const auto& f2 = tf2();                                                    \
                                                                               \
    auto tres =                                                                \
        reuseTmpGeometricField<ReturnType, Type2, PatchField, GeoMesh>::New    \
        (                                                                      \
            f2,                                                                \
            #Func "(" + f1.name() + ',' + f2.name() + ',' + f3.name() + ')',   \
            Func(f1.dimensions(), f2.dimensions())                             \
        );                                                                     \
                                                                               \
    Foam::Func(tres.ref(), f1, f2, f3);                                        \
    tf2.clear();                                                               \
    return tres;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const GeometricField<Type1, PatchField, GeoMesh>& f1,                      \
    const GeometricField<Type2, PatchField, GeoMesh>& f2,                      \
    const tmp<GeometricField<Type3, PatchField, GeoMesh>>& tf3                 \
)                                                                              \
{                                                                              \
    const auto& f3 = tf3();                                                    \
                                                                               \
    auto tres =                                                                \
        reuseTmpGeometricField<ReturnType, Type3, PatchField, GeoMesh>::New    \
        (                                                                      \
            f3,                                                                \
            #Func "(" + f1.name() + ',' + f2.name() + ',' + f3.name() + ')',   \
            Func(f1.dimensions(), f2.dimensions())                             \
        );                                                                     \
                                                                               \
    Foam::Func(tres.ref(), f1, f2, f3);                                        \
    tf3.clear();                                                               \
    return tres;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const tmp<GeometricField<Type1, PatchField, GeoMesh>>& tf1,                \
    const tmp<GeometricField<Type2, PatchField, GeoMesh>>& tf2,                \
    const GeometricField<Type3, PatchField, GeoMesh>& f3                       \
)                                                                              \
{                                                                              \
    const auto& f1 = tf1();                                                    \
    const auto& f2 = tf2();                                                    \
                                                                               \
    auto tres =                                                                \
        reuseTmpTmpGeometricField                                              \
            <ReturnType, Type1, Type1, Type2, PatchField, GeoMesh>             \
        ::New                                                                  \
        (                                                                      \
            tf1,                                                               \
            tf2,                                                               \
            #Func "(" + f1.name() + ',' + f2.name() + ',' + f3.name() + ')',   \
            Func(f1.dimensions(), f2.dimensions())                             \
        );                                                                     \
                                                                               \
    Foam::Func(tres.ref(), f1, f2, f3);                                        \
    tf1.clear();                                                               \
    tf2.clear();                                                               \
    return tres;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const tmp<GeometricField<Type1, PatchField, GeoMesh>>& tf1,                \
    const GeometricField<Type2, PatchField, GeoMesh>& f2,                      \
    const tmp<GeometricField<Type3, PatchField, GeoMesh>>& tf3                 \
)                                                                              \
{                                                                              \
    const auto& f1 = tf1();                                                    \
    const auto& f3 = tf3();                                                    \
                                                                               \
    auto tres =                                                                \
        reuseTmpTmpGeometricField                                              \
            <ReturnType, Type1, Type1, Type3, PatchField, GeoMesh>             \
        ::New                                                                  \
        (                                                                      \
            tf1,                                                               \
            tf3,                                                               \
            #Func "(" + f1.name() + ',' + f2.name() + ',' + f3.name() + ')',   \
            Func(f1.dimensions(), f2.dimensions())                             \
        );                                                                     \
                                                                               \
    Foam::Func(tres.ref(), f1, f2, f3);                                        \
    tf1.clear();                                                               \
    tf3.clear();                                                               \
    return tres;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const GeometricField<Type1, PatchField, GeoMesh>& f1,                      \
    const tmp<GeometricField<Type2, PatchField, GeoMesh>>& tf2,                \
    const tmp<GeometricField<Type3, PatchField, GeoMesh>>& tf3                 \
)                                                                              \
{                                                                              \
    const auto& f2 = tf2();                                                    \
    const auto& f3 = tf3();                                                    \
                                                                               \
    auto tres =                                                                \
        reuseTmpTmpGeometricField                                              \
            <ReturnType, Type2, Type2, Type3, PatchField, GeoMesh>             \
        ::New                                                                  \
        (                                                                      \
            tf2,                                                               \
            tf3,                                                               \
            #Func "(" + f1.name() + ',' + f2.name() + ',' + f3.name() + ')',   \
            Func(f1.dimensions(), f2.dimensions())                             \
        );                                                                     \
                                                                               \
    Foam::Func(tres.ref(), f1, f2, f3);                                        \
    tf2.clear();                                                               \
    tf3.clear();                                                               \
    return tres;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const tmp<GeometricField<Type1, PatchField, GeoMesh>>& tf1,                \
    const tmp<GeometricField<Type2, PatchField, GeoMesh>>& tf2,                \
    const tmp<GeometricField<Type3, PatchField, GeoMesh>>& tf3                 \
)                                                                              \
{                                                                              \
    const auto& f1 = tf1();                                                    \
    const auto& f2 = tf2();                                                    \
    const auto& f3 = tf3();                                                    \
                                                                               \
    auto tres =                                                                \
        reuseTmpTmpGeometricField                                              \
            <ReturnType, Type1, Type1, Type2, PatchField, GeoMesh>             \
        ::New                                                                  \
        (                                                                      \
            tf1,                                                               \
            tf2,                                                               \
            #Func "(" + f1.name() + ',' + f2.name() + ',' + f3.name() + ')',   \
            Func(f1.dimensions(), f2.dimensions())                             \
        );                                                                     \
                                                                               \
    Foam::Func(tres.ref(), f1, f2, f3);                                        \
    tf1.clear();                                                               \
    tf2.clear();                                                               \
    tf3.clear();                                                               \
    return tres;                                                               \
}


#define TERNARY_TYPE_FUNCTION_FFS(ReturnType, Type1, Type2, Type3, Func)       \
                                                                               \
TEMPLATE                                                                       \
void Func                                                                      \
(                                                                              \
    GeometricField<ReturnType, PatchField, GeoMesh>& result,                   \
    const GeometricField<Type1, PatchField, GeoMesh>& f1,                      \
    const GeometricField<Type2, PatchField, GeoMesh>& f2,                      \
    const dimensioned<Type3>& dt3                                              \
)                                                                              \
{                                                                              \
    Foam::Func                                                                 \
    (                                                                          \
        result.primitiveFieldRef(),                                            \
        f1.primitiveField(),                                                   \
        f2.primitiveField(),                                                   \
        dt3.value()                                                            \
    );                                                                         \
    Foam::Func                                                                 \
    (                                                                          \
        result.boundaryFieldRef(),                                             \
        f1.boundaryField(),                                                    \
        f2.boundaryField(),                                                    \
        dt3.value()                                                            \
    );                                                                         \
    result.oriented() = Func(f1.oriented(), f2.oriented());                    \
    result.correctLocalBoundaryConditions();                                   \
    if (GeometricBoundaryField<ReturnType, PatchField, GeoMesh>::debug)        \
    {                                                                          \
        result.boundaryField().check();                                        \
    }                                                                          \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const GeometricField<Type1, PatchField, GeoMesh>& f1,                      \
    const GeometricField<Type2, PatchField, GeoMesh>& f2,                      \
    const dimensioned<Type3>& dt3                                              \
)                                                                              \
{                                                                              \
    auto tres =                                                                \
        reuseTmpGeometricField<ReturnType, Type1, PatchField, GeoMesh>::New    \
        (                                                                      \
            f1,                                                                \
            #Func "(" + f1.name() + ',' + f2.name() + ',' + dt3.name() + ')',  \
            Func(f1.dimensions(), f2.dimensions())                             \
        );                                                                     \
                                                                               \
    Foam::Func(tres.ref(), f1, f2, dt3.value());                               \
    return tres;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const GeometricField<Type1, PatchField, GeoMesh>& f1,                      \
    const GeometricField<Type2, PatchField, GeoMesh>& f2,                      \
    const Type3& s3                                                            \
)                                                                              \
{                                                                              \
    return Foam::Func(f1, f2, dimensioned<Type>(s3));                          \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const tmp<GeometricField<Type1, PatchField, GeoMesh>>& tf1,                \
    const GeometricField<Type2, PatchField, GeoMesh>& f2,                      \
    const dimensioned<Type3>& dt3                                              \
)                                                                              \
{                                                                              \
    const auto& f1 = tf1();                                                    \
                                                                               \
    auto tres =                                                                \
        reuseTmpGeometricField<ReturnType, Type1, PatchField, GeoMesh>::New    \
        (                                                                      \
            f1,                                                                \
            #Func "(" + f1.name() + ',' + f2.name() + ',' + dt3.name() + ')',  \
            Func(f1.dimensions(), f2.dimensions())                             \
        );                                                                     \
                                                                               \
    Foam::Func(tres.ref(), f1, f2, dt3.value());                               \
    tf1.clear();                                                               \
    return tres;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const tmp<GeometricField<Type1, PatchField, GeoMesh>>& tf1,                \
    const GeometricField<Type2, PatchField, GeoMesh>& f2,                      \
    const Type3& s3                                                            \
)                                                                              \
{                                                                              \
    return Foam::Func(tf1, f2, dimensioned<Type>(s3));                         \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const GeometricField<Type1, PatchField, GeoMesh>& f1,                      \
    const tmp<GeometricField<Type2, PatchField, GeoMesh>>& tf2,                \
    const dimensioned<Type3>& dt3                                              \
)                                                                              \
{                                                                              \
    const auto& f2 = tf2();                                                    \
                                                                               \
    auto tres =                                                                \
        reuseTmpGeometricField<ReturnType, Type2, PatchField, GeoMesh>::New    \
        (                                                                      \
            tf2,                                                               \
            #Func "(" + f1.name() + ',' + f2.name() + ',' + dt3.name() + ')',  \
            Func(f1.dimensions(), f2.dimensions())                             \
        );                                                                     \
                                                                               \
    Foam::Func(tres.ref(), f1, f2, dt3.value());                               \
    tf2.clear();                                                               \
    return tres;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const GeometricField<Type1, PatchField, GeoMesh>& f1,                      \
    const tmp<GeometricField<Type2, PatchField, GeoMesh>>& tf2,                \
    const Type3& s3                                                            \
)                                                                              \
{                                                                              \
    return Foam::Func(f1, tf2, dimensioned<Type>(s3));                         \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const tmp<GeometricField<Type1, PatchField, GeoMesh>>& tf1,                \
    const tmp<GeometricField<Type2, PatchField, GeoMesh>>& tf2,                \
    const dimensioned<Type3>& dt3                                              \
)                                                                              \
{                                                                              \
    const auto& f1 = tf1();                                                    \
    const auto& f2 = tf2();                                                    \
                                                                               \
    auto tres =                                                                \
        reuseTmpTmpGeometricField                                              \
            <ReturnType, Type1, Type1, Type2, PatchField, GeoMesh>             \
        ::New                                                                  \
        (                                                                      \
            tf1,                                                               \
            tf2,                                                               \
            #Func "(" + f1.name() + ',' + f2.name() + ',' + dt3.name() + ')',  \
            Func(f1.dimensions(), f2.dimensions())                             \
        );                                                                     \
                                                                               \
    Foam::Func(tres.ref(), f1, f2, dt3.value());                               \
    tf1.clear();                                                               \
    tf2.clear();                                                               \
    return tres;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const tmp<GeometricField<Type1, PatchField, GeoMesh>>& tf1,                \
    const tmp<GeometricField<Type2, PatchField, GeoMesh>>& tf2,                \
    const Type3& s3                                                            \
)                                                                              \
{                                                                              \
    return Foam::Func(tf1, tf2, dimensioned<Type>(s3));                        \
}


// ************************************************************************* //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
