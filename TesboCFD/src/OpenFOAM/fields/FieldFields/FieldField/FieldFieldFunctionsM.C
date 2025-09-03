/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2022-2023 OpenCFD Ltd.
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

#include "FieldM.H"
#include "FieldFieldReuseFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define UNARY_FUNCTION(ReturnType, Type1, Func)                                \
                                                                               \
TEMPLATE                                                                       \
void Func                                                                      \
(                                                                              \
    FieldField<Field, ReturnType>& result,                                     \
    const FieldField<Field, Type1>& f1                                         \
)                                                                              \
{                                                                              \
    const label loopLen = (result).size();                                     \
                                                                               \
    for (label i = 0; i < loopLen; ++i)                                        \
    {                                                                          \
        Func(result[i], f1[i]);                                                \
    }                                                                          \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> Func                                        \
(                                                                              \
    const FieldField<Field, Type1>& f1                                         \
)                                                                              \
{                                                                              \
    auto tres = FieldField<Field, ReturnType>::NewCalculatedType(f1);          \
    Func(tres.ref(), f1);                                                      \
    return tres;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> Func                                        \
(                                                                              \
    const tmp<FieldField<Field, Type1>>& tf1                                   \
)                                                                              \
{                                                                              \
    auto tres = reuseTmpFieldField<Field, ReturnType, Type1>::New(tf1);        \
    Func(tres.ref(), tf1());                                                   \
    tf1.clear();                                                               \
    return tres;                                                               \
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define UNARY_OPERATOR(ReturnType, Type1, Op, OpFunc)                          \
                                                                               \
TEMPLATE                                                                       \
void OpFunc                                                                    \
(                                                                              \
    FieldField<Field, ReturnType>& result,                                     \
    const FieldField<Field, Type1>& f1                                         \
)                                                                              \
{                                                                              \
    const label loopLen = (result).size();                                     \
                                                                               \
    for (label i = 0; i < loopLen; ++i)                                        \
    {                                                                          \
        OpFunc(result[i], f1[i]);                                              \
    }                                                                          \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> operator Op                                 \
(                                                                              \
    const FieldField<Field, Type1>& f1                                         \
)                                                                              \
{                                                                              \
    auto tres = FieldField<Field, ReturnType>::NewCalculatedType(f1);          \
    OpFunc(tres.ref(), f1);                                                    \
    return tres;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> operator Op                                 \
(                                                                              \
    const tmp<FieldField<Field, Type1>>& tf1                                   \
)                                                                              \
{                                                                              \
    auto tres = reuseTmpFieldField<Field, ReturnType, Type1>::New(tf1);        \
    OpFunc(tres.ref(), tf1());                                                 \
    tf1.clear();                                                               \
    return tres;                                                               \
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define BINARY_FUNCTION(ReturnType, Type1, Type2, Func)                        \
                                                                               \
TEMPLATE                                                                       \
void Func                                                                      \
(                                                                              \
    FieldField<Field, ReturnType>& result,                                     \
    const FieldField<Field, Type1>& f1,                                        \
    const FieldField<Field, Type2>& f2                                         \
)                                                                              \
{                                                                              \
    const label loopLen = (result).size();                                     \
                                                                               \
    for (label i = 0; i < loopLen; ++i)                                        \
    {                                                                          \
        Func(result[i], f1[i], f2[i]);                                         \
    }                                                                          \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> Func                                        \
(                                                                              \
    const FieldField<Field, Type1>& f1,                                        \
    const FieldField<Field, Type2>& f2                                         \
)                                                                              \
{                                                                              \
    auto tres = FieldField<Field, ReturnType>::NewCalculatedType(f1);          \
    Func(tres.ref(), f1, f2);                                                  \
    return tres;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> Func                                        \
(                                                                              \
    const FieldField<Field, Type1>& f1,                                        \
    const tmp<FieldField<Field, Type2>>& tf2                                   \
)                                                                              \
{                                                                              \
    auto tres = reuseTmpFieldField<Field, ReturnType, Type2>::New(tf2);        \
    Func(tres.ref(), f1, tf2());                                               \
    tf2.clear();                                                               \
    return tres;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> Func                                        \
(                                                                              \
    const tmp<FieldField<Field, Type1>>& tf1,                                  \
    const FieldField<Field, Type2>& f2                                         \
)                                                                              \
{                                                                              \
    auto tres = reuseTmpFieldField<Field, ReturnType, Type1>::New(tf1);        \
    Func(tres.ref(), tf1(), f2);                                               \
    tf1.clear();                                                               \
    return tres;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> Func                                        \
(                                                                              \
    const tmp<FieldField<Field, Type1>>& tf1,                                  \
    const tmp<FieldField<Field, Type2>>& tf2                                   \
)                                                                              \
{                                                                              \
    auto tres                                                                  \
    (                                                                          \
        reuseTmpTmpFieldField<Field, ReturnType, Type1, Type1, Type2>::        \
            New(tf1, tf2)                                                      \
    );                                                                         \
    Func(tres.ref(), tf1(), tf2());                                            \
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
    FieldField<Field, ReturnType>& result,                                     \
    const Type1& s1,                                                           \
    const FieldField<Field, Type2>& f2                                         \
)                                                                              \
{                                                                              \
    const label loopLen = (result).size();                                     \
                                                                               \
    for (label i = 0; i < loopLen; ++i)                                        \
    {                                                                          \
        Func(result[i], s1, f2[i]);                                            \
    }                                                                          \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> Func                                        \
(                                                                              \
    const Type1& s1,                                                           \
    const FieldField<Field, Type2>& f2                                         \
)                                                                              \
{                                                                              \
    auto tres = FieldField<Field, ReturnType>::NewCalculatedType(f2);          \
    Func(tres.ref(), s1, f2);                                                  \
    return tres;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> Func                                        \
(                                                                              \
    const Type1& s1,                                                           \
    const tmp<FieldField<Field, Type2>>& tf2                                   \
)                                                                              \
{                                                                              \
    auto tres = reuseTmpFieldField<Field, ReturnType, Type2>::New(tf2);        \
    Func(tres.ref(), s1, tf2());                                               \
    tf2.clear();                                                               \
    return tres;                                                               \
}


#define BINARY_TYPE_FUNCTION_FS(ReturnType, Type1, Type2, Func)                \
                                                                               \
TEMPLATE                                                                       \
void Func                                                                      \
(                                                                              \
    FieldField<Field, ReturnType>& result,                                     \
    const FieldField<Field, Type1>& f1,                                        \
    const Type2& s2                                                            \
)                                                                              \
{                                                                              \
    const label loopLen = (result).size();                                     \
                                                                               \
    for (label i = 0; i < loopLen; ++i)                                        \
    {                                                                          \
        Func(result[i], f1[i], s2);                                            \
    }                                                                          \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> Func                                        \
(                                                                              \
    const FieldField<Field, Type1>& f1,                                        \
    const Type2& s2                                                            \
)                                                                              \
{                                                                              \
    auto tres = FieldField<Field, ReturnType>::NewCalculatedType(f1);          \
    Func(tres.ref(), f1, s2);                                                  \
    return tres;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> Func                                        \
(                                                                              \
    const tmp<FieldField<Field, Type1>>& tf1,                                  \
    const Type2& s2                                                            \
)                                                                              \
{                                                                              \
    auto tres = reuseTmpFieldField<Field, ReturnType, Type1>::New(tf1);        \
    Func(tres.ref(), tf1(), s2);                                               \
    tf1.clear();                                                               \
    return tres;                                                               \
}


#define BINARY_TYPE_FUNCTION(ReturnType, Type1, Type2, Func)                   \
    BINARY_TYPE_FUNCTION_SF(ReturnType, Type1, Type2, Func)                    \
    BINARY_TYPE_FUNCTION_FS(ReturnType, Type1, Type2, Func)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define BINARY_OPERATOR(ReturnType, Type1, Type2, Op, OpFunc)                  \
                                                                               \
TEMPLATE                                                                       \
void OpFunc                                                                    \
(                                                                              \
    FieldField<Field, ReturnType>& result,                                     \
    const FieldField<Field, Type1>& f1,                                        \
    const FieldField<Field, Type2>& f2                                         \
)                                                                              \
{                                                                              \
    const label loopLen = (result).size();                                     \
                                                                               \
    for (label i = 0; i < loopLen; ++i)                                        \
    {                                                                          \
        OpFunc(result[i], f1[i], f2[i]);                                       \
    }                                                                          \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> operator Op                                 \
(                                                                              \
    const FieldField<Field, Type1>& f1,                                        \
    const FieldField<Field, Type2>& f2                                         \
)                                                                              \
{                                                                              \
    auto tres = FieldField<Field, ReturnType>::NewCalculatedType(f1);          \
    OpFunc(tres.ref(), f1, f2);                                                \
    return tres;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> operator Op                                 \
(                                                                              \
    const FieldField<Field, Type1>& f1,                                        \
    const tmp<FieldField<Field, Type2>>& tf2                                   \
)                                                                              \
{                                                                              \
    auto tres = reuseTmpFieldField<Field, ReturnType, Type2>::New(tf2);        \
    OpFunc(tres.ref(), f1, tf2());                                             \
    tf2.clear();                                                               \
    return tres;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> operator Op                                 \
(                                                                              \
    const tmp<FieldField<Field, Type1>>& tf1,                                  \
    const FieldField<Field, Type2>& f2                                         \
)                                                                              \
{                                                                              \
    auto tres = reuseTmpFieldField<Field, ReturnType, Type1>::New(tf1);        \
    OpFunc(tres.ref(), tf1(), f2);                                             \
    tf1.clear();                                                               \
    return tres;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> operator Op                                 \
(                                                                              \
    const tmp<FieldField<Field, Type1>>& tf1,                                  \
    const tmp<FieldField<Field, Type2>>& tf2                                   \
)                                                                              \
{                                                                              \
    auto tres                                                                  \
    (                                                                          \
        reuseTmpTmpFieldField<Field, ReturnType, Type1, Type1, Type2>::        \
            New(tf1, tf2)                                                      \
    );                                                                         \
    OpFunc(tres.ref(), tf1(), tf2());                                          \
    tf1.clear();                                                               \
    tf2.clear();                                                               \
    return tres;                                                               \
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define BINARY_TYPE_OPERATOR_SF(ReturnType, Type1, Type2, Op, OpFunc)          \
                                                                               \
TEMPLATE                                                                       \
void OpFunc                                                                    \
(                                                                              \
    FieldField<Field, ReturnType>& result,                                     \
    const Type1& s1,                                                           \
    const FieldField<Field, Type2>& f2                                         \
)                                                                              \
{                                                                              \
    const label loopLen = (result).size();                                     \
                                                                               \
    for (label i = 0; i < loopLen; ++i)                                        \
    {                                                                          \
        OpFunc(result[i], s1, f2[i]);                                          \
    }                                                                          \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> operator Op                                 \
(                                                                              \
    const Type1& s1,                                                           \
    const FieldField<Field, Type2>& f2                                         \
)                                                                              \
{                                                                              \
    auto tres = FieldField<Field, ReturnType>::NewCalculatedType(f2);          \
    OpFunc(tres.ref(), s1, f2);                                                \
    return tres;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> operator Op                                 \
(                                                                              \
    const Type1& s1,                                                           \
    const tmp<FieldField<Field, Type2>>& tf2                                   \
)                                                                              \
{                                                                              \
    auto tres = reuseTmpFieldField<Field, ReturnType, Type2>::New(tf2);        \
    OpFunc(tres.ref(), s1, tf2());                                             \
    tf2.clear();                                                               \
    return tres;                                                               \
}


#define BINARY_TYPE_OPERATOR_FS(ReturnType, Type1, Type2, Op, OpFunc)          \
                                                                               \
TEMPLATE                                                                       \
void OpFunc                                                                    \
(                                                                              \
    FieldField<Field, ReturnType>& result,                                     \
    const FieldField<Field, Type1>& f1,                                        \
    const Type2& s2                                                            \
)                                                                              \
{                                                                              \
    const label loopLen = (result).size();                                     \
                                                                               \
    for (label i = 0; i < loopLen; ++i)                                        \
    {                                                                          \
        OpFunc(result[i], f1[i], s2);                                          \
    }                                                                          \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> operator Op                                 \
(                                                                              \
    const FieldField<Field, Type1>& f1,                                        \
    const Type2& s2                                                            \
)                                                                              \
{                                                                              \
    auto tres = FieldField<Field, ReturnType>::NewCalculatedType(f1);          \
    OpFunc(tres.ref(), f1, s2);                                                \
    return tres;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> operator Op                                 \
(                                                                              \
    const tmp<FieldField<Field, Type1>>& tf1,                                  \
    const Type2& s2                                                            \
)                                                                              \
{                                                                              \
    auto tres = reuseTmpFieldField<Field, ReturnType, Type1>::New(tf1);        \
    OpFunc(tres.ref(), tf1(), s2);                                             \
    tf1.clear();                                                               \
    return tres;                                                               \
}


#define BINARY_TYPE_OPERATOR(ReturnType, Type1, Type2, Op, OpFunc)             \
    BINARY_TYPE_OPERATOR_SF(ReturnType, Type1, Type2, Op, OpFunc)              \
    BINARY_TYPE_OPERATOR_FS(ReturnType, Type1, Type2, Op, OpFunc)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define TERNARY_FUNCTION(ReturnType, Type1, Type2, Type3, Func)                \
                                                                               \
TEMPLATE                                                                       \
void Func                                                                      \
(                                                                              \
    FieldField<Field, ReturnType>& result,                                     \
    const FieldField<Field, Type1>& f1,                                        \
    const FieldField<Field, Type2>& f2,                                        \
    const FieldField<Field, Type3>& f3                                         \
)                                                                              \
{                                                                              \
    const label loopLen = (result).size();                                     \
                                                                               \
    for (label i = 0; i < loopLen; ++i)                                        \
    {                                                                          \
        Func(result[i], f1[i], f2[i], f3[i]);                                  \
    }                                                                          \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> Func                                        \
(                                                                              \
    const FieldField<Field, Type1>& f1,                                        \
    const FieldField<Field, Type2>& f2,                                        \
    const FieldField<Field, Type3>& f3                                         \
)                                                                              \
{                                                                              \
    auto tres = FieldField<Field, ReturnType>::NewCalculatedType(f1);          \
    Func(tres.ref(), f1, f2, f3);                                              \
    return tres;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> Func                                        \
(                                                                              \
    const tmp<FieldField<Field, Type1>>& tf1,                                  \
    const FieldField<Field, Type2>& f2,                                        \
    const FieldField<Field, Type3>& f3                                         \
)                                                                              \
{                                                                              \
    auto tres = reuseTmpFieldField<Field, ReturnType, Type1>::New(tf1);        \
    Func(tres.ref(), tf1(), f2, f3);                                           \
    tf1.clear();                                                               \
    return tres;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> Func                                        \
(                                                                              \
    const FieldField<Field, Type1>& f1,                                        \
    const tmp<FieldField<Field, Type2>>& tf2,                                  \
    const FieldField<Field, Type3>& f3                                         \
)                                                                              \
{                                                                              \
    auto tres = reuseTmpFieldField<Field, ReturnType, Type2>::New(tf2);        \
    Func(tres.ref(), f1, tf2(), f3);                                           \
    tf2.clear();                                                               \
    return tres;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> Func                                        \
(                                                                              \
    const FieldField<Field, Type1>& f1,                                        \
    const FieldField<Field, Type2>& f2,                                        \
    const tmp<FieldField<Field, Type3>>& tf3                                   \
)                                                                              \
{                                                                              \
    auto tres = reuseTmpFieldField<Field, ReturnType, Type3>::New(tf3);        \
    Func(tres.ref(), f1, f2, tf3());                                           \
    tf3.clear();                                                               \
    return tres;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> Func                                        \
(                                                                              \
    const tmp<FieldField<Field, Type1>>& tf1,                                  \
    const tmp<FieldField<Field, Type2>>& tf2,                                  \
    const FieldField<Field, Type3>& f3                                         \
)                                                                              \
{                                                                              \
    tmp<FieldField<Field, ReturnType>> tres                                    \
    (                                                                          \
        reuseTmpTmpFieldField<Field, ReturnType, Type1, Type1, Type2>::        \
            New(tf1, tf2)                                                      \
    );                                                                         \
    Func(tres.ref(), tf1(), tf2(), f3);                                        \
    tf1.clear();                                                               \
    tf2.clear();                                                               \
    return tres;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> Func                                        \
(                                                                              \
    const tmp<FieldField<Field, Type1>>& tf1,                                  \
    const FieldField<Field, Type2>& f2,                                        \
    const tmp<FieldField<Field, Type3>>& tf3                                   \
)                                                                              \
{                                                                              \
    tmp<FieldField<Field, ReturnType>> tres                                    \
    (                                                                          \
        reuseTmpTmpFieldField<Field, ReturnType, Type1, Type1, Type3>::        \
            New(tf1, tf3)                                                      \
    );                                                                         \
    Func(tres.ref(), tf1(), f2, tf3());                                        \
    tf1.clear();                                                               \
    tf3.clear();                                                               \
    return tres;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> Func                                        \
(                                                                              \
    const FieldField<Field, Type1>& f1,                                        \
    const tmp<FieldField<Field, Type2>>& tf2,                                  \
    const tmp<FieldField<Field, Type3>>& tf3                                   \
)                                                                              \
{                                                                              \
    tmp<FieldField<Field, ReturnType>> tres                                    \
    (                                                                          \
        reuseTmpTmpFieldField<Field, ReturnType, Type2, Type2, Type3>::        \
            New(tf2, tf3)                                                      \
    );                                                                         \
    Func(tres.ref(), f1, tf2(), tf3());                                        \
    tf2.clear();                                                               \
    tf3.clear();                                                               \
    return tres;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> Func                                        \
(                                                                              \
    const tmp<FieldField<Field, Type1>>& tf1,                                  \
    const tmp<FieldField<Field, Type2>>& tf2,                                  \
    const tmp<FieldField<Field, Type3>>& tf3                                   \
)                                                                              \
{                                                                              \
    /* TBD: check all three types? */                                          \
    tmp<FieldField<Field, ReturnType>> tres                                    \
    (                                                                          \
        reuseTmpTmpFieldField<Field, ReturnType, Type1, Type1, Type2>::        \
            New(tf1, tf2)                                                      \
    );                                                                         \
    Func(tres.ref(), tf1(), tf2(), tf3());                                     \
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
    FieldField<Field, ReturnType>& result,                                     \
    const FieldField<Field, Type1>& f1,                                        \
    const FieldField<Field, Type2>& f2,                                        \
    const Type3& s3                                                            \
)                                                                              \
{                                                                              \
    const label loopLen = (result).size();                                     \
                                                                               \
    for (label i = 0; i < loopLen; ++i)                                        \
    {                                                                          \
        Func(result[i], f1[i], f2[i], s3);                                     \
    }                                                                          \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> Func                                        \
(                                                                              \
    const FieldField<Field, Type1>& f1,                                        \
    const FieldField<Field, Type2>& f2,                                        \
    const Type3& s3                                                            \
)                                                                              \
{                                                                              \
    auto tres = FieldField<Field, ReturnType>::NewCalculatedType(f1);          \
    Func(tres.ref(), f1, f2, s3);                                              \
    return tres;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> Func                                        \
(                                                                              \
    const tmp<FieldField<Field, Type1>>& tf1,                                  \
    const FieldField<Field, Type2>& f2,                                        \
    const Type3& s3                                                            \
)                                                                              \
{                                                                              \
    auto tres = reuseTmpFieldField<Field, ReturnType, Type1>::New(tf1);        \
    Func(tres.ref(), tf1(), f2, s3);                                           \
    tf1.clear();                                                               \
    return tres;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> Func                                        \
(                                                                              \
    const FieldField<Field, Type1>& f1,                                        \
    const tmp<FieldField<Field, Type2>>& tf2,                                  \
    const Type3& s3                                                            \
)                                                                              \
{                                                                              \
    auto tres = reuseTmpFieldField<Field, ReturnType, Type2>::New(tf2);        \
    Func(tres.ref(), f1, tf2(), s3);                                           \
    tf2.clear();                                                               \
    return tres;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> Func                                        \
(                                                                              \
    const tmp<FieldField<Field, Type1>>& tf1,                                  \
    const tmp<FieldField<Field, Type2>>& tf2,                                  \
    const Type3& s3                                                            \
)                                                                              \
{                                                                              \
    tmp<FieldField<Field, ReturnType>> tres                                    \
    (                                                                          \
        reuseTmpTmpFieldField<Field, ReturnType, Type1, Type1, Type2>::        \
            New(tf1, tf2)                                                      \
    );                                                                         \
    Func(tres.ref(), tf1(), tf2(), s3);                                        \
    tf1.clear();                                                               \
    tf2.clear();                                                               \
    return tres;                                                               \
}


// ************************************************************************* //
