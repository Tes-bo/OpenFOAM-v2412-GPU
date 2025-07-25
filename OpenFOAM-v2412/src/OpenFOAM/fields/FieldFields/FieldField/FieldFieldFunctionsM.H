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

Description
    Macro functions for FieldField\<Type\> algebra.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define UNARY_FUNCTION(ReturnType, Type1, Func)                                \
                                                                               \
TEMPLATE                                                                       \
void Func                                                                      \
(                                                                              \
    FieldField<Field, ReturnType>& result,                                     \
    const FieldField<Field, Type1>& f1                                         \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> Func                                        \
(                                                                              \
    const FieldField<Field, Type1>& f1                                         \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> Func                                        \
(                                                                              \
    const tmp<FieldField<Field, Type1>>& tf1                                   \
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define UNARY_OPERATOR(ReturnType, Type1, Op, OpFunc)                          \
                                                                               \
TEMPLATE                                                                       \
void OpFunc                                                                    \
(                                                                              \
    FieldField<Field, ReturnType>& result,                                     \
    const FieldField<Field, Type1>& f1                                         \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> operator Op                                 \
(                                                                              \
    const FieldField<Field, Type1>& f1                                         \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> operator Op                                 \
(                                                                              \
    const tmp<FieldField<Field, Type1>>& tf1                                   \
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define BINARY_FUNCTION(ReturnType, Type1, Type2, Func)                        \
                                                                               \
TEMPLATE                                                                       \
void Func                                                                      \
(                                                                              \
    FieldField<Field, ReturnType>& result,                                     \
    const FieldField<Field, Type1>& f1,                                        \
    const FieldField<Field, Type2>& f2                                         \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> Func                                        \
(                                                                              \
    const FieldField<Field, Type1>& f1,                                        \
    const FieldField<Field, Type2>& f2                                         \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> Func                                        \
(                                                                              \
    const FieldField<Field, Type1>& f1,                                        \
    const tmp<FieldField<Field, Type2>>& tf2                                   \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> Func                                        \
(                                                                              \
    const tmp<FieldField<Field, Type1>>& tf1,                                  \
    const FieldField<Field, Type2>& f2                                         \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> Func                                        \
(                                                                              \
    const tmp<FieldField<Field, Type1>>& tf1,                                  \
    const tmp<FieldField<Field, Type2>>& tf2                                   \
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define BINARY_TYPE_FUNCTION_SF(ReturnType, Type1, Type2, Func)                \
                                                                               \
TEMPLATE                                                                       \
void Func                                                                      \
(                                                                              \
    FieldField<Field, ReturnType>& result,                                     \
    const Type1& s1,                                                           \
    const FieldField<Field, Type2>& f2                                         \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> Func                                        \
(                                                                              \
    const Type1& s1,                                                           \
    const FieldField<Field, Type1>& f2                                         \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> Func                                        \
(                                                                              \
    const Type1& s1,                                                           \
    const tmp<FieldField<Field, Type1>>& tf2                                   \
);


#define BINARY_TYPE_FUNCTION_FS(ReturnType, Type1, Type2, Func)                \
                                                                               \
TEMPLATE                                                                       \
void Func                                                                      \
(                                                                              \
    FieldField<Field, ReturnType>& result,                                     \
    const FieldField<Field, Type1>& f1,                                        \
    const Type2& s2                                                            \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> Func                                        \
(                                                                              \
    const FieldField<Field, Type1>& f1,                                        \
    const Type2& s2                                                            \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> Func                                        \
(                                                                              \
    const tmp<FieldField<Field, Type1>>& tf1,                                  \
    const Type2& s2                                                            \
);


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
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> operator Op                                 \
(                                                                              \
    const FieldField<Field, Type1>& f1,                                        \
    const FieldField<Field, Type2>& f2                                         \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> operator Op                                 \
(                                                                              \
    const FieldField<Field, Type1>& f1,                                        \
    const tmp<FieldField<Field, Type2>>& tf2                                   \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> operator Op                                 \
(                                                                              \
    const tmp<FieldField<Field, Type1>>& tf1,                                  \
    const FieldField<Field, Type2>& f2                                         \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> operator Op                                 \
(                                                                              \
    const tmp<FieldField<Field, Type1>>& tf1,                                  \
    const tmp<FieldField<Field, Type2>>& tf2                                   \
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define BINARY_TYPE_OPERATOR_SF(ReturnType, Type1, Type2, Op, OpFunc)          \
                                                                               \
TEMPLATE                                                                       \
void OpFunc                                                                    \
(                                                                              \
    FieldField<Field, ReturnType>& result,                                     \
    const Type1& s1,                                                           \
    const FieldField<Field, Type2>& f2                                         \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> operator Op                                 \
(                                                                              \
    const Type1& s1,                                                           \
    const FieldField<Field, Type2>& f2                                         \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> operator Op                                 \
(                                                                              \
    const Type1& s1,                                                           \
    const tmp<FieldField<Field, Type2>>& tf2                                   \
);


#define BINARY_TYPE_OPERATOR_FS(ReturnType, Type1, Type2, Op, OpFunc)          \
                                                                               \
TEMPLATE                                                                       \
void OpFunc                                                                    \
(                                                                              \
    FieldField<Field, ReturnType>& result,                                     \
    const FieldField<Field, Type1>& f1,                                        \
    const Type2& s2                                                            \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> operator Op                                 \
(                                                                              \
    const FieldField<Field, Type1>& f1,                                        \
    const Type2& s2                                                            \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> operator Op                                 \
(                                                                              \
    const tmp<FieldField<Field, Type1>>& tf1,                                  \
    const Type2& s2                                                            \
);


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
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> Func                                        \
(                                                                              \
    const FieldField<Field, Type1>& f1,                                        \
    const FieldField<Field, Type2>& f2,                                        \
    const FieldField<Field, Type3>& f3                                         \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> Func                                        \
(                                                                              \
    const tmp<FieldField<Field, Type1>>& tf1,                                  \
    const FieldField<Field, Type2>& f2,                                        \
    const FieldField<Field, Type3>& f3                                         \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> Func                                        \
(                                                                              \
    const FieldField<Field, Type1>& f1,                                        \
    const tmp<FieldField<Field, Type2>>& tf2,                                  \
    const FieldField<Field, Type3>& f3                                         \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> Func                                        \
(                                                                              \
    const FieldField<Field, Type1>& f1,                                        \
    const FieldField<Field, Type2>& f2,                                        \
    const tmp<FieldField<Field, Type3>>& tf3                                   \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> Func                                        \
(                                                                              \
    const tmp<FieldField<Field, Type1>>& tf1,                                  \
    const tmp<FieldField<Field, Type2>>& tf2,                                  \
    const FieldField<Field, Type3>& f3                                         \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> Func                                        \
(                                                                              \
    const tmp<FieldField<Field, Type1>>& tf1,                                  \
    const FieldField<Field, Type2>& f2,                                        \
    const tmp<FieldField<Field, Type3>>& tf3                                   \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> Func                                        \
(                                                                              \
    const FieldField<Field, Type1>& f1,                                        \
    const tmp<FieldField<Field, Type2>>& tf2,                                  \
    const tmp<FieldField<Field, Type3>>& tf3                                   \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> Func                                        \
(                                                                              \
    const tmp<FieldField<Field, Type1>>& tf1,                                  \
    const tmp<FieldField<Field, Type2>>& tf2,                                  \
    const tmp<FieldField<Field, Type3>>& tf3                                   \
);


#define TERNARY_TYPE_FUNCTION_FFS(ReturnType, Type1, Type2, Type3, Func)       \
                                                                               \
TEMPLATE                                                                       \
void Func                                                                      \
(                                                                              \
    FieldField<Field, ReturnType>& result,                                     \
    const FieldField<Field, Type1>& f1,                                        \
    const FieldField<Field, Type2>& f2,                                        \
    const Type3& s3                                                            \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> Func                                        \
(                                                                              \
    const FieldField<Field, Type1>& f1,                                        \
    const FieldField<Field, Type2>& f2,                                        \
    const Type3& s3                                                            \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> Func                                        \
(                                                                              \
    const tmp<FieldField<Field, Type1>>& tf1,                                  \
    const FieldField<Field, Type2>& f2,                                        \
    const Type3& s3                                                            \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> Func                                        \
(                                                                              \
    const FieldField<Field, Type1>& f1,                                        \
    const tmp<FieldField<Field, Type2>>& tf2,                                  \
    const Type3& s3                                                            \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> Func                                        \
(                                                                              \
    const tmp<FieldField<Field, Type1>>& tf1,                                  \
    const tmp<FieldField<Field, Type2>>& tf2,                                  \
    const Type3& s3                                                            \
);


// ************************************************************************* //
