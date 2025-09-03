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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define UNARY_FUNCTION(ReturnType, Type1, Func, Dfunc)                         \
                                                                               \
TEMPLATE                                                                       \
void Func                                                                      \
(                                                                              \
    GeometricField<ReturnType, PatchField, GeoMesh>& result,                   \
    const GeometricField<Type1, PatchField, GeoMesh>& f1                       \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const GeometricField<Type1, PatchField, GeoMesh>& f1                       \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const tmp<GeometricField<Type1, PatchField, GeoMesh>>& tf1                 \
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define UNARY_OPERATOR(ReturnType, Type1, Op, OpFunc, Dfunc)                   \
                                                                               \
TEMPLATE                                                                       \
void OpFunc                                                                    \
(                                                                              \
    GeometricField<ReturnType, PatchField, GeoMesh>& result,                   \
    const GeometricField<Type1, PatchField, GeoMesh>& f1                       \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> operator Op               \
(                                                                              \
    const GeometricField<Type1, PatchField, GeoMesh>& f1                       \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> operator Op               \
(                                                                              \
    const tmp<GeometricField<Type1, PatchField, GeoMesh>>& tf1                 \
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define BINARY_FUNCTION(ReturnType, Type1, Type2, Func)                        \
                                                                               \
TEMPLATE                                                                       \
void Func                                                                      \
(                                                                              \
    GeometricField<ReturnType, PatchField, GeoMesh>& result,                   \
    const GeometricField<Type1, PatchField, GeoMesh>& f1,                      \
    const GeometricField<Type2, PatchField, GeoMesh>& f2                       \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const GeometricField<Type1, PatchField, GeoMesh>& df1,                     \
    const GeometricField<Type2, PatchField, GeoMesh>& df2                      \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const GeometricField<Type1, PatchField, GeoMesh>& f1,                      \
    const tmp<GeometricField<Type2, PatchField, GeoMesh>>& tf2                 \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const tmp<GeometricField<Type1, PatchField, GeoMesh>>& tf1,                \
    const GeometricField<Type2, PatchField, GeoMesh>& f2                       \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const tmp<GeometricField<Type1, PatchField, GeoMesh>>& tf1,                \
    const tmp<GeometricField<Type2, PatchField, GeoMesh>>& tf2                 \
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define BINARY_TYPE_FUNCTION_SF(ReturnType, Type1, Type2, Func)                \
                                                                               \
TEMPLATE                                                                       \
void Func                                                                      \
(                                                                              \
    GeometricField<ReturnType, PatchField, GeoMesh>& result,                   \
    const dimensioned<Type1>& dt1,                                             \
    const GeometricField<Type2, PatchField, GeoMesh>& f2                       \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const dimensioned<Type1>& dt1,                                             \
    const GeometricField<Type2, PatchField, GeoMesh>& f2                       \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const Type1& s1,                                                           \
    const GeometricField<Type2, PatchField, GeoMesh>& f2                       \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const dimensioned<Type1>& dt1,                                             \
    const tmp<GeometricField<Type2, PatchField, GeoMesh>>& tf2                 \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const Type1& s1,                                                           \
    const tmp<GeometricField<Type2, PatchField, GeoMesh>>& tf2                 \
);


#define BINARY_TYPE_FUNCTION_FS(ReturnType, Type1, Type2, Func)                \
                                                                               \
TEMPLATE                                                                       \
void Func                                                                      \
(                                                                              \
    GeometricField<ReturnType, PatchField, GeoMesh>& result,                   \
    const GeometricField<Type1, PatchField, GeoMesh>& f1,                      \
    const dimensioned<Type2>& dt2                                              \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const GeometricField<Type1, PatchField, GeoMesh>& f1,                      \
    const dimensioned<Type2>& dt2                                              \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const GeometricField<Type1, PatchField, GeoMesh>& f1,                      \
    const Type2& s2                                                            \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const tmp<GeometricField<Type1, PatchField, GeoMesh>>& tf1,                \
    const dimensioned<Type2>& dt2                                              \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const tmp<GeometricField<Type1, PatchField, GeoMesh>>& tf1,                \
    const Type2& s2                                                            \
);


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
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> operator Op               \
(                                                                              \
    const GeometricField<Type1, PatchField, GeoMesh>& f1,                      \
    const GeometricField<Type2, PatchField, GeoMesh>& f2                       \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> operator Op               \
(                                                                              \
    const GeometricField<Type1, PatchField, GeoMesh>& f1,                      \
    const tmp<GeometricField<Type2, PatchField, GeoMesh>>& tf2                 \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> operator Op               \
(                                                                              \
    const tmp<GeometricField<Type1, PatchField, GeoMesh>>& tf1,                \
    const GeometricField<Type2, PatchField, GeoMesh>& f2                       \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> operator Op               \
(                                                                              \
    const tmp<GeometricField<Type1, PatchField, GeoMesh>>& tf1,                \
    const tmp<GeometricField<Type2, PatchField, GeoMesh>>& tf2                 \
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define BINARY_TYPE_OPERATOR_SF(ReturnType, Type1, Type2, Op, OpName, OpFunc)  \
                                                                               \
TEMPLATE                                                                       \
void OpFunc                                                                    \
(                                                                              \
    GeometricField<ReturnType, PatchField, GeoMesh>& result,                   \
    const dimensioned<Type1>& dt1,                                             \
    const GeometricField<Type2, PatchField, GeoMesh>& f2                       \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> operator Op               \
(                                                                              \
    const dimensioned<Type1>& dt1,                                             \
    const GeometricField<Type2, PatchField, GeoMesh>& f2                       \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> operator Op               \
(                                                                              \
    const Type1& s1,                                                           \
    const GeometricField<Type2, PatchField, GeoMesh>& f2                       \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> operator Op               \
(                                                                              \
    const dimensioned<Type1>& dt1,                                             \
    const tmp<GeometricField<Type2, PatchField, GeoMesh>>& tf2                 \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> operator Op               \
(                                                                              \
    const Type1& s1,                                                           \
    const tmp<GeometricField<Type2, PatchField, GeoMesh>>& tf2                 \
);


#define BINARY_TYPE_OPERATOR_FS(ReturnType, Type1, Type2, Op, OpName, OpFunc)  \
                                                                               \
TEMPLATE                                                                       \
void OpFunc                                                                    \
(                                                                              \
    GeometricField<ReturnType, PatchField, GeoMesh>& result,                   \
    const GeometricField<Type1, PatchField, GeoMesh>& f1,                      \
    const dimensioned<Type2>& dt2                                              \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> operator Op               \
(                                                                              \
    const GeometricField<Type1, PatchField, GeoMesh>& f1,                      \
    const dimensioned<Type2>& dt2                                              \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> operator Op               \
(                                                                              \
    const GeometricField<Type1, PatchField, GeoMesh>& f1,                      \
    const Type2& s2                                                            \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> operator Op               \
(                                                                              \
    const tmp<GeometricField<Type1, PatchField, GeoMesh>>& tf1,                \
    const dimensioned<Type2>& dt2                                              \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> operator Op               \
(                                                                              \
    const tmp<GeometricField<Type1, PatchField, GeoMesh>>& tf1,                \
    const Type2& s2                                                            \
);


#define BINARY_TYPE_OPERATOR(ReturnType, Type1, Type2, Op, OpName, OpFunc)     \
    BINARY_TYPE_OPERATOR_SF(ReturnType, Type1, Type2, Op, OpName, OpFunc)      \
    BINARY_TYPE_OPERATOR_FS(ReturnType, Type1, Type2, Op, OpName, OpFunc)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define TERNARY_FUNCTION(ReturnType, Type1, Type2, Type3, Func)                \
                                                                               \
TEMPLATE                                                                       \
void Func                                                                      \
(                                                                              \
    GeometricField<Type1, PatchField, GeoMesh>& result,                        \
    const GeometricField<Type1, PatchField, GeoMesh>& f1,                      \
    const GeometricField<Type2, PatchField, GeoMesh>& f2,                      \
    const GeometricField<Type3, PatchField, GeoMesh>& f3                       \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const GeometricField<Type1, PatchField, GeoMesh>& f1,                      \
    const GeometricField<Type2, PatchField, GeoMesh>& f2,                      \
    const GeometricField<Type3, PatchField, GeoMesh>& f3                       \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const tmp<GeometricField<Type1, PatchField, GeoMesh>>& tf1,                \
    const GeometricField<Type2, PatchField, GeoMesh>& f2,                      \
    const GeometricField<Type3, PatchField, GeoMesh>& f3                       \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const GeometricField<Type1, PatchField, GeoMesh>& f1,                      \
    const tmp<GeometricField<Type2, PatchField, GeoMesh>>& tf2,                \
    const GeometricField<Type3, PatchField, GeoMesh>& f3                       \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const GeometricField<Type1, PatchField, GeoMesh>& f1,                      \
    const GeometricField<Type2, PatchField, GeoMesh>& f2,                      \
    const tmp<GeometricField<Type3, PatchField, GeoMesh>>& tf3                 \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const tmp<GeometricField<Type1, PatchField, GeoMesh>>& tf1,                \
    const tmp<GeometricField<Type2, PatchField, GeoMesh>>& tf2,                \
    const GeometricField<Type3, PatchField, GeoMesh>& f3                       \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const tmp<GeometricField<Type1, PatchField, GeoMesh>>& tf1,                \
    const GeometricField<Type2, PatchField, GeoMesh>& f2,                      \
    const tmp<GeometricField<Type3, PatchField, GeoMesh>>& tf3                 \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const GeometricField<Type1, PatchField, GeoMesh>& f1,                      \
    const tmp<GeometricField<Type2, PatchField, GeoMesh>>& tf2,                \
    const tmp<GeometricField<Type3, PatchField, GeoMesh>>& tf3                 \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const tmp<GeometricField<Type1, PatchField, GeoMesh>>& tf1,                \
    const tmp<GeometricField<Type2, PatchField, GeoMesh>>& tf2,                \
    const tmp<GeometricField<Type3, PatchField, GeoMesh>>& tf3                 \
);


#define TERNARY_TYPE_FUNCTION_FFS(ReturnType, Type1, Type2, Type3, Func)       \
                                                                               \
TEMPLATE                                                                       \
void Func                                                                      \
(                                                                              \
    GeometricField<Type1, PatchField, GeoMesh>& result,                        \
    const GeometricField<Type1, PatchField, GeoMesh>& f1,                      \
    const GeometricField<Type2, PatchField, GeoMesh>& f2,                      \
    const dimensioned<Type3>& dt3                                              \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const GeometricField<Type1, PatchField, GeoMesh>& f1,                      \
    const GeometricField<Type2, PatchField, GeoMesh>& f2,                      \
    const dimensioned<Type3>& dt3                                              \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const GeometricField<Type1, PatchField, GeoMesh>& f1,                      \
    const GeometricField<Type2, PatchField, GeoMesh>& f2,                      \
    const Type3& s3                                                            \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const tmp<GeometricField<Type1, PatchField, GeoMesh>>& tf1,                \
    const GeometricField<Type2, PatchField, GeoMesh>& f2,                      \
    const dimensioned<Type3>& dt3                                              \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const tmp<GeometricField<Type1, PatchField, GeoMesh>>& tf1,                \
    const GeometricField<Type2, PatchField, GeoMesh>& f2,                      \
    const Type3& s3                                                            \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const GeometricField<Type1, PatchField, GeoMesh>& f1,                      \
    const tmp<GeometricField<Type2, PatchField, GeoMesh>>& tf2,                \
    const dimensioned<Type3>& dt3                                              \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const GeometricField<Type1, PatchField, GeoMesh>& f1,                      \
    const tmp<GeometricField<Type2, PatchField, GeoMesh>>& tf2,                \
    const Type3& s3                                                            \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const tmp<GeometricField<Type1, PatchField, GeoMesh>>& tf1,                \
    const tmp<GeometricField<Type2, PatchField, GeoMesh>>& tf2,                \
    const dimensioned<Type3>& dt3                                              \
);                                                                             \
                                                                               \
TEMPLATE                                                                       \
tmp<GeometricField<ReturnType, PatchField, GeoMesh>> Func                      \
(                                                                              \
    const tmp<GeometricField<Type1, PatchField, GeoMesh>>& tf1,                \
    const tmp<GeometricField<Type2, PatchField, GeoMesh>>& tf2,                \
    const Type3& s3                                                            \
);


// ************************************************************************* //
