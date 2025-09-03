/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

Description
    Declaration macros for Field\<Type\> algebra.

\*---------------------------------------------------------------------------*/

#ifndef Foam_FieldM_H
#define Foam_FieldM_H

#include "error.H"
#include "ListLoopM.H"  // For list access macros

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef FULLDEBUG

template<class Type1, class Type2>
void checkFields
(
    const UList<Type1>& f1,
    const UList<Type2>& f2,
    const char* op
)
{
    if (f1.size() != f2.size())
    {
        FatalErrorInFunction
            << " Field<"<<pTraits<Type1>::typeName<<"> f1("<<f1.size()<<')'
            << " and Field<"<<pTraits<Type2>::typeName<<"> f2("<<f2.size()<<')'
            << endl
            << " for operation " << op
            << abort(FatalError);
    }
}

template<class Type1, class Type2, class Type3>
void checkFields
(
    const UList<Type1>& f1,
    const UList<Type2>& f2,
    const UList<Type3>& f3,
    const char* op
)
{
    if (f1.size() != f2.size() || f1.size() != f3.size())
    {
        FatalErrorInFunction
            << " Field<"<<pTraits<Type1>::typeName<<"> f1("<<f1.size()<<')'
            << ", Field<"<<pTraits<Type2>::typeName<<"> f2("<<f2.size()<<')'
            << " and Field<"<<pTraits<Type3>::typeName<<"> f3("<<f3.size()<<')'
            << endl
            << "    for operation " << op
            << abort(FatalError);
    }
}

template<class Type1, class Type2, class Type3, class Type4>
void checkFields
(
    const UList<Type1>& f1,
    const UList<Type2>& f2,
    const UList<Type3>& f3,
    const UList<Type4>& f4,
    const char* op
)
{
    if
    (
        f1.size() != f2.size()
     || f1.size() != f3.size()
     || f1.size() != f4.size()
    )
    {
        FatalErrorInFunction
            << " Field<"<<pTraits<Type1>::typeName<<"> f1("<<f1.size()<<')'
            << ", Field<"<<pTraits<Type2>::typeName<<"> f2("<<f2.size()<<')'
            << ", Field<"<<pTraits<Type3>::typeName<<"> f3("<<f3.size()<<')'
            << " and Field<"<<pTraits<Type4>::typeName<<"> f4("<<f4.size()<<')'
            << endl
            << "    for operation " << op
            << abort(FatalError);
    }
}

#else

template<class Type1, class Type2>
void checkFields
(
    const UList<Type1>&,
    const UList<Type2>&,
    const char*
)
{}

template<class Type1, class Type2, class Type3>
void checkFields
(
    const UList<Type1>&,
    const UList<Type2>&,
    const UList<Type3>&,
    const char*
)
{}

template<class Type1, class Type2, class Type3, class Type4>
void checkFields
(
    const UList<Type1>&,
    const UList<Type2>&,
    const UList<Type3>&,
    const UList<Type4>&,
    const char*
)
{}

#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Unary Free Function : f1 OP Func(f2)

#define TFOR_ALL_F_OP_FUNC_F(typeF1, f1, OP, FUNC, typeF2, f2)                 \
{                                                                              \
    /* Check fields have same size */                                          \
    checkFields(f1, f2, "f1 " #OP " " #FUNC "(f2)");                           \
                                                                               \
    /* Field access */                                                         \
    List_ACCESS(typeF1, f1, f1P);                                              \
    List_CONST_ACCESS(typeF2, f2, f2P);                                        \
                                                                               \
    /* Loop: f1 OP FUNC(f2) */                                                 \
    const label loop_len = (f1).size();                                        \
                                                                               \
    /* pragmas... */                                                           \
    for (label i = 0; i < loop_len; ++i)                                       \
    {                                                                          \
        (f1P[i]) OP FUNC(f2P[i]);                                              \
    }                                                                          \
}


// Nullary Member Function : f1 OP f2.FUNC()

#define TFOR_ALL_F_OP_F_FUNC(typeF1, f1, OP, typeF2, f2, FUNC)                 \
{                                                                              \
    /* Check fields have same size */                                          \
    checkFields(f1, f2, "f1 " #OP " f2" #FUNC);                                \
                                                                               \
    /* Field access */                                                         \
    List_ACCESS(typeF1, f1, f1P);                                              \
    List_CONST_ACCESS(typeF2, f2, f2P);                                        \
                                                                               \
    /* Loop: f1 OP f2.FUNC() */                                                \
    const label loop_len = (f1).size();                                        \
                                                                               \
    /* pragmas... */                                                           \
    for (label i = 0; i < loop_len; ++i)                                       \
    {                                                                          \
        (f1P[i]) OP (f2P[i]).FUNC();                                           \
    }                                                                          \
}


// Binary Free Function : f1 OP FUNC(f2, f3)

#define TFOR_ALL_F_OP_FUNC_F_F(typeF1, f1, OP, FUNC, typeF2, f2, typeF3, f3)   \
{                                                                              \
    /* Check fields have same size */                                          \
    checkFields(f1, f2, f3, "f1 " #OP " " #FUNC "(f2, f3)");                   \
                                                                               \
    /* Field access */                                                         \
    List_ACCESS(typeF1, f1, f1P);                                              \
    List_CONST_ACCESS(typeF2, f2, f2P);                                        \
    List_CONST_ACCESS(typeF3, f3, f3P);                                        \
                                                                               \
    /* Loop: f1 OP FUNC(f2, f3) */                                             \
    const label loop_len = (f1).size();                                        \
                                                                               \
    /* pragmas... */                                                           \
    for (label i = 0; i < loop_len; ++i)                                       \
    {                                                                          \
        (f1P[i]) OP FUNC((f2P[i]), (f3P[i]));                                  \
    }                                                                          \
}


// [reduction] Binary Free Function : s OP FUNC(f1, f2)

#define TFOR_ALL_S_OP_FUNC_F_F(typeS, s, OP, FUNC, typeF1, f1, typeF2, f2)     \
{                                                                              \
    /* Check fields have same size */                                          \
    checkFields(f1, f2, "s " #OP " " #FUNC "(f1, f2)");                        \
                                                                               \
    /* Field access */                                                         \
    List_CONST_ACCESS(typeF1, f1, f1P);                                        \
    List_CONST_ACCESS(typeF2, f2, f2P);                                        \
                                                                               \
    /* Loop: s OP FUNC(f1, f2) */                                              \
    const label loop_len = (f1).size();                                        \
                                                                               \
    /* pragmas, reduction... */                                                \
    for (label i = 0; i < loop_len; ++i)                                       \
    {                                                                          \
        (s) OP FUNC((f1P[i]), (f2P[i]));                                       \
    }                                                                          \
}


// Binary Free Function : f1 OP FUNC(f2, s)

#define TFOR_ALL_F_OP_FUNC_F_S(typeF1, f1, OP, FUNC, typeF2, f2, typeS, s)     \
{                                                                              \
    /* Check fields have same size */                                          \
    checkFields(f1, f2, "f1 " #OP " " #FUNC "(f2, s)");                        \
                                                                               \
    /* Field access */                                                         \
    List_ACCESS(typeF1, f1, f1P);                                              \
    List_CONST_ACCESS(typeF2, f2, f2P);                                        \
                                                                               \
    /* Loop: f1 OP FUNC(f2, s) */                                              \
    const label loop_len = (f1).size();                                        \
                                                                               \
    /* pragmas... */                                                           \
    for (label i = 0; i < loop_len; ++i)                                       \
    {                                                                          \
        (f1P[i]) OP FUNC((f2P[i]), (s));                                       \
    }                                                                          \
}


// [reduction] Binary Free Function : s1 OP FUNC(f, s2)

#define TFOR_ALL_S_OP_FUNC_F_S(typeS1, s1, OP, FUNC, typeF, f, typeS2, s2)     \
{                                                                              \
    /* Field access */                                                         \
    List_CONST_ACCESS(typeF, f, fP);                                           \
                                                                               \
    /* Loop: s1 OP FUNC(f, s2) */                                              \
    const label loop_len = (f).size();                                         \
                                                                               \
    /* pragmas, reduction... */                                                \
    for (label i = 0; i < loop_len; ++i)                                       \
    {                                                                          \
        (s1) OP FUNC((fP[i]), (s2));                                           \
    }                                                                          \
}


// Binary Free Function : f1 OP FUNC(s, f2)

#define TFOR_ALL_F_OP_FUNC_S_F(typeF1, f1, OP, FUNC, typeS, s, typeF2, f2)     \
{                                                                              \
    /* Check fields have same size */                                          \
    checkFields(f1, f2, "f1 " #OP " " #FUNC "(s, f2)");                        \
                                                                               \
    /* Field access */                                                         \
    List_ACCESS(typeF1, f1, f1P);                                              \
    List_CONST_ACCESS(typeF2, f2, f2P);                                        \
                                                                               \
    /* Loop: f1 OP1 f2 OP2 f3 */                                               \
    const label loop_len = (f1).size();                                        \
                                                                               \
    /* pragmas... */                                                           \
    for (label i = 0; i < loop_len; ++i)                                       \
    {                                                                          \
        (f1P[i]) OP FUNC((s), (f2P[i]));                                       \
    }                                                                          \
}


// Binary Free Function : f1 OP FUNC(s1, s2)

#define TFOR_ALL_F_OP_FUNC_S_S(typeF1, f1, OP, FUNC, typeS1, s1, typeS2, s2)   \
{                                                                              \
    /* Field access */                                                         \
    List_ACCESS(typeF1, f1, f1P);                                              \
                                                                               \
    /* Loop: f1 OP FUNC(s1, s2) */                                             \
    const label loop_len = (f1).size();                                        \
                                                                               \
    /* pragmas... */                                                           \
    for (label i = 0; i < loop_len; ++i)                                       \
    {                                                                          \
        (f1P[i]) OP FUNC((s1), (s2));                                          \
    }                                                                          \
}


// Unary Member Function : f1 OP f2 FUNC(s)

#define TFOR_ALL_F_OP_F_FUNC_S(typeF1, f1, OP, typeF2, f2, FUNC, typeS, s)     \
{                                                                              \
    /* Check fields have same size */                                          \
    checkFields(f1, f2, "f1 " #OP " f2 " #FUNC "(s)");                         \
                                                                               \
    /* Field access */                                                         \
    List_ACCESS(typeF1, f1, f1P);                                              \
    List_CONST_ACCESS(typeF2, f2, f2P);                                        \
                                                                               \
    /* Loop: f1 OP f2 FUNC(s) */                                               \
    const label loop_len = (f1).size();                                        \
                                                                               \
    /* pragmas... */                                                           \
    for (label i = 0; i < loop_len; ++i)                                       \
    {                                                                          \
        (f1P[i]) OP (f2P[i]) FUNC((s));                                        \
    }                                                                          \
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Ternary Free Function : f1 OP FUNC(f2, f3, f4)

#define TFOR_ALL_F_OP_FUNC_F_F_F\
(typeF1, f1, OP, FUNC, typeF2, f2, typeF3, f3, typeF4, f4)                     \
{                                                                              \
    /* Check fields have same size */                                          \
    checkFields(f1, f2, f3, f4, "f1 " #OP " " #FUNC "(f2, f3, f4)");           \
                                                                               \
    /* Field access */                                                         \
    List_ACCESS(typeF1, f1, f1P);                                              \
    List_CONST_ACCESS(typeF2, f2, f2P);                                        \
    List_CONST_ACCESS(typeF3, f3, f3P);                                        \
    List_CONST_ACCESS(typeF4, f4, f4P);                                        \
                                                                               \
    /* Loop: f1 OP FUNC(f2, f3, f4) */                                         \
    const label loop_len = (f1).size();                                        \
                                                                               \
    /* pragmas... */                                                           \
    for (label i = 0; i < loop_len; ++i)                                       \
    {                                                                          \
        (f1P[i]) OP FUNC((f2P[i]), (f3P[i]), (f4P[i]));                        \
    }                                                                          \
}

// Ternary Free Function : f1 OP FUNC(f2, f3, s4)

#define TFOR_ALL_F_OP_FUNC_F_F_S\
(typeF1, f1, OP, FUNC, typeF2, f2, typeF3, f3, typeF4, s4)                     \
{                                                                              \
    /* Check fields have same size */                                          \
    checkFields(f1, f2, f3, "f1 " #OP " " #FUNC "(f2, f3, s)");                \
                                                                               \
    /* Field access */                                                         \
    List_ACCESS(typeF1, f1, f1P);                                              \
    List_CONST_ACCESS(typeF2, f2, f2P);                                        \
    List_CONST_ACCESS(typeF3, f3, f3P);                                        \
                                                                               \
    /* Loop: f1 OP FUNC(f2, f3, s4) */                                         \
    const label loop_len = (f1).size();                                        \
                                                                               \
    /* pragmas... */                                                           \
    for (label i = 0; i < loop_len; ++i)                                       \
    {                                                                          \
        (f1P[i]) OP FUNC((f2P[i]), (f3P[i]), (s4));                            \
    }                                                                          \
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Member operator : this field f1 OP1 f2 OP2 f3

#define TFOR_ALL_F_OP_F_OP_F(typeF1, f1, OP1, typeF2, f2, OP2, typeF3, f3)     \
{                                                                              \
    /* Check fields have same size */                                          \
    checkFields(f1, f2, f3, "f1 " #OP1 " f2 " #OP2 " f3");                     \
                                                                               \
    /* Field access */                                                         \
    List_ACCESS(typeF1, f1, f1P);                                              \
    List_CONST_ACCESS(typeF2, f2, f2P);                                        \
    List_CONST_ACCESS(typeF3, f3, f3P);                                        \
                                                                               \
    /* Loop: f1 OP1 f2 OP2 f3 */                                               \
    const label loop_len = (f1).size();                                        \
                                                                               \
    /* pragmas... */                                                           \
    for (label i = 0; i < loop_len; ++i)                                       \
    {                                                                          \
        (f1P[i]) OP1 (f2P[i]) OP2 (f3P[i]);                                    \
    }                                                                          \
}


// Member operator : this field f1 OP1 s OP2 f2

#define TFOR_ALL_F_OP_S_OP_F(typeF1, f1, OP1, typeS, s, OP2, typeF2, f2)       \
{                                                                              \
    /* Check fields have same size */                                          \
    checkFields(f1, f2, "f1 " #OP1 " s " #OP2 " f2");                          \
                                                                               \
    /* Field access */                                                         \
    List_ACCESS(typeF1, f1, f1P);                                              \
    List_CONST_ACCESS(typeF2, f2, f2P);                                        \
                                                                               \
    /* Loop: f1 OP1 s OP2 f2 */                                                \
    const label loop_len = (f1).size();                                        \
                                                                               \
    /* pragmas... */                                                           \
    for (label i = 0; i < loop_len; ++i)                                       \
    {                                                                          \
        (f1P[i]) OP1 (s) OP2 (f2P[i]);                                         \
    }                                                                          \
}


// Member operator : this field f1 OP1 f2 OP2 s

#define TFOR_ALL_F_OP_F_OP_S(typeF1, f1, OP1, typeF2, f2, OP2, typeS, s)       \
{                                                                              \
    /* Check fields have same size */                                          \
    checkFields(f1, f2, "f1 " #OP1 " f2 " #OP2 " s");                          \
                                                                               \
    /* Field access */                                                         \
    List_ACCESS(typeF1, f1, f1P);                                              \
    List_CONST_ACCESS(typeF2, f2, f2P);                                        \
                                                                               \
    /* Loop f1 OP1 s OP2 f2 */                                                 \
    const label loop_len = (f1).size();                                        \
                                                                               \
    /* pragmas... */                                                           \
    for (label i = 0; i < loop_len; ++i)                                       \
    {                                                                          \
        (f1P[i]) OP1 (f2P[i]) OP2 (s);                                         \
    }                                                                          \
}


// Member operator : this field f1 OP f2

#define TFOR_ALL_F_OP_F(typeF1, f1, OP, typeF2, f2)                            \
{                                                                              \
    /* Check fields have same size */                                          \
    checkFields(f1, f2, "f1 " #OP " f2");                                      \
                                                                               \
    /* Field access */                                                         \
    List_ACCESS(typeF1, f1, f1P);                                              \
    List_CONST_ACCESS(typeF2, f2, f2P);                                        \
                                                                               \
    /* Loop: f1 OP f2 */                                                       \
    const label loop_len = (f1).size();                                        \
                                                                               \
    /* pragmas... */                                                           \
    for (label i = 0; i < loop_len; ++i)                                       \
    {                                                                          \
        (f1P[i]) OP (f2P[i]);                                                  \
    }                                                                          \
}


// Member operator : this field f1 OP1 OP2 f2

#define TFOR_ALL_F_OP_OP_F(typeF1, f1, OP1, OP2, typeF2, f2)                   \
{                                                                              \
    /* Check fields have same size */                                          \
    checkFields(f1, f2, #OP1 " " #OP2 " f2");                                  \
                                                                               \
    /* Field access */                                                         \
    List_ACCESS(typeF1, f1, f1P);                                              \
    List_CONST_ACCESS(typeF2, f2, f2P);                                        \
                                                                               \
    /* Loop: f1 OP1 OP2 f2 */                                                  \
    const label loop_len = (f1).size();                                        \
                                                                               \
    /* pragmas... */                                                           \
    for (label i = 0; i < loop_len; ++i)                                       \
    {                                                                          \
        (f1P[i]) OP1 OP2 (f2P[i]);                                             \
    }                                                                          \
}


// Member operator : this field f OP s

#define TFOR_ALL_F_OP_S(typeF, f, OP, typeS, s)                                \
{                                                                              \
    /* Field access */                                                         \
    List_ACCESS(typeF, f, fP);                                                 \
                                                                               \
    /* Loop: f OP s */                                                         \
    const label loop_len = (f).size();                                         \
                                                                               \
    /* pragmas... */                                                           \
    for (label i = 0; i < loop_len; ++i)                                       \
    {                                                                          \
        (fP[i]) OP (s);                                                        \
    }                                                                          \
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Friend operator function : s OP f, allocates storage for s

#define TFOR_ALL_S_OP_F(typeS, s, OP, typeF, f)                                \
{                                                                              \
    /* Field access */                                                         \
    List_CONST_ACCESS(typeF, f, fP);                                           \
                                                                               \
    /* Loop: s OP f */                                                         \
    const label loop_len = (f).size();                                         \
                                                                               \
    /* pragmas, reduction... */                                                \
    for (label i = 0; i < loop_len; ++i)                                       \
    {                                                                          \
        (s) OP (fP[i]);                                                        \
    }                                                                          \
}


// Friend operator function : s OP1 f1 OP2 f2, allocates storage for s

#define TFOR_ALL_S_OP_F_OP_F(typeS, s, OP1, typeF1, f1, OP2, typeF2, f2)       \
{                                                                              \
    /* Field access */                                                         \
    List_CONST_ACCESS(typeF1, f1, f1P);                                        \
    List_CONST_ACCESS(typeF2, f2, f2P);                                        \
                                                                               \
    /* Loop: s OP1 f1 OP2 f2 */                                                \
    const label loop_len = (f1).size();                                        \
                                                                               \
    /* pragmas, reduction... */                                                \
    for (label i = 0; i < loop_len; ++i)                                       \
    {                                                                          \
        (s) OP1 (f1P[i]) OP2 (f2P[i]);                                         \
    }                                                                          \
}


// Friend operator function : s OP FUNC(f), allocates storage for s
#define TFOR_ALL_S_OP_FUNC_F(typeS, s, OP, FUNC, typeF, f)                     \
{                                                                              \
    /* Field access */                                                         \
    List_CONST_ACCESS(typeF, f, fP);                                           \
                                                                               \
    /* Loop: s OP FUNC(f) */                                                   \
    const label loop_len = (f).size();                                         \
                                                                               \
    /* pragmas, reduction... */                                                \
    for (label i = 0; i < loop_len; ++i)                                       \
    {                                                                          \
        (s) OP FUNC(fP[i]);                                                    \
    }                                                                          \
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
