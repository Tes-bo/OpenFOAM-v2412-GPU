/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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
    Macros to ease declaration of member function selection tables.

Note
    Uses macros from runTimeSelectionTables.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_memberFunctionSelectionTables_H
#define Foam_memberFunctionSelectionTables_H

#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//
// Declaration Macros
//

//- Declare a run-time member-function selection (variables and adder classes)
#define declareMemberFunctionSelectionTable\
(returnType,baseType,funcName,argNames,argList,parListUnused)                  \
                                                                               \
    declareRunTimeSelectionTableBase(                                          \
        returnType, funcName##argNames##MemberFunction, argList);              \
                                                                               \
    /* Helper to add funcName to table */                                      \
    template<class baseType##Type>                                             \
    struct add##funcName##argNames##MemberFunctionToTable                      \
    {                                                                          \
        explicit add##funcName##argNames##MemberFunctionToTable                \
        (                                                                      \
            const ::Foam::word& k = baseType##Type::typeName                   \
        )                                                                      \
        {                                                                      \
            funcName##argNames##MemberFunctionTablePtr_construct(true);        \
            if                                                                 \
            (                                                                  \
                !funcName##argNames##MemberFunctionTablePtr_                   \
                   ->insert(k, baseType##Type::funcName)                       \
            )                                                                  \
            {                                                                  \
                std::cerr                                                      \
                    << "Duplicate entry " << k << " in member table "          \
                    << #baseType << std::endl;                                 \
                ::Foam::error::safePrintStack(std::cerr);                      \
            }                                                                  \
        }                                                                      \
                                                                               \
        ~add##funcName##argNames##MemberFunctionToTable()                      \
        {                                                                      \
            funcName##argNames##MemberFunctionTablePtr_construct(false);       \
        }                                                                      \
                                                                               \
        add##funcName##argNames##MemberFunctionToTable                         \
            (const add##funcName##argNames##MemberFunctionToTable&)            \
            = delete;                                                          \
        void operator=                                                         \
            (const add##funcName##argNames##MemberFunctionToTable&)            \
            = delete;                                                          \
    };                                                                         \
                                                                               \
    /* Helper to add funcName to table */                                      \
    template<class baseType##Type>                                             \
    struct addRemovable##funcName##argNames##MemberFunctionToTable             \
    {                                                                          \
        const ::Foam::word name;  /* Retain name for later removal */          \
                                                                               \
        explicit addRemovable##funcName##argNames##MemberFunctionToTable       \
        (                                                                      \
            const ::Foam::word& k = baseType##Type::typeName                   \
        )                                                                      \
        :                                                                      \
            name(k)                                                            \
        {                                                                      \
            funcName##argNames##MemberFunctionTablePtr_construct(true);        \
            funcName##argNames##MemberFunctionTablePtr_                        \
               ->set(k, baseType##Type::funcName);                             \
        }                                                                      \
                                                                               \
        ~addRemovable##funcName##argNames##MemberFunctionToTable()             \
        {                                                                      \
            if (funcName##argNames##MemberFunctionTablePtr_)                   \
            {                                                                  \
                funcName##argNames##MemberFunctionTablePtr_->erase(name);      \
            }                                                                  \
        }                                                                      \
                                                                               \
        addRemovable##funcName##argNames##MemberFunctionToTable                \
            (const addRemovable##funcName##argNames##MemberFunctionToTable&)   \
            = delete;                                                          \
        void operator=                                                         \
            (const addRemovable##funcName##argNames##MemberFunctionToTable&)   \
            = delete;                                                          \
    };                                                                         \
                                                                               \
    /* Helper to add alias for funcName to table */                            \
    template<class baseType##Type>                                             \
    struct addAlias##funcName##argNames##MemberFunctionToTable                 \
    {                                                                          \
        explicit addAlias##funcName##argNames##MemberFunctionToTable           \
        (                                                                      \
            const ::Foam::word& k,                                             \
            const ::Foam::word& alias,                                         \
            const int ver                                                      \
        )                                                                      \
        {                                                                      \
            funcName##argNames##MemberFunctionCompatTable()                    \
                .set(alias, std::pair<::Foam::word,int>(k,ver));               \
        }                                                                      \
    };


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//
// Definition Macros
//

//- Define run-time selection table
#define defineMemberFunctionSelectionTable(baseType,funcName,argNames)         \
                                                                               \
    defineRunTimeSelectionTableBase(                                           \
        baseType,baseType::funcName##argNames##MemberFunction,)


//- Define run-time selection table for template classes
//  use when baseType doesn't need a template argument (eg, is a typedef)
#define defineTemplateMemberFunctionSelectionTable(baseType,funcName,argNames) \
                                                                               \
    defineRunTimeSelectionTableBase(                                           \
        baseType,baseType::funcName##argNames##MemberFunction,template<>)


//- Define run-time selection table for template classes
//  use when baseType requires the Targ template argument
#define defineTemplatedMemberFunctionSelectionTable\
(baseType,funcName,argNames,Targ)                                              \
                                                                               \
    defineRunTimeSelectionTableBase(                                           \
        baseType,baseType<Targ>::funcName##argNames##MemberFunction,template<>)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
