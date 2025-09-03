/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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

InClass
    Foam::memberFunctionSelectionTables

Description
    Macros for easy insertion into member function selection tables

Note
    The helper macro names used here must remain synchronized with
    definitions in memberFunctionSelectionTables.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_addToMemberFunctionSelectionTable_H
#define Foam_addToMemberFunctionSelectionTable_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Add to hash-table of functions with typeName as the key.
#define addToMemberFunctionSelectionTable\
(baseType,thisType,funcName,argNames)                                          \
                                                                               \
    /* Add thisType funcName to the table */                                   \
    baseType::add##funcName##argNames##MemberFunctionToTable<thisType>         \
        add##thisType##funcName##argNames##MemberFunctionTo##baseType##Table_


//- Add to hash-table of functions with 'lookupName' as the key.
#define addNamedToMemberFunctionSelectionTable\
(baseType,thisType,funcName,argNames,lookupName)                               \
                                                                               \
    /* Add thisType funcName to the table, find by lookup name */              \
    baseType::add##funcName##argNames##MemberFunctionToTable<thisType>         \
        add##thisType##funcName##argNames##MemberFunctionTo## \
        baseType##Table_##lookupName##_(#lookupName)


//- Add lookup alias for runTime member selection.
#define addAliasToMemberFunctionSelectionTable\
(baseType,thisType,funcName,argNames,lookup,other,ver)                         \
                                                                               \
    /* Add thisType funcName to the table, find by lookup name */              \
    baseType::addAlias##funcName##argNames##MemberFunctionToTable<thisType>    \
        add##thisType##funcName##argNames##MemberFunctionTo##baseType##Table_  \
        ##lookup##_##other##_(#lookup,#other,ver)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Add to hash-table of functions with typeName as the key.
//  Use when baseType doesn't need a template argument (eg, is a typedef)
#define addTemplateToMemberFunctionSelectionTable\
(baseType,thisType,Targ,funcName,argNames)                                     \
                                                                               \
    /* Add thisType funcName to the table */                                   \
    baseType::add##funcName##argNames##MemberFunctionToTable<thisType<Targ>>   \
        add##thisType##Targ##funcName##argNames##MemberFunctionTo## \
        baseType##Table_


//- Add to hash-table of functions with 'lookupName' as the key.
//  Use when baseType doesn't need a template argument (eg, is a typedef)
#define addNamedTemplateToMemberFunctionSelectionTable\
(baseType,thisType,Targ,funcName,argNames,lookupName)                          \
                                                                               \
    /* Add thisType funcName to the table, find by lookupName */               \
    baseType::add##funcName##argNames##MemberFunctionToTable<thisType<Targ>>   \
        add##thisType##Targ##funcName##argNames##MemberFunctionTo## \
        baseType##Table_##lookupName##_(#lookupName)


//- Add lookup alias for runTime member selection.
#define addAliasTemplateToMemberFunctionSelectionTable\
(baseType,thisType,Targ,funcName,argNames,lookup,other,ver)                    \
                                                                               \
    /* Add thisType funcName to the table, find by lookup name */              \
    baseType::addAlias##funcName##argNames##MemberFunctionToTable              \
        <thisType<Targ>>                                                       \
        add##thisType##Targ##funcName##argNames##MemberFunctionTo##baseType    \
        ##lookup##_##other##_(#lookup,#other,ver)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Add to hash-table of functions with typeName as the key.
//  Use when baseType requires the Targ template argument as well
#define addTemplatedToMemberFunctionSelectionTable\
(baseType,thisType,Targ,funcName,argNames)                                     \
                                                                               \
    /* Add thisType funcName to the table */                                   \
    baseType<Targ>::add##funcName##argNames##MemberFunctionToTable             \
        <thisType<Targ>>                                                       \
        add##thisType##Targ##funcName##argNames##MemberFunctionTo## \
        baseType##Targ##Table_

//- Add to hash-table of functions with 'lookupName' as the key.
//  Use when baseType requires the Targ template argument as well
#define addNamedTemplatedToMemberFunctionSelectionTable\
(baseType,thisType,Targ,funcName,argNames,lookupName)                          \
                                                                               \
    /* Add thisType funcName to the table, find by lookupName */               \
    baseType<Targ>::add##funcName##argNames##MemberFunctionToTable             \
        <thisType<Targ>>                                                       \
        add##thisType##Targ##funcName##argNames##MemberFunctionTo## \
        baseType##Targ##Table_##lookupName##_(#lookupName)

//- Add lookup alias for runTime member selection.
#define addAliasTemplatedToMemberFunctionSelectionTable\
(baseType,thisType,Targ,funcName,argNames,lookup,other,ver)                    \
                                                                               \
    /* Add thisType funcName to the table, find by lookup name */              \
    baseType<Targ>::addAlias##funcName##argNames##MemberFunctionToTable        \
        <thisType<Targ>>                                                       \
        add##thisType##Targ##funcName##argNames##MemberFunctionTo##baseType##  \
        Targ##Table_##lookup##_##other##_(#lookup,#other,ver)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
