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

Description
    Macros for easy insertion into run-time selection tables

Note
    The helper macro names used here must remain synchronized with
    definitions in runTimeSelectionTables.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_addToRunTimeSelectionTable_H
#define Foam_addToRunTimeSelectionTable_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Add to construction table with typeName as the key
#define addToRunTimeSelectionTable\
(baseType,thisType,argNames)                                                   \
                                                                               \
    /* Add thisType factory method to the table */                             \
    baseType::add##argNames##ConstructorToTable<thisType>                      \
        add##thisType##argNames##ConstructorTo##baseType##Table_


//- Add to construction table with typeName as the key
#define addRemovableToRunTimeSelectionTable\
(baseType,thisType,argNames)                                                   \
                                                                               \
    /* Add thisType factory method to the table */                             \
    baseType::addRemovable##argNames##ConstructorToTable<thisType>             \
        addRemovable##thisType##argNames##ConstructorTo##baseType##Table_


//- Add to construction table with 'lookupName' as the key
#define addNamedToRunTimeSelectionTable\
(baseType,thisType,argNames,lookupName)                                        \
                                                                               \
    /* Add thisType factory method to the table, find by lookupName */         \
    baseType::add##argNames##ConstructorToTable<thisType>                      \
        add##thisType##argNames##ConstructorTo## \
        baseType##Table_##lookupName##_(#lookupName)


//- Add to construction table with 'lookupName' as the key
#define addRemovableNamedToRunTimeSelectionTable\
(baseType,thisType,argNames,lookupName)                                        \
                                                                               \
    /* Add thisType factory method to the table, find by lookupName */         \
    baseType::addRemovable##argNames##ConstructorToTable<thisType>             \
        addRemovable##thisType##argNames##ConstructorTo## \
        baseType##Table_##lookupName##_(#lookupName)


//- Add lookup alias for runTime selection
#define addAliasToRunTimeSelectionTable\
(baseType,thisType,argNames,lookup,other,ver)                                  \
                                                                               \
    /* Add thisType constructor function to the table, find by lookup */       \
    baseType::addAlias##argNames##ConstructorToTable<thisType>                 \
        add##thisType##argNames##ConstructorTo##baseType##Table_               \
        ##lookup##_##other##_(#lookup,#other,ver)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Add to construction table with typeName as the key.
//  Use when baseType doesn't need a template argument (eg, is a typedef)
#define addTemplateToRunTimeSelectionTable\
(baseType,thisType,Targ,argNames)                                              \
                                                                               \
    /* Add thisType factory method to the table */                             \
    baseType::add##argNames##ConstructorToTable<thisType<Targ>>                \
        add##thisType##Targ##argNames##ConstructorTo##baseType##Table_


//- Add to construction table with 'lookupName' as the key.
//  Use when baseType doesn't need a template argument (eg, is a typedef)
#define addNamedTemplateToRunTimeSelectionTable\
(baseType,thisType,Targ,argNames,lookupName)                                   \
                                                                               \
    /* Add thisType factory method to the table, find by lookupName */         \
    baseType::add##argNames##ConstructorToTable<thisType<Targ>>                \
        add##thisType##Targ##argNames##ConstructorTo## \
        baseType##Table_##lookupName##_(#lookupName)


//- Add lookup alias for for runTime selection
#define addAliasTemplateToRunTimeSelectionTable\
(baseType,thisType,Targ,argNames,lookup,other,ver)                             \
                                                                               \
    /* Add thisType constructor function to the table, find by lookup */       \
    baseType::addAlias##argNames##ConstructorToTable<thisType<Targ>>           \
        add##thisType##Targs##argNames##ConstructorTo##baseType##Table_        \
        ##lookup##_##other##_(#lookup,#other,ver)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Add to construction table with typeName as the key.
//  Use when baseType requires the Targ template argument as well
#define addTemplatedToRunTimeSelectionTable\
(baseType,thisType,Targ,argNames)                                              \
                                                                               \
    /* Add thisType factory method to the table */                             \
    baseType<Targ>::add##argNames##ConstructorToTable<thisType<Targ>>          \
        add##thisType##Targ##argNames##ConstructorTo##baseType##Targ##Table_


//- Add to construction table with 'lookupName' as the key.
//  Use when baseType requires the Targ template argument as well
#define addNamedTemplatedToRunTimeSelectionTable\
(baseType,thisType,Targ,argNames,lookupName)                                   \
                                                                               \
    /* Add thisType factory method to the table, find by lookupName */         \
    baseType<Targ>::add##argNames##ConstructorToTable<thisType<Targ>>          \
        add##thisType##Targ##argNames##ConstructorTo## \
        baseType##Targ##Table_##lookupName##_(#lookupName)


//- Add lookup alias for for runTime selection
//  Use when baseType requires the Targ template argument as well
#define addAliasTemplatedToRunTimeSelectionTable\
(baseType,thisType,Targ,argNames,lookup,other,ver)                             \
                                                                               \
    /* Add the thisType constructor function to the table, find by lookup */   \
    baseType<Targ>::add##argNames##ConstructorToTable<thisType<Targ>>          \
        add##thisType##Targ##argNames##ConstructorTo##baseType##Targ## \
        Table_##lookup##_##other##_(#lookup,#other,ver)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
