/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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
    Macro definitions for debug switches.

\*---------------------------------------------------------------------------*/

#ifndef Foam_defineDebugSwitch_H
#define Foam_defineDebugSwitch_H

#include "simpleRegIOobject.H"
#include "debug.H"
#include "label.H"  // Also for defining Foam::word etc.

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Define the debug information, lookup as \a name
template<class Type>
class RegisterDebugSwitch
:
    public ::Foam::simpleRegIOobject
{
public:

    //- The unique RegisterDebugSwitch object
    static const RegisterDebugSwitch registerDebugSwitch;

    //- Construct with given registered name
    explicit RegisterDebugSwitch(const char* name)
    :
        ::Foam::simpleRegIOobject(::Foam::debug::addDebugObject, name)
    {}

    virtual ~RegisterDebugSwitch() = default;

    virtual void readData(::Foam::Istream& is)
    {
        is >> Type::debug;
    }

    virtual void writeData(::Foam::Ostream& os) const
    {
        os << Type::debug;
    }

    RegisterDebugSwitch(const RegisterDebugSwitch<Type>&) = delete;
    void operator=(const RegisterDebugSwitch<Type>&) = delete;
};

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Define the debug information, lookup as \a Name
#define registerTemplateDebugSwitchWithName(Type, Name)                        \
                                                                               \
    template<> const Foam::RegisterDebugSwitch<Type>                           \
    Foam::RegisterDebugSwitch<Type>::registerDebugSwitch(Name)


//- Define the debug information, lookup as \a Name
#define registerDebugSwitchWithName(Type, Tag, Name)                           \
                                                                               \
    class add##Tag##ToDebug                                                    \
    :                                                                          \
        public ::Foam::simpleRegIOobject                                       \
    {                                                                          \
    public:                                                                    \
                                                                               \
        /* Construct with given registered name */                             \
        explicit add##Tag##ToDebug(const char* name)                           \
        :                                                                      \
            ::Foam::simpleRegIOobject(::Foam::debug::addDebugObject, name)     \
        {}                                                                     \
                                                                               \
        virtual ~add##Tag##ToDebug() = default;                                \
                                                                               \
        virtual void readData(::Foam::Istream& is)                             \
        {                                                                      \
            is >> Type::debug;                                                 \
        }                                                                      \
                                                                               \
        virtual void writeData(::Foam::Ostream& os) const                      \
        {                                                                      \
            os << Type::debug;                                                 \
        }                                                                      \
                                                                               \
        add##Tag##ToDebug(const add##Tag##ToDebug&) = delete;                  \
        void operator=(const add##Tag##ToDebug&) = delete;                     \
    };                                                                         \
    add##Tag##ToDebug add##Tag##ToDebug_(Name)


//- Define the debug information, lookup as \a Name
#define defineDebugSwitchWithName(Type, Name, Value)                           \
    int Type::debug(::Foam::debug::debugSwitch(Name, Value))

//- Define the debug information
#define defineDebugSwitch(Type, Value)                                         \
    defineDebugSwitchWithName(Type, Type::typeName_(), Value);                 \
    registerDebugSwitchWithName(Type, Type, Type::typeName_())

//- Define the debug information for templates, lookup as \a Name
#define defineTemplateDebugSwitchWithName(Type, Name, Value)                   \
    template<> defineDebugSwitchWithName(Type, Name, Value);                   \
    registerTemplateDebugSwitchWithName(Type, Name)

//- Define the debug information for templates
//  Useful with typedefs
#define defineTemplateDebugSwitch(Type, Value)                                 \
    defineTemplateDebugSwitchWithName(Type, #Type, Value)

//- Define the debug information directly for templates
#define defineNamedTemplateDebugSwitch(Type, Value)                            \
    defineTemplateDebugSwitchWithName(Type, Type::typeName_(), Value)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Old names

// Same as defineTemplateDebugSwitchWithName
#define defineTemplate2DebugSwitchWithName(Type, Name, Value)                  \
    defineTemplateDebugSwitchWithName(Type, Name, Value)

// Same as defineTemplateDebugSwitch
#define defineTemplate2DebugSwitch(Type, Value)                                \
    defineTemplateDebugSwitch(Type, Value)

// Same as defineNamedTemplateDebugSwitch
#define defineNamedTemplate2DebugSwitch(Type, Value)                           \
    defineNamedTemplateDebugSwitch(Type, Value)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
