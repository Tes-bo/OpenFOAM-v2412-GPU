/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2016 OpenFOAM Foundation
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

Class
    Foam::RegisterSwitch

Description
    Class and registration macros for InfoSwitches and OptimisationSwitches
    to support reading from system/controlDict and dynamic update.

\*---------------------------------------------------------------------------*/

#ifndef registerSwitch_H
#define registerSwitch_H

#include "simpleRegIOobject.H"
#include "macros.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class RegisterSwitch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class RegisterSwitch
:
    public simpleRegIOobject
{
public:

    //- Reference to the switch variable that has been registered
    Type& value;

    RegisterSwitch
    (
        void (*registryFn)(const char* name, simpleRegIOobject*),
        const char* name,
        Type& switchVar
    )
    :
        simpleRegIOobject(registryFn, name),
        value(switchVar)
    {}

    virtual ~RegisterSwitch() = default;

    virtual void readData(Istream& is)
    {
        is >> value;
    }

    virtual void writeData(Ostream& os) const
    {
        os << value;
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define registerOptSwitch(Name, Type, SwitchVar)                               \
    static Foam::RegisterSwitch<Type> FILE_UNIQUE(_addToOpt_)                  \
    (Foam::debug::addOptimisationObject, Name, SwitchVar)


#define registerInfoSwitch(Name, Type, SwitchVar)                              \
    static Foam::RegisterSwitch<Type> FILE_UNIQUE(_addToOpt_)                  \
    (Foam::debug::addInfoObject, Name, SwitchVar)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
