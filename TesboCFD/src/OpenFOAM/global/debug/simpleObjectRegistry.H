/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
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

Class
    Foam::simpleObjectRegistry

Description
    Object registry for simpleRegIOobject. Maintains ordering.

SourceFiles
    simpleObjectRegistry.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_simpleObjectRegistry_H
#define Foam_simpleObjectRegistry_H

#include "Dictionary.H"
#include "simpleRegIOobject.H"
#include <string>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class dictionary;

/*---------------------------------------------------------------------------*\
                   Class simpleObjectRegistryEntry Declaration
\*---------------------------------------------------------------------------*/

class simpleObjectRegistryEntry
:
    public Dictionary<simpleObjectRegistryEntry>::link,
    public List<simpleRegIOobject*>
{
public:

    //- Construct with a single object (list size == 1)
    explicit simpleObjectRegistryEntry(simpleRegIOobject* obj)
    :
        List<simpleRegIOobject*>(Foam::one{}, obj)
    {}

    //- Construct with a List of objects
    explicit simpleObjectRegistryEntry(const List<simpleRegIOobject*>& data)
    :
        List<simpleRegIOobject*>(data)
    {}
};


/*---------------------------------------------------------------------------*\
                      Class simpleObjectRegistry Declaration
\*---------------------------------------------------------------------------*/

class simpleObjectRegistry
:
    public Dictionary<simpleObjectRegistryEntry>
{
public:

    // Constructors

        //- Default construct: empty without allocation (capacity=0).
        simpleObjectRegistry() = default;

        //- Construct empty with initial table capacity
        explicit simpleObjectRegistry(const label initialCapacity)
        :
            Dictionary<simpleObjectRegistryEntry>(initialCapacity)
        {}


    // Member Functions

        //- Set values (invoke callbacks) from dictionary entries
        //  Reporting honours the infoDetailLevel
        void setValues(const dictionary& dict, bool report=false);

        //- Set named value, but also handle embedded 'name=value' syntax
        //  Treats 'name=' and 'name' identically.
        //
        //  Will also accept a floating point value, but the called objects
        //  need to handle that.
        //  Reporting honours the infoDetailLevel
        void setNamedValue(std::string name, int val, bool report=false);

        //- Deprecated(2021-09) renamed to setNamedValue
        //  \deprecated(2021-09) renamed to setNamedValue
        void setNamedInt(std::string name, int val, bool report=false)
        {
            this->setNamedValue(name, val, report);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
