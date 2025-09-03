/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenCFD Ltd.
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
    Containers for holding STARCCM interface definitions

\*---------------------------------------------------------------------------*/

#ifndef ccmInterfaceDefinitions_H
#define ccmInterfaceDefinitions_H

#include "Map.H"
#include "Ostream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace ccm
{

class interfaceEntry;
class interfaceDefinitions;

Ostream& operator<<(Ostream& os, const interfaceEntry& entry);
Ostream& operator<<(Ostream& os, const interfaceDefinitions& defs);

/*---------------------------------------------------------------------------*\
                  Class Foam::ccm::interfaceEntry Declaration
\*---------------------------------------------------------------------------*/

//- A STARCCM interface definition is a pair of boundary ids
class interfaceEntry
{
public:
    // Public Data

        //- The internal interface id
        label id;

        //- The first boundary
        label bnd0;

        //- The second boundary
        label bnd1;


    // Constructors

        //- Construct null
        interfaceEntry()
        :
            id(-1),
            bnd0(-1),
            bnd1(-1)
        {}


        //- Construct empty interface definition
        interfaceEntry(const label index)
        :
            id(index),
            bnd0(-1),
            bnd1(-1)
        {}


        //- Construct from components
        interfaceEntry
        (
            const label index,
            const label boundary0,
            const label boundary1
        )
        :
            id(index),
            bnd0(boundary0),
            bnd1(boundary1)
        {}


    // Access

        //- Check for in-place interfaces
        static bool isInPlace(const std::string& configurationType)
        {
            return configurationType == "IN_PLACE";
        }


        //- True if the boundary id is in this interface
        bool inInterface(label bndId) const
        {
            return bndId == bnd0 || bndId == bnd1;
        }

        //- True if all internal ids are non-negative
        bool good() const noexcept
        {
            return (id >= 0 && bnd0 >= 0 && bnd1 >= 0 && bnd0 != bnd1);
        }

        //- Same as good()
        bool valid() const noexcept { return good(); }


        //- Canonical name for boundary 0
        word canonicalName0() const
        {
            return "Interface" + ::Foam::name(id) + "_0";
        }

        //- Canonical name for boundary 1
        word canonicalName1() const
        {
            return "Interface" + ::Foam::name(id) + "_1";
        }

        //- Canonical name for boundary
        word canonicalName(label bndId) const
        {
            if (bndId == bnd0)
            {
                return canonicalName0();
            }
            else if (bndId == bnd1)
            {
                return canonicalName1();
            }
            else
            {
                return word::null;
            }
        }


    // IOstream Operators

        friend Ostream& operator<<
        (
            Ostream& os,
            const interfaceEntry& entry
        )
        {
            os  << "(" << entry.bnd0 << " " << entry.bnd1 << ")";
            return os;
        }

};


/*---------------------------------------------------------------------------*\
                  Class ccm::interfaceDefinitions Declaration
\*---------------------------------------------------------------------------*/

//- A list of available interface definitions
class interfaceDefinitions
:
    public Map<interfaceEntry>
{

    inline Map<interfaceEntry>& map()
    {
        return *this;
    }

    inline const Map<interfaceEntry>& map() const
    {
        return *this;
    }


public:
    // Constructor

        //- Null construct
        interfaceDefinitions()
        {}

        //- Size
        label size() const
        {
            return map().size();
        }

        //- Size
        bool empty() const
        {
            return map().empty();
        }

        //- Clear
        void clear()
        {
            map().clear();
        }


        //- Add (valid) interface entry
        bool add(const interfaceEntry& entry)
        {
            return (entry.good() && map().set(entry.id, entry));
        }


        //- Scan available interface entries for one matching this boundary id
        bool isInterface(label bndId)
        {
            forAllConstIters(map(), iter)
            {
                if (iter.val().inInterface(bndId))
                {
                    return true;
                }
            }
            return false;
        }


        //- Scan interface entries for one matching this boundary id
        //  return the canonical name
        word interfaceName(label bndId)
        {
            word ifname;
            forAllConstIters(map(), iter)
            {
                ifname = iter.val().canonicalName(bndId);
                if (!ifname.empty())
                {
                    break;
                }
            }

            return ifname;
        }


    // IOstream Operators

        friend Ostream& operator<<
        (
            Ostream& os,
            const interfaceDefinitions& defs
        )
        {
            os  << defs.map() << nl;
            return os;
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace ccm
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
