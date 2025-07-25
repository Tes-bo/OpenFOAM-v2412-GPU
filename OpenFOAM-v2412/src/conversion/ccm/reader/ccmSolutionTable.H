/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2024 OpenCFD Ltd.
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
    Containers for holding ccm solution and field listings.

\*---------------------------------------------------------------------------*/

#ifndef ccmSolutionTable_H
#define ccmSolutionTable_H

#include "SLList.H"
#include "Ostream.H"
#include "wordRes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace ccm
{

// Forward Declarations
class fieldEntry;
class fieldTable;
class solutionEntry;

Ostream& operator<<(Ostream& os, const fieldEntry& entry);
Ostream& operator<<(Ostream& os, const fieldTable& entry);
Ostream& operator<<(Ostream& os, const solutionEntry& entry);

/*---------------------------------------------------------------------------*\
                    Class Foam::ccm::namesList Declaration
\*---------------------------------------------------------------------------*/

//- A linked-list that is searchable by the 'name()' of the items
template<class Type>
class namesList
:
    public SLList<Type>
{
public:

    using const_iterator = typename SLList<Type>::const_iterator;
    using iterator = typename SLList<Type>::iterator;


    // Constructors

        //- Default construct
        namesList() = default;


    // Access

        //- True if a list element has a name that matches key
        bool found(const word& key) const
        {
            for (const Type& item : *this)
            {
                if (item.name() == key)
                {
                    return true;
                }
            }

            return false;
        }


        //- Find list element by name
        iterator find(const word& key)
        {
            const auto last = SLList<Type>::end();

            iterator iter = SLList<Type>::begin();

            for (; (iter != last); ++iter)
            {
                if ((*iter).name() == key)
                {
                    break;
                }
            }

            return iter;
        }


        //- Return a list of names in allow-list and not in deny-list
        List<word> findNames
        (
            const wordRes& allow,
            const wordRes& deny = wordRes()
        ) const
        {
            List<word> matched(SLList<Type>::size());

            label count = 0;

            for (const Type& item : *this)
            {
                const word& name = item.name();

                if (allow.match(name) && !deny.match(name))
                {
                    matched[count] = name;
                    ++count;
                }
            }

            matched.resize(count);
            return matched;
        }
};


/*---------------------------------------------------------------------------*\
                    Class Foam::ccm::fieldEntry Declaration
\*---------------------------------------------------------------------------*/

//- A ccm field entry with short name, name, maxId and type
//  shortName => ( fullName, maxId, type );
class fieldEntry
{
    // Private Data

        //- The field name (PROSTAR short name)
        word name_;

        //- The full field name
        string fullName_;

        //- The field units
        string units_;

        //- The max cell id for the field
        label maxCellId_;

        //- The max face id for the field
        label maxFaceId_;

public:

    // Constructors

        //- Construct from components with optional units
        fieldEntry
        (
            const word& shortName,
            const string& fullName,
            const char* units = nullptr
        )
        :
            name_(shortName),
            fullName_(fullName),
            units_(),
            maxCellId_(0),
            maxFaceId_(0)
        {
            if (units && *units)
            {
                units_ = units;
            }
        }


    // Access

        //- The field name (PROSTAR short name)
        const word& name() const noexcept { return name_; }

        //- The full field name
        const string& fullName() const noexcept { return fullName_; }

        //- The field units
        const string& units() const noexcept { return units_; }

        //- The max cell id for the field
        label maxCellId() const noexcept { return maxCellId_; }

        //- The max face id for the field
        label maxFaceId() const noexcept { return maxFaceId_; }


    // Edit

        //- Set the field units
        void units(const std::string& units)
        {
            if (!units.empty())
            {
                units_ = units;
            }
        }

        //- Set the max cell Id for the field
        void maxCellId(const int val)
        {
            if (maxCellId_ < val)
            {
                maxCellId_ = val;
            }
        }


        //- Set the max face Id for the field
        void maxFaceId(const int val)
        {
            if (maxFaceId_ < val)
            {
                maxFaceId_ = val;
            }
        }


    // IOstream Operators

        friend Ostream& operator<<
        (
            Ostream& os,
            const fieldEntry& entry
        )
        {
            os  << entry.name_ << " => " << entry.fullName_
                << " [" << entry.units_.c_str()
                << "] maxCell: " << entry.maxCellId_
                << " maxFace: " << entry.maxFaceId_;

            return os;
        }
};


/*---------------------------------------------------------------------------*\
                  Class Foam::ccm::solutionEntry Declaration
\*---------------------------------------------------------------------------*/

//- A ccm solution entry with name, iteration and time
//  stateName => ( iteration, time );
class solutionEntry
{
    // Private Data

        //- The solution name
        word name_;

        //- The solution iteration/timestep
        label iter_;

        //- The solution time (sec)
        scalar time_;

public:

    // Constructors

        //- Construct from components
        solutionEntry
        (
            const word& name,
            const label iteration,
            const scalar timeValue = 0
        )
        :
            name_(name),
            iter_(iteration),
            time_(timeValue)
        {}

    // Access

        //- The solution name
        const word& name() const noexcept { return name_; }

        //- The solution iteration/timestep
        label iteration() const noexcept { return iter_; }

        //- The solution time (sec)
        scalar timeValue() const noexcept { return time_; }


    // IOstream Operators

        friend Ostream& operator<<
        (
            Ostream& os,
            const solutionEntry& entry
        )
        {
            os  << entry.name_ << " =>"
                << " iter: " << entry.iter_
                << " time: " << entry.time_;

            return os;
        }
};


/*---------------------------------------------------------------------------*\
                  Class Foam::ccm::solutionTable Declaration
\*---------------------------------------------------------------------------*/

// Typedef: ccm::solutionTable
// A list of all the available solutions
typedef namesList<solutionEntry> solutionTable;


/*---------------------------------------------------------------------------*\
                    Class Foam::ccm::fieldTable Declaration
\*---------------------------------------------------------------------------*/

//- A list of the available fields
class fieldTable
:
    public namesList<fieldEntry>
{
public:

    // Constructor

        //- Default construct
        fieldTable() = default;


    // Access

        //- The maximum cell Id referenced in the list
        label maxCellId() const
        {
            label result = 0;
            for (const auto& item : *this)
            {
                const label val = item.maxCellId();

                if (result < val)
                {
                    result = val;
                }
            }

            return result;
        }


        //- The maximum face Id referenced in the list
        label maxFaceId() const
        {
            label result = 0;
            for (const auto& item : *this)
            {
                const label val = item.maxFaceId();

                if (result < val)
                {
                    result = val;
                }
            }

            return result;
        }


    // IOstream Operators

        friend Ostream& operator<<
        (
            Ostream& os,
            const fieldTable& tbl
        )
        {
            os  << static_cast<const namesList<fieldEntry>&>(tbl)
                << nl
                << "maxCell: "  << tbl.maxCellId()
                << " maxFace: " << tbl.maxFaceId();

            return os;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace ccm
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
