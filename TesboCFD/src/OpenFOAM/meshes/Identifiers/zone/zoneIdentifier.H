/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2023 OpenCFD Ltd.
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
    Foam::zoneIdentifier

Description
    Identifies a mesh zone by name and index, with optional physical type
    and group information.

SeeAlso
    patchIdentifier

SourceFiles
    zoneIdentifier.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_zoneIdentifier_H
#define Foam_zoneIdentifier_H

#include "wordList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class dictionary;

/*---------------------------------------------------------------------------*\
                       Class zoneIdentifier Declaration
\*---------------------------------------------------------------------------*/

class zoneIdentifier
{
    // Private Data

        //- Zone name
        word name_;

        //- Zone index in mesh
        label index_;

        //- Zone type (optional)
        word physicalType_;

        //- Groups to which the zone belongs (optional)
        wordList inGroups_;

public:

    // Generated Methods

        //- Copy construct
        zoneIdentifier(const zoneIdentifier&) = default;

        //- Move construct
        zoneIdentifier(zoneIdentifier&&) = default;

        //- Copy assignment
        zoneIdentifier& operator=(const zoneIdentifier&) = default;

        //- Move assignment
        zoneIdentifier& operator=(zoneIdentifier&&) = default;

        //- Destructor
        virtual ~zoneIdentifier() = default;


    // Constructors

        //- Default construct: name="", index=0
        zoneIdentifier();

        //- Construct from mandatory components
        zoneIdentifier(const word& name, const label index);

        //- Construct from components
        zoneIdentifier
        (
            const word& name,
            const label index,
            const word& physicalType,
            const wordList& inGroups = wordList()
        );

        //- Construct from dictionary
        zoneIdentifier
        (
            const word& name,
            const dictionary& dict,
            const label index
        );

        //- Copy construct, resetting the index (if non-negative)
        zoneIdentifier
        (
            const zoneIdentifier& ident,
            const label newIndex
        );

        //- Move construct, resetting the index (if non-negative)
        zoneIdentifier
        (
            zoneIdentifier&& ident,
            const label newIndex
        );


    // Member Functions

        //- The zone name
        const word& name() const noexcept { return name_; }

        //- Modifiable zone name
        word& name() noexcept { return name_; }

        //- The index of this zone in the zone list
        label index() const noexcept { return index_; }

        //- Modifiable index of this zone in the zone list
        label& index() noexcept { return index_; }

        //- The (optional) type of the zone
        const word& physicalType() const noexcept { return physicalType_; }

        //- Modifiable (optional) type of the zone
        word& physicalType() noexcept { return physicalType_; }

        //- The (optional) groups that the zone belongs to
        const wordList& inGroups() const noexcept { return inGroups_; }

        //- Modifiable (optional) groups that the zone belongs to
        wordList& inGroups() noexcept { return inGroups_; }

        //- True if given name is in a group
        bool inGroup(const word& name) const
        {
            return (!name.empty() && inGroups_.contains(name));
        }

        //- Add (unique) group for the zone
        void addGroup(const word& name)
        {
            if (!name.empty() && !inGroups_.contains(name))
            {
                inGroups_.push_back(name);
            }
        }

        //- Remove group for the zone
        void removeGroup(const word& name);

        //- Write (physicalType, inGroups) dictionary entries
        //- (without surrounding braces)
        void write(Ostream& os) const;
};


// Global Operators

//- Write (physicalType, inGroups) dictionary entries
//- (without surrounding braces)
Ostream& operator<<(Ostream& os, const zoneIdentifier& p);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
