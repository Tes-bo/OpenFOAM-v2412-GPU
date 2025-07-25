/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
    Copyright (C) 2020-2023 OpenCFD Ltd.
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
    Foam::patchIdentifier

Description
    Identifies a patch by name and index, with optional physical type
    and group information.

SourceFiles
    patchIdentifier.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_patchIdentifier_H
#define Foam_patchIdentifier_H

#include "wordList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class dictionary;

/*---------------------------------------------------------------------------*\
                       Class patchIdentifier Declaration
\*---------------------------------------------------------------------------*/

class patchIdentifier
{
    // Private Data

        //- Patch name
        word name_;

        //- Patch index in boundary
        label index_;

        //- Patch physical type (optional)
        word physicalType_;

        //- Groups to which the patch belongs (optional)
        wordList inGroups_;

public:

    // Static Member Functions

        //- Default patch name: "patch" or "patchN"
        static word defaultName(const label n = -1)
        {
            return
            (
                n < 0
              ? word("patch", false)
              : word("patch" + std::to_string(n), false)
            );
        }


    // Generated Methods

        //- Copy construct
        patchIdentifier(const patchIdentifier&) = default;

        //- Move construct
        patchIdentifier(patchIdentifier&&) = default;

        //- Copy assignment
        patchIdentifier& operator=(const patchIdentifier&) = default;

        //- Move assignment
        patchIdentifier& operator=(patchIdentifier&&) = default;

        //- Destructor
        virtual ~patchIdentifier() = default;


    // Constructors

        //- Default construct: name="", index=0
        patchIdentifier();

        //- Construct from mandatory components
        patchIdentifier(const word& name, const label index);

        //- Construct from components
        patchIdentifier
        (
            const word& name,
            const label index,
            const word& physicalType,
            const wordList& inGroups = wordList()
        );

        //- Construct from dictionary
        patchIdentifier
        (
            const word& name,
            const dictionary& dict,
            const label index
        );

        //- Copy construct, resetting the index (if non-negative)
        patchIdentifier
        (
            const patchIdentifier& ident,
            const label newIndex
        );

        //- Move construct, resetting the index (if non-negative)
        patchIdentifier
        (
            patchIdentifier&& ident,
            const label newIndex
        );


    // Member Functions

        //- The patch name
        const word& name() const noexcept { return name_; }

        //- Modifiable patch name
        word& name() noexcept { return name_; }

        //- The index of this patch in the boundaryMesh
        label index() const noexcept { return index_; }

        //- Modifiable index of this patch in the boundaryMesh
        label& index() noexcept { return index_; }

        //- The (optional) physical type of the patch
        const word& physicalType() const noexcept { return physicalType_; }

        //- Modifiable (optional) physical type of the patch
        word& physicalType() noexcept { return physicalType_; }

        //- The (optional) groups that the patch belongs to
        const wordList& inGroups() const noexcept { return inGroups_; }

        //- Modifiable (optional) groups that the patch belongs to
        wordList& inGroups() noexcept { return inGroups_; }

        //- True if given name is in a group
        bool inGroup(const word& name) const
        {
            return (!name.empty() && inGroups_.contains(name));
        }

        //- Add (unique) group for the patch
        void addGroup(const word& name)
        {
            if (!name.empty() && !inGroups_.contains(name))
            {
                inGroups_.push_back(name);
            }
        }

        //- Remove group for the patch
        void removeGroup(const word& name);

        //- Write (physicalType, inGroups) dictionary entries
        //- (without surrounding braces)
        void write(Ostream& os) const;
};


// Global Operators

//- Write (physicalType, inGroups) dictionary entries
//- (without surrounding braces)
Ostream& operator<<(Ostream& os, const patchIdentifier& ident);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
