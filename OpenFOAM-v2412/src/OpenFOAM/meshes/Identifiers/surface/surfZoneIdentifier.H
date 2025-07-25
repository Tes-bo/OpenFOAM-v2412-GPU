/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2016-2021 OpenCFD Ltd.
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
    Foam::surfZoneIdentifier

Description
    Identifies a surface patch/zone by name and index,
    with optional geometric type.

SeeAlso
    patchIdentifier

SourceFiles
    surfZoneIdentifier.C

\*---------------------------------------------------------------------------*/

#ifndef surfZoneIdentifier_H
#define surfZoneIdentifier_H

#include "word.H"
#include "label.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class dictionary;

/*---------------------------------------------------------------------------*\
                     Class surfZoneIdentifier Declaration
\*---------------------------------------------------------------------------*/

class surfZoneIdentifier
{
    // Private Data

        //- Patch/zone name
        word name_;

        //- Patch/zone index in meshed surface
        label index_;

        //- Patch/zone type (optional)
        word geometricType_;

public:

    // Public Data

        //- The name for an 'empty' type
        static constexpr const char* const emptyType = "empty";


    // Static Member Functions

        //- Default zone name: "zone" or "zoneN"
        static word defaultName(const label n = -1)
        {
            return
            (
                n < 0
              ? word("zone", false)
              : word("zone" + std::to_string(n), false)
            );
        }


    // Generated Methods

        //- Copy construct
        surfZoneIdentifier(const surfZoneIdentifier&) = default;

        //- Copy assignment
        surfZoneIdentifier& operator=(const surfZoneIdentifier&) = default;


    // Constructors

        //- Default construct. Uses name="", index=0
        surfZoneIdentifier();

        //- Construct null with specified index
        explicit surfZoneIdentifier(const label index);

        //- Construct from mandatory components
        surfZoneIdentifier(const word& name, const label index);

        //- Construct from components
        surfZoneIdentifier
        (
            const word& name,
            const label index,
            const word& geometricType
        );

        //- Construct from dictionary
        surfZoneIdentifier
        (
            const word& name,
            const dictionary& dict,
            const label index
        );

        //- Copy construct, resetting the index
        surfZoneIdentifier
        (
            const surfZoneIdentifier& ident,
            const label index
        );


    // Member Functions

        //- The patch/zone name
        const word& name() const
        {
            return name_;
        }

        //- Modifiable patch/zone name
        word& name()
        {
            return name_;
        }

        //- The geometric type of the patch/zone
        const word& geometricType() const
        {
            return geometricType_;
        }

        //- Modifiable geometric type of the patch/zone
        word& geometricType()
        {
            return geometricType_;
        }

        //- The index of this patch/zone in the surface mesh
        label index() const
        {
            return index_;
        }

        //- Modifiable index of this patch/zone in the surface mesh
        label& index()
        {
            return index_;
        }

        //- Write (geometricType) dictionary entry
        //- (without surrounding braces)
        void write(Ostream& os) const;
};


// Global Operators

//- Compare zone indentifiers for equality
bool operator==(const surfZoneIdentifier& a, const surfZoneIdentifier& b);

//- Compare zone indentifiers for inequality
bool operator!=(const surfZoneIdentifier& a, const surfZoneIdentifier& b);


//- Read name, geometricType
Istream& operator>>(Istream& is, surfZoneIdentifier& obj);

//- Write name, geometricType. Entries are quoted to support empty words.
Ostream& operator<<(Ostream& os, const surfZoneIdentifier& obj);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
