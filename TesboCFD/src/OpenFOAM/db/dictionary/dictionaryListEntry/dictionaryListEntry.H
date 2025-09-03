/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
    Copyright (C) 2017 OpenCFD Ltd.
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
    Foam::dictionaryListEntry

Description
    Read/write List of dictionaries.

    The List entries get stored in a
    dictionary which itself is stored in the parent dictionary with
    keyword 'entryDDD' where DDD is the position in the parent dictionary.
    The printing is again in List format - the keyword is only printed as
    comment. Can be used to e.g. manipulate polyMesh/boundary files.

SourceFiles
    dictionaryListEntryIO.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_dictionaryListEntry_H
#define Foam_dictionaryListEntry_H

#include "dictionaryEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations

class dictionaryListEntry;
Ostream& operator<<(Ostream& os, const dictionaryListEntry& e);

template<>
Ostream& operator<<(Ostream&, const InfoProxy<dictionaryListEntry>&);


/*---------------------------------------------------------------------------*\
                     Class dictionaryListEntry Declaration
\*---------------------------------------------------------------------------*/

class dictionaryListEntry
:
    public dictionaryEntry
{
public:

    // Generated Methods

        //- No copy construct
        dictionaryListEntry(const dictionaryListEntry&) = delete;


    // Constructors

        //- Construct from the parent dictionary and Istream
        dictionaryListEntry(const dictionary& parentDict, Istream& is);

        //- Construct as copy for the given parent dictionary
        dictionaryListEntry
        (
            const dictionary& parentDict,
            const dictionaryListEntry& dictEnt
        );

        autoPtr<entry> clone(const dictionary& parentDict) const
        {
            return autoPtr<entry>(new dictionaryListEntry(parentDict, *this));
        }


    // Member Functions

        //- Write
        virtual void write(Ostream& os) const;

        //- Return info proxy,
        //- used to print token information to a stream
        InfoProxy<dictionaryListEntry> info() const noexcept { return *this; }


    // Ostream Operator

        friend Ostream& operator<<(Ostream& os, const dictionaryListEntry& e);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
