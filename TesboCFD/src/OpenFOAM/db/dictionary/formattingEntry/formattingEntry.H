/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 Sergey Lesnik
    Copyright (C) 2023-2024 OpenCFD Ltd.
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
    Foam::formattingEntry

Description
    A dictionary entry writing only contents by default, used to save
    formatting symbols for the final output

\*---------------------------------------------------------------------------*/

#ifndef Foam_formattingEntry_H
#define Foam_formattingEntry_H

#include "primitiveEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class formattingEntry Declaration
\*---------------------------------------------------------------------------*/

class formattingEntry
:
    public primitiveEntry
{
    // Private Data

        //- The output visibility
        bool active_ = true;

public:

    // Static Member Functions

        //- Generate an entry keyword: "__format-entry__NNN".
        //- The generated names are unlikely to collide with user dictionaries
        static keyType defaultName(label n)
        {
            return keyType
            (
                "__format-entry__" + std::to_string(n),
                keyType::LITERAL
            );
        }


    // Constructors

        //- Inherit all constructors from primitiveEntry
        using primitiveEntry::primitiveEntry;

        //- Construct with character data, using the provided keyword
        formattingEntry(const keyType& key, const char* s, std::streamsize len);

        //- Construct with character data, using the provided keyword
        formattingEntry(const keyType& key, const std::string& content);

        //- Construct with character data, using the provided keyword
        formattingEntry(const keyType& key, std::string&& content);

        //- Construct with character data, using a generated keyword
        formattingEntry(const label n, const char* s, std::streamsize len)
        :
            formattingEntry(defaultName(n), s, len)
        {}

        //- Construct with character data, using a generated keyword
        formattingEntry(const label n, const std::string& content)
        :
            formattingEntry(defaultName(n), content)
        {}

        //- Construct with character data, using a generated keyword
        formattingEntry(const label n, std::string&& content)
        :
            formattingEntry(defaultName(n), std::move(content))
        {}

        //- Construct with token data, using a generated keyword
        formattingEntry(const label n, token&& tok, bool visible=true)
        :
            primitiveEntry(defaultName(n), std::move(tok)),
            active_(visible)
        {}

        //- Clone the entry
        virtual autoPtr<entry> clone(const dictionary&) const
        {
            return autoPtr<entry>(new formattingEntry(*this));
        }


    // Factory Methods


    // Member Functions

        //- Set output visibility on/off.
        //  \return the previous value
        bool active(bool on) noexcept
        {
            bool old(active_);
            active_ = on;
            return old;
        }

        //- Get the output visibility
        bool active() const noexcept
        {
            return active_;
        }

        //- Write content without the keyword.
        //  Special properties:
        //  - ignores any bad tokens on output.
        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
