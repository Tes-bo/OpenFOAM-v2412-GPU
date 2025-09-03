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

\*---------------------------------------------------------------------------*/

#include "formattingEntry.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Write tokens without keyword, suppress/ignore bad tokens.
// Mostly like primitiveEntry::write(os, true);

static void writeTokens(Ostream& os, const tokenList& toks)
{
    bool started = false;  // Separate from previous token?

    for (const token& tok : toks)
    {
        if (!tok.good())  // silently ignore bad tokens
        {
            continue;
        }
        if (started)
        {
            os << token::SPACE;
        }
        else
        {
            started = true;
        }

        // Token output via direct handling in Ostream(s),
        // or normal '<<' output operator
        if (!os.write(tok))
        {
            os << tok;
        }

        if (tok.isCharData())
        {
            // If content appears to be a C++ comment
            // - better make certain that it gets a newline
            //   or later parsing of it will be a problem

            const auto& s = tok.stringToken();
            if (s.starts_with("//") && !s.ends_with('\n'))
            {
                os << '\n';
                started = false;  // Does not need further space separator
            }
        }
    }
}

} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::formattingEntry::formattingEntry
(
    const keyType& key,
    const char* s,
    std::streamsize len
)
:
    primitiveEntry(key, token::undefinedToken)  // Construct with one token
{
    if (s)
    {
        tokenList::front() =
            token(token::tokenType::CHAR_DATA, std::string(s, len));
    }
}


Foam::formattingEntry::formattingEntry
(
    const keyType& key,
    const std::string& content
)
:
    primitiveEntry
    (
        key,
        token(token::tokenType::CHAR_DATA, std::move(content))
    )
{}


Foam::formattingEntry::formattingEntry
(
    const keyType& key,
    std::string&& content
)
:
    primitiveEntry
    (
        key,
        token(token::tokenType::CHAR_DATA, std::move(content))
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::formattingEntry::write(Ostream& os) const
{
    if (active_)
    {
        writeTokens(os, *this);
    }
}


// ************************************************************************* //
