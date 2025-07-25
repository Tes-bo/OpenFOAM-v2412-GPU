/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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

#include "dictionaryListEntry.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

// File-scope: The dictionary size without the "FoamFile" entry
static inline Foam::label realSize(const Foam::dictionary& dict)
{
    return
    (
        (dict.empty() || (dict.front()->keyword() != "FoamFile"))
      ? dict.size()
      : dict.size() - 1
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dictionaryListEntry::dictionaryListEntry
(
    const dictionary& parentDict,
    const dictionaryListEntry& dictEnt
)
:
    dictionaryEntry(parentDict, dictEnt)
{}


Foam::dictionaryListEntry::dictionaryListEntry
(
    const dictionary& parentDict,
    Istream& is
)
:
    dictionaryEntry
    (
        word("entry" + Foam::name(realSize(parentDict))),
        parentDict,
        dictionary::null
    )
{
    token tok(is);
    if (tok.isLabel())
    {
        const label len = tok.labelToken();

        is.readBeginList("List");

        for (label i=0; i<len; ++i)
        {
            entry::New(*this, is);
        }
        is.readEndList("List");
    }
    else if (tok.isPunctuation(token::BEGIN_LIST))
    {
        while (true)
        {
            is >> tok;
            if (tok.error())
            {
                FatalIOErrorInFunction(is)
                    << "parsing error " << tok.info() << nl
                    << exit(FatalIOError);
            }
            else if (tok.isPunctuation(token::END_LIST))
            {
                break;
            }
            is.putBack(tok);
            entry::New(*this, is);
        }
    }
    else
    {
        FatalIOErrorInFunction(is)
            << "incorrect first token, expected <int> or '(', found "
            << tok.info() << nl
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::dictionaryListEntry::write(Ostream& os) const
{
    os  << nl << indent << size()
        << token::SPACE << "// " << keyword() << nl
        << indent << token::BEGIN_LIST << incrIndent << nl;

    // Write contents
    dictionary::write(os, false);

    // Write end delimiter
    os << decrIndent << indent << token::END_LIST << nl;

    os.check(FUNCTION_NAME);
}


// * * * * * * * * * * * * * * Ostream operator  * * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const dictionaryListEntry& e
)
{
    e.write(os);
    return os;
}


template<>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const InfoProxy<dictionaryListEntry>& iproxy
)
{
    const auto& e = *iproxy;

    os  << "    dictionaryListEntry '" << e.keyword() << "'" << endl;

    return os;
}


// ************************************************************************* //
