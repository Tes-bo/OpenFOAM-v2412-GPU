/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2022 OpenCFD Ltd.
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

#include "PackedList.H"
#include "IOstreams.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<unsigned Width>
void Foam::PackedList<Width>::writeEntry(Ostream& os) const
{
    os  << *this;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<unsigned Width>
Foam::Ostream& Foam::PackedList<Width>::printBits
(
    Ostream& os,
    bool debugOutput
) const
{
    os << token::BEGIN_LIST << nl;

    const label nblocks = debugOutput ? blocks_.size() : num_blocks(size());
    for (label blocki = 0; blocki < nblocks; ++blocki)
    {
        BitOps::print(os, blocks_[blocki], '.') << nl;
    }

    os << token::END_LIST << nl;

    return os;
}


template<unsigned Width>
Foam::Istream& Foam::PackedList<Width>::readList(Istream& is)
{
    PackedList<Width>& list = *this;

    // Anull list
    list.clear();

    is.fatalCheck(FUNCTION_NAME);

    token tok(is);

    is.fatalCheck("PackedList::readList(Istream&) : reading first token");

    if (tok.isLabel())
    {
        const label len = tok.labelToken();

        // Set list length to that read
        list.resize(len);

        if (is.format() == IOstreamOption::BINARY)
        {
            // Binary (always contiguous)

            if (len)
            {
                // NOTE: independent of WM_LABEL_SIZE
                is.read(list.data_bytes(), list.size_bytes());

                is.fatalCheck
                (
                    "PackedList::readList(Istream&) : "
                    "reading the binary block"
                );
            }
        }
        else
        {
            // Begin of contents marker
            const char delimiter = is.readBeginList("PackedList");

            if (len)
            {
                if (delimiter == token::BEGIN_LIST)
                {
                    for (label i=0; i<len; ++i)
                    {
                        list.set(i, list.readValue(is));

                        is.fatalCheck
                        (
                            "PackedList::readList(Istream&) : "
                            "reading entry"
                        );
                    }
                }
                else  // token::BEGIN_BLOCK
                {
                    // Assign for all entries
                    list = list.readValue(is);

                    is.fatalCheck
                    (
                        "PackedList::readList(Istream&) : "
                        "reading the single entry"
                    );
                }
            }

            // End of contents marker
            is.readEndList("PackedList");
        }
    }
    else if (tok.isPunctuation(token::BEGIN_LIST))
    {
        is >> tok;
        is.fatalCheck(FUNCTION_NAME);

        while (!tok.isPunctuation(token::END_LIST))
        {
            is.putBack(tok);
            list.push_back(list.readValue(is));

            is >> tok;
            is.fatalCheck(FUNCTION_NAME);
        }
    }
    else if (tok.isPunctuation(token::BEGIN_BLOCK))
    {
        is >> tok;
        is.fatalCheck(FUNCTION_NAME);

        while (!tok.isPunctuation(token::END_BLOCK))
        {
            is.putBack(tok);
            list.setPair(is);

            is >> tok;
            is.fatalCheck(FUNCTION_NAME);
        }
    }
    else
    {
        FatalIOErrorInFunction(is)
            << "incorrect first token, expected <int>, '(' or '{', found "
            << tok.info() << nl
            << exit(FatalIOError);
    }

    return is;
}


template<unsigned Width>
Foam::Ostream& Foam::PackedList<Width>::writeList
(
    Ostream& os,
    const label shortLen
) const
{
    const PackedList<Width>& list = *this;
    const label len = list.size();

    if (os.format() == IOstreamOption::BINARY)
    {
        // Binary (always contiguous)

        os << nl << len << nl;

        if (len)
        {
            // write(...) includes surrounding start/end delimiters
            os.write(list.cdata_bytes(), list.size_bytes());
        }
    }
    else if (len > 1 && list.uniform())
    {
        // Two or more entries, and all entries have identical values.
        os << len << token::BEGIN_BLOCK << list[0] << token::END_BLOCK;
    }
    else if (!shortLen || len <= shortLen)
    {
        // Single-line output

        // Size and start delimiter
        os << len << token::BEGIN_LIST;

        // Contents
        for (label i=0; i < len; ++i)
        {
            if (i) os << token::SPACE;
            os << label(list.get(i));
        }

        // End delimiter
        os << token::END_LIST;
    }
    else
    {
        // Multi-line output

        // Size and start delimiter
        os << nl << len << nl << token::BEGIN_LIST << nl;

        // Contents
        for (label i=0; i < len; ++i)
        {
            os << label(list.get(i)) << nl;
        }

        // End delimiter
        os << token::END_LIST << nl;
    }

    return os;
}


template<unsigned Width>
void Foam::PackedList<Width>::writeEntry
(
    const word& keyword,
    Ostream& os
) const
{
    if (keyword.size())
    {
        os.writeKeyword(keyword);
    }
    writeEntry(os);
    os.endEntry();
}


// * * * * * * * * * * * * * *  Friend Operators * * * * * * * * * * * * * * //

template<unsigned Width>
Foam::Istream& Foam::operator>>(Istream& is, PackedList<Width>& list)
{
    return list.readList(is);
}


template<unsigned Width>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const InfoProxy<PackedList<Width>>& iproxy
)
{
    const auto& list = *iproxy;

    os  << "PackedList<" << Width
        << "> size=" << list.size() << "/" << list.capacity()
        << " (limits: max=" << PackedList<Width>::max_value
        << ", elem_per_block=" << PackedList<Width>::elem_per_block
        << ")"
        << nl;

    return os;
}


// ************************************************************************* //
