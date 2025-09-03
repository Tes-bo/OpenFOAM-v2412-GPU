/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2023 OpenCFD Ltd.
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

#include "HashTable.H"
#include "Istream.H"
#include "Ostream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
Foam::HashTable<T, Key, Hash>::HashTable(Istream& is)
:
    HashTable<T, Key, Hash>()
{
    operator>>(is, *this);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
Foam::Ostream& Foam::HashTable<T, Key, Hash>::printInfo(Ostream& os) const
{
    label used = 0;
    label maxChain = 0;
    unsigned avgChain = 0;

    for (label i=0; i < capacity_; ++i)
    {
        label count = 0;
        for (node_type* ep = table_[i]; ep; ep = ep->next_)
        {
            ++count;
        }

        if (count)
        {
            ++used;
            avgChain += count;

            if (maxChain < count)
            {
                maxChain = count;
            }
        }
    }

    os  << "HashTable<T,Key,Hash>"
        << " elements:" << size() << " slots:" << used << "/" << capacity_
        << " chaining(avg/max):" << (used ? (float(avgChain)/used) : 0)
        << "/" << maxChain << endl;

    return os;
}


template<class T, class Key, class Hash>
Foam::Ostream& Foam::HashTable<T, Key, Hash>::writeKeys
(
    Ostream& os,
    const label shortLen
) const
{
    // Similar to UList::writeList version except the following:
    // - the keys can never be uniform
    // - never write in binary

    label i = this->size();

    if
    (
        (i <= 1 || !shortLen)
     || (i <= shortLen)
    )
    {
        // Write size and start delimiter
        os << i << token::BEGIN_LIST;

        i = 0;
        for (const_iterator iter = this->cbegin(); iter != this->cend(); ++iter)
        {
            if (i++) os << token::SPACE;
            os << iter.key();
        }

        os << token::END_LIST;  // End delimiter
    }
    else
    {
        // Write size and start delimiter
        os << nl << i << nl << token::BEGIN_LIST << nl;

        for (const_iterator iter = this->cbegin(); iter != this->cend(); ++iter)
        {
            os << iter.key() << nl;
        }

        os << token::END_LIST << nl;  // End delimiter
    }

    os.check(FUNCTION_NAME);
    return os;
}


template<class T, class Key, class Hash>
Foam::Istream& Foam::HashTable<T, Key, Hash>::readTable
(
    Istream& is
)
{
    HashTable<T, Key, Hash>& tbl = *this;

    // Anull existing table
    tbl.clear();

    is.fatalCheck(FUNCTION_NAME);

    token tok(is);

    is.fatalCheck
    (
        "operator>>(Istream&, HashTable&) : "
        "reading first token"
    );

    if (tok.isLabel())
    {
        const label len = tok.labelToken();

        // Read beginning of contents
        const char delimiter = is.readBeginList("HashTable");

        if (len)
        {
            if (delimiter != token::BEGIN_LIST)
            {
                FatalIOErrorInFunction(is)
                    << "incorrect first token, '(', found "
                    << tok.info() << nl
                    << exit(FatalIOError);
            }

            tbl.reserve(tbl.size() + len);

            for (label i=0; i<len; ++i)
            {
                Key key;

                is >> key;          // Read the key
                T& val = tbl(key);  // Insert nameless T() into table
                is >> val;          // Read directly into the table value

                is.fatalCheck
                (
                    "operator>>(Istream&, HashTable&) : "
                    "reading entry"
                );
            }
        }

        // Read end of contents
        is.readEndList("HashTable");
    }
    else if (tok.isPunctuation(token::BEGIN_LIST))
    {
        is >> tok;
        while (!tok.isPunctuation(token::END_LIST))
        {
            is.putBack(tok);

            Key key;

            is >> key;          // Read the key
            T& val = tbl(key);  // Insert nameless T() into table
            is >> val;          // Read directly into the table value

            is.fatalCheck
            (
                "operator>>(Istream&, HashTable&) : "
                "reading entry"
            );

            is >> tok;
        }
    }
    else
    {
        FatalIOErrorInFunction(is)
            << "incorrect first token, expected <int> or '(', found "
            << tok.info() << nl
            << exit(FatalIOError);
    }

    is.fatalCheck(FUNCTION_NAME);
    return is;
}


template<class T, class Key, class Hash>
Foam::Ostream& Foam::HashTable<T, Key, Hash>::writeTable
(
    Ostream& os
) const
{
    const HashTable<T, Key, Hash>& tbl = *this;

    const label len = tbl.size();

    if (len)
    {
        // Size and start list delimiter
        os << nl << len << nl << token::BEGIN_LIST << nl;

        // Contents
        for (auto iter = tbl.cbegin(); iter != tbl.cend(); ++iter)
        {
            iter.print(os) << nl;
        }

        os << token::END_LIST;    // End list delimiter
    }
    else
    {
        // Empty hash table
        os << len << token::BEGIN_LIST << token::END_LIST;
    }

    os.check(FUNCTION_NAME);
    return os;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class T, class Key, class Hash>
Foam::Istream& Foam::operator>>
(
    Istream& is,
    HashTable<T, Key, Hash>& tbl
)
{
    return tbl.readTable(is);
}


template<class T, class Key, class Hash>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const HashTable<T, Key, Hash>& tbl
)
{
    return tbl.writeTable(os);
}


// ************************************************************************* //
