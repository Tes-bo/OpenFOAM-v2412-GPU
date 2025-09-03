/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "HashPtrTable.H"
#include "Istream.H"
#include "Ostream.H"
#include "INew.H"
#include "dictionary.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class T, class Key, class Hash>
template<class INew>
void Foam::HashPtrTable<T, Key, Hash>::readIstream
(
    Istream& is,
    const INew& inew
)
{
    is.fatalCheck(FUNCTION_NAME);

    token tok(is);

    is.fatalCheck
    (
        "HashPtrTable::readIstream : "
        "reading first token"
    );

    if (tok.isLabel())
    {
        const label len = tok.labelToken();

        // Read beginning of contents
        const char delimiter = is.readBeginList("HashPtrTable");

        if (len)
        {
            this->reserve(this->size() + len);

            if (delimiter == token::BEGIN_LIST)
            {
                for (label i=0; i<len; ++i)
                {
                    Key key;
                    is >> key;
                    this->set(key, inew(key, is).ptr());

                    is.fatalCheck
                    (
                        "HashPtrTable::readIstream : "
                        "reading entry"
                    );
                }
            }
            else
            {
                FatalIOErrorInFunction(is)
                    << "incorrect first token, '(', found "
                    << tok.info() << nl
                    << exit(FatalIOError);
            }
        }

        // Read end of contents
        is.readEndList("HashPtrTable");
    }
    else if (tok.isPunctuation(token::BEGIN_LIST))
    {
        is >> tok;

        while (!tok.isPunctuation(token::END_LIST))
        {
            is.putBack(tok);
            Key key;
            is >> key;
            this->set(key, inew(key, is).ptr());

            is.fatalCheck
            (
                "HashPtrTable::readIstream : "
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
}


template<class T, class Key, class Hash>
template<class INew>
void Foam::HashPtrTable<T, Key, Hash>::read
(
    const dictionary& dict,
    const INew& inew
)
{
    this->reserve(this->size() + dict.size());
    for (const entry& e : dict)
    {
        this->set(e.keyword(), inew(e.dict()).ptr());
    }
}


template<class T, class Key, class Hash>
void Foam::HashPtrTable<T, Key, Hash>::write(Ostream& os) const
{
    for (const_iterator iter = this->cbegin(); iter != this->cend(); ++iter)
    {
        const T* ptr = iter.val();
        if (ptr)
        {
            ptr->write(os);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
template<class INew>
Foam::HashPtrTable<T, Key, Hash>::HashPtrTable(Istream& is, const INew& inew)
{
    this->readIstream(is, inew);
}


template<class T, class Key, class Hash>
Foam::HashPtrTable<T, Key, Hash>::HashPtrTable(Istream& is)
{
    this->readIstream(is, INew<T>());
}


template<class T, class Key, class Hash>
Foam::HashPtrTable<T, Key, Hash>::HashPtrTable(const dictionary& dict)
{
    this->read(dict, INew<T>());
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class T, class Key, class Hash>
Foam::Istream& Foam::operator>>(Istream& is, HashPtrTable<T, Key, Hash>& tbl)
{
    tbl.clear();
    tbl.readIstream(is, INew<T>());

    return is;
}


// ************************************************************************* //
