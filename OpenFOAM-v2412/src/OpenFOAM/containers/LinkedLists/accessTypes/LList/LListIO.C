/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2022 OpenCFD Ltd.
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

#include "LList.H"
#include "Istream.H"
#include "Ostream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class LListBase, class T>
Foam::LList<LListBase, T>::LList(Istream& is)
{
    operator>>(is, *this);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class LListBase, class T>
Foam::Istream& Foam::LList<LListBase, T>::readList(Istream& is)
{
    LList<LListBase, T>& list = *this;

    // Anull list
    list.clear();

    is.fatalCheck(FUNCTION_NAME);

    token tok(is);

    is.fatalCheck("LList::readList : reading first token");

    if (tok.isLabel())
    {
        const label len = tok.labelToken();

        // Begin of contents marker
        const char delimiter = is.readBeginList("LList");

        if (len)
        {
            if (delimiter == token::BEGIN_LIST)
            {
                for (label i=0; i<len; ++i)
                {
                    T elem;
                    is >> elem;
                    list.push_back(std::move(elem));
                }
            }
            else
            {
                // Uniform content (delimiter == token::BEGIN_BLOCK)

                T elem;
                is >> elem;

                for (label i=0; i<len; ++i)
                {
                    list.push_back(elem);
                }
            }
        }

        // End of contents marker
        is.readEndList("LList");
    }
    else if (tok.isPunctuation(token::BEGIN_LIST))
    {
        is >> tok;
        is.fatalCheck(FUNCTION_NAME);

        while (!tok.isPunctuation(token::END_LIST))
        {
            is.putBack(tok);

            T elem;
            is >> elem;
            list.push_back(std::move(elem));

            is >> tok;
            is.fatalCheck(FUNCTION_NAME);
        }
    }
    else
    {
        FatalIOErrorInFunction(is)
            << "incorrect first token, expected <int> or '(', found "
            << tok.info()
            << exit(FatalIOError);
    }

    is.fatalCheck(FUNCTION_NAME);
    return is;
}


template<class LListBase, class T>
Foam::Ostream& Foam::LList<LListBase, T>::writeList
(
    Ostream& os,
    const label shortLen
) const
{
    // NB: no binary, contiguous output

    const label len = this->size();

    if
    (
        (len <= 1 || !shortLen)
     || (len <= shortLen)
    )
    {
        // Size and start delimiter
        os << len << token::BEGIN_LIST;

        // Contents
        bool space = false;
        for (const T& val : *this)
        {
            if (space) os << token::SPACE;
            os << val;
            space = true;
        }

        // End delimiter
        os << token::END_LIST;
    }
    else
    {
        // Size and start delimiter
        os << nl << len << nl << token::BEGIN_LIST << nl;

        // Contents
        for (const T& val : *this)
        {
            os << val << nl;
        }

        // End delimiter
        os << token::END_LIST;
    }

    os.check(FUNCTION_NAME);
    return os;
}


template<class LListBase, class T>
Foam::Istream& Foam::operator>>(Istream& is, LList<LListBase, T>& list)
{
    return list.readList(is);
}


template<class LListBase, class T>
Foam::Ostream& Foam::operator<<(Ostream& os, const LList<LListBase, T>& list)
{
    return list.writeList(os, -1);  // Always with line breaks
}


// ************************************************************************* //
