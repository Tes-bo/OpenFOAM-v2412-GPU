/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2017-2021 OpenCFD Ltd.
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

#include "LPtrList.H"
#include "Istream.H"
#include "Ostream.H"
#include "INew.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class LListBase, class T>
template<class INew>
void Foam::LPtrList<LListBase, T>::readIstream(Istream& is, const INew& inew)
{
    is.fatalCheck(FUNCTION_NAME);

    token tok(is);

    is.fatalCheck
    (
        "LPtrList::readIstream : "
        "reading first token"
    );

    if (tok.isLabel())
    {
        const label len = tok.labelToken();

        // Read beginning of contents
        const char delimiter = is.readBeginList("LPtrList");

        if (len)
        {
            if (delimiter == token::BEGIN_LIST)
            {
                for (label i=0; i<len; ++i)
                {
                    T* p = inew(is).ptr();
                    this->push_back(p);

                    is.fatalCheck
                    (
                        "LPtrList::readIstream : "
                        "reading entry"
                    );
                }
            }
            else  // Assumed to be token::BEGIN_BLOCK
            {
                T* p = inew(is).ptr();
                this->push_back(p);

                is.fatalCheck
                (
                    "LPtrList::readIstream : "
                    "reading the single entry"
                );

                for (label i=1; i<len; ++i)
                {
                    this->push_back(p->clone().ptr());
                }
            }
        }

        // Read end of contents
        is.readEndList("LPtrList");
    }
    else if (tok.isPunctuation(token::BEGIN_LIST))
    {
        is >> tok;
        is.fatalCheck(FUNCTION_NAME);

        while (!tok.isPunctuation(token::END_LIST))
        {
            is.putBack(tok);
            this->push_back(inew(is).ptr());

            is >> tok;
            is.fatalCheck(FUNCTION_NAME);
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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class LListBase, class T>
template<class INew>
Foam::LPtrList<LListBase, T>::LPtrList(Istream& is, const INew& inew)
{
    this->readIstream(is, inew);
}


template<class LListBase, class T>
Foam::LPtrList<LListBase, T>::LPtrList(Istream& is)
{
    this->readIstream(is, INew<T>());
}


// * * * * * * * * * * * * * * * Istream Operator  * * * * * * * * * * * * * //

template<class LListBase, class T>
Foam::Istream& Foam::operator>>(Istream& is, LPtrList<LListBase, T>& list)
{
    list.clear();
    list.readIstream(is, INew<T>());

    return is;
}


// * * * * * * * * * * * * * * * Ostream Operators * * * * * * * * * * * * * //

template<class LListBase, class T>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const LPtrList<LListBase, T>& list
)
{
    // Size and start delimiter
    os << nl << list.size() << nl << token::BEGIN_LIST << nl;

    // Contents
    for (auto iter = list.cbegin(); iter != list.cend(); ++iter)
    {
        os << *iter << nl;
    }

    // End delimiter
    os << token::END_LIST;

    os.check(FUNCTION_NAME);
    return os;
}

// ************************************************************************* //
