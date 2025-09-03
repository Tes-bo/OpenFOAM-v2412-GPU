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

#include "FixedList.H"
#include "Istream.H"
#include "Ostream.H"
#include "token.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class T, unsigned N>
void Foam::FixedList<T, N>::writeEntry(Ostream& os) const
{
    const word tag("List<" + word(pTraits<T>::typeName) + '>');
    if (token::compound::isCompound(tag))
    {
        os  << tag << token::SPACE;
        if (os.format() == IOstreamOption::BINARY && is_contiguous<T>::value)
        {
            // Need the size too so that List<Type>::readList parses correctly
            os << static_cast<label>(N);
        }
    }
    os << *this;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T, unsigned N>
Foam::FixedList<T, N>::FixedList(Istream& is)
{
    this->readList(is);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class T, unsigned N>
void Foam::FixedList<T, N>::writeEntry
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


template<class T, unsigned N>
Foam::Ostream& Foam::FixedList<T, N>::writeList
(
    Ostream& os,
    const label shortLen
) const
{
    const FixedList<T, N>& list = *this;

    // Unlike UList, no compact ascii output since a FixedList is generally
    // small and we prefer a consistent appearance.
    // Eg, FixedList<T,2> or Pair<T> as "(-1 -1)", never as "2{-1}"

    if (os.format() == IOstreamOption::BINARY && is_contiguous<T>::value)
    {
        // Binary and contiguous. Size is always non-zero

        // write(...) includes surrounding start/end delimiters
        os.write(list.cdata_bytes(), list.size_bytes());
    }
    else if
    (
        (N <= 1 || !shortLen)
     ||
        (
            (N <= unsigned(shortLen))
         &&
            (
                is_contiguous<T>::value
             || Detail::ListPolicy::no_linebreak<T>::value
            )
        )
    )
    {
        // Single-line output

        // Start delimiter
        os << token::BEGIN_LIST;

        // Contents
        for (unsigned i=0; i<N; ++i)
        {
            if (i) os << token::SPACE;
            os << list[i];
        }

        // End delimiter
        os << token::END_LIST;
    }
    else
    {
        // Multi-line output

        // Start delimiter
        os << nl << token::BEGIN_LIST << nl;

        // Contents
        for (unsigned i=0; i<N; ++i)
        {
            os << list[i] << nl;
        }

        // End delimiter
        os << token::END_LIST << nl;
    }

    os.check(FUNCTION_NAME);
    return os;
}


template<class T, unsigned N>
Foam::Istream& Foam::FixedList<T, N>::readList
(
    Istream& is
)
{
    FixedList<T, N>& list = *this;

    is.fatalCheck(FUNCTION_NAME);

    if (is.format() == IOstreamOption::BINARY && is_contiguous<T>::value)
    {
        // Binary and contiguous. Length is non-zero

        Detail::readContiguous<T>
        (
            is,
            list.data_bytes(),
            list.size_bytes()
        );

        is.fatalCheck
        (
            "FixedList<T, N>::readList(Istream&) : "
            "reading the binary block"
        );
        return is;
    }
    else
    {
        token tok(is);

        is.fatalCheck
        (
            "FixedList<T, N>::readList(Istream&) : "
            "reading first token"
        );

        if (tok.isCompound())
        {
            // Compound: transfer contents
            // - in practice probably never reach this branch
            list = tok.transferCompoundToken<List<T>>(is);
            return is;
        }
        else if (tok.isLabel())
        {
            // List lengths must match
            list.checkSize(tok.labelToken());
        }
        else if (!tok.isPunctuation())
        {
            FatalIOErrorInFunction(is)
                << "incorrect first token, expected <label> or '(' , found "
                << tok.info() << nl
                << exit(FatalIOError);
        }
        else
        {
            // Putback the opening bracket
            is.putBack(tok);
        }

        // Begin of contents marker
        const char delimiter = is.readBeginList("FixedList");

        if (delimiter == token::BEGIN_LIST)
        {
            for (unsigned i=0; i<N; ++i)
            {
                is >> list[i];

                is.fatalCheck
                (
                    "FixedList<T, N>::readList(Istream&) : "
                    "reading entry"
                );
            }
        }
        else
        {
            // Uniform content (delimiter == token::BEGIN_BLOCK)
            // - compatibility for v1812 and earlier (see issue #1160)

            T elem;
            is >> elem;

            is.fatalCheck
            (
                "FixedList<T, N>::readList(Istream&) : "
                "reading the single entry"
            );

            // Fill with the value
            this->fill(elem);
        }

        // End of contents marker
        is.readEndList("FixedList");
    }

    return is;
}


// ************************************************************************* //
