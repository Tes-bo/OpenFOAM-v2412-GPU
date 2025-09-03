/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2023 OpenCFD Ltd.
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

#include "UList.H"
#include "Ostream.H"
#include "token.H"
#include "contiguous.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class T>
void Foam::UList<T>::writeEntry(Ostream& os) const
{
    const word tag("List<" + word(pTraits<T>::typeName) + '>');
    if (token::compound::isCompound(tag))
    {
        os  << tag << token::SPACE;
    }

    if (size())
    {
        os << *this;
    }
    else if
    (
        os.format() == IOstreamOption::BINARY
     || std::is_same<char, typename std::remove_cv<T>::type>::value
    )
    {
        // Zero-sized binary - Write size only
        // NB: special treatment for char data (binary I/O only)
        os << label(0);
    }
    else
    {
        // Zero-sized ASCII - Write size and delimiters
        os << label(0) << token::BEGIN_LIST << token::END_LIST;
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class T>
void Foam::UList<T>::writeEntry(const word& keyword, Ostream& os) const
{
    if (keyword.size())
    {
        os.writeKeyword(keyword);
    }
    writeEntry(os);
    os.endEntry();
}


template<class T>
Foam::Ostream& Foam::UList<T>::writeList
(
    Ostream& os,
    const label shortLen
) const
{
    const UList<T>& list = *this;

    const label len = list.size();

    if (os.format() == IOstreamOption::BINARY && is_contiguous<T>::value)
    {
        // Binary and contiguous

        os << nl << len << nl;

        if (len)
        {
            // write(...) includes surrounding start/end delimiters
            os.write(list.cdata_bytes(), list.size_bytes());
        }
    }
    else if (std::is_same<char, typename std::remove_cv<T>::type>::value)
    {
        // Special treatment for char data (binary I/O only)

        const auto oldFmt = os.format(IOstreamOption::BINARY);
        os << nl << len << nl;

        if (len)
        {
            // write(...) includes surrounding start/end delimiters
            os.write(list.cdata_bytes(), list.size_bytes());
        }

        os.format(oldFmt);
    }
    else if (is_contiguous<T>::value && len > 1 && list.uniform())
    {
        // Two or more entries, and all entries have identical values.
        os  << len << token::BEGIN_BLOCK << list[0] << token::END_BLOCK;
    }
    else if
    (
        (len <= 1 || !shortLen)
     ||
        (
            (len <= shortLen)
         &&
            (
                is_contiguous<T>::value
             || Detail::ListPolicy::no_linebreak<T>::value
            )
        )
    )
    {
        // Single-line output

        // Size and start delimiter
        os << len << token::BEGIN_LIST;

        auto iter = list.cbegin();
        const auto last = list.cend();

        // Contents
        if (iter != last)
        {
            os << *iter;

            for (++iter; (iter != last); (void)++iter)
            {
                os << token::SPACE << *iter;
            }
        }

        // End delimiter
        os << token::END_LIST;
    }
    else
    {
        // Multi-line output

        // Size and start delimiter
        os << nl << len << nl << token::BEGIN_LIST;

        auto iter = list.cbegin();
        const auto last = list.cend();

        // Contents
        for (/*nil*/; (iter != last); (void)++iter)
        {
            os << nl << *iter;
        }

        // End delimiter
        os << nl << token::END_LIST << nl;
    }

    os.check(FUNCTION_NAME);
    return os;
}


template<class T>
Foam::Istream& Foam::UList<T>::readList(Istream& is)
{
    UList<T>& list = *this;

    // The target list length - must match with sizes read
    const label len = list.size();

    is.fatalCheck(FUNCTION_NAME);

    token tok(is);

    is.fatalCheck("UList<T>::readList(Istream&) : reading first token");

    if (tok.isCompound())
    {
        // Compound: simply transfer contents

        List<T> elems;
        elems.transfer
        (
            tok.transferCompoundToken<List<T>>(is)
        );

        const label inputLen = elems.size();

        // List lengths must match
        if (inputLen != len)
        {
            FatalIOErrorInFunction(is)
                << "incorrect length for UList. Read "
                << inputLen << " expected " << len
                << exit(FatalIOError);
        }

        std::move(elems.begin(), elems.end(), list.begin());
    }
    else if (tok.isLabel())
    {
        // Label: could be int(..), int{...} or just a plain '0'

        const label inputLen = tok.labelToken();

        // List lengths must match
        if (inputLen != len)
        {
            FatalIOErrorInFunction(is)
                << "incorrect length for UList. Read "
                << inputLen << " expected " << len
                << exit(FatalIOError);
        }

        if (is.format() == IOstreamOption::BINARY && is_contiguous<T>::value)
        {
            // Binary and contiguous

            if (len)
            {
                Detail::readContiguous<T>
                (
                    is,
                    list.data_bytes(),
                    list.size_bytes()
                );

                is.fatalCheck
                (
                    "UList<T>::readList(Istream&) : "
                    "reading binary block"
                );
            }
        }
        else if (std::is_same<char, typename std::remove_cv<T>::type>::value)
        {
            // Special treatment for char data (binary I/O only)
            const auto oldFmt = is.format(IOstreamOption::BINARY);

            if (len)
            {
                // read(...) includes surrounding start/end delimiters
                is.read(list.data_bytes(), list.size_bytes());

                is.fatalCheck
                (
                    "UList<char>::readList(Istream&) : "
                    "reading binary block"
                );
            }

            is.format(oldFmt);
        }
        else
        {
            // Begin of contents marker
            const char delimiter = is.readBeginList("List");

            if (len)
            {
                if (delimiter == token::BEGIN_LIST)
                {
                    for (label i=0; i<len; ++i)
                    {
                        is >> list[i];

                        is.fatalCheck
                        (
                            "UList<T>::readList(Istream&) : "
                            "reading entry"
                        );
                    }
                }
                else
                {
                    // Uniform content (delimiter == token::BEGIN_BLOCK)

                    T elem;
                    is >> elem;

                    is.fatalCheck
                    (
                        "UList<T>::readList(Istream&) : "
                        "reading the single entry"
                    );

                    // Fill with the value
                    this->fill_uniform(elem);
                }
            }

            // End of contents marker
            is.readEndList("List");
        }
    }
    else if (tok.isPunctuation(token::BEGIN_LIST))
    {
        // "(...)" : read into list, handling size-mismatch after

        is >> tok;
        is.fatalCheck(FUNCTION_NAME);

        label inputLen = 0;

        while (!tok.isPunctuation(token::END_LIST))
        {
            is.putBack(tok);
            if (inputLen < len)
            {
                is >> list[inputLen];
            }
            else
            {
                // Read and discard
                T dummy;
                is >> dummy;
            }
            ++inputLen;

            is.fatalCheck
            (
                "UList<T>::readList(Istream&) : "
                "reading entry"
            );

            is >> tok;
            is.fatalCheck(FUNCTION_NAME);
        }

        // List lengths must match
        if (inputLen != len)
        {
            FatalIOErrorInFunction(is)
                << "incorrect length for UList. Read "
                << inputLen << " expected " << len
                << exit(FatalIOError);
        }
    }
    else
    {
        FatalIOErrorInFunction(is)
            << "incorrect first token, expected <int> or '(', found "
            << tok.info() << nl
            << exit(FatalIOError);
    }

    return is;
}


// ************************************************************************* //
