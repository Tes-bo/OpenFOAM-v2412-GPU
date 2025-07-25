/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022-2023 OpenCFD Ltd.
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

#include "DynamicList.H"
#include "Istream.H"
#include "contiguous.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T>
Foam::CircularBuffer<T>::CircularBuffer(Istream& is)
{
    this->readList(is);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class T>
Foam::Ostream& Foam::CircularBuffer<T>::info(Ostream& os) const
{
    os  << "size=" << size() << '/' << capacity()
        << " begin=" << begin_
        << " end=" << end_
        // << " one=" << this->range_one() << this->array_one()
        // << " two=" << this->range_two() << this->array_two()
        << nl;

    return os;
}


template<class T>
Foam::Istream& Foam::CircularBuffer<T>::readList(Istream& is)
{
    // Delegate to DynamicList for reading
    DynamicList<T> elements(std::move(storage_));
    elements.readList(is);

    // Reset the list addressing range
    begin_ = 0;
    end_ = elements.size();

    const label minLen = (end_ + min_size());

    if (!elements.empty() && (elements.capacity() < minLen))
    {
        // Avoid full buffer (beg/end ambiguity)
        // Use setCapacity instead of resize to avoid additional doubling...
        elements.setCapacity(minLen);
    }

    // Use the entire storage
    elements.resize(elements.capacity());

    storage_ = std::move(elements);

    return is;
}


template<class T>
Foam::Ostream& Foam::CircularBuffer<T>::writeList
(
    Ostream& os,
    const label shortLen
) const
{
    const label len = this->size();
    const auto list1 = this->array_one();
    const auto list2 = this->array_two();

    #ifdef FULLDEBUG
    if (len != (list1.size() + list2.size()))
    {
        FatalErrorInFunction
            << "Size check failed"
            << abort(FatalError);
    }
    #endif

    if (os.format() == IOstreamOption::BINARY && is_contiguous<T>::value)
    {
        // Binary and contiguous

        os << nl << len << nl;

        if (len)
        {
            // The TOTAL number of bytes to be written.
            // - possibly add start delimiter
            // This is much like IndirectListBase output

            os.beginRawWrite(len*sizeof(T));

            if (!list1.empty())
            {
                os.writeRaw(list1.cdata_bytes(), list1.size_bytes());
            }
            if (!list2.empty())
            {
                os.writeRaw(list2.cdata_bytes(), list2.size_bytes());
            }

            // End delimiter and/or cleanup.
            os.endRawWrite();
        }
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

        // Contents
        label i = 0;
        for (const T& val : list1)
        {
            if (i++) os << token::SPACE;
            os << val;
        }
        for (const T& val : list2)
        {
            if (i++) os << token::SPACE;
            os << val;
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
        for (const T& val : list1)
        {
            os << val << nl;
        }
        for (const T& val : list2)
        {
            os << val << nl;
        }

        // End delimiter
        os << token::END_LIST << nl;
    }

    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
