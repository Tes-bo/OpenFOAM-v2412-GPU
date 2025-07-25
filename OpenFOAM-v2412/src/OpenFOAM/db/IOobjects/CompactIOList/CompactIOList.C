/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2024 OpenCFD Ltd.
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

#include "CompactIOList.H"
#include "labelList.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class T, class BaseType>
bool Foam::CompactIOList<T, BaseType>::readIOcontents()
{
    if
    (
        readOpt() == IOobject::MUST_READ
     || (isReadOptional() && headerOk())
    )
    {
        Istream& is = readStream(word::null);

        if (headerClassName() == IOList<T>::typeName)
        {
            is >> static_cast<List<T>&>(*this);
            close();
        }
        else if (headerClassName() == typeName)
        {
            is >> *this;
            close();
        }
        else
        {
            FatalIOErrorInFunction(is)
                << "Unexpected class name " << headerClassName()
                << " expected " << typeName
                << " or " << IOList<T>::typeName << endl
                << "    while reading object " << name()
                << exit(FatalIOError);
        }

        return true;
    }

    return false;
}


template<class T, class BaseType>
bool Foam::CompactIOList<T, BaseType>::overflows() const
{
    const List<T>& lists = *this;

    label total = 0;
    for (const auto& sublist : lists)
    {
        const label prev = total;
        total += sublist.size();
        if (total < prev)
        {
            return true;
        }
    }
    return false;
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template<class T, class BaseType>
Foam::CompactIOList<T, BaseType>::CompactIOList(const IOobject& io)
:
    regIOobject(io)
{
    readIOcontents();
}


template<class T, class BaseType>
Foam::CompactIOList<T, BaseType>::CompactIOList
(
    const IOobject& io,
    Foam::zero
)
:
    regIOobject(io)
{
    readIOcontents();
}


template<class T, class BaseType>
Foam::CompactIOList<T, BaseType>::CompactIOList
(
    const IOobject& io,
    const label len
)
:
    regIOobject(io)
{
    if (!readIOcontents())
    {
        List<T>::resize(len);
    }
}


template<class T, class BaseType>
Foam::CompactIOList<T, BaseType>::CompactIOList
(
    const IOobject& io,
    const UList<T>& content
)
:
    regIOobject(io)
{
    if (!readIOcontents())
    {
        List<T>::operator=(content);
    }
}


template<class T, class BaseType>
Foam::CompactIOList<T, BaseType>::CompactIOList
(
    const IOobject& io,
    List<T>&& content
)
:
    regIOobject(io)
{
    List<T>::transfer(content);

    readIOcontents();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T, class BaseType>
bool Foam::CompactIOList<T, BaseType>::writeObject
(
    IOstreamOption streamOpt,
    const bool writeOnProc
) const
{
    if
    (
        streamOpt.format() == IOstreamOption::BINARY
     && overflows()
    )
    {
        streamOpt.format(IOstreamOption::ASCII);

        WarningInFunction
            << "Overall number of elements of CompactIOList of size "
            << this->size() << " overflows the representation of a label"
            << nl << "    Switching to ascii writing" << endl;
    }

    if (streamOpt.format() == IOstreamOption::ASCII)
    {
        // Change type to be non-compact format type
        const word oldTypeName(typeName);

        const_cast<word&>(typeName) = IOList<T>::typeName;

        bool good = regIOobject::writeObject(streamOpt, writeOnProc);

        // Change type back
        const_cast<word&>(typeName) = oldTypeName;

        return good;
    }

    return regIOobject::writeObject(streamOpt, writeOnProc);
}


template<class T, class BaseType>
bool Foam::CompactIOList<T, BaseType>::writeData(Ostream& os) const
{
    return (os << *this).good();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T, class BaseType>
void Foam::CompactIOList<T, BaseType>::operator=
(
    const CompactIOList<T, BaseType>& rhs
)
{
    List<T>::operator=(rhs);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class T, class BaseType>
Foam::Istream& Foam::operator>>
(
    Foam::Istream& is,
    Foam::CompactIOList<T, BaseType>& L
)
{
    // Read compact
    const labelList start(is);
    const List<BaseType> elems(is);

    // Convert
    L.setSize(start.size()-1);

    forAll(L, i)
    {
        T& subList = L[i];

        label index = start[i];
        subList.setSize(start[i+1] - index);

        forAll(subList, j)
        {
            subList[j] = elems[index++];
        }
    }

    return is;
}


template<class T, class BaseType>
Foam::Ostream& Foam::operator<<
(
    Foam::Ostream& os,
    const Foam::CompactIOList<T, BaseType>& L
)
{
    // Keep ASCII writing same
    if (os.format() == IOstreamOption::ASCII)
    {
        os << static_cast<const List<T>&>(L);
    }
    else
    {
        // Convert to compact format
        labelList start(L.size()+1);

        start[0] = 0;
        for (label i = 1; i < start.size(); i++)
        {
            const label prev = start[i-1];
            start[i] = prev+L[i-1].size();

            if (start[i] < prev)
            {
                FatalIOErrorInFunction(os)
                    << "Overall number of elements " << start[i]
                    << " of CompactIOList of size "
                    << L.size() << " overflows the representation of a label"
                    << endl << "Please recompile with a larger representation"
                    << " for label" << exit(FatalIOError);
            }
        }

        List<BaseType> elems(start[start.size()-1]);

        label elemI = 0;
        forAll(L, i)
        {
            const T& subList = L[i];

            forAll(subList, j)
            {
                elems[elemI++] = subList[j];
            }
        }
        os << start << elems;
    }

    return os;
}


// ************************************************************************* //
