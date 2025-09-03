/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018-2024 OpenCFD Ltd.
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

#include "CompactIOField.H"
#include "labelList.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class T, class BaseType>
bool Foam::CompactIOField<T, BaseType>::readIOcontents(bool readOnProc)
{
    if (isReadRequired() || (isReadOptional() && headerOk()))
    {
        // Do reading
        Istream& is = readStream(word::null, readOnProc);

        if (readOnProc)
        {
            if (headerClassName() == IOField<T>::typeName)
            {
                is >> static_cast<Field<T>&>(*this);
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
                    << " or " << IOField<T>::typeName << nl
                    << "    while reading object " << name()
                    << exit(FatalIOError);
            }
        }

        return true;
    }

    return false;
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template<class T, class BaseType>
Foam::CompactIOField<T, BaseType>::CompactIOField(const IOobject& io)
:
    regIOobject(io)
{
    readIOcontents();
}


template<class T, class BaseType>
Foam::CompactIOField<T, BaseType>::CompactIOField
(
    const IOobject& io,
    const bool readOnProc
)
:
    regIOobject(io)
{
    readIOcontents(readOnProc);
}


template<class T, class BaseType>
Foam::CompactIOField<T, BaseType>::CompactIOField
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
Foam::CompactIOField<T, BaseType>::CompactIOField
(
    const IOobject& io,
    const label len
)
:
    regIOobject(io)
{
    if (!readIOcontents())
    {
        Field<T>::resize(len);
    }
}


template<class T, class BaseType>
Foam::CompactIOField<T, BaseType>::CompactIOField
(
    const IOobject& io,
    const UList<T>& content
)
:
    regIOobject(io)
{
    if (!readIOcontents())
    {
        Field<T>::operator=(content);
    }
}


template<class T, class BaseType>
Foam::CompactIOField<T, BaseType>::CompactIOField
(
    const IOobject& io,
    Field<T>&& content
)
:
    regIOobject(io)
{
    Field<T>::transfer(content);

    readIOcontents();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T, class BaseType>
bool Foam::CompactIOField<T, BaseType>::writeObject
(
    IOstreamOption streamOpt,
    const bool writeOnProc
) const
{
    if (streamOpt.format() == IOstreamOption::ASCII)
    {
        // Change type to be non-compact format type
        const word oldTypeName(typeName);

        const_cast<word&>(typeName) = IOField<T>::typeName;

        bool good = regIOobject::writeObject(streamOpt, writeOnProc);

        // Restore type
        const_cast<word&>(typeName) = oldTypeName;

        return good;
    }

    return regIOobject::writeObject(streamOpt, writeOnProc);
}


template<class T, class BaseType>
bool Foam::CompactIOField<T, BaseType>::writeData(Ostream& os) const
{
    return (os << *this).good();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T, class BaseType>
void Foam::CompactIOField<T, BaseType>::operator=
(
    const CompactIOField<T, BaseType>& rhs
)
{
    if (this == &rhs)
    {
        return;  // Self-assigment is a no-op
    }

    Field<T>::operator=(rhs);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class T, class BaseType>
Foam::Istream& Foam::operator>>
(
    Foam::Istream& is,
    Foam::CompactIOField<T, BaseType>& L
)
{
    // Read compact
    const labelList start(is);
    const Field<BaseType> elems(is);

    // Convert
    L.setSize(start.size()-1);

    forAll(L, i)
    {
        T& subField = L[i];

        label index = start[i];
        subField.setSize(start[i+1] - index);

        forAll(subField, j)
        {
            subField[j] = elems[index++];
        }
    }

    return is;
}


template<class T, class BaseType>
Foam::Ostream& Foam::operator<<
(
    Foam::Ostream& os,
    const Foam::CompactIOField<T, BaseType>& L
)
{
    // Keep ASCII writing same
    if (os.format() == IOstreamOption::ASCII)
    {
        os << static_cast<const Field<T>&>(L);
    }
    else
    {
        // Convert to compact format
        labelList start(L.size()+1);

        start[0] = 0;
        for (label i = 1; i < start.size(); i++)
        {
            start[i] = start[i-1]+L[i-1].size();
        }

        Field<BaseType> elems(start[start.size()-1]);

        label elemI = 0;
        forAll(L, i)
        {
            const T& subField = L[i];

            forAll(subField, j)
            {
                elems[elemI++] = subField[j];
            }
        }
        os << start << elems;
    }

    return os;
}


// ************************************************************************* //
