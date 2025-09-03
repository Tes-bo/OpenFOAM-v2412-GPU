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

#include "IOMap.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class T>
bool Foam::IOMap<T>::readIOcontents()
{
    if (isReadRequired() || (isReadOptional() && headerOk()))
    {
        // For if MUST_READ_IF_MODIFIED
        addWatch();

        readStream(typeName) >> *this;
        close();

        return true;
    }

    return false;
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template<class T>
Foam::IOMap<T>::IOMap(const IOobject& io)
:
    regIOobject(io)
{
    readIOcontents();
}


template<class T>
Foam::IOMap<T>::IOMap(const IOobject& io, const label size)
:
    regIOobject(io)
{
    if (!readIOcontents())
    {
        Map<T>::resize(size);
    }
}


template<class T>
Foam::IOMap<T>::IOMap(const IOobject& io, const Map<T>& content)
:
    regIOobject(io)
{
    if (!readIOcontents())
    {
        Map<T>::operator=(content);
    }
}


template<class T>
Foam::IOMap<T>::IOMap(const IOobject& io, Map<T>&& content)
:
    regIOobject(io)
{
    Map<T>::transfer(content);

    readIOcontents();
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class T>
Foam::Map<T> Foam::IOMap<T>::readContents(const IOobject& io)
{
    IOobject rio(io, IOobjectOption::NO_REGISTER);
    if (rio.readOpt() == IOobjectOption::READ_MODIFIED)
    {
        rio.readOpt(IOobjectOption::MUST_READ);
    }

    // The object is global
    rio.globalObject(true);

    IOMap<T> reader(rio);

    return Map<T>(std::move(static_cast<Map<T>&>(reader)));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
bool Foam::IOMap<T>::writeData(Ostream& os) const
{
    os << *this;
    return os.good();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
void Foam::IOMap<T>::operator=(const IOMap<T>& rhs)
{
    Map<T>::operator=(rhs);
}


// ************************************************************************* //
