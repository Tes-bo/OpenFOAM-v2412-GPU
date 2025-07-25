/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2016-2024 OpenCFD Ltd.
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

#include "GlobalIOList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::GlobalIOList<Type>::GlobalIOList(const IOobject& io)
:
    regIOobject(io)
{
    // Check for MUST_READ_IF_MODIFIED
    warnNoRereading<GlobalIOList<Type>>();

    readHeaderOk(IOstreamOption::BINARY, typeName);
}


template<class Type>
Foam::GlobalIOList<Type>::GlobalIOList(const IOobject& io, Foam::zero)
:
    regIOobject(io)
{
    // Check for MUST_READ_IF_MODIFIED
    warnNoRereading<GlobalIOList<Type>>();

    readHeaderOk(IOstreamOption::BINARY, typeName);
}


template<class Type>
Foam::GlobalIOList<Type>::GlobalIOList(const IOobject& io, const label len)
:
    regIOobject(io)
{
    // Check for MUST_READ_IF_MODIFIED
    warnNoRereading<GlobalIOList<Type>>();

    if (!readHeaderOk(IOstreamOption::BINARY, typeName))
    {
        List<Type>::resize(len);
    }
}


template<class Type>
Foam::GlobalIOList<Type>::GlobalIOList
(
    const IOobject& io,
    const UList<Type>& content
)
:
    regIOobject(io)
{
    // Check for MUST_READ_IF_MODIFIED
    warnNoRereading<GlobalIOList<Type>>();

    if (!readHeaderOk(IOstreamOption::BINARY, typeName))
    {
        List<Type>::operator=(content);
    }
}


template<class Type>
Foam::GlobalIOList<Type>::GlobalIOList
(
    const IOobject& io,
    List<Type>&& content
)
:
    regIOobject(io)
{
    // Check for MUST_READ_IF_MODIFIED
    warnNoRereading<GlobalIOList<Type>>();

    List<Type>::transfer(content);

    readHeaderOk(IOstreamOption::BINARY, typeName);
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class Type>
Foam::List<Type> Foam::GlobalIOList<Type>::readContents(const IOobject& io)
{
    IOobject rio(io, IOobjectOption::NO_REGISTER);
    if (rio.readOpt() == IOobjectOption::READ_MODIFIED)
    {
        rio.readOpt(IOobjectOption::MUST_READ);
    }

    // The object is global
    rio.globalObject(true);

    GlobalIOList<Type> reader(rio);

    return List<Type>(std::move(static_cast<List<Type>&>(reader)));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::GlobalIOList<Type>::readData(Istream& is)
{
    is >> *this;
    return is.good();
}


template<class Type>
bool Foam::GlobalIOList<Type>::writeData(Ostream& os) const
{
    os << static_cast<const List<Type>&>(*this);
    return os.good();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::GlobalIOList<Type>::operator=(const GlobalIOList<Type>& rhs)
{
    List<Type>::operator=(rhs);
}


// ************************************************************************* //
