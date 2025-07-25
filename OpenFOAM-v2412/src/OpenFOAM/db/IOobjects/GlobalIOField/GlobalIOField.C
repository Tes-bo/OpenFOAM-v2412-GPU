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

#include "GlobalIOField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::GlobalIOField<Type>::GlobalIOField(const IOobject& io)
:
    regIOobject(io)
{
    // Check for MUST_READ_IF_MODIFIED
    warnNoRereading<GlobalIOField<Type>>();

    readHeaderOk(IOstreamOption::BINARY, typeName);
}


template<class Type>
Foam::GlobalIOField<Type>::GlobalIOField(const IOobject& io, const label len)
:
    regIOobject(io)
{
    // Check for MUST_READ_IF_MODIFIED
    warnNoRereading<GlobalIOField<Type>>();

    if (!readHeaderOk(IOstreamOption::BINARY, typeName))
    {
        Field<Type>::resize(len);
    }
}


template<class Type>
Foam::GlobalIOField<Type>::GlobalIOField
(
    const IOobject& io,
    const UList<Type>& content
)
:
    regIOobject(io)
{
    // Check for MUST_READ_IF_MODIFIED
    warnNoRereading<GlobalIOField<Type>>();

    if (!readHeaderOk(IOstreamOption::BINARY, typeName))
    {
        Field<Type>::operator=(content);
    }
}


template<class Type>
Foam::GlobalIOField<Type>::GlobalIOField
(
    const IOobject& io,
    Field<Type>&& content
)
:
    regIOobject(io)
{
    // Check for MUST_READ_IF_MODIFIED
    warnNoRereading<GlobalIOField<Type>>();

    Field<Type>::transfer(content);

    readHeaderOk(IOstreamOption::BINARY, typeName);
}


template<class Type>
Foam::GlobalIOField<Type>::GlobalIOField
(
    const IOobject& io,
    const tmp<Field<Type>>& tfld
)
:
    regIOobject(io)
{
    // Check for MUST_READ_IF_MODIFIED
    warnNoRereading<GlobalIOField<Type>>();

    const bool reuse = tfld.movable();

    if (reuse)
    {
        Field<Type>::transfer(tfld.ref());
    }

    if (!readHeaderOk(IOstreamOption::BINARY, typeName) && !reuse)
    {
        Field<Type>::operator=(tfld());
    }

    tfld.clear();
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class Type>
Foam::Field<Type> Foam::GlobalIOField<Type>::readContents(const IOobject& io)
{
    IOobject rio(io, IOobjectOption::NO_REGISTER);
    if (rio.readOpt() == IOobjectOption::READ_MODIFIED)
    {
        rio.readOpt(IOobjectOption::MUST_READ);
    }

    // The object is global
    rio.globalObject(true);

    GlobalIOField<Type> reader(rio);

    return Field<Type>(std::move(static_cast<Field<Type>&>(reader)));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::GlobalIOField<Type>::readData(Istream& is)
{
    is >> *this;
    return is.good();
}


template<class Type>
bool Foam::GlobalIOField<Type>::writeData(Ostream& os) const
{
    os << static_cast<const Field<Type>&>(*this);
    return os.good();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::GlobalIOField<Type>::operator=(const GlobalIOField<Type>& rhs)
{
    Field<Type>::operator=(rhs);
}


// ************************************************************************* //
