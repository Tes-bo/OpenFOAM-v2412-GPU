/*---------------------------------------------------------------------------* \
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018-2023 OpenCFD Ltd.
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

#include "IOPtrList.H"

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template<class T>
template<class INew>
Foam::IOPtrList<T>::IOPtrList(const IOobject& io, const INew& inewt)
:
    regIOobject(io)
{
    if (isReadRequired() || (isReadOptional() && headerOk()))
    {
        // For if MUST_READ_IF_MODIFIED
        addWatch();

        PtrList<T>::readIstream(readStream(typeName), inewt);
        close();
    }
}


template<class T>
Foam::IOPtrList<T>::IOPtrList(const IOobject& io)
:
    regIOobject(io)
{
    if (isReadRequired() || (isReadOptional() && headerOk()))
    {
        // For if MUST_READ_IF_MODIFIED
        addWatch();

        PtrList<T>::readIstream(readStream(typeName), INew<T>());
        close();
    }
}


template<class T>
Foam::IOPtrList<T>::IOPtrList(const IOobject& io, const label len)
:
    regIOobject(io),
    PtrList<T>(len)
{
    if (io.isAnyRead())
    {
        FatalErrorInFunction
            << "NO_READ must be set if specifying size" << nl
            << exit(FatalError);
    }
}


template<class T>
Foam::IOPtrList<T>::IOPtrList(const IOobject& io, const PtrList<T>& content)
:
    regIOobject(io)
{
    if (isReadRequired() || (isReadOptional() && headerOk()))
    {
        // For if MUST_READ_IF_MODIFIED
        addWatch();

        PtrList<T>::readIstream(readStream(typeName), INew<T>());
        close();
    }
    else
    {
        PtrList<T>::operator=(content);
    }
}


template<class T>
Foam::IOPtrList<T>::IOPtrList(const IOobject& io, PtrList<T>&& content)
:
    regIOobject(io)
{
    PtrList<T>::transfer(content);

    if (isReadRequired() || (isReadOptional() && headerOk()))
    {
        // For if MUST_READ_IF_MODIFIED
        addWatch();

        PtrList<T>::readIstream(readStream(typeName), INew<T>());
        close();
    }
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class T>
Foam::PtrList<T> Foam::IOPtrList<T>::readContents(const IOobject& io)
{
    IOobject rio(io, IOobjectOption::NO_REGISTER);
    if (rio.readOpt() == IOobjectOption::READ_MODIFIED)
    {
        rio.readOpt(IOobjectOption::MUST_READ);
    }

    IOPtrList<T> reader(rio);

    return PtrList<T>(std::move(static_cast<PtrList<T>&>(reader)));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
bool Foam::IOPtrList<T>::writeData(Ostream& os) const
{
    return (os << *this).good();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
void Foam::IOPtrList<T>::operator=(const IOPtrList<T>& rhs)
{
    PtrList<T>::operator=(rhs);
}


// ************************************************************************* //
