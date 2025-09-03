/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022-2024 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::IOobjectList::IOobjectList(const Foam::zero) noexcept
:
    HashPtrTable<IOobject>(Foam::zero{})
{}


inline Foam::IOobjectList::IOobjectList(const label initialCapacity)
:
    HashPtrTable<IOobject>(initialCapacity)
{}


inline Foam::IOobjectList::IOobjectList(const IOobjectList& list)
:
    HashPtrTable<IOobject>(list)
{}


inline Foam::IOobjectList::IOobjectList(IOobjectList&& list)
:
    HashPtrTable<IOobject>(std::move(list))
{}


inline Foam::IOobjectList::IOobjectList
(
    const objectRegistry& db,
    const fileName& instance,
    IOobjectOption ioOpt
)
:
    IOobjectList(db, instance, fileName::null, ioOpt)
{}


inline Foam::IOobjectList::IOobjectList
(
    const objectRegistry& db,
    const fileName& instance,
    IOobjectOption::registerOption registerObject
)
:
    IOobjectList
    (
        db,
        instance,
        fileName::null,
        IOobjectOption
        (
            IOobjectOption::MUST_READ,
            IOobjectOption::NO_WRITE,
            registerObject
        )
    )
{}


inline Foam::IOobjectList::IOobjectList
(
    const objectRegistry& db,
    const fileName& instance,
    const fileName& local,
    IOobjectOption::registerOption registerObject
)
:
    IOobjectList
    (
        db,
        instance,
        local,
        IOobjectOption
        (
            IOobjectOption::MUST_READ,
            IOobjectOption::NO_WRITE,
            registerObject
        )
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::IOobjectList::add(std::unique_ptr<IOobject>&& objectPtr)
{
    if (objectPtr)
    {
        return insert(objectPtr->name(), std::move(objectPtr));
    }

    return false;
}


inline bool Foam::IOobjectList::add(autoPtr<IOobject>& objectPtr)
{
    if (objectPtr)
    {
        return insert(objectPtr->name(), std::move(objectPtr));
    }

    return false;
}


inline bool Foam::IOobjectList::add(autoPtr<IOobject>&& objectPtr)
{
    if (objectPtr)
    {
        return insert(objectPtr->name(), std::move(objectPtr));
    }

    return false;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline void Foam::IOobjectList::operator=(IOobjectList&& list)
{
    transfer(list);
}


// ************************************************************************* //
