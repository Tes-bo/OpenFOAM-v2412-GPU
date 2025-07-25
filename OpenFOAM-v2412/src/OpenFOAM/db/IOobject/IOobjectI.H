/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2023 OpenCFD Ltd.
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

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class StringType>
inline Foam::word Foam::IOobject::groupName
(
    StringType base,
    const word& group
)
{
    if (group.empty())
    {
        return base;
    }

    return base + ('.' + group);
}


inline Foam::word Foam::IOobject::scopedName
(
    const std::string& scope,
    const word& name
)
{
    if (scope.empty())
    {
        return name;
    }

    word output;
    output.reserve(scope.size() + name.size() + 1);

    output += scope;
    output += IOobject::scopeSeparator;
    output += name;
    return output;
}


inline Foam::word Foam::IOobject::scopedName
(
    const std::string& scope,
    const word& name1,
    const word& name2
)
{
    if (scope.empty())
    {
        return IOobject::scopedName(name1, name2);
    }

    word output;
    output.reserve(scope.size() + name1.size() + name2.size() + 2);

    output += scope;
    output += IOobject::scopeSeparator;
    output += name1;
    if (!name2.empty())
    {
        output += IOobject::scopeSeparator;
        output += name2;
    }
    return output;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::IOobject::IOobject
(
    const word& name,
    const fileName& instance,
    const objectRegistry& registry,
    IOobjectOption::readOption rOpt,
    IOobjectOption::writeOption wOpt,
    bool registerObject,
    bool globalObject
)
:
    IOobject
    (
        name,
        instance,
        registry,
        IOobjectOption(rOpt, wOpt, registerObject, globalObject)
    )
{}


inline Foam::IOobject::IOobject
(
    const word& name,
    const fileName& instance,
    const fileName& local,
    const objectRegistry& registry,
    IOobjectOption::readOption rOpt,
    IOobjectOption::writeOption wOpt,
    bool registerObject,
    bool globalObject
)
:
    IOobject
    (
        name,
        instance,
        local,
        registry,
        IOobjectOption(rOpt, wOpt, registerObject, globalObject)
    )
{}


inline Foam::IOobject::IOobject
(
    const fileName& path,
    const objectRegistry& registry,
    IOobjectOption::readOption rOpt,
    IOobjectOption::writeOption wOpt,
    bool registerObject,
    bool globalObject
)
:
    IOobject
    (
        path,
        registry,
        IOobjectOption(rOpt, wOpt, registerObject, globalObject)
    )
{}


inline Foam::IOobject::IOobject
(
    const IOobject& io,
    const word& name,
    const fileName& local
)
:
    IOobject(io, name)
{
    local_ = local;
}


inline Foam::IOobject::IOobject
(
    const IOobject& io,
    IOobjectOption::readOption rOpt,
    IOobjectOption::writeOption wOpt
)
:
    IOobject(io)
{
    IOobjectOption::readOpt(rOpt);
    IOobjectOption::writeOpt(wOpt);
}


inline Foam::IOobject::IOobject
(
    const IOobject& io,
    IOobjectOption::registerOption regOpt
)
:
    IOobject(io)
{
    IOobjectOption::registerObject(regOpt);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// General access

inline const Foam::word& Foam::IOobject::name() const noexcept
{
    return name_;
}


inline Foam::word Foam::IOobject::group() const
{
    return name_.ext();
}


inline Foam::word Foam::IOobject::member() const
{
    return name_.lessExt();
}


inline const Foam::word& Foam::IOobject::headerClassName() const noexcept
{
    return headerClassName_;
}


inline Foam::word& Foam::IOobject::headerClassName() noexcept
{
    return headerClassName_;
}


inline const Foam::string& Foam::IOobject::note() const noexcept
{
    return note_;
}


inline Foam::string& Foam::IOobject::note() noexcept
{
    return note_;
}


inline unsigned Foam::IOobject::labelByteSize() const noexcept
{
    return static_cast<unsigned>(sizeofLabel_);
}


inline unsigned Foam::IOobject::scalarByteSize() const noexcept
{
    return static_cast<unsigned>(sizeofScalar_);
}


// Checks

inline bool Foam::IOobject::hasHeaderClass() const noexcept
{
    return !headerClassName_.empty();
}


template<class Type>
inline bool Foam::IOobject::isHeaderClass() const
{
    return (Type::typeName == headerClassName_);
}


// Path components

inline const Foam::fileName& Foam::IOobject::instance() const noexcept
{
    return instance_;
}


inline Foam::fileName& Foam::IOobject::instance() noexcept
{
    return instance_;
}


inline const Foam::fileName& Foam::IOobject::local() const noexcept
{
    return local_;
}


inline Foam::fileName Foam::IOobject::objectPath() const
{
    return path()/name();
}


inline Foam::fileName Foam::IOobject::globalObjectPath() const
{
    return globalPath()/name();
}


// Error Handling

inline bool Foam::IOobject::good() const noexcept
{
    return objState_ == objectState::GOOD;
}


inline bool Foam::IOobject::bad() const noexcept
{
    return objState_ == objectState::BAD;
}


// ************************************************************************* //
