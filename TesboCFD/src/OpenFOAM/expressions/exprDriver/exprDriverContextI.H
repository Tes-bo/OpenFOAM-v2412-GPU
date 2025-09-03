/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::expressions::exprDriver::hasContextObjects() const
{
    return !contextObjects_.empty();
}


inline const Foam::regIOobject*
Foam::expressions::exprDriver::cfindContextIOobject
(
    const word& name
) const
{
    // Like objectRegistry::cfindIOobject()
    return contextObjects_.lookup(name, nullptr);
}


inline void Foam::expressions::exprDriver::addContextObject
(
    const word& name,
    const regIOobject* objptr
)
{
    if (objptr)
    {
        contextObjects_.set(name, objptr);
    }
    else
    {
        contextObjects_.erase(name);
    }
}


inline void Foam::expressions::exprDriver::addContextObject
(
    const regIOobject* objptr
)
{
    if (objptr)
    {
        addContextObject(objptr->name(), objptr);
    }
}


inline void Foam::expressions::exprDriver::removeContextObject
(
    const word& name
)
{
    contextObjects_.erase(name);
}


inline void Foam::expressions::exprDriver::removeContextObject
(
    const regIOobject* objptr
)
{
    if (objptr)
    {
        contextObjects_.erase(objptr->name());
    }
}


inline const Foam::HashTable<const Foam::regIOobject*>&
Foam::expressions::exprDriver::contextObjects() const noexcept
{
    return contextObjects_;
}


inline Foam::HashTable<const Foam::regIOobject*>&
Foam::expressions::exprDriver::contextObjects() noexcept
{
    return contextObjects_;
}


template<class ObjType>
const ObjType*
Foam::expressions::exprDriver::cfindContextObject(const word& name) const
{
    // Like objectRegistry::cfindObject()
    return dynamic_cast<const ObjType*>(this->cfindContextIOobject(name));
}


// ************************************************************************* //
