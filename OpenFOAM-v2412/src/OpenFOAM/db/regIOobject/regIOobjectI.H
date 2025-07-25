/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::regIOobject::registered() const noexcept
{
    return registered_;
}


inline bool Foam::regIOobject::ownedByRegistry() const noexcept
{
    return ownedByRegistry_;
}


inline bool Foam::regIOobject::store()
{
    if (checkIn())
    {
        ownedByRegistry_ = true;
    }
    else
    {
        WarningInFunction
            << "Refuse to store unregistered object: " << this->name() << nl;
    }

    return ownedByRegistry_;
}


template<class Type>
inline Type& Foam::regIOobject::store(Type* p)
{
    if (!p)
    {
        FatalErrorInFunction
            << "Object deallocated\n"
            << abort(FatalError);
    }

    const bool ok = p->regIOobject::store();

    if (!ok)
    {
        FatalErrorInFunction
            << "Failed to store pointer: " << p->regIOobject::name()
            << ". Risk of memory leakage\n"
            << abort(FatalError);
    }

    return *p;
}


template<class Type>
inline Type& Foam::regIOobject::store(std::unique_ptr<Type>&& ptr)
{
    // Pass management to objectRegistry
    return store(ptr.release());
}


template<class Type>
inline Type& Foam::regIOobject::store(autoPtr<Type>& ptr)
{
    // Pass management to objectRegistry
    return store(ptr.release());
}


template<class Type>
inline Type& Foam::regIOobject::store(autoPtr<Type>&& ptr)
{
    // Pass management to objectRegistry
    return store(ptr.release());
}


template<class Type>
inline Type& Foam::regIOobject::store(refPtr<Type>& ptr)
{
    Type* p = nullptr;

    if (ptr.is_pointer())
    {
        // Acquire ownership, pass management to objectRegistry
        p = ptr.release();

        store(p);

        // Change parameter to access the stored reference
        ptr.cref(p);
    }
    else
    {
        // Taking ownership of reference does not make much sense.
        // - Storing the object won't actually do so, it will be removed
        //   when the original object goes out of scope.
        // - Storing a clone may not be what we want.

        p = ptr.get();

        WarningInFunction
            << "Refuse to store reference: " << p->name()
            << ". Likely indicates a coding error\n";
    }

    return *p;
}


template<class Type>
inline Type& Foam::regIOobject::store(refPtr<Type>&& ptr)
{
    // Forward as named reference, which also does a move
    return store(ptr);
}


template<class Type>
inline Type& Foam::regIOobject::store(tmp<Type>& ptr)
{
    Type* p = nullptr;

    if (ptr.is_pointer())
    {
        // Acquire ownership, pass management to objectRegistry

        ptr.protect(false);  // Storing (ie, not cached/protected)
        p = ptr.ptr();

        store(p);

        // Change parameter to access the stored reference
        ptr.cref(p);
    }
    else
    {
        // Taking ownership of reference does not make much sense.
        // - Storing the object won't actually do so, it will be removed
        //   when the original object goes out of scope.
        // - Storing a clone may not be what we want.

        p = ptr.get();

        WarningInFunction
            << "Refuse to store reference: " << p->name()
            << ". Likely indicates a coding error\n";
    }

    return *p;
}


template<class Type>
inline Type& Foam::regIOobject::store(tmp<Type>&& ptr)
{
    // Forward as named reference, which also does a move
    return store(ptr);
}


inline void Foam::regIOobject::release(const bool unregister) noexcept
{
    // Note: could also return the old ownedByRegistry_ value
    ownedByRegistry_ = false;
    if (unregister)
    {
        registered_ = false;
    }
}


inline Foam::label Foam::regIOobject::eventNo() const noexcept
{
    return eventNo_;
}

inline Foam::label& Foam::regIOobject::eventNo() noexcept
{
    return eventNo_;
}


inline const Foam::labelList& Foam::regIOobject::watchIndices() const noexcept
{
    return watchIndices_;
}


inline Foam::labelList& Foam::regIOobject::watchIndices() noexcept
{
    return watchIndices_;
}


// ************************************************************************* //
