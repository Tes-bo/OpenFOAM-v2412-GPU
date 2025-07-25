/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2019 OpenFOAM Foundation
    Copyright (C) 2016-2023 OpenCFD Ltd.
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

#include "objectRegistry.H"
#include "predicates.H"
#include <type_traits>

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Templated implementation for classes()
template<class MatchPredicate>
Foam::HashTable<Foam::wordHashSet> Foam::objectRegistry::classesImpl
(
    const objectRegistry& list,
    const MatchPredicate& matchName
)
{
    HashTable<wordHashSet> summary(2*list.size());

    // Summary (key,val) = (class-name, object-names)
    forAllConstIters(list, iter)
    {
        const regIOobject* obj = iter.val();

        if (matchName(obj->name()))
        {
            // Create entry (if needed) and insert
            summary(obj->type()).insert(obj->name());
        }
    }

    return summary;
}


// Templated implementation for count()
template<class MatchPredicate1, class MatchPredicate2>
Foam::label Foam::objectRegistry::countImpl
(
    const objectRegistry& list,
    const MatchPredicate1& matchClass,
    const MatchPredicate2& matchName
)
{
    label count = 0;

    forAllConstIters(list, iter)
    {
        const regIOobject* obj = iter.val();

        if (matchClass(obj->type()) && matchName(obj->name()))
        {
            ++count;
        }
    }

    return count;
}


// Templated implementation for count()
template<class Type, class MatchPredicate>
Foam::label Foam::objectRegistry::countTypeImpl
(
    const objectRegistry& list,
    const MatchPredicate& matchName
)
{
    label count = 0;

    forAllConstIters(list, iter)
    {
        const regIOobject* obj = iter.val();

        if
        (
            (std::is_void<Type>::value || Foam::isA<Type>(*obj))
         && matchName(obj->name())
        )
        {
            ++count;
        }
    }

    return count;
}


// Templated implementation for names(), sortedNames()
template<class MatchPredicate1, class MatchPredicate2>
Foam::wordList Foam::objectRegistry::namesImpl
(
    const objectRegistry& list,
    const MatchPredicate1& matchClass,
    const MatchPredicate2& matchName,
    const bool doSort
)
{
    wordList objNames(list.size());

    label count=0;
    forAllConstIters(list, iter)
    {
        const regIOobject* obj = iter.val();

        if (matchClass(obj->type()) && matchName(obj->name()))
        {
            objNames[count] = obj->name();
            ++count;
        }
    }

    objNames.resize(count);

    if (doSort)
    {
        Foam::sort(objNames);
    }

    return objNames;
}


// Templated implementation for names(), sortedNames()
template<class Type, class MatchPredicate>
Foam::wordList Foam::objectRegistry::namesTypeImpl
(
    const objectRegistry& list,
    const MatchPredicate& matchName,
    const bool doSort
)
{
    wordList objNames(list.size());

    label count = 0;
    forAllConstIters(list, iter)
    {
        const regIOobject* obj = iter.val();

        if
        (
            (std::is_void<Type>::value || Foam::isA<Type>(*obj))
         && matchName(obj->name())
        )
        {
            objNames[count] = obj->name();
            ++count;
        }
    }

    objNames.resize(count);

    if (doSort)
    {
        Foam::sort(objNames);
    }

    return objNames;
}


// Templated implementation for cobjects()/objects(), csorted()/sorted()
template<class Type, class MatchPredicate>
Foam::UPtrList<Type>
Foam::objectRegistry::objectsTypeImpl
(
    const bool strict,
    const objectRegistry& list,
    const MatchPredicate& matchName,
    const bool doSort
)
{
    typedef typename std::remove_cv<Type>::type BaseType;

    UPtrList<Type> result(list.size());

    label count = 0;
    forAllConstIters(list, iter)
    {
        const regIOobject* obj = iter.val();
        const BaseType* ptr = dynamic_cast<const BaseType*>(obj);

        if
        (
            ptr
         && (!strict || Foam::isType<BaseType>(*obj))
         && matchName(obj->name())
        )
        {
            result.set(count, const_cast<BaseType*>(ptr));
            ++count;
        }
    }

    result.resize(count);

    if (doSort)
    {
        Foam::sort(result, nameOp<Type>());  // Sort by object name()
    }

    return result;
}


// Templated implementation for lookupClass()
template<class Type>
Foam::HashTable<Type*>
Foam::objectRegistry::lookupClassTypeImpl
(
    const bool strict,
    const objectRegistry& list
)
{
    typedef typename std::remove_cv<Type>::type BaseType;

    HashTable<Type*> result(list.capacity());

    forAllConstIters(list, iter)
    {
        const regIOobject* obj = iter.val();
        const BaseType* ptr = dynamic_cast<const BaseType*>(obj);

        if
        (
            ptr
         && (!strict || Foam::isType<BaseType>(*obj))
        )
        {
            result.insert(obj->name(), const_cast<BaseType*>(ptr));
        }
    }

    return result;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class MatchPredicate>
Foam::HashTable<Foam::wordHashSet>
Foam::objectRegistry::classes
(
    const MatchPredicate& matchName
) const
{
    return classesImpl(*this, matchName);
}


template<class MatchPredicate>
Foam::label Foam::objectRegistry::count
(
    const MatchPredicate& matchClass
) const
{
    return countImpl(*this, matchClass, predicates::always());
}


template<class MatchPredicate1, class MatchPredicate2>
Foam::label Foam::objectRegistry::count
(
    const MatchPredicate1& matchClass,
    const MatchPredicate2& matchName
) const
{
    return countImpl(*this, matchClass, matchName);
}


template<class Type, class MatchPredicate>
Foam::label Foam::objectRegistry::count
(
    const MatchPredicate& matchName
) const
{
    return countTypeImpl<Type>(*this, matchName);
}


template<class Type>
Foam::label Foam::objectRegistry::count
(
    const bool strict
) const
{
    label nObjects = 0;

    forAllConstIters(*this, iter)
    {
        const regIOobject* obj = iter.val();

        if
        (
            std::is_void<Type>::value
         ||
            (
                strict
              ? bool(Foam::isType<Type>(*obj))
              : bool(Foam::isA<Type>(*obj))
            )
        )
        {
            ++nObjects;
        }
    }

    return nObjects;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, bool Strict>
Foam::UPtrList<const Type>
Foam::objectRegistry::cobjects() const
{
    return objectsTypeImpl<const Type>
    (
        Strict, *this, predicates::always(), false  // doSort = false
    );
}


template<class Type, bool Strict>
Foam::UPtrList<Type>
Foam::objectRegistry::objects()
{
    return objectsTypeImpl<Type>
    (
        Strict, *this, predicates::always(), false  // doSort = false
    );
}


template<class Type, bool Strict>
Foam::UPtrList<const Type>
Foam::objectRegistry::csorted() const
{
    return objectsTypeImpl<const Type>
    (
        Strict, *this, predicates::always(), true  // doSort = true
    );
}


template<class Type, bool Strict>
Foam::UPtrList<Type>
Foam::objectRegistry::sorted()
{
    return objectsTypeImpl<Type>
    (
        Strict, *this, predicates::always(), true  // doSort = false
    );
}


template<class Type, class MatchPredicate>
Foam::UPtrList<const Type>
Foam::objectRegistry::cobjects
(
    const MatchPredicate& matchName
) const
{
    // doSort = false
    return objectsTypeImpl<const Type>(false, *this, matchName, false);
}


template<class Type, class MatchPredicate>
Foam::UPtrList<Type>
Foam::objectRegistry::objects
(
    const MatchPredicate& matchName
)
{
    // doSort = false
    return objectsTypeImpl<Type>(false, *this, matchName, false);
}


template<class Type, class MatchPredicate>
Foam::UPtrList<const Type>
Foam::objectRegistry::csorted
(
    const MatchPredicate& matchName
) const
{
    return objectsTypeImpl<const Type>(false, *this, matchName, true);
}


template<class Type, class MatchPredicate>
Foam::UPtrList<Type>
Foam::objectRegistry::sorted
(
    const MatchPredicate& matchName
)
{
    return objectsTypeImpl<Type>(false, *this, matchName, true);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class MatchPredicate>
Foam::wordList Foam::objectRegistry::names
(
    const MatchPredicate& matchClass
) const
{
    return namesImpl(*this, matchClass, predicates::always(), false);
}


template<class MatchPredicate1, class MatchPredicate2>
Foam::wordList Foam::objectRegistry::names
(
    const MatchPredicate1& matchClass,
    const MatchPredicate2& matchName
) const
{
    return namesImpl(*this, matchClass, matchName, false);
}


template<class Type>
Foam::wordList Foam::objectRegistry::names() const
{
    return namesTypeImpl<Type>(*this, predicates::always(), false);
}


template<class Type, class MatchPredicate>
Foam::wordList Foam::objectRegistry::names
(
    const MatchPredicate& matchName
) const
{
    return namesTypeImpl<Type>(*this, matchName, false);
}


template<class MatchPredicate>
Foam::wordList Foam::objectRegistry::sortedNames
(
    const MatchPredicate& matchClass
) const
{
    return namesImpl(*this, matchClass, predicates::always(), true);
}


template<class MatchPredicate1, class MatchPredicate2>
Foam::wordList Foam::objectRegistry::sortedNames
(
    const MatchPredicate1& matchClass,
    const MatchPredicate2& matchName
) const
{
    return namesImpl(*this, matchClass, matchName, true);
}


template<class Type>
Foam::wordList Foam::objectRegistry::sortedNames() const
{
    return namesTypeImpl<Type>(*this, predicates::always(), true);
}


template<class Type, class MatchPredicate>
Foam::wordList Foam::objectRegistry::sortedNames
(
    const MatchPredicate& matchName
) const
{
    return namesTypeImpl<Type>(*this, matchName, true);
}


template<class Type, bool Strict>
Foam::HashTable<const Type*> Foam::objectRegistry::lookupClass() const
{
    return lookupClassTypeImpl<const Type>(Strict, *this);
}


template<class Type, bool Strict>
Foam::HashTable<Type*> Foam::objectRegistry::lookupClass()
{
    return lookupClassTypeImpl<Type>(Strict, *this);
}


template<class Type>
Foam::HashTable<const Type*> Foam::objectRegistry::lookupClass
(
    const bool strict
) const
{
    return lookupClassTypeImpl<const Type>(strict, *this);
}


template<class Type>
Foam::HashTable<Type*> Foam::objectRegistry::lookupClass
(
    const bool strict
)
{
    return lookupClassTypeImpl<Type>(strict, *this);
}


template<class Type>
bool Foam::objectRegistry::foundObject
(
    const word& name,
    const bool recursive
) const
{
    return this->cfindObject<Type>(name, recursive);
}


template<class Type>
const Type* Foam::objectRegistry::cfindObject
(
    const word& name,
    const bool recursive
) const
{
    return dynamic_cast<const Type*>(this->cfindIOobject(name, recursive));
}


template<class Type>
const Type* Foam::objectRegistry::findObject
(
    const word& name,
    const bool recursive
) const
{
    return this->cfindObject<Type>(name, recursive);
}


template<class Type>
Type* Foam::objectRegistry::findObject
(
    const word& name,
    const bool recursive
)
{
    return const_cast<Type*>(this->cfindObject<Type>(name, recursive));
}


template<class Type>
Type* Foam::objectRegistry::getObjectPtr
(
    const word& name,
    const bool recursive
) const
{
    return const_cast<Type*>(this->cfindObject<Type>(name, recursive));
}


template<class Type>
const Type& Foam::objectRegistry::lookupObject
(
    const word& name,
    const bool recursive
) const
{
    const_iterator iter = cfind(name);

    if (iter.good())
    {
        const Type* ptr = dynamic_cast<const Type*>(iter.val());

        if (ptr)
        {
            return *ptr;
        }

        FatalErrorInFunction
            << nl
            << "    bad lookup of " << name << " (objectRegistry "
            << this->name()
            << ")\n    expected a " << Type::typeName
            << ", found a " << (*iter)->type() << nl
            << exit(FatalError);
    }
    else if (recursive && !this->isTimeDb())
    {
        return parent_.lookupObject<Type>(name, recursive);
    }

    FatalErrorInFunction
        << nl
        << "    failed lookup of " << name << " (objectRegistry "
        << this->name()
        << ")\n    available objects of type " << Type::typeName
        << ':' << nl
        << names<Type>() << nl
        << exit(FatalError);

    return NullObjectRef<Type>();
}


template<class Type>
Type& Foam::objectRegistry::lookupObjectRef
(
    const word& name,
    const bool recursive
) const
{
    const Type& ref = this->lookupObject<Type>(name, recursive);
    // The above will already fail if things didn't work

    return const_cast<Type&>(ref);
}


template<class Type>
bool Foam::objectRegistry::cacheTemporaryObject(Type& obj) const
{
    bool ok = false;

    readCacheTemporaryObjects();

    if (cacheTemporaryObjects_.size())
    {
        temporaryObjects_.insert(obj.name());

        auto iter = cacheTemporaryObjects_.find(obj.name());

        // Cache object if is in the cacheTemporaryObjects list
        // and hasn't been cached yet
        if (iter.good() && iter.val().first() == false)
        {
            iter.val().first() = true;
            iter.val().second() = true;

            Type* cachedPtr = obj.db().template getObjectPtr<Type>(obj.name());

            // Remove any name collisions from the cache
            if (cachedPtr && cachedPtr != &obj && cachedPtr->ownedByRegistry())
            {
                deleteCachedObject(cachedPtr);
            }

            if (debug)
            {
                Info<< "Caching " << obj.name()
                    << " of type " << obj.type() << endl;
            }

            obj.release();
            obj.checkOut();
            regIOobject::store(new Type(std::move(obj)));

            ok = true;
        }
    }

    return ok;
}


// ************************************************************************* //
