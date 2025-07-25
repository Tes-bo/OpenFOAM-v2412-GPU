/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2019 OpenFOAM Foundation
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

Class
    Foam::objectRegistry

Description
    Registry of regIOobjects

SourceFiles
    objectRegistry.C
    objectRegistryTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_objectRegistry_H
#define Foam_objectRegistry_H

#include "HashTable.H"
#include "HashSet.H"
#include "UPtrList.H"
#include "regIOobject.H"
#include "wordRes.H"
#include "Pair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class objectRegistry Declaration
\*---------------------------------------------------------------------------*/

class objectRegistry
:
    public regIOobject,
    public HashTable<regIOobject*>
{
    // Private Data

        //- Master time objectRegistry
        const Time& time_;

        //- Parent objectRegistry
        const objectRegistry& parent_;

        //- Local directory path of this objectRegistry relative to time
        fileName dbDir_;

        //- Current event
        mutable label event_;

        //- State of cacheTemporaryObjects_, set true after reading
        mutable bool cacheTemporaryObjectsActive_;

        //- Names of temporary object with current state
        mutable HashTable<Pair<bool>> cacheTemporaryObjects_;

        //- Accumulated list of temporary objects available to cache
        //  Used to provide diagnostics in case the requested object is not
        //  available
        mutable wordHashSet temporaryObjects_;


    // Private Member Functions

        //- Is the objectRegistry parent_ different from time_
        //  Used to terminate searching within the ancestors
        bool parentNotTime() const noexcept;

        //- Read the cacheTemporaryObjects list from Time controlDict
        void readCacheTemporaryObjects() const;

        //- Delete the cached object. Eg, before caching a new object
        //- A nullptr is ignored.
        void deleteCachedObject(regIOobject* io) const;

        //- Templated implementation for count()
        //  The number of items with a matching class
        template<class MatchPredicate1, class MatchPredicate2>
        static label countImpl
        (
            const objectRegistry& list,
            const MatchPredicate1& matchClass,
            const MatchPredicate2& matchName
        );

        //- Templated implementation for count()
        //  The number of items with a matching class
        template<class Type, class MatchPredicate>
        static label countTypeImpl
        (
            const objectRegistry& list,
            const MatchPredicate& matchName
        );

        //- Templated implementation for classes()
        template<class MatchPredicate>
        static HashTable<wordHashSet> classesImpl
        (
            const objectRegistry& list,
            const MatchPredicate& matchName
        );

        //- Templated implementation for names(), sortedNames()
        template<class MatchPredicate1, class MatchPredicate2>
        static wordList namesImpl
        (
            const objectRegistry& list,
            const MatchPredicate1& matchClass,
            const MatchPredicate2& matchName,
            const bool doSort
        );

        //- Templated implementation for names(), sortedNames()
        template<class Type, class MatchPredicate>
        static wordList namesTypeImpl
        (
            const objectRegistry& list,
            const MatchPredicate& matchName,
            const bool doSort
        );

        //- Templated implementation for csorted()/sorted()
        //  Called with 'Type' or 'const Type'
        template<class Type, class MatchPredicate>
        static UPtrList<Type> objectsTypeImpl
        (
            const bool strict,         // Check with isType<Type>
            const objectRegistry& list,
            const MatchPredicate& matchName,
            const bool doSort          // Sort the list by name
        );

        //- Templated implementation for lookupClass()
        //  Called with 'Type' or 'const Type'
        template<class Type>
        static HashTable<Type*> lookupClassTypeImpl
        (
            const bool strict,  // Check with isType<Type>
            const objectRegistry& list
        );


        //- No copy construct
        objectRegistry(const objectRegistry&) = delete;

        //- No copy assignment
        void operator=(const objectRegistry&) = delete;


public:

    //- Declare type name for this IOobject
    TypeName("objectRegistry");


    // Constructors

        //- Construct the time objectRegistry,
        //- with estimated table capacity (default: 128)
        explicit objectRegistry
        (
            const Time& db,
            const label initialCapacity = 128
        );

        //- Construct sub-registry given an IObject to describe the registry,
        //- with estimated table capacity (default: 128)
        explicit objectRegistry
        (
            const IOobject& io,
            const label initialCapacity = 128
        );


    //- Destructor, with checkOut() for all objects that are ownedByRegistry
    virtual ~objectRegistry();


    // Member Functions

    // Access

        //- Return the object registry
        const objectRegistry& thisDb() const noexcept
        {
            return *this;
        }

        //- Return the parent objectRegistry
        const objectRegistry& parent() const noexcept
        {
            return parent_;
        }

        //- Return time registry
        const Time& time() const noexcept
        {
            return time_;
        }

        //- True if the registry is Time
        bool isTimeDb() const noexcept;

        //- Local directory path of this objectRegistry relative to the time
        virtual const fileName& dbDir() const
        {
            return dbDir_;
        }


    // Helper Functions

        //- Create an IOobject at the current time instance (timeName)
        //- with the specified options
        IOobject newIOobject
        (
            const word& name,       //!< The object name
            IOobjectOption ioOpt    //!< The IOobject options
        ) const;

        //- Create an IOobject at the current time instance (timeName).
        //  By default the object is NO_READ/NO_WRITE/NO_REGISTER
        IOobject newIOobject
        (
            //! The object name
            const word& name,
            //! The read option (default: NO_READ)
            IOobjectOption::readOption rOpt = IOobjectOption::NO_READ,
            //! The write option (default: NO_WRITE)
            IOobjectOption::writeOption wOpt = IOobjectOption::NO_WRITE,
            //! The register option (default: NO_REGISTER)
            IOobjectOption::registerOption regOpt = IOobjectOption::NO_REGISTER
        ) const;


    // Summary of classes

        //- A summary hash of classes used and their associated object names.
        //  Behaviour and usage as per IOobjectList::classes
        HashTable<wordHashSet> classes() const;

        //- A summary hash of classes used and their associated object names,
        //- restricted to objects that have a matching object name.
        template<class MatchPredicate>
        HashTable<wordHashSet> classes(const MatchPredicate& matchName) const;


    // List-wise access (unsorted)

        //- Return unsorted list of objects with a class satisfying
        //- \c isA\<Type\> or \c isType\<Type\> (with Strict)
        //  The lifetime of the returned content cannot exceed the parent!
        template<class Type, bool Strict=false>
        UPtrList<const Type> cobjects() const;

        //- Return unsorted list of objects with a class satisfying
        //- \c isA\<Type\> or \c isType\<Type\> (with Strict)
        //  The lifetime of the returned content cannot exceed the parent!
        template<class Type, bool Strict=false>
        UPtrList<Type> objects();

        //- Return unsorted list of objects with a class satisfying
        //- \c isA\<Type\> that also have a matching object name.
        //  The lifetime of the returned content cannot exceed the parent!
        template<class Type, class MatchPredicate>
        UPtrList<const Type> cobjects(const MatchPredicate& matchName) const;

        //- Return sorted list of objects with a class satisfying \c isA\<Type\>
        //- that also have a matching object name.
        //  The lifetime of the returned content cannot exceed the parent!
        template<class Type, class MatchPredicate>
        UPtrList<Type> objects(const MatchPredicate& matchName);


    // List-wise access (sorted)

        //- Return sorted list of objects with a class satisfying
        //- \c isA\<Type\> or \c isType\<Type\> (with Strict)
        //  The lifetime of the returned content cannot exceed the parent!
        template<class Type, bool Strict=false>
        UPtrList<const Type> csorted() const;

        //- Return sorted list of objects with a class satisfying
        //- \c isA\<Type\> or \c isType\<Type\> (with Strict)
        //  The lifetime of the returned content cannot exceed the parent!
        template<class Type, bool Strict=false>
        UPtrList<Type> sorted();

        //- Return sorted list of objects
        //  The lifetime of the returned content cannot exceed the parent!
        UPtrList<const regIOobject> csorted() const
        {
            return csorted<regIOobject>();
        }

        //- Return sorted list of objects
        //  The lifetime of the returned content cannot exceed the parent!
        UPtrList<regIOobject> sorted()
        {
            return sorted<regIOobject>();
        }

        //- Return sorted list of objects with a class satisfying
        //- \c isA\<Type\> that also have a matching object name.
        //  The lifetime of the returned content cannot exceed the parent!
        template<class Type, class MatchPredicate>
        UPtrList<const Type> csorted(const MatchPredicate& matchName) const;

        //- Return sorted list of objects with a class satisfying
        //- \c isA\<Type\> that also have a matching object name.
        //  The lifetime of the returned content cannot exceed the parent!
        template<class Type, class MatchPredicate>
        UPtrList<Type> sorted(const MatchPredicate& matchName);


    // Number of items

        //- The number of objects of the given class name
        //  \note uses the class type() method
        label count(const char* clsName) const;

        //- The number of objects of the given class name
        //  \note uses the class type() method
        template<class MatchPredicate>
        label count(const MatchPredicate& matchClass) const;

        //- The number of objects of the given class name
        //  \note uses the class type() method
        template<class MatchPredicate1, class MatchPredicate2>
        label count
        (
            const MatchPredicate1& matchClass,
            const MatchPredicate2& matchName
        ) const;

        //- The names of objects with a class satisfying \c isA\<Type\>
        //
        //  \param strict use \c isType\<Type\> instead of \c isA\<Type\>
        //
        //  \note The values of \c count\<Type\>() and \c count(Type::typeName)
        //      may be inconsistent, since they use different mechanisms for
        //      testing the class type.
        //  \note If \a Type is \c void, no isA check is used (always true).
        template<class Type>
        label count(const bool strict = false) const;

        //- The names of objects with a class satisfying \c isA\<Type\>
        //- that also have a matching object name.
        //
        //  \note If \a Type is \c void, no isA check is used (always true).
        template<class Type, class MatchPredicate>
        label count(const MatchPredicate& matchName) const;


    // Summary of names

        //- The unsorted names of all objects
        wordList names() const;

        //- The unsorted names of objects with the given class name.
        //  \note uses the class type() method
        wordList names(const char* clsName) const;

        //- The unsorted names of objects with a matching class name
        //  \note uses the class type() method
        template<class MatchPredicate>
        wordList names(const MatchPredicate& matchClass) const;

        //- The unsorted names of objects with a matching class name
        //- that also have a matching object name.
        //  \note uses the class type() method
        template<class MatchPredicate1, class MatchPredicate2>
        wordList names
        (
            const MatchPredicate1& matchClass,
            const MatchPredicate2& matchName
        ) const;

        //- The unsorted names of objects with a class satisfying \c isA\<Type\>
        //
        //  \note If \a Type is \c void, no isA check is used (always true).
        template<class Type>
        wordList names() const;

        //- The unsorted names of objects with a class satisfying \c isA\<Type\>
        //- that also have a matching object name.
        //
        //  \note If \a Type is \c void, no isA check is used (always true).
        template<class Type, class MatchPredicate>
        wordList names(const MatchPredicate& matchName) const;


    // Summary of names (sorted)

        //- The sorted names of all objects
        wordList sortedNames() const;

        //- The sorted names of objects with the given class name.
        //  \note uses the class type() method
        wordList sortedNames(const char* clsName) const;

        //- The sorted names objects with a matching class name
        //  \note uses the class type() method
        template<class MatchPredicate>
        wordList sortedNames(const MatchPredicate& matchClass) const;

        //- The sorted names of objects with a matching class name
        //- that also have a matching object name.
        //  \note uses the class type() method
        template<class MatchPredicate1, class MatchPredicate2>
        wordList sortedNames
        (
            const MatchPredicate1& matchClass,
            const MatchPredicate2& matchName
        ) const;

        //- The sorted names of objects with a class satisfying \c isA\<Type\>
        //
        //  \note If \a Type is \c void, no isA check is used (always true).
        template<class Type>
        wordList sortedNames() const;

        //- The sorted names of objects with a class satisfying \c isA\<Type\>
        //- that also have a matching object name.
        //
        //  \note If \a Type is \c void, no isA check is used (always true).
        template<class Type, class MatchPredicate>
        wordList sortedNames(const MatchPredicate& matchName) const;


    // Lookup

        //- Lookup and return a const sub-objectRegistry.
        //
        //  \param forceCreate create it if it does not exist.
        //  \param recursive search parent registries.
        const objectRegistry& subRegistry
        (
            const word& name,
            const bool forceCreate = false,
            const bool recursive = false
        ) const;


        //- Return all objects with a class satisfying
        //- \c isA\<Type\> or \c isType\<Type\> (with Strict)
        template<class Type, bool Strict=false>
        HashTable<const Type*> lookupClass() const;

        //- Return all objects with a class satisfying
        //- \c isA\<Type\> or \c isType\<Type\> (with Strict)
        template<class Type, bool Strict=false>
        HashTable<Type*> lookupClass();

        //- Return all objects with a class satisfying \c isA\<Type\>
        //
        //  \param strict use \c isType\<Type\> instead of \c isA\<Type\>
        template<class Type>
        HashTable<const Type*> lookupClass(const bool strict) const;

        //- Return all objects with a class satisfying \c isA\<Type\>
        //
        //  \param strict use \c isType\<Type\> instead of \c isA\<Type\>
        template<class Type>
        HashTable<Type*> lookupClass(const bool strict);

        //- Return const pointer to the regIOobject.
        //
        //  \param recursive search parent registries
        //
        //  \return nullptr if the object was not found.
        const regIOobject* cfindIOobject
        (
            const word& name,
            const bool recursive = false
        ) const;

        //- Does the registry contain the regIOobject object (by name).
        //
        //  \param name the object name
        //  \param recursive search parent registries
        bool contains(const word& name, const bool recursive = false) const;

        //- Is the named Type found?
        //
        //  \param recursive search parent registries
        template<class Type>
        bool foundObject
        (
            const word& name,
            const bool recursive = false
        ) const;

        //- Return const pointer to the object of the given Type.
        //
        //  \param recursive search parent registries
        //
        //  \return nullptr if the object was not found or had incorrect type.
        template<class Type>
        const Type* cfindObject
        (
            const word& name,
            const bool recursive = false
        ) const;

        //- Return const pointer to the object of the given Type.
        //
        //  \param recursive search parent registries
        //
        //  \return nullptr if the object was not found or had incorrect type.
        template<class Type>
        const Type* findObject
        (
            const word& name,
            const bool recursive = false
        ) const;

        //- Return non-const pointer to the object of the given Type.
        //
        //  \param recursive search parent registries
        //
        //  \return nullptr if the object was not found or had incorrect type.
        template<class Type>
        Type* findObject
        (
            const word& name,
            const bool recursive = false
        );

        //- Return non-const pointer to the object of the given Type,
        //- using a const-cast to have it behave like a mutable.
        //  Exercise caution when using.
        //
        //  \param recursive search parent registries.
        //
        //  \return nullptr if the object was not found or had incorrect type.
        template<class Type>
        Type* getObjectPtr
        (
            const word& name,
            const bool recursive = false
        ) const;

        //- Lookup and return const reference to the object
        //- of the given Type. Fatal if not found or the wrong type.
        //
        //  \param recursive search parent registries.
        template<class Type>
        const Type& lookupObject
        (
            const word& name,
            const bool recursive = false
        ) const;

        //- Lookup and return non-const reference to the object
        //- of the given Type. Fatal if not found or the wrong type.
        //
        //  \param recursive search parent registries.
        template<class Type>
        Type& lookupObjectRef
        (
            const word& name,
            const bool recursive = false
        ) const;


    // Events

        //- Return new event number.
        label getEvent() const;


    // Edit

        //- Add a regIOobject to registry. A nullptr is ignored.
        bool checkIn(regIOobject* io) const;

        //- Add a regIOobject to registry
        bool checkIn(regIOobject& io) const;

        //- Remove a regIOobject from registry and free memory if the
        //- object is ownedByRegistry. A nullptr is ignored.
        bool checkOut(regIOobject* io) const;

        //- Remove a regIOobject from registry and free memory if the
        //- object is ownedByRegistry.
        bool checkOut(regIOobject& io) const;

        //- Remove a regIOobject by name from registry and free memory if the
        //- object is ownedByRegistry
        bool checkOut(const word& key) const;

        //- Clear all entries from the registry
        //  Performs a checkOut() for all objects that are ownedByRegistry
        void clear();

        //- Clear all entries from the registry and the table itself.
        void clearStorage();

        //- Erase an entry specified by the given iterator.
        //  Performs a checkOut() if the object was ownedByRegistry.
        //  \return True if the entry existed and was removed
        bool erase(const iterator& iter);

        //- Erase an entry specified by the given key
        //  Performs a checkOut() if the object was ownedByRegistry.
        //  \return True if the entry existed and was removed
        bool erase(const word& key);

        //- Remove entries given by the listed keys
        //  Performs a checkOut() for all objects that are ownedByRegistry.
        //  \return The number of items removed
        label erase(std::initializer_list<word> keys);

        //- Remove entries given by the listed keys
        //  Performs a checkOut() for all objects that are ownedByRegistry.
        //  \return The number of items removed
        label erase(const UList<word>& keys);

        //- Rename
        virtual void rename(const word& newName);


    // Temporaries

        //FUTURE //- Add given name to the set of temporary objects to cache
        //FUTURE void addTemporaryObject(const word& name) const;

        //- True if given name is in the cacheTemporaryObjects set
        bool is_cacheTemporaryObject(const word& name) const;

        //- True if name of object is in the cacheTemporaryObjects set
        bool is_cacheTemporaryObject(const regIOobject* io) const;

        //- True if name of object is in the cacheTemporaryObjects set
        bool is_cacheTemporaryObject(const regIOobject& io) const;

        //- Cache the given object. Moves content and stores
        template<class Type>
        bool cacheTemporaryObject(Type& obj) const;

        //- Reset the cache state of the given object (nullptr is ignored)
        void resetCacheTemporaryObject(const regIOobject* io) const;

        //- Reset the cache state of the given object
        //- in the cacheTemporaryObjects set
        void resetCacheTemporaryObject(const regIOobject& io) const;

        //- Check that all objects specified in the cacheTemporaryObjects
        //- were also cached
        bool checkCacheTemporaryObjects() const;


    // Reading

        //- Return true if any of the object's files have been modified
        virtual bool modified() const;

        //- Read the objects that have been modified
        void readModifiedObjects();

        //- Read object if modified
        virtual bool readIfModified();


    // Writing

        //- The writeData function is required by regIOobject but not used.
        //  For this class, write is used instead
        virtual bool writeData(Ostream&) const
        {
            NotImplemented;
            return false;
        }

        //- Write the objects using stream options
        virtual bool writeObject
        (
            IOstreamOption streamOpt,
            const bool writeOnProc
        ) const;


    // Housekeeping

        //- Same as contains()
        bool found(const word& name, bool recursive = false) const
        {
            return this->contains(name, recursive);
        }

        //- Deprecated(2018-10) find object
        //  \deprecated(2018-10) - use findObject() method
        template<class Type>
        FOAM_DEPRECATED_FOR(2018-10, "findObject / cfindObject() methods")
        const Type* lookupObjectPtr
        (
            const word& name,
            bool recursive = false
        ) const
        {
            return this->cfindObject<Type>(name, recursive);
        }

        //- Deprecated(2018-10) get object pointer, ignoring constness
        //  \deprecated(2018-10) - use getObjectPtr() method
        template<class Type>
        FOAM_DEPRECATED_FOR(2018-10, "getObjectPtr() method")
        Type* lookupObjectRefPtr
        (
            const word& name,
            bool recursive = false
        ) const
        {
            return this->getObjectPtr<Type>(name, recursive);
        }

        //- Deprecated(2023-07) use csorted() method
        //  \deprecated(2023-07) - use csorted() method
        template<class Type>
        FOAM_DEPRECATED_FOR(2023-07, "csorted() method")
        UPtrList<const Type> sorted() const
        {
            return csorted<Type>();
        }

        //- Deprecated(2023-07) use csorted() method
        //  \deprecated(2023-07) - use csorted() method
        FOAM_DEPRECATED_FOR(2023-07, "csorted() method")
        UPtrList<const regIOobject> sorted() const
        {
            return csorted<regIOobject>();
        }

        //- Deprecated(2023-07) use csorted() method
        //  \deprecated(2023-07) - use csorted() method
        template<class Type, class MatchPredicate>
        FOAM_DEPRECATED_FOR(2023-07, "csorted() method")
        UPtrList<const Type> sorted(const MatchPredicate& matchName) const
        {
            return csorted<Type>(matchName);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "objectRegistryTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
