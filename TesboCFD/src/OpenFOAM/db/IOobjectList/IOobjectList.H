/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::IOobjectList

Description
    List of IOobjects with searching and retrieving facilities.
    Implemented as a HashTable, so the various sorted methods should
    be used if traversing in parallel.

SourceFiles
    IOobjectList.C
    IOobjectListI.H
    IOobjectListTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_IOobjectList_H
#define Foam_IOobjectList_H

#include "HashPtrTable.H"
#include "HashSet.H"
#include "IOobject.H"
#include "UPtrList.H"
#include "wordRes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class IOobjectList Declaration
\*---------------------------------------------------------------------------*/

class IOobjectList
:
    public HashPtrTable<IOobject>
{
    // Private Member Functions

        //- Check consistency of names and their order on all processors
        //- (the input list is assumed to be pre-sorted).
        //
        //  With syncPar = true, check that object names are identical
        //  (content and order) on all processors. FatalError if not.
        static void checkNameOrder(const wordList& objectNames, bool syncPar);

        //- Check consistency of object names/order on all processors
        //- (the input list is assumed to be pre-sorted).
        //
        //  With syncPar = true, check that object names are identical
        //  (content and order) on all processors. FatalError if not.
        static void checkObjectOrder
        (
            const UPtrList<const IOobject>& objs,
            bool syncPar
        );

        //- Combine names from all processors and sort
        static void syncNames(wordList& objNames);

        //- Templated implementation for classes()
        template<class MatchPredicate>
        static HashTable<wordHashSet> classesImpl
        (
            const IOobjectList& list,
            const MatchPredicate& matchName
        );

        //- Templated implementation for count()
        //  The number of items with a matching class
        template<class MatchPredicate1, class MatchPredicate2>
        static label countImpl
        (
            const IOobjectList& list,
            const MatchPredicate1& matchClass,
            const MatchPredicate2& matchName
        );

        //- Templated implementation for count()
        template<class Type, class MatchPredicate>
        static label countTypeImpl
        (
            const IOobjectList& list,
            const MatchPredicate& matchName
        );

        //- Templated implementation for names(), sortedNames()
        template<class MatchPredicate1, class MatchPredicate2>
        static wordList namesImpl
        (
            const IOobjectList& list,
            const MatchPredicate1& matchClass,
            const MatchPredicate2& matchName,
            const bool doSort
        );

        //- Templated implementation for names(), sortedNames()
        template<class Type, class MatchPredicate>
        static wordList namesTypeImpl
        (
            const IOobjectList& list,
            const MatchPredicate& matchName,
            const bool doSort
        );

        //- Templated implementation for csorted()/sorted()
        template<class Type, class MatchPredicate>
        static UPtrList<const IOobject> objectsTypeImpl
        (
            const IOobjectList& list,
            const MatchPredicate& matchName,
            const bool doSort
        );

        //- Templated implementation for lookup()
        template<class MatchPredicate>
        static IOobjectList lookupImpl
        (
            const IOobjectList& list,
            const MatchPredicate& matchName
        );

        //- Templated implementation for lookupClass()
        template<class MatchPredicate1, class MatchPredicate2>
        static IOobjectList lookupClassImpl
        (
            const IOobjectList& list,
            const MatchPredicate1& matchClass,
            const MatchPredicate2& matchName
        );

        //- Templated implementation for lookupClass()
        template<class Type, class MatchPredicate>
        static IOobjectList lookupClassTypeImpl
        (
            const IOobjectList& list,
            const MatchPredicate& matchName
        );


public:

    // Constructors

        //- Default construct: empty without allocation (capacity=0).
        IOobjectList() noexcept = default;

        //- Construct empty without allocation (capacity=0)
        inline explicit IOobjectList(const Foam::zero) noexcept;

        //- Construct empty with initial table capacity
        inline explicit IOobjectList(const label initialCapacity);

        //- Copy construct
        inline IOobjectList(const IOobjectList& list);

        //- Move construct
        inline IOobjectList(IOobjectList&& list);

        //- Construct from registry, instance, io options
        inline IOobjectList
        (
            const objectRegistry& db,
            const fileName& instance,
            IOobjectOption ioOpt
        );

        //- Construct from registry, instance, local, io options
        IOobjectList
        (
            const objectRegistry& db,
            const fileName& instance,
            const fileName& local,
            IOobjectOption ioOpt
        );

        //- Construct from registry, instance, registration option
        //- with MUST_READ, NO_WRITE
        inline IOobjectList
        (
            const objectRegistry& db,
            const fileName& instance,
            IOobjectOption::registerOption = IOobjectOption::REGISTER
        );

        //- Construct from registry, instance, local, registration option
        //- with MUST_READ, NO_WRITE
        inline IOobjectList
        (
            const objectRegistry& db,
            const fileName& instance,
            const fileName& local,
            IOobjectOption::registerOption = IOobjectOption::REGISTER
        );


    //- Destructor
    ~IOobjectList() = default;


    // Member Functions

    // Basic methods

        //- Move insert IOobject into the list
        inline bool add(std::unique_ptr<IOobject>&& objectPtr);

        //- Move insert IOobject into the list
        inline bool add(autoPtr<IOobject>& objectPtr);

        //- Move insert IOobject into the list
        inline bool add(autoPtr<IOobject>&& objectPtr);

        //- Remove object from the list, by name or by iterator.
        //
        //  \return autoPtr<IOobject>
        using HashPtrTable<IOobject>::remove;

        //- Remove object from the list by its IOobject::name().
        //
        //  \return autoPtr<IOobject>
        autoPtr<IOobject> remove(const IOobject& io)
        {
            return remove(io.name());
        }


    // Lookup single item

        //- Return const pointer to the object found by name
        //  \return IOobject ptr if found else nullptr
        const IOobject* cfindObject(const word& objName) const;

        //- Return const pointer to the object found by name that also
        //- has headerClassName == Type::typeName
        //  \return IOobject ptr if found and the correct type, else nullptr
        //
        //  \note If \a Type is \c void, no headerClassName check is used
        //      (always true).
        template<class Type>
        const IOobject* cfindObject(const word& objName) const;

        //- Return const pointer to the object found by name
        //  \return IOobject ptr if found else nullptr
        const IOobject* findObject(const word& objName) const;

        //- Return const pointer to the object found by name that also
        //- has headerClassName == Type::typeName
        //  \return IOobject ptr if found and the correct type, else nullptr
        //
        //  \note If \a Type is \c void, no headerClassName check is used
        //      (always true).
        template<class Type>
        const IOobject* findObject(const word& objName) const;

        //- Return non-const pointer to the object found by name
        //  \return IOobject ptr if found else nullptr
        IOobject* findObject(const word& objName);

        //- Return non-const pointer to the object found by name that also
        //- has headerClassName == Type::typeName
        //  \return IOobject ptr if found and the correct type, else nullptr
        //
        //  \note If \a Type is \c void, no headerClassName check is used
        //      (always true).
        template<class Type>
        IOobject* findObject(const word& objName);

        //- Return non-const pointer to the object found by name,
        //- using a const-cast to have it behave like a mutable.
        //  Exercise caution when using.
        //  \return IOobject ptr if found else nullptr
        IOobject* getObject(const word& objName) const;

        //- Return non-const pointer to the object found by name that also
        //- has headerClassName == Type::typeName,
        //- using a const-cast to have it behave like a mutable.
        //  Exercise caution when using.
        //  \return IOobject ptr if found and the correct type, else nullptr
        //
        //  \note If \a Type is \c void, no headerClassName check is used
        //      (always true).
        template<class Type>
        IOobject* getObject(const word& objName) const;


    // Lookup multiple items

        //- The list of IOobjects that have a matching object name.
        template<class MatchPredicate>
        IOobjectList lookup(const MatchPredicate& matchName) const;

        //- The list of IOobjects with the given headerClassName
        IOobjectList lookupClass(const char* clsName) const;

        //- The list of IOobjects with matching headerClassName
        template<class MatchPredicate>
        IOobjectList lookupClass(const MatchPredicate& matchClass) const;

        //- The list of IOobjects with matching headerClassName
        //- that also have a matching object name.
        template<class MatchPredicate1, class MatchPredicate2>
        IOobjectList lookupClass
        (
            const MatchPredicate1& matchClass,
            const MatchPredicate2& matchName
        ) const;

        //- The list of IOobjects with headerClassName == Type::typeName
        //
        //  \note If \a Type is \c void, no headerClassName check is used
        //      (always true).
        template<class Type>
        IOobjectList lookupClass() const;

        //- The list of IOobjects with headerClassName == Type::typeName
        //- that also have a matching object name.
        template<class Type, class MatchPredicate>
        IOobjectList lookupClass(const MatchPredicate& matchName) const;


    // Summary of classes

        //- A summary hash of classes used and their associated object names.
        //  The HashTable representation allows us to leverage various
        //  HashTable methods.
        //  This hashed summary view can be useful when querying particular
        //  aspects. For example,
        //
        //  \code
        //  IOobjectList objects(runTime, runTime.timeName());
        //  HashTable<wordHashSet> classes = objects.classes();
        //
        //  // How many volScalarField?
        //  word checkType = volScalarField::typeName;
        //
        //  Info<< checkType << "="
        //      << (classes.contains(checkType) ? classes[checkType].size() : 0)
        //      << nl;
        //  \endcode
        //  Using the two-parameter HashTable::lookup method lets us avoid
        //  the \c '?' ternary, but still looks fairly ugly:
        //  \code
        //  Info<< checkType << "="
        //      << classes.lookup(checkType, wordHashSet()).size() << nl;
        //  \endcode
        //
        //  If we have non-const access to the hash table, and don't mind
        //  incidentally creating empty entries,
        //  we can use the HashTable::operator() directly:
        //  \code
        //  Info<< checkType << "=" << classes(checkType).size() << nl;
        //  \endcode
        //
        //  Of course, for a single query it would have been easier
        //  and simpler to have used a direct query of the names:
        //  \code
        //  Info<< checkType << "=" << objects.names(checkType).size() << nl;
        //  \endcode
        //
        //  The summary hash, however, becomes most useful when reducing
        //  the objects in consideration to a particular subset. For example,
        //  \code
        //  const wordHashSet interestingTypes
        //  {
        //      volScalarField::typeName,
        //      volVectorField::typeName
        //  };
        //
        //  classes.retain(interestingTypes);
        //  \endcode
        //  Or do just the opposite:
        //  \code
        //  classes.erase(unsupportedTypes);
        //  \endcode
        //  This also works with a hashedWordList, since it provides the
        //  expected '()' operator. But in this case the more general
        //  HashTable::filterKeys is required:
        //  \code
        //  const hashedWordList interestingTypes
        //  {
        //      volScalarField::typeName,
        //      volVectorField::typeName
        //  };
        //
        //  classes.filterKeys(interestingTypes);
        //  \endcode
        //
        //  Of course, there are many other ways to use and manipulate the
        //  summary information.
        //
        //  \return HashTable of class-names for keys and wordHashSet of
        //      of object-names for the values.
        HashTable<wordHashSet> classes() const;

        //- A summary hash of classes used and their associated object names,
        //- restricted to objects that have a matching object name.
        template<class MatchPredicate>
        HashTable<wordHashSet> classes(const MatchPredicate& matchName) const;


    // List-wise access (unsorted)

        //- The unsorted list of IOobjects with
        //- headerClassName == Type::typeName
        //
        //  \note If \a Type is \c void, no headerClassName check is used
        //      (always true).
        //  The lifetime of the returned content cannot exceed the parent!
        template<class Type>
        UPtrList<const IOobject> cobjects() const;

        //- The unsorted list of IOobjects with
        //- headerClassName == Type::typeName
        //- that also have a matching object name.
        //  The lifetime of the returned content cannot exceed the parent!
        template<class Type, class MatchPredicate>
        UPtrList<const IOobject> cobjects
        (
            const MatchPredicate& matchName
        ) const;

    // List-wise access (sorted)

        //- The sorted list of IOobjects with
        //- headerClassName == Type::typeName
        //
        //  \note If \a Type is \c void, no headerClassName check is used
        //      (always true).
        //  The lifetime of the returned content cannot exceed the parent!
        template<class Type>
        UPtrList<const IOobject> csorted() const;

        //- The sorted names of the IOobjects with optional check for
        //- parallel consistency.
        //  FatalError if syncPar = true and names are not consistent on all
        //  processors.
        //  The lifetime of the returned content cannot exceed the parent!
        template<class Type>
        UPtrList<const IOobject> csorted(const bool syncPar) const;

        //- The sorted list of IOobjects
        //  The lifetime of the returned content cannot exceed the parent!
        UPtrList<const IOobject> csorted() const
        {
            return csorted<void>();
        }

        //- The sorted list of IOobjects with optional check for
        //- parallel consistency.
        //  FatalError if syncPar = true and names are not consistent on all
        //  processors.
        //  The lifetime of the returned content cannot exceed the parent!
        UPtrList<const IOobject> csorted(const bool syncPar) const
        {
            return csorted<void>(syncPar);
        }

        //- The sorted list of IOobjects with headerClassName == Type::typeName
        //- that also have a matching object name.
        //  The lifetime of the returned content cannot exceed the parent!
        template<class Type, class MatchPredicate>
        UPtrList<const IOobject> csorted(const MatchPredicate& matchName) const;

        //- The sorted list of IOobjects with headerClassName == Type::typeName
        //- that also have a matching object name.
        //  FatalError if syncPar = true and names are not consistent on all
        //  processors.
        //  The lifetime of the returned content cannot exceed the parent!
        template<class Type, class MatchPredicate>
        UPtrList<const IOobject> csorted
        (
            const MatchPredicate& matchName,
            const bool syncPar
        ) const;


    // Number of items

        //- The number of objects of the given headerClassName
        //  \note uses the class type() method
        label count(const char* clsName) const;

        //- The number of objects of the given headerClassName
        template<class MatchPredicate>
        label count
        (
            const MatchPredicate& matchClass
        ) const;

        //- The number of objects of the given headerClassName
        //- that also have a matching object name.
        template<class MatchPredicate1, class MatchPredicate2>
        label count
        (
            const MatchPredicate1& matchClass,
            const MatchPredicate2& matchName
        ) const;

        //- The number of objects with headerClassName == Type::typeName
        template<class Type>
        label count() const;

        //- The number of objects with headerClassName == Type::typeName
        //- that also have a matching object name.
        //
        //  \note If \a Type is \c void, no headerClassName check is used
        //      (always true).
        template<class Type, class MatchPredicate>
        label count(const MatchPredicate& matchName) const;


    // Summary of names

        //- The unsorted names of the IOobjects
        wordList names() const;

        //- The sorted names of the IOobjects with optional check for
        //- parallel consistency.
        //  FatalError if syncPar = true and names are not consistent on all
        //  processors.
        //  \note Output is always sorted - for consistent serial/parallel
        //  behaviour.
        wordList names(const bool syncPar) const;

        //- The unsorted names of IOobjects with the given headerClassName
        wordList names(const char* clsName) const;

        //- The sorted names of the IOobjects with the given headerClassName
        //  FatalError if syncPar = true and names are not consistent on all
        //  processors.
        //  \note Output is always sorted - for consistent serial/parallel
        //  behaviour.
        wordList names(const char* clsName, const bool syncPar) const;

        //- The unsorted names of IOobjects with the given headerClassName
        template<class MatchPredicate>
        wordList names(const MatchPredicate& matchClass) const;

        //- The sorted names of the IOobjects with the given headerClassName
        //  FatalError if syncPar = true and names are not consistent on all
        //  processors.
        //  \note Output is always sorted - for consistent serial/parallel
        //  behaviour.
        template<class MatchPredicate>
        wordList names
        (
            const MatchPredicate& matchClass,
            const bool syncPar
        ) const;

        //- The unsorted names of IOobjects with the given headerClassName
        //- that also have a matching object name.
        template<class MatchPredicate1, class MatchPredicate2>
        wordList names
        (
            const MatchPredicate1& matchClass,
            const MatchPredicate2& matchName
        ) const;

        //- The sorted names of the IOobjects with the given headerClassName
        //- that also have a matching object name.
        //  FatalError if syncPar = true and names are not consistent on all
        //  processors.
        //  \note Output is always sorted - for consistent serial/parallel
        //  behaviour.
        template<class MatchPredicate1, class MatchPredicate2>
        wordList names
        (
            const MatchPredicate1& matchClass,
            const MatchPredicate2& matchName,
            const bool syncPar
        ) const;


        //- The unsorted names of objects with
        //- headerClassName == Type::typeName
        template<class Type>
        wordList names() const;

        //- The sorted names of objects with
        //- headerClassName == Type::typeName.
        //  FatalError if syncPar = true and names are not consistent on all
        //  processors.
        //  \note Output is always sorted - for consistent serial/parallel
        //  behaviour.
        template<class Type>
        wordList names(bool syncPar) const;

        //- The unsorted names of objects with
        //- headerClassName == Type::typeName and a matching object name.
        template<class Type, class MatchPredicate>
        wordList names(const MatchPredicate& matchName) const;

        //- The sorted names of objects with
        //- headerClassName == Type::typeName and a matching object name.
        //  FatalError if syncPar = true and names are not consistent on all
        //  processors.
        //  \note Output is always sorted - for consistent serial/parallel
        //  behaviour.
        template<class Type, class MatchPredicate>
        wordList names
        (
            const MatchPredicate& matchName,
            const bool syncPar
        ) const;


    // Summary of names (sorted)

        //- The sorted names of the IOobjects
        wordList sortedNames() const;

        //- The sorted names of the IOobjects with optional check for
        //- parallel consistency.
        //  FatalError if syncPar = true and names are not consistent on all
        //  processors.
        wordList sortedNames(const bool syncPar) const;

        //- The sorted names of IOobjects with the given headerClassName
        wordList sortedNames(const char* clsName) const;

        //- The sorted names of the IOobjects with the given headerClassName
        //  FatalError if syncPar = true and names are not consistent on all
        //  processors.
        wordList sortedNames(const char* clsName, const bool syncPar) const;

        //- The sorted names of IOobjects with the given headerClassName
        template<class MatchPredicate>
        wordList sortedNames(const MatchPredicate& matchClass) const;

        //- The sorted names of the IOobjects with the given headerClassName
        //  FatalError if syncPar = true and names are not consistent on all
        //  processors.
        template<class MatchPredicate>
        wordList sortedNames
        (
            const MatchPredicate& matchClass,
            const bool syncPar
        ) const;

        //- The sorted names of IOobjects with the given headerClassName
        //- that also have a matching object name.
        template<class MatchPredicate1, class MatchPredicate2>
        wordList sortedNames
        (
            const MatchPredicate1& matchClass,
            const MatchPredicate2& matchName
        ) const;

        //- The sorted names of the IOobjects with the given headerClassName
        //- that also have a matching object name.
        //  FatalError if syncPar = true and names are not consistent on all
        //  processors.
        template<class MatchPredicate1, class MatchPredicate2>
        wordList sortedNames
        (
            const MatchPredicate1& matchClass,
            const MatchPredicate2& matchName,
            const bool syncPar
        ) const;


        //- The sorted names of objects with headerClassName == Type::typeName
        template<class Type>
        wordList sortedNames() const;

        //- The sorted names of objects with headerClassName == Type::typeName
        //  FatalError if syncPar = true and names are not consistent on all
        //  processors.
        template<class Type>
        wordList sortedNames(bool syncPar) const;

        //- The sorted names of objects with headerClassName == Type::typeName
        //- that also have a matching object name.
        template<class Type, class MatchPredicate>
        wordList sortedNames(const MatchPredicate& matchName) const;

        //- The sorted names of objects with headerClassName == Type::typeName
        //- that also have a matching object name.
        //  FatalError if syncPar = true and names are not consistent on all
        //  processors.
        template<class Type, class MatchPredicate>
        wordList sortedNames
        (
            const MatchPredicate& matchName,
            const bool syncPar
        ) const;


    // Edit

        //- Filter to retain or prune given classes
        //  \return The number of items changed (removed)
        template<class UnaryPredicate>
        label filterClasses
        (
            const UnaryPredicate& pred,
            const bool pruning = false
        );

        //- Filter to retain or prune given object names
        //  \return The number of items changed (removed)
        template<class UnaryPredicate>
        label filterObjects
        (
            const UnaryPredicate& pred,
            const bool pruning = false
        );

        //- Remove objects with names ending with "_0" (restart fields)
        //  \return The number of items changed (removed)
        label prune_0();


    // Parallel

        //- The sorted names of all objects (synchronised across processors)
        wordList allNames() const;

        //- The sorted names of all objects (synchronised across processors)
        //- with headerClassName == Type::typeName
        template<class Type>
        wordList allNames() const;

        //- Verify that object names are synchronised across processors
        //  FatalError if the names are not consistent on all processors.
        void checkNames(const bool syncPar = true) const;


    // Member Operators

        //- No copy assignment
        void operator=(const IOobjectList&) = delete;

        //- Move assignment
        inline void operator=(IOobjectList&& list);


    // Housekeeping

        //- Construct from registry, instance, local, io options
        //  \deprecated(2022-11) prefer IOobjectOption or registerOption
        IOobjectList
        (
            const objectRegistry& db,
            const fileName& instance,
            const fileName& local,
            IOobjectOption::readOption rOpt,
            IOobjectOption::writeOption wOpt = IOobjectOption::NO_WRITE,
            bool registerObject = true
        )
        :
            IOobjectList
            (
                db,
                instance,
                local,
                IOobjectOption(rOpt, wOpt, registerObject)
            )
        {}

        //- Deprecated(2018-11) Locate an object by name (c-string).
        //- Disambiguated from multiple-lookup version by calling parameter.
        //  \deprecated(2018-11) use findObject() for non-ambiguous resolution
        IOobject* lookup(const char* objName) const
        {
            return getObject(objName);
        }

        //- Deprecated(2018-11) Locate an object by name (const word&).
        //- Disambiguated from multiple-lookup version by calling parameter.
        //  \deprecated(2018-11) use findObject() for non-ambiguous resolution
        IOobject* lookup(const word& objName) const
        {
            return getObject(objName);
        }


        //- Deprecated(2023-07) use csorted() method
        //  \deprecated(2023-07) - use csorted() method
        template<class Type>
        FOAM_DEPRECATED_FOR(2023-07, "csorted() method")
        UPtrList<const IOobject> sorted() const
        {
            return csorted<Type>();
        }

        //- Deprecated(2023-07) use csorted() method
        //  \deprecated(2023-07) - use csorted() method
        template<class Type>
        FOAM_DEPRECATED_FOR(2023-07, "csorted() method")
        UPtrList<const IOobject> sorted(const bool syncPar) const
        {
            return csorted<Type>(syncPar);
        }

        //- Deprecated(2023-07) use csorted() method
        //  \deprecated(2023-07) - use csorted() method
        FOAM_DEPRECATED_FOR(2023-07, "csorted() method")
        UPtrList<const IOobject> sorted() const
        {
            return csorted<void>();
        }

        //- Deprecated(2023-07) use csorted() method
        //  \deprecated(2023-07) - use csorted() method
        FOAM_DEPRECATED_FOR(2023-07, "csorted() method")
        UPtrList<const IOobject> sorted(const bool syncPar) const
        {
            return csorted<void>(syncPar);
        }

        //- Deprecated(2023-07) use csorted() method
        //  \deprecated(2023-07) - use csorted() method
        template<class Type, class MatchPredicate>
        FOAM_DEPRECATED_FOR(2023-07, "csorted() method")
        UPtrList<const IOobject> sorted(const MatchPredicate& matchName) const
        {
            return csorted<Type>(matchName);
        }

        //- Deprecated(2023-07) use csorted() method
        //  \deprecated(2023-07) - use csorted() method
        template<class Type, class MatchPredicate>
        FOAM_DEPRECATED_FOR(2023-07, "csorted() method")
        UPtrList<const IOobject> sorted
        (
            const MatchPredicate& matchName,
            const bool syncPar
        ) const
        {
            return csorted<Type>(matchName, syncPar);
        }
};


// Ostream Operator
Ostream& operator<<(Ostream& os, const IOobjectList& list);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "IOobjectListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "IOobjectListTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
