/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
    Copyright (C) 2016-2021 OpenCFD Ltd.
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
    Foam::functionObjects::regionFunctionObject

Description
    Specialization of Foam::functionObject for a region and providing a
    reference to the region Foam::objectRegistry.

    Also provides support for referencing an alternative objectRegistry
    that can hold fields. This may be used, for example, to access
    stored surfaces and fields.

    Dictionary controls
    \table
      Property    | Description                       | Type | Req'd | Dflt
      region      | Name of the mesh region           | word | no    | region0
      subRegion   | Name for alternative objectRegistry | word | no  | ""
    \endtable

See also
    Foam::functionObjects::stateFunctionObject

SourceFiles
    regionFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_regionFunctionObject_H
#define functionObjects_regionFunctionObject_H

#include "stateFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class objectRegistry;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                    Class regionFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class regionFunctionObject
:
    public stateFunctionObject
{
protected:

    // Protected Member Data

        //- Name for alternative object registry
        word subRegistryName_;

        //- Reference to the region objectRegistry
        const objectRegistry& obr_;

        //- Pointer to alternative (eg, sub-region) objectRegistry.
        //  If a sub-region is not in effect, this is a nullptr
        mutable const objectRegistry* obrPtr_;


    // Protected Member Functions

        //- The region or sub-region registry being used
        virtual const objectRegistry& obr() const;


        //- Find object (eg, a field) in the (sub) objectRegistry
        template<class ObjectType>
        bool foundObject(const word& fieldName) const;

        //- Return const pointer to the object (eg, a field) in the
        //- (sub) objectRegistry.
        //
        //  \return nullptr if the object was not found or had incorrect type.
        template<class ObjectType>
        const ObjectType* cfindObject(const word& fieldName) const;

        //- Return const pointer to the object (eg, a field) in the
        //- (sub) objectRegistry.
        //
        //  \return nullptr if the object was not found or had incorrect type.
        template<class ObjectType>
        const ObjectType* findObject(const word& fieldName) const;

        //- Return non-const pointer to the object of the given Type,
        //- (sub) objectRegistry.
        //
        //  \return nullptr if the object was not found or had incorrect type.
        template<class ObjectType>
        ObjectType* findObject(const word& fieldName);

        //- Return non-const pointer to the object of the given Type,
        //- using a const-cast to have it behave like a mutable.
        //  Exercise caution when using.
        //
        //  \return nullptr if the object was not found or had incorrect type.
        template<class ObjectType>
        ObjectType* getObjectPtr(const word& fieldName) const;

        //- Lookup and return object (eg, a field) from the (sub) objectRegistry
        template<class ObjectType>
        const ObjectType& lookupObject(const word& fieldName) const;

        //- Lookup and return object (eg, a field) from the (sub) objectRegistry
        template<class ObjectType>
        ObjectType& lookupObjectRef(const word& fieldName) const;

        //- Store the field in the (sub) objectRegistry under the given name
        //  Note: sets the fieldName to tfield().name() if not already set
        template<class ObjectType>
        bool store
        (
            word& fieldName,
            const tmp<ObjectType>& tfield,
            bool cacheable = false
        );

        //- Store the field in an optional objectRegistry under the given name
        template<class ObjectType>
        bool storeInDb
        (
            const word& fieldName,
            const tmp<ObjectType>& tfield,
            const objectRegistry& obr
        );

        //- Write field if present in the (sub) objectRegistry
        bool writeObject(const word& fieldName);

        //- Clear field from the (sub) objectRegistry if present
        bool clearObject(const word& fieldName);

        //- Clear fields from the (sub) objectRegistry if present
        void clearObjects(const wordList& objNames);


        //- No copy construct
        regionFunctionObject(const regionFunctionObject&) = delete;

        //- No copy assignment
        void operator=(const regionFunctionObject&) = delete;


public:

    //- Runtime type information
    TypeName("regionFunctionObject");


    // Constructors

        //- Construct from Time and dictionary.
        //  The region objectRegistry is looked-up runTime with the name
        //  looked-up from the dictionary (defaults to polyMesh::defaultRegion)
        regionFunctionObject
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Construct from the region objectRegistry and dictionary
        regionFunctionObject
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict
        );


    //- Destructor
    virtual ~regionFunctionObject() = default;


    // Member Functions

        //- Read optional controls
        virtual bool read(const dictionary& dict);


    // Housekeeping

        //- Deprecated(2018-10)
        //  \deprecated(2018-10) - use findObject() method
        template<class ObjectType>
        FOAM_DEPRECATED_FOR(2018-10, "findObject / cfindObject() methods")
        const ObjectType* lookupObjectPtr(const word& fieldName) const
        {
            return this->cfindObject<ObjectType>(fieldName);
        }

        //- Deprecated(2018-10)
        //  \deprecated(2018-10) - use getObjectPtr() method
        template<class ObjectType>
        FOAM_DEPRECATED_FOR(2018-10, "getObjectPtr() method")
        ObjectType* lookupObjectRefPtr(const word& fieldName) const
        {
            return this->getObjectPtr<ObjectType>(fieldName);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "regionFunctionObjectTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
