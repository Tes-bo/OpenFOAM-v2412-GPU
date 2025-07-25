/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2022 OpenCFD Ltd.
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
    Foam::functionObjects::properties

Description
    Storage for function object properties, derived from IOdictionary.
    Provides functionality to read/write state information (data required for
    smooth restart behaviour) and results to/from the state dictionary

    Note: cannot be accessed until after construction of thefunction objects,
    since the owner container functionObjectList is owned by time, and time owns
    [this] i.e. need to wait for time to be fully constructed.

See also
    Foam::functionObject

SourceFiles
    functionObjectProperties.C
    functionObjectPropertiesTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_properties_H
#define functionObjects_properties_H

#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                         Class properties Declaration
\*---------------------------------------------------------------------------*/

class properties
:
    public IOdictionary
{
    // Private Member Data

        //- Name of the results dictionary
        static const word resultsName_;


protected:

    // Protected Member Functions

        //- No copy construct
        properties(const properties&) = delete;

        //- No copy assignment
        void operator=(const properties&) = delete;


public:

    // Constructors

        //- Construct from components
        explicit properties(const IOobject& io);


    //- Destructor
    virtual ~properties() = default;


    // Member Functions

        //- Return list of object names
        wordList objectNames() const;

        //- Return true if the object with objectName exists
        bool hasObjectDict(const word& objectName) const;

        //- Return access to the property dictionary
        dictionary& propertyDict(const word& objectName);


        // Properties

            //- Get dictionary for named object. Creates one if required
            dictionary& getObjectDict(const word& objectName);

            //- Return true if the property exists
            bool foundObjectProperty
            (
                const word& objectName,
                const word& entryName
            ) const;

            //- Remove the trigger index from the properties
            void clearTrigger();

            //- Get the current trigger index
            label getTrigger() const;

            //- Set new trigger index.
            //  \return True if the index changed
            bool setTrigger(const label triggeri);

            //- Set dictionary from named object, return true if set
            bool getObjectDict
            (
                const word& objectName,
                const word& entryName,
                dictionary& dict
            ) const;

            //- Retrieve generic property from named object
            template<class Type>
            Type getObjectProperty
            (
                const word& objectName,
                const word& entryName,
                const Type& defaultValue = Type(Zero)
            ) const;

            //- Set generic property from named object, return true if set
            template<class Type>
            bool getObjectProperty
            (
                const word& objectName,
                const word& entryName,
                Type& value
            ) const;

            //- Add generic property from named object
            template<class Type>
            void setObjectProperty
            (
                const word& objectName,
                const word& entryName,
                const Type& value
            );


        // Results

            //- Get dictionary of object results, return true if set
            bool getObjectResultDict
            (
                const word& objectName,
                dictionary& dict
            ) const;

            //- Add result from named object
            template<class Type>
            void setObjectResult
            (
                const word& objectName,
                const word& entryName,
                const Type& value
            );

            //- Retrieve result from named object
            template<class Type>
            Type getObjectResult
            (
                const word& objectName,
                const word& entryName,
                const Type& defaultValue = Type(Zero)
            ) const;

            //- Set result from named object, return true if set
            template<class Type>
            bool getObjectResult
            (
                const word& objectName,
                const word& entryName,
                Type& value
            ) const;

            //- Return true if the object with objectName exists in results
            bool hasResultObject(const word& objectName) const;

            //- Return list of objects with results
            wordList objectResultNames() const;

            //- Return true if the object with objectName exists and has
            //- entryName in its results
            bool hasResultObjectEntry
            (
                const word& objectName,
                const word& entryName
            ) const;

            //- Return the type of result
            word objectResultType
            (
                const word& objectName,
                const word& entryName
            ) const;

            //- Return result entries for named object
            wordList objectResultEntries(const word& objectName) const;

            //- Write the results entries for all objects to stream
            void writeResultEntries(Ostream& os) const;

            //- Write the results entries for named object to stream
            void writeResultEntries(const word& objectName, Ostream& os) const;

            //- Write the results entries for all objects to stream
            void writeAllResultEntries(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "functionObjectPropertiesTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
