/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2015-2021 OpenCFD Ltd.
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
    Foam::functionObjects::stateFunctionObject

Description
    Base class for function objects, adding functionality to read/write state
    information (data required for smooth restart behaviour) and results
    to/from the state dictionary

    Note: cannot access the state dictionary until after construction of the
    function objects, since the owner container functionObjectList is owned
    by time, and time owns the state dictionary  i.e. need to wait for time
    to be fully constructed.

See also
    Foam::functionObject

SourceFiles
    stateFunctionObject.C
    stateFunctionObjectTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_stateFunctionObject_H
#define functionObjects_stateFunctionObject_H

#include "timeFunctionObject.H"
#include "functionObjectProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class IOdictionary;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                     Class stateFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class stateFunctionObject
:
    public functionObjects::timeFunctionObject
{
    // Private Member Data

        //- Name of the results dictionary
        static const word resultsName_;


protected:

    // Protected Member Functions

        //- Return a const reference to the state dictionary
        const functionObjects::properties& stateDict() const;

        //- Return non-const access to the state dictionary
        functionObjects::properties& stateDict();


        //- No copy construct
        stateFunctionObject(const stateFunctionObject&) = delete;

        //- No copy assignment
        void operator=(const stateFunctionObject&) = delete;


public:

    // Constructors

        //- Construct from components
        stateFunctionObject(const word& name, const Time& runTime);


    //- Destructor
    virtual ~stateFunctionObject() = default;


    // Member Functions

        //- Return access to the property dictionary
        dictionary& propertyDict();


        // Properties

            //- Return true if the property exists
            bool foundProperty(const word& entryName) const;

            //- Remove the trigger index from the properties
            void clearTrigger();

            //- Get the current trigger index
            label getTrigger() const;

            //- Set new trigger index.
            //  \return True if the index changed
            bool setTrigger(const label triggeri);

            //- Set dictionary, return true if set
            bool getDict
            (
                const word& entryName,
                dictionary& dict
            ) const;

            //- Set dictionary from named object, return true if set
            bool getObjectDict
            (
                const word& objectName,
                const word& entryName,
                dictionary& dict
            ) const;

            //- Retrieve generic property
            template<class Type>
            Type getProperty
            (
                const word& entryName,
                const Type& defaultValue = Type(Zero)
            ) const;

            //- Set generic property, return true if set
            template<class Type>
            bool getProperty(const word& entryName, Type& value) const;

            //- Add generic property
            template<class Type>
            void setProperty(const word& entryName, const Type& value);

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

            bool getObjectResultDict
            (
                const word& objectName,
                dictionary& dict
            ) const;

            //- Add result
            template<class Type>
            void setResult
            (
                const word& entryName,
                const Type& value
            );

            //- Add result from named object
            template<class Type>
            void setObjectResult
            (
                const word& objectName,
                const word& entryName,
                const Type& value
            );

            //- Retrieve result
            template<class Type>
            Type getResult
            (
                const word& entryName,
                const Type& defaultValue = Type(Zero)
            ) const;

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

            //- Retrieve the result type
            word resultType(const word& entryName) const;

            //- Return the type of result
            word objectResultType
            (
                const word& objectName,
                const word& entryName
            ) const;

            //- Retrieve the result entries
            wordList objectResultEntries() const;

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
    #include "stateFunctionObjectTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
