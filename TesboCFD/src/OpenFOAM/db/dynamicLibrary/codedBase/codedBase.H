/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

Class
    Foam::codedBase

Description
    Base class for function objects and boundary conditions using dynamic code
    that provides methods for managing loading/unloading/updating
    of a dynamic library. For these purposes, it uses a dynamicCodeContext
    object to maintain information about the state.

    For simple coded objects, the default state management is sufficient.
    When there are more complicated code segments
    (eg, functionObjects::codedFunctionObject), the state management
    must also register these elements as well, starting with an initial
    setCodeContext() call and followed by append() to register each element.

SourceFiles
    codedBase.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_codedBase_H
#define Foam_codedBase_H

#include "dictionary.H"
#include "dynamicCodeContext.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class dynamicCode;
class dlLibraryTable;
class objectRegistry;

/*---------------------------------------------------------------------------*\
                       Class codedBase Declaration
\*---------------------------------------------------------------------------*/

class codedBase
{
    // Private Data

        //- The code context
        dynamicCodeContext context_;

        //- Previously loaded library
        mutable fileName oldLibPath_;


    // Private Member Functions

        //- Load specified library and execute funcName(true)
        void* loadLibrary
        (
            const fileName& libPath,
            const std::string& funcName,
            const dynamicCodeContext& context
        ) const;

        //- Execute funcName(false) and unload specified library
        void unloadLibrary
        (
            const fileName& libPath,
            const std::string& funcName,
            const dynamicCodeContext& context
        ) const;

        //- Create library based on the dynamicCodeContext
        void createLibrary
        (
            dynamicCode& dynCode,
            const dynamicCodeContext& context
        ) const;


protected:

    // Protected Member Functions

        //- Write code-dictionary contents
        static void writeCodeDict(Ostream& os, const dictionary& dict);

        //- Return "codeDict" from objectRegistry or read from disk
        static const dictionary& codeDict
        (
            const objectRegistry& obr,
            const word& dictName = "codeDict"
        );


        //- Access to the dynamic code context
        dynamicCodeContext& codeContext()
        {
            return context_;
        }

        //- Set code context from a dictionary
        void setCodeContext(const dictionary& dict);

        //- Add content to SHA1 hashing
        void append(const std::string& str);


        //- Update library as required, using the given context
        void updateLibrary
        (
            const word& name,
            const dynamicCodeContext& context
        ) const;

        //- Update library as required, using the given code dictionary
        //- to use for the context
        void updateLibrary
        (
            const word& name,
            const dictionary& dict
        ) const;

        //- Update library as required, using the predefined context
        //- or use the codeDict() to generate one
        void updateLibrary(const word& name) const;


        //- Mutable access to the loaded dynamic libraries
        virtual dlLibraryTable& libs() const = 0;

        // Return a description (type + name) for the output
        virtual string description() const = 0;

        // Clear any redirected objects
        virtual void clearRedirect() const = 0;

        // Get the dictionary to initialize the codeContext
        virtual const dictionary& codeDict() const = 0;

        //- Adapt the context for the current object
        virtual void prepare
        (
            dynamicCode& dynCode,
            const dynamicCodeContext& context
        ) const = 0;


public:

    //- Runtime type information
    ClassName("codedBase");


    // Generated Methods

        //- No copy construct
        codedBase(const codedBase&) = delete;

        //- No copy assignment
        void operator=(const codedBase&) = delete;


    // Constructors

        //- Default construct
        codedBase() = default;


    //- Destructor
    virtual ~codedBase() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
