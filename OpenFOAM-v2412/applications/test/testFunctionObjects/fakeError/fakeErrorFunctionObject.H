/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::functionObjects::fakeErrorFunctionObject

Description
    A function object that emits FatalError at different stages
    (or does nothing), which is useful for testing purposes.

    Can request errors from constructor, execute and write methods.

Usage
    Example of function object specification:

    \verbatim
    test
    {
        type    fakeError;
        libs    (testFunctionObjects);
        ...
        errors          warn;
        executeError    true;
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property        | Description                         | Reqd | Deflt
        ioError         | FatalIOError instead of FatalError  | no  | no
        constructError  | generate error on construct         | no  | no
        executeError    | generate error on execute()         | no  | no
        writeError      | generate error on write()           | no  | no
    \endtable

SourceFiles
    fakeErrorFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_fakeErrorFunctionObject_H
#define functionObjects_fakeErrorFunctionObject_H

#include "timeFunctionObject.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                   Class fakeErrorFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class fakeErrorFunctionObject
:
    public timeFunctionObject
{
    // Private Data

        Switch ioError_;
        Switch constructError_;
        Switch executeError_;
        Switch writeError_;


    // Private Member Functions

        //- Emit error for whatever
        void emitError(const char* what) const;

        //- No copy construct
        fakeErrorFunctionObject(const fakeErrorFunctionObject&) = delete;

        //- No copy assignment
        void operator=(const fakeErrorFunctionObject&) = delete;


public:

    //- Runtime type information
    TypeName("fakeError");


    // Constructors

        //- Construct from Time and dictionary
        fakeErrorFunctionObject
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~fakeErrorFunctionObject() = default;


    // Member Functions

        //- Read and report settings
        virtual bool read(const dictionary& dict);

        //- Execute. Emit error or do nothing.
        virtual bool execute();

        //- Write. Emit error or do nothing.
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
