/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2016-2022 OpenCFD Ltd.
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
    Foam::functionObjects::writeDictionary

Group
    grpUtilitiesFunctionObjects

Description
    Reports dictionary contents on change.

SourceFiles
    writeDictionary.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_functionObjects_writeDictionary_H
#define Foam_functionObjects_writeDictionary_H

#include "regionFunctionObject.H"
#include "wordList.H"
#include "SHA1Digest.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class writeDictionary Declaration
\*---------------------------------------------------------------------------*/

class writeDictionary
:
    public functionObjects::regionFunctionObject
{
    // Private Data

        //- Names of dictionaries to monitor
        wordList dictNames_;

        //- Digest of dictionary contents
        List<SHA1Digest> digests_;

        //- First time that checking has been called (per execution cycle)
        bool firstChange_;


    // Private Member Functions

        //- Write the output header
        void writeHeader();

        //- Helper to check and write the dictionary if its sha1 has changed
        void checkDictionary(const dictionary& dict, const label dicti);

        //- Attempt read dictionary at location
        bool tryDirectory(const word& location, const label dicti);

        //- Check all dictionaries
        bool performCheck();

        //- No copy construct
        writeDictionary(const writeDictionary&) = delete;

        //- No copy assignment
        void operator=(const writeDictionary&) = delete;


public:

    //- Runtime type information
    TypeName("writeDictionary");


    // Constructors

        //- Construct from Time and dictionary
        writeDictionary
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~writeDictionary() = default;


    // Member Functions

        //- Read the writeDictionary data
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the selected dictionaries
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
