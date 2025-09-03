/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020-2021 OpenCFD Ltd.
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
    Foam::baseIOdictionary

Description
    baseIOdictionary is derived from dictionary and IOobject to give the
    dictionary automatic IO functionality via the objectRegistry.

    To facilitate IO, baseIOdictionary is provided with a constructor from
    IOobject and with readData/writeData functions.

SourceFiles
    baseIOdictionary.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_baseIOdictionary_H
#define Foam_baseIOdictionary_H

#include "dictionary.H"
#include "regIOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class baseIOdictionary Declaration
\*---------------------------------------------------------------------------*/

class baseIOdictionary
:
    public regIOobject,
    public dictionary
{
    // Private Data

        static bool writeDictionaries;


public:

    //- Declare type-name, virtual type (with debug switch)
    TypeName("dictionary");


    // Generated Methods

        //- Copy construct
        baseIOdictionary(const baseIOdictionary&) = default;

        //- Move construct
        baseIOdictionary(baseIOdictionary&&) = default;

        //- Destructor
        virtual ~baseIOdictionary() = default;


    // Constructors

        //- Construct given an IOobject
        //- and optional fallback dictionary content (ununsed)
        //  A null dictionary pointer is treated like an empty dictionary.
        explicit baseIOdictionary
        (
            const IOobject& io,
            const dictionary* fallback = nullptr
        );

        //- Construct given an IOobject and fallback dictionary
        //- content (ununsed)
        baseIOdictionary(const IOobject& io, const dictionary& dict);

        //- Construct given an IOobject and Istream (ununsed)
        baseIOdictionary(const IOobject& io, Istream& is);


    // Member Functions

        //- Return complete path + object name if the file exists
        //- either in the case/processor or case otherwise null
        virtual fileName filePath() const = 0;

        //- Name function is needed to disambiguate those inherited
        //- from regIOobject and dictionary
        const word& name() const;

        //- The readData function required by regIOobject read operation
        virtual bool readData(Istream&);

        //- The writeData function required by regIOobject write operation
        virtual bool writeData(Ostream&) const;

        //- Is object global
        virtual bool global() const = 0;


    // Member Operators

        //- Copy assignment of dictionary entries (leave regIOobject untouched)
        void operator=(const baseIOdictionary& rhs);

        //- Copy assignment of dictionary entries
        using dictionary::operator=;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
