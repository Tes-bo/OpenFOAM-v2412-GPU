/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2021-2024 OpenCFD Ltd.
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
    Foam::IOdictionary

Description
    IOdictionary is derived from dictionary and IOobject to give the dictionary
    automatic IO functionality via the objectRegistry.  To facilitate IO,
    IOdictionary is provided with a constructor from IOobject and with
    readData/writeData functions.

SourceFiles
    IOdictionary.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_IOdictionary_H
#define Foam_IOdictionary_H

#include "baseIOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class IOdictionary Declaration
\*---------------------------------------------------------------------------*/

class IOdictionary
:
    public baseIOdictionary
{
public:

    // Constructors

        //- Construct given an IOobject
        //- and optional fallback dictionary content
        //  A null dictionary pointer is treated like an empty dictionary.
        explicit IOdictionary
        (
            const IOobject& io,
            const dictionary* fallback = nullptr
        );

        //- Construct given an IOobject
        //- and fallback dictionary content
        IOdictionary(const IOobject& io, const dictionary& dict);

        //- Construct given an IOobject, wanted typeName
        //- and optional fallback dictionary content
        //  A null dictionary pointer is treated like an empty dictionary.
        IOdictionary
        (
            const IOobject& io,
            const word& wantedType,
            const dictionary* fallback = nullptr
        );

        //- Construct given an IOobject and Istream
        IOdictionary(const IOobject& io, Istream& is);


    // Factory Methods

        //- Read and return contents, testing for "dictionary" type.
        //- The IOobject will not be registered
        static dictionary readContents(const IOobject& io);

        //- Read and return contents, testing for expected type.
        //- The IOobject will not be registered
        static dictionary readContents
        (
            const IOobject& io,
            const word& wantedType
        );


    //- Destructor
    virtual ~IOdictionary() = default;


    // Member Functions

        //- The object is global
        virtual bool global() const
        {
            return true;
        }

        //- Return complete path + object name if the file exists
        //- either in the case/processor or case otherwise null
        virtual fileName filePath() const
        {
            return globalFilePath(type());
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Global file type for IOdictionary
template<>
struct is_globalIOobject<IOdictionary> : std::true_type {};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
