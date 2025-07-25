/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2017 OpenFOAM Foundation
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
    Foam::localIOdictionary

Description
    localIOdictionary is derived from IOdictionary but excludes parallel
    master reading.

SourceFiles
    localIOdictionary.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_localIOdictionary_H
#define Foam_localIOdictionary_H

#include "baseIOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class localIOdictionary Declaration
\*---------------------------------------------------------------------------*/

class localIOdictionary
:
    public baseIOdictionary
{
public:

    // Constructors

        //- Construct given an IOobject
        //- and optional fallback dictionary content
        //  A null dictionary pointer is treated like an empty dictionary.
        explicit localIOdictionary
        (
            const IOobject& io,
            const dictionary* fallback = nullptr
        );

        //- Construct given an IOobject and fallback dictionary content
        localIOdictionary(const IOobject& io, const dictionary& dict);

        //- Construct given an IOobject, wanted typeName
        //- and optional fallback dictionary content
        //  A null dictionary pointer is treated like an empty dictionary.
        localIOdictionary
        (
            const IOobject& io,
            const word& wantedType,
            const dictionary* fallback = nullptr
        );

        //- Construct given an IOobject and Istream
        localIOdictionary(const IOobject& io, Istream& is);


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
    virtual ~localIOdictionary() = default;


    // Member Functions

        //- The object is not global
        virtual bool global() const
        {
            return false;
        }

        //- Return complete path + object name if the file exists
        //- in the case otherwise null
        virtual fileName filePath() const
        {
            // Use default (local only) search strategy
            return localFilePath(type());
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
