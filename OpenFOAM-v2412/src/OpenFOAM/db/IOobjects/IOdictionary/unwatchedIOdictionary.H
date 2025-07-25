/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2017 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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
    Foam::unwatchedIOdictionary

Description
    unwatchedIOdictionary is like IOdictionary but stores
    dependencies as files instead of fileMonitor watchIndices.
    Used to read controlDict since there fileMonitor not yet setup.

SourceFiles
    unwatchedIOdictionary.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_unwatchedIOdictionary_H
#define Foam_unwatchedIOdictionary_H

#include "baseIOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class unwatchedIOdictionary Declaration
\*---------------------------------------------------------------------------*/

class unwatchedIOdictionary
:
    public baseIOdictionary
{
    // Private Data

        //- The files that would be watched
        fileNameList files_;


public:

    // Constructors

        //- Construct given an IOobject
        //- and optional fallback dictionary content
        //  A null dictionary pointer is treated like an empty dictionary.
        explicit unwatchedIOdictionary
        (
            const IOobject& io,
            const dictionary* fallback = nullptr
        );

        //- Construct given an IOobject and fallback dictionary content
        unwatchedIOdictionary
        (
            const IOobject& io,
            const dictionary& dict
        );

        //- Construct given an IOobject, wanted typeName
        //- and optional fallback dictionary content
        //  A null dictionary pointer is treated like an empty dictionary.
        unwatchedIOdictionary
        (
            const IOobject& io,
            const word& wantedType,
            const dictionary* fallback = nullptr
        );

        //- Construct given an IOobject and Istream
        unwatchedIOdictionary(const IOobject& io, Istream& is);


    //- Destructor
    virtual ~unwatchedIOdictionary() = default;


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

        //- Add file watch on object (READ_IF_MODIFIED)
        virtual void addWatch();

        //- Add file watch for fileName on object if not yet watched.
        //  \return index of watch
        virtual label addWatch(const fileName&);

        //- Return the files that would be watched
        const fileNameList& files() const noexcept
        {
            return files_;
        }

        //- Access to the files that would be watched
        fileNameList& files() noexcept
        {
            return files_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Global file type for unwatchedIOdictionary
template<>
struct is_globalIOobject<unwatchedIOdictionary> : std::true_type {};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
