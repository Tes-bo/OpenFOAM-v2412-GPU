/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2024 OpenCFD Ltd.
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
    Foam::coordinateSystems

Description
    A centralized collection of named coordinate systems.

Note
    Mixing normal constructors and the coordinateSystems::New constructor
    may yield unexpected results.

    \verbatim
    cat1
    {
        coordinateSystem
        {
            type   indirect;
            name   _10;
        }
        porosity    0.781;
        Darcy
        {
            d   d [0 -2 0 0 0]  (-1000 -1000 0.50753e+08);
            f   f [0 -1 0 0 0]  (-1000 -1000 12.83);
        }
    }
    \endverbatim

    For this to work correctly, the coordinateSystem constructor must be
    supplied with an objectRegistry as well as the dictionary.

SourceFiles
    coordinateSystems.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_coordinateSystems_H
#define Foam_coordinateSystems_H

#include "regIOobject.H"
#include "PtrList.H"
#include "coordinateSystem.H"
#include "wordRes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class coordinateSystems Declaration
\*---------------------------------------------------------------------------*/

class coordinateSystems
:
    public regIOobject,
    public PtrList<coordinateSystem>
{
    // Private Member Functions

        //- Read if IOobject flags set.
        //- Reads "coordinateSystems" or older "IOPtrList<coordinateSystem>"
        //  Return true if read.
        bool readIOcontents();

public:

    //- Declare type-name, virtual type (without debug switch)
    TypeNameNoDebug("coordinateSystems");


    // Generated Methods

        //- No copy construct
        coordinateSystems(const coordinateSystems&) = delete;

        //- No copy assignment
        void operator=(const coordinateSystems&) = delete;


    // Constructors

        //- Read construct from IOobject
        explicit coordinateSystems(const IOobject& io);

        //- Read construct "coordinateSystems" from "constant"
        //- using given registry
        explicit coordinateSystems(const objectRegistry& obr);

        //- Construct from IOobject and PtrList content
        coordinateSystems
        (
            const IOobject& io,
            const PtrList<coordinateSystem>& content
        );

        //- Construct from IOobject and transferring PtrList content
        coordinateSystems
        (
            const IOobject& io,
            PtrList<coordinateSystem>&& content
        );


    // Selectors

        //- Return previously registered or read construct from "constant"
        static const coordinateSystems& New(const objectRegistry& obr);


    // Member Functions

        //- Find and return indices for all matches
        //  A no-op (returns empty list) for an empty key
        labelList indices(const wordRe& key) const;

        //- Find and return indices for all matches
        //  A no-op (returns empty list) for an empty matcher
        labelList indices(const wordRes& matcher) const;

        //- Find and return index for the first match, return -1 if not found
        //  A no-op (returns -1) for an empty key
        label findIndex(const wordRe& key) const;

        //- Find and return index for the first match, return -1 if not found
        //  A no-op (returns -1) for an empty matcher
        label findIndex(const wordRes& matcher) const;

        //- Search if given key exists
        bool found(const wordRe& key) const;

        //- Return pointer to named coordinateSystem or nullptr on error
        const coordinateSystem* cfind(const word& name) const;

        //- Return reference to named coordinateSystem or FatalErrror
        const coordinateSystem& lookup(const word& name) const;

        //- A list of the coordinate-system names
        wordList names() const;

        //- A list of the coordinate-system names satisfying the input matcher
        wordList names(const wordRe& key) const;

        //- A list of the coordinate-system names satisfying the input matcher
        wordList names(const wordRes& matcher) const;

        //- Identical to names()
        inline wordList toc() const
        {
            return names();
        }


    // IO

        //- Write data
        bool writeData(Ostream& os) const;

        //- Write using stream options
        virtual bool writeObject
        (
            IOstreamOption streamOpt,
            const bool writeOnProc = true
        ) const;


    // Housekeeping

        //- Identical to the indices() method (AUG-2018)
        FOAM_DEPRECATED_FOR(2018-08, "indices() method")
        labelList findIndices(const wordRe& key) const
        {
            return this->indices(key);
        }

        //- Deprecated(2020-03) find named coordinateSystem or nullptr
        //
        //  \deprecated(2020-03) - use cfind() method
        FOAM_DEPRECATED_FOR(2020-03, "cfind() method")
        const coordinateSystem* lookupPtr(const word& name) const
        {
            return this->cfind(name);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Global file type for coordinateSystems - same content for all ranks
template<>
struct is_globalIOobject<coordinateSystems> : std::true_type {};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
